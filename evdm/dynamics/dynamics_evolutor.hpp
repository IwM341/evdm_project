#pragma once 
#include <grob/grid.hpp>
#include <grob/linear_interpolator.hpp>
#include "../body_potential.hpp"
#include "../form_factors.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>

#include "../traj_pool.hpp"
#include "../utils/quadtree.hpp"
#include "../utils/progress_bar.hpp"
#include "../dynamics/dynamics.hpp"
#include "../utils/discret_generator.hpp"
#include <list>
#include <algorithm>
#include "../utils/pool_allocator.hpp"
namespace evdm {

    template <typename ReducerFunc_t,typename Arg0,typename Arg1>
    inline constexpr auto arg_reducer(ReducerFunc_t && reducer,Arg0 const & arg0,Arg1 const & arg1){
        return reducer(arg0,arg1);
    }

    template <typename ReducerFunc_t,typename Arg0,typename...Args>
    inline constexpr auto arg_reducer(ReducerFunc_t && reducer,Arg0 const & arg0,Args const&... args){
        return arg_reducer(reducer,arg0,arg_reducer(reducer,args...));
    }
    
    template <typename Tuple_t>
    auto tuple_min( Tuple_t const& _Tp){
        return std::apply([](auto const&... args){
            return arg_reducer([](auto const & a,auto const & b){return std::min(a,b);},args...);
        },_Tp);
    }
    template <typename Tuple_t>
    auto tuple_max( Tuple_t const& _Tp){
        return std::apply([](auto const&... args){
            return arg_reducer([](auto const & a,auto const & b){return std::max(a,b);},args...);
        },_Tp);
    }
    template <typename Tuple_t>
    auto tuple_sum( Tuple_t const& _Tp){
        return std::apply([](auto const&... args){
            return arg_reducer([](auto const & a,auto const & b){return a+b;},args...);
        },_Tp);
    }

    
    /// @brief struct with form factor and kinematics
    struct ElementInfo_t {
		size_t A,Z;
		float mN;
		float mX;
		float dmX;
		ScatterEvent ff;
	};

	template <typename FormFactorConstructor_impl>
	std::vector<ElementInfo_t> make_element_info(
		FormFactorConstructor_impl const & impl,
		size_t ptype_in,size_t ptype_out
	){
		return impl.construct(ptype_in,ptype_out);
	}

	template <typename T>
	using VecFunc_t = grob::GridFunction<
		grob::linear_interpolator,grob::GridVector<T>, std::vector<T>
	>;

	/// @brief function cfrac{2x}{ln(\frac{1+x}{1-x})}
	template <typename T>
	T WLogLinear(T x) {
		x = std::abs(x);
		if (x >= (T)1) {
			return 0;
		}
		
		// Для малых x используем разложение в ряд для лучшей точности
		if (x < 0.1) {
			double x2 = x * x;
			constexpr T cf[5] = {1./3,4./45,44./945,428./14175,10196./467775}; 
			return 1 - x2*(cf[0] + x2*(cf[1] + x2*(cf[2] + x2*(cf[3] + x2*cf[4]))));
		}
		else {
			return (2*x)/std::log((T(1) + x) / (T(1) - x));
		}
	}

	/// @brief struct to sample theta from random value 
	template <typename T>
	struct ScatterInfoTheta {
		T full_prob;
		VecFunc_t<T> SampleTheta;
		ScatterInfoTheta(){}

		template <typename Allocator>
		ScatterInfoTheta(std::list<grob::Point<T, T>,Allocator>  m_prob)
		{
			using namespace grob::literals;
			{
				//delete intervals with zero probability
				auto AvP= [](auto p0,auto p1){
					auto Pav = (p0+p1)/2;
					if(Pav > 0){
						return Pav*WLogLinear( 2*(p0-p1)/Pav );
					} 
					else {
						return Pav;
					}
				};
				
				if(m_prob.size() < 2){
					full_prob = 0;
					return;
					throw std::runtime_error("ScatterInfoTheta: can't construct from empty list");
				}
				//integrate probability	
				{
					T p_sum = 0;
					auto it = m_prob.begin();
					grob::Point<T,T> &_it = *it;
					T p0 = _it[0_c];
					T x0 = _it[1_c];
					_it[0_c] = 0;
					++it;
					for (; it != m_prob.end();++it) {
						auto [x1,p1] = *it;
						p_sum += AvP(p0,p1)*(x1-x0);
						(*it)[1_c] = p_sum;
						p0 = p1;
						x0 = x1;
					}
					full_prob = p_sum;
					if(full_prob == 0){
						return;
					}
					for(auto & [X,CDF] : m_prob){
						CDF *= (1/p_sum);
					}
					m_prob.back()[1_c] = 1;
				}
				
				//remove duplicates
				{
					auto it = m_prob.begin();
					while (it != m_prob.end()) {
						auto next = it;
						++next;
						
						if (next != m_prob.end() && (*it)[1_c] == (*next)[1_c]) {
							m_prob.erase(next);
						} else {
							++it;
						}
					}
				}
				
			}
			if(m_prob.size()<2){
				full_prob =0;
				#ifndef NDEBUG
				throw std::runtime_error("m_prob.size()<2");
				#endif
				return;
			}
			SampleTheta.Grid.resize(m_prob.size());
			SampleTheta.Values.resize(m_prob.size());
			size_t i=0;
			for(auto [Xpoint,Cdf] : m_prob){
				SampleTheta.Values[i] = Xpoint;
				SampleTheta.Grid[i] = Cdf;
				++i;
			}
		}
		ScatterInfoTheta(T full_prob, VecFunc_t<T> SampleTheta) :
			full_prob(full_prob), SampleTheta(std::move(SampleTheta)) {}
		SERIALIZATOR_FUNCTION(PROPERTY_NAMES("full_prob", "SampleTheta"), PROPERTIES(full_prob, SampleTheta));
		DESERIALIZATOR_FUNCTION(ScatterInfoTheta, PROPERTY_NAMES("full_prob", "SampleTheta"), PROPERTY_TYPES(full_prob, SampleTheta));

	};

	/// @brief struct with info: to sample the element toscatter + 
	/// trajectory distribution for each element 
	template <typename T>
	struct ScatterInfoElements {
		T full_prob;
		size_t first_nonzero;
		size_t last_nonzero;
		std::vector<ScatterInfoTheta<T>> ThetaTrajs;
		std::vector<T> ProbsElements;
		
		///indexes: i0_lmax,i1_lmax,i0_rmin,i1_rmin,i0_rmax,i1_rmax
		typedef std::tuple<size_t,size_t,size_t,size_t,size_t,size_t> indexes_t;
		indexes_t indexes;

		inline constexpr size_t i0_lmax() const{
			return std::get<0>(indexes);
		}
		inline constexpr size_t i1_lmax() const{
			return std::get<1>(indexes);
		}
		inline constexpr size_t i0_rmin() const{
			return std::get<2>(indexes);
		}
		inline constexpr size_t i1_rmin() const{
			return std::get<3>(indexes);
		}
		inline constexpr size_t i0_rmax() const{
			return std::get<4>(indexes);
		}
		inline constexpr size_t i1_rmax() const{
			return std::get<5>(indexes);
		}
		
		ScatterInfoElements(){}
		ScatterInfoElements(
			std::vector<ScatterInfoTheta<T>> ThetaTrajs,
			indexes_t indexes
		):
			ThetaTrajs(std::move(ThetaTrajs)),
			indexes(indexes)
		{
			full_prob = std::accumulate(ThetaTrajs.begin(), ThetaTrajs.end(),T(0),
				[](const auto& x, const auto& y) {return x + y.full_prob; });
			{
				size_t i=0;
				for( auto const & th :  ThetaTrajs){
					if(th.full_prob > 0){
						first_nonzero = i;
						break;
					}
				}
			}
			// find last nonzero
			{
				for( size_t i=1;i<=ThetaTrajs.size();++i){
					size_t j = ThetaTrajs.size()-i;
					if(ThetaTrajs[j].full_prob > 0){
						last_nonzero = j;
						break;
					}
				}
			}
			//count ProbsElements
			if(full_prob == 0){
				return;
			}
			ProbsElements.resize(ThetaTrajs.size());
			ProbsElements[0] = ThetaTrajs[0].full_prob/full_prob;
			for(size_t i=1;i<ThetaTrajs.size();++i){
				ProbsElements[i] = ProbsElements[i-1] + ThetaTrajs[i].full_prob/full_prob;
			}
			ProbsElements.back() = 1;
		}
		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("ThetaTrajs","indexes" ), 
			PROPERTIES(ThetaTrajs, indexes)
		);
		DESERIALIZATOR_FUNCTION(
			ScatterInfoElements, 
			PROPERTY_NAMES("ThetaTrajs", "indexes" ), 
			PROPERTY_TYPES(ThetaTrajs, indexes)
		);
		inline int correct_index(int i) const{
			for(;i>=0;--i){
				if(ThetaTrajs[i].full_prob > 0){
					return i;
				}
			}
			return -1;
		}
		size_t gen_index(T xi) const{
			size_t index = evdm::IndexGenBigHelper::gen(
				ProbsElements.begin(), ProbsElements.end(), xi
			);
			if(ProbsElements[index].full_prob == 0) [[unlikely]]{
				while(index < first_nonzero){
					++index;
					if(ProbsElements[index].full_prob > 0){
						return index;
					}
				}
				while(index > last_nonzero){
					--index;
					if(ProbsElements[index].full_prob > 0){
						return index;
					}
				}
			}else{
				return index;
			}
		}
	};

	template <typename T>
	static T interpolate_log(T p0,T p1,T alpha){
		if(p0 <= 0 || p1 <= 0){
			return 0;
		} else {
			return std::pow(p0,1-alpha)*std::pow(p1,alpha);
		}
	}
	template <typename T>
	static T interpolate_log(T p00,T p01,T p10,T p11,T alpha,T beta){
		return interpolate_log(
			interpolate_log(p00,p10,alpha),
			interpolate_log(p01,p11,alpha),
			beta
		);
	}

};