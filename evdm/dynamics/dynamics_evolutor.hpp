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
#include <ranges>
#include <span>
namespace evdm {

	//template <typename T>
	//using PoolAllocator = std::allocator<T>;

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

	template <typename Tuple_t,typename Functor_t,std::size_t...I>
	inline constexpr auto tuple_map_impl(Tuple_t && _Tp,Functor_t && F,std::index_sequence<I...>){
		return std::make_tuple(F(std::get<I>(std::forward<Tuple_t >(_Tp)))...);
	}

	template <typename Tuple_t,typename Functor_t>
	inline constexpr auto tuple_map(Tuple_t && _Tp,Functor_t && F){
		return tuple_map_impl(std::forward<Tuple_t >(_Tp),std::forward<Functor_t>(F),
			std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple_t>>>{});	
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
		auto Elements = impl.construct(ptype_in, ptype_out);
		std::sort(
			Elements.begin(), Elements.end(),
			[](const auto& a, const auto& b) {
				if (a.A > b.A) { return true; };
				if (a.A < b.A) { return false; };
				if (a.Z > b.Z) { return true; };
				if (a.Z < b.Z) { return false; };
				return  true;
			}
		);
		return Elements;
	}

	template <typename T>
	using VecFunc_t = grob::GridFunction<
		grob::linear_interpolator,grob::GridVector<T>, std::vector<T>
	>;

	/// @brief function cfrac{2x}{ln(\frac{1+x}{1-x})}
	template <typename T>
	T WLogLinear(T x) {
		x = std::abs(x);
		if (!(x < (T)1)) {
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
	

	template <typename T>
	inline T interpolate_log(T p0,T p1,T alpha){
		if(p0 <= 0 || p1 <= 0){
			return 0;
		} else {
			return std::pow(p0,alpha)*std::pow(p1,1-alpha);
		}
	}
	template <typename T>
	inline T interpolate_lin(T p0,T p1,T alpha){
		return p0*alpha + p1*(1-alpha);
	}
	template <typename T>
	inline T interpolate_lin(T p00,T p01,T p10,T p11,T alpha,T beta){
		return interpolate_lin(
			interpolate_lin(p00,p10,alpha),
			interpolate_lin(p01,p11,alpha),
			beta
		);
	}
	template <typename T>
	inline T interpolate_log(T p00,T p01,T p10,T p11,T alpha,T beta){
		return interpolate_log(
			interpolate_log(p00,p10,alpha),
			interpolate_log(p01,p11,alpha),
			beta
		);
	}

	template <typename T>
	using corner_tuple = std::tuple<T,T,T,T>;

	struct NoExceptionHandler{
		inline void operator () () const{}
	};


	template <typename T>
	VecFunc_t<T> List_to_Sample(std::vector<grob::Point<T, T>>  m_prob){
		VecFunc_t<T> SampleTheta;
		VecFunc_t<T> NullFunc(std::vector<T>({0,1}),std::vector<T>({0,0}));
		using namespace grob::literals;
		{
			//delete intervals with zero probability
			auto AvP= [](auto p0,auto p1){
				auto Pav = (p0+p1)/2;
				if(Pav > 0){
					return Pav*WLogLinear( (p0-p1)/Pav );
				} 
				else {
					return Pav;
				}
			};
			
			if(m_prob.size() < 2){
				return NullFunc;
				throw std::runtime_error("ScatterInfoTheta: can't construct from empty list");
			}
			//integrate probability	
			{
				T p_sum = 0;
				auto it = m_prob.begin();
				grob::Point<T,T> &_it = *it;
				T p0 = _it[1_c];
				T x0 = _it[0_c];
				_it[1_c] = 0;
				++it;
				for (; it != m_prob.end();++it) {
					auto [x1,p1] = *it;
					p_sum += AvP(p0,p1)*(x1-x0);
					(*it)[1_c] = p_sum;
					p0 = p1;
					x0 = x1;
				}
				if(p_sum == 0){
					return NullFunc;
				}
				for(auto & [X,CDF] : m_prob){
					CDF *= (1/p_sum);
				}
				m_prob.back()[1_c] = 1;
			}
			
		}
		//remove duplicates
		size_t input_i = 0;
		{
			size_t output_i = 1;
			for (; output_i < m_prob.size(); ++output_i) {
				if (m_prob[output_i] != m_prob[input_i]) [[likely]] {
					input_i++;
					if (input_i != output_i) [[unlikely]]
						m_prob[input_i] = m_prob[output_i];
				}
			}
		}
		std::span<grob::Point<T, T>> m_prob_slice(
			m_prob.begin(),
			input_i + 1
		);

		if(m_prob_slice.size()<2){
			#ifndef NDEBUG
			throw std::runtime_error("m_prob.size()<2");
			#endif
			return NullFunc;
		}
		
		size_t i=0;
		auto Xpoints = std::views::transform(m_prob_slice,
			[](const auto& P) {return P[0_c]; });
		auto Cdfs = std::views::transform(m_prob_slice,
			[](const auto& P) {return P[1_c]; });

		SampleTheta.Grid = std::vector<T>(Cdfs.begin(), Cdfs.end());
		SampleTheta.Values = std::vector<T>(Xpoints.begin(), Xpoints.end());
		return SampleTheta;
	}
	
	/// @brief struct to sample theta from random value 
	template <typename T>
	struct ScatterInfoTheta {
		T full_prob;
		T Tin;
		//VecFunc_t<T> SampleTheta;
		VecFunc_t<T> LogThetaProbs; // x - dimentionless theta, y - density
		ScatterInfoTheta():full_prob(0),Tin(0){}

		ScatterInfoTheta(std::span<grob::Point<T, T>>  m_prob,T Tin):Tin(Tin)
		{	
			using namespace grob::literals;
			std::vector<grob::Point<T, T>> m_reserve;
			if(m_prob.size() < 2){
				m_reserve = { {0,0},{0,1} };
				m_prob = m_reserve;
			}
			auto Xp = std::views::transform(m_prob,[](auto const & P){return P[0_c];});
			auto Yp = std::views::transform(m_prob,[](auto const & P){return std::log(P[1_c]);});

			LogThetaProbs.Grid = std::vector<T>(Xp.begin(), Xp.end());
			LogThetaProbs.Values = std::vector<T>(Yp.begin(),Yp.end());
			T sum = 0;
			const auto & X = LogThetaProbs.Grid;
			const auto & Y = LogThetaProbs.Values;
			for(size_t i=0;i< LogThetaProbs.size()-1;++i){
				sum += (X[i+1]- X[i]) * std::exp( (Y[i]+Y[i+1])/2 );
			}
			full_prob = sum;
		}
		inline T log_prob(T x)const{
			if(x >= LogThetaProbs.Grid.front() && x < LogThetaProbs.Grid.back())
				return LogThetaProbs(x);
			else
				return std::numeric_limits<T>::max();
		}

		inline T prob(T x)const{
			if(x >= LogThetaProbs.Grid.front() && x < LogThetaProbs.Grid.back())
				return std::exp( LogThetaProbs(x) );
			else
				return 0;
		}
		ScatterInfoTheta(T full_prob, T Tin,VecFunc_t<T> LogThetaProbs) :
			full_prob(full_prob),Tin(Tin), LogThetaProbs(std::move(LogThetaProbs)) {}
		SERIALIZATOR_FUNCTION(PROPERTY_NAMES("full_prob","Tin", "LogThetaProbs"), PROPERTIES(full_prob,Tin, LogThetaProbs));
		DESERIALIZATOR_FUNCTION(ScatterInfoTheta, PROPERTY_NAMES("full_prob","Tin", "LogThetaProbs"), PROPERTY_TYPES(full_prob,Tin, LogThetaProbs));

	};

	/// @brief struct to sample theta
	template <typename T>
	struct ScatterSampleTheta{
		VecFunc_t<T> SampleTheta;
		VecFunc_t<T> MajorantPDF;
		ScatterSampleTheta(){}
		
		template <typename X>
		using array4 = std::array<X,4>;
		
		inline grob::Point<T,T> sample(T xi)const{
			T x = SampleTheta(xi);
			return {x,MajorantPDF(x)};
		}

		ScatterSampleTheta(
			VecFunc_t<T> const & logP0,
			VecFunc_t<T> const & logP1,
			VecFunc_t<T> const & logP2,
			VecFunc_t<T> const & logP3)
		{
			using namespace grob::literals;
			size_t m_size = logP0.size()+ logP1.size()+ logP2.size()+ logP3.size();
			std::vector<grob::Point<T,T> > m_points;
			m_points.reserve(m_size);
			constexpr T max_val = std::numeric_limits<T>::max();
			
			VecFunc_t<T> null_prob(std::vector<T>({0,1}),std::vector<T>({ max_val,max_val}));
			std::array<std::reference_wrapper<const VecFunc_t<T>> ,4> Ps = { logP0,logP1,logP2,logP3};
			for(size_t i=0;i>4;++i){
				const VecFunc_t<T> & F = Ps[i];
				if(F.size() < 2){
					Ps[i] = null_prob;
				}
			}
			
			array4<size_t> m_index = {0,0,0,0};
			array4<size_t> m_end = grob::map(Ps,[](const VecFunc_t<T> & P){return P.size();});

			while(m_index != m_end){
				T min_next_x = 1;
				for(size_t j=0;j<4;++j){
					if(m_index[j] != m_end[j]){
						min_next_x = std::min(min_next_x,Ps[j].get().Grid[ m_index[j] ]);
					}
				}
				T max_value = 0;
				for(size_t j=0;j<4;++j){
					if(m_index[j] != m_end[j]) {
						if(Ps[j].get().Grid[ m_index[j] ] == min_next_x){
							max_value = std::max(max_value,std::exp(Ps[j].get().Values[ m_index[j] ]) );
							++m_index[j];
						} else {
							const VecFunc_t<T> & F = Ps[j];
							max_value = std::max(max_value,std::exp(F(min_next_x)));
						}
					}
				}
				m_points.push_back({min_next_x,max_value});
			}
			auto Xp = std::views::transform(m_points,[](auto const & P){return P[0_c];});
			auto Yp = std::views::transform(m_points,[](auto const & P){return P[1_c];});
			MajorantPDF.Grid = std::vector<T>(Xp.begin(),Xp.end());
			MajorantPDF.Values = std::vector<T>(Yp.begin(),Yp.end());
			SampleTheta = List_to_Sample(std::move(m_points));
		}

		template <typename Gen_t,typename dP_dth_functor_t>
		T gen_theta(dP_dth_functor_t ProbFunc,
			Gen_t G,size_t attempts,
			std::span<T> buff) const
		{
			constexpr size_t Nt = 80;

			for (size_t _a = 0; _a < attempts;++_a){
				T xi_prob = G();
				auto [m_th,majorant] = sample(xi_prob);
				T maj_xi = G();
				if(majorant * maj_xi <= ProbFunc(m_th)){
					return m_th;
				}
			}
			//no return: make by honest integration
			T xi_prob = G();
			T max_theta_undim = 1;
			T prob_func_0 = ProbFunc(0);
			for(size_t i=0;i<20;++i){
				if(ProbFunc(max_theta_undim/(Nt-1)) < prob_func_0 * 1e-10){
					max_theta_undim/=2;
				}
			}
			std::vector<T> m_bf;
			
			if(buff.size() < Nt){
				m_bf.resize(Nt);
				buff = m_bf;
			} else {
				buff = std::span{buff.begin(),Nt};
			}
			buff[0] = 0;
			T last_prob = prob_func_0;
			for(size_t i=1;i<buff.size();++i){
				T next_prob = ProbFunc(i*max_theta_undim/(Nt-1));
				T sum_prob = (last_prob+next_prob);
				T diff_prob = (last_prob-next_prob);
				buff[i] = buff[i-1] + (sum_prob)*WLogLinear(diff_prob/sum_prob)/2;
				last_prob = next_prob;
			}
			T target_prob = buff.back()*xi_prob;
			auto i = IndexGenBigHelper::gen(buff.begin(),buff.end(),target_prob);
			auto th0 = i*max_theta_undim/(Nt-1);
			auto dth = max_theta_undim/(Nt - 1);
			T alpha = upbound( (buff[i+1] - target_prob)/(buff[i+1]-buff[i]),1);
			return th0 + dth*(1-alpha);
		}

		ScatterSampleTheta(VecFunc_t<T> SampleTheta, VecFunc_t<T> MajorantPDF) :
			SampleTheta(std::move(SampleTheta)), MajorantPDF(std::move(MajorantPDF)) {}
		SERIALIZATOR_FUNCTION(PROPERTY_NAMES("SampleTheta", "MajorantPDF"), PROPERTIES(SampleTheta, MajorantPDF));
		DESERIALIZATOR_FUNCTION(ScatterSampleTheta, PROPERTY_NAMES("SampleTheta", "MajorantPDF"), PROPERTY_TYPES(SampleTheta, MajorantPDF));
	};

	/// @brief struct with info: to sample the element toscatter + 
	/// trajectory distribution for each element 
	template <typename T>
	struct ScatterInfoElements {
		T full_prob;
		std::vector<ScatterInfoTheta<T>> ThetaTrajs;
		//std::vector<ScatterSampleTheta<T>> Samples;
		std::vector<T> ProbsElements;
		
		///indexes: i0_lmax,i1_lmax,i0_rmin,i1_rmin,i0_rmax,i1_rmax
		typedef std::tuple<size_t,size_t,size_t> indexes_t;
		indexes_t indexes;

		inline constexpr size_t i_lmax() const{
			return std::get<0>(indexes);
		}
		inline constexpr size_t i_rmin() const{
			return std::get<1>(indexes);
		}
		inline constexpr size_t i_rmax() const{
			return std::get<2>(indexes);
		}
		
		ScatterInfoElements(){}
		ScatterInfoElements(
			std::vector<ScatterInfoTheta<T>> _ThetaTrajs,
			indexes_t indexes
		):
			ThetaTrajs(std::move(_ThetaTrajs)),
			indexes(indexes)
		{
			auto prob_getter = std::views::transform(ThetaTrajs,[](auto const & x){return x.full_prob;});
			full_prob = std::accumulate(prob_getter.begin(), prob_getter.end(),T(0));
			//count ProbsElements
			if(full_prob == 0){
				return;
			}
			ProbsElements = std::vector<T>(prob_getter.begin(),prob_getter.end());
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
		
	};
	/// @brief struct with info: to sample the element toscatter + 
	/// to sample theta
	template <typename T>
	struct ScatterSampleElements {
		std::vector<ScatterSampleTheta<T>> Samples;
		//std::vector<ScatterSampleTheta<T>> Samples;
		//std::vector<T> SamplesElements;
		
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
		
		ScatterSampleElements(){}
		
		using sci_cref_t = const ScatterInfoElements<T> & ;
		ScatterSampleElements( sci_cref_t If0,sci_cref_t If1,sci_cref_t If2,sci_cref_t If3){
			size_t N = If0.ProbsElements.size();
			Samples.reserve(N);
			for(size_t i=0;i< N;++i){
				auto mget = [i](auto const& If) {
					return If.ThetaTrajs[i].LogThetaProbs;
				};

				Samples.push_back(ScatterSampleTheta<T>(
					mget(If0), mget(If1), mget(If2),mget(If3)
				));
			}
			auto Ift = std::tie(If0,If1,If2,If3);
			
			indexes = std::make_tuple( 
				tuple_min(tuple_map(Ift,[](sci_cref_t X){return X.i_lmax();})),
				tuple_max(tuple_map(Ift,[](sci_cref_t X){return X.i_lmax() + 1;})),
				tuple_min(tuple_map(Ift,[](sci_cref_t X){return X.i_rmin();})),
				tuple_max(tuple_map(Ift,[](sci_cref_t X){return X.i_rmin()+1;})),
				tuple_min(tuple_map(Ift,[](sci_cref_t X){return X.i_rmax();})),
				tuple_max(tuple_map(Ift,[](sci_cref_t X){return X.i_rmax()+1;}))
			);
		}

		ScatterSampleElements(
			std::vector<ScatterSampleTheta<T>> Samples,
			//std::vector<T> SamplesElements,
			indexes_t indexes
		):
			Samples(std::move(Samples)),
			//SamplesElements(std::move(SamplesElements)),
			indexes(indexes)
		{
			
		}
		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("Samples",
				//"SamplesElements",
				"indexes"), 
			PROPERTIES(Samples,
				 //SamplesElements,
				 indexes)
		);
		DESERIALIZATOR_FUNCTION(
			ScatterSampleElements, 
			PROPERTY_NAMES("Samples", 
				//"SamplesElements",
				 "indexes"), 
			PROPERTY_TYPES(Samples,
				// SamplesElements,
				indexes)
		);

		template <typename ExceptionFunc_t = NoExceptionHandler>
		size_t gen_index(
			T xi,std::span<T> buffer,
			corner_tuple<ScatterInfoElements<T> const &> m_els,
			T alpha,T beta,ExceptionFunc_t && exception_null_prob = NoExceptionHandler{}) const{
			std::vector<T> M_buf;
			auto ProbFunc = [m_els,alpha,beta](size_t index){
				return interpolate_log( 
					std::get<0>(m_els).ProbsElements[index],
					std::get<1>(m_els).ProbsElements[index],
					std::get<2>(m_els).ProbsElements[index],
					std::get<3>(m_els).ProbsElements[index],
					alpha,beta
				);
			};
			size_t N_elements = std::get<0>(m_els).ProbsElements.size();
			if(buffer.size() < N_elements)[[unlikely]]{
				M_buf.resize(N_elements);
				buffer = M_buf;
			}
			else {
				buffer = std::span{ buffer.begin(),N_elements };
			}
			buffer[0] = ProbFunc(0);
			
			for(size_t i=1;i<N_elements;++i){
				buffer[i] = buffer[i-1] + ProbFunc(i);
			}

			size_t index = evdm::IndexGenBigHelper::gen(
				buffer.begin(), buffer.end(), xi*buffer.back()
			);

			if(buffer[index] == 0) [[unlikely]]{
				size_t i = index;
				while(i > 0){
					i--;
					if(buffer[i] > 0){
						return i;
					}
				}
				i = index;
				while(i < buffer.size()){
					++i;
					if(buffer[i] > 0){
						return i;
					}
				}
				exception_null_prob();
			}
			return index;
		}
	};

	

	

	template <typename T>
	struct PlotInfo{
		std::vector<T> X,Y,Values;
		std::vector<int> Indicies;
	};
	template <typename T, typename Vt, typename U, size_t N, typename ValueAdaptor>
	inline PlotInfo<T> quadtree_plotinfo_xyz(
		Quadtree<T, Vt, U, N> const& qt, ValueAdaptor&& F
	){
		size_t leafs = qt.leafs();
		size_t points = qt.point_num();
		size_t triangles = leafs * 2;
		std::vector<T> X(4 * leafs), Y(leafs * 4), Values(leafs * 4);
		std::vector<int> TriangleIndicies(3 * triangles);

		int i = 0;
		for (const QuadtreeNode<T, Vt, U>& Nd : qt) {
			// fill X and Y
			X[4 * i] = Nd.x0();
			X[4 * i + 1] = Nd.x0();
			X[4 * i + 2] = Nd.x1();
			X[4 * i + 3] = Nd.x1();

			Y[4 * i] = Nd.y0();
			Y[4 * i + 1] = Nd.y1();
			Y[4 * i + 2] = Nd.y0();
			Y[4 * i + 3] = Nd.y1();
			// fill Values
			Values[4 * i] = F(Nd.x0(), Nd.y0(), Nd(0, 0));
			Values[4 * i + 1] = F(Nd.x0(), Nd.y1(), Nd(0, 1));
			Values[4 * i + 2] = F(Nd.x1(), Nd.y0(), Nd(1, 0));
			Values[4 * i + 3] = F(Nd.x1(), Nd.y1(), Nd(1, 1));

			std::tie(
				TriangleIndicies[6 * i + 0],
				TriangleIndicies[6 * i + 1],
				TriangleIndicies[6 * i + 2]
			) = std::make_tuple(4 * i, 4 * i + 1, 4 * i + 2);
			std::tie(
				TriangleIndicies[6 * i + 3 + 0],
				TriangleIndicies[6 * i + 3 + 1],
				TriangleIndicies[6 * i + 3 + 2]
			) = std::make_tuple(4 * i + 1, 4 * i + 2, 4 * i + 3);
			++i;
		}
		return { X,Y,Values,TriangleIndicies };
	}

	template <typename T,typename Vt,typename U,size_t N,typename ValueAdaptor>
	PlotInfo<T> quadtree_plotinfo(Quadtree<T,Vt,U,N> const & qt,ValueAdaptor && F){
		return quadtree_plotinfo_xyz(
			qt, [&F](auto const& x, auto const& y, auto const& z) {
				return F(z);
			}
		);
	}


};