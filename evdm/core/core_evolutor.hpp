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
#include "../dynamics/dynamics_evolutor.hpp"

#if defined(__GNUC__) || defined(__clang__)
    #define PURE_FUNCTION [[gnu::pure]]
#elif defined(_MSC_VER)
    #define PURE_FUNCTION [[msvc::pure]]
#else
    #define PURE_FUNCTION
#endif
namespace evdm {


	template <typename T>
	struct scatter_pre_result{
		T cos_theta_in;
		T cos_theta_out;
		T gen_pre_form_factor;
	};
	template <typename T>
	inline PURE_FUNCTION scatter_pre_result<T> gen_scatter(
		T xi_theta_in,same_t<T> xi_theta_out,
		same_t<T> V_in,same_t<T> v1_in, same_t<T> u_in,
		same_t<T> DeltaV2){
		T A = (v1_in*v1_in + V_in*V_in + DeltaV2);
		T B = 2*V_in*v1_in;
		
		if(A > B){
			T t = B/A;
			//evaluate (1 - (1+x)^{2/3})/x
			auto Wfunc_2_3 = [](T one_plus_x){
				T q = std::pow(one_plus_x,T(1)/T(3));
				return -(1+q)/(1+q+q*q);
			};
			//evaluate ((1+x)^{3/2} - 1)/x
			auto Wfunc_3_2_p = [](T x){
				T x1 = 1+x;
				T q = std::pow(x1,T(3)/T(2));
				return (1 + x1 + x1*x1)/(1 + q*q*q);
			};
			//evaluate (1 - (1-x)^{3/2})/x
			auto Wfunc_3_2_m = [](T x){
				T x1 = 1-x;
				T q = std::pow(x1,T(3)/T(2));
				return (1 + x1 + x1*x1)/(1 + q*q*q);
			};


			// evaluate [ 1 - ( (1+t)^{3/2} - xi*t*{ ((1+t)^{3/2} - (1-t)^{3/2} )/t} )^{2/3}]/t 
			// =  A * [ (1-(1-A*t)^{2/3})/(At) ] = A*Wfunc_2_3(A*t)
			// 1 + A*t = (1+t)^{3/2} - xi*t*{ ((1+t)^{3/2} - (1-t)^{3/2} )/t}
			// A = Wfunc_3_2_p(t) - xi*(Wfunc_3_2_p(t)+Wfunc_3_2_m(t))
			T a = Wfunc_3_2_p(t) - xi_theta_in*(Wfunc_3_2_p(t)+Wfunc_3_2_m(t));
			T cos_theta_in = a * Wfunc_2_3(1+a*t);
			T cos_theta_out = -1 + 2*xi_theta_out;

			T m_factor = (Wfunc_3_2_p(t)+Wfunc_3_2_m(t))*(T(2)/T(3));
			//T G1 = (T(2)/T(3))*std::sqrt(A)*(Wfunc_3_2_p(t)+Wfunc_3_2_m(t));

			return {cos_theta_in,cos_theta_out, 
				m_factor*std::sqrt(A * u_in * std::numbers::inv_pi_v<T>)};
		}
		else if(B > A && (A+B)>=0 ){
			T t = A/B;

			T t_plus_1_pow_3_2 = std::pow(1+t,T(1.5));
			T cos_theta_in = t - std::pow( t_plus_1_pow_3_2*(1-xi_theta_in),T(2)/T(3) );
			T cos_theta_out = -1 + 2*xi_theta_out;
			T m_factor = t_plus_1_pow_3_2*(T(2)/T(3));
			return {
				cos_theta_in,cos_theta_out, 
				m_factor*std::sqrt(B * u_in * std::numbers::inv_pi_v<T>)};
		}
		else{
			return {
				-1 + 2*xi_theta_in,-1 + 2*xi_theta_out,\
				T(0)};
		}
	}


	template <typename T,typename FF_t>
	inline T get_ff_out(
		T Vu_in,T Vu_out,T mu,T delta,T VescMin,T cos_theta_out,
		FormFactor_t const & FF){
		T m_factor = mu*VescMin*VescMin;
		
		T delta_Vu2 = Vu_in*Vu_in+ Vu_out*Vu_out - 2*Vu_in*Vu_out*cos_theta_out;
		T t = delta_Vu2 >0 ? delta/(m_factor*delta_Vu2)	: T(0);
		
		T Vu_in_t = Vu_in *( T(0.5)-t);
		T Vu_out_t = Vu_out *( T(0.5)+t);

		T V2_perp = Vu_in_t*Vu_in_t + Vu_out_t*Vu_out_t - 2*Vu_in_t*Vu_out_t*cos_theta_out;
		
		return FF.ScatterFactor(
			m_factor*delta_Vu2,
			V2_perp,
			delta
		);
	}

	
	template <typename T>
	std::pair<T,T> get_el_out(
		T mu_xi,T mu_i,
		T r, T phi_r, T vr,T vtau,T v,T v1,T Vu_in,T Vu_out,
		T cos_theta_in,T phi_in,
		T cos_theta_out,T phi_out)
	{
		vec3<T> V_in(vtau,0,vr);
		T vinv = v > 0 ? 1/v : T(0);
		vec3<T> N_V(vr*vinv,0,vtau*vinv);

		vec3<T> N1(-N_V[2],0,N_V[0]);
		vec3<T> N2(0,1,0);

		T sin_theta_in = std::sqrt(1-cos_theta_in*cos_theta_in);
		vec3<T> v1_vec = N_V*(v1*cos_theta_in) +
			N1*(v1*sin_theta_in*std::cos(phi_in)) + 
			N2*(v1*sin_theta_in*std::sin(phi_in)); 
		
		vec3<T> Vu = V_in-v1_vec;
		auto [e1,e2,e_vu] = Basis3(Vu);
		
		T sin_theta_out = std::sqrt(1-cos_theta_out*cos_theta_out);
		vec3<T> V_cm = mu_xi*V_in+mu_i*v1_vec;
		vec3<T> Vu_out_vec = Vu_out*(
			e_vu*cos_theta_out + 
			e1*sin_theta_out*std::cos(phi_out),
			e2*sin_theta_out*std::sin(phi_out),
		);
		vec3<T> Vfinal = V_cm + mu_i*Vu_out_vec;
		
		T Vout2 = Vfinal.squaredNorm();
		T L = r*std::sqrt(Vfinal[0]*Vfinal[0]+Vfinal[1]*Vfinal[1]);
		return {Vout2-phi_r,L};
	}



	template <typename T>
	struct ScatterRVExpInfo_t{
		constexpr static T m_factor = 1.2;
		constexpr static T m_factor_inv = 1/m_factor;
		constexpr static size_t size = 48;
		constexpr static T min_value = 1e-37;
		std::vector<T> Majorants;
		T full_prob;
		T u_tres;

		ScatterRVExpInfo_t(){}
		ScatterRVExpInfo_t( std::vector<T> Majorants,T full_prob,T u_tres):
			Majorants(std::move(Majorants)), 
			full_prob(full_prob),
			u_tres(u_tres)
		{
			if(Majorants.size() != size){
				throw std::runtime_error("ScatterExpInfo_t: invalid Majorants size");
			}
		}
	
		inline static grob::Rect<T> box(size_t i){
			return {
				i !=0 ? std::pow(m_factor_inv, size-i) : 0, 
				std::pow(m_factor_inv, size-i-1)};
		}
		

		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("Majorants","full_prob","u_tres"), 
			PROPERTIES(Majorants, full_prob,u_tres)
		);
		DESERIALIZATOR_FUNCTION(
			ScatterRVExpInfo_t, 
			PROPERTY_NAMES("Majorants","full_prob","u_tres"), 
			PROPERTY_TYPES(Majorants, full_prob,u_tres)
		);
	};
	
	template <typename T>
	struct ScatterRVExp : ScatterRVExpInfo_t<T>{
		using ScatterRVExpInfo_t<T>::ScatterRVExpInfo_t;

		std::vector<T> interval_samples;

		ScatterRVExp(){}
		ScatterRVExp(std::vector<T> majorants std::vector<T> interval_samples):
			interval_samples(interval_samples){}
		//"inherited" box function
		inline static grob::Rect<T> box(size_t i){
			return ScatterExpInfo_t<T>::box(i);
		}
		//construct from ScatterExpInfo_t
		template <typename iter_t>
		ScatterRVExp(iter_t majorant_begin,iter_t majorant_end):
		{
			if( majorant_end - majorant_begin != size){
				throw std::runtime_error("majorant_end - majorant_begin != size");
			}
			interval_samples[0] = (*majorant_begin)*box(0).volume();
			++majorant_begin;
			for(size_t i=1;i<size;++i,++majorant_begin){
				interval_samples[i] = 
					interval_samples[i] + 
					SEI.Majorants[i]*box(i).volume();
			}
			T max_prob = interval_samples.back();
			for(size_t i=0;i<size;++i){
				interval_samples[i] /= max_prob;
			}
		}
		size_t gen_index(T xi)const{
			return evdm::IndexGenBigHelper::gen(
				interval_samples.begin(), interval_samples.end(), xi
			);
		}

		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("interval_samples"), 
			PROPERTIES(interval_samples)
		);
		DESERIALIZATOR_FUNCTION(
			ScatterRVExp, 
			PROPERTY_NAMES("interval_samples"), 
			PROPERTY_TYPES(interval_samples)
		);
	};


	template <typename T>
	using QT_RV_t = Quadtree<T,ScatterRVExpInfo_t<T>,ScatterRVExp<T>>;

	template <typename T>
	struct ElementsRVQuadTree{
		std::vector<QT_RV_t> QTs_per_element;
	};
	
	template <typename T>
	struct StateEL{
		size_t ptype;
		T e,l;
	};

	template <typename T,typename FormFactorConstructor_impl>
	struct ScatterInfoEL_in_t {
		typedef ScatterInfoElements<T> ScInfoEls_t;
		using QT_t =  Quadtree<T,ScatterInfoElements<T>>;
		using Node_t = typename QT_t::Node;

		typedef std::vector<QT_t>  QTrees_t;
		typedef std::vector<ElementsRVQuadTree<T>>  RV_Info_t;

		size_t ptypes;
		FormFactorConstructor_impl FFC_impl;
		std::vector<int> ValueIndexes;
			//-1 if no value
		std::vector<std::vector<ElementInfo_t>> FormFactors;
			//index of vectors - is ptype_in*ptypes + ptype_out
		QTrees_t QuadTrees;
		RV_Info_t RV_Info;
		ScatterInfoEL_in_t(
			size_t ptypes,
			FormFactorConstructor_impl FFC_impl,				
			std::vector<int> ValueIndexes = {},
			QTrees_t QuadTrees = {},
			RV_Info_t RV_Info = {} 
		):
			ptypes(ptypes),
			FFC_impl(std::move(FFC_impl)),
			ValueIndexes(std::move(ValueIndexes)),
			QuadTrees(std::move(QuadTrees)),
			RV_Info(std::move(RV_Info))
		{
			if(ValueIndexes.size() < ptypes*ptypes){
				ValueIndexes.resize(ptypes*ptypes,-1);
			}
			FormFactors.resize(QuadTrees.size());
			size_t m_size = QuadTrees.size();
			for(size_t i=0;i<ptypes;++i){
				for(size_t j=0;j<ptypes;++j){
					int index = ValueIndexes[i*ptypes + j];
					if(index >= m_size){
						throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t:"
							" invalid QuadTrees entry");
					}
					if(index>=0){
						FormFactors[index] = make_element_info(
							FFC_impl,i,j
						);
						//sort: first with bigger mN
						FormFactors[index].sort(
							[](const ElementInfo_t & a,const ElementInfo_t & b){
								return a.mN > b.mN;
							}
						);
					} 
				}
			}
		}

		template <typename Gen_t>
		std::pair<StateEL<T>,T> gen_next(
			StateEL<T> current,T max_time,
			T * p_in_p_out_ProbMatrix,const Body<T> & B,Gen_t && G
		)const{
			//calculate scatter probability
			//generate ptype_out
			std::array<T,4> ptype_out_probs;
			std::array<size_t ,4> ptype_out_indexes;
			if(ptypes > ptype_out_probs.size()){
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::gen_next:"
					" too many ptypes for fixed size arrays");
			}
			size_t array_index = 0;
			for(size_t ptype_out=0;ptype_out<ptypes;++ptype_out){
				if(hasValue(current.ptype,ptype_out)){
					const QT_t & qt = getQT(current.ptype,ptype_out);
					typename QT_t::Node const & Node = qt.find_node(current.e,current.l);
					auto [alpha,betha ] = Node.coeffs(current.e,current.l);
					auto [p00,p01,p10,p11] = Node.values_view(  
						[](const ScatterInfoElements<T> & siot){
							return siot.full_prob;
						}
					);
					T p_elem = interpolate_log(p00,p01,p10,p11,alpha,betha);
					ptype_out_probs[array_index] = 
						p_elem * p_in_p_out_ProbMatrix[current.ptype*ptypes + ptype_out];
					ptype_out_indexes[array_index] = ptype_out;
					++array_index;
				}
			}
			T scat_prob = std::accumulate(
				ptype_out_probs.begin(),
				ptype_out_probs.begin()+array_index,
				(T)0
			);

			if(scat_prob <= 0){
				return {current, max_time};
			}
			//generate scatter time
			auto xi = E0I1_G(G);
			T dt_scat = -std::log(xi) / scat_prob;
			if (dt_scat >= max_time) {
				return { current, max_time };
			}
			//scatter
			//generate ptype_out
			for(size_t i=1;i<array_index;++i){
				ptype_out_probs[i] += ptype_out_probs[i-1];
			}
			auto xi_pt = E0I1_G(G)*ptype_out_probs[array_index-1];
			size_t sel_index = evdm::IndexGenBigHelper::gen(
				ptype_out_probs.begin(),
				ptype_out_probs.begin()+array_index,
				xi_pt
			);
			size_t ptype_out = ptype_out_indexes[sel_index];

			//generate e,l after scatter
			//generate element
			const QT_t & qt = getQT(current.ptype,ptype_out);
			auto const & Node = qt.find_node(current.e,current.l);

			//to generate element we need to:
			//1) make interpolated range from corners
			T xi_elem = G();
			int i00 = Node(0,0).gen_index(xi_elem);
			int i01 = Node(0,1).gen_index(xi_elem);
			int i10 = Node(1,0).gen_index(xi_elem);
			int i11 = Node(1,1).gen_index(xi_elem);

			int sel_elem_index = std::min(std:min(i00,i01),std::min(i10,i11));
			if(sel_elem_index < 0){
#ifndef NDEBUG
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::gen_next:"
					" can't generate element index after scatter");
#else
				return {current, max_time};
#endif
			}
			//3) get theta distribution for selected element
			// we will use min from corners
			auto xi_th = G();
			grob::Point<T,T,T,T> thetas = Node.values_view(
				[sel_elem_index,xi_th](const ScatterInfoElements<T> & siot)->T{
					const auto & theta_sample = siot.ThetaTrajs[sel_elem_index].SampleTheta;
					return evdm::IndexGenBigHelper::gen(
						theta_sample.begin(),
						theta_sample.end(),
						xi_th
					);
				}
			);
			using namespace grib::literals;
			T theta = tuple_min(thetas);
			
			grob::Point<size_t,size_t,size_t,size_t> il0_hints = Node.indices_view(
				[sel_elem_index](const ScatterInfoElements<T> & siot)->size_t{
					return siot.i0_lmax();
				}
			);
			//find minimum of tuple with boost
			size_t il0_hint =  tuple_min(il0_hints); //minimum of tuple
			

			grob::Point<size_t,size_t,size_t,size_t> il1_hints = Node.indices_view(
				[sel_elem_index](const ScatterInfoElements<T> & siot)->size_t{
					return siot.i1_lmax();
				}
			);
			size_t il1_hint =  tuple_max(il1_hints);
			auto [L2,r2_max ] = B.maxL2(-current.e,il0_hint,il1_hint);

			
			auto [rmin,rmax] = B.find_rmin_rmax(-current.e,L2,r2_max);

			T theta_max = B.get_theta_max(rmin,rmax);
			T theta_tmp = theta_max*theta;
			auto cth2 = rmin * std::cos(theta / 2);
			auto sth2 = rmax * std::sin(theta / 2);
			T r = std::sqrt( cth2*cth2+ sth2*sth2);
			T phi_r = B.Phi(r);
			T v2 = downbound(phi_r + current.e, 0);
			auto v = std::sqrt(v2);
			T vt = upbound(std::sqrt(L2) / r, v);
			T vr = ssqrt(v2 - vt * vt);
			//vec3<T> V_in = V_in_nd * VescMin;
			T vesc = B.Vesc;
			auto const & ff =  FormFactors[getIndex( current.ptype ,ptype_out)][sel_elem_index];
			T mu_chi = ff.mX*(1/(ff.mN+ff.mX));
			T mu_i = ff.mN*(1/(ff.mN+ff.mX));
			T Mu_chi_i = ff.mN*mu_chi;
			T Delta_M_chi = ff.dmX;
			T v2_delta = - 2*ff.dmX/Mu_chi_i;
			T v1_temp_sq = B.Temp(r)/ff.mN/(vesc*vesc);
			T v1_tres = v2_delta > v2 ? std::sqrt(v2_delta) - v : T(0);
			T u_tresh = v2_delta > 0 ? (v1_tres*v1_tres)/(2*v1_temp_sq) : T(0);

			const QT_RV_t<T> & m_RV_Qt = RV_Info[
				getIndex(current.ptype,ptype_out)
			][sel_elem_index];

			typename QT_RV_t<T>::Node const & NodeRV =  
				m_RV_Qt.find_node(r,v/std::sqrt(B.Phi(r)));

			
			for(size_t attempt=0;attempt<100;++attempt){
				auto xi_rv = G();
				ScatterRVExp<T> const & m_scatter_rv = NodeRV.uvalue();

				size_t rv_index =  m_scatter_rv.gen_index(G());
				T majorant = m_scatter_rv.Majorants[rv_index];
				auto [xi0,xi1] = m_scatter_rv.box(rv_index);
				T xi = m_scatter_rv.min_value + ((xi0) + G()*(xi1 - xi0));
				T s = std::log(xi);
				T u_in = (s + u_tresh);
				T v1 = std::sqrt(2*u_in*v1_temp_sq);
				
				scatter_pre_result<T> m_s = 
					gen_scatter(G(),G(),v,v1,u_in,v2_delta);
				T Vu_in2 = v*v+v1*v1+2*v*v1*m_s.cos_theta_in;
				T Vu_in = std::sqrt(Vu_in2);
				T Vu_out = std::sqrt(downbound(Vu_in2 + v2_delta,0));
				
				T m_factor = get_ff_out(
					Vu_in,Vu_out,Mu_chi_i,Delta_M_chi,vesc,m_s.cos_theta_out,ff.ff);
				
				T from_factor = m_factor*m_s.gen_pre_form_factor;

				if(from_factor <= majorant){
					
				}
				#ifndef NDEBUG
				if(from_factor > majorant){
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::gen_next:"
					" can't generate element index after scatter");
				}
				#endif
				if(from_factor > majorant*G()){
					auto [newE,newL] = get_el_out(
						mu_chi,mu_i,r,phi_r,vr,vt,v,v1,Vu_in,Vu_out,
						m_s.cos_theta_in,2*std::numbers::pi_v<T>*G(),
						m_s.cos_theta_out,2*std::numbers::pi_v<T>*G()
					);
					if(newE < 0){
						T Lmax = std::sqrt(B.maxL2(-newE));
						T l = Lmax > 0 ? downbound(newL/Lmax,1) : T(0);
						return { {ptype_out,newE,l}, max_time-dt_scat};
					} else{
						return { {ptype_out,newE,0}, 0};
					}
					return  
				}
			}
			
		}
		
		T ScatterProb(T e,T l,size_t ptype_in,T * ptype_out_mlt_factors) const{
			T prob = 0;
			for(size_t ptype_out=0;ptype_out<ptypes;++ptype_out){
				if(hasValue(ptype_in,ptype_out)){
					const QT_t & qt = getQT(ptype_in,ptype_out);
					typename QT_t::Node const & Node = qt.find_node(e,l);
					auto [alpha,betha ] = Node.coeffs(e,l);
					auto [p00,p01,p10,p11] = Node.values_view(  
						[](const ScatterInfoElements & siot){
							return siot.full_prob;
						}
					);
					T p_elem = interpolate_log(p00,p01,p10,p11,alpha,betha);
					prob += p_elem * ptype_out_mlt_factors[ptype_out];
				}
			}
			return prob;
		}

		SERIALIZATOR_FUNCTION(
		PROPERTY_NAMES("ptypes","FFC_impl","ValueIndexes","QuadTrees","RV_Info"), 
		PROPERTIES(ptypes,FFC_impl,ValueIndexes,QuadTreesRV_Info) 
		);
		DESERIALIZATOR_FUNCTION(
			ScatterInfoEL_in_t, 
			PROPERTY_NAMES("ptypes","FFC_impl","ValueIndexes","QuadTrees","RV_Info"), 
			PROPERTY_TYPES(ptypes,FFC_impl,ValueIndexes,QuadTrees,RV_Info)
		);

		void addQT(
			size_t ptype_in,
			size_t ptype_out,
			QT_t qt,
			ElementsRVQuadTree<T> rv_info
		){
			ValueIndexes[ptype_in*ptypes + ptype_out] = QuadTrees.size();
			QuadTrees.push_back(std::move(qt));

			FormFactors.push_back(
				make_element_info(
					FFC_impl,ptype_in,ptype_out
				)
			);
			RV_Info.push_back(
				std::move(rv_info)
			);
		}

		bool hasValue(size_t ptype_in,size_t ptype_out) const{
			return HasValueMask[ptype_in*ptypes + ptype_out]>0;
		}
		int getIndex(size_t ptype_in,size_t ptype_out) const{
			return ValueIndexes[ptype_in*ptypes + ptype_out];
		}
		const ElementsRVQuadTree<T> &getRVInfo(size_t ptype_in,size_t ptype_out) const{
			if(!hasValue(ptype_in,ptype_out)){
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::getRVInfo:"
					" no quadtree for given ptype_in,ptype_out");
			}
			return RV_Info[getIndex(ptype_in,ptype_out)];
		}
		ElementsRVQuadTree<T> &getRVInfo(size_t ptype_in,size_t ptype_out){
			return const_cast<ElementsRVQuadTree<T> &>(
				static_cast<const ScatterInfoEL_in_t*>(this)->getRVInfo(ptype_in,ptype_out)
			);
		}
		const QT_t & getQT(size_t ptype_in,size_t ptype_out) const{ 
			if(!hasValue(ptype_in,ptype_out)){
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::getRVInfo:"
					" no quadtree for given ptype_in,ptype_out");
			}
			return QuadTrees[getIndex(ptype_in,ptype_out)];
		}
		QT_t & getQT(size_t ptype_in,size_t ptype_out){ 
			return const_cast<Quadtree<T,ScatterInfoElements>&>(
				static_cast<const ScatterInfoEL_in_t*>(this)->getQT(ptype_in,ptype_out)
			);
		}
	}; 


	template <typename T,typename FFC_impl_t>
	struct MCEvolutor {
		//content
		evdm::Body<T> B;
		//typedef decltype (B.get_le(0)) EL_Func_t;
		//EL_Func_t LE;
		ScatterInfoEL_in_t<T,FFC_impl_t> ScatterInfoEL;

		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("body","ScatterInfoEL"), 
			PROPERTIES(B,ScatterInfoEL)
		);
		DESERIALIZATOR_FUNCTION(
			MCEvolutor, 
			PROPERTY_NAMES("body","ScatterInfoEL"), 
			PROPERTY_TYPES(B,ScatterInfoEL)
		);
		
		MCEvolutor(
			evdm::Body<T> B,
			//EL_Func_t LE,
			ScatterInfoEL_in_t ScatterInfoEL
		):
			B(std::move(B)),
			//LE(std::move(LE)),
			ScatterInfoEL(std::move(ScatterInfoEL))
		{}

		MCEvolutor(
			evdm::Body<T> B,
			size_t ptypes,
			FormFactorConstructor_impl FFC_impl
		):B(std::move(B)),ScatterInfoEL(ptypes,FFC_impl)
		{}
		void AddProcess(
			size_t ptype_in,size_t ptype_out,
			size_t initial_rv_level, size_t max_rv_level,size_t max_rv_size,
			size_t initial_el_level, size_t max_el_level,size_t max_el_size,
			){
			// TODO 1) for each element fill RV info quadtree
			// critrium to split: 
			//		bad interpolation of probs
			//		big difference in treshold
			//		big difference in majorants
			// 2) using RV info fill elinfo quadtree
			//	critrium to split: 
			//		bad interpolation of probs for any element
			//ScatterInfoEL.addQT(ptype_in,ptype_out,);
			throw std::runtime_error("unimplemented");
		}
		void AddSameProcess(size_t ptype_in_dsc,size_t ptype_out_dsc,
		size_t ptype_in_src,size_t ptype_out_src){
			ScatterInfoEL.addQT(ptype_in_dsc,ptype_out_dsc,
					ScatterInfoEL.getQT(ptype_in_src,ptype_out_src),
					ScatterInfoEL.getRVInfo(ptype_in_src,ptype_out_src)
			)
		}

		template <typename Gen_t>
		std::vector<StateEL<T>> evolute(
			std::vector<StateEL<T>> init,
			Gen_t && _G,
			std::array<T,16> prob_matrix,
			T evolve_time,
			size_t max_scatter,
		){
			int N = init.size();
			const auto _seed = _G.state;
			#pragma omp parallel for
			for(int i=0;i<N;++i){
				auto G = _G;
				G.set_seed(_seed ^ i + 1);
				max_scatter += std::floor( G()*ScatterInfoEL.ptypes);
				StateEL<T> tmp_state = init[i];
				T tmp_time = evolve_time;
				for(size_t j=0;j<max_scatter;++j){
					std::tie(tmp_state,tmp_time) = 
						ScatterInfoEL.gen_next(in,tmp_time,prob_matrix.data(),B,G);
					if(tmp_time == 0){
						break;
					}
				}
				init[i] = tmp_time;
			}
		}

		private:
		ScatterRVExpInfo_t<T> genRVinfo(
			size_t ptype_in, size_t ptype_out,T r,T v,T n_r, ElementInfo_t const & El,size_t Nmk)
		{
			T mu_chi = ff.mX*(1/(ff.mN+ff.mX));
			T mu_i = ff.mN*(1/(ff.mN+ff.mX));
			T Mu_chi_i = ff.mN*mu_chi;
			T Delta_M_chi = ff.dmX;
			T v2_delta = - 2*ff.dmX/Mu_chi_i;
			T vesc = B.Vesc;

			T v2 = v*v;
			T v1_temp_sq = B.Temp(r)/ff.mN/(vesc*vesc);
			T v1_tres = v2_delta > v2 ? std::sqrt(v2_delta) - v : T(0);
			T u_tresh = v2_delta > 0 ? (v1_tres*v1_tres)/(2*v1_temp_sq) : T(0);

			using RVPi = ScatterRVExpInfo_t<T>;
			T sum = 0;
			std::vector<T> max_ffs[RVPi::size];
			for (int i=0;i<RVPi::size){
				auto [xiu_0,xiu_1] = RVPi::box(i);
				T max_ff = 0;
				for(int _k : std::iota(Nmk)){
					T xi = RVPi::min_value+(xiu_0 + (xiu_1-xiu_0)*G());
					T s = std::log(xi);
					T u_in = (s + u_tresh);
					T v1 = std::sqrt(2*u_in*v1_temp_sq);
					scatter_pre_result<T> m_s = 
						gen_scatter(G(),G(),v,v1,u_in,v2_delta);
				
					T Vu_in2 = v*v+v1*v1+2*v*v1*m_s.cos_theta_in;
					T Vu_in = std::sqrt(Vu_in2);
					T Vu_out = std::sqrt(downbound(Vu_in2 + v2_delta,0));
					
					T m_factor = get_ff_out(
						Vu_in,Vu_out,Mu_chi_i,Delta_M_chi,vesc,m_s.cos_theta_out,ff.ff);
					T from_factor = m_factor*m_s.gen_pre_form_factor;
					max_ff = std::max(max_ff,from_factor);
					sum += from_factor*n_r;
				}
				max_ffs[i] = max_ff;
			}
			return ScatterRVExpInfo_t<T>(std::move(max_ffs),sum/Nmk,u_tresh);
		}

		template <typename Func_t>
		ScatterInfoTheta<T> genTrajProbs(
			T e, T l, T rmin,T rmax,T theta_max,T Tin,T Tout,QT_RV_t<T> const & rv_probs,
			size_t max_steps,size_t max_lvl,T interpol_accept,
		){
			struct PDF_t{
				T th;
				T p;
				//T max_error=0;
				//T total=0;
			};
			std::list<PDF_t,PoolAllocator<T> > 
				m_points(PoolAllocator<T>(max_steps));

			T _1 = 1;

			T dP_dth = [&](auto theta_undim){
				auto [tau, dtau_dth] = tau_theta(theta);
				T theta_full = theta_undim*theta_max;
				T ct = std::cos(theta_full);
				T st = std::sin(theta_full);
				
				// Check that rmax > 1 in outer trajectories?
				T r = sqrt(powint(rmin*ct)+powint(rmax*st)); 
				T dtau_dth = std::sqrt(B.S(r,rmin,rmax));//S or 1/S?
				T v = ssqrt(e - B.Phi(r));
				T vmax = ssqrt(B.Phi(r));
				if(vmax == 0){
					return 0;
				}
				auto const & N = 
					rv_probs.find_node(r,v/vmax);
				auto [P00,P01,P10,P11] = N.values_view([](auto const & x){return x.full_prob;});
				auto [alph,bth]= N.coeffs(r,v/vmax);
				
				return dtau_dth*interpolate_log(
					interpolate_log(P00,P10,alph),
					interpolate_log(P01,P11,alph),
					bth
				);
			}

			m_points.push_back({0,dP_dth(0)});
			m_points.push_back({_1,dP_dth(_1)});

			T FullProbLow = m_points.back().p * 1;

			T PFD_max = m_points.front().p;

			if(max_steps < 5){
				throw std::runtime_error("genTrajProbs: max_steps < 4");
			}
			auto intr_log(PDF_t P0,PDF_t P1){
				return interpolate_log(P0.p,P1.p,T(0.5));
			};
			auto intr_lin(PDF_t P0,PDF_t P1){
				return (P0.p + P1.p)/2;
			};
			
			T m_zero = 1e-9*m_points.front().second;

			max_steps -= 5;
			
			for(auto it = m_points.begin();it != std::prev(m_points.end());++it){
				T tmp = it->p;
				T nxt = std::next(it)->p;
				if( nxt < tmp/2){
					m_points.
				}
			}

			//first step: minimize TV:
			while (max_steps >0)
			{
				for(auto it = m_points.begin();it != std::prev(m_points.end());++it){
					//
					m_points.insert(it,);
				}
			}
			
			

			
		}
	};

};
