#pragma once

#include "dynamics_evolutor.hpp"

namespace evdm{

    template <typename T>
	struct scatter_pre_result{
		T cos_theta_in;
		T cos_theta_out;
		T gen_pre_form_factor;
	};
	/// @brief DeltaV2 is positive  in exotermic process
	template <typename T>
	inline scatter_pre_result<T> gen_scatter(
		T xi_theta_in,same_t<T> xi_theta_out,
		same_t<T> V_in,same_t<T> v1_in, same_t<T> u_in,
		same_t<T> DeltaV2){
		T A = (v1_in*v1_in + V_in*V_in + DeltaV2);
		T B = 2*V_in*v1_in;
		
		if(A > B){
			T t = B/A;
			//evaluate (1 - (1+x)^{2/3})/x
			auto Wfunc_2_3 = [](T one_plus_x){
				T q = std::cbrt(one_plus_x);
				return -(1+q)/(1+q+q*q);
			};
			//evaluate ((1+x)^{3/2} - 1)/x
			auto Wfunc_3_2_p = [](T x){
				T x1 = 1+x;
				T q = x1 * std::sqrt(x1);
				return (1 + x1 + x1*x1)/(1 + q);
			};
			//evaluate (1 - (1-x)^{3/2})/x
			auto Wfunc_3_2_m = [](T x){
				T x1 = 1-x;
				T q = x1 * std::sqrt(x1);
				return (1 + x1 + x1*x1)/(1 + q);
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

			T t_plus_1_pow_3_2 = (1+t)*std::sqrt(1+t);
			T cos_theta_in = t - std::pow( t_plus_1_pow_3_2*(1-xi_theta_in),T(2)/T(3) );
			T cos_theta_out = -1 + 2*xi_theta_out;
			T m_factor = t_plus_1_pow_3_2*(T(2)/T(3));
			return {
				cos_theta_in,cos_theta_out, 
				m_factor*std::sqrt(B * u_in * std::numbers::inv_pi_v<T>)};
		}
		else{
			return {
				-1 + 2*xi_theta_in,-1 + 2*xi_theta_out,
				T(0)};
		}
	}

	///@brief v2delta positive in exotermic process
	template <typename T>
	inline T get_ff_out(
		T Vu_in,T Vu_out,T mu,T v2delta,T VescMin,T cos_theta_out,
		FormFactor_t const & FF){
		T m_factor = mu*mu*(VescMin*VescMin);
		
		T delta_Vu2 = Vu_in*Vu_in+ Vu_out*Vu_out - 2*Vu_in*Vu_out*cos_theta_out;
		T t = delta_Vu2 >0 ? v2delta/(delta_Vu2) : T(0);
		
		T Vu_in_t = Vu_in *( T(0.5)-t);
		T Vu_out_t = Vu_out *( T(0.5)+t);

		T V2_perp = Vu_in_t*Vu_in_t + Vu_out_t*Vu_out_t - 2*Vu_in_t*Vu_out_t*cos_theta_out;
		
		return FF.ScatterFactor(
			m_factor*delta_Vu2,
			(VescMin*VescMin)*V2_perp,
			0
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
		T cth_V = v > 0 ? vr * (1/v) : 1;
		T sth_V = v > 0 ? vtau * (1/v) : 0;

		vec3<T> N_V(sth_V,0, cth_V);

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
			e1*sin_theta_out*std::cos(phi_out) +
			e2*sin_theta_out*std::sin(phi_out)
		);
		vec3<T> Vfinal = V_cm + mu_i*Vu_out_vec;
		
		T Vout2 = Vfinal.squaredNorm();
		T L = r*std::sqrt(Vfinal[0]*Vfinal[0]+Vfinal[1]*Vfinal[1]);
		return {Vout2-phi_r,L};
	}


    /// @brief struct with majorants for particular element, r, v 
	template <typename T>
	struct ScatterRVExpInfo_t{
		constexpr static T m_factor = 3;
		constexpr static T m_factor_inv = 1/m_factor;
		constexpr static size_t size = 12;
		constexpr static T min_value = 1e-20;
		std::array<T,size> Majorants;
		T full_prob;
		T u_tres;

		ScatterRVExpInfo_t():full_prob(0),u_tres(0){}
		ScatterRVExpInfo_t( std::array<T, size> Majorants,T full_prob,T u_tres):
			Majorants(Majorants), 
			full_prob(full_prob),
			u_tres(u_tres)
		{
			if(this->Majorants.size() != size){
				throw std::runtime_error("ScatterExpInfo_t: invalid Majorants size");
			}
		}
	
		inline static grob::Rect<T> box(size_t i){
			return {
				i !=0 ? powint(m_factor_inv, size-i) : (T)0, 
				powint(m_factor_inv, size-i-1)};
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
	
	/// @brief struct to sample intervals from majorants
	template <typename T>
	struct ScatterRVExp{
		using Base = ScatterRVExpInfo_t<T>;
		constexpr static T m_factor = Base::m_factor;
		constexpr static T m_factor_inv = Base::m_factor_inv;
		constexpr static size_t size = Base::size;
		constexpr static T min_value = Base::min_value;

		std::array<T,size> Majorants;
		std::array<T,size> interval_samples;

		ScatterRVExp(){}
		ScatterRVExp(std::array<T, size> Majorants,std::array<T, size> interval_samples):
			Majorants(Majorants),interval_samples(interval_samples){}
		//"inherited" box function
		inline static grob::Rect<T> box(size_t i){
			return Base::box(i);
		}
		//construct from ScatterExpInfo_t
		template <typename iter_t>
		ScatterRVExp(iter_t majorant_begin,iter_t majorant_end)
		{
			//interval_samples.resize(std::distance(majorant_begin,majorant_end));
			if (std::distance(majorant_begin, majorant_end) != size) {
				throw std::runtime_error("majorant_end - majorant_begin != size");
			}
			size_t i = 0;
			for(auto it =majorant_begin;it!=majorant_end;++it,++i){
				Majorants[i] = *it;
			}
			
			interval_samples[0] = (*majorant_begin)*box(0).volume();
			++majorant_begin;
			for(size_t i=1;i<size;++i,++majorant_begin){
				interval_samples[i] = 
					interval_samples[i-1] + 
					*majorant_begin*box(i).volume();
			}
			T max_prob = interval_samples.back();
			for(size_t i=0;i<size;++i){
				interval_samples[i] /= max_prob;
			}
		}
		size_t gen_index(T xi)const{
			return IndexGenBigHelper::gen(
				interval_samples.begin(), interval_samples.end(), xi
			);
		}

		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("Majorants","interval_samples"), 
			PROPERTIES(Majorants,interval_samples)
		);
		DESERIALIZATOR_FUNCTION(
			ScatterRVExp, 
			PROPERTY_NAMES("Majorants","interval_samples"), 
			PROPERTY_TYPES(Majorants,interval_samples)
		);
	};


	/// @brief quad tree for each element in space (r,v/v_max(r)) 
	template <typename T>
	using QT_RV_t = Quadtree<T,ScatterRVExpInfo_t<T>,ScatterRVExp<T>>;

	
	/// @brief vector of qt for each element
	template <typename T>
	using ElementsRVQuadTree = 
		std::vector<QT_RV_t<T>> ; //QTs_per_element
	
	

    template <typename T>
    struct dP_dth_functor{
        evdm::Body<T> const & B;
        QT_RV_t<T> const & rv_probs;
		T e, rmin, rmax;
        T theta_max;
        T FullPeriod;
		
        dP_dth_functor(
            evdm::Body<T> const & B,QT_RV_t<T> const & rv_probs,
            T e, T rmin, T rmax,
            T theta_max,T FullPeriod
        ):B(B),rv_probs(rv_probs),
        e(e),rmin(rmin),rmax(rmax),
        theta_max(theta_max),FullPeriod(FullPeriod){}
        
        inline T operator ()(T theta_undim)const{
            T theta_full = theta_undim*theta_max;
            T ct = std::cos(theta_full/2);
            T st = std::sin(theta_full/2);
            
            // Check that rmax > 1 in outer trajectories?
            T r = sqrt(powint(rmin*ct,2)+powint(rmax*st,2)); 
            if(r>=1){
                return 0;
            }
            T dtau_dth = 1/(2*std::sqrt(B.S(r,rmin,rmax)) );//S or 1/S?
            T v = ssqrt(e + B.Phi(r));
            T vmax = ssqrt(B.Phi(r));
            if(vmax == 0){
                return 0;
            }
            auto const & N = 
                rv_probs.find_node(r,v/vmax);
            auto [P00,P01,P10,P11] = N.values_view([](auto const & x){return x.full_prob;});
            auto [alph,bth]= N.coeffs(r,v/vmax);
            
            return dtau_dth* theta_max *interpolate_log(
                P00,P01,P10,P11,
                alph,bth
            )/FullPeriod;
        }

    };
};