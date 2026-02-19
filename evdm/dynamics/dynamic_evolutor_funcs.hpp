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

	struct GenTrajProbs_return_full {};
	template <typename T,typename U = GenTrajProbs_return_full>
	T GenTrajProbs(
		const Body<T>& B, T e, T l, T rmin, T rmax, T theta_max,
		T Tin, T Tout, QT_RV_t<T> const& rv_probs,
		std::span<grob::Point<T, T>> m_points, 
		std::span<grob::Point<T, T>> m_buffer,
		size_t max_lvl, T MaxProbDiff,
		U xi_or_return_full = GenTrajProbs_return_full{}
	) {

		using namespace grob::literals;
		typedef grob::Point<T, T> PDF_t;


		if (m_points.size() < 17) {
			throw std::runtime_error("GenTrajProbs: max_steps < 16");
		}
		size_t max_steps = m_points.size() - 1;
		if (m_buffer.size() < max_steps + 1) {
			throw std::runtime_error("GenTrajProbs: m_buffer.size() < m_points.size()");
		}

		size_t tmp_size = 0;

		T _1 = 1;

		//B.get_internal_period(rmin,rmax,100);genTrajProbs
		dP_dth_functor<T> dP_dth(B, rv_probs, e, rmin, rmax, theta_max, Tin + Tout);

		m_points[0] = { 0,dP_dth(0) };
		m_points[1] = { _1 / 2,dP_dth(_1 / 2) };
		m_points[2] = { _1,dP_dth(_1) };
		tmp_size += 3;

		T PFD_max = m_points.front()[1_c];


		max_steps -= 2;
		auto intr_log = [](PDF_t P0, PDF_t P1) {
			return interpolate_log(P0[1_c], P1[1_c], T(0.5));
			};
		auto intr_lin = [](PDF_t P0, PDF_t P1) {
			return (P0[1_c] + P1[1_c]) / 2;
			};

		auto new_point = [&dP_dth](PDF_t P0, PDF_t P1)->PDF_t {
			T x = (P0[0_c] + P1[0_c]) / 2;
			return PDF_t(x, dP_dth(x));
			};
		//T m_zero = 1e-9*m_points.front().second;


		{// get rid of zero elements until there is positive element

			PDF_t& p0 = m_points.front();
			PDF_t& p1_2 = *std::next(m_points.begin());
			PDF_t& p1 = m_points.back();
			T remain_prob = std::max(p1_2[1_c], p1[1_c]) * p1_2[0_c];
			for (size_t i = 0; i < max_lvl; ++i) {
				if (p0[1_c] * p1_2[0_c] * 1e-10 > remain_prob) {
					PDF_t np1 = new_point(p0, p1_2);
					p1 = p1_2;
					p1_2 = np1;
					remain_prob += p1_2[1_c] * p1_2[0_c];
				}
				else {
					break;
				}
			}
		}

		auto PointComp = [](auto const& x, auto const& y) {
			return x[0_c] < y[0_c];
			};
		auto push_points = [PointComp](
			std::span<PDF_t>& m_points, size_t& tmp_size,
			size_t new_points_num, std::span<PDF_t>& buff)
			{
				if (new_points_num) {
					auto beg = m_points.begin();
					std::merge(beg, beg + tmp_size, beg + tmp_size,
						beg + tmp_size + new_points_num, buff.begin(),
						PointComp);
					std::swap(m_points, buff);
					tmp_size += new_points_num;
				}
			};

		for (size_t i = 0; i < 3; ++i) {
			size_t num = 0;
			for (size_t i = 0; i < tmp_size - 1; ++i) {
				T x0 = m_points[i][0_c];
				T x1 = m_points[i + 1][0_c];
				T x = (x0 + x1) / 2;
				m_points[tmp_size + i] = { x,dP_dth(x) };
				++num;
				--max_steps;
			}
			push_points(m_points, tmp_size, num, m_buffer);
		}//refine grid to 17 points

		auto return_value = [&](size_t size) -> T {
			
			// calclulate full probability
			//log avarage
			auto AvP = [](auto p0, auto p1) {
				auto Pav = (p0 + p1) / 2;
				if (Pav > 0) {
					return Pav * WLogLinear((p0 - p1) / Pav);
				}
				else {
					return Pav;
				}
			};
			//integrate probability	
			std::span<PDF_t> m_p{ m_points.begin(),size };
			T p_sum = 0;
			{
				auto it = m_p.begin();
				grob::Point<T, T>& _it = *it;
				T p0 = _it[1_c];
				T x0 = _it[0_c];
				_it[1_c] = 0;
				++it;
				for (; it != m_p.end(); ++it) {
					auto [x1, p1] = *it;
					p_sum += AvP(p0, p1) * (x1 - x0);
					(*it)[1_c] = p_sum;
					p0 = p1;
					x0 = x1;
				}
			}

			if constexpr (std::is_same_v< U, GenTrajProbs_return_full> ) {
				return p_sum;
			}
			else {
				T target = p_sum * xi_or_return_full;
				auto m_cdf = m_p | std::views::transform([](auto x) {return x[1_c];});
				auto m_theta = m_p | std::views::transform([](auto x) {return x[0_c]; });
				size_t index = IndexGenBigHelper::gen(m_cdf.begin(), m_cdf.end(), target);
				index = (index > 0) ? index - 1 : 0;
				if (m_cdf[index] != m_cdf[index + 1]) {
					T dth = m_theta[index + 1] - m_theta[index];
					T alpha = ( target- m_cdf[index]) / (m_cdf[index + 1] - m_cdf[index]);
					return m_theta[index] + dth * alpha;
				}
				else {
					return m_theta[index];
				}
			}
		};
		// we want to make prob less than MaxProbInPin in each bin
		while (max_steps) {
			size_t new_points = 0;
			auto inplace_iter = m_points.begin() + tmp_size;
			for (size_t i = 0; i < tmp_size - 1; ++i) {
				PDF_t const& p = m_points[i];
				PDF_t const& nxt = m_points[i + 1];
				T h = nxt[0_c] - p[0_c];
				if (2 * std::abs(p[1_c] - nxt[1_c]) > MaxProbDiff * (p[1_c] + nxt[1_c])) {
					if (max_steps) [[likely]] {
						new_points++;
						max_steps--;
						PDF_t ph = new_point(p, nxt);
						*inplace_iter = ph;
						++inplace_iter;
					}
					else {
						push_points(m_points, tmp_size, new_points, m_buffer);
						return return_value(tmp_size);
					}
				}
			}
			if (new_points) {
				push_points(m_points, tmp_size, new_points, m_buffer);
			}
			else {
				return return_value(tmp_size);
			}
		}
		return return_value(tmp_size);
	}
};