#pragma once 

#include "dynamics.hpp"
#include "../utils/progress_bar.hpp"
#include "../measure.hpp"
namespace evdm {

	template <class Gen_t>
	/*MK generator of input velocity*/
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> HaloVelocity(
		Gen_t&& G,Gen_vt<Gen_t> VescTmp,
		Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> mU0) 
	{
		using T = Gen_vt<Gen_t>;
		auto ksi = sqrt(-2 * std::log( E0I1_G(G)) );
		auto phi = RandomPhi(G)/2;


		auto sinPhi = sin(phi);
		auto cosPhi = cos(phi);

		auto u0 = mU0 / Vdisp;
		auto ve = VescTmp / Vdisp;

		auto u = sqrt(u0 * u0 + ksi * ksi + 2 * u0 * ksi * cosPhi);


		auto sinTheta = (u != 0.0 ? ksi * sinPhi / u : 0);
		auto v = sqrt(u * u + ve * ve);
		
		auto n = RandomNvec(G);
		
		return MCResult<vec3<T>,T>(n * (v * Vdisp), 
			sinTheta * v * sqrt((T)std::numbers::pi / 2));
	}
	template <class Gen_t>
	/*MK generator of input velocity*/
	inline MCResult<vec3< Gen_vt<Gen_t>>, Gen_vt<Gen_t>> 
		HaloVelocityConstrained(
		Gen_t&& G, Gen_vt<Gen_t> VescTmp, Gen_vt<Gen_t> Vmin,
		Gen_vt<Gen_t> Vmax,Gen_vt<Gen_t>  Vdisp, Gen_vt<Gen_t>  mU0)
	{
		typedef Gen_vt<Gen_t> number_t;
		number_t umin = (Vmin > VescTmp ? sqrt(Vmin * Vmin - VescTmp * VescTmp) / Vdisp : 0.0);
		number_t umax = ((Vmax > Vmin && Vmax > VescTmp) ? sqrt(Vmax * Vmax - VescTmp * VescTmp) / Vdisp : umin);

		number_t u0 = mU0 / Vdisp;

		number_t ksi_max = umax + u0;
		number_t ksi_min = (mU0 - umax > 0 ? u0 - umax : (umin - u0 > 0 ? umin - u0 : 0));

		number_t ksi_rd = exp(-ksi_min * ksi_min / 2) - exp(-ksi_max * ksi_max / 2);

		number_t ksi = sqrt(-2 * log(exp(-ksi_max * ksi_max / 2) + G() * ksi_rd));

		number_t cosThMin = std::min(1.0, std::max(-1.0, (umin * umin - ksi * ksi - u0 * u0) / (2 * ksi * u0)));
		number_t cosThMax = std::max(-1.0, std::min(1.0, (umax * umax - ksi * ksi - u0 * u0) / (2 * ksi * u0)));

		number_t max_theta = acos(cosThMin);
		number_t min_theta = acos(cosThMax);

		number_t rd_th = (max_theta - min_theta) / (number_t)std::numbers::pi;

		number_t theta = min_theta + (max_theta - min_theta) * G();

		number_t ve = VescTmp / Vdisp;

		number_t u = sqrt(u0 * u0 + ksi * ksi + 2 * u0 * ksi * cos(theta));

		number_t sinTheta = (u != 0.0 ? ksi * sin(theta) / u : 0);

		auto v = sqrt(u * u + ve * ve);
		auto n = RandomNvec(G);

		return { n * (v * Vdisp), 
			ksi_rd * rd_th * sinTheta * v * sqrt((number_t)std::numbers::pi/ 2) 
		};
	}

	template <class Gen_t,typename VectorT>
	/*MK generator of output nu'*/
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> VuOut(
		Gen_t && G, const  VectorT& Vcm, const VectorT& Vu,
		Gen_vt<Gen_t> Vesc, Gen_vt<Gen_t> VescMin, 
		Gen_vt<Gen_t>  mp, Gen_vt<Gen_t>  mk, 
		Gen_vt<Gen_t>  mp_frac, Gen_vt<Gen_t>  mk_frac,
		Gen_vt<Gen_t>  m_cm, Gen_vt<Gen_t> deltaE = 0) {

		typedef Gen_vt<Gen_t> number_t;
		number_t VcmN = Vcm.norm();
		typedef vec3<number_t> vec3_t;
		vec3_t  n_v = (VcmN > 0 ? Vcm / VcmN : vec3_t({0,0,1}) );

		number_t cosThetaVcm = n_v[2];
		number_t sinThetaVcm = sqrt(n_v[0] * n_v[0] + n_v[1] * n_v[1]);

		number_t cosPhiVcm = 1.0;
		number_t sinPhiVcm = 0.0;

		if (sinThetaVcm > 1e-10) {
			cosPhiVcm = n_v[0] / sinThetaVcm;
			sinPhiVcm = n_v[1] / sinThetaVcm;
		}

		vec3_t n_1(cosThetaVcm * cosPhiVcm, cosThetaVcm * sinPhiVcm, -sinThetaVcm);
		vec3_t n_2(-sinPhiVcm, cosPhiVcm, 0);

		/*
		std::cout << n_1*n_1 << "\t" << n_1*n_2 << "\t" << n_1*n_v << std::endl;
		std::cout << n_2*n_1 << "\t" << n_2*n_2 << "\t" << n_2*n_v << std::endl;
		std::cout << n_v*n_1 << "\t" << n_v*n_2 << "\t" << n_v*n_v << std::endl<< std::endl;
		*/

		number_t Vu1_squared = Vu.squaredNorm() - deltaE * 2/ m_cm;
		if (Vu1_squared <= 0.0)
			return MCResult<vec3_t, number_t>(vec3_t(0, 0, 0), 0);

		number_t Vu1 = sqrt(Vu1_squared);

		number_t cosTh1max = (
			Vesc * Vesc - Vu1_squared * mp_frac * mp_frac - VcmN * VcmN
			) / (2 * VcmN * Vu1* mp_frac);

		if (!(cosTh1max > -1))
			return MCResult<vec3_t, number_t>(vec3_t(0, 0, 0), 0);
		else if (cosTh1max >= 1) {
			cosTh1max = 1;
		}

		number_t cosTh1 = (1 + cosTh1max) * G() - 1;
		number_t sinTh1 = sqrt(1.0 - cosTh1 * cosTh1);
		number_t phi1 = RandomPhi(G);

		const vec3_t Vu1_vec = Vu1 * (n_v * cosTh1 + n_1 * sinTh1 * cos(phi1) + n_2 * sinTh1 * sin(phi1));
		return MCResult<vec3_t, number_t>(Vu1_vec, (1 + cosTh1max)/2 * Vu1 / VescMin);
	}


	/// <summary>
	/// generate out velocity (with weight) in capture.
	/// </summary>
	/// <param name="mp">nuclei mass</param>
	/// <param name="mk">DM mass</param>
	/// <param name="delta_mk">DM delta mass </param>
	/// <param name="dF">form factor function</param>
	/// <param name="VescR"></param>
	/// <param name="VescMin"></param>
	/// <param name="nR"></param>
	/// <param name="TempR"></param>
	/// <param name="G"></param>
	/// <param name="Vdisp"></param>
	/// <param name="mU0"></param>
	/// <param name="pow_r"></param>
	/// <returns>tuple(\vec{V}_{out},r - location, Vdisp - disp velocity at r)</returns>
	template <
		typename Gen_t,
		typename dF_Type,
		typename VescRFuncType,
		typename N_FuncType,
		typename TempRFuncType
		>
	inline MCResult<
		std::tuple<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>, Gen_vt<Gen_t>>,
		Gen_vt<Gen_t>
	> Vout1(
		Gen_t&& G,
		Gen_vt<Gen_t> mi, Gen_vt<Gen_t> mk, 
		Gen_vt<Gen_t> mi_frac, Gen_vt<Gen_t> mk_frac,
		Gen_vt<Gen_t> m_cm, Gen_vt<Gen_t> delta_mk,
		dF_Type dF,
		VescRFuncType const& VescR, 
		Gen_vt<Gen_t> VescMin,
		N_FuncType const& nR, 
		TempRFuncType const& TempR,
		Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> mU0, 
		Gen_vt<Gen_t> pow_r = 1) 
	{
		using T = Gen_vt<Gen_t>;
		T factor = 1;

		//generate radius
		T r_nd = pow(G(), pow_r);//pow(G(),1.0/3.0);
		factor *= (3 * pow_r * pow(r_nd, (3 * pow_r - 1) / pow_r));
		//gain escape velocity from redius
		//if(r_nd < 0.3){
		//    r_nd+=0.0;
		//}
		auto Vesc = VescR(r_nd)*VescMin;

		//random input velocity
		auto VelocityMk = HaloVelocity(G, Vesc, Vdisp, mU0);
		auto V_wimp = VelocityMk.Result;//vec3::PolarCos(sqrt(VelocityMk.Result*VelocityMk.Result+Vesc*Vesc),
		//RandomCos(G),RandomPhi(G));
		factor *= VelocityMk.RemainDensity;


		auto n_nd = nR(r_nd);//TODO n_nd as a function of radius
		factor *= n_nd;
		
		//thermal generation
		auto Input_Vel_gen = Gauss3_BeyondD(G, std::sqrt(TempR(r_nd) / mi),0.8,8);
		
		factor *= Input_Vel_gen.RemainDensity;
		vec3<T> V1 = Input_Vel_gen.Result;

		//Vcm - is a vrlocity of momentum center
		vec3<T> Vcm = (V_wimp * mk_frac + V1 * mi_frac);
		//Nu is input velocity of WIMP in cm coordinatesd
		vec3<T> Vdelta =  (V_wimp - V1);

		//Ecm - is kinetic energy in cm
		auto E_cm = m_cm * Vdelta.squaredNorm() / 2;


		//Energyloss is used in inelastic processes
		//auto EnLoss_mc = dF.EnergyLoss(E_cm - delta_mk);
		T inel_enloss = 0;
		//factor *= EnLoss_mc.RemainDensity;

		// Generating out velocity
		auto Vumk = VuOut(G, Vcm, Vdelta, Vesc, VescMin, 
			mi, mk,mi_frac,mk_frac,m_cm, inel_enloss + delta_mk
		);
		vec3<T> Vu1 = Vumk.Result;
		factor *= Vumk.RemainDensity;

		// q - exchange momentum
		vec3<T> vec_Q = m_cm * (Vdelta - Vu1);
		auto q_2 = vec_Q.squaredNorm();

		//
		vec3<T> vec_V_T_inel = (Vdelta + Vu1) / 2;
		if(q_2 != 0) 
			vec_V_T_inel -= delta_mk * vec_Q / q_2;
		
		T v_2 = vec_V_T_inel.squaredNorm();
		/*if ((Vdelta + Vu1).squaredNorm() != 0.0 && v_2 == 0.0) {
			throw std::runtime_error("Eigen bug");
		}*/
		
		//////////
		// v_2 = \vec_{v}_{inel T}^{\perp 2}
		// \vec_{v}_{inel T}^{\perp 2} = 1/2( v_x_1 + v_x_2 - v_N_1 - v_N_2)+delta/q^2*q
		// v_x_1 + v_x_2 - v_N_1 - v_N_2 = (nu+nu1)*(1-mk/mi)
		////////// 

		factor *= dF.ScatterFactor(q_2, v_2,inel_enloss);

		using Tuple_t = std::tuple<vec3<T>, T, T>;
		return MCResult<Tuple_t,T>(Tuple_t(Vu1*mi_frac + Vcm, r_nd, Vesc), factor);
	}


	/// <summary>
	/// calculates capture.
	/// </summary>
	/// <param name="G"></param>
	/// <param name="se"></param>
	/// <param name="TempR"></param>
	/// <param name="H_le_sp_type"></param>
	/// <param name="mk"></param>
	/// <param name="dm"></param>
	/// <param name="mi"></param>
	/// <param name="_mp"></param>
	/// <param name="_3p_r"></param>
	/// <param name="mU0"></param>
	/// <param name="Vdisp"></param>
	/// <param name="VescMin"></param>
	/// <param name="Nmk"></param>
	template <
		typename Histo_t,
		//typename Histo_opt_t,
		typename LEFunc_t,
		typename Rm_type_t,
		typename Gen_t,
		typename TempFunc_t,
		typename VescFunc_t>
	inline std::pair<double,double> CaptureImpl(
		Gen_t&& G,
		ScatterEvent const & se,
		TempFunc_t const& TempR,
		Histo_t &H_le_sp_type,
		LEFunc_t const & LEf,
		Gen_vt<Gen_t> mk, Gen_vt<Gen_t> dm,
		Gen_vt<Gen_t> mi,
		Rm_type_t _3p_r,
		Gen_vt<Gen_t> mU0,
		Gen_vt<Gen_t> Vdisp,
		VescFunc_t const &VescR,
		Gen_vt<Gen_t> VescMin,
		size_t Nmk,
		Gen_vt<Gen_t> weight = 1)
	{
		using T = Gen_vt<Gen_t>;
		T sum = 0;
		T sum2 = 0;
		auto const_fact_rd = weight / Nmk;
		auto m_cm = (mk * mi) / (mk + mi);
		auto mi_frac = (mi) / (mk + mi);
		auto mk_frac = (mk) / (mk + mi);
		auto action = [&](auto&& dF) {
			for (size_t i = 0; i < Nmk; ++i) {
				auto mk_res = Vout1(
					G,mi, mk, mi_frac, mk_frac,m_cm, dm, 
					dF, VescR, VescMin, se.n_e, TempR, Vdisp, mU0, _3p_r);
				vec3<T> v_nd = std::get<0>(mk_res.Result) / VescMin;
				

				T r_nd = std::get<1>(mk_res.Result);
				T  v_esc_nd = std::get<2>(mk_res.Result) / VescMin;

				T E_nd = (v_nd.squaredNorm() - v_esc_nd * v_esc_nd);
				
				if (E_nd < 0) {
					T L_nd = r_nd * sqrt(v_nd.x() * v_nd.x() + v_nd.y() * v_nd.y());
					T l = L_nd / LEf(-E_nd);
					T dens = mk_res.RemainDensity * const_fact_rd;;
					H_le_sp_type.put_force(dens, E_nd, l);
					sum += dens;
					sum2 += dens * mk_res.RemainDensity;
				}/*
				else{
					auto l = L_nd/H.LE_func(E_nd);
				}*/
				/*
			   if(Ewas != H.values[1].values[0]){
				   PVAR(dens);
				   PVAR(v_nd);
				   PVAR(r_nd);
				   PVAR(v_esc_nd);
				   PVAR(E_nd);
				   PVAR(L_nd);
				   PVAR(H.values[1].values[1]);
				   print();
			   }*/
			   /*
		   if(std::isnan(E_nd) or std::isnan(L_nd)){
			   PVAR(r_nd);
			   PVAR(v_nd);
			   PVAR(v_esc_nd);
			   PVAR(v_esc_nd);
		   }*/
			}
		};
		std::visit(action,se.sf);

		return {sum,std::sqrt(sum2* weight - sum * sum )};
	}
};

