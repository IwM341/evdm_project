#pragma once
#include "dynamics.hpp"
#include "../utils/progress_bar.hpp"
namespace evdm {
	/*MK generator of output vu', considering evaporation */
	template <class Gen_t, typename VectorT>
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> VuOut_evap(
		Gen_t&& G, const  VectorT& Vcm, const VectorT& Vu,
		Gen_vt<Gen_t> Vesc, Gen_vt<Gen_t> VescMin,
		Gen_vt<Gen_t>  mp, Gen_vt<Gen_t>  mk,
		Gen_vt<Gen_t>  mp_frac, Gen_vt<Gen_t>  mk_frac,
		Gen_vt<Gen_t>  m_cm, Gen_vt<Gen_t> deltaE = 0) {

		typedef Gen_vt<Gen_t> number_t;
		number_t VcmN = Vcm.norm();
		typedef vec3<number_t> vec3_t;
		vec3_t  n_v = (VcmN > 0 ? Vcm / VcmN : vec3_t({ 0,0,1 }));

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

		number_t Vu1_squared = Vu.squaredNorm() - deltaE * 2 / m_cm;
		if (Vu1_squared <= 0.0)
			return MCResult<vec3_t, number_t>(vec3_t(0, 0, 0), 0);

		number_t Vu1 = sqrt(Vu1_squared);

		number_t cosTh1 = 2 * G() - 1;
		number_t sinTh1 = sqrt(1.0 - cosTh1 * cosTh1);
		number_t phi1 = RandomPhi(G);

		const vec3_t Vu1_vec = Vu1 *
			(n_v * cosTh1 + n_1 * sinTh1 * cos(phi1) + n_2 * sinTh1 * sin(phi1));
		return MCResult<vec3_t, number_t>(Vu1_vec, Vu1 / VescMin);
	}

	template <
		typename Gen_t,
		typename dF_Type,
		typename VescR_t,
		typename NR_t,
		typename Temp_t
	>
	inline MCResult<
		vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>
	> Vout1_Scatter(
		vec3<Gen_vt<Gen_t>> V_wimp,
		Gen_t&& G,
		Gen_vt<Gen_t> mi, Gen_vt<Gen_t> mk,
		Gen_vt<Gen_t> mi_frac, Gen_vt<Gen_t> mk_frac,
		Gen_vt<Gen_t> m_cm, Gen_vt<Gen_t> delta_mk,
		dF_Type dF,
		VescR_t const VescR,
		Gen_vt<Gen_t> VescMin,
		NR_t const n_nd,
		Temp_t const TempR
	) {
		using T = Gen_vt<Gen_t>;
		T factor = 1;
		factor *= n_nd;

		//thermal generation
		auto Input_Vel_gen = Gauss3_BeyondD(G, std::sqrt(TempR / mi), 0.8, 8);
		factor *= Input_Vel_gen.RemainDensity;
		vec3<T> V1 = Input_Vel_gen.Result;
		//Vcm - is a vrlocity of momentum center
		vec3<T> Vcm = (V_wimp * mk_frac + V1 * mi_frac);
		//Nu is input velocity of WIMP in cm coordinatesd
		vec3<T> Vu = (V_wimp - V1);

		//Ecm - is kinetic energy in cm
		auto E_cm = m_cm * Vu.squaredNorm() / 2;

		//Energyloss is used in inelastic processes
		//auto EnLoss_mc = dF.EnergyLoss(E_cm - delta_mk);
		T inel_enloss = 0;
		//factor *= EnLoss_mc.RemainDensity;

		// Generating out velocity
		auto Vumk = VuOut_evap(
			G, Vcm, Vu, VescR, VescMin, mi, mk, mi_frac, mk_frac, m_cm,
			inel_enloss + delta_mk);
		vec3<T> Vu1 = Vumk.Result;
		factor *= Vumk.RemainDensity;

		// q - exchange momentum
		vec3<T> vec_Q = m_cm * (Vu - Vu1);
		auto q_2 = vec_Q.squaredNorm();

		//
		vec3<T> vec_V_T_inel = (Vu + Vu1) * (mi_frac - mk_frac) / 2;
		if (q_2 != 0)
			vec_V_T_inel += delta_mk * vec_Q / q_2;
		T v_2 = vec_V_T_inel.squaredNorm();
		//////////
		// v_2 = \vec_{v}_{inel T}^{\perp 2}
		// \vec_{v}_{inel T}^{\perp 2} = 1/2( v_x_1 + v_x_2 - v_N_1 - v_N_2)+delta/q^2*q
		// v_x_1 + v_x_2 - v_N_1 - v_N_2 = (nu+nu1)*(1-mk/mi)
		////////// 

		factor *= dF.ScatterFactor(q_2, v_2, inel_enloss);

		/**/
		return { vec3<T>(Vu1 * mi_frac + Vcm), factor };
	}


	template <
		typename Gen_t,
		typename HistoBank_t,
		typename HistoEvap_t,
		typename LEFunc_t,
		typename BinTrajPool_t,

		typename BinMes_t,

		typename PhiFunc_t,
		typename SFunc_t,
		typename TempFunc_t,
		typename VescFunc_t
	>
	inline void ScatterImpl(
		Gen_t&& G,
		Gen_vt<Gen_t> mk, Gen_vt<Gen_t> dm,
		Gen_vt<Gen_t> mi,
		ScatterEvent const& se,
		HistoBank_t& ScatMat_bank,
		HistoEvap_t&& EvapDistrib,
		bool count_evap,
		LEFunc_t const& LEf,
		Gen_vt<Gen_t> _F2_value,
		const BinTrajPool_t& TrajPoolVec,
		BinMes_t m_mes,
		PhiFunc_t const& Phi,
		SFunc_t const& S_func,

		TempFunc_t const& TempR,
		VescFunc_t const& VescR,
		Gen_vt<Gen_t> VescMin,
		size_t Nmk,
		size_t Nmk_per_traj,
		Gen_vt<Gen_t> weight,
		progress_omp_function<>& m_progress_func
	) {
		typedef
			bin_traj_p_vt_t<std::decay_t<decltype(TrajPoolVec[0])>>
			Traj_t;
		auto FromFactorVar_Action = [&, &Nir = se.n_e](auto&& dF) {
			const size_t N_in = ScatMat_bank.Grid.size();
			progress_omp_bar<> m_bar(
				m_progress_func, N_in, std::max((int)(N_in / 1000), 1)
			);
#pragma	omp parallel for
			auto m_cm = mk * mi / (mk + mi);
			auto mi_frac = mi / (mk + mi);
			auto mk_frac = mk / (mk + mi);
			for (size_t i = 0; i < N_in; ++i) {
				auto IJ = ScatMat_bank.Grid.FromLinear(i);

				auto el_bin = ScatMat_bank.Grid[IJ];

				auto OutHisto = ScatMat_bank.Values[i];

				const bin_traj_pool_t<Traj_t>& el_trajpool = TrajPoolVec[i];

				auto TinFunc = make_bin_traj_pool_Tin_func(
					el_trajpool, LEf, _F2_value);
				auto ToutFunc = make_bin_traj_pool_Tout_func(
					el_trajpool, LEf);

				auto m_bin_el_gen = gen_EL(m_mes, el_bin, G, LEf);


				for (size_t nm = 0; nm < Nmk; ++nm) {
					auto [e, l] = m_bin_el_gen();
					auto Lmax = LEf(-e);
					auto Ltmp = l * Lmax;

					auto [u0, u1, theta_max] =
						TinFunc.u0_u1_theta1(e, l, Lmax);
					u0 = std::clamp(u0, (decltype(u0))0, (decltype(u0))1);
					u1 = std::clamp(u1, (decltype(u1))0, (decltype(u1))1);
					auto r0 = std::sqrt(u0);
					auto r1 = std::sqrt(u1);
					auto [tp_i, tp_j] = el_trajpool.Grid.pos(e, l);

					auto const& th00 = el_trajpool[{tp_i, tp_j}].theta_tau;
					//prelimenary trajectory

					auto Tin_Teheta = TinFunc.tin_theta(e, l);

					auto Tin = Tin_Teheta * theta_max;
					auto Tout = ToutFunc(e, l, Ltmp);

					auto mk_factor =
						weight * Tin / ((Tin + Tout) * Nmk * Nmk_per_traj);

					for (size_t nt = 0; nt < Nmk_per_traj; ++nt) {
						auto tau = G();
						auto [theta_undim, d_theta_undim] = th00(tau);
						auto theta = theta_undim * theta_max;
						auto d_theta = d_theta_undim * theta_max;

						auto cth2 = std::cos(theta / 2);
						auto sth2 = std::sin(theta / 2);
						Traj_t r = std::sqrt(u0 * cth2 * cth2 + u1 * sth2 * sth2);

						auto renorm_factor =
							d_theta_undim /
							(2 * Tin_Teheta * std::sqrt(S_func(r, r0, r1)));

						Traj_t v2 = downbound(Phi(r) + e, 0);
						auto v = std::sqrt(v2);
						Traj_t vt = upbound(Ltmp / r, v);
						Traj_t vr = ssqrt(v2 - vt * vt);

						vec3<Traj_t > V_in(vt * VescMin, 0, vr * VescMin);
						Traj_t v_esc_nd = VescR(r);
						auto [Vout, factor] =
							Vout1_Scatter(
								V_in, G, mi, mk, mi_frac, mk_frac, m_cm, dm, dF,
								v_esc_nd * VescMin, VescMin,
								Nir(r), TempR(r)
							);

						vec3<Traj_t > v_nd = Vout / VescMin;//OK
						Traj_t e_out = (v_nd.squaredNorm() - v_esc_nd * v_esc_nd);


						auto final_factor = factor * mk_factor * renorm_factor;
						if (e_out < 0) {
							Traj_t L_nd = r * std::sqrt(v_nd.x() * v_nd.x() + v_nd.y() * v_nd.y());
							Traj_t l_out = (Lmax > 0 ? L_nd / Lmax : 0);
							OutHisto.put_force(final_factor, e_out, l_out);
						}
						else if (count_evap) {
							EvapDistrib.Values[i] += final_factor;
						}
					}
				}
				m_bar.next();
			}
		};
		std::visit(FromFactorVar_Action, se.sf);
	}


}; //namespace evdm 