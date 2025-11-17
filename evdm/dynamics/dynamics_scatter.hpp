#pragma once
#include "dynamics.hpp"
#include "../utils/progress_bar.hpp"
#include <ranges>
#include "../core/core_matrix.hpp"
#include <deque>
#include <numeric>
//#include <>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_MSC_VER)
#define NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define NOINLINE __attribute__((noinline))
#else
#define NOINLINE
#endif

namespace evdm {

	template <typename T, typename Func_t>
	inline T critical_value(T left, T right, Func_t&& F, bool expect_first_true, size_t Nstep = 20) {
		if (expect_first_true) {
			if (!F(left)) {
				return left;
			} if (F(right)) {
				return right;
			}
			for (size_t i = 0; i < Nstep; ++i) {
				T midl = (left + right) / 2;
				if (F(midl)) {
					left = midl;
				}
				else {
					right = midl;
				}
			}
			return right;
		}
		else {
			if (F(left)) {
				return left;
			} if (!F(right)) {
				return right;
			}
			for (size_t i = 0; i < Nstep; ++i) {
				T midl = (left + right) / 2;
				if (F(midl)) {
					right = midl;
				}
				else {
					left = midl;
				}
			}
			return left;
		}
	}

	template <typename inner_iterator_t>
	struct TripletView {
		int RowShift, ColShift;
		inner_iterator_t m_begin, m_end;
		typedef inner_iterator_t inner_iterator;
		typedef typename inner_iterator::value_type SpVec_t;
		typedef typename std::decay_t<SpVec_t>::InnerIterator sp_iterator;
		typedef typename std::decay_t<SpVec_t>::Scalar inner_value_type;
		typedef Eigen::Triplet<inner_value_type> value_type;

		struct iterator : std::iterator<std::output_iterator_tag, value_type> {

			friend struct TripletView;

			inner_iterator vec_iter;
			inner_iterator vec_end;

			sp_iterator sp_iter;
			int RowShift, ColShift;
			int j;

		private:
			iterator(inner_iterator vec_iter,
				inner_iterator vec_end,
				sp_iterator sp_iter,
				int RowShift, int ColShift, int j
			) :vec_iter(vec_iter), vec_end(vec_end), sp_iter(sp_iter),
				RowShift(RowShift), ColShift(ColShift), j(j) {}

		public:
			inline bool operator !=(iterator other)const {
				return (vec_iter != other.vec_iter);
			}
			inline bool operator ==(iterator other) const {
				return (vec_iter == other.vec_iter);
			}
			iterator& operator++() {
				if (sp_iter) {
					++sp_iter;
				}
				if (sp_iter) {
					return *this;
				}
				else {
					++vec_iter;
					for (; vec_iter != vec_end; ++vec_iter) {
						if ((*vec_iter).nonZeros()) {
							break;
						}
						++j;
					}

					if (vec_iter != vec_end) {
						sp_iter = sp_iterator(*vec_iter);
						++j;
					}
				}
				return *this;
			}
			inline value_type operator*()const {
				return { sp_iter.index() + ColShift,j + RowShift,sp_iter.value() };
			}
			struct proxy {
				value_type m_value;
				inline value_type* operator->() noexcept{
					return &m_value;
				}
				inline value_type& operator *() noexcept {
					return m_value;
				}
			};
			inline proxy operator->()const noexcept {
				return { operator*() };
			}
		};
		TripletView(inner_iterator_t _begin, inner_iterator_t _end,
			int ColShift, int RowShift) :m_begin(_begin), m_end(_end), 
			ColShift(ColShift), RowShift(RowShift) {
		}

		iterator begin() {
			sp_iterator it;
			int j = 0;
			for (; m_begin != m_end; ++m_begin) {
				if ((*m_begin).nonZeros()) {
					it = sp_iterator(*m_begin);
					break;
				}
				++j;
			}
			return { m_begin,m_end,it,ColShift,RowShift,j };
		}
		iterator end() {
			sp_iterator it;
			return { m_end,m_end,it,ColShift,RowShift,0 };
		}
	};

	template <typename iter_t>
	TripletView<iter_t> triplet_view(
		iter_t begin, iter_t end,  int ColShift = 0, int RowShift = 0) {
		return { begin, end,RowShift, ColShift };
	}


	/*MK generator of output vu', considering evaporation */
	template <class Gen_t, typename VectorT>
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> VuOut_evap(
		Gen_t&& G, const  VectorT& Vcm, const VectorT& Vu,
		Gen_vt<Gen_t> Vesc, Gen_vt<Gen_t> VescMin,
		Gen_vt<Gen_t>  mp, Gen_vt<Gen_t>  mk,
		Gen_vt<Gen_t>  mp_frac, Gen_vt<Gen_t>  mk_frac,
		Gen_vt<Gen_t>  m_cm, Gen_vt<Gen_t> deltaE_div_m_cm = 0) {

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

		number_t Vu1_squared = Vu.squaredNorm() - deltaE_div_m_cm;
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
		typename Temp_t,
		typename ThermGen_t
	>
	inline MCResult<
		vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>
	> DeltaVout1_Scatter(
		vec3<Gen_vt<Gen_t>> V_wimp,
		Gen_vt<Gen_t> V_wimp_norm,
		Gen_t& G,
		ThermGen_t ThrmGen,
		Gen_vt<Gen_t> mi, Gen_vt<Gen_t> mk,
		Gen_vt<Gen_t> mi_frac, Gen_vt<Gen_t> mk_frac,
		Gen_vt<Gen_t> m_cm, Gen_vt<Gen_t> delta_mk,
		Gen_vt<Gen_t> deltaE_div_m_cm,
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
		
		MCResult< T, T> Input_Vel_gen_abs = 
			ThrmGen.gen_abs(G, deltaE_div_m_cm, V_wimp_norm, TempR, mi);
		factor *= Input_Vel_gen_abs.RemainDensity;
		T V1_abs = Input_Vel_gen_abs.Result;

		T max_cos = upbound(
			downbound(
				(V_wimp_norm * V_wimp_norm + 
					V1_abs * V1_abs - 
					deltaE_div_m_cm
				) / (2 * V_wimp_norm * V1_abs),
			-1), 1);

		
		factor *= (max_cos + 1) / 2;
		if (factor == (T)0) {
			return { vec3<T>(0, 0, 0), 0};
		}
		vec3<T> V1 = GenVecCos(G, V_wimp, -1, max_cos) * V1_abs;
		auto norm = V1.norm();

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
			deltaE_div_m_cm);
		vec3<T> Vu1 = Vumk.Result;
		factor *= Vumk.RemainDensity;

		// q - exchange momentum
		vec3<T> dVu = Vu1 - Vu;
		auto vec_minusQ = m_cm * dVu;
		auto q_2 = vec_minusQ.squaredNorm();

		//
		vec3<T> vec_V_T_inel = (Vu + Vu1) / 2;
		if (delta_mk != 0 && q_2 != 0)
			vec_V_T_inel += delta_mk * vec_minusQ / q_2;
		T v_2 = vec_V_T_inel.squaredNorm();
		//////////
		// v_2 = \vec_{v}_{inel T}^{\perp 2}
		// \vec_{v}_{inel T}^{\perp 2} = 1/2( v_x_1 + v_x_2 - v_N_1 - v_N_2)+delta/q^2*q
		// v_x_1 + v_x_2 - v_N_1 - v_N_2 = (nu+nu1)*(1-mk/mi)
		////////// 

		factor *= dF.ScatterFactor(q_2, v_2, inel_enloss);

		/**/
		return { vec3<T>(dVu * mi_frac), factor };
	}

	template <typename T,typename Func_Phi_t, typename Func_Tr_t,typename TinTraj_t, typename LEf_t>
	bin_dedl_t<T> PrepairBin(bin_dedl_t<T> m_bin,
		evdm::same_t<T> m, evdm::same_t<T> mi, evdm::same_t<T> deltaE_div_m_cm,
		Func_Phi_t &&Phi, Func_Tr_t &&TempR, TinTraj_t && TinFunc, LEf_t&&LEf)
	{
		auto reaction_goes_e = [&](auto e) {
			auto [u0, u1, theta_max] =
				TinFunc.u0_u1_theta1(e, std::get<1>(m_bin).left, LEf(-e));
			T r = std::sqrt(u0);
			T v2 = downbound(Phi(r) + e, 0);
			auto Tr = TempR(r);
			auto Vmax = std::sqrt(v2) + 8 * std::sqrt(Tr / mi);
			return Vmax >= ssqrt(deltaE_div_m_cm);
		};
		auto& e0e1 = std::get<0>(m_bin);
		e0e1.left = critical_value(e0e1.left, e0e1.right, reaction_goes_e, false);
		
		auto &l0l1 = std::get<1>(m_bin);
		T critical_l_worst = std::get<1>(m_bin).left;
		for (size_t lb = 0; lb < 2; ++lb) {
			auto e = lb ? e0e1.left : e0e1.right;
			auto reaction_goes_l = [&](auto l) {
				auto Le = LEf(-e);
				auto [u0, u1, theta_max] =
					TinFunc.u0_u1_theta1(e, l, Le);
				T r = std::sqrt(u0);
				T v2 = downbound(Phi(r) + e, 0);
				auto Tr = TempR(r);
				auto Vmax = std::sqrt(v2) + 8 * std::sqrt(Tr / mi);
				return Vmax >= ssqrt(deltaE_div_m_cm);
			};
			critical_l_worst = std::max(
				critical_l_worst, 
				critical_value(l0l1.left, l0l1.right, reaction_goes_l, true)
			);
		}
		l0l1.right = critical_l_worst;
		return m_bin;
	}

	template < typename Te_t,typename U, typename V,typename T,
		typename Func_Phi_t, typename Func_Tr_t, typename TinTraj_t>
	auto find_tau_max(Te_t e,U r0,U r1, V theta_max,
		T m, T mi, evdm::same_t<T> deltaE_div_m_cm,
		Func_Phi_t&& Phi, Func_Tr_t&& TempR, TinTraj_t&& th00) {
		
		typedef  decltype(th00.Grid[0]) Tau_t;
		auto u0 = r0 * r0;
		auto u1 = r1 * r1;
		auto reaction_goes = [&](auto tau) {
			auto [theta_undim, d_theta_undim] = th00(tau);
			auto theta = theta_undim * theta_max;

			auto cth2 = std::cos(theta / 2);
			auto sth2 = std::sin(theta / 2);
			auto r = std::sqrt(u0 * cth2 * cth2 + u1 * sth2 * sth2);
			auto v2 = downbound(Phi(r) + e, 0);
			auto Tr = TempR(r);
			auto Vmax = std::sqrt(v2) + 8 * std::sqrt(Tr / mi);
			return Vmax >= ssqrt(deltaE_div_m_cm);
		};
		return critical_value(Tau_t(0), Tau_t(1), reaction_goes, true);
	}

	
	struct ScatterImplExtra{
		bool cut_el = false;
		bool cut_tau = false;
		size_t cut_factor = 8;
	};

	template <
		typename Grid_vt,
		typename Mat_vt,
		typename Gen_t,
		typename ScatterFactor_t,
		typename NDensityEvent_t,
		typename ThermGenerator_t,
		typename HistoEvap_t,
		typename LEFunc_t,
		typename BinTrajPool_t,


		typename PhiFunc_t,
		typename SFunc_t,
		typename TempFunc_t,
		typename VescFunc_t,
		typename Nmk_Vec_t
	>
	NOINLINE SpMatrix_t<Mat_vt> ScatterImpl(
		std::type_identity<Grid_vt>,
		std::type_identity<Mat_vt>,
		size_t MatSize,
		Gen_t&& _G,
		ThermGenerator_t ThermGen,
		Gen_vt<Gen_t> mk, Gen_vt<Gen_t> dm,
		Gen_vt<Gen_t> mi,
		ScatterFactor_t const & dF,
		NDensityEvent_t const & Nir,
		size_t MatrixIndexShift_In,
		size_t MatrixIndexShift_Out,
		HistoEvap_t&& EvapDistrib,
		Gen_vt<Gen_t> ZeroValue,
		LEFunc_t const& LEf,
		Gen_vt<Gen_t> _F2_value,
		const BinTrajPool_t& TrajPoolVec,
		measure_dEpdlq<Gen_vt< Gen_t>> m_mes,
		PhiFunc_t const& Phi,
		SFunc_t const& S_func,

		TempFunc_t const& TempR,
		VescFunc_t const& VescR,
		Gen_vt<Gen_t> VescMin,
		Nmk_Vec_t const & Nmk_v,
		size_t Nmk_per_traj,
		Gen_vt<Gen_t> weight,
		ScatterImplExtra extra_params = {}
	) {
		typedef Grid_vt Traj_t;

		
		const size_t N_in = EvapDistrib.Grid.size();
		auto&& grid = EvapDistrib.Grid;
		
		auto m_cm = mk * mi / (mk + mi);

		auto deltaE_div_m_cm = 2 * dm / m_cm;

		auto mi_frac = mi / (mk + mi);
		auto mk_frac = mk / (mk + mi);
		const auto _seed = _G.state;
		std::vector<Mat_vt> OutIBuffer(N_in);
		std::vector<std::vector<SpTriplet_t<Mat_vt>>> triplets(N_in);
		#pragma	omp parallel for  firstprivate(OutIBuffer)
		for (int i = 0; i < N_in; ++i) {
			std::fill(OutIBuffer.begin(), OutIBuffer.end(), 0);
			auto G = _G;
			G.set_seed(_seed ^ i + 1);
			auto IJ = grid.FromLinear(i);

			auto el_bin = grid[IJ];

			

			auto OutHisto = grob::make_histo_view(grid,OutIBuffer);
			const bin_traj_pool_t<Traj_t>& el_trajpool = TrajPoolVec[i];

			auto TinFunc = make_bin_traj_pool_Tin_func(
				el_trajpool, LEf, _F2_value);
			auto ToutFunc = make_bin_traj_pool_Tout_func(
				el_trajpool, LEf);
			
			Traj_t bin_cut_factor = 1;
			if(extra_params.cut_el){
				Traj_t old_mes = m_mes(el_bin);
				el_bin = PrepairBin(el_bin, mk, mi, deltaE_div_m_cm, Phi, TempR, TinFunc, LEf);
				Traj_t new_mes = m_mes(el_bin);
				bin_cut_factor = new_mes/old_mes;
			}

			auto m_bin_el_gen = m_mes.get_gen(el_bin);

			size_t Nmk = 0;
			if constexpr (std::is_same_v<Nmk_Vec_t, size_t>) {
				Nmk = Nmk_v;
			} else {
				Nmk = Nmk_v[i];
			}
			for (size_t nm = 0; nm < Nmk; ++nm) {
				
				auto [e, l] = m_bin_el_gen(G);
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


				Traj_t tau_max = (Nmk_per_traj <= 10 || !extra_params.cut_tau) ? Traj_t(1) :
					find_tau_max(e,r0, r1, theta_max, mk, mi, deltaE_div_m_cm, Phi, TempR, th00);

				auto Tin = Tin_Teheta * theta_max;
				auto Tout = ToutFunc(e, l, Ltmp);

				auto mk_factor =
					bin_cut_factor*weight * tau_max* Tin / ((Tin + Tout) * Nmk * Nmk_per_traj);

				if (tau_max > 0) {
					for (size_t nt = 0; nt < Nmk_per_traj; ++nt) {
						auto tau = G() * tau_max;
						auto [theta_undim, d_theta_undim] = th00(tau);
						auto theta = theta_undim * theta_max;
						auto d_theta = d_theta_undim * theta_max;

						auto cth2 = std::cos(theta / 2);
						auto sth2 = std::sin(theta / 2);
						Traj_t r = std::sqrt(u0 * cth2 * cth2 + u1 * sth2 * sth2);

#ifndef NDEBUG
						auto _C0 = Phi(std::sqrt(u0)) + e - Ltmp * Ltmp / u0;
						auto _C0_p = Phi(std::sqrt(u0 * (1 + 1e-4))) + e - Ltmp * Ltmp / (u0 * (1 + 1e-4));
						auto _C0_m = Phi(std::sqrt(u0 * (1 - 1e-4))) + e - Ltmp * Ltmp / (u0 * (1 - 1e-4));
						auto _C1 = Phi(std::sqrt(u1)) + e - Ltmp * Ltmp / u1;
						auto _C1_p = Phi(std::sqrt(u1 * (1 + 1e-4))) + e - Ltmp * Ltmp / (u1 * (1 + 1e-4));
						auto _C1_m = Phi(std::sqrt(u1 * (1 - 1e-4))) + e - Ltmp * Ltmp / (u1 * (1 - 1e-4));
#endif
						auto renorm_factor =
							d_theta_undim /
							(2 * Tin_Teheta * std::sqrt(S_func(r, r0, r1)));

						Traj_t v2 = downbound(Phi(r) + e, 0);
						auto v = std::sqrt(v2);
						Traj_t vt = upbound(Ltmp / r, v);
						Traj_t vr = ssqrt(v2 - vt * vt);

						vec3<Traj_t > V_in_nd(vt, 0, vr);
						vec3<Traj_t > V_in = V_in_nd * VescMin;
						Traj_t v_esc_nd = VescR(r);
						auto [dVout, factor] =
							DeltaVout1_Scatter(
								V_in, v, G, ThermGen,
								mi, mk, mi_frac, mk_frac, m_cm, dm, deltaE_div_m_cm,
								dF,
								v_esc_nd * VescMin, VescMin,
								Nir(r), TempR(r)
							);
						vec3<Traj_t > dVout_nd = dVout * (1 / VescMin);
						auto VplusV1 = 2 * V_in_nd + dVout_nd;
						auto dVtau2 = dVout_nd[0] * VplusV1[0] +
							dVout_nd[1] * VplusV1[1];

						Traj_t deltaE = dVtau2 + dVout_nd[2] * VplusV1[2];
						//vec3<Traj_t > v_nd = Vout / VescMin;//OK
						Traj_t e_out = e + deltaE;
						//Traj_t e_out1 = ((V_in_nd+ dVout_nd).squaredNorm() - v_esc_nd * v_esc_nd);

						auto final_factor = factor * mk_factor * renorm_factor;
						auto m_zero = ZeroValue / Nmk;
						if (e_out < 0) {
							if (final_factor > m_zero) {
								Traj_t deltaL2 = (r * r) * dVtau2;
								Traj_t L_nd2 = Ltmp * Ltmp + deltaL2;
								Traj_t L_nd = ssqrt(L_nd2);

								auto Lcutter = [](Traj_t Lmax) ->Traj_t {
									return Lmax > 0 ? 1 / Lmax : 0;
								};

								Traj_t Lm_inv = Lcutter(Lmax);
								Traj_t Lmax_out = LEf(-e_out);
								Traj_t Lm_inv1 = Lcutter(Lmax_out);

								Traj_t l_out = upbound(L_nd * Lm_inv1, 1);
								//Traj_t l_out1= (Lmax_out > 0 ? L_nd / Lmax_out : 0);
								OutHisto.put_force(final_factor, e_out, l_out);
							}
						}
						else {
							EvapDistrib.Values[i] += final_factor;
						}
					}
				}
			}

			size_t NonZeros = std::count_if(
				OutIBuffer.begin(), OutIBuffer.end(),
				[ZeroValue](auto x) {return x > ZeroValue; }
			);
			triplets[i].reserve(NonZeros);
			for (size_t k = 0; k < OutIBuffer.size(); ++k) {
				if (OutIBuffer[k] > 0) {
					triplets[i].push_back(SpTriplet_t<Mat_vt>(
						k + MatrixIndexShift_Out ,
						i + MatrixIndexShift_In,
						OutIBuffer[k] 
					));
				}
			}
		}
		std::vector<SpTriplet_t<Mat_vt>>
			AllTriplets(flatten(triplets));

		SpMatrix_t<Mat_vt> Ret(MatSize, MatSize);
		Ret.setFromTriplets(AllTriplets.begin(),AllTriplets.end());

		return Ret;
	}


	template <
		typename Grid_vt,
		typename Mat_vt,
		typename Gen_t,
		typename ScatterFactor_t,
		typename NDensityEvent_t,
		typename ThermGenerator_t,
		typename HistoEvap_t,
		typename LEFunc_t,
		typename BinTrajPool_t,

		typename PhiFunc_t,
		typename SFunc_t,
		typename TempFunc_t,
		typename VescFunc_t,
		typename Nmk_Vec_t
	>
	NOINLINE SpMatrix_t<Mat_vt> ScatterImplShift (
		std::type_identity<Grid_vt>,
		std::type_identity<Mat_vt>,
		size_t MatSize,
		Gen_t&& _G,
		ThermGenerator_t ThermGen,
		Gen_vt<Gen_t> mk, Gen_vt<Gen_t> dm,
		Gen_vt<Gen_t> mi,
		ScatterFactor_t const& dF,
		NDensityEvent_t const& Nir,
		size_t MatrixIndexShift_In,
		size_t MatrixIndexShift_Out,
		HistoEvap_t&& EvapDistrib,
		Gen_vt<Gen_t> ZeroValue,
		LEFunc_t const& LEf,
		Gen_vt<Gen_t> _F2_value,
		const BinTrajPool_t& TrajPoolVec,
		measure_dEpdlq<Grid_vt> m_mes,
		PhiFunc_t const& Phi,
		SFunc_t const& S_func,

		TempFunc_t const& TempR,
		VescFunc_t const& VescR,
		Gen_vt<Gen_t> VescMin,
		Nmk_Vec_t const& Nmk_v,
		size_t Nmk_per_traj,
		Gen_vt<Gen_t> weight,
		ScatterImplExtra extra_params = {}
	) {
		typedef Grid_vt Traj_t;


		const size_t N_in = EvapDistrib.Grid.size();
		auto&& grid = EvapDistrib.Grid;

		auto m_cm = mk * mi / (mk + mi);

		auto deltaE_div_m_cm = 2 * dm / m_cm;

		auto mi_frac = mi / (mk + mi);
		auto mk_frac = mk / (mk + mi);
		const auto _seed = _G.state;


		std::vector<Eigen::SparseVector<Mat_vt>> Columns(N_in);
		
#pragma	omp parallel for
		for (int i = 0; i < N_in; ++i) {
			
			Columns[i] = Eigen::SparseVector<Mat_vt>(N_in);

			auto G = _G;
			G.set_seed(_seed ^ i + 1);
			auto IJ = grid.FromLinear(i);

			auto el_bin = grid[IJ];
			auto [X0X1, Y0Y1] = m_mes.toXY(el_bin);

			//auto OutHisto = grob::make_histo_view(grid, OutIBuffer);
			const bin_traj_pool_t<Traj_t>& el_trajpool = TrajPoolVec[i];

			auto TinFunc = make_bin_traj_pool_Tin_func(
				el_trajpool, LEf, _F2_value);
			auto ToutFunc = make_bin_traj_pool_Tout_func(
				el_trajpool, LEf);
			
			Traj_t bin_cut_factor = 1;
			if(extra_params.cut_el){
				Traj_t old_mes = m_mes(el_bin);
				el_bin = PrepairBin(el_bin, mk, mi, deltaE_div_m_cm, Phi, TempR, TinFunc, LEf);
				Traj_t new_mes = m_mes(el_bin);
				bin_cut_factor = new_mes/old_mes;
			}
			//auto m_bin_el_gen = gen_EL(m_mes, el_bin, G, LEf);

			size_t Nmk = 0;
			if constexpr (std::is_same_v<Nmk_Vec_t, size_t>) {
				Nmk = Nmk_v;
			}
			else {
				Nmk = Nmk_v[i];
			}

			for (size_t nm = 0; nm < Nmk; ++nm) {
				auto X0 = X0X1.left + X0X1.volume() * G();
				auto Y0 = Y0Y1.left + Y0Y1.volume() * G();
				
				auto [e, l] = m_mes.fromXY(X0, Y0);
				Traj_t Lmax = LEf(-e);
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

				Traj_t tau_max = (Nmk_per_traj <= 10|| !extra_params.cut_tau) ? Traj_t(1) :
					find_tau_max(e,r0, r1, theta_max, mk, mi, deltaE_div_m_cm, Phi, TempR, th00);

				auto mk_factor =
					weight *bin_cut_factor * tau_max*Tin / ((Tin + Tout) * Nmk * Nmk_per_traj);


				if(tau_max > 0)
				{
					for (size_t nt = 0; nt < Nmk_per_traj; ++nt) {
						auto tau = G() * tau_max;
						auto [theta_undim, d_theta_undim] = th00(tau);
						auto theta = theta_undim * theta_max;
						auto d_theta = d_theta_undim * theta_max;

						auto cth2 = std::cos(theta / 2);
						auto sth2 = std::sin(theta / 2);
						Traj_t r = std::sqrt(u0 * cth2 * cth2 + u1 * sth2 * sth2);

#ifndef NDEBUG
						auto _C0 = Phi(std::sqrt(u0)) + e - Ltmp * Ltmp / u0;
						auto _C0_p = Phi(std::sqrt(u0 * (1 + 1e-4))) + e - Ltmp * Ltmp / (u0 * (1 + 1e-4));
						auto _C0_m = Phi(std::sqrt(u0 * (1 - 1e-4))) + e - Ltmp * Ltmp / (u0 * (1 - 1e-4));
						auto _C1 = Phi(std::sqrt(u1)) + e - Ltmp * Ltmp / u1;
						auto _C1_p = Phi(std::sqrt(u1 * (1 + 1e-4))) + e - Ltmp * Ltmp / (u1 * (1 + 1e-4));
						auto _C1_m = Phi(std::sqrt(u1 * (1 - 1e-4))) + e - Ltmp * Ltmp / (u1 * (1 - 1e-4));
#endif
						auto renorm_factor =
							d_theta_undim /
							(2 * Tin_Teheta * std::sqrt(S_func(r, r0, r1)));

						Traj_t v2 = downbound(Phi(r) + e, 0);
						auto v = std::sqrt(v2);
						Traj_t vt = upbound(Ltmp / r, v);
						Traj_t vr = ssqrt(v2 - vt * vt);

						vec3<Traj_t > V_in_nd(vt, 0, vr);
						vec3<Traj_t > V_in = V_in_nd * VescMin;
						Traj_t v_esc_nd = VescR(r);
						auto [dVout, factor] =
							DeltaVout1_Scatter(
								V_in, v* VescMin, G, ThermGen,
								mi, mk, mi_frac, mk_frac, m_cm, dm, deltaE_div_m_cm,
								dF,
								v_esc_nd * VescMin, VescMin,
								Nir(r), TempR(r)
							);
						vec3<Traj_t > dVout_nd = dVout * (1 / VescMin);
						auto VplusV1 = 2 * V_in_nd + dVout_nd;
						auto dVtau2 = dVout_nd[0] * VplusV1[0] +
							dVout_nd[1] * VplusV1[1];

						Traj_t deltaE = dVtau2 + dVout_nd[2] * VplusV1[2];
						//vec3<Traj_t > v_nd = Vout / VescMin;//OK
						Traj_t e_out = e + deltaE;
						//Traj_t e_out = ((V_in_nd + dVout_nd).squaredNorm() - v_esc_nd * v_esc_nd);

						auto final_factor = factor * mk_factor * renorm_factor;
						if (e_out < 0 && final_factor > ZeroValue / Nmk) {
							Traj_t deltaL2 = (r * r) * dVtau2;
							Traj_t L_nd2 = Ltmp * Ltmp + deltaL2;
							Traj_t L_nd = ssqrt(L_nd2);

							auto Lcutter = [](Traj_t Lmax) ->Traj_t {
								return Lmax > 0 ? 1 / Lmax : 0;
							};

							Traj_t Lm_inv = Lcutter(Lmax);
							Traj_t Lmax_out = LEf(-e_out);
							Traj_t Lm_inv1 = Lcutter(Lmax_out);

							Traj_t l_out = upbound(L_nd * Lm_inv1, 1);
							//Traj_t l_out1= (Lmax_out > 0 ? L_nd / Lmax_out : 0);
							//Traj_t deltal = l_out - l;
							//auto [x0, y0] = m_mes.toXY(e, l);
							auto [x1, y1] = m_mes.toXY(e_out, l_out);
							Traj_t dx = x1 - X0;
							Traj_t dy = y1 - Y0;


							Traj_t Xmin = X0X1.left + dx;
							Traj_t Xmax = X0X1.right + dx;
							Traj_t Ymin = Y0Y1.left + dy;
							Traj_t Ymax = Y0Y1.right + dy;
							bin_dedl_t<Traj_t> binout(
								grob::Rect<Traj_t>(Xmin, Xmax),
								grob::Rect<Traj_t>(Ymin, Ymax)
							);

							auto [emin, lmin] = m_mes.fromXY(Xmin, Ymin);
							auto [emax, lmax] = m_mes.fromXY(Xmax, Ymax);

							size_t i2_0 = grid.grid().pos(emin);
							size_t i2_1 = grid.grid().pos(emax) + 1;
							for (size_t i1 = i2_0; i1 < i2_1; ++i1) {
								auto gridl = grid.inner(i1);

								size_t j2_0 = gridl.pos(lmin);
								size_t j2_1 = gridl.pos(lmax) + 1;

								for (size_t j = j2_0; j < j2_1; ++j) {
									auto mbin = m_mes.toXY(grid[{i1, j}]);

									auto [binOut, isint] = grob::intersect(binout, mbin);
									auto VolOut = binOut.volume();
									if (isint) {
										auto vol_factor = VolOut * Lcutter(binout.volume());
										size_t k = grid.LinearIndex({ i1, j });
										Columns[i].coeffRef(k) += final_factor * vol_factor;
									}
								}
							}
						}
						else if (e_out >= 0) {
							EvapDistrib.Values[i] += final_factor;
						}
						else if (std::isnan(final_factor) || std::isnan(e_out) ){
							std::ostringstream error;
							error << "i = " << i << " nm = " << nm << " nt = " << nt << std::endl;
							throw std::runtime_error(error.str());
						}
					}
				}
			}
		}
		//std::vector<SpTriplet_t<Mat_vt>>
		//	AllTriplets(flatten(triplets));
		auto nonZeros = std::accumulate(
			Columns.begin(),
			Columns.end(),
			0ull,
			[](size_t acc, const auto& vec) {
				return acc + vec.nonZeros();
			}
		);
		auto m_triplets = triplet_view(Columns.begin(), Columns.end(), 
			MatrixIndexShift_Out, MatrixIndexShift_In);

		SpMatrix_t<Mat_vt> Ret(MatSize, MatSize);


		

		Ret.setFromTriplets(m_triplets.begin(), m_triplets.end());

		Ret.makeCompressed();

		return Ret;
	}



	template <
		typename Grid_vt,
		typename Mat_vt,
		typename Gen_t,
		typename ScatterFactor_t,
		typename NDensityEvent_t,
		typename ThermGenerator_t,
		typename HistoEvap_t,
		typename LEFunc_t,
		typename BinTrajPool_t,

		typename PhiFunc_t,
		typename SFunc_t,
		typename TempFunc_t,
		typename VescFunc_t,
		typename Nmk_Vec_t
	>
	NOINLINE SpMatrix_t<Mat_vt> ScatterImplShift1(
		std::type_identity<Grid_vt>,
		std::type_identity<Mat_vt>,
		size_t MatSize,
		Gen_t&& _G,
		ThermGenerator_t ThermGen,
		Gen_vt<Gen_t> mk, Gen_vt<Gen_t> dm,
		Gen_vt<Gen_t> mi,
		ScatterFactor_t const& dF,
		NDensityEvent_t const& Nir,
		size_t MatrixIndexShift_In,
		size_t MatrixIndexShift_Out,
		HistoEvap_t&& EvapDistrib,
		Gen_vt<Gen_t> ZeroValue,
		LEFunc_t const& LEf,
		Gen_vt<Gen_t> _F2_value,
		const BinTrajPool_t& TrajPoolVec,
		measure_dEpdlq<Grid_vt> m_mes,
		PhiFunc_t const& Phi,
		SFunc_t const& S_func,

		TempFunc_t const& TempR,
		VescFunc_t const& VescR,
		Gen_vt<Gen_t> VescMin,
		Nmk_Vec_t const& Nmk_v,
		size_t Nmk_per_traj,
		Gen_vt<Gen_t> weight,
		ScatterImplExtra extra_params = {}
	) {
		typedef Grid_vt Traj_t;


		const size_t N_in = EvapDistrib.Grid.size();
		auto&& grid = EvapDistrib.Grid;

		auto m_cm = mk * mi / (mk + mi);

		auto deltaE_div_m_cm = 2 * dm / m_cm;

		auto mi_frac = mi / (mk + mi);
		auto mk_frac = mk / (mk + mi);
		const auto _seed = _G.state;


		std::vector<Eigen::SparseVector<Mat_vt>> Columns(N_in);

#pragma	omp parallel for
		for (int i = 0; i < N_in; ++i) {

			Columns[i] = Eigen::SparseVector<Mat_vt>(N_in);

			auto G = _G;
			G.set_seed(_seed ^ i + 1);
			auto IJ = grid.FromLinear(i);

			auto el_bin = grid[IJ];
			auto [X0X1, Y0Y1] = m_mes.toXY(el_bin);

			//auto OutHisto = grob::make_histo_view(grid, OutIBuffer);
			const bin_traj_pool_t<Traj_t>& el_trajpool = TrajPoolVec[i];

			auto TinFunc = make_bin_traj_pool_Tin_func(
				el_trajpool, LEf, _F2_value);
			auto ToutFunc = make_bin_traj_pool_Tout_func(
				el_trajpool, LEf);

			//auto m_bin_el_gen = gen_EL(m_mes, el_bin, G, LEf);
			Traj_t bin_cut_factor = 1;
			if(extra_params.cut_el){
				Traj_t old_mes = m_mes(el_bin);
				el_bin = PrepairBin(el_bin, mk, mi, deltaE_div_m_cm, Phi, TempR, TinFunc, LEf);
				Traj_t new_mes = m_mes(el_bin);
				bin_cut_factor = new_mes/old_mes;
			}

			size_t Nmk = 0;
			if constexpr (std::is_same_v<Nmk_Vec_t, size_t>) {
				Nmk = Nmk_v;
			}
			else {
				Nmk = Nmk_v[i];
			}

			for (size_t corner_i = 0; corner_i < 4; ++corner_i) {
				auto X0 = corner_i < 2 ? X0X1.left : X0X1.right;
				auto Y0 = corner_i % 2 ? Y0Y1.left : Y0Y1.right;
				auto [e, l] = m_mes.fromXY(X0, Y0);
				Traj_t Lmax = LEf(-e);
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

				Traj_t tau_max = !extra_params.cut_tau ? Traj_t(1) : 
					find_tau_max(e,r0, r1, theta_max, mk, mi, deltaE_div_m_cm, Phi, TempR, th00);

				auto mk_factor =
					tau_max* bin_cut_factor *weight * Tin / ((Tin + Tout) * Nmk) / 4;
				

				if (tau_max > 0) {
					for (size_t nt = 0; nt < Nmk; ++nt) {
						auto tau = G()*tau_max;
						auto [theta_undim, d_theta_undim] = th00(tau);
						auto theta = theta_undim * theta_max;
						
						auto cth2 = std::cos(theta / 2);
						auto sth2 = std::sin(theta / 2);
						Traj_t r = std::sqrt(u0 * cth2 * cth2 + u1 * sth2 * sth2);
						auto renorm_factor = 1;

						Traj_t v2 = downbound(Phi(r) + e, 0);
						auto v = std::sqrt(v2);
						Traj_t vt = upbound(Ltmp / r, v);
						Traj_t vr = ssqrt(v2 - vt * vt);

						vec3<Traj_t > V_in_nd(vt, 0, vr);
						vec3<Traj_t > V_in = V_in_nd * VescMin;
						Traj_t v_esc_nd = VescR(r);
						auto [dVout, factor] =
							DeltaVout1_Scatter(
								V_in, v* VescMin, G, ThermGen,
								mi, mk, mi_frac, mk_frac, m_cm, dm, deltaE_div_m_cm,
								dF,
								v_esc_nd * VescMin, VescMin,
								Nir(r), TempR(r)
							);
						vec3<Traj_t > dVout_nd = dVout * (1 / VescMin);
						auto VplusV1 = 2 * V_in_nd + dVout_nd;
						auto dVtau2 = dVout_nd[0] * VplusV1[0] +
							dVout_nd[1] * VplusV1[1];

						Traj_t deltaE = dVtau2 + dVout_nd[2] * VplusV1[2];
						//vec3<Traj_t > v_nd = Vout / VescMin;//OK
						Traj_t e_out = e + deltaE;
						//Traj_t e_out = ((V_in_nd + dVout_nd).squaredNorm() - v_esc_nd * v_esc_nd);

						auto final_factor = factor * mk_factor * renorm_factor;
						if (e_out < 0 && final_factor > ZeroValue / Nmk) {
							Traj_t deltaL2 = (r * r) * dVtau2;
							Traj_t L_nd2 = Ltmp * Ltmp + deltaL2;
							Traj_t L_nd = ssqrt(L_nd2);

							auto Lcutter = [](Traj_t Lmax) ->Traj_t {
								return Lmax > 0 ? 1 / Lmax : 0;
							};

							Traj_t Lm_inv = Lcutter(Lmax);
							Traj_t Lmax_out = LEf(-e_out);
							Traj_t Lm_inv1 = Lcutter(Lmax_out);

							Traj_t l_out = upbound(L_nd * Lm_inv1, 1);
							//Traj_t l_out1= (Lmax_out > 0 ? L_nd / Lmax_out : 0);
							//Traj_t deltal = l_out - l;
							//auto [x0, y0] = m_mes.toXY(e, l);
							auto [x1, y1] = m_mes.toXY(e_out, l_out);
							Traj_t dx = x1 - X0;
							Traj_t dy = y1 - Y0;


							Traj_t Xmin = X0X1.left + dx;
							Traj_t Xmax = X0X1.right + dx;
							Traj_t Ymin = Y0Y1.left + dy;
							Traj_t Ymax = Y0Y1.right + dy;
							bin_dedl_t<Traj_t> binout(
								grob::Rect<Traj_t>(Xmin, Xmax),
								grob::Rect<Traj_t>(Ymin, Ymax)
							);

							auto [emin, lmin] = m_mes.fromXY(Xmin, Ymin);
							auto [emax, lmax] = m_mes.fromXY(Xmax, Ymax);

							size_t i2_0 = grid.grid().pos(emin);
							size_t i2_1 = grid.grid().pos(emax) + 1;
							for (size_t i1 = i2_0; i1 < i2_1; ++i1) {
								auto gridl = grid.inner(i1);

								size_t j2_0 = gridl.pos(lmin);
								size_t j2_1 = gridl.pos(lmax) + 1;

								for (size_t j = j2_0; j < j2_1; ++j) {
									auto mbin = m_mes.toXY(grid[{i1, j}]);

									auto [binOut, isint] = grob::intersect(binout, mbin);
									auto VolOut = binOut.volume();
									if (isint) {
										auto vol_factor = VolOut * Lcutter(binout.volume());
										size_t k = grid.LinearIndex({ i1, j });
										Columns[i].coeffRef(k) += final_factor * vol_factor;
									}
								}
							}
						}
						else if (e_out >= 0) {
							EvapDistrib.Values[i] += final_factor;
						}

					}
				}
				
			}
		}
		//std::vector<SpTriplet_t<Mat_vt>>
		//	AllTriplets(flatten(triplets));
		auto nonZeros = std::accumulate(
			Columns.begin(),
			Columns.end(),
			0ull,
			[](size_t acc, const auto& vec) {
				return acc + vec.nonZeros();
			}
		);
		auto m_triplets = triplet_view(Columns.begin(), Columns.end(),
			MatrixIndexShift_Out, MatrixIndexShift_In);

		SpMatrix_t<Mat_vt> Ret(MatSize, MatSize);




		Ret.setFromTriplets(m_triplets.begin(), m_triplets.end());

		Ret.makeCompressed();

		return Ret;
	}
}; //namespace evdm 