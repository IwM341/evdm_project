#pragma once
#include "../measure.hpp"
#include <grob/grid.hpp>
#include <grob/grid_objects.hpp>
#include <numbers>
#include "../utils/mc.hpp"
#include "../utils/progress_bar.hpp"

namespace evdm {
	template <typename Vec_t>
	concept VectorC = requires(Vec_t t, size_t i) {
		{t.size()}->std::convertible_to<size_t>;
		{t.operator [](i)};
	};

	template <typename Vec_t>
	concept MatrixC = requires (Vec_t t,size_t i) {
		{t(i, i)};
		{t.col(i)}->VectorC;
		{t.row(i)}->VectorC;
	};

	template <
		MatrixC pre_mat_t,
		typename GridEL_t,
		typename EL_Functype,
		typename Phi_Func_t,
		typename TrajPoolArray_t,
		typename _F2_t,
		typename Gen_t
	>
	void AnnImpl(
		pre_mat_t &A0,
		pre_mat_t &Av,
		GridEL_t const & grid,
		Gen_vt<Gen_t> Rmin, Gen_vt<Gen_t> Rmax,
		TrajPoolArray_t Trajects,
		EL_Functype const& LEf,
		Phi_Func_t const& Phi,
		_F2_t _F2,
		Gen_t _G,
		size_t Nmk_r,
		progress_omp_function<>& m_progress_func
	) {
		const size_t N_in = grid.size();
		progress_omp_bar<> m_bar(
			m_progress_func, N_in, std::max((int)(N_in / 1000), 1)
		);
		typedef std::decay_t<decltype(grid.grid()[0].left)> T;
		const auto _seed = _G.state;
		#pragma	omp parallel for
		for (int I_in = 0; I_in < N_in; ++I_in) {
			auto G = _G;
			G.set_seed(_seed ^ (size_t)(I_in)+1);

			auto IJ_in = grid.FromLinear(I_in);


			bin_traj_pool_t<T> const& pool_in = Trajects[I_in];
			size_t p_i_e = pool_in.Grid.grid().size();
			size_t p_i_l = pool_in.Grid.inner(0).size();
			using TI_ref = traj_info<T>const&;
			auto r_min_in = std::max(
				((TI_ref)pool_in[{0, 0}]).rmin,
				(T)Rmin
			);
			auto r_max_in = std::min(
				((TI_ref)pool_in[{p_i_e - 1, 0}]).rmax,
				(T)Rmax
			);

			bin_dedl_t<T> m_bin_in = grid[IJ_in];
			if (r_min_in <= r_max_in) {
				for (size_t i_out = 0; i_out < grid.grid().size(); ++i_out)
				{
					for (size_t j_out = 0; j_out < grid.inner(i_out).size(); ++j_out) {
						auto push_seed = G.state;
						auto I_out = grid.LinearIndex({ i_out,j_out });
						bin_traj_pool_t<T> const& pool_out = Trajects[I_out];
						size_t p_o_e = pool_out.Grid.grid().size();
						size_t p_o_l = pool_out.Grid.inner(0).size();
						auto r_min_out = ((TI_ref)pool_out[{0, 0}]).rmin;
						auto r_max_out = ((TI_ref)pool_out[{p_o_e - 1, 0}]).rmax;
						auto r_min = std::max(r_min_out, r_min_in);
						auto r_max = std::min(r_max_out, r_max_in);
						if (r_min <= r_max) {
							auto m_bin_out = grid[{i_out, j_out}];

							auto dens_in = 1 / dEdL2(
								m_bin_in, LEf(-std::get<0>(m_bin_in).left),
								LEf(-std::get<0>(m_bin_in).right)
							);
							auto dens_out = 1 / dEdL2(
								m_bin_out, LEf(-std::get<0>(m_bin_out).left),
								LEf(-std::get<0>(m_bin_out).right)
							);

							auto TinFunc_in = make_bin_traj_pool_Tin_func(
								pool_in, LEf, _F2);
							auto ToutFunc_in = make_bin_traj_pool_Tout_func(
								pool_in, LEf);
							auto PeriodFunc_in = [&](auto e, auto l, auto Lmax) {
								return TinFunc_in(e, l, Lmax) + ToutFunc_in(e, l, Lmax);
							};

							auto TinFunc_out = make_bin_traj_pool_Tin_func(
								pool_out, LEf, _F2);
							auto ToutFunc_out = make_bin_traj_pool_Tout_func(
								pool_out, LEf);
							auto PeriodFunc_out = [&](auto e, auto l, auto Lmax) {
								return TinFunc_out(e, l, Lmax) + ToutFunc_out(e, l, Lmax);
							};

							auto gen_esl = [](
								bin_dedl_t<T> const& bin, auto r, auto phi_r,
								auto&& LEf, auto&& G
								) {
									T e1 = std::get<0>(bin).right;
									T e0 = std::clamp(-(T)phi_r, std::get<0>(bin).left, e1);
									T e = e0 + (e1 - e0) * G();
									T Lmax_e = LEf(-e);
									T v2_max = downbound(phi_r + e, 0);
									T L2m = v2_max * r * r;
									T L0 = Lmax_e * std::get<1>(bin).left;
									auto l1 = std::get<1>(bin).right;
									T L1 = Lmax_e * l1;
									auto sqr0 = ssqrt(L2m - L1 * L1);
									auto sqr1 = ssqrt(L2m - L0 * L0);
									auto weight_s = (sqr1 - sqr0);
									auto sqr = sqr0 + G() * weight_s;
									auto weight = (e1 - e0) * weight_s;
									auto l = upbound(std::sqrt(L2m - sqr * sqr) / Lmax_e, l1);
									return std::make_tuple(
										e, l, v2_max,
										Lmax_e, weight
									);
							};

							auto& A0_ij = A0(I_in, I_out);
							auto& Av_ij = Av(I_in, I_out);
							auto factor_r = r_max - r_min;
							auto dens_weight = factor_r * dens_in * dens_out * 4 / (3 * Nmk_r);

							for (size_t _mki = 0; _mki < Nmk_r; ++_mki) {
								auto r = r_min + factor_r * G();

								auto phi_r = (r < 1 ? Phi(r) : 1 / r);


								auto [e1, l1, v21, Lm1, W1] =
									gen_esl(m_bin_in, r, phi_r, LEf, G);

								auto [e2, l2, v22, Lm2, W2] =
									gen_esl(m_bin_out, r, phi_r, LEf, G);

								auto Tl1 = PeriodFunc_in(e1, l1, Lm1);
								auto Tl2 = PeriodFunc_out(e2, l2, Lm2);

								auto weight_s = dens_weight * W1 * W2 / (Tl1 * Tl2);
								if (std::isnan(weight_s)) {
									std::string message =
										"nan at annihilation, seed = " +
										std::to_string(push_seed) +
										", (I_in,I_out) = (" +
										std::to_string(I_in) + ", " +
										std::to_string(I_out) +
										"), mki = " + std::to_string(_mki);

									throw std::runtime_error(message);
								}
								A0_ij += weight_s;
								Av_ij += weight_s * (v21 + v22) / 2;
							}
						}
					}
				}
			}

			m_bar.next();
		}
	}

};//END NAMESPACE