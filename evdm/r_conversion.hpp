#pragma once
#include "measure.hpp"
#include <grob/grid_objects.hpp>
#include "utils/progress_bar.hpp"
#include "utils/math_external.hpp"
#include <numbers>
#include <utility>
#include <numbers>
#include "r_conv_util.hpp"
/*
* Conversion from EL density to r density
* Scatter integral with sigma
*/

namespace evdm{

    

    /// @brief return grid function dN/d^3r(r)
    /// @param H EL distribution fistogramm
    /// @param Trajects vector of trajectory pool, corresponding to EL grid
    /// ptypes vector of considered particle types
    /// @param LEf function of L(E)
    /// @param mu_el measur mu(dEdL)
    /// @param G PRNG
    /// @param r_grid grid of r
    /// @param Nmk_per_bin number of mk samples integrating through LE bin
    /// @return GridFunction r->dN/d3r(r)
    template <
        typename HistoType,
        typename TrajPoolArray_t,
        typename EL_Functype,
        typename Phi_Func_t,
        typename _F2_t,
        typename Gen_t,
        typename RGrid_t,
        typename N_mk_t
    >
    auto convert_r_density(
        HistoType const &H_vec,
        TrajPoolArray_t Trajects,
        std::vector<size_t> const &ptypes,
        EL_Functype const & LEf,
        Phi_Func_t const & Phi,
        _F2_t _F2,
        Gen_t _G,
        RGrid_t r_grid, N_mk_t const & Nmk_per_bin_distrib,
        progress_omp_function<> m_prog_bar_func = progress_omp_function<>())
    {
        const size_t Nrg = r_grid.size();
        auto RDens = grob::make_function(std::move(r_grid),
            std::vector<double>(Nrg,0)/*Eigen::VectorXd(Nrg)*/);

        auto& ELGrid = H_vec.Grid.inner(0);
        typedef std::decay_t<decltype(ELGrid.grid()[0].left)> T;
        


        progress_omp_bar<> m_prog_bar(m_prog_bar_func, Nrg,
            std::max((int)Nrg / 100, (int)1));
        auto seed = _G.state;

        #pragma omp parallel for
        for (int ird = 0; ird < Nrg; ++ird) {
            auto G = _G;
            G.set_seed(seed ^ (ird + 1));
            auto& dense_accum = RDens[ird];
            T r = RDens.Grid[ird];
            T phi_r = (r < 1 ? Phi(r) : 1/r);

            const size_t i_E_min = ELGrid.grid().pos(-phi_r);
            for (size_t i = i_E_min; i < ELGrid.grid().size(); ++i) {
                auto max_L_e1 = LEf(-ELGrid.grid()[i].left);
                auto max_L_bound = ssqrt(phi_r + ELGrid.grid()[i].right) * r;
                const size_t j_L_max= (ELGrid.inner(i).pos(
                    max_L_e1 > 0 ? max_L_bound/ max_L_e1 : 1));

                for (size_t j = 0; j <= j_L_max; ++j) {
                    auto IJ_index = grob::make_MI(i, j);
                    auto IJ_LinIndex = ELGrid.LinearIndex(IJ_index);
                    auto m_bin = ELGrid[IJ_index];
                    /*FORMULA DEBUG,REMOVE*/
                    //std::get<0>(m_bin) = { -1.21,-1.2 };
                    //std::get<1>(m_bin) = { 0.7,0.701 };
                    /*FORMULA DEBUG,REMOVE*/
                    auto c_bin = _bin_cut(LEf, m_bin,phi_r, r);
                    if (std::get<0>(c_bin).left != std::get<0>(c_bin).right) {
                        auto Vol = dEdL2(
                            m_bin, LEf(-std::get<0>(m_bin).left),
                            LEf(-std::get<0>(m_bin).right)
                        );
                        auto N_paticles = [&]() {
                            std::decay_t<decltype(H_vec.Values[0])> sum = 0;
                            for (auto p : ptypes) {
                                sum += H_vec[{p, IJ_index}];
                            }
                            return sum;
                        }();
                        /*FORMULA DEBUG,REMOVE*/
                        //N_paticles = 1;
                        /*FORMULA DEBUG,REMOVE*/
                        auto Bin_Dens = N_paticles / Vol;



                        auto TinFunc = make_bin_traj_pool_Tin_func(
                            Trajects[IJ_LinIndex], LEf, _F2);
                        auto ToutFunc = make_bin_traj_pool_Tout_func(
                            Trajects[IJ_LinIndex], LEf);

                        auto PeriodFunc = [&](auto e, auto l, auto Lmax) {
                            return TinFunc(e, l, Lmax) + ToutFunc(e, l, Lmax);
                        };

                        /*checks {
                            auto [e0, e1] = std::get<0>(c_bin);
                            auto [l0, l1] = std::get<1>(c_bin);
                            if (e0 < -phi_r) {
                                throw std::runtime_error("error e0 < -phi_r");
                            }
                            if (e1 < -phi_r) {
                                throw std::runtime_error("error e1 < -phi_r");
                            }
                            if (e0 > e1) {
                                throw std::runtime_error("error e1 < e0");
                            }
                        }*/

                        double RDensValue = 0;
                        size_t Nmk_per_bin = 0;
                        if constexpr (
                            std::is_arithmetic_v<N_mk_t>
                        ) {
                            Nmk_per_bin = Nmk_per_bin_distrib;
                        } else {
                            Nmk_per_bin = static_cast<size_t>(
                                Nmk_per_bin_distrib[IJ_LinIndex]
                            );
                        }

                        double factor = Bin_Dens*2 / (Nmk_per_bin * 3);
                        auto bin_gen = mc_d3v(r, phi_r, c_bin, G, LEf);

                        for (size_t n = 0; n < Nmk_per_bin; ++n) {
                            auto [EL, dense] = bin_gen();
                            auto [tE, tl, Lmax] = EL;
                            auto Tl = PeriodFunc(tE, tl, Lmax);
                            RDensValue += factor * dense / Tl;
                        }
                        dense_accum += RDensValue;
                        /*FORMULA DEBUG,REMOVE*/
                        //goto r_next_label;
                        /*FORMULA DEBUG,REMOVE*/
                    }
                }
            }
            m_prog_bar.next();
            /*FORMULA DEBUG,REMOVE*/
        //r_next_label:;
            /*FORMULA DEBUG,REMOVE*/
        }
        return RDens;
    }


    
};
