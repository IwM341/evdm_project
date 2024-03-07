#ifndef R_CONVERSION_HPP
#define R_CONVERSION_HPP
#include "measure.hpp"
#include <grob/grid_objects.hpp>
#include "utils/math_external.hpp"
#include <numbers>
#include <utility>
#include <numbers>
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
        typename RGrid_t
    >
    auto convert_r_density(
        HistoType const &H_vec,
        TrajPoolArray_t Trajects,
        std::vector<size_t> const &ptypes,
        EL_Functype const & LEf,
        Phi_Func_t const & Phi,
        _F2_t _F2,
        Gen_t G,
        RGrid_t r_grid,size_t Nmk_per_bin = 10)
    {
        const size_t Nrg = r_grid.size();
        auto RDens = grob::make_function(std::move(r_grid),
            std::vector<double>(Nrg,0)/*Eigen::VectorXd(Nrg)*/);

        auto& ELGrid = H_vec.Grid.inner(0);
        typedef std::decay_t<decltype(ELGrid.grid()[0].left)> T;
        
        auto bin_cut = [&LEf](auto const& m_bin, auto phi, auto r) {
            auto bin1 = m_bin;
            auto& [e0, e1] = std::get<0>(bin1);
            const auto& l0 = std::get<1>(bin1).left;
            auto line_parabole_solver = [](T x0, T x1, T y0, T y1,T x_p_0,T C_p)->
                std::pair<T,T>{
                auto TanL = (y1 - y0) / (x1 - x0);
                auto T2 = TanL * TanL;
                auto _sqrt_s = std::sqrt(C_p*(C_p-4* TanL * y0 - 4*T2* x_p_0+4*T2*x0));
                auto _b_s = 2 * T2 * x0 - 2 * TanL * y0 + C_p;
                auto delta_sqrt_div_2t2 =  
                    2*T2*x0*x0 - 4*TanL*x0*y0 + 2*C_p*x_p_0+2*y0*y0;
                if (TanL == 0) {
                    return { x_p_0,std::numeric_limits<T>::infinity() };
                }
                if (!(_sqrt_s >= 0)) {
                    return {std::numeric_limits<T>::infinity(),
                            std::numeric_limits<T>::infinity() };
                }
                else {
                    if (_b_s > 0) {
                        return { delta_sqrt_div_2t2 / (_sqrt_s + _b_s),
                            (_b_s + _sqrt_s) / (2 * T2) };
                    }
                    else if(_b_s < 0) {
                        return { (_b_s - _sqrt_s) / (2 * T2),
                            delta_sqrt_div_2t2 / (_b_s - _sqrt_s) };
                    }
                    else {
                        return { -_sqrt_s / (2 * T2) ,_sqrt_s / (2 * T2) };
                        return { -_sqrt_s / (2 * T2) ,_sqrt_s / (2 * T2) };
                    }
                }
            };
            auto [em0, em1] = 
                line_parabole_solver(e0,e1, LEf(-e0)*l0, LEf(-e1) * l0,-phi,r*r);
            em0 = std::max(em0, -(T)phi);
            e0 = std::max((em0 < e1 ? em0 : e1), e0);
            e1 = std::min((em1 > e0 ? em1 : e0), e1);
            return bin1;
        };
        
        //#pragma omp omp_in_parallel for private(G)
        for (size_t ird = 0; ird < Nrg; ++ird) {
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
                    auto c_bin = bin_cut(m_bin, phi_r, r);
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
                            return TinFunc(e, l, Lmax) + ToutFunc(e, l);
                        };

                        /*checks*/ {
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
                        }

                        double RDensValue = 0;
                        double factor = Bin_Dens / (Nmk_per_bin * std::numbers::pi * 2);
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
            /*FORMULA DEBUG,REMOVE*/
        //r_next_label:;
            /*FORMULA DEBUG,REMOVE*/
        }
        return RDens;
    }


    /// @brief calculates annihilation matrix
    /// @param Grid EL Grid
    /// @param m_matrix lambda matrix accessor m_matrix(i,j)
    /// @param LE L(E) func
    /// @param d3v measure, d3v(dEdL)
    /// @param mu_el measur mu(dEdL)
    /// @param G random number generator
    /// @param VSigma function, equal to |v1-v2|*sigma(|v1-v2|), where 
    /// sigma --- annihilation cross section, |v1-v2| --- velosity difference  
    /// @return 
    template <typename GridType,
                typename MatrixAcsessor_t,
                typename LE_Functype,
                typename d3v_measure_t,
                typename el_measure_t,
                typename Gen_t,
                typename SigmaFuc_t>
    auto FillAnnMatrix(GridType const & Grid,
                        MatrixAcsessor_t && m_matrix,
                        LE_Functype const & LE,
                        d3v_measure_t d3v,
                        el_measure_t mu_el,
                        Gen_t && G,
                        SigmaFuc_t VSigma)
    {
        //TODO
    }
    
};
#endif//R_CONVERSION_HPP