#pragma once

#include "../dynamics/dynamics_scatter.hpp"

#include "../core/core_matrix.hpp"
#include "../core/core_distrib.hpp"

namespace evdm {
    template <typename DistrT, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type,
        typename Gen_t,typename ThermGaussGenerator_t, 
        typename BinMes_t,
        typename Nmk_vec_t>
    void Scatter(
        GridMatrix<DistrT, Body_vt, GridEL_vt, grid_type>& mMatrix,
        Distribution<DistrT, Body_vt, GridEL_vt, grid_type>& mEvap,
        bool count_evap,
        size_t ptype_in, size_t ptype_out, ScatterEvent const& se,
        Gen_t&& G, ThermGaussGenerator_t ThermGen, BinMes_t m_measure,
        Gen_vt< Gen_t> mk, Gen_vt< Gen_t> dm, Gen_vt< Gen_t> mi,
        Nmk_vec_t  const &Nmk_v, size_t Nmk_per_traj,
        Gen_vt< Gen_t> weight,
        progress_omp_function<>& m_progress_func)
    {
        auto const& m_body = mMatrix.body();
        auto const& Temp = mMatrix.body().Temp;
        auto const& Phi = mMatrix.body().Phi;
        auto mHisto_full = mMatrix.as_histo_bank(ptype_in, ptype_out);
        auto EvapVector = mEvap.as_histo();
        auto LEf = mMatrix.Grid.LE();
        auto _F2 = mMatrix.body()._DD_F_C(1);
        auto Sfunc = [&](auto r, auto rmin, auto rmax) {
            return m_body.S(r, rmin, rmax);
        };
        std::visit(
            [&](auto const& dF) {
                ScatterImpl(
                    std::type_identity< GridEL_vt>{}, G, ThermGen,
                    mk, dm, mi, dF,se.n_e,
                    mHisto_full, EvapVector, count_evap,
                    LEf, _F2, mMatrix.Grid.TrajPools(),
                    m_measure, Phi, Sfunc, Temp,
                    mEvap.body().VescFunc, mEvap.body().Vesc,
                    Nmk_v, Nmk_per_traj, weight, m_progress_func
                );
            },
            se.sf
        );
    }
};//namespace evdm 