#pragma once

#include "../dynamics/dynamics_scatter.hpp"

#include "../core/core_matrix.hpp"
#include "../core/core_distrib.hpp"

namespace evdm {
    template <typename DistrT, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type,
        typename Gen_t, typename BinMes_t>
    void Scatter(
        GridMatrix<DistrT, Body_vt, GridEL_vt, grid_type>& mMatrix,
        Distribution<DistrT, Body_vt, GridEL_vt, grid_type>& mEvap,
        bool count_evap,
        size_t ptype_in, size_t ptype_out, ScatterEvent const& se,
        Gen_t&& G, BinMes_t m_measure,
        Gen_vt< Gen_t> mk, Gen_vt< Gen_t> dm, Gen_vt< Gen_t> mi,
        size_t  Nmk, size_t Nmk_per_traj,
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
        ScatterImpl(
            G,
            mk, dm, mi, se,
            mHisto_full, EvapVector, count_evap,
            LEf, _F2, mMatrix.Grid.TrajPools(),
            m_measure, Phi, Sfunc, Temp,
            mEvap.body().VescFunc, mEvap.body().Vesc,
            Nmk, Nmk_per_traj, weight, m_progress_func
        );
    }
};//namespace evdm 