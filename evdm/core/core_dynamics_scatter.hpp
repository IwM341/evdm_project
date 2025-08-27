#pragma once

#include "../dynamics/dynamics_scatter.hpp"

#include "../core/core_matrix.hpp"
#include "../core/core_distrib.hpp"

namespace evdm {
    template <typename DistrT, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type,
        typename Gen_t,typename ThermGaussGenerator_t, 
        typename Nmk_vec_t>
    void Scatter(
        GridMatrix<DistrT, Body_vt, GridEL_vt, grid_type>& mMatrix,
        Distribution<DistrT, Body_vt, GridEL_vt, grid_type>& mEvap,
        bool count_evap,
        size_t ptype_in, size_t ptype_out, ScatterEvent const& se,
        Gen_t&& G, ThermGaussGenerator_t ThermGen, measure_dEpdlq<Gen_vt< Gen_t>> m_measure,
        Gen_vt< Gen_t> mk, Gen_vt< Gen_t> dm, Gen_vt< Gen_t> mi,
        Nmk_vec_t  const &Nmk_v, size_t Nmk_per_traj,
        Gen_vt< Gen_t> weight,
        progress_omp_function<>& m_progress_func)
    {
        auto const& m_body = mMatrix.body();
        auto const& Temp = mMatrix.body().Temp;
        auto const& Phi = mMatrix.body().Phi;
        
        //auto mHisto_full = mMatrix.as_histo_bank(ptype_in, ptype_out);
        
        size_t GridSize = mMatrix.grid().size1();
        size_t ptypes = mMatrix.grid().ptypes();

        auto EvapVector = mEvap.as_histo();
        auto EvapSlice = EvapVector.inner_slice(ptype_in);
        auto LEf = mMatrix.Grid.LE();
        auto _F2 = mMatrix.body()._DD_F_C(1);
        auto Sfunc = [&](auto r, auto rmin, auto rmax) {
            return m_body.S(r, rmin, rmax);
        };
        size_t InShift = ptype_in* GridSize;
        size_t OutShift = ptype_out * GridSize;
        typedef decltype(mMatrix.get_vtype()) Mat_t;
        SpMatrix_t< Mat_t> MmatNew = std::visit(
            [&](auto const& dF) -> SpMatrix_t< Mat_t> {
                return ScatterImpl(
                    std::type_identity< GridEL_vt>{},
                    std::type_identity< Mat_t>{}, GridSize* ptypes,G, ThermGen,
                    mk, dm, mi, dF,se.n_e,
                    InShift,OutShift, EvapSlice, count_evap,
                    LEf, _F2, mMatrix.Grid.TrajPools(),
                    m_measure, Phi, Sfunc, Temp,
                    mEvap.body().VescFunc, mEvap.body().Vesc,
                    Nmk_v, Nmk_per_traj, weight, m_progress_func
                );
            },
            se.sf
        );
        mMatrix.values() += MmatNew;
    }
};//namespace evdm 