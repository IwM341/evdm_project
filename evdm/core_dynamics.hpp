#ifndef CORE_DYNAMICS_HPP
#define CORE_DYNAMICS_HPP
#include "core.hpp"
#include "traj_pool.hpp"
#include "trajectory.hpp"


namespace evdm {
    
    /// <summary>
    /// Calculates Capture
    /// </summary>
    /// <param name="mDistrib">Vector of LE distribution of particles</param>
    /// <param name="ptype">index of paricle type</param>
    /// <param name="se">Scatter Event (e.g. nucleus concentration. + form factor)</param>
    /// <param name="G">random generator</param>
    /// <param name="mk">WIMP mass</param>
    /// <param name="dm">out wimp mass - in wimp mass</param>
    /// <param name="mi">nuclei mass</param>
    /// <param name="mU0">speed of body relative to halo</param>
    /// <param name="Vdisp">dispersion of DM speed in halo</param>
    /// <param name="Nmk">number of monte-carlo samples</param>
    /// <param name="weight">additional weight to scale result, default - 1</param>
    template <typename DistrT,typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type,typename Gen_t>
    std::pair<double, double> Capture(Distribution<DistrT, Body_vt, GridEL_vt, grid_type> & mDistrib,
        size_t ptype, ScatterEvent const & se, Gen_t && G,
        Gen_vt< Gen_t> mk, Gen_vt< Gen_t> dm, Gen_vt< Gen_t> mi,
        Gen_vt< Gen_t> mU0, Gen_vt< Gen_t> Vdisp,float pow3_r,
        size_t  Nmk, Gen_vt< Gen_t> weight = 1 ) 
    {
        
        auto const& Temp = mDistrib.body().Temp;
        auto mHisto_full = mDistrib.as_histo();
        auto mHisto = mHisto_full.inner_slice(ptype);
        auto LEf = mDistrib.Grid.LE();
        return CaptureImpl(G, se, Temp, mHisto, LEf, mk, dm, mi, pow3_r, mU0, Vdisp,
            mDistrib.body().VescFunc,mDistrib.body().Vesc, Nmk, weight);
    }

    template <typename DistrT, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type, 
        typename Gen_t,typename BinMes_t>
    void Scatter(
        GridMatrix<DistrT, Body_vt, GridEL_vt, grid_type>& mMatrix,
        Distribution<DistrT, Body_vt, GridEL_vt, grid_type>& mEvap,
        bool count_evap,
        size_t ptype_in, size_t ptype_out, ScatterEvent const& se,
        Gen_t&& G, BinMes_t m_measure,
        Gen_vt< Gen_t> mk, Gen_vt< Gen_t> dm, Gen_vt< Gen_t> mi,
        size_t  Nmk,size_t Nmk_per_traj, 
        Gen_vt< Gen_t> weight,
        progress_omp_function<>& m_progress_func)
    {
        auto const & m_body = mMatrix.body();
        auto const& Temp = mMatrix.body().Temp;
        auto const& Phi = mMatrix.body().Phi;
        auto mHisto_full = mMatrix.as_histo_bank(ptype_in, ptype_out);
        auto EvapVector = mEvap.as_histo();
        auto LEf = mMatrix.Grid.LE();
        auto _F2 = mMatrix.body()._DD_F_C(1);
        auto Sfunc = [&](auto r, auto rmin, auto rmax) {
            return m_body.S(r,rmin,rmax);
        };
        ScatterImpl(
            G,
            mk, dm, mi, se,
            mHisto_full, EvapVector, count_evap,
            LEf, _F2, mMatrix.Grid.TrajPools(),
            m_measure, Phi, Sfunc,Temp,
            mEvap.body().VescFunc, mEvap.body().Vesc,
            Nmk, Nmk_per_traj, weight, m_progress_func
        );
    }

    /*
    template <typename Body_vt, typename GridEL_vt, GridEL_type grid_type>
    struct TrajPool {
        using TrajT = GridEL_vt;

        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt, GEL_vt, grid_type>;
        GridEL_t Grid;

        using bin_pull_t = bin_traj_pool_t<TrajT>;
        std::shared_ptr<std::vector<bin_pull_t>> Values;


        static constexpr Body_vt get_body_vtype() {
            static_assert("TrajPool::get_body_vtype() is not callable");
        };
        static constexpr TrajT get_vtype()
        {
            static_assert("TrajPool::get_vtype() is not callable");
        }
        static constexpr GridEL_vt get_grid_vtype() {
            static_assert("TrajPool::get_grid_vtype() is not callable");
        }
        static constexpr GridEL_type get_grid_type() {
            return grid_type;
        }
        static constexpr bin_traj_pool_t<TrajT> get_pull_t() {
            static_assert("TrajPool::bin_pull_t() is not callable");
        }

        inline auto const& grid()const {
            return *(Grid.Grid);
        }
        inline auto const& values() const {
            return *Values;
        }
        inline auto& values() {
            return *Values;
        }
        TrajPool(
            GridEL_t const& Grid,
            size_t traj_bins = 64, TrajT _err = 0.05,
            size_t n_e_max = 2, size_t n_l_max = 2) :
            Grid(Grid), Values(std::make_shared<std::vector<bin_pull_t>>() )
        {
            n_e_max = std::max(n_e_max, (size_t)2);
            n_l_max = std::max(n_l_max, (size_t)2);
            bool is_static = (n_e_max <= 2) && (n_l_max <= 2);
            Values->reserve(Grid.Grid->inner(0).size());
            auto LE = Grid.LE();
            auto const& B = *Grid.body;
            for (auto el_bin : Grid.Grid->inner(0)) {
                Values->push_back(
                    GetTrajPool<TrajT>(
                        is_static,
                        el_bin,B,LE,
                        _err, traj_bins, 
                        n_e_max,n_l_max)
                );
            }
        }
        auto as_grid_object()const {
            return grob::make_grid_object_ref(grid(), values());
        }
    };
    */
}

#endif//CORE_DYNAMICS_HPP