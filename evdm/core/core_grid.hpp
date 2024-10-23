#pragma once
#include <Eigen/Dense>
#include "../grid_variants.hpp"
#include "../traj_pool.hpp"
#include "core_body.hpp"
#include <memory>
namespace evdm {

    struct TrajPoolInitParams_t {
        bool is_static;
        float T_error;
        size_t traj_bins;
        size_t n_e_max;
        size_t n_l_max;
        TrajPoolInitParams_t(
            bool is_static = true,
            float T_error = 0.05,
            size_t traj_bins = 100,
            size_t n_e_max = 2,
            size_t n_l_max = 2
        ) :
            is_static(is_static),
            T_error(T_error),
            traj_bins(traj_bins),
            n_e_max(n_e_max),
            n_l_max(n_l_max) {}
    };




    template <typename BodyModel_vt, typename GridEL_vt, GridEL_type _grid_type>
    struct EL_Grid {
        using BM_vt = BodyModel_vt;
        using GEL_vt = GridEL_vt;
        using GridEL_t = GridEL<GridEL_vt, _grid_type>;
        using Body_t = Body<BM_vt>;
        using LE_func_t = typename Body<BM_vt>::LE_func_t;
        using TrajPool_t = std::vector<bin_traj_pool_t<GridEL_vt>>;

        static constexpr BodyModel_vt get_bm_vtype() {
            static_assert("EL_Grid::get_bm_vtype() is not callable");
        }
        static constexpr GEL_vt get_grid_vtype() {
            static_assert("EL_Grid::get_grid_vtype() is not callable");
        }
        static constexpr GridEL_type get_grid_type() {
            return _grid_type;
        }

        BodyModel<BM_vt> body; ///shared ptr to Body Model
        std::shared_ptr<GridEL_t> Grid; /// shared ptr to grid variants
        std::shared_ptr<LE_func_t > _LE; /// shared ptr to L(E) struct function
        std::shared_ptr<TrajPool_t> _TrajPools;

    private:
        void _init_trajs(TrajPoolInitParams_t TrajPoolInitParams) {
            auto const& mEL_Grid = Grid->inner(0);
            //_TrajPools = std::make_shared<>
            std::vector<bin_traj_pool_t<GridEL_vt>> mTP;
            mTP.reserve(mEL_Grid.size());
            for (auto m_bin : mEL_Grid) {
                mTP.push_back(
                    GetTrajPool<GridEL_vt>(
                        TrajPoolInitParams.is_static,
                        m_bin, *body, *_LE,
                        TrajPoolInitParams.T_error,
                        TrajPoolInitParams.traj_bins,
                        TrajPoolInitParams.n_e_max,
                        TrajPoolInitParams.n_l_max
                        )
                );
            }
            _TrajPools = std::make_shared<
                std::vector<bin_traj_pool_t<GridEL_vt>>
            >(std::move(mTP));
        }
    public:
        EL_Grid(
            BodyModel<BM_vt> const& body,
            GridEL_t _Grid, LE_func_t _LE_, TrajPool_t _TrajPools_
        ): body(body),
            Grid(std::make_shared<GridEL_t>(std::move(_Grid))),
            _LE(std::make_shared<LE_func_t>(std::move(_LE_))),
            _TrajPools(std::make_shared<TrajPool_t>(std::move(_TrajPools_))){}

        EL_Grid(
            BodyModel<BM_vt> const& body,
            GridEL_t _Grid, size_t NE_e_size,
            TrajPool_t _TrajPools_
        ) : body(body),
            Grid(std::make_shared<GridEL_t>(std::move(_Grid))),
            _LE(std::make_shared<typename Body<BM_vt>::LE_func_t>(
                body->get_le(NE_e_size)
                )
            ),
            _TrajPools(std::make_shared<TrajPool_t>(std::move(_TrajPools_)))
        {}

        EL_Grid(
            BodyModel<BM_vt> const& body,
            GridEL_t _Grid, LE_func_t _LE_,
            TrajPoolInitParams_t TrajPoolInitParams =
            TrajPoolInitParams_t()
        ) : body(body),
            Grid(std::make_shared<GridEL_t>(std::move(_Grid))),
            _LE(std::make_shared<LE_func_t>(std::move(_LE_)))
        {
            _init_trajs(TrajPoolInitParams);
        }

        template <typename GridType>
        EL_Grid(
            BodyModel<BM_vt> const& body, GridType&& EL_Grid,
            size_t NE_e_size = 100, 
            TrajPoolInitParams_t TrajPoolInitParams = 
            TrajPoolInitParams_t()
        ) :
            body(body),
            Grid(std::make_shared<GridEL_t>(std::forward<GridType>(EL_Grid))),
            _LE(std::make_shared<typename Body<BM_vt>::LE_func_t>(
                    body->get_le(NE_e_size)
                )
            )
        {
            _init_trajs(TrajPoolInitParams);
        }

        EL_Grid refine(size_t Er, size_t Lr) const {
            size_t NE_e_size = this->_LE->Internal_lm.size() - 1;
            return EL_Grid(
                body, 
                this->Grid->refine(Er,Lr),
                NE_e_size* Er
            );
        }

        inline EL_Func<BM_vt> getLE()const {
            return EL_Func<BM_vt>(body, _LE);
        }
        inline auto const& getBodyLE()const {
            return *_LE;
        }
        auto LE()const {
            return _LE->i_lm();
        }
        auto const& getLE_inner_grid() const {
            return Grid->inner(0);
        }

        using TrajPools_t = std::vector<bin_traj_pool_t<GridEL_vt>>;
        TrajPools_t const& TrajPools() const {
            return *_TrajPools;
        }

    };


    template <GridEL_type grid_type, typename BM_vt, typename GridType>
    auto make_core_GridEL(
        BodyModel<BM_vt> const& body, 
        GridType&& in_EL_Grid, 
        size_t NE_e_size = 100
    ) {
        using EL_t = decltype(in_EL_Grid.inner(0).grid()[0].left);
        //in_EL_Grid.inner(0) -- is EL grid
        //in_EL_Grid.inner(0).grid() -- is E grid
        //in_EL_Grid.inner(0).grid()[0] -- is Rect<T>
        //Rect<T>.left is T
        // we need only T
        return EL_Grid<BM_vt, EL_t, grid_type>(
            body, std::forward<GridType>(in_EL_Grid), 
            NE_e_size
        );
    }


}; //namespace evdm
