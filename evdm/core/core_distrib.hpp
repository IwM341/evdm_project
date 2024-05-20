#pragma once

#include "../histo_norm.hpp"
#include "../utils/progress_bar.hpp"
#include "../r_conversion.hpp"

#include "core_grid.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>

namespace evdm {
    template <typename T, typename Body_vt, typename GridEL_vt, GridEL_type grid_type>
    struct Distribution {
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt, GEL_vt, grid_type>;

        GridEL_t Grid;
        typedef Eigen::VectorX<T> _Vec_t;
        typedef Eigen::VectorBlock<Eigen::VectorX<T>> Vec_t;
        //using const_Vec_t = Eigen::VectorX<T>::ConstSegmentReturnType;
        typedef Eigen::VectorBlock<const Eigen::VectorX<T>> const_Vec_t;
    private:
        std::shared_ptr<Eigen::VectorX<T>> _Values;
        Vec_t _ValuesView;
        const_Vec_t _const_ValuesView;
        size_t padding;

        static Eigen::VectorX<T> process_vector(Eigen::VectorX<T> X, size_t _size, size_t padding) {
            size_t delta = _size % padding;
            if (X.size() >= _size && X.size() % padding == 0) {
                return std::move(X);
            }
            else {
                Eigen::VectorX<T> _X(delta ? _size + padding - delta : _size);
                _X.setZero();
                _X.segment(0, X.size()).noalias() = X;
                return std::move(_X);
            }
        }

    public:
        size_t get_padding()const {
            return padding;
        }

        void check_ptype(int ptype) const {
            if (!(ptype < grid().grid().size())) {
                throw std::runtime_error("ptype should be less than max ptypes");
            }
        }
        static constexpr Body_vt get_body_vtype() {
            static_assert("Distribution::get_body_vtype() is not callable");
        };
        static constexpr T get_vtype()
        {
            static_assert("Distribution::get_vtype() is not callable");
        }
        static constexpr GridEL_vt get_grid_vtype() {
            static_assert("Distribution::get_grid_vtype() is not callable");
        }
        static constexpr GridEL_type get_grid_type() {
            return grid_type;
        }

        inline _Vec_t& raw_vector() {
            return *_Values;
        }
        inline _Vec_t const& raw_vector()const {
            return *_Values;
        }

        inline auto const& grid()const {
            return *(Grid.Grid);
        }
        inline Body<Body_vt> const& body() const {
            return *Grid.body;
        }
        inline auto body_ptr() const {
            return Grid.body;
        }

        inline auto const& TrajPools() const {
            return Grid.TrajPools();
        }

        inline Vec_t values() {
            return _ValuesView;
        }
        inline const_Vec_t values() const {
            return _const_ValuesView;
        }
        inline auto as_histo() {
            return grob::make_histo_view(grid(), values());
        }
        inline auto as_histo() const {
            return grob::make_histo_view(grid(), values());
        }

        Distribution(GridEL_t const& Grid, Eigen::VectorX<T> Values = {}, size_t padding = 128) :
            Grid(Grid), _Values(std::make_shared< Eigen::VectorX<T>>(process_vector(std::move(Values), grid().size(), padding))),
            _ValuesView(_Values->segment(0, grid().size())),
            _const_ValuesView(static_cast<const Eigen::VectorX<T>&>(*_Values).segment(0, grid().size())),
            padding(padding)
        {}

        template <typename U>
        Distribution<U, Body_vt, GridEL_vt, grid_type> as_type() const {
            return Distribution<U, Body_vt, GridEL_vt, grid_type>(
                Grid, _Values->template cast<U>()
                );
        }

        T count(int ptype)const {
            if (ptype < 0) {
                return _const_ValuesView.sum();
            }
            else {
                check_ptype(ptype);
                size_t N_size = grid().inner(0).size();
                //Eigen::VectorXf Values;
                return _const_ValuesView.segment(ptype * N_size, N_size).sum();
            }
        }
        void CopyBroadcast(size_t ptype_src, size_t ptype_dst, T weight = 1) {
            check_ptype(ptype_src);
            check_ptype(ptype_dst);
            size_t N_size = grid().inner(0).size();
            _ValuesView.segment(ptype_dst * N_size, N_size) =
                _ValuesView.segment(ptype_src * N_size, N_size);
        }
        void SummBroadcast(size_t ptype_src, size_t ptype_dst, T weight = 1) {
            check_ptype(ptype_src);
            check_ptype(ptype_dst);
            size_t N_size = grid().inner(0).size();
            _ValuesView.segment(ptype_dst * N_size, N_size) +=
                _ValuesView.segment(ptype_src * N_size, N_size) * weight;
        }
        Eigen::VectorBlock<Eigen::VectorX<T>> block(int ptype) {
            size_t N_size = grid().inner(0).size();
            if (ptype >= 0) {
                check_ptype(ptype);
                return _Values->segment(ptype * N_size, N_size);
            }
            else {
                return _ValuesView;
            }

        }
        const Eigen::VectorBlock<const Eigen::VectorX<T>> block(int ptype) const {
            size_t N_size = grid().inner(0).size();
            if (ptype >= 0) {
                check_ptype(ptype);
                return static_cast<const Eigen::VectorX<T> &>(*_Values).
                    segment(ptype * N_size, N_size);
            }
            else {
                return _const_ValuesView;
            }
        }

        template <typename Gen_t, typename RGrid_t>
        auto r_dens(
            std::vector<size_t> ptypes, Gen_t&& G, RGrid_t&& RGrid,
            size_t Nmk_per_bin = 1000,
            progress_omp_function<> m_prog_bar_func = progress_omp_function<>()
        ) const {
            for (auto ptype : ptypes) {
                check_ptype(ptype);
            }
            return convert_r_density(
                as_histo(), TrajPools(),
                ptypes, Grid.LE(), body().Phi, body()._DD_F_C(1),
                G, std::forward<RGrid_t>(RGrid), Nmk_per_bin, m_prog_bar_func
            );
        }
    };
    template <
        typename T, typename Init_t,
        typename Body_vt, typename GridEL_vt,
        GridEL_type grid_type
    >
    auto make_Distribution(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid, Init_t&& init)
    {
        Distribution<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid);
        auto LE = Dstr.Grid.LE();
        if constexpr (!std::is_same_v<bool, std::decay_t<Init_t>>) {
            auto H1 = Dstr.as_histo();
            for (auto [bin, val] : H1) {
                auto type = std::get<0>(bin);
                auto [e0, e1] = std::get<0>(bin.tail());
                auto [l0, l1] = std::get<1>(bin.tail());
                auto Lm0 = LE(-e0);
                auto Lm1 = LE(-e1);
                val = evdm::integrateAB2([&](auto e) {
                    auto Lm = LE(-e);
                return evdm::integrateAB2([&](auto l) {
                    return init(e, l);
                    }, Lm * l0, Lm * l1, 1);
                    }, e0, e1, 1);
            }
        }
        return Dstr;
    }

    template <typename V1, typename BT1, typename T1, GridEL_type GT1,
        typename V2, typename BT2, typename T2, GridEL_type GT2,
        typename Deg_t, typename BinMeasure_t>
    auto distrib_norm(
        Distribution<V1, BT1, T1, GT1> const& D1,
        Distribution<V2, BT2, T2, GT2> const& D2,
        Deg_t deg,
        BinMeasure_t mu
    ) {
        auto m_bin_norm = get_bin_mes(mu, D1.Grid.LE());
        return histo_norm(
            D1.as_histo(),
            D2.as_histo(),
            deg, m_bin_norm
        );
    }

    template <
        typename T,
        typename Body_vt, typename GridEL_vt,
        GridEL_type grid_type
    >
    auto make_Distribution_data(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid,
        T* _data, size_t _size, size_t padding
    ) {
        size_t _size_grid = Grid.Grid->size();
        size_t _new_size = _size_grid + padding - _size_grid % padding;
        Eigen::VectorX<T> X(_new_size);
        X.setZero();
        size_t _copy_size = std::min(_size, _new_size);
        Eigen::Map<Eigen::VectorX<T>, Eigen::Unaligned> Data_Map(_data, _copy_size);
        X.segment(0, _copy_size) = Data_Map;
        return Distribution<T, Body_vt, GridEL_vt, grid_type>(
            Grid,
            std::move(X)
            );
    }

    template <
        typename Array_t,
        typename Body_vt, typename GridEL_vt,
        GridEL_type grid_type
    >
    auto make_Distribution_array(EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid,
        Array_t&& _data, size_t padding) {
        
        using T = std::decay_t<decltype(_data[0])>;
        size_t _size = _data.size();

        size_t _size_grid = Grid.Grid->size();
        size_t _new_size = _size_grid + padding - _size_grid % padding;
        Eigen::VectorX<T> X(_new_size);
        X.setZero();
        size_t _copy_size = std::min(_size, _new_size);
        
        for (size_t i = 0; i < _copy_size; ++i) {
            X[i] = _data[i];
        }

        return Distribution<T, Body_vt, GridEL_vt, grid_type>(
            Grid,
            std::move(X)
            );
    }
}; //namespace evdm 
