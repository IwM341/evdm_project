#pragma once

#include "../histo_norm.hpp"
#include "../utils/progress_bar.hpp"
#include "../r_conversion.hpp"

#include "core_grid.hpp"
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>
//#include <pybind11/pybind11.h>

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
        Eigen::VectorX<T> _Values;
        Vec_t _ValuesView;
        const_Vec_t _const_ValuesView;

        static Eigen::VectorX<T> process_vector(
            Eigen::VectorX<T> X, size_t _size
        ) {
            if (X.size() == _size) {
                return std::move(X);
            }
            else {
                Eigen::VectorX<T> _X(_size);
                _X.setZero();
                _X.segment(0, std::min((size_t)X.size(),_size) ).noalias() = X;
                return std::move(_X);
            }
        }

    public:
        size_t get_padding()const {
            return 1;
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
            return _Values;
        }
        inline _Vec_t const& raw_vector()const {
            return _Values;
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
            return grob::make_histo_view(grid(), grob::vector_view<T>(values().data(), grid().size()));
        }
        inline auto as_histo() const {
            return grob::make_histo_view(grid(), grob::vector_view<const T>(values().data(), grid().size()));
        }

        Distribution(GridEL_t const& Grid, Eigen::VectorX<T> Values = {}, size_t padding = 1) :
            Grid(Grid), _Values(process_vector(std::move(Values), grid().size())),
            _ValuesView(_Values.segment(0, grid().size())),
            _const_ValuesView(static_cast<const Eigen::VectorX<T>&>(_Values).segment(0, grid().size()))
        {}

        template <
            typename U, typename Body_vt1,
            typename GridEL_vt1, GridEL_type grid_type1
        >
        friend struct Distribution;

        template <typename U>
        Distribution(
            Distribution<U, Body_vt, GridEL_vt, grid_type> const& _original
        ): Distribution(
            _original.Grid, 
            _original._Values.template cast<T>()
           ){}

        Distribution(Distribution  const& _original ) : Distribution(
            _original.Grid,
            _original._Values
        ) {}

        template <typename U>
        Distribution<U, Body_vt, GridEL_vt, grid_type> as_type() const {
            return Distribution<U, Body_vt, GridEL_vt, grid_type>(*this);
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
                return _Values.segment(ptype * N_size, N_size);
            }
            else {
                return _ValuesView;
            }

        }
        const Eigen::VectorBlock<const Eigen::VectorX<T>> block(int ptype) const {
            size_t N_size = grid().inner(0).size();
            if (ptype >= 0) {
                check_ptype(ptype);
                return static_cast<const Eigen::VectorX<T> &>(_Values).
                    segment(ptype * N_size, N_size);
            }
            else {
                return _const_ValuesView;
            }
        }

        template <typename Gen_t, typename RGrid_t,typename Nmk_distr_t>
        auto r_dens(
            std::vector<size_t> ptypes, Gen_t&& G, RGrid_t&& RGrid,
            Nmk_distr_t const & Nmk_distrib = 1000,
            progress_omp_function<> m_prog_bar_func = progress_omp_function<>()
        ) const {
            for (auto ptype : ptypes) {
                check_ptype(ptype);
            }
            return convert_r_density(
                as_histo(), TrajPools(),
                ptypes, Grid.LE(), body().Phi, body()._DD_F_C(1),
                G, std::forward<RGrid_t>(RGrid), Nmk_distrib, m_prog_bar_func
            );
        }
    };
    template <
        typename T, typename Init_t,
        typename Body_vt, typename GridEL_vt,
        GridEL_type grid_type
    >
    auto make_Distribution(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid, Init_t&& init,
        size_t padding = 1)
    {
        Distribution<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid, {}, padding);
        auto LE = Dstr.Grid.LE();
        if constexpr (!std::is_same_v<bool, std::decay_t<Init_t>>) {
            auto H1 = Dstr.as_histo();
            for (auto [bin, val] : H1) {
                auto type = std::get<0>(bin);
                auto [e0, e1] = std::get<0>(bin.tail());
                auto [l0, l1] = std::get<1>(bin.tail());
                val = evdm::integrateAB2([&](auto e) {
                    auto Lm = LE(-e);
                return evdm::integrateAB2([&](auto l) {
                    return init(e, l);
                    }, l0, l1, 1);
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
        const T* _data, size_t data_step, size_t _size, size_t padding
    ) {
        size_t _size_grid = Grid.Grid->size();
        size_t _new_size = _size_grid;
        Eigen::VectorX<T> X(_new_size);
        X.setZero();
        size_t _copy_size = std::min(_size, _new_size);
        if (data_step == 1) {
            Eigen::Map<const Eigen::VectorX<T>>
                Data_Map(_data, _copy_size);
            X.segment(0, _copy_size) = Data_Map;
            return Distribution<T, Body_vt, GridEL_vt, grid_type>(
                Grid,
                std::move(X)
                );
        }
        else {
            Eigen::Map<const Eigen::VectorX<T>, Eigen::Unaligned,
                Eigen::InnerStride<Eigen::Dynamic>
            > DataMap(_data, _copy_size, Eigen::InnerStride<Eigen::Dynamic>(data_step));
            X.segment(0, _copy_size) = DataMap;
            return Distribution<T, Body_vt, GridEL_vt, grid_type>(
                Grid,
                std::move(X)
                );
        }
        
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
        size_t _new_size = _size_grid;
        Eigen::VectorX<T> X(_size_grid);
        X.setZero();
        size_t _copy_size = std::min(_size, _size_grid);
        
        for (size_t i = 0; i < _copy_size; ++i) {
            X[i] = _data[i];
        }

        return Distribution<T, Body_vt, GridEL_vt, grid_type>(
            Grid,
            std::move(X)
            );
    }

    template <
        typename T,typename B_vt,typename G_vt,GridEL_type grid_type,
        typename Functor_t
    >
    T avarage_by_grid(
        Distribution<T,B_vt,G_vt, grid_type> const & dstr,
        int ptype,
        Functor_t && F_EL
    ) {
        auto get_sum_count = [&](size_t m_ptype) {
            auto H = dstr.as_histo();
            auto H1 = H.inner_slice(m_ptype);
            T sum = 0;
            T avarage = 0;

            for (auto [mbin, value] : H1) {
                avarage += value;
                auto [e, l] = mbin.center();
                sum += value * F_EL(e, l);
            }
            return std::tuple<T, T>(sum, avarage);
        };
        if (ptype >= 0) {
            auto [sum, count] = get_sum_count(ptype);
            return sum / count;
        }
        else {
            T full_sum = 0;
            T full_count = 0;
            size_t ptypemax = dstr.grid().ptypes();
            for (size_t ptype = 0; ptype < ptypemax; ++ptype) {
                auto [sum, count] = get_sum_count(ptype);
                full_sum += sum;
                full_count += count;
            }
            return full_sum / full_count;
        }
        
    }

    template <
        typename T, typename B_vt, typename G_vt, GridEL_type grid_type
    >
    std::vector<T> get_E_distrib(
        Distribution<T, B_vt, G_vt, grid_type> const& dstr,
        int ptype = -1
    ) {
        size_t gridEsize = dstr.grid().inner(0).grid().size();
        //pybind11::print("grid size = ", gridEsize);
        std::vector<T> CountValues(gridEsize,0);
        //pybind11::print("CountValues.size()", CountValues.size());
        //CountValues.resize(gridEsize);
        //for (auto& x : CountValues) { x = 0; }
        auto put_E_distrib = [&](size_t m_ptype,T * data) {
            auto H = dstr.as_histo();
            auto H1 = H.inner_slice(m_ptype);
            for (size_t i = 0; i < gridEsize; ++i) {
                data[i] += H1.inner_slice(i).count();
            }
        };

        
        
        if (ptype >= 0) {
            put_E_distrib(ptype, CountValues.data());
        }
        else {
            size_t ptypemax = dstr.grid().ptypes();
            for (size_t pt = 0; pt < ptypemax; ++pt) {
                put_E_distrib(pt, CountValues.data());
            }
        }
        return CountValues;
    }
}; //namespace evdm 
