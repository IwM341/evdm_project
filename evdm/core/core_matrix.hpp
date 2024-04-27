#pragma once

#include "core_grid.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>

namespace evdm {

    template <typename T>
    using Matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    template <typename T, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type>
    struct GridMatrix {
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt, GEL_vt, grid_type>;

        GridEL_t Grid;

        typedef Matrix_t<T> _Mat_t;
        typedef Eigen::Block<Matrix_t<T>> Mat_t;
        //using const_Vec_t = Eigen::VectorX<T>::ConstSegmentReturnType;
        typedef Eigen::Block<const Matrix_t<T>> const_Mat_t;
    private:
        std::shared_ptr<_Mat_t> _Values;
        Mat_t _ValuesView;
        const_Mat_t _const_ValuesView;
        size_t padding;


        static _Mat_t process_matrix(
            _Mat_t X, size_t _N,
            size_t padding) {
            size_t delta = _N % padding;
            if (X.cols() >= _N && X.size() % padding == 0) {
                return std::move(X);
            }
            else {
                size_t _N_true = delta ? _N + padding - delta : _N;
                _Mat_t _X(_N_true, _N_true);
                _X.setZero();
                _X.block(0, 0, X.cols(), X.rows()).noalias() = X;
                return std::move(_X);
            }
        }

    public:

        size_t get_padding()const {
            return padding;
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


        void check_ptype(int ptype) const {
            if (!(ptype < grid().grid().size())) {
                throw std::runtime_error("ptype should be less than max ptypes");
            }
        }


        inline Body<Body_vt> const& body() const {
            return *Grid.body;
        }
        inline auto body_ptr() const {
            return Grid.body;
        }
        inline auto const& grid()const {
            return *(Grid.Grid);
        }
        inline Mat_t values() {
            return _ValuesView;
        }
        inline const_Mat_t values() const {
            return _const_ValuesView;
        }
        inline _Mat_t& raw_matrix() {
            return *_Values;
        }
        inline _Mat_t const& raw_matrix() const {
            return *_Values;
        }

        template <typename Index>
        inline auto as_histo(Index const& _I) {
            size_t _shift = grid().LinearIndex(_I);
            return grob::make_histo_view(grid(), _ValuesView.col(_shift));
        }

        template <typename Index>
        inline auto as_histo(Index const& _I) const {
            size_t _shift = grid().LinearIndex(_I);
            return grob::make_histo_view(grid(), _ValuesView.col(_shift));
        }

        inline auto as_histo_i(size_t i) {
            return grob::make_histo_view(grid(), _ValuesView.col(i));
        }
        inline auto as_histo_i(size_t i) const {
            return grob::make_histo_view(grid(), _ValuesView.col(i));
        }

        inline auto as_histo_bank(size_t ptype_in, size_t ptype_out) {
            check_ptype(ptype_in);
            check_ptype(ptype_out);
            auto const& _grid = grid().inner(0); //geom grid without ptypes
            size_t inner_size = _grid.size(); //size of geom grid
            size_t ptypes = grid().grid().size(); //number of ptypes
            return grob::make_grid_object_ref(
                _grid,
                grob::as_container(
                    [&_grid, m_block = block(ptype_in, ptype_out)](size_t i)mutable {
                        auto m_col = m_block.col(i);
            return grob::make_histo_view(_grid,
                grob::vector_view(m_col.data(), m_col.size())
            );
                    }, _grid.size()
                        )
            );
        }

        GridMatrix(GridEL_t const& Grid,
            Matrix_t<T> Values = {}, size_t padding = 128) :
            Grid(Grid), _Values(
                std::make_shared< Matrix_t<T>>(
                    process_matrix(std::move(Values), Grid.Grid->size(), padding)
                    )),
            _ValuesView(_Values->block(0, 0, grid().size(), grid().size())),
            _const_ValuesView(
                static_cast<_Mat_t const&>(*_Values)
                .block(0, 0, grid().size(), grid().size())
            ), padding(padding) {}

        template <typename U>
        GridMatrix<U, Body_vt, GridEL_vt, grid_type> as_type() const {
            return GridMatrix<U, Body_vt, GridEL_vt, grid_type>(
                Grid, _Values->template cast<U>()
                );
        }

        void CopyBroadcast(
            size_t ptype_in_src, size_t ptype_out_src,
            size_t ptype_in_dst, size_t ptype_out_dst, T weight = 1) 
        {
            block(ptype_in_dst, ptype_out_dst) = block(ptype_in_src, ptype_out_src) * weight;
        }
        void SummBroadcast(
            size_t ptype_in_src, size_t ptype_out_src,
            size_t ptype_in_dst, size_t ptype_out_dst, T weight = 1) {

            block(ptype_in_dst, ptype_out_dst) += block(ptype_in_src, ptype_out_src) * weight;
        }

        inline Eigen::Block<Matrix_t<T>> block(int ptype_in, int ptype_out) {



            size_t N_size = grid().inner(0).size();
            size_t N_size_full = grid().size();
            if (ptype_in < 0 || ptype_out < 0) {
                return _ValuesView;
            }
            else {
                check_ptype(ptype_in);
                check_ptype(ptype_out);
                return raw_matrix().block(
                    ptype_out * N_size, ptype_in * N_size,
                    N_size, N_size
                );
            }
        }
        inline auto diag(int ptype_in, int ptype_out) {
            return block(ptype_in, ptype_out).diagonal();
        }
        inline auto diag(int ptype_in, int ptype_out) const{
            return block(ptype_in, ptype_out).diagonal();
        }
        inline auto raw_diag() {
            return raw_matrix().diagonal();
        }
        inline auto raw_diag() const{
            return raw_matrix().diagonal();
        }

        inline Eigen::Block<const Matrix_t<T>> block(int ptype_in, int ptype_out)const {



            size_t N_size = grid().inner(0).size();
            size_t N_size_full = grid().size();
            if (ptype_in < 0 || ptype_out < 0) {
                return _const_ValuesView;
            }
            else {
                check_ptype(ptype_in);
                check_ptype(ptype_out);
                return raw_matrix().block(
                    ptype_out * N_size, ptype_in * N_size,
                    N_size, N_size
                );
            }
        }

        /// @brief makes diagonal elements to scatter
        template <typename  mVector_t>
        void to_scatter(mVector_t const & _EvapValues) {
            for (size_t i = 0; i < _ValuesView.rows(); ++i) {
                _ValuesView(i, i) = 0;
                auto _sum = _ValuesView.col(i).sum();
                _ValuesView(i, i) = -_sum - _EvapValues[i];
            }
        }
        /// @brief makes diagonal elements to scatter
        void to_scatter() {
            for (size_t i = 0; i < _ValuesView.rows(); ++i) {
                _ValuesView(i, i) = 0;
                auto _sum = _ValuesView.col(i).sum();
                _ValuesView(i, i) = -_sum;
            }
        }
    };

    template <typename T, typename Body_vt, typename GridEL_vt, GridEL_type grid_type>
    auto make_Matrix(EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid) {
        GridMatrix<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid);
        return Dstr;
    }

}; //namespace evdm