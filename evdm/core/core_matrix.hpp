#pragma once

#include "core_grid.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>
#include <ranges>
namespace evdm {

    
    template <typename T>
    using Matrix_t =
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;


    template <typename T>
    using SpMatrix_t =
        Eigen::SparseMatrix<T,Eigen::ColMajor>;

    template <typename T>
    using SpTriplet_t = Eigen::Triplet<T>;

    template<typename T>
    std::vector<T> flatten(const std::vector<std::vector<T>>& nested) {
        auto joined = nested | std::views::join;
        return { joined.begin(), joined.end() };
    }


    template <typename T>
    std::vector<SpTriplet_t<T>> getTriplets(const SpMatrix_t<T> &M) {
        std::vector<SpTriplet_t<T>> outTriplets;
        outTriplets.reserve(M.nonZeros());
        for (int k = 0; k < M.outerSize(); ++k) {
            for (typename SpMatrix_t<T>::InnerIterator it(M, k); it; ++it) {
                outTriplets.push_back(
                    SpTriplet_t<T>( it.row(), it.col(),it.value() )
                );
            }
        }
        return outTriplets;
    }
    template <typename T, typename U>
    std::vector<SpTriplet_t<U>> getTriplets(
        const SpMatrix_t<T>& M,std::type_identity<U> 
    ) {
        std::vector<SpTriplet_t<U>> outTriplets;
        outTriplets.reserve(M.nonZeros());
        for (int k = 0; k < M.outerSize(); ++k) {
            for (typename SpMatrix_t<T>::InnerIterator it(M, k); it; ++it) {
                outTriplets.push_back(
                    SpTriplet_t<U>(it.row(), it.col(), it.value())
                );
            }
        }
        return outTriplets;
    }

    template <typename T>
    bool isDiagonalFullyStored(const SpMatrix_t<T>& mat) {
        int n = std::min(mat.rows(), mat.cols());
        for (int i = 0; i < n; ++i) {
            bool found = false;
            for (typename SpMatrix_t<T>::InnerIterator it(mat, i); it; ++it) {
                if (it.row() == i) { // Found Diagonal Element
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        return true;
    }

    template <typename T, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type>
    struct GridMatrix {
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt, GEL_vt, grid_type>;

        GridEL_t Grid;

        typedef SpMatrix_t<T> Mat_t;
        
    private:
        Mat_t Values;

    public:

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
        inline Mat_t & values() {
            return Values;
        }
        inline Mat_t const & values() const {
            return Values;
        }


        template <typename Index>
        auto as_histo(Index const& _I) const{
            size_t _shift = grid().LinearIndex(_I);
            auto m_grid = grid();
            return grob::make_histo(
                m_grid,
                Eigen::VectorXd(Values.col(_shift))
            );
        }

        inline auto as_histo_i(size_t i)  const {
            auto m_grid = grid();
            return grob::make_histo(
                m_grid,
                Eigen::VectorXd(Values.col(i))
            );
        }

        template <typename MatrixArray_t = Mat_t>
        GridMatrix(GridEL_t const& Grid,
            Mat_t _Values = {}) :
            Grid(Grid), Values(std::move(_Values)) 
        {
            size_t N = Grid.Grid->size();
            if (Values.rows() == 0 || Values.cols() == 0) {
                Values.resize(N, N);
            }
            if (Values.rows() != N || Values.cols() != N) {
                std::ostringstream ErrorText;
                ErrorText << "evdm::GridMatrix constructor "
                    "error: matrix size(" << 
                    Values.rows() << "*" << Values.cols() <<
                    ") doesn't match grid size (" << N << "*" << N << ")";
                throw std::runtime_error(
                    ErrorText.str()
                );
            }
            std::vector<SpTriplet_t<T>> mOtherTriplets;
            int n = Values.rows();
            for (int i = 0; i < n; ++i) {
                bool found = false;
                for (typename SpMatrix_t<T>::InnerIterator it(Values, i); it; ++it) {
                    if (it.row() == i) { // Found Diagonal Element
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    mOtherTriplets.push_back({ i,i,0 });
                }
            }
            if (mOtherTriplets.size()) {
                auto mTriplets = getTriplets(Values);
                mOtherTriplets.insert(
                    mOtherTriplets.end(), 
                    mTriplets.begin(),
                    mTriplets.end()
                );
                Values.setFromTriplets(
                    mOtherTriplets.begin(), mOtherTriplets.end()
                );
            }
        }


        template <typename U,typename T1,typename T2,GridEL_type grt1>
        friend struct GridMatrix;

        template <typename U>
        GridMatrix(
            GridMatrix<U, Body_vt, GridEL_vt, grid_type> const& original
        ): 
            GridMatrix(
                original.Grid, 
                original.Values.template cast<T>()
            )
        {}

        GridMatrix(GridMatrix const& original) :
            GridMatrix(
                original.Grid,
                original.Values) {}

        template <typename U>
        GridMatrix<U, Body_vt, GridEL_vt, grid_type> 
            as_type() const {
            return GridMatrix<
                U, Body_vt, GridEL_vt, grid_type
            >(*this);
        }

        auto _block(
            int ptype_in, int ptype_out
        ) const {
            size_t N_size = grid().inner(0).size();
            size_t N_size_full = grid().size();
            if (ptype_in < 0 || ptype_out < 0) {
                return Values.block(
                    0, 0,
                    N_size_full, N_size_full
                );
            }
            else {
                check_ptype(ptype_in);
                check_ptype(ptype_out);
                return Values.block(
                    ptype_out * N_size, ptype_in * N_size,
                    N_size, N_size
                );
            }
        }

        template <typename U>
        void add_to_block(SpMatrix_t<U> const& OtherMat,
            int ptype_in, int ptype_out) {
            auto mTriplets = getTriplets(OtherMat, std::type_identity<T>{});
            ptype_in = std::max(ptype_in, int(0));
            ptype_out = std::max(ptype_in, int(0));
            check_ptype(ptype_in);
            check_ptype(ptype_out);
            size_t N_size = grid().inner(0).size();
            
            int InShift = ptype_in * N_size;
            int OutShift = ptype_out * N_size;

            for (auto & tr : mTriplets) {
                tr = SpTriplet_t<T>(
                    tr.row() + OutShift, InShift + tr.col(), tr.value()
                );
            }
            SpMatrix_t<T> Other(Values.rows(), Values.cols());
            Other.setFromTriplets(mTriplets.begin(), mTriplets.end());
            Values += Other;
        }
        template <typename U>
        void add_to_block(Matrix_t<U> const& OtherMat,
            int ptype_in, int ptype_out) {
            return add_to_block(SpMatrix_t<U>(OtherMat.sparseView()), ptype_in, ptype_out);
        }
        Mat_t block(
            int ptype_in, int ptype_out
        ) const{
            size_t N_size = grid().inner(0).size();
            size_t N_size_full = grid().size();
            if (ptype_in < 0 || ptype_out < 0) {
                return Values;
            }
            else {
                check_ptype(ptype_in);
                check_ptype(ptype_out);
                return Values.block(
                    ptype_out * N_size, ptype_in * N_size,
                    N_size, N_size
                );
            }
        }
        inline Eigen::VectorX<T> diag(int ptype_in, int ptype_out) const{
            size_t N_size = grid().inner(0).size();
            size_t N_size_full = grid().size();
            if (ptype_in < 0 or ptype_out < 0) {
                return Values.diagonal().eval();
            }
            check_ptype(ptype_in);
            check_ptype(ptype_out);
            Eigen::VectorX<T> mDiag(N_size);
            mDiag.setZero();

            for (int i = N_size* ptype_in; i < N_size * (ptype_in +1); ++i) {
                for (typename SpMatrix_t<T>::InnerIterator it(Values, i); it; ++it) {
                    if (it.row() == i + N_size*(ptype_out- ptype_in)) {
                        // Found Diagonal Element
                        mDiag[i] = it.value();
                    }
                }
            }
            return mDiag;
        }

        inline auto raw_diag() const {
            return diag(-1,-1);
        }
        /// @brief makes diagonal elements to scatter
        template <typename  mVector_t>
        void to_scatter(mVector_t const & _EvapValues) {
            for (size_t i = 0; i < Values.rows(); ++i) {
                Values.coeffRef(i, i) = 0;
                auto _sum = Values.col(i).sum();
                Values.coeffRef(i, i) = -_sum - _EvapValues[i];
            }
        }
        /// @brief makes diagonal elements to scatter
        void to_scatter() {
            for (size_t i = 0; i < Values.rows(); ++i) {
                Values.coeffRef(i, i) = 0;
                auto _sum = Values.col(i).sum();
                Values.coeffRef(i, i) = -_sum;
            }
        }
    };

    template <
        typename T, typename Body_vt,     
        typename GridEL_vt, GridEL_type grid_type
    >
    auto make_Matrix(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid
    ) {
        GridMatrix<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid, {});
        return Dstr;
    }

    template <
        typename T, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type
    >
    auto make_Matrix(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid,
        SpMatrix_t<T> mSpMat
    ) {
        GridMatrix<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid, std::move(mSpMat));
        return Dstr;
    }

    /*
    template <
        typename T, typename Body_vt, 
        typename GridEL_vt, GridEL_type grid_type
    >   
    auto make_Matrix(
            EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid,
            const T * _mdata,size_t N, size_t padding) 
    {
        Eigen::Map<const Matrix_t<T>> data_view(_mdata, N, N);
        return GridMatrix<T, Body_vt, GridEL_vt, grid_type> 
            (Grid, data_view, padding);
    }

    template <
        typename T, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type
    >
    auto make_Matrix(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid,
        const T* _mdata, 
        size_t N,
        size_t _strid_outer, size_t _strid_inner, 
        size_t padding)
    {
        typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> stride_t;
        Eigen::Map<
            const Matrix_t<T>,
            Eigen::Unaligned,
            Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
        > data_view(_mdata, N, N, stride_t(_strid_inner, _strid_outer));
        return GridMatrix<T, Body_vt, GridEL_vt, grid_type>
            (Grid, data_view, padding);
    }*/

}; //namespace evdm