#include "core_matrix.hpp"
#include "../dynamics/dynamics_annihilation.hpp"

namespace evdm {
	
    template <typename T, typename Body_vt,
        typename GridEL_vt, GridEL_type grid_type>
    struct GridAnnPreMatrix {
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt, GEL_vt, grid_type>;
        GridEL_t Grid;
        Matrix_t<T> A0, Av;

        GridAnnPreMatrix(
            GridEL_t const& Grid,
            Matrix_t<T> _A0,
            Matrix_t<T> _Av
        ):Grid(Grid),A0(std::move(_A0)),Av(std::move(_Av)){
            size_t N = Grid.getLE_inner_grid().size();
            if (A0.cols() != N || A0.rows() != N) {
                throw std::out_of_range("GridAnnPreMatrix: A0.cols() != N || A0.rows() != N");
            } 
            if (Av.cols() != N || Av.rows() != N) {
                throw std::out_of_range("GridAnnPreMatrix: Av.cols() != N || Av.rows() != N");
            }
        }  


        template <typename Gen_t>
        GridAnnPreMatrix(
            GridEL_t const &Grid,
            GridEL_vt Rmin, GridEL_vt Rmax,
            size_t Nmk_traj, 
            Gen_t && G,progress_omp_function<> m_progress):
            Grid(Grid){
            const auto& grid = Grid.getLE_inner_grid();
            const size_t N = grid.size();
            A0.resize(N, N);
            A0.setZero();
            Av.resize(N, N);
            Av.setZero();
            auto _F2 = Grid.body->_DD_F_C(1);

            AnnImpl(A0, Av,
                grid, Rmin, Rmax, *Grid._TrajPools,
                Grid.LE(), Grid.body->Phi, _F2,
                G, Nmk_traj, m_progress);
            Av = 0.5 * (Av + Av.transpose());
            A0 = 0.5 * (A0 + A0.transpose());
        }

        template <typename U>
        void add_to(GridMatrix<U, Body_vt, GridEL_vt, grid_type>& M, size_t ptype0, size_t ptype1,T a_0, T a_v) const {
            auto pre_result = (a_0 * A0 + a_v * Av).template cast<U>();
            M.block(ptype0, ptype1).noalias() += pre_result;
            if(ptype0== ptype1)
                M.block(ptype0, ptype1).diagonal() += pre_result.diagonal();
        }
        //A0 - annihilation when sigma v ~ 1
        //Av - annihilation when sigma v ~ v^2

        template <typename U>
        GridAnnPreMatrix(
            GridAnnPreMatrix<U, Body_vt, GridEL_vt, grid_type> const& original
        ) : Grid(original.Grid),
            A0(original.A0.template cast<T>()),
            Av(original.Av.template cast<T>()) {}

        GridAnnPreMatrix(GridAnnPreMatrix const& original) :
            Grid(original.Grid),
            A0(original.A0),
            Av(original.Av)
        {}

        template <typename U>
        GridAnnPreMatrix<U, Body_vt, GridEL_vt, grid_type>
            as_type() const {
            return GridAnnPreMatrix<
                U, Body_vt, GridEL_vt, grid_type
            >(*this);
        }

    };

};
 
