#ifndef CORE_HPP
#define CORE_HPP
#include "body_potential.hpp"
#include "grid_variants.hpp"
#include "histo_norm.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>
#include "dynamics/dynamics.hpp"

namespace evdm{

    
    template <typename T>
    using BodyModel = std::shared_ptr<Body<T>>;

    struct forward_value{
        template <typename T>
        static decltype(auto) forward(T && x){
            return std::forward<T>(x);
        } 
    };
    struct forward_shared{
        template <typename T>
        static decltype(auto) forward(T && x){
            return std::make_shared<std::decay_t<T>>(std::forward<T>(x));
        } 
    };

    template <typename T,typename forwarder = forward_value,
                typename IterableType,typename U>
    auto load_body_model(IterableType const &_array,size_t Rpoints,U Velocity){
        std::vector<T> Rho_values(Rpoints);
        for(size_t i=0;i<Rho_values.size();++i){
            Rho_values[i] = _array[i];
        }
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0,(T)1,Rpoints),
            std::move(Rho_values)
        );
        return forwarder::forward(
                    Body<T>::FromRho(Velocity,std::move(Rho))
                );
    }

    template <typename forwarder = forward_value, typename T,typename U>
    auto load_body_model(std::vector<T> X,U Velocity) {
        std::vector<T> Rho_values(std::move(X));
        size_t size = Rho_values.size();
        
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0, (T)1, size),
            std::move(Rho_values)
        );
       
        
        return forwarder::forward(
            Body<T>::FromRho(Velocity,std::move(Rho))
        );
    }

    template <typename T,typename forwarder = forward_value,
                typename Callable_t, typename U>
    auto make_body_model(Callable_t &&F,size_t Rpoints, U Velocity){
        std::vector<T> Rho_values(Rpoints);
        T h_inv = (T)1.0 / (Rpoints - 1);
        for(size_t i=0;i<Rho_values.size();++i){
            Rho_values[i] = F(i* h_inv);
        }
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0,(T)1,Rpoints),
            std::move(Rho_values)
        );
        return forwarder::forward(
                    Body<T>::FromRho(Velocity,std::move(Rho))
                );
    }



    template <typename T>
    struct EL_Func{
        BodyModel<T> body;
        std::shared_ptr<typename Body<T>::LE_func_t> LE;

        static T get_vtype() {
            static_assert("EL_Func::get_vtype() is not callable");
        }
        EL_Func(BodyModel<T> const& body, decltype(LE) const& LE):
            body(body),LE(LE){}
        EL_Func(BodyModel<T> const& body,size_t Ne_bins):
            body(body),
            LE(
                std::make_shared<
                    typename Body<T>::LE_func_t
                >(body->get_le(Ne_bins))
            ){}
    };

    template <typename BodyModel_vt,typename GridEL_vt, GridEL_type _grid_type>
    struct EL_Grid {
        using BM_vt = BodyModel_vt;
        using GEL_vt = GridEL_vt;
        using GridEL_t = GridEL<GridEL_vt, _grid_type>;
        using Body_t = Body<BM_vt>;



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
        std::shared_ptr<typename Body<BM_vt>::LE_func_t > _LE; /// shared ptr to L(E) struct function

        template <typename GridType>
        EL_Grid(BodyModel<BM_vt> const& body, GridType&& EL_Grid, size_t NE_e_size = 100) :
            body(body),
            Grid(std::make_shared<GridEL_t>(std::forward<GridType>(EL_Grid))),
            _LE(std::make_shared<typename Body<BM_vt>::LE_func_t>(body->get_le(NE_e_size)))
        {}

        inline EL_Func<BM_vt> getLE()const {
            return EL_Func<BM_vt>(body, _LE);
        }
        auto LE()const {
            return _LE->i_lm();
        }

    };


    template <GridEL_type grid_type,typename BM_vt,typename GridType>
    auto make_core_GridEL(BodyModel<BM_vt> const& body,GridType && in_EL_Grid, size_t NE_e_size = 100){
        using EL_t = decltype(in_EL_Grid.inner(0).grid()[0].left);
        //in_EL_Grid.inner(0) -- is EL grid
        //in_EL_Grid.inner(0).grid() -- is E grid
        //in_EL_Grid.inner(0).grid()[0] -- is Rect<T>
        //Rect<T>.left is T
        // we need only T
        return EL_Grid<BM_vt,EL_t, grid_type>(body,std::forward<GridType>( in_EL_Grid));
    }



    template <typename T,typename Body_vt,typename GridEL_vt, GridEL_type grid_type>
    struct Distribution{
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt,GEL_vt, grid_type>;
        
        GridEL_t Grid;
        std::shared_ptr<Eigen::VectorX<T>> Values;
        

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

        inline auto const &grid()const{
            return *(Grid.Grid);
        }
        inline auto& body() const {
            return Grid.body;
        }

        inline auto const& values() const {
            return *Values;
        } 
        inline auto & values() {
            return *Values;
        }
        inline auto as_histo() {
            return grob::make_histo_view(grid(), values());
        }
        inline auto as_histo() const {
            return grob::make_histo_view(grid(), values());
        }

        Distribution(GridEL_t const &Grid,
                     Eigen::VectorX<T> Values = {}) : 
            Grid(Grid),Values(std::make_shared< Eigen::VectorX<T>>(std::move(Values))){
                if(this->Values->size() == 0){
                    this->Values->resize(Grid.Grid->size());
                    this->Values->setZero();
                }
            }
        
        template <typename U>
        Distribution<U, Body_vt, GridEL_vt, grid_type> as_type() const {
            return Distribution<U, Body_vt, GridEL_vt, grid_type>(Grid,Values->template cast<U>());
        }

        T count(int ptype)const {
            if (ptype < 0) {
                return Values->sum();
            }
            else {
                size_t N_size = grid().inner(0).size();
                //Eigen::VectorXf Values;
                return Values->block(ptype * N_size, 0, N_size, 1).sum();
            }
        }
    };
    template <typename T,typename Init_t, typename Body_vt, typename GridEL_vt, GridEL_type grid_type>
    auto make_Distribution(EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid, Init_t&& init) {
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
    
    template <typename V1,typename BT1,typename T1, GridEL_type GT1,
              typename V2,typename BT2,typename T2, GridEL_type GT2,
                typename Deg_t,typename BinMeasure_t>
    auto distrib_norm(
        Distribution<V1, BT1,T1, GT1> const& D1,
        Distribution<V2, BT2,T2, GT2> const& D2,
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

    template <typename T>
    using Matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    template <typename T, typename Body_vt, typename GridEL_vt, GridEL_type grid_type>
    struct GridMatrix{
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt, GEL_vt, grid_type>;

        GridEL_t Grid;
        std::shared_ptr<Matrix_t<T>> Values;


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

        inline auto const& grid()const {
            return *(Grid.Grid);
        }
        inline auto const& values() const {
            return *Values;
        }
        inline auto& values() {
            return *Values;
        }

        template <typename Index>
        inline auto as_histo(Index const & _I) {
            size_t _size = grid().size();
            size_t _shift = grid().LinearIndex(_I);
            return grob::make_histo_view(grid(), values().data() + _shift);
        }
        template <typename Index>
        inline auto as_histo(Index const& _I) const {
            size_t _shift = grid().LinearIndex(_I);
            return grob::make_histo_view(grid(), values().data() + _shift);
        }

        GridMatrix(GridEL_t const& Grid,
            Matrix_t<T> Values = {}) :
            Grid(Grid), Values(std::make_shared< Matrix_t<T>>(std::move(Values))) {
            if (this->Values->size() == 0) {
                Values = Matrix_t<T>(Grid.Grid->size(),Grid.Grid->size());
                this->Values->setZero();
            }
        }

        template <typename U>
        GridMatrix<U, Body_vt, GridEL_vt, grid_type> as_type() const {
            return GridMatrix<U, Body_vt, GridEL_vt, grid_type>(Grid, Values->template cast<U>());
        }
    };

    template <typename T, typename Body_vt, typename GridEL_vt, GridEL_type grid_type>
    auto make_Matrix(EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid) {
        GridMatrix<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid);
        return Dstr;
    }
};
#endif//CORE_HPP