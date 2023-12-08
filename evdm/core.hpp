#ifndef CORE_HPP
#define CORE_HPP
#include "body_potential.hpp"
#include "grid_variants.hpp"
#include "histo_norm.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>

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
                typename IterableType>
    auto load_body_model(IterableType const &_array,size_t Rpoints){
        std::vector<T> Rho_values(Rpoints);
        for(size_t i=0;i<Rho_values.size();++i){
            Rho_values[i] = _array[i];
        }
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0,(T)1,Rpoints),
            std::move(Rho_values)
        );
        return forwarder::forward(
                    Body<T>::FromRho(std::move(Rho))
                );
    }

    template <typename forwarder = forward_value, typename T>
    auto load_body_model(std::vector<T> X) {
        std::vector<T> Rho_values(std::move(X));
        size_t size = Rho_values.size();
        
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0, (T)1, size),
            std::move(Rho_values)
        );
       
        
        return forwarder::forward(
            Body<T>::FromRho(std::move(Rho))
        );
    }

    template <typename T,typename forwarder = forward_value,
                typename Callable_t>
    auto make_body_model(Callable_t &&F,size_t Rpoints){
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
                    Body<T>::FromRho(std::move(Rho))
                );
    }



    template <typename T>
    struct EL_Func{
        BodyModel<T> body;
        std::shared_ptr<typename Body<T>::LE_func_t> LE;

        static T get_vtype() {
            static_assert("EL_Func::get_vtype() is not callable");
        }
        EL_Func(BodyModel<T> const& body,size_t Ne_bins):
            body(body),
            LE(
                std::make_shared<
                    typename Body<T>::LE_func_t
                >(body->get_le(Ne_bins))
            ){}
    };

    template <typename BodyModel_vt,typename GridEL_vt>
    struct EL_Grid{
        using BM_vt = BodyModel_vt;
        using GEL_vt = GridEL_vt;
        using GridEL_t = GridEL<GridEL_vt>;
        using Body_t = Body<BM_vt>;



        static constexpr BodyModel_vt get_bm_vtype() {
            static_assert("EL_Grid::get_bm_vtype() is not callable");
        }
        static constexpr GEL_vt get_grid_vtype() {
            static_assert("EL_Grid::get_grid_vtype() is not callable");
        }

        BodyModel<BM_vt> body; ///shared ptr to Body Model
        std::shared_ptr<GridEL_t> Grid; /// shared ptr to grid variants
        std::shared_ptr<typename Body<BM_vt>::LE_func_t > LE; /// shared ptr to L(E) struct function
        
        template <typename GridType>
        EL_Grid(BodyModel<BM_vt> const& body,GridType && EL_Grid,size_t NE_e_size=100):
            body(body),
            Grid(std::make_shared<GridEL_t>(std::forward<GridType>(EL_Grid))),
            LE(std::make_shared<typename Body<BM_vt>::LE_func_t>(body->get_le(NE_e_size)))
            {}
    };

    template <typename BM_vt,typename GridType>
    auto make_core_GridEL(BodyModel<BM_vt> const& body,GridType && in_EL_Grid, size_t NE_e_size = 100){
        using EL_t = decltype(in_EL_Grid.inner(0).grid()[0].left);
        //in_EL_Grid.inner(0) -- is EL grid
        //in_EL_Grid.inner(0).grid() -- is E grid
        //in_EL_Grid.inner(0).grid()[0] -- is Rect<T>
        //Rect<T>.left is T
        // we need only T
        return EL_Grid<BM_vt,EL_t>(body,std::forward<GridType>( in_EL_Grid));
    }

    template <typename T>
    using VectorX = std::conditional_t<
            std::is_same_v<T,double>,
            Eigen::VectorXd ,
            Eigen::VectorXf >;
    
    template <typename T>
    using MatrixX = std::conditional_t<
            std::is_same_v<T,double>,
            Eigen::MatrixXd ,
            Eigen::MatrixXf >;

    template <typename GridEL_vt,typename T>
    struct Distribution{
        using GEL_vt = GridEL_vt;
        using GridEL_t = GridEL<GEL_vt>;
        std::shared_ptr<GridEL_t> __Grid;
        VectorX<T> Values;


        static constexpr T get_vtype() {
            static_assert("Distribution::get_vtype() is not callable");
        }
        static constexpr GridEL_vt get_grid_vtype() {
            static_assert("Distribution::get_grid_vtype() is not callable");
        }

        auto const &Grid()const{
            return *__Grid;
        }
        Distribution(std::shared_ptr<GridEL_t> const &__Grid,
                     VectorX<T> Values = {}) : 
            __Grid(__Grid),Values(std::move(Values)){
                if(this->Values.size() == 0){
                    this->Values.resize(__Grid->size());
                    this->Values *= 0;
                }
            }
        
        template <typename U>
        Distribution<GridEL_vt,U> as_type() const {
            return Distribution<GridEL_vt,U>(__Grid,Values);
        }
    };
    
    template <typename V1,typename T1,typename V2,typename T2,
                typename Deg_t,typename BinNorm_t>
    auto distrib_norm(Distribution<V1,T1> const& D1,
                        Distribution<V2,T2> const& D2,
                        Deg_t deg,
                        BinNorm_t const & mu)
    {
        return std::visit([&](auto const & Grid1,auto const & Grid2){
                return histo_norm(grob::make_histo_view(Grid1,D1.Values),
                              grob::make_histo_view(Grid2,D2.Values),
                              deg,mu);
            },D1.Grid(),D2.Grid());
    } 

};
#endif//CORE_HPP