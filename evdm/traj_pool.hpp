#ifndef TRAJ_POOL_HPP
#define TRAJ_POOL_HPP
#include <grob/grid_objects.hpp>
#include "grid_variants.hpp"
#include "body_potential.hpp"
#include "trajectory.hpp"
/*
*   In this headre we define array of
*   trajectories, and periods
*/

namespace evdm{

    template <typename T>
    struct traj_info{
        T rmin,rmax;
        T T_in_theta; /// regularized by theta period
        T theta_max;
        typedef grob::GridFunction<grob::splineD1D,grob::GridUniform<T>,std::vector<T>> 
            GridFunction_t;
        GridFunction_t theta_tau;
    };
    
    struct hidden_reg_l{
        template <typename T>
        inline auto operator ()(T const & l) const{
            return std::sqrt(1-l*l);
        }
    };
    struct unhidden_reg_l:public hidden_reg_l{};
    
    template <typename Tr_Type>
    using Bin_e_grid_t = grob::GridUniform<Tr_Type>;
    
    template <typename Tr_Type>
    using Bin_l_grid_t = grob::GridFunctional<Tr_Type,hidden_reg_l,unhidden_reg_l>;

    template <typename Tr_Type>
    using Bin_el_grid_t = decltype(
        grob::mesh_grids(
            std::declval<Bin_e_grid_t<Tr_Type>>(),
            std::declval<Bin_l_grid_t<Tr_Type>>()
        ) 
    );

    template <typename Tr_Type>
    using bin_traj_pool_t = grob::GridObject<
            Bin_el_grid_t<Tr_Type>,
            std::vector<traj_info<Tr_Type>
        >;

    //template <typename Tr_Type>
    //

    template <typename Tr_Type,typename U,typename T>
    bin_traj_pool_t<Tr_Type> GetTrajPool(
                            bool is_static,
                            const grob::Point<grob::Rect<U>,grob::Rect<U>> & dEdL,
                            Body<T> const & B,
                            typename Body<T>::LE_func_t const & LE,
                            const Tr_Type T_error = 0.05,
                            const size_t traj_bins = 100,
                            const size_t n_e_max= 2,const size_t n_l_max = 2)
    {
        auto Bin_e_grid = grob::GridUniform<Tr_Type>(std::get<0>(dEdL).left,std::get<0>(dEdL).right,2);
        auto Bin_l_p_grid = Bin_l_grid_t<Tr_Type>(std::get<1>(dEdL).left,
                                          std::get<1>(dEdL).right,2,
                                          hidden_reg_l{},
                                          unhidden_reg_l{});

        typedef traj_info<Tr_Type> traj_info_t;
        auto get_traj = [&LE,&B,traj_bins](U e, U l){
            auto L2 = LE.i_lmq(e)*l*l;
            auto [rm,rp] = B.find_rmin_rmax(e,L2,LE);
            auto m_traj = B.get_internal_traj(rm,rp,traj_bins);
            return traj_info_t {
                rm,rp,
                std::get<1>(m_traj),
                std::get<0>(m_traj),
                InverseTrajectory<
                    typename traj_info<Tr_Type>::GridFunction_t::interpolator_t
                >(std::get<2>(m_traj),traj_bins)
            };
        };
        if (is_static){
            auto mGrid = grob::mesh_grids(Bin_e_grid,Bin_l_p_grid);
            std::vector<traj_info_t> ti_values = {
                get_traj(std::get<0>(dEdL).left,std::get<1>(dEdL).left),
                get_traj(std::get<0>(dEdL).left,std::get<1>(dEdL).right),
                get_traj(std::get<0>(dEdL).right,std::get<1>(dEdL).left),
                get_traj(std::get<0>(dEdL).right,std::get<1>(dEdL).right),
            };
            return grob::make_grid_object(mGrid,std::move(ti_values));
        } else {
            std::vector<traj_info_t> TI = {
                    get_traj(std::get<0>(dEdL).left,std::get<1>(dEdL).left),
                    get_traj(std::get<0>(dEdL).left,std::get<1>(dEdL).right),
                    get_traj(std::get<0>(dEdL).right,std::get<1>(dEdL).left),
                    get_traj(std::get<0>(dEdL).right,std::get<1>(dEdL).right)
                };
            auto norm1 = [](auto x,auto y){return std::abs(x-y)/(x+y);};


            auto Calculate_dT = [&](){
                Tr_Type dT_e_max = 0;
                Tr_Type dT_l_max = 0;
                const size_t ne = Bin_e_grid.size();
                const size_t nl = Bin_l_p_grid.size();
                for(size_t il=0;il<nl;++il){
                    for(size_t ie=0;ie<ne-1;++ie){
                        Tr_Type dT = norm1(TI[(ie+1)*nl + il].T_in_theta,TI[ie*nl + il].T_in_theta);
                        dT_e_max = std::max(dT_e_max,dT);
                    }
                }
                for(size_t il=0;il<nl-1;++il){
                    for(size_t ie=0;ie<ne;++ie){
                        Tr_Type dT = norm1(TI[ie*nl + il].T_in_theta,TI[ie*nl + il+1].T_in_theta);
                        dT_l_max = std::max(dT_l_max,dT);
                    }
                }
                return std::make_pair(dT_e_max,dT_l_max);
            };
            auto Increase_NE = [&](){
                auto Bin_e_grid_tmp = Bin_e_grid;
                Bin_e_grid_tmp.resize(Bin_e_grid.size()*2-1);
                //   e0     e1     e2     ...
                //[*****][*****][*****][*****]
                //   |        |      |
                //   |        ----   -----------
                //   |           |             |
                //[*****][?????][*****][?????][*****][?????][*****]
                std::vector<traj_info_t> TI_tmp(Bin_e_grid_tmp.size()*Bin_l_p_grid.size());
                for(size_t i=0;i<Bin_e_grid_tmp.size();++i){
                    for(size_t j=0;j<Bin_l_p_grid.size();++j){
                        if(i%2){
                            TI_tmp[i*Bin_l_p_grid.size() + j] = 
                                get_traj(Bin_e_grid_tmp[i],Bin_l_p_grid[j]);
                        } else {
                            TI_tmp[i*Bin_l_p_grid.size() + j] = 
                                std::move(TI[(i/2)*Bin_l_p_grid.size() + j]);
                        }
                    }
                }
                Bin_e_grid=Bin_e_grid_tmp;
                TI = std::move(TI_tmp);
            };
            auto Increase_NL = [&](){
                auto Bin_l_grid_tmp = Bin_l_p_grid;
                Bin_l_grid_tmp.resize(Bin_l_p_grid.size()*2-1);
                // 12345 l12345 ...
                //[*****][*****]
                //[*?*?*?*?*][*?*?*?*?*]
                std::vector<traj_info_t> TI_tmp(Bin_e_grid.size()*Bin_l_grid_tmp.size());
                for(size_t i=0;i<Bin_e_grid.size();++i){
                    for(size_t j=0;j<Bin_l_grid_tmp.size();++j){
                        if(j%2){
                            TI_tmp[i*Bin_l_grid_tmp.size() + j] = 
                                get_traj(Bin_e_grid[i],Bin_l_grid_tmp[j]);
                        } else {
                            TI_tmp[i*Bin_l_grid_tmp.size() + j] = 
                                std::move(TI[i*Bin_l_p_grid.size() + j/2]);
                        }
                    }
                }
                Bin_l_p_grid.hidden() = Bin_l_grid_tmp.hidden();
                TI = std::move(TI_tmp);
            };
            
            for(size_t i=0;i< n_e_max*n_l_max;++i){
                auto [eps_e,eps_l] = Calculate_dT();
                if( (Bin_e_grid.size()*2-1 > n_e_max || eps_e <= T_error) &&
                    (Bin_l_p_grid.size()*2-1 > n_l_max || eps_l < T_error))
                    {
                        break;   
                    }
                if(eps_e > eps_l){
                    Increase_NE();
                    if(eps_l > T_error && Bin_l_p_grid.size()*2-1 <= n_l_max){
                        Increase_NL();
                    }
                }
                else {
                    Increase_NL();
                    if(eps_e > T_error && Bin_e_grid.size()*2-1 <= n_e_max){
                        Increase_NE();
                    }
                }
            }

            auto mGrid = grob::mesh_grids(Bin_e_grid,Bin_l_p_grid);
            return grob::make_grid_object(mGrid,std::move(TI));
        }
    }
};

#endif//TRAJ_POOL_HPP