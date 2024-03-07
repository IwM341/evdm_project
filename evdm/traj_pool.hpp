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
        T rmin,rmax;  ///min and max r of trajectory
        T dul, su, su_minus_2_q;/// (umax-umin)/sqrt(1-l^2), umax+umin, (umax+umin-2)/q_e
        T u_delta_0z, u_delta_1z;
        T T_in_theta; /// regularized by theta period
        T theta_max;  /// max theta parametr of trajectory
        
        typedef grob::GridFunction<grob::splineD1D,grob::GridUniform<T>,std::vector<T>> 
            GridFunction_t;
        GridFunction_t theta_tau; /// function of theta of tau
    };
    
    struct hidden_reg_l{
        template <typename T>
        inline auto operator ()(T const & l) const{
            return l*l;
        }
    };
    struct unhidden_reg_l {
        template <typename T>
        inline auto operator ()(T const& _hid_l) const {
            return std::sqrt(_hid_l);
        }
    };
    
    template <typename Tr_Type>
    using Bin_e_grid_t = grob::GridUniform<Tr_Type>;
    
    template <typename Tr_Type>
    using Bin_l_grid_t = grob::GridUniform<Tr_Type>; 
    //grob::GridFunctional<Tr_Type,hidden_reg_l,unhidden_reg_l>;

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
            std::vector<traj_info<Tr_Type>>
        >;
    
    template <typename bin_traj_pool_t_inst_t>
    struct bin_traj_pool_value_type;

    template <typename Tr_Type>
    struct bin_traj_pool_value_type<
        bin_traj_pool_t<Tr_Type>>
    {
        typedef Tr_Type type;
    };

    template <typename bin_traj_pool_t_inst_t>
    using bin_traj_p_vt_t = 
        typename bin_traj_pool_value_type<
            bin_traj_pool_t_inst_t
        >::type;


    /// <summary>
    /// makes inbin grid function (with interpolator) of some
    /// quantity, depending on traj_info(e,l) and, optionally (e,undim l)
    /// </summary>
    /// <param name="TrajPool">trajectory infos of particular bin</param>
    /// <param name="traj_info_transformator">sum function of trajectory info (and e, l)</param>
    /// <returns>GridFunction, which interpolates traj_info_transformator function</returns>
    template <typename Tr_Type,typename FuncType>
    auto bin_traj_pool_func(
        bin_traj_pool_t<Tr_Type> const& TrajPool,
        FuncType &&traj_info_transformator
        )
    {
        using Lin2Interpol = grob::interpolator_product<
            grob::linear_interpolator,
            grob::linear_interpolator
        > ;

        return grob::make_function_ref<Lin2Interpol>(
            TrajPool.Grid,
            grob::as_container([&](size_t i) {
                if constexpr (
                    std::is_invocable_v<
                        FuncType,
                        traj_info<Tr_Type>,
                        Tr_Type,
                        Tr_Type
                    >) 
                {
                    auto [e, l] = TrajPool.Grid[TrajPool.Grid.FromLinear(i)];
                    return traj_info_transformator(TrajPool.Values[i], e, l);
                }
                else {
                    return traj_info_transformator(TrajPool.Values[i]);
                }
                
            }, TrajPool.size())
        );
    }

    template <typename T>
    inline constexpr T _id_func(T x) { return x; }

    template <typename T>
    inline constexpr T _quad_func(T x) { return x; }

#define DECLARE_TRAJ_POOL_PARAMFUNCTION(final_param_name,transform,init_param)\
        template <typename Tr_Type, typename LE_Func_t> \
        inline auto traj_pool_##final_param_name##_func(bin_traj_pool_t<Tr_Type> const& TrajPool,\
        LE_Func_t const& LE_i) { \
        using Lin2Interpol = grob::interpolator_product< \
            grob::linear_interpolator,\
            grob::linear_interpolator \
        >; \
        return grob::make_function_ref<Lin2Interpol>( \
            TrajPool.Grid, \
            grob::as_container([&](size_t i) { \
                return transform(TrajPool.Values[i].init_param); \
                }, TrajPool.size()) \
            ); \
    }
    DECLARE_TRAJ_POOL_PARAMFUNCTION(dul, _id_func, dul)
    DECLARE_TRAJ_POOL_PARAMFUNCTION(su, _id_func, su)
    DECLARE_TRAJ_POOL_PARAMFUNCTION(tinth, _id_func, T_in_theta)
    DECLARE_TRAJ_POOL_PARAMFUNCTION(umin, _quad_func, rmin)
    DECLARE_TRAJ_POOL_PARAMFUNCTION(umax, _quad_func, rmax)
    DECLARE_TRAJ_POOL_PARAMFUNCTION(u0z, _id_func, u_delta_0z)
    DECLARE_TRAJ_POOL_PARAMFUNCTION(u1z, _id_func, u_delta_1z)
    DECLARE_TRAJ_POOL_PARAMFUNCTION(dusq, _id_func, su_minus_2_q)

//#undef DECLARE_TRAJ_POOL_PARAMFUNCTION


#define DECLARE_TRAJ_POOL_PARAM_MEMBER(param)\
    using param##_f_t = decltype(traj_pool_##param##_func(TrajPool, LE_i));\
    param##_f_t param;

#define CONSTRUCT_TRAJ_POOL_PARAM_MEMBER(param)\
    param(traj_pool_##param##_func(TrajPool, LE_i))

    template <typename Tr_Type, typename LE_Func_t>
    struct bin_traj_pool_Tin_func {
        bin_traj_pool_t<Tr_Type> const& TrajPool;
        LE_Func_t const& LE_i;
        DECLARE_TRAJ_POOL_PARAM_MEMBER(dul);
        DECLARE_TRAJ_POOL_PARAM_MEMBER(su);
        DECLARE_TRAJ_POOL_PARAM_MEMBER(dusq);
        DECLARE_TRAJ_POOL_PARAM_MEMBER(tinth);
        DECLARE_TRAJ_POOL_PARAM_MEMBER(u0z);
        DECLARE_TRAJ_POOL_PARAM_MEMBER(u1z);


        Tr_Type _minus_F2_inv;

        
        inline bin_traj_pool_Tin_func(
            bin_traj_pool_t<Tr_Type> const& TrajPool,
            LE_Func_t const& LE_i,
            Tr_Type _F2 = -0.25
        ) :
            TrajPool(TrajPool), LE_i(LE_i),
            CONSTRUCT_TRAJ_POOL_PARAM_MEMBER(dul),
            CONSTRUCT_TRAJ_POOL_PARAM_MEMBER(su),
            CONSTRUCT_TRAJ_POOL_PARAM_MEMBER(dusq),
            CONSTRUCT_TRAJ_POOL_PARAM_MEMBER(tinth),
            CONSTRUCT_TRAJ_POOL_PARAM_MEMBER(u0z),
            CONSTRUCT_TRAJ_POOL_PARAM_MEMBER(u1z),
            _minus_F2_inv(-1/_F2) {}

        
        inline auto theta1(Tr_Type e, Tr_Type l_undim, Tr_Type Lme) const {
            using T = Tr_Type;
            auto L = l_undim * Lme;
            auto L2 = L * L;
            
            auto z = ((1 + e)-L2) * _minus_F2_inv;
            bool small_z = z < 1e-3;


            T pot_cos;
            auto q_e = (1 + 2 * e) * _minus_F2_inv / 2;
            if (q_e <=0) {
                if (z < 1e-6) {
                    return pi<T>;
                }
                auto sqr_qe = std::sqrt(q_e * q_e + 2 * z);
                auto sqr_safe = Lme * 
                    std::sqrt(2* _minus_F2_inv * (1 - l_undim * l_undim));
                auto delta1 = u1z(e, l_undim) * (2 * z / (sqr_safe - q_e));
                auto delta0 = u0z(e, l_undim) * (sqr_safe - q_e);
                

                if (q_e < -1. / 1024 || l_undim < 1-1e-2) {
                    auto x = delta1 / (delta0);
                    bool bad_delta = delta0 <= 0;
                    auto m_asin = bad_delta?
                        0 :
                        std::asin(2 * std::sqrt(x) / (1 + x));
                    return (bad_delta || x > 1 ? m_asin : pi<T> -m_asin);
                }
                else {
                    auto d_delta = dusq(e, l_undim) * q_e;
                    pot_cos = d_delta / (delta0+ delta1);
                }
            }
            else {
                auto sqr_qe = std::sqrt(q_e * q_e + 2 * z);
                auto delta1 = (q_e + sqr_qe);
                auto delta0 = u0z(e, l_undim) * (2 * z / (sqr_qe + q_e));
                if ( q_e > 1. / 1024 || l_undim < 1 - 1e-2) {
                    auto x = delta1 / (delta0);
                    bool bad_delta = delta0 <= 0;
                    auto m_asin = bad_delta?
                        0 :
                        std::asin(2 * std::sqrt(x) / (1 + x));
                    return ( bad_delta || x > 1 ? m_asin : pi<T> -m_asin);
                }
                else {
                    auto d_delta = dusq(e, l_undim) * q_e;
                    pot_cos = d_delta / (delta0 + delta1);
                }
            }
            
            /*else if (l_undim >(T)0.96 || small_z) {
                
                if (q_e < -(T)1. / 16) {
                    auto delta_u = dul(e, l_undim) * std::sqrt(1 - l_undim * l_undim);
                    auto delta_1 = (q_e + sqr_qe);
                    auto x = delta_1 / ( delta_u - delta_1);
                    return pi<T> -std::asin(2 * std::sqrt(x) / (1 + x));
                }
                pot_cos = q_e / std::sqrt(q_e * q_e + 2 * z);
            }
            else {
                auto delta_u = dul(e, l_undim) * std::sqrt(1 - l_undim * l_undim);
                auto sum_u = su(e, l_undim);
                auto u02 = (sum_u - delta_u) ;
                if (u02 > 2)
                    return pi<T>;
                pot_cos = (sum_u-2) / delta_u;
            }*/
            auto ThetaM = std::acos(
                pot_cos < -1 ? -1 : (
                    pot_cos <= 1 ? pot_cos : 1
                    )
            );
            return ThetaM;
        }
        inline auto operator()(Tr_Type e, Tr_Type l_undim, Tr_Type Lme) const {
            auto ThetaM = theta1(e, l_undim, Lme);
            return tinth(e, l_undim) * ThetaM;
        }
    };

    template <typename Tr_Type, typename LE_Func_t>
    struct bin_traj_pool_Tout_func {
        bin_traj_pool_t<Tr_Type> const& TrajPool;
        LE_Func_t const& LE_i;

        inline bin_traj_pool_Tout_func(
            bin_traj_pool_t<Tr_Type> const& TrajPool,
            LE_Func_t const& LE_i
        ) :
            TrajPool(TrajPool), LE_i(LE_i) {}

        template <typename T>
        inline auto operator()(T e, T l_undim) const {
            T L = l_undim * LE_i(-e);
            return eT(-e, L * L); 
        }
    };
    template <typename T,typename LEF_t,typename U>
    inline bin_traj_pool_Tin_func<T, LEF_t> 
        make_bin_traj_pool_Tin_func(
            bin_traj_pool_t<T> const& TrajPool,
            LEF_t const& LE_i,U _F2 = -0.25) 
    {
        return bin_traj_pool_Tin_func<T, LEF_t>(TrajPool, LE_i, _F2);
    }
    template <typename T, typename LEF_t>
    inline bin_traj_pool_Tout_func<T, LEF_t>
        make_bin_traj_pool_Tout_func(
            bin_traj_pool_t<T> const& TrajPool,
            LEF_t const& LE_i)
    {
        return bin_traj_pool_Tout_func<T, LEF_t>(TrajPool, LE_i);
    }


    struct T_Full_Id_t {
        template <typename T>
        T operator()(T x) {
            return x;
        }
    };

    /// <summary>
    /// class of transform function of T
    /// </summary>
    template <typename LE_Func_t,typename T_Full_Func_t>
    struct traj_info_transform_T_full {
        LE_Func_t const& LE_i;
        T_Full_Func_t TFunc;
        inline traj_info_transform_T_full(
            LE_Func_t const& LE_i,
            T_Full_Func_t TFunc
        ) : LE_i(LE_i), TFunc(std::move(TFunc)){}

        template <typename T,typename U>
        inline T operator ()(const traj_info<T>& TI, U e, U l) const{
            auto T_in = TI.T_in_theta*TI.theta_max;
            U L = l * LE_i(-e);
            auto T_out = eT(-e, L * L);
            return TFunc(T_in, T_out);
        }
    };

    /// <summary>
    /// enum class, what to extract in traj_info_transform_t
    /// </summary>
    enum class tit_ext_v {
        rmin, rmax, theta
    };

    template <tit_ext_v ext_val>
    struct traj_info_transform_t {

        template <typename T>
        inline T operator ()(const traj_info<T>& TI) const {
            if constexpr (ext_val == tit_ext_v::rmin)
                return TI.rmin;
            else if (ext_val == tit_ext_v::rmax)
                return TI.rmax;
            else
                return TI.theta_max;
        }
    };


    /// <summary>
    /// return tansform for makeing f(tin,tout) function in bin traj pool
    /// </summary>
    /// <typeparam name="T">type of e and l</typeparam>
    /// <param name="LE_i">LE function L(E)</param>
    /// <param name="TFunc">Function of Tin and Tout</param>
    /// <returns>transformator object</returns>
    template <typename LE_Func_t,typename T_full_func_t = T_Full_Id_t>
    traj_info_transform_T_full<LE_Func_t,T_full_func_t>
        make_traj_func(
            LE_Func_t  const & LE_i,
            T_full_func_t TFunc = T_Full_Id_t{}
        )
    {
        return traj_info_transform_T_full( LE_i, std::move(TFunc));
    }
    template <typename LE_Func_t>
    auto make_traj_tin_getter(LE_Func_t const& LE_i) {
        return make_traj_func(LE_i, [](auto t_in, auto t_out) {return t_in; });
    }
    template <typename LE_Func_t>
    auto make_traj_tout_getter(LE_Func_t const& LE_i) {
        return make_traj_func(LE_i, [](auto t_in, auto t_out) {return t_out; });
    }
    template <typename LE_Func_t>
    auto make_traj_tfull_getter(LE_Func_t const& LE_i) {
        return make_traj_func(LE_i, [](auto t_in, auto t_out) {return t_in+t_out; });
    }


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
        auto Bin_e_grid = grob::GridUniform<Tr_Type>(
            std::get<0>(dEdL).left,std::get<0>(dEdL).right,2);
        auto Bin_l_p_grid = Bin_l_grid_t<Tr_Type>(
            std::get<1>(dEdL).left,
            std::get<1>(dEdL).right, 2//,
            //hidden_reg_l, unhidden_reg_l{}
        );
        typedef traj_info<Tr_Type> traj_info_t;
        auto get_traj = [&LE,&B,traj_bins](U e, U l)->traj_info_t {
            auto Lmax = LE.i_lm(-e);
            auto L = (Lmax * l);
            auto L2 = L * L;
            auto [rm,rp] = B.find_rmin_rmax(-e,L2,LE);
            auto m_traj = B.get_internal_traj<Tr_Type>(rm,rp,traj_bins);
            auto up = rp * rp;
            auto um = rm * rm;
            auto du = up - um;
            auto us = up + um;
            auto F2 = B._DD_F_C(1);
            
            auto z = (L2 - (1 + e)) / F2;
            auto q_e = (1 + 2 * e) / (-2 * F2);

            auto z_qe_ratio = std::abs(2 * z) / (q_e * q_e);

            auto sqr_e = std::sqrt(q_e * q_e + 2 * z);
            auto sqr_safe = (
                q_e < 0 ? 
                Lmax *std::sqrt(-2/ F2 * (1-l*l)) :
                sqr_e
            );
            auto dZp = (q_e < 0 ? (2 * z) / (sqr_safe - q_e) : sqr_safe + q_e);
            auto dZm = (q_e < 0 ? sqr_safe - q_e : (2 * z) / (sqr_safe + q_e));
            auto u_delta_0z = (
                (l > 1 - 1e-4) ?
                1 :
                (1 - um) / dZm
            );

            auto u_delta_1z = (
                (z_qe_ratio < 1e-3 || l > 1 - 1e-4) ?
                1 :
                (up - 1) / dZp
            );

            auto dus = us - 2;
            auto dus_q = (std::abs(q_e) > 1.0/32 ? dus / q_e : 2);
            if (dus < 0)
                dus_q = 2;
            traj_info_t Ret = {
                rm,rp,
                (l > 1-1e-4 ? 2* Lmax*std::sqrt(2/ -F2) : du/std::sqrt(1-l*l)),
                us,dus_q,
                u_delta_0z,u_delta_1z,
                std::get<1>(m_traj),
                std::get<0>(m_traj),
                InverseTrajectory<
                    typename traj_info<Tr_Type>::GridFunction_t::interpolator_t
                >(std::get<2>(m_traj),traj_bins)
            };
            return Ret;
        };
        if (is_static){
            auto mGrid = grob::mesh_grids(Bin_e_grid,Bin_l_p_grid);
            std::vector<traj_info_t> ti_values = {
                get_traj(std::get<0>(dEdL).left,std::get<1>(dEdL).left),
                get_traj(std::get<0>(dEdL).left,std::get<1>(dEdL).right),
                get_traj(std::get<0>(dEdL).right,std::get<1>(dEdL).left),
                get_traj(std::get<0>(dEdL).right,std::get<1>(dEdL).right)
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
            //bool e_1 = dEdL.template x<0>().left > -0.9;

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
                //Bin_l_p_grid.hidden() = Bin_l_grid_tmp.hidden();
                Bin_l_p_grid = Bin_l_grid_tmp;
                TI = std::move(TI_tmp);
            };
            
            for(size_t i=0;i< n_e_max*n_l_max;++i){
                auto [eps_e,eps_l] = Calculate_dT();
                
                bool eps_e_incr = (eps_e > T_error && Bin_e_grid.size() * 2 - 1 <= n_e_max);
                bool eps_l_incr = (eps_l > T_error && Bin_l_p_grid.size() * 2 - 1 <= n_l_max);

                if( !eps_e_incr && !eps_l_incr)
                {
                    break;   
                }
                if(eps_e_incr){
                    Increase_NE();
                }
                if(eps_l_incr){
                    Increase_NL();
                }
               
            }
            //if (Bin_e_grid.size() * Bin_l_p_grid.size() != TI.size()) {
            //    throw std::runtime_error("error in traj pool size");
            //}
            auto mGrid = grob::mesh_grids(Bin_e_grid,Bin_l_p_grid);
            return grob::make_grid_object(mGrid,std::move(TI));
        }
    }
};

#endif//TRAJ_POOL_HPP