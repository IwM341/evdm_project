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
        T umin_l2, umax;
        T dul, su, su_minus_2_q;/// (umax-umin)/sqrt(1-l^2), umax+umin, (umax+umin-2)/q_e
        T u_delta_0z, u_delta_1z;
        T T_in_theta; /// regularized by theta period
        T theta_max;  /// max theta parametr of trajectory
        

        typedef grob::GridFunction<grob::splineD1D,grob::GridUniform<T>,std::vector<T>> 
            GridFunction_t;
        GridFunction_t theta_tau; /// function of theta of tau

        template <typename...Args>
        inline static traj_info constructor(Args&&...args) {
            return traj_info{ std::forward<Args>(args)... };
        }

        SERIALIZATOR_FUNCTION(
            PROPERTY_NAMES("rmin","rmax", "umin_l2", "umax",
                "dul", "su", "su_minus_2_q",
                "u_delta_0z", "u_delta_1z",
                "T_in_theta","theta_max",
                "theta_tau"
            ),
            PROPERTIES(rmin, rmax, umin_l2, umax,
                dul, su, su_minus_2_q, 
                u_delta_0z,  u_delta_1z,
                T_in_theta,theta_max, theta_tau
            )
        )
        DESERIALIZATOR_FUNCTION(constructor,
            PROPERTY_NAMES("rmin", "rmax", "umin_l2", "umax",
                "dul", "su", "su_minus_2_q",
                "u_delta_0z", "u_delta_1z",
                "T_in_theta", "theta_max",
                "theta_tau"
            ),
            PROPERTY_TYPES(rmin, rmax, umin_l2, umax,
                dul, su, su_minus_2_q,
                u_delta_0z, u_delta_1z,
                T_in_theta, theta_max,
                theta_tau
            )
        )
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


    namespace trajpool_param_arrays{
    
        #define DECLARE_TRAJ_POOL_PARAM_STRUCT(param,inner_param) \
        template <typename Array_t> \
        struct g_##param{ \
            Array_t const& traj_pools;\
            inline g_##param(Array_t const& traj_pools):traj_pools(traj_pools){}\
            inline size_t size()const {\
                return traj_pools.size();\
            }\
            inline auto operator [](size_t i) const{\
                return traj_pools[i].inner_param;\
            }\
        };
        DECLARE_TRAJ_POOL_PARAM_STRUCT(min_l2, min_l2);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(dul, dul);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(su, su);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(tinth, T_in_theta);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(umin_l2, umin_l2);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(umax, umax);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(u0z, u_delta_0z);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(u1z, u_delta_1z);
        DECLARE_TRAJ_POOL_PARAM_STRUCT(dusq, su_minus_2_q);
        #undef DECLARE_TRAJ_POOL_PARAM_STRUCT
    };

    template <typename Tr_Type, typename LE_Func_t>
    struct bin_traj_pool_Tin_func {
        bin_traj_pool_t<Tr_Type> const& TrajPool;
        LE_Func_t const& LE_i;

        Tr_Type _minus_F2_inv;

        
        inline bin_traj_pool_Tin_func(
            bin_traj_pool_t<Tr_Type> const& TrajPool,
            LE_Func_t const& LE_i,
            Tr_Type _F2 = -0.25
        ) :
            TrajPool(TrajPool), LE_i(LE_i),
            _minus_F2_inv(-1/_F2) {}

        using Lin2Interpol = grob::interpolator_product< 
            grob::linear_interpolator, 
            grob::linear_interpolator 
        >;
        #define DECLARE_PARAM_FUNC(param) \
            auto param = grob::make_function_ref<Lin2Interpol>\
            (TrajPool.Grid,trajpool_param_arrays::g_##param{TrajPool.Values});
        #define DECLARE_ALL_FUNCS\
            DECLARE_PARAM_FUNC(umin_l2);\
            DECLARE_PARAM_FUNC(umax);\
            DECLARE_PARAM_FUNC(dul);\
            DECLARE_PARAM_FUNC(su);\
            DECLARE_PARAM_FUNC(dusq);\
            DECLARE_PARAM_FUNC(tinth);\
            DECLARE_PARAM_FUNC(u0z);\
            DECLARE_PARAM_FUNC(u1z);

        inline std::tuple<Tr_Type, Tr_Type, Tr_Type> 
            u0_u1_theta1(Tr_Type e, Tr_Type l_undim, Tr_Type Lme) const 
        {
            DECLARE_ALL_FUNCS;
            using T = Tr_Type;
            auto L = l_undim * Lme;
            auto L2 = L * L;
            auto l2 = l_undim * l_undim;
            if (L2 >= 1 + e) { //Left Zone
                if(l_undim < 0.8) { // Down Zone
                    auto u0 = umin_l2(e, l_undim)* l2/(1 + std::sqrt(1 - l2));
                    auto u1 = std::max(su(e, l_undim) - u0, u0);
                    return { u0,u1 ,pi<T> };
                } else {

                    auto u1_min_u0 = dul(e, l_undim) * std::sqrt(1 - l2);
                    auto u1_plus_u0 = su(e, l_undim);
                    auto u0 = downbound((u1_plus_u0 - u1_min_u0) / 2,0);
                    auto u1 = upbound((u1_plus_u0 + u1_min_u0) / 2, 1);  
                    return {u0,u1 ,pi<T>};
                }
            }
            else { // Right Zone
                auto z = ((1 + e) - L2) * _minus_F2_inv;
                auto q_e = (1 + 2 * e) * _minus_F2_inv / 2;
                auto sqr_qe = std::sqrt(q_e * q_e + 2 * z);
                auto delta1 = (q_e + sqr_qe);
                auto u1 = 1 + delta1;
                if (l_undim < 0.9) { // Down Zone
                    auto u0 = umin_l2(e, l_undim)* l2 / (1 + std::sqrt(1 - l2));
                    auto delta0 = 1 - u0;
                    auto x = delta1 / (delta0);
                    auto bad_delta = delta0 <= 0 || delta1 <= 0;
                    auto m_asin = (!bad_delta) ? std::asin(2 * std::sqrt(x) / (1 + x)) : 0;
                    auto RetTh = (bad_delta || x > 1 ? m_asin : pi<T> -m_asin);
                    return { u0 ,u1, RetTh }; 
                }
                else { // Up Zone
                    auto delta0 = u0z(e, l_undim) * (2 * z / (sqr_qe + q_e));
                    auto u0 = 1 - delta0;
                    auto x = delta1 / (delta0);
                    auto bad_delta = delta0 <= 0 || delta1 <= 0;
                    auto m_asin = (!bad_delta) ? std::asin(2 * std::sqrt(x) / (1 + x)) : 0;
                    auto RetTh = (bad_delta || x > 1 ? m_asin : pi<T> -m_asin);
                    return { u0 ,u1, RetTh };
                }
            }
        }
        inline Tr_Type
            theta1(Tr_Type e, Tr_Type l_undim, Tr_Type Lme) const
        {
            DECLARE_ALL_FUNCS;

            using T = Tr_Type;
            auto L = l_undim * Lme;
            auto L2 = L * L;

            if (L2 >= 1 + e) { //Left Zone
                return pi<T>;
            }
            else { // Right Zone

                auto z = ((1 + e) - L2) * _minus_F2_inv;
                auto q_e = (1 + 2 * e) * _minus_F2_inv / 2;
                auto sqr_qe = std::sqrt(q_e * q_e + 2 * z);
                auto delta1 = (q_e + sqr_qe);
                if (l_undim < 0.9) { // Down Zone
                    auto l2 = l_undim * l_undim;
                    auto u0 = umin_l2(e, l_undim)*l2 / (1 + std::sqrt(1 - l2));
                    auto delta0 = 1 - u0;
                    auto x = delta1 / (delta0);
                    auto bad_delta = delta0 <= 0;
                    auto m_asin = (!bad_delta) ? std::asin(2 * std::sqrt(x) / (1 + x)) : 0;
                    auto RetTh = (bad_delta || x > 1 ? m_asin : pi<T> -m_asin);
                    return RetTh;
                }
                else { // Up Zone
                    auto delta0 = u0z(e, l_undim) * (2 * z / (sqr_qe + q_e));
                    auto x = delta1 / (delta0);
                    auto bad_delta = delta0 <= 0;
                    auto m_asin = (!bad_delta) ? std::asin(2 * std::sqrt(x) / (1 + x)) : 0;
                    auto RetTh = (bad_delta || x > 1 ? m_asin : pi<T> -m_asin);
                    return RetTh;
                }
            }
        }
        
        inline auto operator()(Tr_Type e, Tr_Type l_undim, Tr_Type Lme) const {
            auto ThetaM = theta1(e, l_undim, Lme);
            DECLARE_PARAM_FUNC(tinth);
            return tinth(e, l_undim) * ThetaM;
        }
        inline auto tin_theta(Tr_Type e, Tr_Type l_undim) const {
            DECLARE_PARAM_FUNC(tinth);
            return tinth(e, l_undim);
        }
        inline auto tin_full(Tr_Type e, Tr_Type l_undim, Tr_Type ThetaM) const {
            DECLARE_PARAM_FUNC(tinth);
            return tinth(e, l_undim) * ThetaM;
        }
        #undef DECLARE_ALL_FUNCS
        #undef DECLARE_PARAM_FUNC
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

        inline auto operator()(Tr_Type e, Tr_Type l_undim) const {
            Tr_Type L = l_undim * LE_i(-e);
            return eT(-e, L * L); 
        }

        inline auto operator()(Tr_Type e, Tr_Type l_undim, Tr_Type  Lmax) const {
            Tr_Type L = l_undim * Lmax;
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
        auto e0_positive = B.Phi[0];
        auto C0 = -B._DD_F_r(0);//-F''(0)
        typedef traj_info<Tr_Type> traj_info_t;
        auto get_traj = [&LE,&B, e0_positive, C0,traj_bins](U e, U l)->traj_info_t {
            auto Lmax = LE.i_lm(-e);
            auto rm_e = LE.i_rm(-e);
            auto L = (Lmax * l);
            auto L2 = L * L;
            auto [rm,rp] = B.find_rmin_rmax(-e,L2,LE);
            auto m_traj = B.get_internal_traj<Tr_Type>(rm,rp,traj_bins);
            auto um = rm * rm;
            auto up = rp * rp;
            auto du = up - um;
            auto us = up + um;
            auto F2_1 = B._DD_F_C(1);
            auto F2_r = B._DD_F_r(rm_e);

            auto z = (L2 - (1 + e)) / F2_1;
            auto q_e = (1 + 2 * e) / (-2 * F2_1);

            auto z_qe_ratio = std::abs(2 * z) / (q_e * q_e);

            auto q_e_sqrt_expr = q_e * q_e + 2 * z;
            auto sqr_e = std::sqrt(q_e_sqrt_expr);
            auto sqr_safe = (
                q_e_sqrt_expr < 0 ? 
                0 :
                sqr_e
            );
            auto pot_sqr = Lmax * std::sqrt(-2 / F2_1 * (1 - l * l));
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

            auto _l2 = l * l;

            auto _de = (e0_positive + e);
            auto Lm_de_lim =( 
                _de < 1e-5 ?
                1 / std::sqrt(2 * C0) : 
                Lmax/_de
            );
            auto u0_l2_limit = Lm_de_lim * Lmax;
            auto _l2_s = _l2 / (1 + std::sqrt(1 - _l2));
            auto u0_l2 = (l > 1e-3 ? um / _l2_s : 2*u0_l2_limit);
            traj_info_t Ret = {
                rm,rp,
                u0_l2 ,up,
                (l > 1-1e-4 ? 2* Lmax*std::sqrt(2/ -F2_r) : du/std::sqrt(1-l*l)), //dul
                us,dus_q,
                u_delta_0z,u_delta_1z,
                std::get<1>(m_traj),
                std::get<0>(m_traj),
                InverseTrajectory<
                    typename traj_info<Tr_Type>::GridFunction_t::interpolator_t
                >(std::get<2>(m_traj),traj_bins)
            };
            return std::move(Ret);
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