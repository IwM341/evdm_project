#ifndef MEASURE_HPP
#define MEASURE_HPP

#include <grob/grid_objects.hpp>
#include "utils/prng.hpp"
#include "utils/math_external.hpp"
namespace evdm{

    /**
    * the full square measure, where l = L/Lmax(e)
    */
    struct measure_dEdl {};

    /**
    * measure, \int{dEdL}
    */
    struct measure_dEdL {};

    /**
    * measure \int{dEdL^2}
    */
    struct measure_dEdL2 {};
    
    template <typename T>
    using bin_dedl_t = grob::Point<grob::Rect<T>, grob::Rect<T>>;

    template <typename T,typename U>
    inline T dEdl(bin_dedl_t<T> const & m_bin,
                    U const& Lm0, U const& Lm1){
        return m_bin.volume();
    }

    template <typename T,typename U>
    inline T dEdL(bin_dedl_t<T> const & m_bin,
                    U const & Lm0,U const & Lm1)
    {
        return m_bin.volume()*(Lm0+Lm1)/2;
    }
    template <typename T,typename U>
    inline T dEdL2(bin_dedl_t<T> const & m_bin,
                    U const & Lm0,U const & Lm1){
        auto [l0, l1] = std::get<1>(m_bin);
        return std::get<0>(m_bin).volume()*(l1*l1-l0*l0)*(Lm0*Lm0+ Lm0*Lm1 + Lm1*Lm1)/3;
    }
    template <typename T>
    inline T dEd1_up(bin_dedl_t<T> const& m_bin,
        T const& Lm0, T const& Lm1) {
        return m_bin.volume();
    }

    template <typename T>
    inline T dEdL_up(bin_dedl_t<T> const& m_bin,
        T const& Lm0, T const& Lm1)
    {
        return m_bin.volume() * 2*(Lm1* Lm1* Lm1- Lm0* Lm0* Lm0)/3;
    }
    template <typename T>
    inline T dEdL2_up(bin_dedl_t<T> const& m_bin,
        T const& Lm0, T const& Lm1) {
        return m_bin.volume() * std::get<1>(m_bin).center() * (Lm0 * Lm0 + Lm1 *Lm1) / 2;
    }

    /// <summary>
    /// measure of el bin assuming e < -1/2
    /// </summary>
    /// <param name="mes_t">value of measure indicator measure_dEd*</param>
    /// <param name="m_bin">bin in el grid</param>
    /// <param name="Lm0">Lmax(e left)</param>
    /// <param name="Lm1">Lmax(e right)</param>
    /// <returns>real measure</returns>
    template <typename T,typename U,typename measure_t>
    inline T mes_down(measure_t mes_t,bin_dedl_t<T> const& m_bin, U const& Lm0, U const& Lm1) {
        if constexpr (std::is_same_v<measure_t, measure_dEdl>) {
            return dEdl(m_bin, Lm0, Lm1);
        }
        else if (std::is_same_v < measure_t, measure_dEdL>) {
            return dEdL(m_bin, Lm0, Lm1);
        }
        else if (std::is_same_v < measure_t, measure_dEdL2>) {
            return dEdL2(m_bin, Lm0, Lm1);
        }
        else {
            static_assert(true,"unexpected mesure_t type");
        }
    }

    /// <summary>
    /// measure of el bin assuming e < -1/2
    /// </summary>
    /// <param name="mes_t">value of measure indicator type measure_dEd*</param>
    /// <param name="LE">L(E) function</param>
    /// <returns>measure function of bin</returns>
    template <typename measure_t,typename LE_Functype>
    inline auto get_bin_mes(measure_t mes_t,LE_Functype && LE){
        return [&LE](auto const &m_bin){
            typedef decltype(m_bin.volume()) T;
            T Lm0 = (T)LE(-std::get<0>(m_bin).left);
            T Lm1 = (T)LE(-std::get<0>(m_bin).right);
            return mes_down(measure_t{},m_bin,Lm0,Lm1);
        };
    }
    /// <summary>
    /// measure of el bin assuming e > -1/2
    /// </summary>
    /// <param name="mes_t">value of measure indicator type</param>
    /// <param name="m_bin">bin in el grid</param>
    /// <param name="Lm0">Lmax(e left)</param>
    /// <param name="Lm1">Lmax(e right)</param>
    /// <returns>real measure</returns>
    template <typename T, typename measure_t>
    inline T mes_EL(measure_t mes_t, bin_dedl_t<T> const& m_bin, T const& Lm0, T const& Lm1) {
        if constexpr (std::is_same_v<measure_t, measure_dEdl>) {
            return dEdl_up(m_bin, Lm0, Lm1);
        }
        else if (std::is_same_v < measure_t, measure_dEdL>) {
            return dEdL_up(m_bin, Lm0, Lm1);
        }
        else if (std::is_same_v < measure_t, measure_dEdL2>) {
            return dEdL2_up(m_bin, Lm0, Lm1);
        }
        else {
            static_assert(true, "unexpected mesure_t type");
        }
    }


    /// <summary>
    /// struct, which indicates no addtional boundaries in EL generation
    /// </summary>
    struct _empt_bounds_{};

    /// @brief uniform generatin of (E,L) in bin for dEd1 measure 
    /// @param m_bin 
    /// @param b 
    /// @param G 
    /// @param LE 
    /// @return bin generetor lambda: WARNING LE function captured by reference
    template <
        typename T,typename U,
        typename GenType,typename LE_FuncType>
    inline auto gen_EL_dEdl(bin_dedl_t<T> const & m_bin,
                                U const & b,
                                GenType && G,LE_FuncType && LE)
    {   
        return [G,b,m_bin,&LE](){
            auto xi_e = G();
            auto e = std::get<0>(m_bin).reduction(xi_e);
            return grob::Point<T,T>{e,
                LE(-e)*std::get<1>(m_bin).reduction(G())};
        };
    }

    template <
        typename T, typename U,
        typename GenType, 
        typename LE_FuncType,
        typename LE_Bound_t>
    inline auto gen_EL_dEdl_up_bound(bin_dedl_t<T> const& m_bin,
        U const& b,
        GenType&& G, LE_FuncType&& LE, LE_Bound_t &&Bound)
    {
        return [G, b, m_bin, &LE,&Bound]() {
            auto e = std::get<0>(m_bin).reduction(G());
            
            auto LE_tmp = LE(-e);
            T Lmax_ordinary = std::get<1>(m_bin).right * LE_tmp;
            T Lmin_ordinary = std::get<1>(m_bin).left * LE_tmp;
            T mLmax = std::min(Lmax_ordinary, (T)Bound(-e));

            return MCResult<grob::Point<T, T>, T>{ 
                { e, Lmin_ordinary + (mLmax- Lmin_ordinary)*G() },
                (mLmax - Lmin_ordinary) / (Lmax_ordinary - Lmin_ordinary)
            };
        };
    }


    /// @brief uniform generatin of (E,L) in bin for dEdL measure 
    /// @param m_bin bin [e0,e1]*[l0,l1]
    /// @param b (L1-L0)/(L1+L0)
    /// @param G random generator
    /// @param LE L(e)
    /// @return bin generetor lambda: WARNING LE function captured by reference
    template <typename T, typename U, typename GenType,typename LE_FuncType>
    inline auto gen_EL_dEdL(bin_dedl_t<T> const & m_bin,
                                U const & b,
                                GenType && G,LE_FuncType && LE)
    {   
        auto b1 = (1-b)/2;
        return [m_bin,&LE,G,b,b1](){
            auto xi = G();
            auto ue =(b >= 1 ? std::sqrt(xi) : xi/(b1 + std::sqrt(b1*b1+b*xi)));
            auto e = std::get<0>(m_bin).reduction(ue);
            return grob::Point<T,T> {e,
                    std::get<1>(m_bin).reduction(G())*LE(-e)};
        };
    }

    template <typename T, typename U, 
        typename GenType, typename LE_FuncType,
        typename LE_Bound_t>
    inline auto gen_EL_dEdL_up_bound(bin_dedl_t<T> const& m_bin,
        U const& b,
        GenType&& G, LE_FuncType&& LE, LE_Bound_t && Bound)
    {
        auto b1 = (1 - b) / 2;
        return [m_bin, &LE, &Bound,G, b, b1]() {
            auto xi = G();
            auto ue = (b >= 1 ? std::sqrt(xi) : xi / (b1 + std::sqrt(b1 * b1 + b * xi)));
            auto e = std::get<0>(m_bin).reduction(ue);

            auto LE_tmp = LE(-e);
            T Lmax_ordinary = std::get<1>(m_bin).right * LE_tmp;
            T Lmin_ordinary = std::get<1>(m_bin).left * LE_tmp;
            T mLmax = std::min(Lmax_ordinary, (T)Bound(-e));

            return MCResult<grob::Point<T, T>, T>{ 
                {e, Lmin_ordinary + (mLmax - Lmin_ordinary) * G() },
                (mLmax - Lmin_ordinary) / (Lmax_ordinary - Lmin_ordinary)
            };
        };
    }

    /// @brief uniform generatin of (E,L) in bin for dEdL2 measure 
    /// @param m_bin bin [e0,e1]*[l0,l1]
    /// @param b (L1-L0)/(L1+L0)
    /// @param G random generator
    /// @param LE L(e)
    /// @return 
    template <typename T,typename U,typename GenType,typename LE_FuncType>
    inline auto gen_EL_dEdL2(
        bin_dedl_t<T> const &  m_bin,
        U const & b,
        GenType && G,LE_FuncType && LE)
    {
       
        auto b1 = (1-b); 
        auto b12 = b1 * b1; //(1-b)^2
        auto b13 = b1 * b12; // (1-b)^3

        auto b23 = (3 + b * b);
        auto l0 = std::get<1>(m_bin).left;
        auto l1 = std::get<1>(m_bin).right;
        auto l02 =l0*l0;
        auto l12 =l1*l1; 
        


        return [G,m_bin,l02,l12,b1, b12, b13, b23, &LE](){
            auto xi = G();
            auto _sqr = std::cbrt(b13 + 2 * xi * b23);
            auto ue = xi* b23/(b12+b1* _sqr+ _sqr* _sqr);

            auto e = std::get<0>(m_bin).reduction(ue);
            auto lx = G();
            return grob::Point<T,T>{e,
                        std::sqrt(l02*(1-lx)+l12*lx)*LE(-e)
                    };
        };
    }

    template <typename T, typename U, 
        typename GenType, typename LE_FuncType,
        typename LE_Bound_t >
    inline auto gen_EL_dEdL2_up_bound(
        bin_dedl_t<T> const& m_bin,
        U const& b,
        GenType&& G, LE_FuncType&& LE, LE_Bound_t && Bound)
    {

        auto b1 = (1 - b);
        auto b12 = b1 * b1; //(1-b)^2
        auto b13 = b1 * b12; // (1-b)^3

        auto b23 = (3 + b * b);
        auto l0 = std::get<1>(m_bin).left;
        auto l1 = std::get<1>(m_bin).right;
        auto l02 = l0 * l0;
        auto l12 = l1 * l1;



        return [G, m_bin, l02, l12, b1, b12, b13, b23, &LE,&Bound]() {
            auto xi = G();
            auto _sqr = std::cbrt(b13 + 2 * xi * b23);
            auto ue = xi * b23 / (b12 + b1 * _sqr + _sqr * _sqr);

            auto e = std::get<0>(m_bin).reduction(ue);
            
            auto LE_tmp = LE(-e);
            T Lmax_ordinary = l02 * LE_tmp;
            T Lmin_ordinary = l12 * LE_tmp;
            T mLmax = std::min(Lmax_ordinary, (T)Bound(-e));

            T Lmax2 = Lmax_ordinary * Lmax_ordinary;
            T Lmin2 = Lmin_ordinary * Lmin_ordinary;
            T Lm2 = mLmax * mLmax;


            return MCResult<grob::Point<T, T>, T>{ 
                {e, std::sqrt(Lmin2 + (Lm2 - Lmin2) *G() )},
                    (Lm2 - Lmin2)/(Lmax2-Lmin2)
            };
        };
    }

    /// <summary>
    /// returns generator function of bin, according to measure mes_t;
    /// </summary>
    /// if auto f = gen_EL(...) then f() = {e,l}, distributed by measure mes_t
    /// <param name="mes_t">value of measure indicator type</param>
    /// <param name="m_bin">bin [e0,e1]*[l0,l1]</param>
    /// <param name="b">bin in el grid</param>
    /// <param name="G">generator from 0 to 1, captured by value</param>
    /// <param name="LE">LE L(e)</param>
    /// <returns>bin generetor lambda: WARNING LE function captured by reference</returns>
    template <
        typename measure_t, typename T, 
        typename GenType, typename LE_FuncType,
        typename LE_Bound_t = _empt_bounds_>
    inline auto gen_EL(
        measure_t mes_t,bin_dedl_t<T> const& m_bin, 
        GenType&& G, LE_FuncType&& LE, LE_Bound_t&& Bound = _empt_bounds_{}) {
        auto Lm0 = LE(-std::get<0>(m_bin).left);
        auto Lm1 = LE(-std::get<0>(m_bin).right);
        auto b = (Lm1 - Lm0) / (Lm1 + Lm0);

        if constexpr (std::is_same_v<LE_Bound_t, _empt_bounds_>) {
            if constexpr (std::is_same_v<measure_t, measure_dEdl>) {
                return gen_EL_dEdl(m_bin, b, G, LE);
            }
            else if constexpr (std::is_same_v < measure_t, measure_dEdL>) {
                return gen_EL_dEdL(m_bin, b, G, LE);
            }
            else if constexpr (std::is_same_v < measure_t, measure_dEdL2>) {
                return gen_EL_dEdL2(m_bin, b, G, LE);
            }
            else {
                static_assert(true, "unexpected mesure_t type");
            }
        }
        else {
            if constexpr (std::is_same_v<measure_t, measure_dEdl>) {
                return gen_EL_dEdl_up_bound(m_bin, b, G, LE,Bound);
            }
            else if constexpr (std::is_same_v < measure_t, measure_dEdL>) {
                return gen_EL_dEdL_up_bound(m_bin, b, G, LE, Bound);
            }
            else if constexpr (std::is_same_v < measure_t, measure_dEdL2>) {
                return gen_EL_dEdL2_up_bound(m_bin, b, G, LE, Bound);
            }
            else {
                static_assert(true, "unexpected mesure_t type");
            }
        }
    }

    /// <summary>
    /// creating generator function of e,l, with measure d^3v.
    /// </summary>
    /// if f = mc_d3v(...) then f() is MCResult(Point(e,l),weight), where weight, corresponding to d^3v measure
    /// therefore the average of weight is integral d^3v over el bin
    /// <param name="m_bin">bin [e0,e1]*[l0,l1]</param>
    /// <param name="b">bin in el grid. NOTE: the bin should be cut in e, so that l could be from l0 to lmax(r,e)</param>
    /// <param name="G">generator from 0 to 1, captured by value</param>
    /// <param name="LE">LE L(e)</param>
    /// <param name="is_Bound">any type in case of considering kinematic restriction</param>
    /// <returns>mk result, grob::point<T,T,T>(e,L/Le,Le) </returns>
    template < typename U,typename T, typename GenType, typename LE_FuncType>
    inline auto mc_d3v(U r, U Phi_r, bin_dedl_t<T> const& m_bin, GenType&& G, LE_FuncType&& LE) {
            
        return [m_bin, r, Phi_r, G, &LE]() {
            auto xi = G();
            T e = std::get<0>(m_bin).reduction(xi);
            auto W_e = std::get<0>(m_bin).volume();
            auto LE_max = LE(-e);

            auto Vmax = std::abs(std::sqrt(Phi_r + e));
            //if (std::isnan(Vmax)) {
            //    throw std::runtime_error("v is nan!!!");
            //  }
            auto Lmax = r * Vmax;

            auto Lmax_undim = (LE_max > 0 ? Lmax / LE_max : 1);

            //if Lmax_undim != 0, then _d0 is just l0*Le/Lmax <=1
            //if Lmax_undim == 0, then _d0 = 0, _d1 = 1
            auto [_l0, _l1] = std::get<1>(m_bin);

            auto _d0 = upbound(downbound(_l0 / Lmax_undim, 0), 1);
            auto _d1 = upbound(_l1 / Lmax_undim, 1);

            auto _sqr_l0 = std::sqrt(1 - _d0 * _d0);
            auto _sqr_l1 = std::sqrt(1 - _d1 * _d1);
            auto d_sqr = ((_d1 < 0.5 || _d0 < 0.5) ?
                (_d1 - _d0) * (_d1 + _d0) / (_sqr_l1 + _sqr_l0) : 
                _sqr_l0 - _sqr_l1);

            auto _L_sqr = _sqr_l1 + (_sqr_l0 - _sqr_l1) * G();
            auto W_L = d_sqr * Vmax;

            T l = Lmax_undim * std::sqrt(1- _L_sqr* _L_sqr);
            return MCResult<grob::Point<T, T,T>, T>(
                grob::Point<T, T,T>(e, l, LE_max),
                W_e* W_L);
        };
    }

}

#endif//MEASURE_HPP