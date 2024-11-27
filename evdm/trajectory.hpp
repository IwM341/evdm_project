#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include "body_potential.hpp"
#include <grob/grid_objects.hpp>
#include <cmath>
#include <numbers>
#include "utils/math_external.hpp"
namespace evdm{


    namespace __detail{
        
        /// @brief solve equation x-y*sin(x) = z, z \in [0,\pi] 
        /// @return 
        template <typename T>
        T solve_xsin(T y,T z){
            T u = 1 - y;
            auto Tailor_o_Solution = [u](T z1){
                T _sqrt_part = std::cbrt(std::sqrt(8*u*u*u+9*z1*z1)+3*z1);
                return (_sqrt_part*_sqrt_part - 2*u)/_sqrt_part;
            };
            T c = Tailor_o_Solution(pi<T>)*(1/pi<T>);
            T x0 = z + (Tailor_o_Solution(z)-c*z)+u*z*(pi<T>-z)/8;
            auto Step = [y,z](T x){return x - (x - y*sin(x) - z)/(1 - y*cos(x));};
            return (Step(Step(x0)));
        }
        template <typename T>
        T xsin(T y,T x){
            return x - y*sin(x);
        }
    }

    /// @brief external reqularized trajectory duration
    /// (i.e. T(e,l) - pi/(2*e^1.5) )
    /// @param e energy
    /// @param L2 squared momentum
    /// @return 
    template <typename T>
    inline T reT(T const & e,T const & L2){
        T s_sqre = ssqrt(1 - e - L2);
        if(e < 1e-2){
            auto z = 2*(s_sqre)/(1-2*e);
            auto z2 = z*z;
            return -z + z*z2*(
                [](auto x){
                    return ((T)1./6 + x*(
                        -(T)1./10 + x*(
                            (T)1./14 + x*(
                                -(T)1./18
                            )
                        )
                    ));
                }(e*z*z)
            );
        } else {
            auto _e = std::sqrt(e);
            auto sqr = 2*_e*(s_sqre);
            return (sqr+std::atan((1-2*e)/sqr) - pi<T>/2 )/(2*e*_e);
        }
    }
    /// @brief external trajectory duration
    /// @param e energy
    /// @param L2 squared momentum
    /// @return 
    template <typename T>
    inline T eT(T const & e,T const & L2){
        if (e >= 0.5 && L2 >= (1 - e) - 1e-7) {
            return 0;
        }
        else if (e <= 0) {
            return std::numeric_limits<T>::max();
        } else {
            return reT(e, L2) + pi<T> / (2 * e * std::sqrt(e));
        }
        
    }
    /// @brief full duration of trajectory in 1/r potential
    /// @param e energy
    /// @return 
    template <typename T>
    inline T fT(T const & e){
        return pi<T>()/(2*e*std::sqrt(e));
    }

    /// @brief T_extranal(e,l)/T_full(e)
    /// @param e energy
    /// @param L2 squared momentum
    /// @return 
    template <typename T>
    inline T e_frac_T(T const & e,T const & L2){
        return 1 +  reT(e,L2)*(2/pi<T>()*e*std::sqrt(e));
    }

    template <typename T>
    T exTrajCoord(T const & e,T const & L2,T const & tau){
        auto r_av = 1/(2*e);
        auto y = unnan(std::sqrt(1-4*e*L2));
        auto r_df = y*r_av;
        return r_av - r_df*std::cos(
            __detail::solve_xsin(y,pi<T>()*(1 - (1-tau)*e_frac_T(e,L2)))
        );
    }

    template <typename Interpol,typename TauFunc>
    auto InverseTrajectory(TauFunc const & F,size_t Ntau){
        typedef std::decay_t<decltype(F.Grid)> Grid_t;
        typedef typename Grid_t::value_type T;

        auto InversedFunction = grob::make_function
            <grob::spline1D>(
            grob::Grid1<decltype(F.Values),
                    grob::vector_array_grid_helper
                >(F.Values),
            F.Grid
        );
        //T fT = F.Values[F.size()-1]; = 1
        auto theta_tau = [&](auto const &tau)->T{
            return InversedFunction(tau);
        };
        
        return grob::make_function_f<Interpol>(
            Grid_t(0,1,Ntau+1),
            theta_tau
        );
    }
}

#endif//TRAJECTORY_HPP