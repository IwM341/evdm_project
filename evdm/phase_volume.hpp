#ifndef PHASE_VOLUME_HPP
#define PHASE_VOLUME_HPP
#include <cmath>
#include <utility>
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
            T c = Tailor_o_Solution(M_PI)*(1/M_PI);
            T x0 = z + (Tailor_o_Solution(z)-c*z)+u*z*(M_PI-z)/8;
            auto Step = [y,z](T x){return x - (x - y*sin(x) - z)/(1 - y*cos(x));};
            return (Step(Step(x0)));
        }
        template <typename T>
        T xsin(T y,T x){
            return x - y*sin(x);
        }

        template <typename T>
        struct alignas(alignof(T)*4) Array4{
            T data[4];
            inline T & operator [](size_t i){
                return data[i];
            }
            inline T const & operator [](size_t i)const{
                return data[i];
            }
            inline T accum()const {
                return data[0]+data[1]+data[2]+data[3];
            }
        };
        /// @brief eval b spline
        /// @param x_h parametr, equal to (x-x_i)/h
        /// @return F(x) so that F(h) = f0, F(2h) = f1, F'(h) = (f2-f0)/h, F'(2h) = (f3-f1)/h
        template <typename T>
        std::pair<T,T>   eval_B_spline(T x_h,T f0,T f1,T f2,T f3){
            Array4<T> X = {x_h,x_h,x_h,x_h};
            Array4<T> F = {f0,f1,f2,f3};
            Array4<T> C0 = {2.,-3.,3.,-1.};
            Array4<T> C1 = {-4.,19.0/2,-8.,5.0/2};
            Array4<T> C2 = {5.0/2,-7.,13./2,-2.};
            Array4<T> C3 = {-1./2,3./2,-3./2,1./2};

            Array4<T> Sum = C3;
            
            Array4<T> DSum = C3;
            for(size_t i=0;i<4;++i){
                DSum[i] = C3[i]*3;    
            }
            
            for(size_t i=0;i<4;++i){
                Sum[i] =Sum[i]*X[i] + C2[i];
            }
            for(size_t i=0;i<4;++i){
                DSum[i] = DSum[i]*X[i] + C2[i]*2;
            }
            for(size_t i=0;i<4;++i){
                Sum[i] =Sum[i]*X[i] + C1[i];
            }
            for(size_t i=0;i<4;++i){
                DSum[i] = DSum[i]*X[i] + C1[i];
            }
            for(size_t i=0;i<4;++i){
                Sum[i] =Sum[i]*X[i] + C0[i];
            }
            for(size_t i=0;i<4;++i){
                Sum[i]*=F[i];
            }
            for(size_t i=0;i<4;++i){
                DSum[i] *=F[i];
            }
            return {Sum.accum(),DSum.accum()};

        }
    };

};

#endif//PHASE_VOLUME_HPP