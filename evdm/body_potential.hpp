#ifndef BODY_POTENTIAL_HPP
#define BODY_POTENTIAL_HPP

#include <grob/grid_objects.hpp>
#include "utils/integration.hpp"
#include <cmath>
#include <algorithm>
#include <numbers>
#include <limits>
namespace evdm{

    namespace __detail{
        template <typename Functype,typename T>
        auto GaussQuad(Functype && f,T const & x,T const & h){
            constexpr T Cl = (T) 0.2113248653;
            constexpr T Cr = (T)0.7886751347;
            return (f(x + h*Cl)+f(x + h*Cr))*h/2;
        }
        template <typename T>
        T NewtonStep(T const &x,T const &F,T const &D_F,T const &DD_F){
            int S = (D_F >=0 ? 1 : -1);
            return x - 2*F/(
                D_F + S*
                    sqrt(
                        std::max(D_F*D_F-2*F*DD_F,(T)0)
                    )
               );
        }
        template <typename T>
        T NewtonStep(T const &x,T const &F,T const &D_F){
            return x - F/D_F;
        }

    };


    template <typename T = float>
    struct Body{
        typedef T value_type;
        typedef grob::GridUniform<T> GridR;
        typedef std::vector<T> Values_t;
        typedef grob::GridFunction<
            grob::interpolator_spline1D,
            GridR,Values_t> RFunc2_t;
        typedef grob::GridFunction<
            grob::linear_interpolator,
            GridR,Values_t> RFunc1_t;

        RFunc1_t Rho;
        RFunc2_t Q;
        RFunc2_t Phi;
        
        RFunc2_t VescFunc; /// sqrt(phi(r))
        RFunc1_t Temp;


        size_t low_steps_num;
        T Rho_max;
        T Vesc;
        
        
        static value_type get_vtype(){
            static_assert("Body::get_vtype() is not callable");
        }
        SERIALIZATOR_FUNCTION(PROPERTY_NAMES("Q","rho","phi","Temp","vesc"), PROPERTIES(Q, Rho, Phi, Temp, Vesc))
        WRITE_FUNCTION(Q,Rho,Phi,Temp,Vesc)
        DESERIALIZATOR_FUNCTION(
            Body,PROPERTY_NAMES("Q","rho","phi","Temp","vesc"), PROPERTY_TYPES(Q, Rho, Phi, Temp,Vesc))
        READ_FUNCTION(Body,PROPERTY_TYPES(Q,Rho,Phi,Temp,Vesc))
    
        auto const& getVesc() const{
            return VescFunc;
        }
        auto const& getTemp() const {
            return Temp;
        }
        template <typename Array_t>
        void setTemp(Array_t&& values) {
            if (Temp.size() != values.size()) {
                throw std::range_error(
                    "error trying to set temperature of body model,"
                    " incorrect data size");
            }
            for (size_t i = 0; i < Temp.size(); ++i) {
                Temp[i] = values[i];
            }
        }
        auto M() const{
            return grob::make_function_ref(
                Rho.Grid,
                grob::as_container([&](size_t i) {
                        return std::pow(Rho.Grid[i], 3) * Q[i];
                    },Rho.Grid.size())
                );
        }
        //CHECKED, OK
        template <typename...ArgsToCreateRho>
        static Body FromRho(T VescMax,ArgsToCreateRho&&...args){
            using namespace __detail;
            RFunc1_t Rho(std::forward<ArgsToCreateRho>(args)...);
            auto Grid = Rho.Grid;
            T h = Grid.h();
            Body B{RFunc1_t(std::move(Rho)),
                    RFunc2_t(Grid,std::vector<T>(Grid.size(),0)),
                    RFunc2_t(Grid,std::vector<T>(Grid.size(),0)),
                    RFunc2_t(Grid,std::vector<T>(Grid.size(),0)),
                    RFunc2_t(Grid,std::vector<T>(Grid.size(),0))
                };
            B.low_steps_num = 2;
            B.Q[0] = B.Rho[0];
            B.Phi[0] = 0;
            auto cube = [](auto x){return x*x*x;};
            //constexpr float euler_points =  
            for(size_t i=1;i<Grid.size();++i){
                B.Q[i] = (
                    B.Q[i-1]*cube(Grid[i-1]) + 
                    GaussQuad([&](T x){
                        return 3*B.Rho(x)*x*x;
                    },Grid[i-1],h))/(cube(Grid[i])
                );
            }
            T fac = 1/B.Q.Values.back();
            for(size_t i=0;i<Grid.size();++i){
                B.Q[i] *= fac;
                B.Rho[i] *= fac; 
            }
            for(size_t i=1;i<Grid.size();++i){
                B.Phi[i] = (
                    B.Phi[i-1] - 
                    GaussQuad([&](T x){return B.Q(x)*x;},Grid[i-1],h)
                );
            }
            T phi_shift = 1 - B.Phi.Values.back();
            for(auto & phi : B.Phi.Values){
                phi += phi_shift;
            }
            B.Rho_max = B.Rho[B.Rho.size()-1];
            
            B.Vesc = VescMax;
            for (size_t i = 0; i < B.Phi.size(); ++i) {
                B.VescFunc[i] = std::sqrt(B.Phi[i]);
            }

            return B;
        }

        struct LE_func_t{
            RFunc2_t Internal_lm;
            RFunc2_t Internal_um;

            /// @brief Lmax^2(e>0, internal
            inline auto i_lmq(T const &e)const{
                return (e > 0.5 ? std::pow(Internal_lm(e),2) : 1-e);
            }
            /// @brief Lmax^2(e>0), external trajs
            inline auto lmq(T const &e)const{
                return (
                    e > 0.5 ? 
                    std::pow(Internal_lm(e),2) : 
                    (e > 0 ?
                        1 / (4 * e) : 
                        (std::numeric_limits<T>::max())
                    )
                );
            }
            /// @brief Lmax(e>0) internal
            inline auto i_lm(T const &e)const{
                return (e > 0.5 ? Internal_lm(e) : sqrt(1-e));;
            }
            /// @brief Lmax^2(e>0) external
            inline auto lm(T const &e)const{
                return (
                    e > 0.5 ? 
                    Internal_lm(e) : (
                        e > 0 ? 
                        std::sqrt(1/(4*e)) : 
                        std::sqrt((std::numeric_limits<T>::max()))
                        )
                    );
            }
            /// @brief emax(e>0) internal
            inline auto i_rm(T const &e)const{
                return (e > 0.5 ? sqrt(std::abs(Internal_um(e))) : 1);
            }
            /// @brief emax(e>0) external
            inline auto rm(T const &e)const{
                return (
                    e > 0.5 ? 
                    sqrt(std::abs(Internal_um(e))) : 
                    (
                        e > 0 ? 
                        1 / (2 * e) : 
                        (std::numeric_limits<T>::max())
                    )
                );
            }

            #define MAKE_FUNC_WRAPPER(func) inline auto func()const{ \
                return [&](T const &e){return func(e);};\
            } 
            MAKE_FUNC_WRAPPER(i_lmq)
            MAKE_FUNC_WRAPPER(lmq)
            MAKE_FUNC_WRAPPER(i_lm)
            MAKE_FUNC_WRAPPER(lm)
            MAKE_FUNC_WRAPPER(i_rm)
            MAKE_FUNC_WRAPPER(rm)

        };

        //NOT_CHECKED
        /// @brief makes objects, containing Lmax(e)
        /// @param Ne number of steps in uniform e grid
        /// @return 
        LE_func_t get_le(size_t Ne)const{
            GridR Grid(0.5,Phi[0],Ne+1);
            RFunc2_t lm(Grid,std::vector<T>(Grid.size()));
            RFunc2_t um(Grid,std::vector<T>(Grid.size()));
            for(size_t i=0;i<Grid.size();++i){
                auto p = maxL2(Grid[i]);
                lm[i] = std::sqrt(p.first);
                um[i] = p.second*p.second;
            }
            return LE_func_t{std::move(lm),std::move(um)};
        }

        /// @brief Lmax_squared(e) and rm(e) assuming trajectory intersect body
        /// @param e positive dimentionless energy
        /// @return pair(Lmax^2(e),rm(e)), where rm(e) --- where Lm(e) reaches 
        /// as Lm(e) = argmax(r)(rv_t)
        inline std::pair<T,T> maxL2(T e)const{
            if(e <= 0.5){
                return {1 - e,1};
            }
            else {
                auto FuncArray = 
                grob::as_container([&](size_t i){
                    //want to get F'(u) = phi(r) - r^2/2Q(r)
                    auto r = Phi.Grid[i];
                    return Phi[i] - r*r*Q[i]/2;
                },Phi.Grid.size());
                size_t i = grob::find_last_index(FuncArray,e,std::greater<>{});
                auto r0 = Phi.Grid[i];
                auto r1 = Phi.Grid[i+1];
                auto f0 = FuncArray[i] - e;
                auto f1 = FuncArray[i+1] - e;
                auto r = sqrt((r1*r1*f0-r0*r0*f1)/(f0-f1));
                return {r*r*(Phi(r)-e),r};
            }
        }

        /// @brief Lmax^2(e) and rm(e) 
        /// @param e positive dimentionless energy
        /// @return pair(Lmax^2(e),rm(e)), where rm(e) --- where Lm(e) reaches 
        /// as Lm(e) = argmax(r)(rv_t)
        inline std::pair<T,T> maxL2d(T e)const{
            if(e <= 0.5){
                return {1/(4*e),1/(2*e)};
            }
            else {
                return maxL2(e);
            }
        }

        ////CHECKED, OK
        /// @brief return S function
        /// S[phi](r,rm,rp) = [u = r^2, F(u) = r^2*phi(r)] 
        /// = -( (F(up)-F(u))/(up-u) - (F(u-F(um))/(u-um) )/(up-um) 
        /// @param r 
        /// @param rm 
        /// @param rp 
        /// @return
        inline T S(T r,T rm,T rp)const{
            T h = Rho.Grid.h();
            T u = r*r;
            T um = rm*rm;
            T up = rp*rp;
            auto F = [&](T r){return r*r*Phi(r);};
            auto F1 = [&](T r){return Phi(r) - Q(r)*r*r/2;};
            auto mF2 = [&](T r){return (3*Rho(r)+Q(r))/8;};

            auto _dr_max = low_steps_num*h*(1+2*h/(rp+rm));

            if(rp-rm <= _dr_max*2){
                T u_s = (u+um+up)/3;
                T r_s = sqrt((u+um+up)/3);
                auto F2_s = (r_s <= 1 ? mF2(r_s) : -_DD_F_C(u_s)/2); 
                return F2_s;
            } 
            else if (r - rm <=_dr_max){
                T r_s = sqrt((u+um)/2);
                auto Frp = (rp <= 1 ? F(rp) : __F_C(up)); 
                return (F1(r_s) - (Frp-F(r))/(up-u))/(up-um);
            }
            else if (rp - r <= _dr_max){
                T u_s = (u+up)/2;
                T r_s = sqrt(u_s);
                auto F1_s = (r_s <= 1 ? F1(r_s) : _D_F_C(u_s)); 
                return ((F(r)-F(rm))/(u-um) - F1_s)/(up-um);
            }
            else {
                auto Frp = (rp <= 1 ? F(rp) : __F_C(up)); 
                return ((F(r)-F(rm))/(u-um) - (Frp-F(r))/(up-u))/(up-um);
            }
        }


        inline T __F_C(T const &u)const{
            auto u1 = u - 1;
            return 1 + u1/2 - u1*u1*(3*Rho_max+1)/8;
        }
        inline T _D_F_C(T const &u)const{
            return 1./2 - (u-1)*(3*Rho_max+1)/4;
        }
        inline T _DD_F_C(T const &u)const{
            return -(3*Rho_max+1)/4;
        }
        inline T _DD_F(T const& u) const {
            return -(3 * Rho(std::sqrt(u)) + 1) / 4;
        }
        inline T _DD_F_r(T const& r) const {
            return -(3 * Rho(r) + 1) / 4;
        }
        //NOT_CHECKED
        /// @brief finds u, that more than one in quadratic potential
        /// @param e positive energy, e<1/2
        /// @param L2 is squared of dimless momentum l^2
        /// @return up(e,l) > 1
        inline T _find_u_more_1(T const & e,T const & L2)const{
            auto Rho_max3 = 1 + 3*Rho_max; 
            auto e2 = 1-2*e;  
            return (
                        2*e2+Rho_max3 + 
                        2*std::sqrt((1 + e2 - 2*L2)*Rho_max3 + e2*e2)
                    )/Rho_max3;
        }


        //NOT_CHECKED
        /// @brief finds rmin or rmax (e,l) between indexes i0,i1
        /// @tparam b if 0 than rmin if 1 than rmax
        /// @param e positive energy
        /// @param L2 positive squared momentum
        /// @param i0 
        /// @param i1 
        /// @return 
        template <size_t b>
        inline T find_r(T const & e,T const & L2,
        size_t i0,size_t i1)const{
            auto F = [&](size_t i){
                auto r = Phi.Grid[i];
                return r*r*(Phi[i] - e);
            };

            size_t i = grob::find_index(
                grob::as_container(F,Phi.Grid.size()),
                L2,i0,i1,
                typename std::conditional<
                    b == 0,std::less<>,std::greater<>
                >::type {});
            
            auto [_i0,_i1] = (b==0 ? std::make_pair(i, i+1) :
                             std::make_pair(i+1, i)); 
            auto r0 = Phi.Grid[_i0];
            auto r1 = Phi.Grid[_i1];
            auto u0 = r0*r0;
            auto u1 = r1*r1;
            
            auto u = std::clamp(
                __detail::NewtonStep(u0,u0*(Phi[_i0]-e)-L2,Phi[_i0]-u0*Q[_i0]/2-e,-(3*Rho[_i0]+Q[_i0])/4)
                ,( b == 0 ? u0 : u1),( b == 0 ? u1 : u0));
            return sqrt(u);
        }

        inline std::pair<T,T> find_close_rmin_rmax(T const& e,T const& L2,size_t im) const{
            auto r0 = Phi.Grid[im];
            // f0 = r^2(phi(r) - e) - l^2, r = r0
            auto u0 = r0*r0;
            T f0 = u0*(Phi[im] - e) - L2;
            // f1 = df/du, r = r0
            T f1 = Phi[im] - Q[im]*u0/2 - e;
            T f2 = (3*Rho[im]+Q[im])/4;
            auto msqrt = std::sqrt(std::max((T)0,f1*f1+2*f0*f2))/f2;
            auto avar = u0 + f1/f2;
            return {std::sqrt(avar-msqrt),std::sqrt(avar+msqrt)};
        }

        inline std::pair<T,T> find_rmin_rmax(T const& e,T const& L2,LE_func_t const &LE,
        size_t imin0_g = 0,size_t imin1_g = std::numeric_limits<size_t>::max(),
        size_t imax0_g = 0,size_t imax1_g = std::numeric_limits<size_t>::max())const{
            //more e -> less rm -> less im
            auto rm = LE.rm(e);
            size_t im = Phi.Grid.pos(rm);
            auto F = [&](size_t i){
                auto r = Phi.Grid[i];
                return r*r*(Phi[i] - e);
            }; 
            if(F(im) > L2 && F(im+1) > L2){
                T r_min = find_r<0>
                (
                    e,L2,
                    std::max((size_t)0,imin0_g),
                    std::min(im,imin1_g)
                );
                // if im == size()-2 than this condition failed
                // im can't be size()-1
                size_t i_max = Phi.Grid.size()-1;
                if(F(i_max) < L2){
                    T r_max = find_r<1>(e,L2,
                        std::max(im+1,imax0_g),
                        std::min(i_max,imax1_g));
                    return {r_min,r_max};
                }
                else {
                    return {r_min,std::sqrt(_find_u_more_1(e,L2))};
                }
            } else {
                //auto F = [&](T r){return r*r*Phi(r);};
                //auto F1 = [&](T r){return Phi(r) - Q(r)*r*r/2;};
                //auto mF2 = [&](T r){return (3*Rho(r)+Q(r))/8;};
                return find_close_rmin_rmax(e,L2,im);
            }
        }

        /// @brief returns i0, i1, so that r[i0] <= rmin <= rmax <= r[i1]
        /// @param rmin 
        /// @param rmax 
        /// @return 
        inline std::pair<size_t,size_t> _i0i1r(T const & rmin,T const & rmax){
            return {Phi.Grid.pos(rmin),Phi.Grid.pos(rmax)+1};
        }
        
        /// @brief get grid function of tau(theta/theta1)
        /// tau_max = real Tin/theta_max
        /// @tparam U type of trajectory function
        /// @param rmin 
        /// @param rmax 
        /// @param Nbins 
        /// @return tuple(theta_max,tau_max,Taus)
        template <typename U = T>
        inline auto get_internal_traj(T const & rmin,T const & rmax,size_t Nbins)const{
            auto u0 = rmin*rmin;
            auto u1 = rmax*rmax;
            
            auto u_av = (u0+u1)/2;
            auto u_d = (u0-u1)/2;
            
            

            T theta_max = std::numbers::pi;
            if (u1 > 1) {
                auto pot_cos = (1 - u_av) / u_d;
                theta_max = std::acos(
                    pot_cos < -1 ? -1 : (
                        pot_cos <= 1 ? pot_cos : 1
                    )
                );
            }

            auto R = [&](auto theta){
                return std::sqrt(u_av + u_d*std::cos(theta));
            };
            auto dTau_dTheta = [&](auto theta){
                return 1/(2*std::sqrt(S(R(theta*theta_max),rmin,rmax)));
            };

            using GridR_mod = grob::GridUniform<U>;
            using RFunc2_t_mod = grob::GridFunction<
                grob::interpolator_spline1D,
                GridR_mod, std::vector<U>>;
            RFunc2_t_mod Taus(GridR_mod(0,1,Nbins+1),std::vector<U>(Nbins+1));
            Taus[0] = 0;
            for(size_t i=1;i<Taus.size();++i){
                Taus[i] = Taus[i-1] + integrateAB2(dTau_dTheta,Taus.Grid[i-1],Taus.Grid[i],1);
            }
            const auto tau_max = Taus[Taus.size()-1];
            for(auto & tau : Taus.Values){
                tau *= (1/tau_max);
            }
            return std::make_tuple(theta_max,tau_max,Taus);
        }
        inline auto get_internal_period(T const & rmin,T const & rmax,size_t Nbins){
            auto u0 = rmin*rmin;
            auto u1 = rmax*rmax;
            
            auto u_av = (u0+u1)/2;
            auto u_d = (u0-u1)/2;
            auto pot_cos = (1-u_av)/u_d;
            auto theta_max = std::acos(
                pot_cos < -1 ? -1 : (
                    pot_cos <= 1 ? pot_cos : 1
                )  
            );

            auto R = [&](auto theta){
                return std::sqrt(u_av + u_d*std::cos(theta));
            };
            auto dTau_dTheta = [&](auto theta){
                return 1/(2*std::sqrt(S(R(theta),rmin,rmax)));
            };
            GridR ThetaGrid(0,theta_max,Nbins+1);
            T tau = 0;
            for(size_t i=1;i<ThetaGrid.size();++i){
                tau += integrateAB2(dTau_dTheta,ThetaGrid[i-1],ThetaGrid[i],1);
            }
            return tau;
        }
    };
};

#endif//BODY_POTENTIAL_HPP