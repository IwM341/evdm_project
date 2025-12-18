#ifndef MK_HPP
#define MK_HPP
#include "prng.hpp"
#include "cmath"
#include <numbers>
#include <utility>
namespace evdm{
    template <typename ResultType,typename DensityType = ResultType>
    struct MCResult{
        ResultType Result;
        DensityType RemainDensity;
        MCResult(ResultType Result, DensityType RemainDensity = 1.0):Result(Result),RemainDensity(RemainDensity){}
        MCResult():RemainDensity(0){}
    };


    template <typename T>
    struct MCIntegral{
        T Result;
        T Sigma;
        MCIntegral(T Result,T Sigma = 0):Result(Result),Sigma(Sigma){}

        operator T (){
            return Result;
        }
        operator std::pair<T,T> (){
            return std::pair<T,T>(Result,Sigma);
        }
    };

    template<typename GenType>
    auto MCIntegrate(GenType G,size_t N){
        typedef typename std::decay<decltype(G())>::type value_type;
        value_type sum = 0;
        value_type sum2 = 0;
        for (size_t i=0;i<N;++i){
            auto f = G();
            sum += f;
            sum2 += f*f;
        }
        sum = sum/N;
        sum2 = sum2/N;
        return MCIntegral<value_type>(sum,sqrt(sum2-sum*sum));
    }


    template <typename Gen_t>
    using Gen_vt = std::decay_t<Gen_t>::value_type;

    template <class Gen_t>
    inline Gen_vt<Gen_t> RandomPhi(Gen_t&& G) {
        return G() * 2 * std::numbers::pi;
    }

    template <class Gen_t>
    inline Gen_vt<Gen_t> RandomCos(Gen_t&& G) {
        return G() * 2 - 1;
    }

    template <class Gen_t>
    inline MCResult<Gen_vt<Gen_t>> RandomCos_c(Gen_t&& G, Gen_vt<Gen_t> cosTheta_min, Gen_vt<Gen_t> cosTheta_max) {
        return CResult<Gen_vt<Gen_t>>(
            cosTheta_min + G() * (cosTheta_max - cosTheta_min), 
            (cosTheta_max - cosTheta_min)/2
            );
    }

    template <class Gen_t>
    inline Gen_vt<Gen_t> RandomPhi(Gen_t&& G, Gen_vt<Gen_t> phi0, Gen_vt<Gen_t> phi1) {
        return phi0 + G() * (phi1 - phi0);
    }

    template <class Gen_t>
    inline Gen_vt<Gen_t> Gauss2Norm(Gen_t&& G, Gen_vt<Gen_t> Vdisp) {
        if constexpr (G.bounds == genpol::E0I1)
            return Vdisp * std::sqrt(-2 * std::log(G()));
        else if (G.bounds == genpol::I0E1)
            return Vdisp * std::sqrt(-2 * std::log(1 - G()));
        else
            static_assert(true, "incorrect generator bounds");
    }
    template <class Gen_t>
    inline Gen_vt<Gen_t> Gauss2Norm2(Gen_t&& G, Gen_vt<Gen_t> Vdisp) {
        if constexpr (G.bounds == genpol::E0I1)
            return Vdisp* Vdisp * ( - 2 * std::log(G()));
        else if (G.bounds == genpol::I0E1)
            return Vdisp * Vdisp * (-2 * std::log(1 - G()));
        else
            static_assert(true, "incorrect generator bounds");
    }

    template <class Gen_t>
    inline Gen_vt<Gen_t> Gauss(Gen_t&& G, Gen_vt<Gen_t> Vdisp) {
        Gen_vt<Gen_t> V = Gauss2Norm(G, Vdisp);
        return V * cos(RandomPhi(G));
    }

    template <class Gen_t>
    inline MCResult<Gen_vt<Gen_t>> Gauss3_min_abs(Gen_t&& G, Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> Vmin) {
        Vmin = (Vmin > 0 ? Vmin : 0);
        Gen_vt<Gen_t> a0 = exp(-Vmin * Vmin / (2 * Vdisp * Vdisp));

        Gen_vt<Gen_t> v_nd; 
        if constexpr (G.bounds == genpol::E0I1)
            v_nd = std::sqrt(-2 * std::log(a0 * (G())));
        else if (G.bounds == genpol::I0E1)
            v_nd = std::sqrt(-2 * std::log(a0 * (1 - G())  )  );
        else
            static_assert(true, "incorrect generator bounds");

        return MCResult<Gen_vt<Gen_t>>(v_nd * Vdisp, sqrt(2.0 / std::numbers::pi) * v_nd * a0);
    }

};
#endif//MK_HPP