#ifndef POLYNOM_HPP
#define POLYNOM_HPP
#include <iostream>


namespace evdm{
    template <typename T,size_t _size,size_t _align_el = 1>
    struct PolynomHorner{
        alignas(_align_el*sizeof(T)) T coeffs[_size];

        constexpr static size_t size(){return _size;}

        inline constexpr T const & operator [](size_t i) const {
            return coeffs[_size-1-i];
        }
        inline constexpr T  & operator [](size_t i) {
            return coeffs[_size-1-i];
        }

        inline constexpr T operator () (T x) const{
            T sum = coeffs[0];
            for(size_t i=1;i<_size;++i){
                sum = coeffs[i] + sum*x;
            }
            return sum;
        }

        PolynomHorner(){}
        template <typename Array_t>
        PolynomHorner(Array_t && cfs,size_t _array_size = _size){
            int N = _array_size - 1;
            for(auto it = cfs.begin();it != cfs.end() && N >= 0;++it,--N){
                coeffs[N] = *it;
            }
            for(;N>=0;--N){
                coeffs[N] = 0;
            }
        }
        template <typename Q>
        PolynomHorner(std::initializer_list<Q> cfs,size_t _array_size = _size){
            int N = _array_size - 1;
            for(auto it = cfs.begin();it != cfs.end() && N >= 0;++it,--N){
                coeffs[N] = *it;
            }
            for(;N>=0;--N){
                coeffs[N] = 0;
            }
        }
        static inline constexpr PolynomHorner x_plus_by(PolynomHorner const & X,PolynomHorner const & Y,T b){
            PolynomHorner Z;
            for(size_t i=0;i<_size;++i){
                Z.coeffs[i] = X.coeffs[i] + b*Y.coeffs[i];
            }
            return Z;
        }
    };
};

#endif//POLYNOM_HPP