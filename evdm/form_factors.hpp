#ifndef FORM_FACTORS_HPP
#define FORM_FACTORS_HPP
#include <cmath>
#include <array>
#include "utils/mc.hpp"
#include "utils/variant_tools.hpp"
#include "utils/polynom.hpp"
#include <tuple>
#include <string>
#include <sstream>
namespace evdm{

    namespace detail{
        template <size_t N>
        struct coeff_evaulator{
            constexpr static size_t N_2 = ( N % 2 ? (N+1)/2 : N/2);

            template <typename T>
            static inline constexpr T eval(const T * data,T const &x){
                return coeff_evaulator<N_2>(data,x) + 
                        coeff_evaulator<N-N_2>(data + N_2,x)*x;
            }
        };
        template <>
        struct coeff_evaulator<0>{
            template <typename T>
            static inline constexpr T eval(const T * data,T const &x){
                return data[0];
            }
        };
    };

    template <size_t _size>
    using Polynom = PolynomHorner<float,_size,16>;

    struct ElasticFactorBase{
        inline MCResult<float,size_t> EnergyLoss(float Emax) const noexcept{
            return 0;
        }
        inline float ScatterFactor(float q_2,float v_2,float Eloss)const noexcept{
            return 1;
        }
    };

    inline std::string degree_print(int deg){
        if (deg < 0)
            return "*y^(" + std::to_string(deg) + ")";
        else if (deg == 0) {
            return "";
        }
        else if (deg == 1) {
            return "*y";
        }
        else {
            return "*y^" + std::to_string(deg);
        }
    }

    template <size_t _poly_size,bool y_inv = false>
    struct QexpFactor : ElasticFactorBase{
        using Poly = Polynom<_poly_size>;
        Poly _Pol;
        float b_2;

        QexpFactor(){}
        template <typename Array_t>
        QexpFactor(Array_t const & coeffs,float b):_Pol(coeffs),b_2(b*b){}
        inline float ScatterFactor(float q_2,float v_2,float Eloss)const noexcept{
            float y = (b_2*q_2/4);
            if constexpr(y_inv)
                return _Pol(y)*exp(-2*y)/y;
            else {
                return _Pol(y)*exp(-2*y);
            }
        }
        inline std::string repr() const{
            std::ostringstream S;
            S << "exp(-2y)*(";
            for (size_t i = 0; i < _poly_size; ++i) {
                if (i != 0 && _Pol[i] >= 0) {
                    S << " + ";
                }
                S << _Pol[i] << degree_print((int)i - (y_inv ? 1 : 0));
            }
            S << ")";
            return S.str();
        }
        template <typename T>
        inline void rescale(T x) {
            _Pol.rescale(x);
        }
    };

    template <size_t _poly_size,bool y_inv  = false>
    struct QexpFactor_v : ElasticFactorBase{
        using Poly = Polynom<_poly_size>;
        Poly _P_0;
        Poly _P_V;
        float b_2;

        QexpFactor_v(){}
        template <typename Array1_t,typename Array2_t>
        QexpFactor_v(Array1_t const & coeffs_v0,Array2_t const & coeffs_v1,float b):
            _P_0(coeffs_v0),_P_V(coeffs_v1),b_2(b*b){}

        inline float ScatterFactor(float q_2,float v_2,float Eloss)const noexcept{
            float y = (b_2*q_2/4);
            Poly eff_poly = Poly::x_plus_by(_P_0,_P_V,v_2);
            if constexpr(y_inv)
                return eff_poly (y)*exp(-2*y)/y;
            else {
                return eff_poly (y)*exp(-2*y);
            }
        }
        inline std::string repr()const {
            std::ostringstream S;
            S << "exp(-2y)*(";
            S << "{";
            for (size_t i = 0; i < _poly_size; ++i) {
                if (i != 0 && (_P_0[i] >= 0)) {
                    S << " + ";
                }
                S << _P_0[i] << degree_print((int)i - (y_inv ? 1 : 0));
            }
            S << "} + v^2*{";
            for (size_t i = 0; i < _poly_size; ++i) {
                if (i != 0 && _P_V[i] >= 0) {
                    S << " + ";
                }
                S << _P_V[i] << degree_print((int)i - (y_inv ? 1 : 0));
            }
            S << "}";
            S << ")";
            return S.str();
        }

        template <typename T>
        inline void rescale(T x) {
            _P_0.rescale(x);
            _P_V.rescale(x);
        }
    };

    namespace __index_detail{
        template<size_t i,typename IndexSeq>
        struct get_index;
        
        
        template<size_t i,size_t I0,size_t...I>
        struct get_index<i,std::index_sequence<I0,I...>>{
            constexpr static size_t value = (i==0 ? I0 : get_index<i-1,std::index_sequence<I...>>::value );
        };
        template<size_t i,size_t I0>
        struct get_index<i,std::index_sequence<I0>>{
            constexpr static size_t value = I0;
        };

        

        template <typename IndexSeq>
        struct index_seq_last;

        template <size_t I0>
        struct index_seq_last<std::index_sequence<I0>>{
            static constexpr size_t value = I0;
        };

        template <size_t I0,size_t...I>
        struct index_seq_last<std::index_sequence<I0,I...>>{
            static constexpr size_t value = index_seq_last<std::index_sequence<I...>>::value;
        };

        template <typename IndexSeq>
        struct index_seq_first;

        template <size_t I0,size_t...I>
        struct index_seq_first<std::index_sequence<I0,I...>>{
            static constexpr size_t value = I0;
        };

        template <typename IS1,size_t I>
        struct index_seq_push_back;
        
        template <size_t I,size_t...I1>
        struct index_seq_push_back<std::index_sequence<I1...>,I>{
            typedef std::index_sequence<I1...,I> type;
        };

        template <typename IS1,size_t I>
        struct index_seq_push_front;
        
        template <size_t I,size_t...I1>
        struct index_seq_push_front<std::index_sequence<I1...>,I>{
            typedef std::index_sequence<I,I1...> type;
        };


        template <typename IS>
        struct index_seq_pop_front;
        
        template <size_t I,size_t...I1>
        struct index_seq_pop_front<std::index_sequence<I,I1...>>{
            typedef std::index_sequence<I1...> type;
        };
        
        template <typename T1,typename T2>
        struct type_pair{
            typedef T1 first;
            typedef T2 second;
        };

        template <typename IS1,typename IS2>
        constexpr auto index_seq_split_rearrange(type_pair<IS1,IS2>){
            if constexpr(IS1::size() < IS2::size()){
                typedef typename index_seq_push_back<IS1,index_seq_first<IS2>::value>::type IS1_n;
                typedef typename index_seq_pop_front<IS2>::type IS2_n;
                typedef type_pair<IS1_n,IS2_n> tuple_pair_new;
                return index_seq_split_rearrange(tuple_pair_new{});
            } else {
                return type_pair<IS1,IS2>{};
            }
        }

        template <typename IndexSequence>
        struct index_seq_split : 
            decltype(
                index_seq_split_rearrange(
                        type_pair<std::index_sequence<>,IndexSequence>{}
                    )
                ){};

        template <size_t I1>
        size_t find_first_moreq_index(size_t i, std::index_sequence<I1>) {
            return 0;
        }

        template <size_t I1,size_t... I>
        inline size_t find_first_moreq_index(size_t i, std::index_sequence<I1,I...>) {
            typedef std::index_sequence<I1,I...> IS;
            typedef index_seq_split<IS> type;
            if( i <= index_seq_last<typename type::first>::value){
                return find_first_moreq_index(i,typename type::first{});
            } else {
                return type::first::size() + find_first_moreq_index(i,typename type::second{});
            }
        }
    };

    template <typename IndexSeq>
    struct QexpFactors;

    template <size_t...I>
    struct QexpFactors<std::index_sequence<I...>> :public std::variant<
            QexpFactor<I,false>...,QexpFactor<I,true>...,
            QexpFactor_v<I,false>...,QexpFactor_v<I,true>...
        >{
        typedef std::index_sequence<I...> PolySizes;

        constexpr static bool Y_INV = true; 

        typedef std::variant<
            QexpFactor<I,false>...,QexpFactor<I,true>...,
            QexpFactor_v<I,false>...,QexpFactor_v<I,true>...
        > Base;

        typedef std::variant<QexpFactor<I>...> Base_v0;
        typedef std::variant<QexpFactor<I,Y_INV>...> Base_yinv_v0;
        typedef std::variant<QexpFactor_v<I>...> Base_v1;
        typedef std::variant<QexpFactor_v<I,Y_INV>...> Base_yinv_v1;

        //typedef QexpFactors<std::index_sequence<I...>> this_t;
        QexpFactors(){}
        
        template <typename Array_t>
        inline static auto MakeBase_v0(bool y_inv,double b,Array_t const & coeffs){
            if(!y_inv){
                return variant_cast<Base>(
                    make_variant<Base_v0>(
                        __index_detail::find_first_moreq_index(coeffs.size(),PolySizes{}),
                        [b,&coeffs](auto t){
                            return typename decltype(t)::type (coeffs,b);
                        }
                    )
                );
            } else {
                return variant_cast<Base>(
                    make_variant<Base_yinv_v0>(
                        __index_detail::find_first_moreq_index(coeffs.size(),PolySizes{}),
                        [b,&coeffs](auto t){
                            return typename decltype(t)::type (coeffs,b);
                        }
                    )
                );
            }
        }

        template <typename Array1_t,typename Array2_t>
        inline static auto MakeBase_v1(bool y_inv,double b,Array1_t const & coeffs_v0,Array2_t const & coeffs_v1){
            if(!y_inv){
                return variant_cast<Base>(
                    make_variant<Base_v1>(
                        __index_detail::find_first_moreq_index(std::max(coeffs_v0.size(),coeffs_v1.size()),PolySizes{}),
                        [&](auto t){
                            return typename decltype(t)::type (coeffs_v0,coeffs_v1,b);
                        }
                    )
                );
            } else {
                return variant_cast<Base>(
                    make_variant<Base_yinv_v1>(
                        __index_detail::find_first_moreq_index(std::max(coeffs_v0.size(),coeffs_v1.size()),PolySizes{}),
                        [&](auto t){
                            return typename decltype(t)::type (coeffs_v0,coeffs_v1,b);
                        }
                    )
                );
            }
        }

        template <typename Array_t>
        QexpFactors(bool y_inv,double b,Array_t const & coeffs): 
        Base(vmove(MakeBase_v0(y_inv,b,coeffs))){}

        template <typename Array1_t,typename Array2_t>
        QexpFactors(bool y_inv,double b,Array1_t const & coeffs_v0,Array2_t const & coeffs_v1): 
        Base(vmove(MakeBase_v1(y_inv,b,coeffs_v0,coeffs_v1))){}

        inline Base & as_variant(){
            return static_cast<Base &>(*this);
        }
        inline Base const & as_variant() const{
            return static_cast<Base const &>(*this);
        }
        std::string repr() const {
            return std::visit([](auto const& Poly) ->std::string{
                return Poly.repr();
            }, as_variant());
        }

        template <typename T>
        void rescale(T sc) {
            std::visit([](auto const& Poly) ->std::string {
                return Poly;
            }, as_variant());
        }
    };

};

#endif//FORM_FACTORS_HPP