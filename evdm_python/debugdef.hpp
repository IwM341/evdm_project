#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <type_traits>
#include <stdio.h>
#include <list>
#include <set>
#include <regex>
#include <mutex>
#include <random>



template <typename U,typename V>
std::ostream & operator << (std::ostream & os, const std::pair<U,V> & P){
    os << "pair(" << P.first << ", " << P.second << ")";
    return os;
}


struct _ff{
    struct _format_args{
        size_t prec;
        bool is_exp;
    };
    template <typename T>
    struct printable{
        T value;
        _format_args _format;
        operator std::string () const{
            std::stringstream oss;
            oss << std::setprecision(_format.prec);
            if(_format.is_exp){
                oss << std::scientific;
            }
            oss << value;
            return oss.str();
        }
        printable(T value,_format_args _format) : value(value), _format(_format){}
        inline friend std::ostream & operator << 
            (std::ostream & os,printable const & __val)
        {
            return  os << static_cast<std::string>(__val);
        }
    };

    constexpr inline _ff(size_t prec = 6,bool is_exp = false):
        _fmt{prec,is_exp}{}

    template <typename T>
    inline friend auto operator *(T const &value,_ff _format){
        return printable<T>(value,_format._fmt);
    }
    template <typename T>
    inline friend auto operator *(_ff _format,T const &value){
        return printable<T>(value,_format._fmt);
    }

    private:
    _format_args _fmt;
};

namespace debugdefs{




    template <typename T>
    std::string to_debug_string(const T &x){
        std::ostringstream os;
        os <<  x;
        return os.str();
    }

    template <typename T,size_t N>
    std::string to_debug_string(const std::array<T,N> &x){
        std::ostringstream os;
        os << "array[";
        if(N >= 1)
            os <<  x[0];
        for(size_t i=1;i<N;++i){
            os << ", " << x[i];
        }
        os << "]";
        return os.str();
    }

    template<typename T, typename U>
    std::string to_debug_string(const std::pair<T,U> &x){
        return "{" + to_debug_string(x.first)+", "+to_debug_string(x.second) + "}";
    }

    template <typename Tuple,int remain>
    struct __tuple_elements_string_list{
        constexpr static int i = std::tuple_size<Tuple>::value - remain;
        static std::string get_list(Tuple const & _Tp,const std::string &delim){
            return to_debug_string(std::get<i>(_Tp)) + delim +
                    __tuple_elements_string_list<Tuple,remain-1>::get_list(_Tp,delim);
        }
    };

    template <class Tuple>
    struct __tuple_elements_string_list<Tuple,0>{
        static std::string get_list(Tuple const & _Tp,const std::string &delim){
            return "";
        }
    };

    template<typename...Args>
    std::string to_debug_string(const std::tuple<Args...> &Tp){
        return "(" +
                __tuple_elements_string_list<
                        std::tuple<Args...>,
                        std::tuple_size<std::tuple<Args...>>::value
                    >::get_list(Tp,", ") + ")";
    }
};
#define SVAR(x) (std::string(#x) + std::string(" = ") + debugdefs::to_debug_string(x))
#define PVAR(x) print(SVAR(x))

#ifndef _NDEBUG
#define VERBOSE_VAR(x)
#else
#define VERBOSE_VAR(x)
#endif


#ifdef _VERBOSE
#ifndef _NDEBUG
#define VERBOSE_VAR1(x)
#endif
#else
#define VERBOSE_VAR1(x)
#endif
