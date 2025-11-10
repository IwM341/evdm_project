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
#include <chrono>
#include <span>

#ifdef _MSC_VER
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif




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
    concept Printable = requires (T const & t,std::ostream & os){ {os<<t};};
    template <Printable T>
    std::string to_debug_string(const T &x){
        std::ostringstream os;
        os <<  x;
        return os.str();
    }

    template <typename T>
    concept VectorSimilar = requires(T const& x){
        {x.size()};
        {x[size_t{}]};
    };

    template <VectorSimilar T>
    std::string to_debug_string_iter(T const & x){
        std::ostringstream os;
        os << "array[";
        if(x.size() >= 1)
            os <<  x[0];
        for(size_t i=1;i<x.size();++i){
            os << ", " << to_debug_string(x[i]);
        }
        os << "]";
        return os.str();
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
    std::string to_debug_string_tuple(const std::tuple<Args...> &Tp){
        return "(" +
                __tuple_elements_string_list<
                        std::tuple<Args...>,
                        std::tuple_size<std::tuple<Args...>>::value
                    >::get_list(Tp,", ") + ")";
    }
};

template <typename U,typename V>
std::ostream & operator << (std::ostream & os, const std::pair<U,V> & P){
    os << "pair(" << P.first << ", " << P.second << ")";
    return os;
}

template <typename...U>
std::ostream & operator << (std::ostream & os, const std::tuple<U...> & TP){
    os << debugdefs::to_debug_string_tuple(TP);
    return os;
}

template <typename T,typename...Args>
std::ostream & operator << (std::ostream & os,std::vector<T,Args...> const & V){
    return os << debugdefs::to_debug_string_iter(V);
} 
template <typename T,std::size_t N>
std::ostream & operator << (std::ostream & os,std::span<T,N> const & V){
    return os << debugdefs::to_debug_string_iter(V);
} 

template <typename T,size_t N>
std::ostream & operator << (std::ostream & os,std::array<T, N> const & V){
    return os << debugdefs::to_debug_string_iter(V);
} 
#define SVAR(x) (std::string(#x) + std::string(" = ") + debugdefs::to_debug_string(x))
#define PVAR(x) println(SVAR(x))

#ifndef _NDEBUG
#define DEBUG(x) println(SVAR(x))
#else
#define DEBUG(x)
#endif


#ifdef _VERBOSE
#ifndef _NDEBUG
#define DEBUG1(x) println(SVAR(x))
#endif
#else
#define DEBUG1(x)
#endif


#define PDEL() (std::cout << "-----------------------------------------------" <<std::endl)
#define SCOMPARE(x,y) (SVAR(x) + " vs " + SVAR(y))
#define COMPARE(x,y) std::cout << SCOMPARE(x,y) <<std::endl;

template <typename T>
inline std::string TypeString(){
    std::string fname = __PRETTY_FUNCTION__ ;
    return fname.substr(35,fname.size() - 49-34);
}
template <typename T>
inline std::string TypeString(T &&x){
    return TypeString<std::decay_t<T>>();
}
#define TypeToString(type) TypeString<type>()

template <typename EnumClass>
struct EnNameStr{
	template <EnumClass em>
    inline static std::string _str(){
		std::string enum_name = __PRETTY_FUNCTION__ ;
		return  enum_name.substr(69,enum_name.size()-135);
	}
};
#define ValueToString(V) EnNameStr<decltype(V)>::_str<V>()

inline void println(){
    std::cout <<std::endl;
}
template <typename T>
inline void println(T &&data){
    std::cout << data <<std::endl;
}

template <typename ...Args, typename T>
inline void println(T &&data,Args&&...args){
    std::cout << data;
    println(args...);
}

template <typename Delimtype,typename T>
inline void printd(Delimtype delim,T data){
    std::cout << data <<std::endl;
}

template <typename Delimtype,typename ...Args, typename T>
inline void printd(Delimtype delim,T data,Args...args){
    std::cout << data << delim;
    printd(delim,args...);
}


inline std::string make_path(const std::string & str){
	return std::string("\"") + std::regex_replace(str, std::regex("\\\\"), "\\\\") + "\"";
}


inline double rnd(){
    static std::random_device __rd;
    static std::mt19937 gen(__rd());
    static std::uniform_real_distribution<> dis;
    return dis(gen);
}
inline float rndf(){
    static std::random_device __rd;
    static std::mt19937 gen(__rd());
    static std::uniform_real_distribution<float> dis;
    return dis(gen);
}


template <typename T>
inline T norm1(T const &x, T const & y){
    return std::abs( (x-y)/(x+y));
}
template <typename T>
inline T norm2(T const &x, T const & y){
    return std::abs( (x*x-y*y)/(x*x+y*y));
}

auto time_now(){
    return std::chrono::high_resolution_clock::now();
}

template <typename T>
auto time_count_micros(T start,T stop){
    return std::chrono::duration<double, std::micro>(stop - start).count();
}
template <typename T>
auto time_count_millis(T start,T stop){
    return std::chrono::duration<double, std::milli>(stop - start).count();
}
template <typename T>
auto time_count_sec(T start,T stop){
    return std::chrono::duration<double>(stop - start).count();
}


#define debug_return(expr) PVAR((expr)); return expr;
