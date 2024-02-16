#ifndef DEBUGDEF_H
#define DEBUGDEF_H

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

#ifndef __GNUC__
#define __PRETTY_FUNCTION__ "__PRETTY_FUNCTION__ macro is not supported"
#endif

template <typename U,typename V>
std::ostream & operator << (std::ostream & os, const std::pair<U,V> & P){
    os << "pair(" << P.first << ", " << P.second << ")";
    return os;
}

/// @brief struct for printing numbers
/// use: auto format = _ff(precision,is_exponental)
/// then cout << format*4 << endl;
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

#define PDEL() (std::cout << "-----------------------------------------------" <<std::endl)
#define SCOMPARE(x,y) (SVAR(x) + " vs " + SVAR(y))
#define COMPARE(x,y) std::cout << SCOMPARE(x,y) <<std::endl;

#define _P std::make_pair
#define _T std::make_tuple

template <typename T>
std::string TypeString(){
    std::string fname = __PRETTY_FUNCTION__ ;
    return fname.substr(35,fname.size()-84);
}
#define TypeToString(type) TypeString<type>()

template <typename EnumClass>
struct EnNameStr{
	template <EnumClass em>
	static std::string _str(){
		std::string enum_name = __PRETTY_FUNCTION__ ;
		return  enum_name.substr(69,enum_name.size()-135);
	}
};
#define ValueToString(V) EnNameStr<decltype(V)>::_str<V>()

void print(){
    std::cout <<std::endl;
}
template <typename T>
void print(T data){
    std::cout << debugdefs::to_debug_string(data) <<std::endl;
}

template <typename ...Args, typename T>
void print(T data,Args...args){
    std::cout << debugdefs::to_debug_string(data);
    print(args...);
}

template <typename Delimtype,typename T>
void printd(Delimtype delim,T data){
    std::cout << data <<std::endl;
}

template <typename Delimtype,typename ...Args, typename T>
void printd(Delimtype delim,T data,Args...args){
    std::cout << data << delim;
    printd(delim,args...);
}


std::string make_path(const std::string & str){
	return std::string("\"") + std::regex_replace(str, std::regex("\\\\"), "\\\\") + "\"";
}


double rnd(){
    static std::random_device __rd;
    static std::mt19937 gen(__rd());
    static std::uniform_real_distribution<> dis;
    return dis(gen);
}
float rndf(){
    static std::random_device __rd;
    static std::mt19937 gen(__rd());
    static std::uniform_real_distribution<float> dis;
    return dis(gen);
}


template <typename T>
T norm1(T const &x, T const & y){
    return std::abs( (x-y)/(x+y));
}
template <typename T>
T norm2(T const &x, T const & y){
    return std::abs( (x*x-y*y)/(x*x+y*y));
}

//plot with lines
// gp << "plot" << gp.file1d(...) << "with lines"
template <typename GridFunc_t>
auto make_cols(GridFunc_t const & V){
    std::vector<std::decay_t<decltype(V.Grid[0])>> X(V.size());
    for(size_t i=0;i<V.size();++i){
        X[i] = V.Grid[i];
    }
    return std::make_pair(std::move(X),V.Values);
}


//plot with boxes
// gp << "plot" << gp.file1d(...) << "with boxes"
template <typename Histogramm_t>
auto make_cols_h(Histogramm_t const & H){
    auto const & _Grid = H.Grid.unhisto();
    std::vector<double> Grid(_Grid.begin(),_Grid.end());
    Grid.resize(Grid.size()-1);
    for(size_t i=0;i<_Grid.size()-1;++i){
        Grid[i] = (_Grid[i]+_Grid[i+1])/2;
    }
    return std::make_pair(std::move(Grid),H.Values);
}



#define debug_return(expr) PVAR((expr)); return expr;

#endif
