#ifndef PRNG_HPP
#define PRNG_HPP

#include <stdint.h>
#include <type_traits>
#include <limits>

namespace evdm{

    enum class genpol{
        //I0I1, //include 0, include 1
        I0E1, //include 0, exclude 1 
        E0I1 //exclude 0, include 1 
        //E0E1  //exclude 0, exclude 1 
    };

    template <genpol P>
    struct MetaPol{
        constexpr static genpol value = P; 
    };
    
    template <typename T>
    struct MetaType{
        typedef T type;
    };

    template <size_t _mantiss,typename OutType,typename T>
    inline constexpr OutType shift_value(T value){
        constexpr size_t _len = std::numeric_limits<T>::digits;
        if constexpr(_len >= _mantiss){
            return ((OutType)value) >> (_len-_mantiss);
        } else {
            return ((OutType)value) << (_mantiss - _len);
        }
    }
    template <typename Value_t>
    inline float int_to_float(Value_t value,MetaPol<genpol::I0E1>,MetaType<float>){
        union {
            uint_least32_t _uint;
            float _float;
        } m_val;
        m_val._uint = shift_value<23,uint_least32_t>(value);
        m_val._uint |= (127<<23);
        return m_val._float-1;
    }
    template <typename Value_t>
    inline float int_to_float(Value_t value,MetaPol<genpol::E0I1>,MetaType<float>){
        union {
            uint_least32_t _uint;
            float _float;
        } m_val;
        m_val._uint = shift_value<23,uint_least32_t>(value);
        m_val._uint |= (1<<31 | 127<<23);
        return m_val._float+2;
    }
    template <typename Value_t>
    inline double int_to_float(Value_t value,MetaPol<genpol::I0E1>,MetaType<double>){
        union {
            uint_least64_t _uint;
            double _float;
        } m_val;
        m_val._uint = shift_value<52,uint_least64_t>(value);
        m_val._uint |= (1023ULL<<52);
        return m_val._float-1;
    }
    template <typename Value_t>
    inline double int_to_float(Value_t value,MetaPol<genpol::E0I1>,MetaType<double>){
        union {
            uint_least64_t _uint;
            double _float;
        } m_val;
        m_val._uint = shift_value<52,uint_least64_t>(value);
        m_val._uint |= (1ULL<<63 | 1023ULL<<52);
        return m_val._float+2;
    }

    template <typename T,genpol policy>
    struct _bounds_base {
        typedef T value_type;
        constexpr static genpol bounds = policy;
    };

    template <typename T = float,genpol policy = genpol::I0E1>
    struct  xorshift32f : public _bounds_base <T,policy>{
        mutable uint_least32_t state;

        xorshift32f(uint_least32_t seed = 2121885558):state(seed){}
        static inline uint_least32_t next(uint_least32_t x){
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            return x;
        }
        inline T operator ()() const{
            return int_to_float(state = next(state),
                    MetaPol<policy>{},
                    MetaType<T>{});
        }
        void thread_seed(size_t thread_num = 0) const{
            state = next(next(thread_num + 1));
        }
    };

    template <typename T = float,genpol policy = genpol::I0E1>
    struct  xorshift64f : public _bounds_base <T,policy> {
        mutable uint_least64_t state;
        xorshift64f(uint_least64_t seed = 818855585854798547):state(seed){}
        static inline uint_least64_t next(uint_least64_t x){
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            return x;
        }
        inline T operator ()() const{
            return int_to_float(state = next(state),
                    MetaPol<policy>{},
                    MetaType<T>{});
        }
        void thread_seed(size_t thread_num = 0) const{
            state = next(next(thread_num + 1));
        }
    };


    template <typename T, genpol policy = genpol::I0E1>
    struct universal_xorshift_t : 
    public std::conditional_t<
        std::is_same_v<float, T>,
        xorshift32f<T, policy>,
        xorshift64f<T, policy>
    > {};

    template <typename T, genpol policy = genpol::I0E1>
    using xorshift = universal_xorshift_t<T, policy>;

    template <typename Gen_t>
    inline auto E0I1_G(Gen_t&& G) {
        if constexpr (G.bounds == genpol::E0I1)
            return G();
        else if (G.bounds == genpol::I0E1)
            return 1-G();
        else
            static_assert(true, "incorrect generator bounds"); 
    }
};

#endif//PRNG_HPP