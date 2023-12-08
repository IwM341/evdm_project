#ifndef PRNG_HPP
#define PRNG_HPP
#include <omp.h>
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
            uint32_t _uint;
            float _float;
        } m_val;
        m_val._uint = shift_value<23,uint32_t>(value);
        m_val._uint |= (127<<23);
        return m_val._float-1;
    }
    template <typename Value_t>
    inline float int_to_float(Value_t value,MetaPol<genpol::E0I1>,MetaType<float>){
        union {
            uint32_t _uint;
            float _float;
        } m_val;
        m_val._uint = shift_value<23,uint32_t>(value);
        m_val._uint |= (1<<31 | 127<<23);
        return m_val._float+2;
    }
    template <typename Value_t>
    inline double int_to_float(Value_t value,MetaPol<genpol::I0E1>,MetaType<double>){
        union {
            uint64_t _uint;
            double _float;
        } m_val;
        m_val._uint = shift_value<52,uint64_t>(value);
        m_val._uint |= (1023ULL<<52);
        return m_val._float-1;
    }
    template <typename Value_t>
    inline double int_to_float(Value_t value,MetaPol<genpol::E0I1>,MetaType<double>){
        union {
            uint64_t _uint;
            double _float;
        } m_val;
        m_val._uint = shift_value<52,uint64_t>(value);
        m_val._uint |= (1ULL<<63 | 1023ULL<<52);
        return m_val._float+2;
    }

    template <typename T = float,genpol policy = genpol::I0E1>
    struct  xorshift32f{
        mutable uint32_t state;
        xorshift32f(uint32_t seed = 2121885558):state(seed){}
        static inline uint32_t next(uint32_t x){
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
    struct  xorshift64f{
        mutable uint64_t state;
        xorshift64f(uint64_t seed = 818855585854798547):state(seed){}
        static inline uint64_t next(uint64_t x){
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
};

#endif//PRNG_HPP