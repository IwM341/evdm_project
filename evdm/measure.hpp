#ifndef MEASURE_HPP
#define MEASURE_HPP

#include <grob/grid_objects.hpp>

namespace evdm{
    template <typename T>
    inline T dEd1(grob::Point<grob::Rect<T>,grob::Rect<T>> const & m_bin){
        return m_bin.volume();
    }

    template <typename T,typename LE_FuncType>
    inline T dEdL(grob::Point<grob::Rect<T>,grob::Rect<T>> const & m_bin,
                    T const & L0,T const & L1)
    {
        return m_bin.volume()*(L0+L1)/2;
    }
    template <typename T>
    inline T dEdL2(grob::Point<grob::Rect<T>,grob::Rect<T>> const & m_bin,
                    T const & L0,T const & L1){
        auto L0 = LE(m_bin.x<0>().left)
        auto L1 = LE(m_bin.x<0>().right)
        return m_bin.volume()*m_bin.x<1>().center()*(L0*L0+L0*L1+L1*L1)/3;
    }

    /// @brief uniform generatin of (E,L) in bin for dEd1 measure 
    /// @param m_bin 
    /// @param b 
    /// @param G 
    /// @param LE 
    /// @return 
    template <typename T,typename U,typename GenType,typename LE_FuncType && LE>
    auto gen_EL_dEd1(grob::Point<grob::Rect<T>,grob::Rect<T>> const & m_bin,
                                U const & b,
                                GenType && G,LE_FuncType && LE)
    {   
        return [G,b,m_bin,&LE](){
            auto e = m_bin.x<0>().interpol_coeff(G());
            return grob::Point<T,T>{e,
                    m_bin.x<0>().interpol_coeff(G())*LE(e)};
        }
    }

    /// @brief uniform generatin of (E,L) in bin for dEdL measure 
    /// @param m_bin bin [e0,e1]*[l0,l1]
    /// @param b (L1-L0)/(L1+L0)
    /// @param G random generator
    /// @param LE L(e)
    /// @return 
    template <typename T,typename GenType,typename LE_FuncType && LE>
    auto gen_EL_dEdL(grob::Point<grob::Rect<T>,grob::Rect<T>> const & m_bin,
                                U const & b,
                                GenType && G,LE_FuncType && LE)
    {   
        auto b1 = (1-b)/2;
        return [m_bin,&LE,G,b,b1](){
            auto xi = G();
            auto ue =(b >= 1 ? std::sqrt(xi) : xi/(b1 + std::sqrt(b1*b1+b*xi)));
            auto e = m_bin.x<0>().interpol_coeff(ue);
            return grob::Point<T,T> {e,
                    m_bin.x<0>().interpol_coeff(G())*LE(e)};
        }
    }

    /// @brief uniform generatin of (E,L) in bin for dEdL2 measure 
    /// @param m_bin bin [e0,e1]*[l0,l1]
    /// @param b (L1-L0)/(L1+L0)
    /// @param G random generator
    /// @param LE L(e)
    /// @return 
    template <typename T,typename GenType,typename LE_FuncType && LE>
    auto gen_EL_dEdL2(
        grob::Point<grob::Rect<T>,grob::Rect<T>> const &  m_bin,
        U const & b,
        GenType && G,LE_FuncType && LE)
    {
       
        auto b1 = (1-b)/2;
        auto l0 = m_bin.x<1>().left;
        auto l1 = m_bin.x<1>().right;
        auto l02 =l0*l0;
        auto l12 =l1*l1; 
        auto weight = 1/(4*b*b/6+b1);
        return [G,b1,b,&LE,weight](){
            auto xi = G();
            auto ue =(b >= 1 ? std::sqrt(xi) : xi/(b1 + std::sqrt(b1*b1+b*xi)));
            auto e = m_bin.x<0>().interpol_coeff(ue);
            auto lx = G();
            return MCResult<grob::Point<T,T>,T>{{e,
                    std::sqrt(l02*(1-lx)+l12*x)*LE(e)},
                        (b*ue+b1)*weight
                    };
        }
    }
}

#endif//MEASURE_HPP