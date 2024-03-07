#ifndef HISTO_NORM_HPP
#define HISTO_NORM_HPP
/*
*   Comparing two distributions with different grid,
*   with different norm
*/
#include "utils/mc.hpp"
namespace evdm{

    /// @brief calculates L_p norm of difference between two distribs
    /// @param H1 histo 1
    /// @param H2 histo 2
    /// @param p norm param
    /// @param mu measure of dEdL mu : Point<Rect<U>,Rect<U>>
    /// @return 
    template <typename HistoType1,typename HistoType2,typename BinNorm_t,typename T>
    auto histo_norm(HistoType1 const & H1,HistoType2 const & H2,T const & p,BinNorm_t const & mu){
        double sum = 0;

        if(H1.Grid.grid().size() != H2.Grid.grid().size())
        {
            throw std::range_error(
                "error comparing 2 distributions:"
                "different number of particle types");
        }
        size_t p_types = H1.Grid.grid().size();
        auto & EL_Grid1 = H1.Grid.inner(0);
        auto & EL_Grid2 = H2.Grid.inner(0);

        for(size_t p_type =0;p_type <p_types;++p_type){
            for(auto it1 = EL_Grid1.begin();it1 != EL_Grid2.end();++it1){
                auto dEdL1 = *it1;
                auto H1_v = H1[grob::make_MI_rec(p_type,it1.index())];
                auto mes1 = mu(dEdL1);
                auto rho1 = H1_v/mes1;
                size_t i2_0 = EL_Grid2.grid().pos(std::get<0>(dEdL1).left);
                size_t i2_1 = EL_Grid2.grid().pos(std::get<0>(dEdL1).right) + 1;
                for(size_t i2 = i2_0; i2 < i2_1;++i2){
                    size_t j2_0 = EL_Grid2.inner(i2).pos(std::get<1>(dEdL1).left);
                    size_t j2_1 = EL_Grid2.inner(i2).pos(std::get<1>(dEdL1).right)+1;
                    for(size_t j2 = j2_0; j2 < j2_1;++j2){
                        auto m_dEdL2 = EL_Grid2[grob::make_MI(i2,j2)];
                        auto H2_v = H2[grob::make_MI(p_type,i2,j2)];
                        auto mes2 = mu(m_dEdL2);

                        auto rho2 = H2_v/mes2;

                        auto [_m_cube, is_intersec] = grob::intersect(dEdL1, m_dEdL2);
                        if(is_intersec){
                            sum += std::pow(std::abs(rho1 - rho2), p) *
                                mu(_m_cube);
                        }
                    }        
                }
            }
        }
        return std::pow(sum,1.0/p);
    }

    

};

#endif//HISTO_NORM_HPP