#ifndef R_CONVERSION_HPP
#define R_CONVERSION_HPP

/*
* Conversion from EL density to r density
* Scatter integral with sigma
*/

namespace evdm{
    /// @brief return grid function dN/d^3r(r)
    /// @param H EL distribution fistogramm
    /// @param d3v measure, d3v(dEdL)
    /// @param mu_el measur mu(dEdL)
    /// @return GridFunction r->dN/d3r(r)
    template <typename HistoType,
                typename EL_Functype,
                typename d3v_measure_t,
                typename el_measure_t>
    auto convert_r_density(HistoType const & H,
                            EL_Functype ELf,
                            d3v_measure_t const & d3v,
                            el_measure_t mu_el,double r0 = 0,double r1 = 1)
    {
        //TODO
    }


    /// @brief calculates annihilation matrix
    /// @param Grid EL Grid
    /// @param m_matrix lambda matrix accessor m_matrix(i,j)
    /// @param LE L(E) func
    /// @param d3v measure, d3v(dEdL)
    /// @param mu_el measur mu(dEdL)
    /// @param G random number generator
    /// @param VSigma function, equal to |v1-v2|*sigma(|v1-v2|), where 
    /// sigma --- annihilation cross section, |v1-v2| --- velosity difference  
    /// @return 
    template <typename GridType,
                typename MatrixAcsessor_t,
                typename LE_Functype,
                typename d3v_measure_t,
                typeename el_measure_t,
                typenaem Gen_t,
                typename SigmaFuc_t>
    auto FillAnnMatrix(GridType const & Grid,
                        MatrixAcsessor_t && m_matrix,
                        LE_Functype const & LE,
                        d3v_measure_t d3v,
                        el_measure_t mu_el,
                        Gen_t && G,
                        SigmaFuc_t VSigma)
    {
        //TODO
    }
    
};
#endif//R_CONVERSION_HPP