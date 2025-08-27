#pragma once
#include "core_matrix.hpp"
#include "../measure.hpp"


namespace evdm {
	template <typename T>
	struct GridProjector {
		measure_dEpdlq<T> m_mes;
		/// <summary>
		/// 
		/// </summary>
		/// <param name="p">in measure dxE^p</param>
		/// <param name="q">in measure dl^q</param>
		GridProjector(T p, T q, T Emin) : m_mes(p, q, Emin) {}

		template <typename GridEL_vt1, typename GridEL_vt2,
			GridEL_type grid_type1, GridEL_type grid_type2
		> SpMatrix_t<T> operator () (
			const GridEL< GridEL_vt1, grid_type1>& G_out,
			const GridEL< GridEL_vt2, grid_type2>& G_in) const {
			

            auto p_types = G_out.ptypes();
            auto& grid_out = G_out.inner(0);
            auto& grid_in = G_in.inner(0);

            
            size_t NGrid_in = grid_in.size();

            size_t NGrid_out = grid_out.size();

            std::vector<std::vector<SpTriplet_t<T>>> Triplets(p_types* NGrid_in);

            #pragma omp parallel for
            for (int Iin = 0; Iin < NGrid_in; ++Iin) {
                auto [i_in,j_in] = grid_in.FromLinear(Iin);
                auto m_bin_in = grid_in[{i_in, j_in}];
                auto VolumeIn = m_mes(m_bin_in);

                size_t i2_0 = grid_out.grid().pos(std::get<0>(m_bin_in).left);
                size_t i2_1 = grid_out.grid().pos(std::get<0>(m_bin_in).right) + 1;

                for (size_t i2 = i2_0; i2 < i2_1; ++i2) {
                    size_t j2_0 = grid_out.inner(i2).pos(std::get<1>(m_bin_in).left);
                    size_t j2_1 = grid_out.inner(i2).pos(std::get<1>(m_bin_in).right) + 1;

                    for (size_t j2 = j2_0; j2 < j2_1; ++j2) {
                        size_t Iout = grid_out.LinearIndex({ i2,j2 });
                        auto m_bin_out = grid_out[{i2,j2}];
                        auto [bin_intersec,is_intersec] = grob::intersect(
                            m_bin_in, m_bin_out
                        );
                        
                        if (is_intersec) {
                            auto weight = 
                                VolumeIn > 0 ? 
                                (T)m_mes(bin_intersec) / VolumeIn : 
                                (T)0;
                            for (size_t ptype = 0; ptype < p_types; ++ptype) {
                                size_t p_fac_in = ptype * NGrid_in;
                                size_t p_fac_out = ptype * NGrid_out;
                                Triplets[Iin + p_fac_in].push_back(
                                    SpTriplet_t<T>(
                                        Iout + p_fac_out, 
                                        Iin + p_fac_in, 
                                        weight
                                    )
                                );
                            }
                        }
                        
                    }
                }
            }
            std::vector<SpTriplet_t<T>>
                AllTriplets(flatten(Triplets));

            SpMatrix_t<T> Ret(p_types * NGrid_out, p_types * NGrid_in);
            Ret.setFromTriplets(AllTriplets.begin(), AllTriplets.end());
            return Ret;
		}
	};
};
