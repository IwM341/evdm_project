#include "core_distrib.hpp"
#include "../dynamics/dynamics.hpp"
#include "../utils/discret_generator.hpp"
#include "../grid_variants.hpp"
#include <span>


namespace evdm {

	template <typename Bt, typename Gt, typename St, evdm::GridEL_type gtype>
	Distribution<St, Bt, Gt, gtype> put_points(
		EL_Grid<Bt, Gt, gtype> const& _m_grid, 
		std::span<const StateEL<St>> m_states, 
		same_t<St> total_count
	) {
		GridEL<Gt, gtype> const & m_grid  = *_m_grid.Grid;
		St factor = total_count / m_states.size();
		Eigen::VectorX<St> values(m_grid.size());
		values.setZero();
		size_t ptypes = m_grid.ptypes();
		for (int i = 0; i < m_states.size(); ++i) {
			if (m_states[i].ptype < ptypes, m_states[i].e < 0){
				auto [ptype, e, l] = m_states[i];
				auto m_index = m_grid.pos(ptype, e, l);
				int index = m_grid.LinearIndex(
					m_index
				);
				values[index] += factor;
			}
		}
		return Distribution<St, Bt, Gt, gtype>(_m_grid, std::move(values));
	}

	template<  typename T,typename Bt, typename Gt, evdm::GridEL_type gtype, typename St,typename Gen_t>
	void sample_distrib(Distribution<T, Bt, Gt, gtype> const& distrib,std::span< StateEL<St>> states,Gen_t && G,same_t<Gt> mes_p, same_t<Gt> mes_q) {
		auto const & V = distrib.values();
		std::vector<T> m_sample(V.begin(),V.end());
		for(int i=1;i<m_sample.size();++i){
			m_sample[i] += m_sample[i-1];
		}
		 
		evdm::measure_dEpdlq<Gt> m_mes(mes_p,mes_q,distrib.body().Phi.Values[0]);
		GridEL<Gt, gtype> const & m_grid = distrib.grid();
		for(StateEL<St> & state : states){
			T xi = G()*m_sample.back();
			size_t idx = IndexGenBigHelper::gen(m_sample.begin(),m_sample.end(),xi);
			auto [ptype,e0e1,l0l1] = m_grid[m_grid.FromLinear(idx)];
			using namespace grob::literals;
			bin_dedl_t<Gt> dxdy = m_mes.toXY(grob::make_point(e0e1,l0l1));
			auto [e,l] = m_mes.fromXY(dxdy[0_c].reduction(G()),dxdy[1_c].reduction(G()));
			state.ptype = ptype;
			state.e = e;
			state.l = l;
		}
	}
};


