
#include "dist_sampler.hpp"

pybind11::array init_samples_impl(type_fd_var, Py_Distribution const & m_distrib, size_t N, size_t seed,double p,double q) {
	return std::visit(
		[&]<class Vt, _DISTRIB_TMPL_>(
			std::type_identity<Vt>, 
			evdm::Distribution<_DISTRIB_PARS_> const& m_d
		)->pybind11::array {
		pybind11::array_t<evdm::StateEL<Vt>> m_array(N);
		evdm::sample_distrib(m_d, std::span< evdm::StateEL<Vt>>(m_array.mutable_data(),m_array.size()), evdm::xorshift<_T3>(seed),p,q);
		return m_array;
	},type_fd_var{}, m_distrib.m_distrib);
}

Py_Distribution put_points_impl(type_fd_var, Py_EL_Grid const& Grid, pybind11::array m_states, double total_count) {
	
	return std::visit(
		[&]<class T>(std::type_identity<T>, auto const& m_grid)->Py_Distribution 
		{
		pybind11::array_t<evdm::StateEL<T>> m_states_t = m_states;
		return evdm::put_points(m_grid, std::span<const evdm::StateEL<T>>(m_states_t.data(), m_states_t.size()), total_count);
		}, type_fd_var{}, Grid.m_grid
	);
}