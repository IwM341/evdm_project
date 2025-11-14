#pragma once

#include <evdm/core/core_dynamics_scatter.hpp>
#include <evdm/utils/prng.hpp>
#include <evdm/measure.hpp>
#include "../matrix_python.hpp"

/*template <
	class Nmk_t,
	class ThermGen_t,
	class Vt, class Bt, class Gvt,
	evdm::GridEL_type G_tp,
	int algol
>
void ScatterImpl_inner(
	std::integral_constant<int, algol>,
	evdm::GridMatrix< Vt,Bt,Gvt, G_tp> & m_mat,
	evdm::Distribution< Vt, Bt, Gvt, G_tp>& m_evp,
	double zero_val,int ptype_in, int ptype_out,size_t seed,
	evdm::ScatterEvent const& sc_event, ThermGen_t m_therm_gen,
	evdm::measure_dEpdlq< Gvt> m_el_measure,
	float M_DM,float deltaM, float NucleiM,
	Nmk_t const & Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
) {
	evdm::xorshift<Gvt> G(seed);
	evdm::Scatter(std::integral_constant<int, algol>{},
		m_mat, m_evp, zero_val, ptype_in, ptype_out,
		sc_event, G, m_therm_gen, m_el_measure,
		M_DM, deltaM, NucleiM, Nmk, Nmk_per_traj, weight, ProgFunc
	);

}*/

struct Nmk_vector {
	std::variant<size_t, std::vector<size_t>> Nmks;
	size_t operator [](size_t i) const {
		return std::visit([i](auto const& Nmks)->size_t {
			if constexpr (std::is_same_v <std::decay_t<decltype(Nmks)>, size_t>) {
				return Nmks;
			}
			else {
				if (i > Nmks.size()) {
					throw std::range_error("i > Nmks.size()");
				}
				return Nmks[i];
			}
		}, Nmks);
	}
	template <typename T>
	Nmk_vector(T&& value):Nmks(std::forward<T>(value)){}

	Nmk_vector(){}
};


struct ScatterImpl_Args {
	Py_Matrix& m_mat;
	double zero_val; int ptype_in; int ptype_out; size_t seed;
	evdm::ScatterEvent const& sc_event; ScatterMethodVariant_t therm_gen_var;
	std::tuple<double, double> m_el_measure_pair;
	double M_DM; double deltaM; double NucleiM;
	Nmk_vector const& Nmk; size_t Nmk_per_traj;
	double weight; evdm::ScatterImplExtra extra_args;
	
	template <int alg_i>
	void apply_impl(std::integral_constant<int, alg_i> algol) {
		
		std::visit([&]<
			class Vt, class Bt, class Gvt,
			evdm::GridEL_type G_tp
		>(Matrix_Pair_Inst<Vt, Bt, Gvt, G_tp> &m_matrix) {
			evdm::Distribution<Vt,Bt,Gvt, G_tp>& m_evp = m_matrix.second;
			evdm::GridMatrix<Vt, Bt, Gvt, G_tp>& m_mat = m_matrix.first;
			evdm::xorshift<Gvt> G(seed);
			auto Emin = m_evp.grid().inner(0).grid()[0].left;
			auto [p, q] = m_el_measure_pair;
			evdm::measure_dEpdlq<Gvt> m_el_measure(p, q, Emin);
			evdm::ThermGaussGenerator_Variant<
				ScatterMethodVariant_t> therm_gen{ therm_gen_var };
			evdm::Scatter(std::integral_constant<int, algol>{},
				m_mat, m_evp, zero_val, ptype_in, ptype_out,
				sc_event, G, therm_gen, m_el_measure,
				M_DM, deltaM, NucleiM, Nmk, Nmk_per_traj, weight, extra_args
			);
		}, m_mat.m_matrix);
		
	}
	void apply(std::integral_constant<int, evdm::AlgolNaive>);
	void apply(std::integral_constant<int, evdm::AlgolShift>);
	void apply(std::integral_constant<int, evdm::AlgolDiffuse>);
};

