#include "dynamic_python.hpp"
#include <evdm/core/core_dynamics_scatter.hpp>
#include <evdm/utils/prng.hpp>
#include "progress_log.hpp"


void Py_ScatterProcess(
	Py_Matrix& ScatterAccum,
	int ptype_in, int ptype_out,
	float M_DM,
	float deltaM,
	float NucleiM,
	const Py_ScatterEvent& sc_event,
	size_t Nmk,
	size_t Nmk_per_traj,
	bool count_evap,
	float weight,
	pybind11::handle update_function
) {
	auto m_compare =
		make_compare_sc_event(sc_event, ptype_in, ptype_out, ScatterAccum.events);
	if (!m_compare.empty()) {
		throw pybind11::value_error(m_compare);
	}

	auto ProgFunc = make_progress_func(update_function);
	std::visit([&]
		<class Vt, class Bt, class Gvt, evdm::GridEL_type G_tp>
		(Matrix_Pair_Inst<Vt, Bt, Gvt, G_tp> &m_matrix)->void {
		typedef Gvt T;
		evdm::Distribution<Vt, Bt, Gvt, G_tp> & m_evp = m_matrix.second;
		evdm::GridMatrix<Vt, Bt, Gvt, G_tp> & m_mat = m_matrix.first;

		evdm::Scatter(
			m_mat, m_evp, count_evap, ptype_in, ptype_out,
			sc_event, evdm::xorshift<T>{}, evdm::measure_dEdL{},
			M_DM, deltaM, NucleiM, Nmk, Nmk_per_traj, weight, ProgFunc
		);
	}, ScatterAccum.m_matrix);


	auto sce_info = sc_event.unique ?
		scatter_event_info(sc_event.name, ptype_in, ptype_out, 0, Nmk) :
		scatter_event_info();
	ScatterAccum.events.push_back(sce_info);
}



void add_pyscatter_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;

	m.def("CalcScatterImpl", &Py_ScatterProcess,
		"Calculates scatter matrix part, add event to matrix class\n\n"
		"Parameters:\n"
		"___________\n"
		"matrix : Matrix\n\tScatter matrix histo.\n"
		"ptype_in : int\n\tindex of input particle type.\n"
		"ptype_out : int\n\tindex of output particle type.\n"
		"m_wimp : float\n\tdm particle mass, GeV.\n"
		"delta : float\n\tdelta mass GeV: output mass - input mass.\n"
		"m_nuc : float\n\tnuclear mass, GeV.\n"
		"sc_event : ScatterEvent.\n"
		"Nmk : int\n\tnumber of monte-carle steps.\n"
		"Nmk_traj: int\n\tnumber of monte-carle steps on each trajectory.\n"
		"weight : float\n\t[optional] scale factor, default - 1.\n"
		"bar : object\n\t[optional] progress bar update function.",
		py::arg("matrix"),
		py::arg("ptype_in"),
		py::arg("ptype_out"),
		py::arg("m_wimp"),
		py::arg("delta"),
		py::arg("m_nuc"),
		py::arg("sc_event"),
		py::arg("Nmk"),
		py::arg_v("Nmk_traj",10),
		py::arg_v("count_evap",false),
		py::arg_v("weight",1),
		py::arg_v("bar",py::none())
	);
}

