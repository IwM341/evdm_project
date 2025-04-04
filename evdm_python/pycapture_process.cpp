#include "dynamic_python.hpp"
#include <evdm/core/core_dynamics_capture.hpp>
#include <evdm/utils/prng.hpp>
#include "capture_python.hpp"

pybind11::tuple Py_CaptureProcess(
	Py_Capture& CaptAccum,
	int ptype_in,
	int ptype_out,
	float M_DM,
	float deltaM,
	float NucleiM,
	const Py_ScatterEvent& sc_event,
	float body_halo_v,
	float dm_v_disp,
	float Vhalomax,
	size_t Nmk,
	float r_pow, // the r sistribution
	float weight,// result will be multiplyed by weught
	size_t seed,
	bool ConstrainV
) {
	auto m_compare =
		make_compare_sc_event(sc_event, ptype_in, ptype_out, CaptAccum.events);
	if (!m_compare.empty()) {
		throw pybind11::value_error(m_compare);
	}
	if (seed == 0) {
		seed = std::numeric_limits<size_t>::max();
	}
	auto [sum, dsum] = std::visit([&]<_DISTRIB_TMPL_>(evdm::Distribution<_DISTRIB_PARS_> &mDistrib) {
		typedef decltype(mDistrib.get_grid_vtype()) T;
		pybind11::gil_scoped_release m_gil_release;
		evdm::xorshift<T> G(seed%std::numeric_limits<
			std::conditional_t<std::is_same_v<float,T>,uint_least32_t, uint_least64_t>
		>::max());
		return evdm::Capture(mDistrib, ptype_out, sc_event, G,
			M_DM, deltaM, NucleiM, body_halo_v, dm_v_disp, ConstrainV, Vhalomax, r_pow, Nmk, weight);

	}, CaptAccum.m_distrib);

	auto sce_info = sc_event.unique ?
		scatter_event_info(sc_event.name, ptype_in, ptype_out, sum, Nmk) :
		scatter_event_info();
	{
		pybind11::gil_scoped_acquire m_lock;
		CaptAccum.events.push_back(sce_info);
	}
	
	return pybind11::make_tuple(sum, dsum);
}

void add_pycapture_to_python_module(pybind11::module_& m) {
	namespace py = pybind11;
	m.def("CalcCaptureImpl", Py_CaptureProcess,
		"Calculates capture, add event to capture vector,"
		"returns tuple (capture, sigma)\n\n"
		"Parameters:\n"
		"___________\n"
		"capt_vector : Capture\n\tCapture histogramm.\n"
		"ptype_in : int\n\tindex of input particle type.\n"
		"ptype_out : int\n\tindex of output particle type.\n"
		"m_wimp : float\n\tdm particle mass, GeV.\n"
		"delta : float\n\tdelta mass GeV: output mass - input mass.\n"
		"m_nuc : float\n\tnuclear mass, GeV.\n"
		"sc_event : ScatterEvent.\n"
		"Vbody : float\n\tspeed of body relative to halo.\n"
		"Vdisp : float\n\tdispersion of DM speed in halo.\n"
		"Vmax : float\n\t max DM speed in halo.\n"
		"Nmk : int\n\tnumber of monte-carle steps.\n"
		"r_pow :float\n\t"
		"impact on r distribution: r = (xi)^(r_pow), where xi uniforemly distributed.\n"
		"weight :float\n\t[optional] scale factor, default - 1.\n\t"
		"seed : int\n\t[optional] for random generator \n\t"
		"constrain : bool \n\t if false, DM speed will not truncated",
		py::arg("capt_vector"),
		py::arg("ptype_in"),
		py::arg("ptype_out"),
		py::arg("m_wimp"),
		py::arg("delta"),
		py::arg("m_nuc"),
		py::arg("sc_event"),
		py::arg("Vbody"),
		py::arg("Vdisp"),
		py::arg("Vmax"),
		py::arg("Nmk"),
		py::arg_v("r_pow", 2.0),
		py::arg_v("weight", 1.0),
		py::arg_v("seed", evdm::default_seed),
		py::arg_v("constrain", true)
	);
}