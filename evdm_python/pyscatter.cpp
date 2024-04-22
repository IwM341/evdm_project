#include "dynamic_python.hpp"
#include <evdm/core_dynamics.hpp>
#include <evdm/utils/prng.hpp>
#include "progress_log.hpp"

Py_ScatterFactor qexp_factor(float b, bool y_inv,
	pybind11::array_t<float> P_0_coeffs, pybind11::array_t<float> P_V_coeffs)
{
	auto P_0_arr = grob::vector_view<const float>(
			P_0_coeffs.data(), pybind11::len(P_0_coeffs)
		);
	auto P_V_arr = grob::vector_view<const float>(
			P_V_coeffs.data(), pybind11::len(P_V_coeffs)
		);
	size_t _msize = std::max(P_0_arr.size(), P_V_arr.size());
	return Py_ScatterFactor(y_inv, b, P_0_arr, P_V_arr);
}
typedef Py_ScatterFactor(*Qexp_1_t)(float, bool, pybind11::array_t<float>, pybind11::array_t<float>);

Py_ScatterFactor qexp_factor(float b, bool y_inv,
	pybind11::array_t<float> P_0_coeffs)
{
	auto P_0_arr = grob::vector_view<const float>(
		P_0_coeffs.data(), pybind11::len(P_0_coeffs)
		);

	size_t _msize = P_0_arr.size();
	return Py_ScatterFactor(y_inv, b, P_0_arr);
}

Py_ScatterFactor func_factor(pybind11::object _function) {
	float (*func)(float, float);
	func = reinterpret_cast<float (*)(float, float)> (
		_function.attr("address").cast<size_t>()
		);
	Py_ScatterFactor sf(func);
	sf._opt_functur = _function;
	return sf;
}


std::string make_compare_sc_event(
	const Py_ScatterEvent& sc_event,
	int ptype_in,int ptype_out,
	const std::vector< scatter_event_info> & events
) {
	for (auto const& ev : events) {
		if (ev.unique && sc_event.unique &&
			ev.name == sc_event.name &&
			ev.ptype_in == ptype_in && ev.ptype_out == ptype_out) {
			std::string exc_str = "error: scatter event of name: '" + ev.name +
				"' and ptype_in/out = " +
				std::to_string(ptype_in) + "/" + std::to_string(ptype_in) +
				" is alredy considered, avoid double summation";
			return exc_str;
		}
	}
	return "";
}

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
	size_t Nmk,
	float r_pow, // the r sistribution
	float weight// result will be multiplyed by weught
	)
{
	auto m_compare = 
		make_compare_sc_event(sc_event, ptype_in, ptype_out, CaptAccum.events);
	if (!m_compare.empty()) {
		throw pybind11::value_error(m_compare);
	}
	auto [sum,dsum] = std::visit([&]_DISTRIB_TMPL_(evdm::Distribution<_DISTRIB_PARS_> & mDistrib) {
			typedef decltype(mDistrib.get_grid_vtype()) T;
			return evdm::Capture(mDistrib, ptype_out, sc_event, evdm::xorshift<T>{},
				M_DM, deltaM, NucleiM, body_halo_v, dm_v_disp, r_pow, Nmk, weight);
			
		}, CaptAccum.m_distrib);

	auto sce_info = sc_event.unique ?
		scatter_event_info(sc_event.name, ptype_in, ptype_out, sum, Nmk) :
		scatter_event_info();
	CaptAccum.events.push_back(sce_info);
	return pybind11::make_tuple(sum, dsum);
}

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
		<class Vt, class Bt, class Gvt,evdm::GridEL_type G_tp>
		(Matrix_Pair_Inst<Vt,Bt,Gvt,G_tp> & m_matrix)->void {
		typedef Gvt T;
		evdm::Scatter(
			m_matrix.first, m_matrix.second,count_evap, ptype_in, ptype_out,
			sc_event, evdm::xorshift<T>{}, evdm::measure_dEdL{},
			M_DM, deltaM, NucleiM, Nmk, Nmk_per_traj, weight, ProgFunc
		);
	}, ScatterAccum.m_matrix);


	auto sce_info = sc_event.unique ?
		scatter_event_info(sc_event.name, ptype_in, ptype_out, 0, Nmk) :
		scatter_event_info();
	ScatterAccum.events.push_back(sce_info);
}


void add_scatter_funcs_to_python_module(pybind11::module_& m)
{
	Py_ScatterFactor::add_to_python_module(m);
	Py_ScatterEvent::add_to_python_module(m);
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
		"Nmk : int\n\tnumber of monte-carle steps.\n"
		"r_pow :float\n\t"
		"impact on r distribution: r = (xi)^(r_pow), where xi uniforemly distributed.\n"
		"weight :float\n\t[optional] scale factor, default - 1.",
		py::arg("capt_vector"),
		py::arg("ptype_in"),
		py::arg("ptype_out"),
		py::arg("m_wimp"),
		py::arg("delta"),
		py::arg("m_nuc"),
		py::arg("sc_event"),
		py::arg("Vbody"),
		py::arg("Vdisp"),
		py::arg("Nmk"),
		py::arg_v("r_pow",2.0),
		py::arg_v("weight",1.0)
	);
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


typedef Py_ScatterFactor(*Qexp_0_t)(float, bool, pybind11::array_t<float>);

void Py_ScatterFactor::add_to_python_module(pybind11::module_& m) {
	namespace py = pybind11;
	py::class_<Py_ScatterFactor>(m, "ScatterFactor")
		.def("__repr__", [](const Py_ScatterFactor& sf) {return sf.repr(); })
		.def("__str__", [](const Py_ScatterFactor& sf) {return sf.repr(); })
		.def("__call__", 
			[](const Py_ScatterFactor& sf, float y, float v2T)
			{return sf.eval(y, v2T); },
			"evaluates form factor for demonstration",
			py::arg("y"),
			py::arg_v("v2T",0)
			);
	m.def("qexp_factor", static_cast<Qexp_1_t>( &qexp_factor),
		"create exponential form factor exp(-2y)(p0_i y^i + pv_i y^i*v_{perp}^2  )/y^t\n"
		"where y = b^2*q^2/4\n"
		"Parameters:\n"
		"___________\n"
		"b : float\n\t size of nuclei in GeVn\n"
		"y_inv : bool\n\tif false than t = 0 and y^t = 1, else t = 1\n"
		"P_0 : array\n\tcoefficients of y^i\n"
		"P_V : array\n\toptional coefficints of y^i v_perp^2",
		py::arg("b"),
		py::arg("y_inv"),
		py::arg("P_0"),
		py::arg("P_V")
	)
	.def("qexp_factor", static_cast<Qexp_0_t>(&qexp_factor),
		"create exponential form factor exp(-2y)(p0_i y^i)/y^t\n"
		"where y = b^2*q^2/4\n"
		"Parameters:\n"
		"___________\n"
		"b : float\n\t size of nuclei in GeVn\n"
		"y_inv : bool\n\tif false than t = 0 and y^t = 1, else t = 1\n"
		"P_0 : array\n\tcoefficients of y^i",
		py::arg("b"),
		py::arg("y_inv"),
		py::arg("P_0"))
	.def("func_factor", func_factor,
		"create elastic form fractor from function\n"
		"Input should contain pointer to funcion with signature float ScatterFunc(float q_2,float v2T)\n"
		"where q_2 = q^2 - transferred momentum in GeV,v2T - squared norm of inelastic transfer velocity\n"
		"Parameters:\n"
		"___________\n"
		"func : function\n\tfloat ScatterFunc(float q_2,float v2T)",
		py::arg("func"))
	.def("__repr__", [](const Py_ScatterFactor& sf) {
			return sf.repr();
		})
	.def("__str__", [](const Py_ScatterFactor& sf) {
			return sf.repr();
		});
}

evdm::ScatterEvent Py_MakeScatterEvent(
	pybind11::handle n_e, 
	Py_ScatterFactor const& _factor) 
{
	pybind11::array_t<float> Conc = n_e.cast<pybind11::array_t<float>>();
	if (Conc.size() < 2) {
		Conc.resize({2});
		Conc.mutable_at(1) = Conc.at(0);
	}
	return evdm::ScatterEvent(Conc.size(), (float*)Conc.data(), _factor.to_ff());
}

Py_ScatterEvent::Py_ScatterEvent(
	pybind11::handle n_e, 
	Py_ScatterFactor const &_factor, 
	const char* name,
	bool unique
): evdm::ScatterEvent(Py_MakeScatterEvent(n_e,_factor)),
	name(name), unique(unique)
{

}

void Py_ScatterEvent::add_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;
	py::class_<Py_ScatterEvent>(m, "ScatterEvent")
		.def(
			py::init<
				pybind11::handle, 
				Py_ScatterFactor const&,
				const char*,
				bool
			>(),
			"creating ScatterEvent.\n\n"
			"Parameters:\n"
			"___________\n"
			"n_e : array \n\t relative to n_p concentration of targets,"
			" where n_p - avarage concentration of nuclons.\n"
			"sf : ScatterFactor\n\tscatter factor class instance.\n"
			"name : string\n\toptional name of event.\n"
			"unique : bool\n\tif true, sum of capture with same event names will be blocked.",
			py::arg("n_e"),
			py::arg("sf"),
			py::arg_v("name","__unnamed__"),
			py::arg_v("unique", false)
		)
		.def("__repr__", [](const Py_ScatterEvent& event) {
				return event.repr();
			})
		.def("__str__", [](const Py_ScatterEvent& event) {
				return event.repr();
			});
}
