#include "dynamic_python.hpp"


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

Py_ScatterFactor helm_factor(float R,float s2,float cns_fac) {
	return Py_ScatterFactor(evdm::BesselFormFactor(R, s2, cns_fac));
}


std::string make_compare_sc_event(
	const Py_ScatterEvent& sc_event,
	int ptype_in, int ptype_out,
	const std::vector< scatter_event_info>& events
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
			py::arg_v("v2T", 0)
		);
	m.def("qexp_factor", static_cast<Qexp_1_t>(&qexp_factor),
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
	.def("helm_factor", helm_factor,
		"create helm form fractor\n"
		"FF(q^2) = CF*(Bessels[qR])^2*exp(-q^2*s2)"
		"Parameters:\n",
		py::arg("R"), py::arg("s2"), py::arg("CF"))
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
		Conc.resize({ 2 });
		Conc.mutable_at(1) = Conc.at(0);
	}
	return evdm::ScatterEvent(Conc.size(), (float*)Conc.data(), _factor.to_ff());
}

Py_ScatterEvent::Py_ScatterEvent(
	pybind11::handle n_e,
	Py_ScatterFactor const& _factor,
	const char* name,
	bool unique
) : evdm::ScatterEvent(Py_MakeScatterEvent(n_e, _factor)),
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
			py::arg_v("name", "__unnamed__"),
			py::arg_v("unique", false)
		)
		.def("__repr__", [](const Py_ScatterEvent& event) {
		return event.repr();
			})
		.def("__str__", [](const Py_ScatterEvent& event) {
				return event.repr();
			});
}

