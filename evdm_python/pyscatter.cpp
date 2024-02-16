#include "dynamic_python.hpp"
#include <evdm/core_dynamics.hpp>
#include <evdm/utils/prng.hpp>
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


pybind11::tuple Py_CaptureProccess(
	Py_Capture& CaptAccum,
	int ptype,
	float M_DM,
	float deltaM,
	float NucleiM,
	const Py_ScatterEvent& sc_event,
	float body_halo_v,
	float dm_v_disp,
	size_t Nmk,
	std::string gen_dtype,
	float r_pow, // the r sistribution
	float weight// result will be multiplyed by weught
	)
{
	auto [sum,dsum] = std::visit([&](auto& mDistrib) {
			if(gen_dtype == "float"){
				return evdm::Capture(mDistrib, ptype, sc_event, evdm::xorshift32f<float>{},
					M_DM, deltaM, NucleiM, body_halo_v, dm_v_disp, r_pow, Nmk, weight);
			}
			else if (gen_dtype == "double") {
				return evdm::Capture(mDistrib, ptype, sc_event, evdm::xorshift64f<double>{},
					M_DM, deltaM, NucleiM, body_halo_v, dm_v_disp, r_pow, Nmk, weight);
			}
		}, CaptAccum.m_distrib);
	return pybind11::make_tuple(sum, dsum);
}

void add_scatter_funcs_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;
	m.def("Capture", Py_CaptureProccess,
		"Calculates capture, add event to capture vector,"
		" returns tuple (capture,mk sigma)\n"
		"capt_vector -- Capture histogramm\n"
		"ptype -- index of particle type\n"
		"m_wimp -- dm particle mass, GeV\n"
		"delta -- delta mass GeV: output mass - input mass\n"
		"m_nuc -- nuclear mass, GeV\n"
		"sc_event -- ScatterEvent instance\n"
		"Vbody -- speed of body relative to halo\n"
		"Vdisp -- dispersion of DM speed in halo\n"
		"Nmk -- number of monte-carle steps\n"
		"gen_vtype -- value type of prng: 'float' or 'double'\n"
		"r_pow -- impact on r distribution: r = (xi)^(r_pow), where xi uniforemly distributed\n"
		"weight -- additional scale factor, default - 1",
		py::arg("capt_vector"),
		py::arg("ptype"),
		py::arg("m_wimp"),
		py::arg("delta"),
		py::arg("m_nuc"),
		py::arg("sc_event"),
		py::arg("Vbody"),
		py::arg("Vdisp"),
		py::arg("Nmk"),
		py::arg_v("gen_vtype","float"),
		py::arg_v("r_pow",2.0),
		py::arg_v("weight",1.0)
	);
}


typedef Py_ScatterFactor(*Qexp_0_t)(float, bool, pybind11::array_t<float>);

void Py_ScatterFactor::add_to_python_module(pybind11::module_& m) {
	namespace py = pybind11;
	py::class_<Py_ScatterFactor>(m, "ScatterFactor")
		.def("__repr__", [](const Py_ScatterFactor& sf) {return sf.repr(); });
	m.def("qexp_factor", static_cast<Qexp_1_t>( &qexp_factor),
		"create exponential form factor exp(-2y)(p0_i y^i + pv_i y^i*v_{perp}^2  )/y^t\n"
		"where y = b^2*q^2/4\n"
		"if y_inv is false than t = 0 and y^t = 1, else t = 1\n"
		"P_0 -- array like with coefficients of y^i\n"
		"P_V -- optional, array like with coefficints of y^i v_perp^2",
		py::arg("b"),
		py::arg("y_inv"),
		py::arg("P_0"),
		py::arg("P_V")
	)
	.def("qexp_factor", static_cast<Qexp_0_t>(&qexp_factor),"",
		py::arg("b"),
		py::arg("y_inv"),
		py::arg("P_0"))
	.def("__repr__", [](const Py_ScatterFactor& sf) {
			return sf.repr();
		})
	.def("__str__", [](const Py_ScatterFactor& sf) {
			return sf.repr();
		});
}

Py_ScatterEvent::Py_ScatterEvent(
	pybind11::array n_e, 
	Py_ScatterFactor const &_factor, 
	std::string_view name,
	bool unique
): evdm::ScatterEvent(n_e.size(),(float * ) n_e.data(),_factor.to_ff()),
	_name(name), unique(unique)
{

}

void Py_ScatterEvent::add_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;
	py::class_<Py_ScatterEvent>(m, "ScatterEvent")
		.def(
			py::init<
				pybind11::array, 
				Py_ScatterFactor const&,
				std::string_view
			>(),
			"creating ScatterEvent\n"
			"n_e -- relative to n_p concentration of event,"
			" where n_p - concentration of nuclons\n"
			"sf -- ScatterFactor instance\n"
			"unique -- if true, sum of capture with same event names wil be blocked",
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
