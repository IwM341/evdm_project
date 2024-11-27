#include "capture_python.hpp"
#include <pybind11/stl.h>



void scatter_event_info::add_to_python_module(pybind11::module_& m) {
	namespace py = pybind11;
	py::class_<scatter_event_info>(m, "scatter_event_info")
		.def("__repr__", &scatter_event_info::repr)
		.def("__str__", &scatter_event_info::repr);
}


inline Py_Capture Py_Capture::as_type(const char* dtype) const
{
	std::string_view _dtype = dtype;
	if (_dtype == "float")
#ifdef DISTRIB_USE_FLOAT
		return as_type_t<float>();
#else
		throw pybind11::type_error("unsupported distrib type 'float'");
#endif
	else if (_dtype == "double") {
#ifdef DISTRIB_USE_DOUBLE
		return as_type_t<double>();
#else
		throw pybind11::type_error("unsupported distrib type 'double'");
#endif
	}
	else {
		throw pybind11::type_error("wrong data type: " + std::string(dtype) + ", expect float or double");
	}
}
Py_Capture& Py_Capture::add(Py_Capture const& _another)
{
	for (auto const& event : events) {
		if (event.unique) {
			for (const auto & _event2 : _another.events) {
				if (event == _event2) {
					throw pybind11::value_error("trying to combine"
						"captures with intersecting events");
				}
			}
		}
	}
	std::visit([](auto& m_d1, const auto& m_d2)->void {
		using T = decltype(m_d1.get_vtype());
		m_d1.values() += m_d2.values().template cast<T>();
	}, m_distrib, _another.m_distrib);
	return *this;
}
Py_Capture Py_Capture::operator_plus(Py_Capture const& _another) const
{
	Py_Capture Ret = copy();
	Ret.add(_another);
	return Ret;
}
Py_Capture::Py_Capture(
	Py_EL_Grid const & m_grid,
	pybind11::dict m_py_object
)
	:Py_Distribution(CreateDistribFromDict(m_grid,
		m_py_object)),
	events(
		m_py_object["events"].cast<decltype(events)>()
	) {}

Py_Capture Py_Capture::from_dict(pybind11::dict const & distr_dict) {
	pybind11::handle G = distr_dict["grid"];
	if (pybind11::isinstance<pybind11::dict>(G)) {
		return Py_Capture(
			Py_EL_Grid::from_dict1(
				G.cast<pybind11::dict>()
			), distr_dict
		);
	}
	else {
		return Py_Capture(
			G.cast<Py_EL_Grid>(), distr_dict
		);
	}
}

Py_Capture Py_Capture::copy() const
{
	return Py_Capture(static_cast<Py_Distribution const&>(*this));
}
pybind11::dict Py_Capture::get_object(pybind11::handle self)
{
	using namespace pybind11::literals;
	return pybind11::dict(
		"type"_a = "evdm.Capture",
		"values"_a = get_array(self, -1, true),
		"grid"_a = getGrid(),
		"events"_a = grob::map(
			events,
			[](scatter_event_info const& event) {
				return event.to_object();
			}),
		"padding"_a = get_padding()
	);
}

void Py_Capture::add_to_python_module(pybind11::module_& m) {
	namespace py = pybind11;
	py::class_<Py_Capture, Py_Distribution>(m, "Capture")
		.def(py::init([](
				Py_EL_Grid const& mGridEL,
				const char* dtype,
				pybind11::handle Init) ->Py_Capture
			{ return CreatePyDistrib(mGridEL, dtype, Init); }),
			"constructor of Capture(Distrib) class\n\n"
			"Parameters:\n"
			"___________\n"
			"ELGrid : GridEL\n\tgrid, where distribution is created.\n"
			"dtype : string\n\t'float' or 'double'.\n"
			"init : lambda\n\tinitialiser fuinction (density)"
			"of 3 args: (ptype,e,l), "
			"default - None, meaning zero dirtribution.",
			py::arg("ELGrid"),
			py::arg_v("dtype", "float"),
			py::arg_v("init", py::none()))
		.def(py::init(&CreateDistribFromDict),
			"constructor of Capture(Distrib) class\n\n"
			"Parameters:\n"
			"___________\n"
			"ELGrid : GridEL\n\tgrid, where distribution is created.\n"
			"values : array\n\tnumpy array of values.\n",
			py::arg("ELGrid"),
			py::arg("values")
		).def(py::init(&Py_Capture::from_dict))
		.def("__getstate__", [](py::handle self)->py::dict {
			return py::cast<Py_Capture&>(self).get_object(self);
			}
		)
		.def("__setstate__", [](py::handle self, py::dict state) {
				self.attr("__init__")(state);
			}
		)
		.def(py::init(&Py_Distribution::CreatePyDistribFromArray))
		.def("to_object", [](py::handle self)->py::dict {
			return py::cast<Py_Capture&>(self).get_object(self);
		})
		.def("as_type", &Py_Distribution::as_type, 
			"creating Capture with another dtype",
			py::arg("dtype")
		)
		.def("copy",&Py_Capture::copy)
		.def("append",&Py_Capture::add,
			"adds to capture extra events capture", 
			py::arg("second_capture"))
		.def("__add__", &Py_Capture::operator_plus)
		.def_property_readonly("events", 
			[](const Py_Capture& m_c) {return m_c.get_events(); });
}


