#include "core_python.hpp"

Py_EL_Grid Py_Matrix::getGrid() const
{
	return std::visit([](const auto& distrib) {
		return Py_EL_Grid(distrib.Grid);
	}, m_matrix);
}

std::vector<scatter_event_info> const& Py_Matrix::get_events() const
{
	return events;
}

std::string Py_Matrix::repr() const
{
	return std::visit([](auto const& _distrib) {
		return std::string("Matrix( dtype = ") +
		type_name<decltype(_distrib.get_vtype())>() + ")";
	}, m_matrix);
}

Py_Matrix Py_Matrix::as_type(const char* type_n) const
{
	std::string_view _dtype = type_n;
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
		throw pybind11::type_error("wrong data type: " + std::string(type_n) + ", expect float or double");
	}
}

Py_Distribution Py_Matrix::total_probs() const
{
	throw 1;
}


Py_Matrix CreatePyMatrix(
	Py_EL_Grid const& mGridEL,
	const char* dtype
)
{
	auto Creator = [&](auto _type_marker) {
		return Py_Matrix(mGridEL, _type_marker);
	};

	std::string_view _dtype = dtype;
	if (_dtype == "float") {
#ifdef DISTRIB_USE_FLOAT
		return Creator(Py_EL_Grid::m_type_marker<float>{});
#else
		throw pybind11::type_error("unsupported distrib type 'float'");
#endif
		
	}
	else if (_dtype == "double") {
#ifdef DISTRIB_USE_DOUBLE
		return Creator(Py_EL_Grid::m_type_marker<double>{});
#else
		throw pybind11::type_error("unsupported distrib type 'double'");
#endif
	}
	else {
		throw pybind11::type_error("wrong data type: " + std::string(_dtype) + ", expect float or double");
	}
}

void Py_Matrix::add_to_python_module(pybind11::module_& m) {
	namespace py = pybind11;
	py::class_<Py_Matrix>(m, "Matrix")
		.def(py::init(&CreatePyMatrix),
			"constructor of Matrix class\n"
			"ELGrid -- Created EL Grid\n"
			"dtype -- float or double\n",
			py::arg("ELGrid"),
			py::arg_v("dtype", "float"))
		.def("__repr__", &Py_Matrix::repr)
		.def("__str__", &Py_Matrix::repr)
		.def("as_type", &Py_Matrix::as_type,
			"creating Matrix with another dtype", 
			py::arg("dtype")
		)
		.def_property_readonly("grid", &Py_Matrix::getGrid)
		.def_property_readonly("events", &Py_Matrix::get_events);
}