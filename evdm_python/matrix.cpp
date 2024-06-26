#include "core_python.hpp"
#include <pybind11/stl.h>

Py_EL_Grid Py_Matrix::getGrid() const
{
	return std::visit([](const auto& distrib) {
		return Py_EL_Grid(distrib.first.Grid);
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
		type_name<decltype(_distrib.first.get_vtype())>() + ")";
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

void Py_Matrix::calc_diag() {
	std::visit([]<_MATRIX_TMPL_>(Matrix_Pair_Inst<_MATRIX_PARS_> &matr)->void {
		matr.first.to_scatter(matr.second.values());
	}, m_matrix);
}

Py_Distribution Py_Matrix::evap_distrib()
{
	return std::visit(
		[]<_MATRIX_TMPL_>(Matrix_Pair_Inst<_MATRIX_PARS_> &m_p)->Py_Distribution {
			return m_p.second;
		},m_matrix
	);
}

template <typename BlockType>
pybind11::array make_array_from_eigen(BlockType && m_block, pybind11::handle self) {
	using _T1 = std::remove_pointer_t<decltype(m_block.data())>;
	constexpr size_t _v_s_ = sizeof(_T1);
	auto shape = pybind11::array::ShapeContainer(
		{ m_block.rows(), m_block.cols() }
	);
	auto strides = pybind11::array::ShapeContainer(
		{ _v_s_ * m_block.rowStride(),_v_s_ * m_block.colStride() }
	);
	return pybind11::array_t<_T1>(
		pybind11::buffer_info(m_block.data(), shape, strides), self
	);
}

pybind11::array Py_Matrix::get_matrix(
	pybind11::handle self,int p_in , int p_out,bool raw
)
{
	return std::visit(
		[p_in, p_out, self,raw]<_MATRIX_TMPL_>(
			Matrix_Pair_Inst<_MATRIX_PARS_> & mat
		)->pybind11::array {
			evdm::GridMatrix<_MATRIX_PARS_>& m_mat = mat.first;
			if (!raw) {
				return make_array_from_eigen(m_mat.block(p_in, p_out),self);
			}
			else {
				return make_array_from_eigen(m_mat.raw_matrix(), self);
			}
			
		},
		m_matrix
	);
}

pybind11::array Py_Matrix::get_diag(pybind11::handle self, int p_in, int p_out, bool raw) {
	return std::visit(
		[p_in, p_out, self, raw]<_MATRIX_TMPL_>(
			Matrix_Pair_Inst<_MATRIX_PARS_> &mat
		)->pybind11::array {
		evdm::GridMatrix<_MATRIX_PARS_>& m_mat = mat.first;
		constexpr size_t _v_s_ = sizeof(_T1);
		if (!raw) {
			return make_array_from_eigen(m_mat.diag(p_in, p_out), self);
		}
		else {
			return make_array_from_eigen(m_mat.raw_diag(), self);
		}
	},m_matrix
	);
}
pybind11::array Py_Matrix::get_evap(pybind11::handle self, int p_type, bool raw) {
	return std::visit(
		[p_type, self, raw]<_MATRIX_TMPL_>(
			Matrix_Pair_Inst<_MATRIX_PARS_> &mat
		)->pybind11::array {
		evdm::Distribution<_MATRIX_PARS_>& m_ev = mat.second;
		constexpr size_t _v_s_ = sizeof(_T1);
		if (!raw) {
			return make_array_from_eigen(m_ev.block(p_type), self);
		}
		else {
			return make_array_from_eigen(m_ev.raw_vector(), self);
		}
	}, m_matrix
	);
}

Py_Distribution Py_Matrix::get_diag_distrib() {
	return std::visit(
		[]<_MATRIX_TMPL_>(
			Matrix_Pair_Inst<_MATRIX_PARS_> &mat
		)-> Py_Distribution{
		evdm::GridMatrix<_MATRIX_PARS_>& m_mat = mat.first;
		auto arr = m_mat.values().diagonal();
		return Py_Distribution(m_mat.Grid, arr, m_mat.get_padding());
		
	}, m_matrix
	);
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
	py::class_<Py_Matrix>(m, "Matrix", 
		"class containing scattering matrix and evaporation vector")
		.def(py::init(&CreatePyMatrix),
			"constructor of Matrix class\n\n"
			"Parameters:\n"
			"___________\n"
			"ELGrid : GridEL\n\tCreated EL Grid.\n"
			"dtype : string\n\tfloat or double.\n",
			py::arg("ELGrid"),
			py::arg_v("dtype", "float"))
		.def("__repr__", &Py_Matrix::repr)
		.def("__str__", &Py_Matrix::repr)
		.def("as_type", &Py_Matrix::as_type,
			"creating Matrix with another dtype",
			py::arg("dtype"))
		.def("calc_diag", &Py_Matrix::calc_diag,
			"make diag values -summ of scatter probabilities")
		.def("to_numpy", [](
			py::handle self, int ptype_in, int ptype_out,bool raw)->py::array {
				return py::cast<Py_Matrix&>(self).get_matrix(
					self, ptype_in, ptype_out, raw
				);
			},
			"gives numpy array view to scatter matrix.\n\n"
			"Parameters:\n"
			"___________\n"
			"ptype_in : int\n\tin wimp type, default -1, meaning all.\n"
			"ptype_out : int\n\tout wimp type, default -1, meaning all.\n"
			"is_raw : bool\n\tgives raw array with padding.",
			py::arg_v("ptype_in",-1), py::arg_v("ptype_out", -1),
				py::arg_v("is_raw", false))
		.def("to_numpy_diag", [](
			py::handle self, int ptype, bool raw)->py::array {
				return py::cast<Py_Matrix&>(self).get_diag(
					self, ptype, ptype, raw
				);
			},
			"gives numpy array view to diag scatter matrix.\n\n"
			"Parameters:\n"
			"___________\n"
			"ptype : int\n\twimp type, default -1, meaning all.\n"
			"is_raw : bool\n\tgives raw array with padding.",
			py::arg_v("ptype", -1), 
			py::arg_v("is_raw", false))
		.def("diag_distrib", &Py_Matrix::get_diag_distrib,
			"givesdistribution of diag scatter matrix.\n\n"
				"Parameters:\n"
				"___________\n"
				"ptype : int\n\twimp type, default -1, meaning all.")
		.def_property_readonly("evap_histo", &Py_Matrix::evap_distrib,
			"gives evaporation vector of matrix\n")
		.def_property_readonly("grid", &Py_Matrix::getGrid)
		.def_property_readonly("events", &Py_Matrix::get_events);
}