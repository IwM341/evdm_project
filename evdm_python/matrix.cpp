#include "core_python.hpp"
#include <pybind11/stl.h>
#include "grob_python.hpp"

Py_Matrix make_Py_Matrix_from_arrays(
	Py_EL_Grid const& m_grid,
	pybind11::array const& mat,
	pybind11::array const& evap,
	size_t padding,std::vector<scatter_event_info> m_events) 
{
	size_t N = mat.shape()[0];
	if (mat.ndim() != 2 ||
		!(mat.shape()[1] == N)
		) {
		throw pybind11::value_error(
			"input matrix have incorrect shape,"
			"expect N*N"
		);
	}
	typedef std::variant <
#ifdef DISTRIB_USE_FLOAT
		std::type_identity<float>
#endif
#ifdef DISTRIB_USE_DOUBLE
#ifdef DISTRIB_USE_FLOAT
		, 
#endif
		std::type_identity<double>
#endif
	> Float_variants_t;
	
	Float_variants_t m_type;
	bool type_setted = false;
#ifdef DISTRIB_USE_FLOAT
	if (check_np_type<float>(mat)) {
		m_type = std::type_identity<float>{};
		type_setted = true;
	}
#endif
#ifdef DISTRIB_USE_DOUBLE
	if (check_np_type<double>(mat)) {
		m_type = std::type_identity<double>{};
		type_setted = true;
	}
#endif
	if(!type_setted){
		throw pybind11::type_error("not supported type of array");
	}
	
		return std::visit([&]<typename T, _GRID_EL_TMPL_>(
					std::type_identity<T> m_type,
					evdm::EL_Grid<_GRID_EL_PARS_> const & _grid
				) ->Py_Matrix
		{
			pybind11::array_t<T> mat_ = mat;
			pybind11::array_t<T> evap_ = evap;
			const T* mat_ptr = mat_.data();
			const T* evap_ptr = evap_.data();
			size_t evap_N = evap_.size();
			size_t evap_stride = evap_.strides()[0];
			auto Ret = Py_Matrix(
				evdm::make_Matrix(_grid, mat_ptr, N, padding),
				evdm::make_Distribution_data(
					_grid,evap_ptr, evap_stride, evap_N,padding
				)
			);
			Ret.events = m_events;
			return Ret;
		}, m_type,m_grid.m_grid);
}
Py_Matrix::Py_Matrix(
	Py_EL_Grid const& mGridEL, 
	pybind11::dict const& saved_obj):
	Py_Matrix(make_Py_Matrix_from_arrays(
		mGridEL, 
		saved_obj["matrix"].cast<pybind11::array>(),
		saved_obj["evap"].cast<pybind11::array>(),
		saved_obj["padding"].cast<size_t>(),
		saved_obj["events"].cast<decltype(events)>()
	)){}

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

Py_Matrix Py_Matrix::copy() const
{
	return std::visit([]<_MATRIX_TMPL_>(
		Matrix_Pair_Inst<_MATRIX_PARS_> const& mat_dstrb
	){
		return Py_Matrix(mat_dstrb.first, mat_dstrb.second);
	},m_matrix);
}

void Py_Matrix::calc_diag() {
	std::visit([]<_MATRIX_TMPL_>(Matrix_Pair_Inst<_MATRIX_PARS_> &matr)->void {
		matr.first.to_scatter(matr.second.values());
	}, m_matrix);
}

Py_Matrix Py_Matrix::evolve_matrix() const
{
	Py_Matrix m_copy = *this;
	m_copy.calc_diag();
	return m_copy;
}



void Py_Matrix::combine(Py_Matrix const& _another)
{
	for (auto const& event : events) {
		if (event.unique) {
			for (const auto& _event2 : _another.events) {
				if (event == _event2) {
					throw pybind11::value_error("trying to combine"
						"matrixes with intersecting events");
				}
			}
		}
	}
	std::visit([](auto& matr, auto const& matr2)->void {
		typedef decltype(matr.first.get_vtype()) T1;
		typedef decltype(matr2.first.get_vtype()) T2;
		if constexpr (std::is_same_v<T1,T2>) 
		{
			matr.first.values() += matr2.first.values();
			matr.second.values() += matr2.second.values();
		}
		else {
			matr.first.values() += matr2.first.values().template cast<T1>();
			matr.second.values() += matr2.second.values().template cast<T1>();
		}
		
	}, m_matrix,_another.m_matrix);
}

Py_Matrix& Py_Matrix::add(Py_Matrix const& _another)
{
	combine(_another);
	return *this;
}

Py_Matrix Py_Matrix::operator_plus(Py_Matrix const& _another) const
{
	Py_Matrix m_copy = *this;
	return m_copy.add(_another);
}

size_t Py_Matrix::get_padding() const
{
	return std::visit([](const auto& m_mat_pair) {
		return m_mat_pair.second.get_padding();
	},m_matrix);
}

Py_Distribution Py_Matrix::evap_distrib()
{
	return std::visit(
		[]<_MATRIX_TMPL_>(
			Matrix_Pair_Inst<_MATRIX_PARS_> &m_p
		)->Py_Distribution {
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

pybind11::dict Py_Matrix::get_object(pybind11::handle self)
{
	using namespace pybind11::literals;
	return pybind11::dict(
		"type"_a = "evdm.Matrix",
		"matrix"_a = get_matrix(self, -1, -1, true),
		"evap"_a = get_evap(self, -1, true),
		"padding"_a = get_padding(),
		"events"_a = grob::map(
			events,
			[](scatter_event_info const& event) {
				return event.to_object();
			})
	);
}



Py_Matrix Py_Matrix::from_object(Py_EL_Grid const& grid, pybind11::dict const& m_dict)
{
	if (m_dict["type"].cast<std::string_view>() != "evdm.Matrix") {
		throw pybind11::value_error("expect evdm.Matrix type");
	}
	pybind11::array mat_arr = m_dict["matrix"];
	pybind11::array evap_arr = m_dict["evap"];
	size_t padding = m_dict["padding"].cast<size_t>();
	auto _events = m_dict["events"];
	std::vector<pybind11::handle>
		m_events(_events.begin(), _events.end());

	auto events = grob::map(
		m_events,
		[](pybind11::handle D_ev) {
			return scatter_event_info(D_ev.cast< pybind11::dict>());
		}
	);
	Py_Matrix Matr = Py_Matrix::fromArray_impl(grid, mat_arr, evap_arr, padding);
	Matr.events = events;
	return Matr;
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
		.def(py::init(&Py_Matrix::fromArray1))
		.def(py::init(&Py_Matrix::fromArray2))
		.def(py::init(&Py_Matrix::from_object))
		.def("__repr__", &Py_Matrix::repr)
		.def("__str__", &Py_Matrix::repr)
		.def("copy",&Py_Matrix::copy)
		.def("as_type", &Py_Matrix::as_type,
			"creating Matrix with another dtype",
			py::arg("dtype"))
		.def("append",&Py_Matrix::add,"adds to capture extra events capture")
		.def("__add__",&Py_Matrix::operator_plus)
		.def("calc_diag", &Py_Matrix::calc_diag,
			"make diag values -summ of scatter probabilities")
		.def("",&Py_Matrix::evap_distrib,
			"make new evolution matrix (S_{ii} = -\\sum_j{S_{ji}})")
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
			py::arg_v("is_raw", false)
		)
		.def("to_object", [](py::handle self) {
				return py::cast<Py_Matrix&>(self).get_object(self);
		},"serialization into object dict")
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

Py_Matrix Py_Matrix::fromArray_impl(
	Py_EL_Grid const& mGridEL,
	pybind11::array const& Mat,
	pybind11::array const& Evap,
	int padding)
{
	auto array_var = array_variant<DISTRIB_TYPE_LIST>(Mat);

	return std::visit([&]<class T>(pybind11::array_t<T> const& mMat) {
		pybind11::array_t<T> const& mEvap = Evap;
		size_t gr_size = mGridEL.size() * mGridEL.ptypes();
		if (mMat.ndim() != 2) {
			throw pybind11::index_error("evdm.Matrix: ndim of matrix != 2");
		}
		if (mMat.shape()[0] < gr_size || mMat.shape()[1] < gr_size) {
			throw pybind11::index_error(
				"evdm.Matrix: matrix shouldn't have less size than grid"
			);
		}
		if (mMat.shape()[0] != mMat.shape()[1]) {
			throw pybind11::index_error(
				"expect square matrix"
			);
		}
		if (mEvap.shape()[0] < gr_size) {
			throw pybind11::index_error(
				"evdm.Matrix: evap shouldn't have less size than grid"
			);
		}
		size_t m_stride_outer = mMat.strides()[0] / sizeof(T);
		size_t m_stride_inner = mMat.strides()[1] / sizeof(T);

		size_t m_stride_evap = mEvap.strides()[0] / sizeof(T);

		return Py_Matrix(mGridEL, mMat.data(), gr_size,
			m_stride_outer, m_stride_inner,
			mEvap.data(), mEvap.size(), m_stride_evap, padding);
	}, array_var);
}

Py_Matrix Py_Matrix::fromArray1(
	Py_EL_Grid const& mGridEL,
	pybind11::array const& Mat, 
	pybind11::array const& Evap)
{
	return Py_Matrix::fromArray_impl(mGridEL, Mat, Evap,-1);
}

Py_Matrix Py_Matrix::fromArray2(
	Py_EL_Grid const& mGridEL,
	pybind11::array const& Mat
) 
{
	auto array_var = array_variant<DISTRIB_TYPE_LIST>(Mat);

	return std::visit([&]<class T>(pybind11::array_t<T> const& mMat) {
		size_t gr_size = mGridEL.size() * mGridEL.ptypes();
		if (mMat.ndim() != 2) {
			throw pybind11::index_error("evdm.Matrix: ndim of matrix != 2");
		}
		if (mMat.shape()[0] < gr_size || mMat.shape()[1] < gr_size) {
			throw pybind11::index_error(
				"evdm.Matrix: matrix shouldn't have less size than grid"
			);
		}
		if (mMat.shape()[0] != mMat.shape()[1] ) {
			throw pybind11::index_error(
				"expect square matrix"
			);
		}
		
		size_t m_stride_outer = mMat.strides()[0] / sizeof(T);
		size_t m_stride_inner = mMat.strides()[1] / sizeof(T);

		return Py_Matrix(mGridEL, mMat.data(), gr_size,
			m_stride_outer, m_stride_inner,
			(T*)nullptr, 0, 1,-1);
	}, array_var);
}

