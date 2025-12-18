#include "matrix_python.hpp"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "grob_python.hpp"



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
		return Creator(std::type_identity<float>{});
#else
		throw pybind11::type_error("unsupported distrib type 'float'");
#endif
		
	}
	else if (_dtype == "double") {
#ifdef DISTRIB_USE_DOUBLE
		return Creator(std::type_identity<double>{});
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
			"ELGrid : GridEL\n\tEL Grid.\n"
			"dtype : string\n\tfloat or double.\n",
			py::arg("ELGrid"),
			py::arg_v("dtype", "float"))
		.def(py::init(&Py_Matrix::fromArray_impl))
		.def(py::init(&Py_Matrix::fromArray2))
		.def(py::init(&Py_Matrix::from_object))
		.def(py::init(&Py_Matrix::from_dict))
		.def(py::pickle([](py::handle self)->py::dict {
				return py::cast<Py_Matrix&>(self).get_object(self);
			}, &Py_Matrix::from_dict)
		)
		.def("__repr__", &Py_Matrix::repr)
		.def("__str__", &Py_Matrix::repr)
		.def("copy",&Py_Matrix::copy)
		.def("as_type", &Py_Matrix::as_type,
			"creating Matrix with another dtype",
			py::arg("dtype"))
		.def("append",&Py_Matrix::add,"adds to capture extra events capture")
		.def("__add__",&Py_Matrix::operator_plus)
		.def("calc_diag", &Py_Matrix::calc_diag,
			"make diag values -sum of scatter probabilities",
			py::arg_v("count_evap",false))
		.def("",&Py_Matrix::evap_distrib,
			"make new evolution matrix (S_{ii} = -\\sum_j{S_{ji}})")
		.def("to_numpy", [](
			py::handle self, int ptype_in, int ptype_out,bool raw)->py::object {
				return py::cast<Py_Matrix&>(self).get_matrix(
					self, ptype_in, ptype_out, raw
				);
			},
			"gives numpy array view to scatter matrix.\n\n"
			"Parameters:\n"
			"___________\n"
			"ptype_in : int\n\tin wimp type, default -1, meaning all.\n"
			"ptype_out : int\n\tout wimp type, default -1, meaning all.\n"
			/*"is_raw : bool\n\tgives raw array with padding."*/,
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
			"gives distribution of diag scatter matrix.\n\n"
				"Parameters:\n"
				"___________\n"
				"ptype : int\n\twimp type, default -1, meaning all.")
		.def_property_readonly("evap_histo", &Py_Matrix::evap_distrib,
			"gives evaporation vector of matrix\n")
		.def_property_readonly("grid", &Py_Matrix::getGrid)
		.def_property_readonly("events", &Py_Matrix::get_events);
}


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
		) {
		return Py_Matrix(mat_dstrb.first, mat_dstrb.second);
	}, m_matrix);
}

void Py_Matrix::calc_diag(bool CountEvap) {
	std::visit([CountEvap]<_MATRIX_TMPL_>(Matrix_Pair_Inst<_MATRIX_PARS_> &matr)->void {
		if(CountEvap)
			matr.first.to_scatter(matr.second.values());
		else
			matr.first.to_scatter();
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
	if constexpr (std::is_same_v<T1, T2>)
	{
		matr.first.values() += matr2.first.values();
		matr.second.values() += matr2.second.values();
	}
	else {
		matr.first.values() += matr2.first.values().template cast<T1>();
		matr.second.values() += matr2.second.values().template cast<T1>();
	}

		}, m_matrix, _another.m_matrix);
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
		}, m_matrix);
}

Py_Distribution Py_Matrix::evap_distrib()
{
	return std::visit(
		[]<_MATRIX_TMPL_>(
			Matrix_Pair_Inst<_MATRIX_PARS_> &m_p
			)->Py_Distribution {
		return m_p.second;
	}, m_matrix
	);
}



pybind11::object Py_Matrix::get_matrix(
	pybind11::handle self, int p_in, int p_out, bool raw
)
{
	return std::visit(
		[p_in, p_out, self, raw]<_MATRIX_TMPL_>(
			Matrix_Pair_Inst<_MATRIX_PARS_> &mat
			)->pybind11::object {
		evdm::GridMatrix<_MATRIX_PARS_>& m_mat = mat.first;
		if (!raw) {
			return pybind11::cast(m_mat.block(p_in, p_out));
		}
		else {
			return pybind11::cast(m_mat.block(p_in, p_out));
		}

	},m_matrix);
}

pybind11::dict Py_Matrix::get_object(pybind11::handle self)
{
	using namespace pybind11::literals;
	return pybind11::dict(
		"type"_a = "evdm.Matrix",
		"matrix"_a = get_matrix(self, -1, -1, true),
		"evap"_a = get_evap(self, -1, true),
		"grid"_a = getGrid(),
		/*"padding"_a = get_padding(),*/
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
	pybind11::object mat_arr = m_dict["matrix"];
	pybind11::array evap_arr = m_dict["evap"];
	/*size_t padding = m_dict["padding"].cast<size_t>();*/
	auto _events = m_dict["events"];
	std::vector<pybind11::handle>
		m_events(_events.begin(), _events.end());

	auto events = grob::map(
		m_events,
		[](pybind11::handle D_ev) {
			return scatter_event_info(D_ev.cast< pybind11::dict>());
		}
	);
	Py_Matrix Matr = Py_Matrix::fromArray_impl(grid, mat_arr, evap_arr/*, padding*/);
	Matr.events = events;
	return Matr;
}

Py_Matrix Py_Matrix::from_dict(pybind11::dict const& distr_dict) {
	pybind11::handle G = distr_dict["grid"];
	if (pybind11::isinstance<pybind11::dict>(G)) {
		return Py_Matrix::from_object(
			Py_EL_Grid::from_dict1(
				G.cast<pybind11::dict>()
			), distr_dict
		);
	}
	else {
		return Py_Matrix::from_object(
			G.cast<Py_EL_Grid>(), distr_dict
		);
	}
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
	}, m_matrix
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
			)-> Py_Distribution {
		evdm::GridMatrix<_MATRIX_PARS_>& m_mat = mat.first;
		auto arr = m_mat.values().diagonal();
		return Py_Distribution(m_mat.Grid, arr/*, m_mat.get_padding()*/);

	}, m_matrix
	);
}

Py_Matrix Py_Matrix::fromArray_impl(
	Py_EL_Grid const& mGridEL,
	pybind11::handle Mat,
	pybind11::array const& Evap/*,
	int padding*/)
{
	auto array_var = array_variant<DISTRIB_TYPE_LIST>(Evap);

	return std::visit([&]<class T>(pybind11::array_t<T> const& mEvap) {
		size_t gr_size = mGridEL.size() * mGridEL.ptypes();

		return Py_Matrix(
			mGridEL, 
			Mat.cast<evdm::SpMatrix_t<T>>(),
			mEvap.template cast<Eigen::VectorX<T>>() );
	}, array_var);
}

Py_Matrix Py_Matrix::fromArray1(
	Py_EL_Grid const& mGridEL,
	pybind11::array const& Mat, 
	pybind11::array const& Evap)
{
	return Py_Matrix::fromArray_impl(mGridEL, Mat, Evap);
}

Py_Matrix Py_Matrix::fromArray2(
	Py_EL_Grid const& mGridEL,
	pybind11::handle Mat
) 
{
	static const pybind11::object SpM_float = 
		pybind11::cast(evdm::SpMatrix_t<float>(1, 1));
	static const pybind11::object SpM_double =
		pybind11::cast(evdm::SpMatrix_t<double>(1, 1));

	auto array_var = getType(
		[Mat](auto X)->bool {
			typedef decltype(X) T;
			if (pybind11::isinstance<pybind11::array_t<T>>(Mat)) {
				return true;
			}
			if (pybind11::isinstance(Mat, pybind11::type::of(SpM_float))) {
				auto mdtype = Mat.attr("dtype");
				if constexpr (std::is_same_v<T, float>) {
					return mdtype.is(SpM_float.attr("dtype"));
				}
				else if constexpr (std::is_same_v<T, double>) {
					return mdtype.is(SpM_double.attr("dtype"));
				}
			}
			if (pybind11::isinstance<pybind11::list>(Mat)) {
				return std::is_same_v<T, float>;
			}
			return false;
		}, std::type_identity<std::variant<DISTRIB_TYPE_LIST>>{},
		std::type_identity<float>{}
	);

	return std::visit([&]<class T>(std::type_identity<T>) {
		Eigen::VectorX<T> Evap = {};
		return Py_Matrix(
			mGridEL, pybind11::cast<evdm::SpMatrix_t<T>>(Mat), Evap
		);
	}, array_var);
}


/*
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
			size_t evap_stride = evap_.strides()[0]/sizeof(T);
			/
			auto Ret = Py_Matrix(
				evdm::make_Matrix(_grid, mat_ptr, N, padding),
				evdm::make_Distribution_data(
					_grid,evap_ptr, evap_stride, evap_N,padding
				)
			);
			Py_Matrix Ret(_grid, mat_ptr,_grid.Grid->size(),)
			Ret.events = m_events;
			return Ret;
		}, m_type,m_grid.m_grid);
}*/
/*
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
	*/