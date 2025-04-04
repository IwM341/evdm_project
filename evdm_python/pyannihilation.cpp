#include "annihilation.hpp"
#include "progress_log.hpp"
#include <evdm/utils/variant_tools.hpp>
#include "grid_python.hpp"
#include "matrix_python.hpp"
#include "grob_python.hpp"
#include <type_traits>
#include <pybind11/eigen.h>


void Py_Pre_Ann::add_to_python_module(pybind11::module& m)
{
	namespace py = pybind11;
	py::class_<Py_Pre_Ann>(m, "Ann", "python class for pre annihilation matrix")
		.def(py::init<
				Py_EL_Grid const&, size_t, 
				std::string_view, double , double , 
				size_t, py::handle
			>(),
			"constructor.\n\n"
			"Parameters:\n\t"
			"grid : el grid.\n\t"
			"Nmk : number of MK integratons per bin.\n\t"
			"dtype : type of values.\n\t"
			"bar : update progress bar function.\n",
			py::arg("grid"), py::arg("Nmk"),
			py::arg_v("dtype", "float"), 
			py::arg_v("rmin", 0),
			py::arg_v("rmax", 1),
			py::arg_v("seed", 1),
			py::arg_v("bar", py::none())
		).def("__repr__", &Py_Pre_Ann::repr)
		.def(py::init([](
			Py_EL_Grid const& G,
			py::array const& A0,
			py::array const& Av)
			{
				using namespace pybind11::literals;
				return Py_Pre_Ann::from_object(G, py::dict("A0"_a = A0, "Av"_a = Av));
			}), 
			py::arg("Grid"),py::arg("A0"), py::arg("Av")
		).def_property_readonly("A0",
			[](pybind11::handle self) {
				return self.cast<Py_Pre_Ann &>().A0(self);
			}
		).def_property_readonly("Av",
			[](pybind11::handle self) {
				return self.cast<Py_Pre_Ann &>().Av(self);
			}
		)
		.def("__str__", &Py_Pre_Ann::repr)
		.def(py::init(&Py_Pre_Ann::from_object))
		.def(py::init(&Py_Pre_Ann::from_dict))
		.def(py::pickle([](py::handle self)->py::dict {
				return py::cast<Py_Pre_Ann&>(self).get_object(self);
			}, &Py_Pre_Ann::from_dict)
		).def("copy", &Py_Pre_Ann::copy)
		.def("as_type", &Py_Pre_Ann::as_type,
			"creating Ann with another dtype",
			py::arg("dtype"))
		.def("to_object", [](py::handle self) {
				return py::cast<Py_Pre_Ann&>(self).get_object(self);
			},"serialization into object dict"
		).def_property_readonly("grid", &Py_Pre_Ann::getGrid)
		.def("add_to_matrix", &Py_Pre_Ann::add_ann,
			"build part of annihilation matrix of ptype0, ptype1, the matrix would by Aij += a0*a_0_ij+av*a_v_ij"
			"where annihilation determines by dN_i/dt = (n_{\\infty} \\sigma_{ann} v_esc) (A_{ij}  N_j) N_i"
			"Parameters:\n\t"
			"matrix : input matrix to build (modify)\n\t"
			"ptype0, ptype1: scattering types\n\t"
			"a0 : coefficient of \\sigma_{ann} v_esc = const annihilation\n\t"
			"av : coefficient of \\sigma_{ann} v_esc = const * v^2 annihilation\n",
			py::arg("matrix"), py::arg("ptype0"), py::arg("ptype1"), py::arg("a0"), py::arg("av")
		);
}


Py_Pre_Ann::Py_Pre_Ann(
	Py_EL_Grid const& mGridEL, size_t Nmk_bin,
	std::string_view dtype, double Rmin, double Rmax, size_t seed, 
	pybind11::handle update_function): 
	m_preann(
		std::visit([Nmk_bin, dtype, Rmin , Rmax, &update_function,_seed=seed]
			<_GRID_EL_TMPL_,typename T>
			(const  evdm::EL_Grid<_GRID_EL_PARS_>& m_gr,
				evdm::self_type<T>)->PreAnn_Variant_t
	{
		evdm::progress_omp_function<> m_prog = make_progress_func(update_function);
		if constexpr (std::is_same_v<T, void>) {
			throw pybind11::type_error("unrecongnesed dtype");
		}
		size_t seed = _seed;
		if (seed == 0) {
			seed = std::numeric_limits<size_t>::max();
		}
		
		if constexpr (evdm::variant_consist_of_v<
			evdm::GridAnnPreMatrix<T, _T1, _T2, _m_grid_t>, 
			PreAnn_Variant_t
		>) {
			evdm::xorshift<T> G(seed);
			pybind11::gil_scoped_release m_lock;
			return evdm::GridAnnPreMatrix<T, _T1, _T2, _m_grid_t>(
				m_gr, Rmin, Rmax, Nmk_bin,G,
				m_prog
			);
		}
		else {
			throw pybind11::type_error("unsupported dtype of annihilation");
		}
		
	}, mGridEL.m_grid, 
		evdm::make_variant_alt(dtype == "float" ? 0 : (dtype == "double" ? 1 : 2), 
			evdm::self_type<float>{},
			evdm::self_type<double>{}, 
			evdm::self_type<void>{})
		)
	)
{}

void Py_Pre_Ann::add_ann(Py_Matrix& AnnMatrix, size_t ptype0, size_t ptype1, double a0, double av) const
{
	std::visit([ptype0, ptype1, a0, av]<_ANN_MATRIX_TMPL_,_MATRIX_TMPL_>
		(evdm::GridAnnPreMatrix<_ANN_MATRIX_PARS_> const& m_preann, 
			Matrix_Pair_Inst<_MATRIX_PARS_> & m_mat) {
		if constexpr (std::is_same_v<_T2, _T2a> &&
			std::is_same_v<_T3, _T3a> &&
			std::is_same_v<_T2, _T2a> &&
			_ma_grid_t == _m_grid_t) 
		{
			m_preann.add_to(m_mat.first, ptype0, ptype1, a0, av);
		}
		else {
			throw pybind11::type_error("incompatible types of annihilation pre matrix and matrix");
		}
		
	}, m_preann, AnnMatrix.m_matrix);
}

pybind11::array Py_Pre_Ann::A0(pybind11::handle self) {
	return std::visit(
		[self]<_DISTRIB_TMPL_>(
			evdm::GridAnnPreMatrix<_DISTRIB_PARS_> & m_ann
			) {
		return make_array_from_eigen(m_ann.A0, self);
	}, m_preann);
}
pybind11::array Py_Pre_Ann::Av(pybind11::handle self) {
	return std::visit(
		[self]<_DISTRIB_TMPL_>(
			evdm::GridAnnPreMatrix<_DISTRIB_PARS_> & m_ann
			) {
		return make_array_from_eigen(m_ann.Av, self);
	}, m_preann);
}

Py_EL_Grid Py_Pre_Ann::getGrid() const
{
	return std::visit([](const auto& item) {
		return Py_EL_Grid(item.Grid);
	}, m_preann);
}

Py_Pre_Ann Py_Pre_Ann::as_type(const char* type_name) const {
	auto m_type = type_from_str(
		type_name,
		std::type_identity<std::tuple< DISTRIB_TYPE_LIST >> {}
	);
	return std::visit([this]<class U>(std::type_identity<U>) {
		return as_type_t<U>();
	}, m_type);
}



Py_Pre_Ann Py_Pre_Ann::from_object(
	Py_EL_Grid const& Grid,
	pybind11::dict const& m_dict) 
{
	pybind11::array A0 = m_dict["A0"].cast< pybind11::array>();
	pybind11::array Av = m_dict["Av"].cast< pybind11::array>();
	auto A0_var = array_variant<DISTRIB_TYPE_LIST>(A0);
	return std::visit(
		[&Av]<_DISTRIB_TMPL_>(
			evdm::EL_Grid<_T2,_T3,_m_grid_t> const& Grid,
			pybind11::array_t<_T1> const & A0
		){
		pybind11::array_t<_T1> const& _Av =
			Av.cast<pybind11::array_t<_T1>>();

		size_t N = Grid.getLE_inner_grid().size();
		if (A0.ndim() != 2 || A0.shape()[0] != N || A0.shape()[1] != N) {
			throw pybind11::index_error("Ann: A0.ndim != or A0.shape != (N,N)");
		}
		if (_Av.ndim() != 2 || _Av.shape()[0] != N || _Av.shape()[1] != N) {
			throw pybind11::index_error("Ann: A0.ndim != or A0.shape != (N,N)");
		}
		//
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> stride_t;

		size_t m_stride_outer_0 = A0.strides()[0] / sizeof(_T1);
		size_t m_stride_inner_0 = A0.strides()[1] / sizeof(_T1);
		size_t m_stride_outer_v = _Av.strides()[0] / sizeof(_T1);
		size_t m_stride_inner_v = _Av.strides()[1] / sizeof(_T1);
		Eigen::Map<
			const evdm::Matrix_t<_T1>,
			Eigen::Unaligned,
			Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
		> A0_map(A0.data(), N, N, stride_t(m_stride_inner_0, m_stride_outer_0));

		Eigen::Map<
			const evdm::Matrix_t<_T1>,
			Eigen::Unaligned,
			Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
		> Av_map(_Av.data(), N, N, stride_t(m_stride_inner_v, m_stride_outer_v));

		return Py_Pre_Ann(
			evdm::GridAnnPreMatrix<_DISTRIB_PARS_>(
				Grid,
				A0_map,
				A0_map
				)
		);
		/*
		return Py_Pre_Ann(
			evdm::GridAnnPreMatrix<_DISTRIB_PARS_>(
				Grid,
					A0.template cast<evdm::Matrix_t<_T1>>(),
					_Av.template cast<evdm::Matrix_t<_T1>>()
				)
		);*/
	}, Grid.m_grid, A0_var);
}
Py_Pre_Ann Py_Pre_Ann::from_dict(pybind11::dict const& distr_dict) {
	pybind11::handle G = distr_dict["grid"];
	if (pybind11::isinstance<pybind11::dict>(G)) {
		return Py_Pre_Ann::from_object(
			Py_EL_Grid::from_dict1(
				G.cast<pybind11::dict>()
			), distr_dict
		);
	}
	else {
		return Py_Pre_Ann::from_object(
			G.cast<Py_EL_Grid>(), distr_dict
		);
	}
}
pybind11::dict Py_Pre_Ann::get_object(pybind11::handle self){
	using namespace pybind11::literals;
	return pybind11::dict(
		"type"_a = "evdm.Ann",
		"A0"_a = A0(self),
		"Av"_a = Av(self),
		"grid"_a = getGrid()
	);
}
std::string Py_Pre_Ann::repr() const {
	return std::visit([]<_DISTRIB_TMPL_>(
		evdm::GridAnnPreMatrix<_DISTRIB_PARS_> const& m_pr) 
		{
			return std::string("Ann( dtype = ") + 
				type_name<_T1>() + ")";
		}, m_preann
	);
}
Py_Pre_Ann Py_Pre_Ann::copy() const {
	return std::visit([](auto const& m_pr) {
			return Py_Pre_Ann(m_pr);
		}, m_preann
	);
} 
