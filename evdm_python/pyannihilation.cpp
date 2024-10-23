#include "annihilation.hpp"
#include "progress_log.hpp"
#include <evdm/utils/variant_tools.hpp>
Py_Pre_Ann::Py_Pre_Ann(Py_EL_Grid const& mGridEL, size_t Nmk_bin, std::string_view dtype, pybind11::handle update_function):
	m_preann(std::visit([Nmk_bin, dtype,&update_function]<_GRID_EL_TMPL_,typename T>(const  evdm::EL_Grid<_GRID_EL_PARS_>& m_gr,evdm::self_type<T>)->PreAnn_Variant_t
	{
		evdm::progress_omp_function<> m_prog = make_progress_func(update_function);
		if constexpr (std::is_same_v<T, void>) {
			throw pybind11::type_error("unrecongnesed dtype");
		}
		else {
			if constexpr (evdm::variant_consist_of_v<evdm::GridAnnPreMatrix<T, _T1, _T2, _m_grid_t>, PreAnn_Variant_t>) {
				return evdm::GridAnnPreMatrix<T, _T1, _T2, _m_grid_t>(m_gr, Nmk_bin, evdm::xorshift<_T2>{}, m_prog);
			}
			else {
				throw pybind11::type_error("unsupported dtype");
			}
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

Py_EL_Grid Py_Pre_Ann::getGrid() const
{
	return std::visit([](const auto& item) {
		return Py_EL_Grid(item.Grid);
	}, m_preann);
}

void Py_Pre_Ann::add_to_python_module(pybind11::module& m)
{
	namespace py = pybind11;
	py::class_<Py_Pre_Ann>(m, "AnnPre","python class for pre annihilation matrix")
		.def(py::init<Py_EL_Grid const& ,size_t,std::string_view,py::handle>(),
			"constructor.\n\n"
			"Parameters:\n\t"
			"grid : el grid.\n\t"
			"Nmk : number of MK integratons per bin.\n\t"
			"dtype : type of values.\n\t"
			"bar : update progress bar function.\n",
			py::arg("grid"), py::arg("Nmk"),
			py::arg_v("dtype","float"), py::arg_v("bar", py::none()))
		.def("add_to_matrix",&Py_Pre_Ann::add_ann,
			"build part of annihilation matrix of ptype0, ptype1, the matrix would by Aij += a0*a_0_ij+av*a_v_ij"
			"where annihilation determines by dN_i/dt = (n_{\\infty} \\sigma_{ann} v_esc) (A_{ij}  N_j) N_i"
			"Parameters:\n\t"
			"matrix : input matrix to build (modify)\n\t"
			"ptype0, ptype1: scattering types\n\t"
			"a0 : coefficient of \\sigma_{ann} v_esc = const annihilation\n\t"
			"av : coefficient of \\sigma_{ann} v_esc = const * v^2 annihilation\n",
			py::arg("matrix"), py::arg("ptype0"), py::arg("ptype1"), py::arg("a0"), py::arg("av"));
}
