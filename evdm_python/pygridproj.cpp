#include "pygridproj.hpp"
#include "grid_python.hpp"
#include "grob_python.hpp"
#include <evdm/core/core_project.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <variant>

pybind11::object GridProjectionMatrix(
	Py_EL_Grid OutGrid, Py_EL_Grid inGrid, double p , double q 
) 
{
	return std::visit(
		[p, q](
			auto const & Gout, 
			auto const& Gin
		){
			typedef decltype(Gout.get_grid_vtype()) T;
			
			T Emin = Gout.Grid->inner(0).grid()[0].left;
			evdm::GridProjector<T> mproj(p,q, Emin);
			evdm::SpMatrix_t<T> spmat = mproj(*Gout.Grid, *Gin.Grid);
			return pybind11::cast(spmat);
		},
		OutGrid.m_grid, inGrid.m_grid
	);
}

void add_to_python_module_GridProjectionMatrix(pybind11::module_& m) {
	namespace py = pybind11;
	m.def("GridProjectionMatrix", GridProjectionMatrix,
		"Make projection matrix from Grid_in to Grid_out distribution",
		py::arg("Grid_out"), py::arg("Grid_in"),
		py::arg_v("p",1,"measure of bin proportional dE^p"),
		py::arg_v("q", 2, "measure of bin proportional dl^q")
	);
}