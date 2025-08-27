#pragma once
#include "grid_python.hpp"
pybind11::object GridProjectionMatrix(
	Py_EL_Grid OutGrid, Py_EL_Grid inGrid,double p = 1, double q = 2
);


void add_to_python_module_GridProjectionMatrix(pybind11::module_& m);