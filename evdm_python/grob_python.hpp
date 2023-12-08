#ifndef GROB_PYTHON_HPP
#define GROB_PYTHON_HPP

#include <grob/grid_objects.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "debugdef.hpp"

template <typename GridFunction_t>
auto make_python_function_1D(GridFunction_t const& F) {
	namespace py = pybind11;
	typedef typename GridFunction_t::value_type T;
	py::array_t<T> Grid(F.Grid.size());
	py::array_t<T> Values(F.Grid.size());

	T* GridData = Grid.mutable_data();
	T* ValueData = Values.mutable_data();
	for (size_t i = 0; i < F.Grid.size(); ++i) {
		GridData[i] = F.Grid[i];
		ValueData[i] = F.Values[i];
	}
	return py::make_tuple( Grid ,Values);
}


#endif//GROB_PYTHON_HPP