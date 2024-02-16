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

template <typename T>
auto make_py_array_slice(pybind11::array_t<T>& mvector) {
	return grob::make_slice(mvector.mutable_data(),0,mvector.size());
}
template <typename T>
auto make_py_array_slice(pybind11::array_t<T>const& mvector) {
	return grob::make_slice(mvector.data(), 0, mvector.size());
}


#endif//GROB_PYTHON_HPP