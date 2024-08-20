#ifndef GROB_PYTHON_HPP
#define GROB_PYTHON_HPP

#include <grob/grid_objects.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "debugdef.hpp"

template <typename GridFunction_t>
pybind11::tuple make_python_function_1D(GridFunction_t const& F) {
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

template <typename Vect_t>
auto make_py_array(Vect_t const & arr) {
	namespace py = pybind11;
	typedef typename std::decay_t<decltype(arr[0])> T;
	py::array_t<T> Array(arr.size());
	T* _data = Array.mutable_data();
	for (size_t i = 0; i < arr.size(); ++i) {
		_data[i] = arr[i];
	}
	return Array;
}

template <typename T>
auto make_py_array_slice(pybind11::array_t<T>const& mvector) {
	return grob::make_slice(mvector.data(), 0, mvector.size());
}

template <typename T>
inline size_t tuple_check_array(
	pybind11::array const& X,std::type_identity<T>,
) {
	if (X.dtype() == pybind11::dtype::of<T>()) {
		return 0;
	}else {
		return 1;
	}
}

template <typename T1,typename...Ts>
inline size_t tuple_check_array(
	pybind11::array const& X, 
	std::type_identity<T1> self_T1,
	std::type_identity<Ts>...self_Ts
) {
	if (tuple_check_array(X, self_T1) == 0){
		return 0;
	}
	else {
		return 1 + tuple_check_array(X, self_Ts...);
	}
}

template <typename...Ts>
std::variant<std::type_identity<Ts>...> array_type(pybind11::array const &X) {
	size_t i = tuple_check_array(X, std::type_identity<Ts>{}...);
	if (i < sizeof...(Ts)) {
	
	}
	else {
		throw pybind11::type_error("got array of unsupported type");
	}
	return evdm::make_variant_alt(i, std::type_identity<Ts>{}...);
}

#endif//GROB_PYTHON_HPP