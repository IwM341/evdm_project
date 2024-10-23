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

template <typename Vect_t>
auto make_py_array_view(Vect_t const& arr,pybind11::handle owner) {
	namespace py = pybind11;
	typedef typename std::decay_t<decltype(arr[0])> T;
	T const* data = arr.data();
	size_t size = arr.size();
	return py::array_t<T>(size,data,owner);
}



template <typename T>
auto make_py_array_slice(pybind11::array_t<T>const& mvector) {
	return grob::make_slice(mvector.data(), 0, mvector.size());
}

template <typename T>
inline size_t tuple_check_array(
	pybind11::array const& X,std::type_identity<T>
) {
	if (X.dtype().equal(pybind11::dtype::of<T>())) {
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

template <typename...Ts>
std::variant<pybind11::array_t<Ts>...> array_variant(
	pybind11::array const& X
) {
	typedef std::variant<pybind11::array_t<Ts>...> Ret_t;
	size_t i = tuple_check_array(X, std::type_identity<Ts>{}...);
	if (i < sizeof...(Ts)) {

	}
	else {
		throw pybind11::type_error("got array of unsupported type");
	}
	return evdm::make_variant<Ret_t>(i, [&X](auto type_id) {
		return typename decltype(type_id)::type(X);
	});
}

template <typename T>
bool py_type_check(pybind11::handle h, std::type_identity<std::tuple<T>>) {
	return pybind11::type(h).equal(pybind11::dtype::of<T>());
}
template <typename T, typename...Ts>
bool py_type_check(pybind11::handle h, std::type_identity<std::tuple<T, Ts...>>) {
	return py_type_check(h, std::type_identity<std::tuple<T>>{}) ||
		py_type_check(h, std::type_identity<std::tuple<Ts...>>{});
}

template <typename T>
const char* type_name();

template <typename T1>
inline size_t tuple_check_type(
	std::string_view S,
	std::type_identity<T1> self_T1
) {
	return !(S == type_name<T1>());
}

template <typename T1, typename...Ts>
inline size_t tuple_check_type(
	std::string_view S,
	std::type_identity<T1> self_T1,
	std::type_identity<Ts>...self_Ts
) {
	if (tuple_check_type(S, self_T1) == 0) {
		return 0;
	}
	else {
		return 1 + tuple_check_type(S, self_Ts...);
	}
}

template <typename...Ts>
std::variant<
	std::type_identity<Ts>...
> type_from_str(
	std::string_view S,
	std::type_identity<std::tuple<Ts...>>
) {
	size_t i = tuple_check_type(S, std::type_identity<Ts>{}...);

	if ( i >= sizeof...(Ts)) {
		throw pybind11::type_error("fail in attempt to determine type from string");
	}
	return evdm::make_variant_alt(i, std::type_identity<Ts>{}...);
}


template <typename Validator_t, typename T1>
inline size_t tuple_get_type(
	Validator_t validator,
	std::string_view S,
	std::type_identity<T1> self_T1
) {
	return !validator(self_T1,S);
}

template <
	typename Validator_t,
	typename T1, typename...Ts
> inline size_t tuple_get_type(
	Validator_t validator,
	std::string_view S,
	std::type_identity<T1> self_T1,
	std::type_identity<Ts>...self_Ts
) {
	if (tuple_get_type(validator,S, self_T1) == 0) {
		return 0;
	}
	else {
		return 1 + tuple_get_type(validator,S, self_Ts...);
	}
}


template <typename Validator_t, typename...Ts>
std::variant<
	std::type_identity<Ts>...
> VariantFromString(
	Validator_t Validator, 
	std::string_view S, 
	std::type_identity<std::tuple<Ts...>>
) {
	size_t i = tuple_get_type(S, std::type_identity<Ts>{}...);

	if (i >= sizeof...(Ts)) {
		throw std::runtime_error(
			"fail in attempt to determine type from string"
		);
	}
	return evdm::make_variant_alt(i, std::type_identity<Ts>{}...);
}

template <typename Validator_t, typename...Ts>
std::variant<
	std::type_identity<Ts>...
> VariantFromString(
	Validator_t Validator,
	std::string_view S,
	std::type_identity<std::variant<Ts...>>
) {
	size_t i = tuple_get_type(Validator,S, std::type_identity<Ts>{}...);

	if (i >= sizeof...(Ts)) {
		throw std::runtime_error(
			"fail in attempt to determine type from string"
		);
	}
	return evdm::make_variant_alt(i, std::type_identity<Ts>{}...);
}

/// @brief get value of type T from dict
template <typename T>
inline T pyget(
	T default_value,
	pybind11::dict const& D,
	const char *key,
	bool warning = false) {
	if (D.contains(key)) {
		return D[key].cast<T>();
	}
	else {
		if (warning) {
			pybind11::print(
				"'", key,
				"' is not specified, set default value: ",
				default_value
			);
		}
		return default_value;
	}
}

#endif//GROB_PYTHON_HPP