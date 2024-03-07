#include <iostream>
#include "core_python.hpp"

template <>
const char* type_name<float>() {
	return "float";
}
template <>
const char* type_name<double>() {
	return "double";
}
template <>
const char* type_name<int>() {
	return "int";
}
template <>
const char* type_name<size_t>() {
	return "size_t";
}

