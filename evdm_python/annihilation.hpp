#pragma once
#include "core_python.hpp"

struct Py_Pre_Ann {
	PreAnn_Variant_t m_preann;
	Py_Pre_Ann(Py_EL_Grid const& mGridEL, size_t Nmk_bin, std::string_view dtype, pybind11::handle update_function);
	void add_ann(Py_Matrix & AnnMatrix,size_t ptype0, size_t ptype1, double a0,double av) const;

	static void add_to_python_module(pybind11::module& m);
};