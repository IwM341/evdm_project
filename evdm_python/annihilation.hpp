#pragma once
#include "matrix_python.hpp"


struct Py_Pre_Ann {
	PreAnn_Variant_t m_preann;
	
	template <typename T, typename B_vt,
		typename G_vt, evdm::GridEL_type grid_type>
	Py_Pre_Ann(
		evdm::GridAnnPreMatrix<T, B_vt, G_vt, grid_type> m_ann
	) : m_preann(std::move(m_ann)){}


	Py_Pre_Ann(
		Py_EL_Grid const& mGridEL,
		size_t Nmk_bin, std::string_view dtype,
		double Rmin, double Rmax, size_t seed,
		pybind11::handle update_function);
	void add_ann(
		Py_Matrix & AnnMatrix,
		size_t ptype0, size_t ptype1, 
		double a0,double av
	) const;

	static Py_Pre_Ann from_object(Py_EL_Grid const & Grid,pybind11::dict const&);
	static Py_Pre_Ann from_dict(pybind11::dict const&);
	pybind11::dict get_object(pybind11::handle self);
	std::string repr() const;
	Py_Pre_Ann copy() const;

	template <typename T>
	Py_Pre_Ann as_type_t() const {
		return std::visit([this](auto const& m_ann)->Py_Pre_Ann {
			return m_ann.template as_type<T>();
		},m_preann);
	}

	Py_Pre_Ann as_type(const char * type_name) const;

	pybind11::array A0(pybind11::handle self);
	pybind11::array Av(pybind11::handle self);
	Py_EL_Grid getGrid() const;

	static void add_to_python_module(pybind11::module& m);
};