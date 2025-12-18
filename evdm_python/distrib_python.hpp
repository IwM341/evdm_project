#pragma once
#include "grid_python.hpp"


struct Py_Distribution {
	Distrib_Variant_t m_distrib;

	template <typename T, typename Initializer_t>
	Py_Distribution(Py_EL_Grid const& mGridEL,
		std::type_identity<T> _type,
		Initializer_t&& init) :
		m_distrib(
			std::visit([&init](auto const& _grid)->Distrib_Variant_t {
				return evdm::make_Distribution<T>(_grid, init,1);
				}, mGridEL.m_grid)
		)
	{}

	template <typename ArrayType>
	Py_Distribution(Py_EL_Grid const& mGridEL,
		ArrayType&& _array/*, size_t padding*/) :
		m_distrib(
			std::visit([&_array](auto const& _grid)->Distrib_Variant_t {
				return evdm::make_Distribution_array(
					_grid, std::forward<ArrayType>(_array), 1
				);
				}, mGridEL.m_grid)
		)
	{}


	template <typename T>
	Py_Distribution(Py_EL_Grid const& mGridEL,
	const T* _data, size_t _stride, size_t _size/*, size_t padding*/) :
	m_distrib(
		std::visit([&](auto const& _grid)->Distrib_Variant_t {
			return evdm::make_Distribution_data<T>(_grid, _data, _stride, _size, 1/*padding*/);
			}, mGridEL.m_grid)
	)
	{}

	template <typename T, typename Bt, typename Gt, evdm::GridEL_type g_type>
	Py_Distribution(
		evdm::Distribution<T, Bt, Gt, g_type> const& distrib) :
		m_distrib(distrib) {}

	static Py_Distribution CreatePyDistribFromArray(
		Py_EL_Grid const& mGridEL,
		pybind11::array values
	);
	static Py_Distribution from_dict(pybind11::dict const&);

	Py_EL_Grid getGrid()const;

	template <typename T>
	Py_Distribution as_type_t() const {
		return std::visit([](auto const& distrib) {
			return Py_Distribution(distrib.template  as_type<T>());
			}, m_distrib);
	}
	Py_Distribution copy() const;

	std::string repr()const;
	Py_Distribution as_type(const char* type_n)const;

	pybind11::dict get_object(pybind11::handle self);

	pybind11::array get_array(
		pybind11::handle self, int ptype, bool raw = false);

	pybind11::tuple get_E_distrib(int ptype) const;
	double get_avarage(pybind11::handle Functor, int ptype) const;

	double count(int ptype = -1)const;
	pybind11::tuple get_r_dense(
		pybind11::handle ptypes,
		double r_min, double  r_max, size_t Nr,
		pybind11::handle Nperbin,
		pybind11::handle opt_dense_funct,
		pybind11::handle update_function) const;

	pybind11::tuple plot2o(size_t ptype)const;
	pybind11::tuple plot1o(
		size_t ptype,
		std::string_view m_measure,
		std::string_view LE_Space)const;

	static void add_to_python_module(pybind11::module_& m);
	size_t get_padding() const;
};
Py_Distribution CreatePyDistrib(
	Py_EL_Grid const& mGridEL,
	const char* dtype,
	pybind11::handle Init
);
Py_Distribution CreateDistribFromDict(
	Py_EL_Grid const& mGridEL, pybind11::dict X
);

struct Py_DistribMeasure {
	double p_deg;
	int ptype;
	std::string measure;
	Py_DistribMeasure(std::string measure, double p_deg, int ptype = -1);
	double call(Py_Distribution const& D1, Py_Distribution const& D2);
	std::string repr()const;
};

double compare_distribs(
	Py_Distribution const& D1,
	Py_Distribution const& D2,
	std::string const& measure,
	double p_deg,
	int ptype);