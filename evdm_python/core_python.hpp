#ifndef CORE_PYTHON_HPP
#define CORE_PYTHON_HPP
#include <evdm/core.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


template <typename T>
const char* type_name();


struct Py_BodyModel {
	std::variant<
		evdm::BodyModel<float>,
		evdm::BodyModel<double>
	> m_body;

	Py_BodyModel(pybind11::array_t<double> const & RhoValues);
	Py_BodyModel(pybind11::array_t<float> const& RhoValues);
	Py_BodyModel(pybind11::handle  const& RhoObject,
				std::string const& dtype,
				std::optional<size_t> _size);


	const char* dtype()const;
	std::string repr() const;
	size_t size()const;
	pybind11::tuple getRho()const;
	pybind11::tuple getPhi()const;
	pybind11::tuple getM()const;
	pybind11::tuple getQ()const;
	
	static Py_BodyModel Create(pybind11::handle  RhoObject,
		std::string const& dtype = "float",
		std::optional<size_t> _size = std::nullopt);

	static void add_to_python_module(pybind11::module_& m);
};

struct Py_EL_Grid {
	std::variant<
		//           B vtype, GEL vtype
		evdm::EL_Grid<float,float>,
		evdm::EL_Grid<float, double>,
		evdm::EL_Grid<double, float>,
		evdm::EL_Grid<double, double>
	> m_grid;
	
	size_t size()const;
	size_t ptypes()const;
	const char* dtype()const;
	std::string repr() const;
	std::string printd1() const;
	
	template <typename T> struct m_type_marker 
	{
		typedef T type;
	};

	template <typename T,typename FuncType_NL,typename GEL_t, typename RhoE_ft, typename RhoL_ft>
	inline Py_EL_Grid(evdm::BodyModel<T> const& BM,
		size_t ptypes, size_t Ne,
		FuncType_NL&& Nl_func,
		RhoE_ft && RhoE_func_or_false,
		RhoL_ft && RhoL_func_or_false,
		m_type_marker<GEL_t>,
		std::integral_constant<bool,true> is_CUU):
			m_grid(
				evdm::EL_Grid <T, GEL_t>(BM,
					evdm::EL_Grid <T, GEL_t>::GridEL_t::grid_var_CUU(
						ptypes, (T) - BM->Phi[0], Ne, Nl_func
					),
					2 * Ne
				)
			){}
	template <typename T, typename FuncType_NL, typename GEL_t, typename RhoE_ft, typename RhoL_ft>
	inline Py_EL_Grid(evdm::BodyModel<T> const& BM,
		size_t ptypes, size_t Ne,
		FuncType_NL&& Nl_func,
		RhoE_ft&& RhoE_func_or_false,
		RhoL_ft&& RhoL_func_or_false,
		m_type_marker<GEL_t>,
		std::integral_constant<bool, false> is_CUU) :
			m_grid(
				evdm::EL_Grid <T, GEL_t>(BM,
					evdm::EL_Grid <T, GEL_t>::GridEL_t::grid_var_CVV(
						ptypes, (T) - BM->Phi[0], Ne, RhoE_func_or_false,
						Nl_func, RhoL_func_or_false
					),2*Ne
				)
			) {}

	pybind11::tuple getPlot(bool is_internal) const;
	pybind11::tuple getLE() const;
	static void add_to_python_module(pybind11::module_& m);

};

#endif//CORE_PYTHON_HPP