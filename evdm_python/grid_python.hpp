#pragma once
#include "body_python.hpp"


struct Py_EL_Grid {
	ELGrid_Variant_t m_grid;

	/// @brief size of inner part of grid
	size_t size()const;
	size_t size_e() const;
	size_t ptypes()const;
	const char* dtype()const;
	std::string repr() const;
	std::string printd1() const;



	template <typename Bt, typename Gt, evdm::GridEL_type _m_type>
	Py_EL_Grid(evdm::EL_Grid<Bt, Gt, _m_type> const& _m_grid) :m_grid(_m_grid) {}
	template <typename T, typename FuncType_NL, typename GEL_t, typename RhoE_ft, typename RhoL_ft>
	inline Py_EL_Grid(evdm::BodyModel<T> const& BM,
		size_t ptypes, size_t Ne,
		FuncType_NL&& Nl_func,
		RhoE_ft&& RhoE_func_or_false,
		RhoL_ft&& RhoL_func_or_false,
		std::type_identity<GEL_t>,
		std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>,
		evdm::TrajPoolInitParams_t TrajPoolInit = evdm::TrajPoolInitParams_t()) :
		m_grid(
			evdm::EL_Grid <T, GEL_t, evdm::GridEL_type::GridCUU>(BM,
				evdm::Grid_types<GEL_t, evdm::GridEL_type::GridCUU>::construct(
					ptypes, (T)-BM->Phi[0], Ne, Nl_func
				),
				2 * Ne, TrajPoolInit
				)
		) {}
	template <typename T, typename FuncType_NL, typename GEL_t,
		typename RhoE_ft, typename RhoL_ft>
	inline Py_EL_Grid(evdm::BodyModel<T> const& BM,
		size_t ptypes, size_t Ne,
		FuncType_NL&& Nl_func,
		RhoE_ft&& RhoE_func_or_false,
		RhoL_ft&& RhoL_func_or_false,
		std::type_identity<GEL_t>,
		std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>,
		evdm::TrajPoolInitParams_t TrajPoolInit = evdm::TrajPoolInitParams_t()) :
		m_grid(
			evdm::EL_Grid <T, GEL_t, evdm::GridEL_type::GridCVV>(BM,
				evdm::Grid_types <GEL_t, evdm::GridEL_type::GridCVV>::construct(
					ptypes, (T)-BM->Phi[0], Ne, RhoE_func_or_false,
					Nl_func, RhoL_func_or_false
				), 4 * Ne, TrajPoolInit
				)
		) {}
	Py_EL_Grid refine(size_t Ne, size_t Nl) const;

	static Py_EL_Grid from_dict(Py_BodyModel BM, pybind11::dict const& grid_dict);
	static Py_EL_Grid from_dict1(pybind11::dict const& grid_dict);

	pybind11::dict get_object(pybind11::handle self);

	Py_BodyModel getBody()const;
	pybind11::array getPlot(bool is_internal) const;
	pybind11::tuple getLE() const;

	pybind11::tuple getRminRmax(double e, double l_undim) const;

	pybind11::list getTrajTFuncs(pybind11::handle output_param_name, bool LE_mult = false) const;
	std::function<double(double, double)> getTraj_callFunc(pybind11::handle output_param_name) const;
	std::function<double(double)> getLE_callFunc() const;

	static void add_to_python_module(pybind11::module_& m);
	std::variant<size_t, std::pair<size_t, size_t>> get_index(
		double e, double L, bool linear, bool hidden)const;
	size_t get_e_index(double e)const;
	pybind11::tuple get_el(
		std::variant<size_t, std::pair<size_t, size_t>> index)const;
	pybind11::dict get_el_all(
		std::variant<size_t, std::pair<size_t, size_t>> index)const;
	pybind11::array_t<size_t> indexes(pybind11::handle) const;

	pybind11::array get_E_array(pybind11::handle self)const;
	pybind11::array get_L_array(
		pybind11::handle, int index, bool hidden
	)const;



	pybind11::dict _get_traj_all(double e, double L);

	pybind11::dict _get_traj_all_inter(double e, double L);
};

