#ifndef DYNAMIC_PYTHON_HPP
#define DYNAMIC_PYTHON_HPP
#include "core_python.hpp"
#include <evdm/dynamics/dynamics.hpp>

struct Py_ScatterFactor : public evdm::FormFactor_t {
	typedef evdm::FormFactor_t Base;
	using Base::Base;
	static void add_to_python_module(pybind11::module_& m);
	inline  evdm::FormFactor_t& to_ff() {
		return *this;
	}
	inline  const evdm::FormFactor_t& to_ff()const {
		return *this;
	}
};
Py_ScatterFactor qexp_factor(float b, bool y_inv,
	pybind11::array_t<float> P_0_coeffs, pybind11::array_t<float> P_V_coeffs);
Py_ScatterFactor qexp_factor(float b, bool y_inv,
	pybind11::array_t<float> P_0_coeffs);

struct Py_ScatterEvent : public evdm::ScatterEvent {
	std::string _name;
	bool unique;
	Py_ScatterEvent(pybind11::array n_e, Py_ScatterFactor const & _factor, std::string_view name,bool unique);
	static void add_to_python_module(pybind11::module_& m);
};


struct Py_TrajectoryInfo {
	
	static void add_to_python_module(pybind11::module_& m);
};

struct Py_AnnihilationInfo {
	static void add_to_python_module(pybind11::module_& m);
};



pybind11::tuple Py_CaptureProccess(
	Py_Capture& CaptAccum,
	int ptype,
	float M_DM,
	float deltaM,
	float NucleiM,
	const Py_ScatterEvent& sc_event,
	float body_halo_v,
	float dm_v_disp,
	size_t Nmk,
	std::string gen_dtype = "float",
	float r_pow = 2, // the r sistribution
	float weight = 1 // result will be multiplyed by weught
	);
void Scatter(std::vector<std::pair<int,int>> ptypes_in_out, 
	const Py_ScatterEvent& sc_event, 
	Py_TrajectoryInfo const& TrajInfo , 
	Py_Matrix& CaptAccum);

void add_scatter_funcs_to_python_module(pybind11::module_& m);


#endif//DYNAMIC_PYTHON_HPP