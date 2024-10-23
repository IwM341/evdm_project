#ifndef DYNAMIC_PYTHON_HPP
#define DYNAMIC_PYTHON_HPP
#include "core_python.hpp"
#include <evdm/dynamics/dynamics.hpp>
#include <evdm/utils/prng.hpp>

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
	double eval(double y,double v2_T = 0)const {
		return eval_slow(y, v2_T);
	}
	pybind11::object _opt_functur;
};
struct Py_ScatterEvent : public evdm::ScatterEvent {
	std::string name;
	bool unique;
	Py_ScatterEvent(
		pybind11::handle n_e, Py_ScatterFactor const & _factor, 
		const char * name,bool unique);
	static void add_to_python_module(pybind11::module_& m);
};





pybind11::tuple Py_CaptureProcess(
	Py_Capture& CaptAccum,
	int ptype_in,
	int ptype_out,
	float M_DM,
	float deltaM,
	float NucleiM,
	const Py_ScatterEvent& sc_event,
	float body_halo_v,
	float dm_v_disp,
	size_t Nmk,
	float r_pow = 2, // the r sistribution
	float weight = 1, // result will be multiplyed by weight
	size_t seed = evdm::default_seed
	);
void Py_ScatterProcess(
	Py_Matrix & ScatterAccum,
	int ptype_in, int ptype_out, 
	float M_DM,
	float deltaM,
	float NucleiM,
	const Py_ScatterEvent& sc_event,
	pybind11::handle Nmk_v,
	pybind11::kwargs ExtraArgs
);

std::string make_compare_sc_event(
	const Py_ScatterEvent& sc_event,
	int ptype_in, int ptype_out,
	const std::vector< scatter_event_info>& events
);

void add_pycapture_to_python_module(pybind11::module_& m);
void add_pyscatter_to_python_module(pybind11::module_& m);


inline void add_scatter_funcs_to_python_module(pybind11::module_& m) {
	Py_ScatterFactor::add_to_python_module(m);
	Py_ScatterEvent::add_to_python_module(m);
	add_pycapture_to_python_module(m);
	add_pyscatter_to_python_module(m);
}


#endif//DYNAMIC_PYTHON_HPP