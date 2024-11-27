#pragma once
#include "distrib_python.hpp"

struct Py_Capture : public Py_Distribution {
	using Py_Distribution::Py_Distribution;
	inline Py_Capture(Py_Distribution _distrib) :
		Py_Distribution(std::move(_distrib)) {}

	static Py_Capture from_dict(pybind11::dict const&);

	std::vector<scatter_event_info> events;
	static void add_to_python_module(pybind11::module_& m);
	inline std::vector<scatter_event_info>const& get_events()const { return events; }

	Py_Capture(Py_EL_Grid const& m_grid,
		pybind11::dict m_py_object);
	pybind11::dict get_object(pybind11::handle self);
	Py_Capture copy() const;

	template <typename T>
	inline Py_Capture as_type_t() const {
		return Py_Capture(Py_Distribution::as_type_t<T>());
	}

	inline Py_Capture as_type(const char* dtype) const;

	/// <summary>
	/// summ two captures WITH DIFFERENT EVENTS
	/// </summary>
	/// <param name="_another"></param>
	Py_Capture& add(Py_Capture const& _another);
	Py_Capture operator_plus(Py_Capture const& _another) const;
};