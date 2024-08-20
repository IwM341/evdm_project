#include "core_python.hpp"

scatter_event_info::scatter_event_info(
	pybind11::dict const& _object) :
	scatter_event_info(
		_object["name"].cast<std::string>(),
		_object["ptype_in"].cast<size_t>(),
		_object["ptype_out"].cast<size_t>(),
		_object["amount"].cast<double>(),
		_object["Nmk"].cast<double>(),
		_object["unique"].cast<bool>()
	) {}

pybind11::dict scatter_event_info::to_object() const
{
	using namespace pybind11::literals;
	return pybind11::dict(
		"name"_a = name,
		"ptype_in"_a = ptype_in,
		"ptype_out"_a = ptype_out,
		"amount"_a = amount,
		"Nmk"_a = Nmk,
		"unique"_a = unique
	);
}
std::string scatter_event_info::repr()const
{
	return "(" + name + " + W_m_" + std::to_string(ptype_in)
		+ "->" + std::to_string(ptype_in) + ": " +
		std::to_string(Nmk) +
		" [" + std::to_string(amount) + "]" + ")";
}
