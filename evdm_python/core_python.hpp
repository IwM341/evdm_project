#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <tuple>
#include <set>
#include "value_types.hpp"
#include <ctype.h>

template <typename T>
const char* type_name();

template <typename T, typename Tuple_t>
struct tuple_contains;

template <typename T, typename... Args>
struct tuple_contains<T, std::tuple<Args...>> :
	std::disjunction<std::is_same<T, Args>...> {};

template <typename IndexSequenceType,typename IndexSequenceType::value_type I>
struct index_sequence_contain;

template <typename _T,_T...Is, _T I>
struct index_sequence_contain<std::integer_sequence<_T, Is...>, I> :
	std::disjunction<std::integral_constant<bool,(Is==I)>...> {};

template <typename T>
bool check_np_type(pybind11::array const& np_array) {
	if constexpr (std::is_floating_point_v<T>) {
		if (np_array.dtype().kind() == 'f') {
			return np_array.dtype().itemsize() == sizeof(T);
		} else {
			return false;
		}
	} if constexpr (std::is_integral_v<T>) {
		if (np_array.dtype().kind() == 'i') {
			return np_array.dtype().itemsize() == sizeof(T);
		}
		else {
			return false;
		}
	}
}

	




struct scatter_event_info {
	std::string name;
	size_t ptype_in;
	size_t ptype_out;
	double amount;
	double Nmk;
	bool unique;

	scatter_event_info(pybind11::dict const & _object);
	inline scatter_event_info() :
		name(""), ptype_in(0), ptype_out(0), amount(0), Nmk(0), unique(false) {}
	pybind11::dict to_object() const;
	inline scatter_event_info(
		std::string name, size_t ptype_in, size_t ptype_out,
		double amount, double Nmk,
		bool unique = true):
		name(name), ptype_in(ptype_in), ptype_out(ptype_out),
		amount(amount), Nmk(Nmk),unique(unique){}

	inline bool operator == (const scatter_event_info& event) const {
		return name == event.name && ptype_in == event.ptype_in 
			&& ptype_out == event.ptype_out;
	}
	inline bool operator < (const scatter_event_info& event) const {
		return (name < event.name) || (ptype_in < event.ptype_in)
			|| ((ptype_out < event.ptype_out));
	}
	inline bool operator > (const scatter_event_info& event) const {
		return event < *this;
	}
	inline bool operator <= (const scatter_event_info& event) const {
		return (event < *this) || (event == *this);
	}
	inline bool operator >= (const scatter_event_info& event) const {
		return event <= *this;
	}
	friend inline scatter_event_info linear_combination(
		const scatter_event_info & E1, 
		const scatter_event_info & E2,
		double a1,double a2) {
		return scatter_event_info(
			E1.name,
			E1.ptype_in,
			E1.ptype_out,
			a1 * E1.amount + a2 * E2.amount,
			1 / (a1 * a1 / E1.Nmk + a2 * a2 / E2.Nmk),
			E1.unique || E2.unique
		);
	}
	std::string repr()const;
	static void add_to_python_module(pybind11::module_& m);
};

template <typename Grid_t>
std::vector<size_t> get_N_distrib_from_handle(Grid_t const& Grid, pybind11::handle h) {
	if (pybind11::isinstance<pybind11::array>(h)) {
		auto A = h.cast< pybind11::array_t<size_t>>();
		if (Grid.size() != A.size()) {
			pybind11::index_error("Array of values not match grid size");
		}
		return std::vector<size_t>((size_t*)A.data(), (size_t*)A.data() + A.size());
	}
	if (pybind11::isinstance<pybind11::function>(h)) {
		typedef decltype(Grid.grid()[0].left) T;
		try {
			auto F_el = h.cast < std::function<size_t(T, T)>>();
			(void)(F_el(0, 0));
			std::vector<size_t> NRet(Grid.size());
			auto it = NRet.begin();
			for (auto [de, dl] : Grid) {
				*it = F_el(de.center(), dl.center());
				++it;
			}
			return NRet;
		}
		catch (std::exception& e) {
			auto F_el = h.cast < std::function<size_t(T, T, T, T)>>();
			(void)(F_el(0, 0, 0, 0));
			std::vector<size_t> NRet(Grid.size());
			auto it = NRet.begin();
			for (auto [de, dl] : Grid) {
				auto [e0, e1] = de;
				auto [l0, l1] = dl;
				*it = F_el(e0, e1, l0, l1);
				++it;
			}
			return NRet;
		}
	}
	else {
		pybind11::value_error("expect function or array to make distrib of N");
		return {};
	}
}

