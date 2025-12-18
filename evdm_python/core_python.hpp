#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <tuple>
#include <set>
#include "value_types.hpp"
#include "grob_python.hpp"
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

struct NpDictSerializator {

	template <typename T>
	struct serializable_s : std::false_type {};

	template <typename T>
	struct deserializable_s : std::false_type {};

	template <typename T, typename...Args>
	struct serializable_s<
		std::vector<T, Args...>
	> : std::is_arithmetic<T> {};

	template <typename T, typename...Args>
	struct deserializable_s<
		std::vector<T, Args...>
	> : std::is_arithmetic<T> {};


	template <typename T>
	pybind11::object MakePrimitive(T const& value) {
		return pybind11::cast(
			value,
			pybind11::return_value_policy::automatic_reference
		);
	}
	template <typename Keys_t, typename Values_t>
	pybind11::object MakeDict(Keys_t&& keys, Values_t&& values) {
		pybind11::dict m_dict;
		for (size_t i = 0; i < keys.size(); ++i) {
			if constexpr (
				std::is_same_v<
				std::decay_t<decltype(keys[i])>,
				std::string
				>
				) {
				m_dict[keys[i].c_str()] = values[i];
			}
			else if constexpr (
				std::is_same_v<
				std::decay_t<decltype(keys[i])>,
				std::string_view
				>
				) {
				std::string X(keys[i].begin(), keys[i].end());
				m_dict[X.c_str()] = values[i];
			}
			else {
				m_dict[keys[i]] = values[i];
			}
		}
		return m_dict;
	}

	template <typename Values_t>
	pybind11::object MakeArray(Values_t&& values) {
		typedef std::decay_t<decltype(values[0])> value_type;
		if constexpr (std::is_fundamental_v<value_type>) {
			return pybind11::cast(make_py_array(values));
		}
		else {
			pybind11::list L;
			size_t N = values.size();
			for (size_t i = 0; i < N; ++i) {
				L.append(values[i]);
			}
			return L;
		}
	}

	template <typename...Ts>
	pybind11::object Make(std::vector<Ts...> const& value) {
		return make_py_array(value);
	}

	template <typename T, typename...Args>
	std::vector<T, Args...> get_impl(
		std::type_identity<std::vector<T, Args...>>,
		pybind11::handle Obj)
	{
		try {
			pybind11::array_t<T> m_array = Obj.cast<pybind11::array_t<T>>();
			if (m_array.ndim() != 1) {
				throw pybind11::index_error("number of dimentions in array should be 1 to cast to std::vector");
			}
			const T* _ptr = m_array.data();
			size_t _stride = m_array.strides()[0] / sizeof(T);
			size_t _size = m_array.shape()[0];
			std::vector<T, Args...> Ret(_size);
			for (size_t i = 0; i < _size; ++i) {
				Ret[i] = *(_ptr + i * _stride);
			}
			return Ret;
		}
		catch (pybind11::cast_error&) {

		}

		std::vector<T, Args...> Ret;
		for (auto it = Obj.begin(); it != Obj.end(); ++it) {
			Ret.push_back(it->cast<T>());
		}

		return Ret;
	}


	template <typename T>
	T GetPrimitive(pybind11::handle Obj) {
		return pybind11::cast<T>(Obj);
	}
	template <typename T>
	T Get(pybind11::handle Obj) {
		return get_impl(std::type_identity<T>{}, Obj);
	}


	auto BeginArray(pybind11::handle Obj) {
		return Obj.begin();
	}
	auto EndArray(pybind11::handle Obj) {
		return Obj.end();
	}

	auto BeginDict(pybind11::handle Obj) {
		return pybind11::cast<pybind11::dict>(Obj).begin();
	}
	auto EndDict(pybind11::handle Obj) {
		return pybind11::cast<pybind11::dict>(Obj).end();
	}
	template <typename Iter_type>
	auto GetKey(Iter_type it) {
		return it->first.template cast<std::string_view>();
	}
	template <typename Iter_type>
	auto GetValue(Iter_type it) {
		return it->second;
	}
	template <typename Iter_type>
	auto GetItem(Iter_type it) {
		return *it;
	}
	template <typename T>
	pybind11::handle GetPrimitive(pybind11::handle Obj, T& value) {
		return value = Obj.cast<T>();
		grob::GridUniform<T>::Serialize;
	}
};