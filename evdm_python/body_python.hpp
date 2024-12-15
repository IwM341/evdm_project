#pragma once
#include "core_python.hpp"


struct Py_BodyModel {
	BodyModel_Variant_t	m_body;
#ifdef BODY_MODEL_USE_DOUBLE
	using max_vtype = double;
#else
	using max_vtype = float;
#endif // DISTRIB_USE_DOUBLE
#ifdef BODY_MODEL_USE_FLOAT
	using min_vtype = float;
#else
	using min_vtype = double;
#endif // DISTRIB_USE_DOUBLE

	template <typename T>
	Py_BodyModel(evdm::BodyModel<T> m_body) :m_body(std::move(m_body)) {}

#ifdef BODY_MODEL_USE_DOUBLE
	Py_BodyModel(pybind11::array_t<double> const& RhoValues, double Velocity,
		pybind11::handle Temp = pybind11::none());
#endif
#ifdef BODY_MODEL_USE_FLOAT
	Py_BodyModel(pybind11::array_t<float> const& RhoValues, float Velocity,
		pybind11::handle Temp = pybind11::none());
#endif
	Py_BodyModel(pybind11::handle  const& RhoObject, double Velocity,
		std::string_view dtype,
		pybind11::handle Temp = pybind11::none());

	void setTemp(pybind11::handle mTemp);

	double get_vesc() const;
	const char* dtype()const;

	template <typename U>
	Py_BodyModel as_type_t() const {
		return std::visit([](auto const& Bptr) {
			return Py_BodyModel(
				std::make_shared<evdm::Body<U>>(Bptr->template as_type<U>())
			);
			}, m_body
		);
	}
	Py_BodyModel as_type(const char* type_name) const;

	std::string repr() const;
	size_t size()const;
	pybind11::tuple getRho()const;
	pybind11::tuple getPhi()const;
	pybind11::tuple getM()const;
	pybind11::tuple getQ()const;


	pybind11::dict get_object(pybind11::handle self) const;
	static Py_BodyModel from_dict(pybind11::dict const& Object);

	static void add_to_python_module(pybind11::module_& m);

};