#include "core_python.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <evdm/core.hpp>
#include <pybind11/stl.h>
#include "debugdef.hpp"
#include "grob_python.hpp"


size_t Py_BodyModel::size()const {
	return std::visit([](auto const& B) {return B->Phi.Grid.size(); }, m_body);
}

void Py_BodyModel::setTemp(pybind11::handle mTemp){
	std::visit([&mTemp](auto & mb) {
		if (pybind11::array_t<double>::check_(mTemp)) {
			mb->setTemp(
				make_py_array_slice(
					mTemp.cast< pybind11::array_t<double>>()
				)
			);
		}
		else if (pybind11::array_t<float>::check_(mTemp)) {
			mb->setTemp(
				make_py_array_slice(
					mTemp.cast< pybind11::array_t<float >>()
				)
			);
		}
		}, m_body);
	
}
Py_BodyModel::Py_BodyModel(pybind11::array_t<double> const& RhoValues, double Velocity, pybind11::handle Temp)
:m_body(evdm::load_body_model<double,evdm::forward_shared>(RhoValues.data(),RhoValues.size(), Velocity))
{
	setTemp(Temp);
}

Py_BodyModel::Py_BodyModel(pybind11::array_t<float> const& RhoValues, float Velocity, pybind11::handle Temp)
	:m_body(evdm::load_body_model<double, evdm::forward_shared>(RhoValues.data(), RhoValues.size(), Velocity))
{
	setTemp(Temp);
}

template <typename T>
struct T_t{
	typedef T type;
};
Py_BodyModel::Py_BodyModel(pybind11::handle const& Rho,
	std::string_view dtype,
	std::optional<size_t> _size, 
	double Velocity, pybind11::handle Temp)
{
	//DEBUG1(__PRETTY_FUNCTION__);

	size_t m_size = 0;
	bool is_function = false;

	try {
		m_size = pybind11::len(Rho);
	} catch(...){
		is_function = true;
	}

	
	if (!is_function) {
		auto ExtractList = [&](auto var_type) {
			typedef typename decltype(var_type)::type T;
			std::vector<T> values(m_size);
			size_t i = 0;
			for (auto item : Rho) {
				values[i] = item.cast<T>();
				++i;
			}
			return evdm::load_body_model<evdm::forward_shared>(std::move(values), Velocity);
		};
		if (dtype == "double") {
			m_body = ExtractList(T_t<double>{});
			return;
		}
		else if (dtype == "float") {
			m_body = ExtractList(T_t<float>{});
			return;
		}
		else {
			throw pybind11::type_error("wrong data type: " + std::string(dtype) + ", expect float or double");
		}
	}
	else {
		if (_size.has_value()) {
			m_size = _size.value();
		}
		else {
			throw pybind11::value_error("for callable Rho function need size");
		}
		auto Fd = [&](auto x) {
			return Rho(x).cast<double>();
		};
		auto Ff = [&](auto x) {
			return Rho(x).cast<float>();
		};
		if (dtype == "double") {
			m_body = evdm::make_body_model<double,evdm::forward_shared>(Fd, m_size, Velocity);
			return;
		}
		else if (dtype == "float") {
			m_body = evdm::make_body_model<float, evdm::forward_shared>(Ff, m_size, Velocity);
			return;
		}
		else {
			throw pybind11::type_error("wrong data type: " + std::string(dtype) + ", expect float or double");
		}

	}
	setTemp(Temp);
}

const char* Py_BodyModel::dtype() const
{
	if(m_body.index() == 0)
		return "float";
	else if (m_body.index() == 1) {
		return "double";
	}
	else {
		return "notype";
	}
}

std::string Py_BodyModel::repr() const
{
	return std::visit([](const auto & B) {
			const char* _tn = type_name<decltype(B->get_vtype())>();
			return "Body( size = " + std::to_string(B->Phi.size()) + ", dtype = " + _tn + ")";
		}, m_body);
}

pybind11::tuple Py_BodyModel::getRho() const
{
	return std::visit([](auto const& B) {
		return make_python_function_1D(B->Rho);
	},m_body);
}

pybind11::tuple Py_BodyModel::getPhi() const
{
	return std::visit([](auto const& B) {
			return make_python_function_1D(B->Phi);
		}, m_body);
}

pybind11::tuple Py_BodyModel::getM() const
{
	return std::visit([](auto const& B) {
		return make_python_function_1D(B->M());
	}, m_body);
}

pybind11::tuple Py_BodyModel::getQ() const
{
	return std::visit([](auto const& B) {
		return make_python_function_1D(B->Q);
		}, m_body);
}

Py_BodyModel Create(pybind11::handle Rho, std::optional<size_t> _size, 
	std::string_view dtype,double Velocity, pybind11::handle mTemp = pybind11::none())
{	
	if (pybind11::array_t<double>::check_(Rho)) {
		return Py_BodyModel(Rho.cast<pybind11::array_t<double>>(), Velocity, mTemp);
	}
	else if (pybind11::array_t<float>::check_(Rho)) {
		return Py_BodyModel(Rho.cast<pybind11::array_t<float>>(), (float)Velocity, mTemp);
	}
	else {
		return Py_BodyModel(Rho, dtype, _size, Velocity, mTemp);
	}
	
}

void Py_BodyModel::add_to_python_module(pybind11::module_& m)
{ 
	
	namespace py = pybind11;
	using pyobj_ref = pybind11::object const&;
	py::class_<Py_BodyModel>(m, "Body")
		.def(
			py::init([](pyobj_ref _o, std::optional<size_t> const& _size,
				double Velocity, std::string_view dtype, pybind11::handle mTemp)
				{return Create(_o, _size, dtype, Velocity, mTemp); }),
			"constructs Body from Rho\n"
			"Rho -- array of values of mass density or function [0,1]->real\n"
			"size -- if Rho is func then size is the number of points\n"
			"velocity -- velocity of body inside halo\n"
			"dtype -- string, float or double\n"
			"Temp -- optional temperature array, default - 0\n",
			py::arg("Rho"),
			py::arg_v("size", std::nullopt),
			py::arg_v("velocity", 0.7e-3),
			py::arg_v("dtype", std::string("")),
			py::arg_v("Temp", pybind11::none())
		)
		.def("setTemp", &Py_BodyModel::setTemp, py::arg_v("Temp", pybind11::none()))
		.def("__repr__", &Py_BodyModel::repr)
		.def_property_readonly("size", &Py_BodyModel::size)
		.def_property_readonly("dtype", &Py_BodyModel::dtype)
		.def_property_readonly("phi", &Py_BodyModel::getPhi)
		.def_property_readonly("rho", &Py_BodyModel::getRho)
		.def_property_readonly("M", &Py_BodyModel::getM)
		.def_property_readonly("Q", &Py_BodyModel::getQ);
}

