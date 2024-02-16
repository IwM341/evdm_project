#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Eigen>
#include <evdm/core.hpp>
#include "core_python.hpp"
#include "debugdef.hpp"
#include "dynamic_python.hpp"

namespace py = pybind11;

struct Pet {
    Pet(const std::string& name) : name(name) { }
    std::string name;
};

struct Dog : Pet {
    using Pet::Pet;
    Dog(const Pet& _m_pet) :Pet(_m_pet) {}
    std::string bark() const { return "woof!"; }
};

PYBIND11_MODULE(pyevdm, m)
{
    m.doc() = "pyevdm module with classes:\n"
        "Body -- representation of gravity potential of body\n"
        "GridEL -- grid in ptypes, E,L variables\n"
        "Distrib -- class of distribution of particles\n"
        "Capture -- Distribution of captured particles with event info\n"
        "Matrix -- class of Scatter Matrix\n"
        "ScatterFactor -- form factor appered in collisions\n"
        "ScatterEvent -- concentration of scatter target + form factor (ScatterFactor)";

    Py_BodyModel::add_to_python_module(m);
    Py_EL_Grid::add_to_python_module(m);
    Py_Distribution::add_to_python_module(m);
    scatter_event_info::add_to_python_module(m);
    Py_Capture::add_to_python_module(m);
    Py_Matrix::add_to_python_module(m);
    Py_ScatterFactor::add_to_python_module(m);
    Py_ScatterEvent::add_to_python_module(m);

    py::class_<Pet>(m, "Pet")
        .def(py::init<const std::string&>(),
            "Documentation for class Pet")
        .def_readwrite("name", &Pet::name);

    // Method 1: template parameter:
    py::class_<Dog, Pet /* <- specify C++ parent type */>(m, "Dog")
        .def(py::init<const Pet&>(),
            "Pet Constructor"
        )
        .def("bark", &Dog::bark);
}
