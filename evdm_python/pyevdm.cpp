#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Eigen>
#include <evdm/core.hpp>
#include "core_python.hpp"
#include "debugdef.hpp"
#include "dynamic_python.hpp"

namespace py = pybind11;


PYBIND11_MODULE(pyevdm, m)
{
    m.doc() = "pyevdm module with classes:\n"
        "Body : representation of gravity potential of body\n"
        "GridEL : grid in ptypes, E,L variables\n"
        "Distrib : class of distribution of particles\n"
        "Capture : Distribution of captured particles with event info\n"
        "Matrix : class of Scatter Matrix\n"
        "ScatterFactor : form factor appered in collisions\n"
        "ScatterEvent : concentration of scatter target + form factor (ScatterFactor)";

    Py_BodyModel::add_to_python_module(m);
    Py_EL_Grid::add_to_python_module(m);
    Py_Distribution::add_to_python_module(m);
    scatter_event_info::add_to_python_module(m);
    Py_Capture::add_to_python_module(m);
    Py_Matrix::add_to_python_module(m);
 

    add_scatter_funcs_to_python_module(m);
}
