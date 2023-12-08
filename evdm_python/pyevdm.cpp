#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Eigen>
#include <evdm/core.hpp>
#include "core_python.hpp"
#include "debugdef.hpp"


namespace py = pybind11;

double summ_obj(py::object const& m_list){
    std::cout << __PRETTY_FUNCTION__<<std::endl;
    double summ = 0;
    for(auto item: m_list){
        summ+= item.cast<double>();
    }
    return summ;
}

struct Some {
    Some(std::string S) :S(S) {}

    std::string S;
};

Some factory(pybind11::object const& _o) {
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout << _o.get_type().cast<std::string>() << std::endl;
    return Some("123");
    pybind11::array_t<float> G;
    G.mutable_data();
}

PYBIND11_MODULE(pyevdm, m)
{
    m.doc() = "pyevdm module";
    //std::cout << __PRETTY_FUNCTION__<<std::endl;
    //np::initialize();
    m.def("add", [](double a,double b){return a+b;},
            "A function that adds two numbers",
            py::arg("i"), py::arg("j"));

    //pybind11::class_<Some>(m, "Some")
    //    .def(pybind11::init([](pybind11::object const& _o) {return factory(_o); }));

    Py_BodyModel::add_to_python_module(m);
    Py_EL_Grid::add_to_python_module(m);
    //m.def("factory", factory);
}
