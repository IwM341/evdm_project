#include <evdm/core/core_distrib_sampler.hpp>
#include "../distrib_python.hpp"



typedef  std::variant<std::type_identity<float>, std::type_identity<double>> type_fd_var;

pybind11::array init_samples_impl(type_fd_var, Py_Distribution const& m_distrib, size_t N, size_t seed, double p, double q);
Py_Distribution put_points_impl(type_fd_var, Py_EL_Grid const& Grid, pybind11::array m_states, double total_count);