#include "../scatter_impl.hpp"

void ScatterImpl_Args::apply(std::integral_constant<int, evdm::AlgolDiffuse> arg) {
	apply_impl(arg);
}