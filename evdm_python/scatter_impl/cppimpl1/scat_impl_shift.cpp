#include "../scatter_impl.hpp"

void ScatterImpl_Args::apply(std::integral_constant<int, evdm::AlgolShift> arg) {
	apply_impl(arg);
}