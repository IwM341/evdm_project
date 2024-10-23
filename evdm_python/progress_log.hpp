#pragma once

#include <pybind11/pybind11.h>
#include <evdm/utils/progress_bar.hpp>

inline evdm::progress_omp_function<> make_progress_func(
	pybind11::handle update_function
)
{
	evdm::progress_omp_function<> ProgF;
	if (!update_function.is_none()) {
		if (pybind11::hasattr(update_function, "update")) {
			ProgF = evdm::progress_omp_function(
				[update_function](int m_curr, int m_full) {
					update_function.attr("update")(m_curr, m_full);
				}
			);
		}
		else if (pybind11::hasattr(update_function, "__call__")) {
			ProgF = evdm::progress_omp_function(
				[update_function](int m_curr, int m_full) {
					update_function.operator()(m_curr, m_full);
				}
			);
		}
		else {
			throw pybind11::value_error(
				"provided progress bar class don't contain"
				"'update' functon or not callable");
		}
	}
	return ProgF;
}