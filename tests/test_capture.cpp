
//#define HISTO_OPTIM
#include <evdm/utils/mc.hpp>
#include <iostream>
#include <evdm/core/core_distrib.hpp>
#include <evdm/core/core_dynamics_capture.hpp>
#include <grob/const_container.hpp>
#include <chrono>

auto time_now() {
	return std::chrono::high_resolution_clock::now();
}

template <typename T>
auto time_count_micros(T start, T stop) {
	return std::chrono::duration<double, std::micro>(stop - start).count();
}
template <typename T>
auto time_count_millis(T start, T stop) {
	return std::chrono::duration<double, std::milli>(stop - start).count();
}
template <typename T>
auto time_count_sec(T start, T stop) {
	return std::chrono::duration<double>(stop - start).count();
}


int main() {
	auto r_dens = [](float x) {return 1; };
	evdm::BodyModel<float> B = evdm::make_body_model<float,evdm::forward_shared>(
		r_dens, 1001, 0.002);
	auto Grid_1 = evdm::Grid_types<float, evdm::GridEL_type::GridCVV>::construct(
		1, -B->Phi[0], 100, [](auto t_e) {return 0.1+t_e; }, [](auto t_e) {return 100; }, [](auto t_e,auto t_l) {return 1; }
	);
	auto gridel = evdm::EL_Grid<float, float, evdm::GridEL_type::GridCVV>(
		B, Grid_1,200
	);
	evdm::Distribution<float,float,float, evdm::GridEL_type::GridCVV> m_distrib(
		evdm::make_Distribution<float>(gridel, false)
	);
	
	auto n_p = grob::as_container([](size_t i) {return 1; }, 101);

	auto event = evdm::ScatterEvent(
		n_p.size(), n_p,
		evdm::FormFactor_t(false, 1, std::vector<float>({1,0.1,0.1,0.1}))
	);

	size_t N = 1000000;
	auto G = evdm::xorshift32f<float>{};
	
	for (size_t j = 0; j < 20; ++j) {	
		auto t0 = time_now();
		G.set_seed(123);
		auto [C, dC] = evdm::Capture(
			m_distrib, 0, event, G, 100, 0, 10, 0.001, 0.001, false, 0.001 * 4, 2, N, 1
		);
		auto t1 = time_now();
		std::cout << "time: " << time_count_millis(t0, t1) << std::endl;
		std::cout << "Capt = " << m_distrib.count(0) << std::endl;
		std::cout << "dCapt = " << dC / sqrt(N) << std::endl;
		m_distrib.values() *= 0;
	}


	

	return 0;
}