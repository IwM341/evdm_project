#include <evdm/utils/mc.hpp>
#include <iostream>
#include <evdm/core/core_distrib.hpp>
#include <evdm/core/core_dynamics_capture.hpp>
#include <grob/const_container.hpp>
int main() {
	auto r_dens = [](float x) {return 1; };
	evdm::BodyModel<float> B = evdm::make_body_model<float,evdm::forward_shared>(
		r_dens, 1001, 0.002);
	auto Grid_1 = evdm::Grid_types<float, evdm::GridEL_type::GridCUU>::construct(
		1, -B->Phi[0], 100, [](auto t_e) {return 100; }
	);
	auto gridel = evdm::EL_Grid<float, float, evdm::GridEL_type::GridCUU>(
		B, Grid_1,200
	);
	auto m_distrib = evdm::make_Distribution<float>(gridel, false);
	
	auto n_p = grob::as_container([](size_t i) {return 1; }, 101);

	auto event = evdm::ScatterEvent(
		n_p.size(), n_p,
		evdm::FormFactor_t(false, 1, std::vector<float>({1,0.1,0.1,0.1}))
	);

	size_t N = 1000000;
	auto G = evdm::xorshift32f<float>{};
	auto [C, dC] = evdm::Capture(
		m_distrib, 0, event, G, 100, 0, 10, 0.001, 0.001, 2, N, 1
	);
	std::cout << "Capt = " << m_distrib.count(0) << std::endl;
	std::cout << "dCapt = " << dC/sqrt(N) << std::endl;

	return 0;
}