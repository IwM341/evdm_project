// Combined template specializations
// Generated automatically - do not edit

#pragma once

#include <evdm/core/core_dynamics_scatter.hpp>

#include <evdm/utils/prng.hpp>

#include <evdm/measure.hpp>

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 0>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 1>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< float, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< float, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, float, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, float, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, float, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, float, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< float> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCUU>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCUU>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	size_t const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Naive m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Full m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8 m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_NoTherm m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

void ScatterImpl(
	std::integral_constant<int, 2>,
	evdm::GridMatrix< double, double, double, evdm::GridEL_type::GridCVV>& m_mat,
	evdm::Distribution< double, double, double, evdm::GridEL_type::GridCVV>& m_evp,
	double zero_val, int ptype_in, int ptype_out, size_t seed,
	evdm::ScatterEvent const& sc_event, evdm::ThermGaussGenerator_Soft8_Treshold m_therm_gen,
	evdm::measure_dEpdlq< double> m_el_measure,
	float M_DM, float deltaM, float NucleiM,
	std::vector<size_t> const& Nmk, size_t Nmk_per_traj,
	double weight, evdm::progress_omp_function<> ProgFunc
);
// ======================

