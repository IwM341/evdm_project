#include "dynamic_python.hpp"
#include <evdm/core/core_dynamics_scatter.hpp>
#include <evdm/utils/prng.hpp>
#include "grob_python.hpp"
#include "progress_log.hpp"
#include "value_types.hpp"
#include <evdm/measure.hpp>
#include "matrix_python.hpp"

void Py_ScatterProcess(
	Py_Matrix& ScatterAccum,
	int ptype_in, int ptype_out,
	float M_DM,
	float deltaM,
	float NucleiM,
	const Py_ScatterEvent& sc_event,
	pybind11::handle Nmk_v,
	pybind11::kwargs ExtraArgs
) {
	size_t Nmk_per_traj = pyget<size_t>(1, ExtraArgs, "Nmk_traj");
	bool count_evap = true;// pyget<bool>(false, ExtraArgs, "count_evap");
	size_t seed = pyget<size_t>(evdm::default_seed, ExtraArgs, "seed");
	pybind11::handle update_function = pyget<pybind11::handle>(
		pybind11::none(), ExtraArgs,"bar"
	);

	std::string methodS = pyget<std::string>("notherm", ExtraArgs, "method");
	double weight = pyget<double>(1, ExtraArgs, "weight");
	
	auto m_compare =
		make_compare_sc_event(sc_event, ptype_in, ptype_out, ScatterAccum.events);
	if (!m_compare.empty()) {
		throw pybind11::value_error(m_compare);
	}
	if (seed == 0) {
		seed = std::numeric_limits<size_t>::max();
	}
	auto ProgFunc = make_progress_func(update_function);
	size_t Nmk = 0;
	ScatterMethodVariant_t ThermGenVar = std::visit(
		[]<class T>(std::type_identity<T>)->
			ScatterMethodVariant_t {
			return T{};
		}, VariantFromString (
			[]<class T>(
				std::type_identity<T>, 
				std::string_view S
			) {
				return T::detect(S);
			}, methodS,
			std::type_identity<ScatterMethodVariant_t>{}
		)
	);
	auto mes_consr =
		[]<class T>(std::type_identity<T>)->
		ScatterMeasureVariant_t
	{
		return T{};
	};
	std::tuple < double , double > measure_tp = pyget<std::tuple<double, double>>(
		std::make_tuple((double)1,(double)2), ExtraArgs, "measure");

	
	std::visit([&]
		<class ThermGen_t,
			class Vt, class Bt, class Gvt, 
			evdm::GridEL_type G_tp
		>
		(Matrix_Pair_Inst<Vt, Bt, Gvt, G_tp> &m_matrix, 
			ThermGen_t m_therm_gen)->void {
		typedef Gvt T;


		evdm::Distribution<Vt, Bt, Gvt, G_tp>& m_evp = m_matrix.second;
		evdm::GridMatrix<Vt, Bt, Gvt, G_tp>& m_mat = m_matrix.first;
		evdm::xorshift<T> G(seed);
		
		try {
			Nmk = Nmk_v.cast<size_t>();
		}
		catch (pybind11::cast_error&) {}

		auto [p, q] = measure_tp;
		auto Emin = m_evp.grid().inner(0).grid()[0].left;
		evdm::measure_dEpdlq<Gvt> m_el_measure(p, q, Emin);


		pybind11::gil_scoped_release m_unlock;
		if (Nmk) {
			evdm::Scatter(
				m_mat, m_evp, count_evap, ptype_in, ptype_out,
				sc_event, G, m_therm_gen, m_el_measure,
				M_DM, deltaM, NucleiM, Nmk, Nmk_per_traj, weight, ProgFunc
			);
		} else {
			auto Nmk_vc = get_N_distrib_from_handle(m_evp.grid().inner(0), Nmk_v);
			Nmk = std::accumulate(Nmk_vc.begin(), Nmk_vc.end(),0) /
				 Nmk_vc.size();
			evdm::Scatter(
				m_mat, m_evp, count_evap, ptype_in, ptype_out,
				sc_event, G, m_therm_gen, m_el_measure,
				M_DM, deltaM, NucleiM, Nmk_vc, Nmk_per_traj, weight, ProgFunc
			);
		}
		
	}, ScatterAccum.m_matrix,ThermGenVar);


	auto sce_info = sc_event.unique ?
		scatter_event_info(sc_event.name, ptype_in, ptype_out, 0, Nmk) :
		scatter_event_info();

	{
		pybind11::gil_scoped_acquire m_lock;
		ScatterAccum.events.push_back(sce_info);
	}

}



void add_pyscatter_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;

	m.def("CalcScatterImpl", &Py_ScatterProcess,
		"Calculates scatter matrix part, add event to matrix class\n\n"
		"Parameters:\n"
		"___________\n"
		"matrix : Matrix\n\tScatter matrix histo.\n"
		"ptype_in : int\n\tindex of input particle type.\n"
		"ptype_out : int\n\tindex of output particle type.\n"
		"m_wimp : float\n\tdm particle mass, GeV.\n"
		"delta : float\n\tdelta mass GeV: output mass - input mass.\n"
		"m_nuc : float\n\tnuclear mass, GeV.\n"
		"sc_event : ScatterEvent.\n"
		"Nmk : int | function | vector\n\tnumber of monte-carle steps.\n\t"
		"May depends on (e,l) or on (e0,e1,l0,l1)\n"
		"measure: how E,L distributed in bin.\n\t"
		"should be a tuple (p,q). E and l would be uniformly distributed by dE^pdl^q. default (1,2)\n\t"
		"method : str\n\t method of generating therm velocity of nuclei."
		"can be: \n\t 'notherm', 'naive',"
		"'soft' (more probability of high velocities)"
		", 'soft_tresh' (same as soft, but considering inelastic treshold)\n"
		"Nmk_traj: int\n\tnumber of monte-carle steps on each trajectory.\n"
		"weight : float\n\t[optional] scale factor, default - 1.\n"
		"bar : object\n\t[optional] progress bar update function.",
		py::arg("matrix"),
		py::arg("ptype_in"),
		py::arg("ptype_out"),
		py::arg("m_wimp"),
		py::arg("delta"),
		py::arg("m_nuc"),
		py::arg("sc_event"),
		py::arg("Nmk")
	);
}

