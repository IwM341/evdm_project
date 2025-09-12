#include "dynamic_python.hpp"
//#include <evdm/core/core_dynamics_scatter.hpp>
#include "scatter_impl/scatter_impl.hpp"
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
	std::string algolS = pyget<std::string>("naive", ExtraArgs, "algol");

	double zero_val = pyget<double>(0.0, ExtraArgs, "zero");

	std::variant<
		std::integral_constant<int,0>,
		std::integral_constant<int, 1>,
		std::integral_constant<int, 2>
	> algol = std::integral_constant<int, 0>{};
	if (algolS == "naive") {
		algol = std::integral_constant <int, evdm::AlgolNaive>{};
	}
	else if (algolS == "shift") {
		algol = std::integral_constant <int, evdm::AlgolShift>{};
	}
	else if (algolS == "diffuse") {
		algol = std::integral_constant <int, evdm::AlgolDiffuse>{};
	}
	else {
		std::ostringstream error;
		error << "algol should be 'naive', 'shift' or 'diffuse', but got " << algolS;
		throw std::runtime_error(error.str());
	}


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
	Nmk_vector Nmk_var;

	try {
		Nmk_var = Nmk_v.cast<size_t>();
	}
	catch (pybind11::cast_error&) {
		Nmk_var =
			std::visit([&](auto && m_matevap)->Nmk_vector {
				return get_N_distrib_from_handle(
					m_matevap.second.grid().inner(0), Nmk_v
				);
			}, ScatterAccum.m_matrix);
	}

	
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


	ScatterImpl_Args m_scatter_args{ ScatterAccum,
	 zero_val,ptype_in, ptype_out, seed,
	sc_event,  ThermGenVar,
	measure_tp, M_DM, deltaM, NucleiM,
	Nmk_var, Nmk_per_traj,
	weight,  ProgFunc };
	
	{
		pybind11::gil_scoped_release m_unlock;
		std::visit([&](auto m_algol) {
			m_scatter_args.apply(m_algol);
		}, algol);
	}

	auto sce_info = sc_event.unique ?
		scatter_event_info(sc_event.name, ptype_in, ptype_out, 0, Nmk) :
		scatter_event_info();

	{
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
		", 'soft_tresh' (same as soft, but considering inelastic treshold), full\n"
		"algol : str\n\t method of calc matrix: could be \n\t"
		"naive, shift, diffuse"
		"zero : float\n\t if matrix element Sij < zero, then it assumed to be zero"
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

