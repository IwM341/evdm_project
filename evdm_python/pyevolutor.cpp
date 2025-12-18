#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <evdm/utils/variant_tools.hpp>
#include "grob_python.hpp"
#include <evdm/utils/prng.hpp>

#include <evdm/core/core_evolutor.hpp>
#include "dynamic_python.hpp"
#include "distrib_python.hpp"
#include "scatter_impl/dist_sampler.hpp"
#include <format>
struct Py_FFC_Impl {
	pybind11::object impl;

	static pybind11::object __getstate__(Py_FFC_Impl const & _this) {
		return _this.impl;
	}
	Py_FFC_Impl(pybind11::object impl) :impl(impl) {}

	static Py_FFC_Impl __setstate__(pybind11::object _impl) {
		return Py_FFC_Impl{ _impl };
	}
	std::vector<evdm::ElementInfo_t>  construct(size_t pin, size_t pout) const{
		pybind11::list m_infos = impl.attr("construct").call(pin, pout);
		std::vector<evdm::ElementInfo_t> m_infos_ret;
		m_infos_ret.reserve(m_infos.size());

		for(auto it : m_infos) 
		{
			pybind11::dict m_info = it.cast<pybind11::dict>();
			size_t A = m_info["A"].cast<size_t>();
			size_t Z = m_info["Z"].cast<size_t>();
			float mN = m_info["mN"].cast<float>();
			float mX = m_info["mX"].cast<float>();
			float dmX = m_info["delta"].cast<float>();
			Py_ScatterEvent info = m_info["event"].cast<Py_ScatterEvent>();
			m_infos_ret.push_back( evdm::ElementInfo_t{
			   .A = A,.Z = Z,.mN = mN,.mX = mX,.dmX = dmX,
			   .ff = std::move(info) }
			);
		}
		return m_infos_ret;
	}
};

struct Py_Evolutor {
	template <typename T>
	using MCE_t = evdm::MCEvolutor<T, Py_FFC_Impl>;
	std::variant<
		MCE_t<float>, MCE_t<double>
	> m_evolutor;

	size_t ptypes()const {
		return std::visit([]<class T>(MCE_t<T> const& ev) {
			return ev.ptypes();
		}, m_evolutor);
	}
	Py_Evolutor(std::variant<
		MCE_t<float>, MCE_t<double>
	> m_evolutor) : m_evolutor(std::move(m_evolutor)) {}
	Py_Evolutor(Py_BodyModel B, size_t ptypes,pybind11::object ffc_impl): 
		m_evolutor(std::visit([&]<class T>(evdm::BodyModel<T> const& B)->decltype(m_evolutor){
		return evdm::MCEvolutor<T, Py_FFC_Impl>(*B, ptypes, Py_FFC_Impl{ ffc_impl });
	}, B.m_body)) {}

	pybind11::array evolute(
		pybind11::array States,  double evolve_time, 
		pybind11::dict m_probs, size_t max_scatter, size_t seed,std::string errors) const
	{
		if (!seed) {
			seed = 1;
		}
		std::list<std::string> m_errors;
		{
			std::replace(errors.begin(), errors.end(), ';', '|');
			std::replace(errors.begin(), errors.end(), ',', '|');
			std::istringstream f(errors);
			std::string s;
			
			while (std::getline(f, s, '|')) {
				m_errors.push_back(s);
			}
		}
		return std::visit([&]<class T>(MCE_t<T> const& ev)->pybind11::array {
			std::array<T, 16> m_probs_arr;
			size_t pt = ptypes();
			for (auto& a : m_probs_arr) a = 1;
			for (auto m_info : m_probs) {
				auto indexes = m_info.first.cast<std::pair<size_t, size_t> >();
				T value = evdm::downbound(m_info.second.cast<T>(), (T)0);
				if (indexes.first < pt && indexes.second < pt) {
					m_probs_arr[indexes.first * pt + indexes.second] = value;
				}
				else {
					throw std::out_of_range("evolute: probs factors are out of range ptype");
				}
			}
			pybind11::array_t<evdm::StateEL<T>> m_states = States;
			ev.evolute(
				std::span < evdm::StateEL<T>>(m_states.mutable_data(), m_states.size()),
				evdm::xorshift<T>(seed), m_probs_arr, evolve_time, max_scatter, m_errors
			);
			return m_states;
		}, m_evolutor);
	}
	void addScatterProcess(size_t ptypein, size_t ptypeout,pybind11::kwargs params) {
		std::visit([&]<class T>(MCE_t<T> & ev) {
			size_t seed = pyget<size_t>(123234324, params, "seed");
			size_t initial_rv_level = pyget<size_t>(6, params, "rv_linit");
			size_t max_rv_level = pyget<size_t>(16, params, "rv_lmax");
			size_t max_rv_size = pyget<size_t>(10000, params, "rv_size");
			T rv_cmp_accept = pyget<T>(0.05, params, "rv_acc");
			size_t initial_el_level = pyget<size_t>(6, params, "el_linit");
			size_t max_el_level = pyget<size_t>(16, params, "el_lmax");
			size_t max_el_size = pyget<size_t>(10000, params, "el_size"); 
			T el_cmp_accept = pyget<T>(0.05, params, "el_acc");
			size_t ThetaMaxSteps = pyget<size_t>(100, params, "max_theta_steps");
			T maxProbTheta =  pyget<T>(0.05, params, "max_theta_prob"); 
			T zeroProb = pyget<T>(1e-16, params, "zero_prob");
			size_t Nmk = pyget<size_t>(200, params, "Nmk");
			ev.AddProcess(ptypein, ptypeout,
				initial_rv_level, max_rv_level, max_rv_size, rv_cmp_accept, 
				initial_el_level, max_el_level, max_el_size, el_cmp_accept,
				ThetaMaxSteps, maxProbTheta, zeroProb,evdm::xorshift<T>(seed),Nmk
			);

		}, m_evolutor);
	}
	pybind11::tuple plot_info(size_t pin,size_t pout,std::string what_to_plot, pybind11::kwargs extra){
		return std::visit(
			[&]<class T>(MCE_t<T>&ev)->pybind11::tuple {
			evdm::PlotInfo<T> m_plt;
			if (what_to_plot == "el") {
				m_plt = ev.plotinfo_el(pin, pout, extra["param"].cast<std::string_view>());
			} else if (what_to_plot == "rv") {
				m_plt = ev.plotinfo_rv(pin, pout, extra["A"].cast<size_t>(), extra["Z"].cast<size_t>());
			}
			else {
				throw std::runtime_error(std::format("unexpected key {}", what_to_plot));
			}
			
			using shape = pybind11::array::ShapeContainer;
			pybind11::array_t<T> XArray(shape({ (int)m_plt.X.size() }));
			pybind11::array_t<T> YArray(shape({ (int)m_plt.Y.size() }));

			
			pybind11::array_t<int> TrIndexes(shape({ (int)m_plt.Indicies.size() / 3, 3 }));
			pybind11::array_t<T> Values(m_plt.Values.size());
			std::copy(m_plt.X.begin(), m_plt.X.end(), XArray.mutable_data());
			std::copy(m_plt.Y.begin(), m_plt.Y.end(), YArray.mutable_data());
			std::copy(m_plt.Indicies.begin(), m_plt.Indicies.end(), TrIndexes.mutable_data());

			std::copy(m_plt.Values.begin(), m_plt.Values.end(), Values.mutable_data());
			return pybind11::make_tuple(XArray, YArray, TrIndexes, Values);

		}, m_evolutor);
	}



	pybind11::array initial_states(
		Py_Distribution const &init_distrib,
		size_t N, pybind11::tuple m_mes, size_t seed)
	{
		return std::visit([&]<class T>(MCE_t<T> const& ev)->pybind11::array
		{
			return init_samples_impl(
				std::type_identity<T>{}, init_distrib, N, seed , 
				m_mes[0].cast<double>(), m_mes[1].cast<double>()
			);
		},m_evolutor);
	}
	Py_Distribution to_distrib(Py_EL_Grid grid, pybind11::array states,double total_count) {
		return std::visit([&]<class T>(MCE_t<T> const& ev) ->
			Py_Distribution 
		{
			return put_points_impl(std::type_identity<T>{}, grid, states, total_count);

		}, m_evolutor);
	}
	static pybind11::dict __getstate__(Py_Evolutor const& _this) {
		using namespace pybind11::literals;
		return std::visit([&]<class T>(MCE_t<T> const& ev) ->pybind11::dict {

			NpDictSerializator  S;	
			return pybind11::dict(
				"type"_a = type_name<T>(),
				"evolutor"_a = grob::Serialize(ev, S)
			);
		}, _this.m_evolutor);
	}

	static Py_Evolutor  __setstate__(pybind11::dict m_val) {

		auto type_var_t = type_from_str(
			m_val["type"].cast<std::string_view>(),
			std::type_identity<std::tuple<float, double>>{}
		);
		return std::visit(
			[&]<class T>(std::type_identity<T>)-> Py_Evolutor {
			NpDictSerializator  S;
			return Py_Evolutor{ grob::DeSerialize<MCE_t<T>>(m_val["evolutor"],S) };
		},type_var_t);
	}
};
void PyEvolutor_add_to_python_module(pybind11::module& m) {
	namespace py = pybind11;

	auto add_State = [&m]<class T>(std::type_identity<T>, const char* name) {
		using mst_t = evdm::StateEL<T>;
		mst_t (*init)(size_t,T,T) =  [](size_t pt, T e, T l) {return mst_t{ pt,e,l }; };
		py::class_<mst_t>(m, name)
			.def(py::init<size_t, T, T>())
			.def_readwrite("e", &mst_t::e)
			.def_readwrite("l", &mst_t::l)
			.def_readwrite("ptype", &mst_t::ptype);
		PYBIND11_NUMPY_DTYPE(evdm::StateEL<T>, ptype,e, l);
	};
	add_State(std::type_identity<float>{}, "Statef");
	add_State(std::type_identity<double>{}, "Stated");
	
	py::class_<Py_FFC_Impl>(m, "FF_Provider_Impl")
		.def(py::init<pybind11::object>())
		.def(py::pickle(&Py_FFC_Impl::__getstate__, &Py_FFC_Impl::__setstate__));

	py::class_<Py_Evolutor>(m, "MCEvolutor")
		.def(py::init<Py_BodyModel, size_t, pybind11::object>(),
			"evolutor class\n\n"
			"Parameters:\n"
			"___________\n"
			"body : evdm.Body\n"
			"ptypes : int\n"
			"ffprovider : evdm.FF_Provider",
			py::arg("body"), py::arg("ptypes"), py::arg("ffprovider")
		)
		.def(py::pickle(&Py_Evolutor::__getstate__, &Py_Evolutor::__setstate__))
		.def("addScatterProcess", &Py_Evolutor::addScatterProcess,
			"add all scatter elements for process pin->pout\n\n"
			"Parameters:\n"
			"___________\n"
			"pin : int\n"
			"pout : \n"
			"rv_linit\n"
			"rv_lmax\n"
			"rv_size\n"
			"rv_acc\n"
			"el_linit\n"
			"el_lmax\n"
			"el_size\n"
			"el_acc\n"
			"max_theta_steps\n"
			"max_theta_prob\n"
			"zero_prob\n"
			"Nmk\n"
			"seed\n",
			py::arg("pin"), py::arg("pout"))
		.def("gen_initial",
			&Py_Evolutor::initial_states,
			"generates initial states according to distribution\n\n"
			"Parameters:\n"
			"___________\n"
			"distrib : evdm.Distrib\n\t initial distribution\n"
			"N : int\n\tnumber of states\n"
			"mes : Tuple[float,float]\n\tbin measure\n"
			"seed : int",
			py::arg("distrib"), py::arg("N"), 
			py::arg_v("mes",py::make_tuple(1,2)), py::arg_v("seed",123))
		.def("evolute", &Py_Evolutor::evolute,
			"calculate evolution of initial states\n"
			"Parameters:\n"
			"___________\n"
			"states : array\n\t initial state to evolute\n"
			"probs : dict [optional]\n\t "
			"time : float\n\t time of evolution\n"
			"dict of additional mult factor, for example:\n\t"
			"{(0,0):0,(0,1):2} (by default they are 1)\n"
			"max_scat : int\n\t max number of scatter in evolution (to avoid infinit loop)\n"
			"seed : int",
			py::arg("states"), py::arg("time"),
			py::arg_v("probs",py::dict()), 
			py::arg_v("max_scat", 1000000000), py::arg_v("seed", 123), 
			py::arg_v("errors", ""))
		.def("to_distrib", &Py_Evolutor::to_distrib,
			"make grid distribution from states\n"
			"Parameters:\n"
			"___________\n"
			"grid : evd.GridEL\n\t grid where to put points\n"
			"states : array\n"
			"count: float\n\ttotal amount of particles ()",
			py::arg("grid"),
			py::arg("states"),
			py::arg("count"))
		.def("plot_info",&Py_Evolutor::plot_info,
			"plot probabilities\n"
			"Parameters:\n"
			"___________\n"
			"pin, pout \n\t initial and final ptype\n"
			"data: str\n\t 'rv' or 'el' (what to plot)\n\t"
			"if plot 'rv', need to point A and Z of element\n",
			py::arg("pin"), py::arg("pout"), py::arg("data")
		);
}