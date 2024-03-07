#include "core_python.hpp"
#include <evdm/print_distrib.hpp>
#include "debugdef.hpp"
#include <evdm/measure.hpp>
#include "grob_python.hpp"
Py_EL_Grid Py_Distribution::getGrid()const {
	return std::visit([](const auto& distrib) {
		return Py_EL_Grid(distrib.Grid);
	}, m_distrib);
}

Py_Distribution Py_Distribution::copy() const
{
	return std::visit([](const auto& m_distr) {
		return Py_Distribution(m_distr);
	}, m_distrib);
}

std::string Py_Distribution::repr() const
{
	double n_particles = this->count(-1);
	return std::visit([n_particles](auto const& _distrib) {
		return std::string("Distrib( dtype = ") +
			type_name<decltype(_distrib.get_vtype())>() +
			", count = " + std::to_string(n_particles) + ")";
	},m_distrib);
}

Py_Distribution Py_Distribution::as_type(const char* type_n)const
{
	std::string_view _dtype = type_n;
	if (_dtype == "float")
#ifdef DISTRIB_USE_FLOAT
		return as_type_t<float>();
#else
		throw pybind11::type_error("unsupported distrib type 'float'");
#endif
	else if (_dtype == "double") {
#ifdef DISTRIB_USE_DOUBLE
		return as_type_t<double>();
#else
		throw pybind11::type_error("unsupported distrib type 'double'");
#endif
	}
	else {
		throw pybind11::type_error("wrong data type: " + std::string(type_n) + ", expect float or double");
	}
}

double Py_Distribution::count(int ptype)const
{

	return std::visit([ptype](const auto& distrib)->double {
		if (ptype  > (int)distrib.grid().grid().size()) {
			throw pybind11::index_error("ptype is more than ptypes number");
		}
		auto H1 = distrib.as_histo();
		DEBUG1(H1.Grid.size());
		DEBUG1(H1.Values.size());
		return distrib.count(ptype);
			
	},m_distrib);
}

pybind11::tuple Py_Distribution::get_r_dense(
	pybind11::handle ptypes,
	double r_min, double  r_max, size_t Nr,
	size_t Nperbin,
	pybind11::handle opt_dense_funct) const
{
	int Ntypes = getGrid().ptypes();
	auto check_ptype = [Ntypes](size_t p) {
		if (p < Ntypes)
			return p;
		else {
			throw pybind11::value_error("ptypes contain incorrect index");
		}
	};
	std::vector<size_t> m_ptypes;
	try {
		for (auto it : ptypes) {
			m_ptypes.push_back(check_ptype(it.cast<size_t>()));
		}
	}
	catch (...) {
		m_ptypes.clear();
		int p = ptypes.cast<int>();
		if (p == -1) {
			m_ptypes.resize(Ntypes);
			for (size_t i = 0; i < Ntypes; ++i)
				m_ptypes[i] = i;
		}
		else {
			m_ptypes.push_back(check_ptype(p));
		}
	}
	auto dense_function = [opt_dense_funct](double t_r) {
		return opt_dense_funct.call(t_r).cast<double>();
	};
	bool is_given_dense = !opt_dense_funct.is_none();

	grob::GridVector<double> rGrid =
		is_given_dense ?
		grob::density_vector_nt(r_min, r_max, dense_function, Nr) :
		grob::GridVector<double>(r_min, r_max, Nr);

	return  std::visit(
		[&m_ptypes,Nperbin,rgrid = std::move(rGrid)]
		(auto const& distrib)->pybind11::tuple {
			using Gen_t = decltype(distrib.get_grid_vtype());
			auto r_dence_func = distrib.r_dens(
				m_ptypes, evdm::xorshift<Gen_t>{},
				std::move(rgrid),
				Nperbin
			);
			return make_python_function_1D(r_dence_func);
		},
		m_distrib
	);
	
}

pybind11::tuple Py_Distribution::plot(size_t ptype) const
{
	return std::visit([ptype](const auto& _distrib) {
		if (ptype > _distrib.grid().grid().size()) {
			throw pybind11::index_error("ptype is more than ptypes number");
		}
		typedef decltype(_distrib.get_vtype()) T;
		typedef decltype(_distrib.get_grid_vtype()) GT;

		evdm::tri_sizes m_sizes = evdm::GridEl_TriSizes(*_distrib.Grid.Grid);

		using shape = pybind11::array::ShapeContainer;

		pybind11::array_t<GT> Vertexes(shape({ (int)m_sizes.pts, 2 }));
		pybind11::array_t<int> TrIndexes(shape({ (int)m_sizes.trs, 3 }));
		pybind11::array_t<T> Values(m_sizes.pts);

		evdm::DirstributionPrinting(
			_distrib.as_histo(), ptype, _distrib.Grid.LE(), 
			Values.mutable_data(), 
			Vertexes.mutable_data(), 
			TrIndexes.mutable_data(), m_sizes
		);
		return pybind11::make_tuple(Vertexes, TrIndexes, Values);
	}, m_distrib);
}

Py_Distribution CreatePyDistrib(
	Py_EL_Grid const& mGridEL,
	const char * dtype,
	pybind11::handle Init
) 
{
	auto Creator = [&](auto _type_marker) {
		typedef typename decltype(_type_marker)::type T;

		auto InitFunction = [Init](T e, T l) {
			return pybind11::cast<T>(Init(e, l));
		};
		if (Init.is_none()) {
			return Py_Distribution(mGridEL, _type_marker,false);
		}else {
			return Py_Distribution(mGridEL, _type_marker, InitFunction);
		}
	};
	std::string_view _dtype = dtype;
	if (_dtype == "float") {
#ifdef DISTRIB_USE_FLOAT
		return Creator(Py_EL_Grid::m_type_marker<float>{});
#else
		throw pybind11::type_error("unsupported distrib value type 'float'");
#endif
	} else if(_dtype == "double"){
#ifdef DISTRIB_USE_DOUBLE
		return Creator(Py_EL_Grid::m_type_marker<double>{});
#else
		throw pybind11::type_error("unsupported distrib value type 'double'");
#endif
	} else {
		throw pybind11::type_error("wrong data type: " + std::string(_dtype) + ", expect float or double");
	}
}

double compare_distribs(
	Py_Distribution const& D1,
	Py_Distribution const& D2,
	std::string const& measure,
	double p_deg,
	int ptype) 
{
	return std::visit(
		[p_deg, &measure, ptype]
	(const auto& D1, const auto& D2)->double {
			if (
				static_cast<void*>(D1.Grid.body.get()) != 
				static_cast<void*>(D2.Grid.body.get())
			) 
			{
				throw pybind11::value_error(
					"compared distributions should have grids of the same Body"
				);
			}
			if (measure == "dEdl") {
				return evdm::distrib_norm(D1, D2, p_deg, evdm::measure_dEdl{});
			}
			else if (measure == "dEdL") {
				return evdm::distrib_norm(D1, D2, p_deg, evdm::measure_dEdL{});
			}
			else if (measure == "dEdL2")
			{
				evdm::distrib_norm(D1, D2, p_deg, evdm::measure_dEdL2{});
			}
			else {
				throw pybind11::type_error(
					"undefined meusure type '" +
					measure + "', expect dEdl, dEdL, dEdL2");
			}
			return -1;
		}, D1.m_distrib,D2.m_distrib
	);
}
	

Py_DistribMeasure::Py_DistribMeasure(std::string measure, double p_deg, int ptype):
	measure(std::move(measure)), p_deg(p_deg), ptype(ptype){}

double Py_DistribMeasure::call(Py_Distribution const& D1, Py_Distribution const& D2) {
	return compare_distribs(D1, D2, measure, p_deg, ptype);
}

std::string Py_DistribMeasure::repr() const
{
	return "Metric( dmu = " + measure + 
		", p = " + std::to_string(p_deg) + 
		", ptype = " + 
		(ptype >= 0 ? std::to_string(ptype) : std::string("all"));
}

void Py_Distribution::add_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;
	using pyobj_ref = py::object const&;
	py::class_<Py_Distribution>(m, "Distrib")
		.def(py::init(&CreatePyDistrib),
			"constructor of Distrib class\n" 
			"ELGrid -- Created EL Grid\n"
			"dtype -- float or double\n"
			"init -- initialiser fuinction (density)"
			" of 3 args: (ptype,e,l),\n"
			"default -- None, meaning zero dirtribution",
			py::arg("ELGrid"),
			py::arg_v("dtype", "float"),
			py::arg_v("init" , py::none()))
		.def("plot",&Py_Distribution::plot,
			 "returns tuple: (vertexes,triangles,values)\n"
			"where vertexes -- (N,2) shape array of vert coords,\n"
			"triangles -- triangles (N,3) shape array of indexes\n"
			"values -- values corresponding to vertexes", py::arg("ptype"))
		.def("count", &Py_Distribution::count,
			 "calculates number of particles",
				py::arg_v("ptype", -1))
		.def("as_type", &Py_Distribution::as_type, py::arg("dtype"),
			"creating Distrib with another dtype")
		.def("copy", &Py_Distribution::copy)
		.def_property_readonly("grid",&Py_Distribution::getGrid)
		.def("rdens",&Py_Distribution::get_r_dense,
			"return density function, i.e. d^N/d^3r\n"
			"ptypes -- array or int: considered ptypes\n"
			"rmin, rmax -- min and max radius in r distribution\n"
			"Nr -- number of points in r grid from rmin to rmax\n"
			"(Nb) -- number of MK iterations per bin for integration\n"
			"(rden) -- optional density function of r points location",
			py::arg_v("ptypes",-1), 
			py::arg_v("rmin", 0), py::arg_v("rmax", 1), py::arg_v("Nr", 100),
			py::arg_v("Nb", 10000), py::arg_v("rden", py::none())
		);
	
	py::class_<Py_DistribMeasure>(m, "Metric")
		.def(
			py::init<std::string, double, int >(),
			"Creating Lp Metric of distributions\n"
			"using: M(D1,D2), constructor arguments:\n"
			"Mes --\n"
			"\t'dEdl': measure differential is dE and dl, where l = L/Lmax(E)\n"
			"\t'dEdL': measure differential is dE and dL\n"
			"\t'dEdL2': measure differential is dE and dL^2 = 2LdL\n"
			"p_deg -- value of p parameter of Lp norm\n"
			"ptype -- undex of particles to compare, if -1, than summ of all particles",
			py::arg_v("Mes", "dEdL"),
			py::arg_v("p_deg", 1),
			py::arg_v("ptype", -1)
		)
		.def("__call__", &Py_DistribMeasure::call,"compare two distribs");
}
