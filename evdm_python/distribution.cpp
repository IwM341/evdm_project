#include "distrib_python.hpp"
#include <evdm/print_distrib.hpp>
#include "debugdef.hpp"
#include <evdm/measure.hpp>
#include "grob_python.hpp"
#include <pybind11/functional.h>
#include "progress_log.hpp"
#include <pybind11/stl.h>
//#include <format>


Py_Distribution Py_Distribution::CreatePyDistribFromArray(
	Py_EL_Grid const& mGridEL, 
	pybind11::array values
) {
	auto array_var = array_variant<DISTRIB_TYPE_LIST>(values);
	return std::visit([&]<class T>( pybind11::array_t<T> const& m_array) {
		size_t _stride = m_array.strides()[0] / sizeof(T);
		size_t extra_size = m_array.size() - mGridEL.size()*mGridEL.ptypes();
		size_t padding = evdm::min_deg_2(extra_size);
		return Py_Distribution(mGridEL, m_array.data(), _stride, m_array.size());
	}, array_var);
}

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

pybind11::dict Py_Distribution::get_object(pybind11::handle self)
{
	using namespace pybind11::literals;
	return pybind11::dict(
		"type"_a = "evdm.Distribution",
		"values"_a = get_array(
				self, -1, true
		),
		"grid"_a = getGrid(),
		"padding"_a = get_padding()
	);
}

double Py_Distribution::count(int ptype)const
{

	return std::visit([ptype](const auto& distrib)->double {
		if (ptype  > (int)distrib.grid().grid().size()) {
			throw pybind11::index_error("ptype is more than ptypes number");
		}
		auto H1 = distrib.as_histo();
		VERBOSE_VAR1(H1.Grid.size());
		VERBOSE_VAR1(H1.Values.size());
		return distrib.count(ptype);
			
	},m_distrib);
}

pybind11::tuple Py_Distribution::get_r_dense(
	pybind11::handle ptypes,
	double r_min, double  r_max, size_t Nr,
	pybind11::handle Nperbin,
	pybind11::handle opt_dense_funct,
	pybind11::handle update_function) const
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
		return opt_dense_funct.operator()(t_r).cast<double>();
	};
	bool is_given_dense = !opt_dense_funct.is_none();

	grob::GridVector<double> rGrid =
		is_given_dense ?
		grob::density_vector_nt(r_min, r_max, dense_function, Nr) :
		grob::GridVector<double>(r_min, r_max, Nr);

	evdm::progress_omp_function ProgF = 
		make_progress_func(update_function);
	
	return  std::visit(
		[&m_ptypes,Nperbin,rgrid = std::move(rGrid), &ProgF]<_DISTRIB_TMPL_>(
			evdm::Distribution<_DISTRIB_PARS_> const& distrib
		)->pybind11::tuple 
		{
			using Gen_t = decltype(distrib.get_grid_vtype());
			
			int NperBin_i = -1;
			try {
				NperBin_i = Nperbin.cast<int>();
			}catch(pybind11::cast_error &){
				NperBin_i = -1;
			}
			if (NperBin_i > 0) {
				auto r_dence_func = distrib.r_dens(
					m_ptypes, evdm::xorshift<Gen_t>{},
					std::move(rgrid),
					NperBin_i, ProgF
				);
				return make_python_function_1D(r_dence_func);
			}
			else {
				auto Nperbin_v = 
					get_N_distrib_from_handle(
						distrib.grid().inner(0), Nperbin
					);
				return make_python_function_1D(
					distrib.r_dens(
						m_ptypes, evdm::xorshift<Gen_t>{},
						std::move(rgrid),
						Nperbin_v, ProgF
					)
				);
			}
			
			
		},
		m_distrib
	);
	
}

pybind11::tuple Py_Distribution::plot1o(
	size_t ptype,
	std::string_view m_measure,
	std::string_view LE_Space) const
{
	return std::visit([ptype, m_measure, LE_Space]<_DISTRIB_TMPL_>(
	const evdm::Distribution<_DISTRIB_PARS_>& _distrib) {

	if (ptype > _distrib.grid().grid().size()) {
		throw pybind11::index_error("ptype is more than ptypes number");
	}
	typedef decltype(_distrib.get_vtype()) T;
	typedef decltype(_distrib.get_grid_vtype()) GT;

	evdm::tri_sizes m_sizes = evdm::GridEl_TriSizes_1order(*_distrib.Grid.Grid);

	using shape = pybind11::array::ShapeContainer;

	pybind11::array_t<GT> XArray(shape({ (int)m_sizes.pts }));
	pybind11::array_t<GT> YArray(shape({ (int)m_sizes.pts }));

	pybind11::array_t<int> TrIndexes(shape({ (int)m_sizes.trs, 3 }));
	pybind11::array_t<T> Values(m_sizes.pts);

	size_t measure_var = 0;
	if (m_measure == "1" || m_measure == "" || m_measure == "") {
		measure_var = 0;
	} else if (m_measure == "dEdl" || m_measure == "El" || m_measure == "E-l") {
		measure_var = 1;
	}else if (m_measure == "dEdL" || m_measure == "EL" || m_measure == "default") {
		measure_var = 2;
	} else if (m_measure == "dEdL2" || m_measure == "dEdL^2" || m_measure == "EL2") {
		measure_var = 3;
	} else {
		using namespace std::string_literals;
		throw pybind11::value_error(
			"unexpected measure type '" + std::string(m_measure) + "', "
				"expect '', 'dEdl', 'dEdL' or 'dEdL2'"
		);
	}
	auto m_mes_variant = evdm::make_variant_alt(
		measure_var,
		evdm::measure_1{}, evdm::measure_dEdl{}, evdm::measure_dEdL{}, evdm::measure_dEdL2{}
	);
	auto m_histo_full = _distrib.as_histo();
	auto m_histo_slice = m_histo_full.inner_slice(ptype);

	size_t le_space_var_num = 0;
	if (LE_Space == "latent" || LE_Space == "inner") {
		le_space_var_num = 1;
	}
	auto LE_func_variant = evdm::make_variant_alt(
		le_space_var_num,
		_distrib.Grid.LE(),
		[](auto x) {return 1; }
	);

	std::visit([&](auto m_variant,auto && le_func_plot) {
		evdm::DirstributionPrinting_1order(
			m_histo_slice, _distrib.Grid.LE(), le_func_plot,
			Values.mutable_data(),
			XArray.mutable_data(),
			YArray.mutable_data(),
			TrIndexes.mutable_data(), m_sizes,
			m_variant
		);
	}, m_mes_variant, LE_func_variant);
	
	return pybind11::make_tuple(XArray, YArray, TrIndexes, Values);
	
	}, m_distrib);
}

pybind11::tuple Py_Distribution::plot2o(size_t ptype) const
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

pybind11::array Py_Distribution::get_array(pybind11::handle  self, int ptype,bool raw) 
{
	return std::visit(
		[ptype, self, raw]<_DISTRIB_TMPL_>(
			evdm::Distribution<_DISTRIB_PARS_> &distr
		)->pybind11::array {
		using T = decltype(distr.get_vtype());
			if (!raw) {
			Eigen::VectorBlock<Eigen::VectorX<T>> m_block = distr.block(ptype);
			auto shape = pybind11::array::ShapeContainer({ m_block.size() });
			auto strides = pybind11::array::ShapeContainer({ sizeof(T) });
			return pybind11::array_t<T>(
				pybind11::buffer_info(m_block.data(), shape, strides), self
				);
		}
		else {
			Eigen::VectorX<T> & m_block = distr.raw_vector();
			auto shape = pybind11::array::ShapeContainer({ m_block.size() });
			auto strides = pybind11::array::ShapeContainer({ sizeof(T) });
			return pybind11::array_t<T>(
				pybind11::buffer_info(m_block.data(), shape, strides), self
				);
		}
		
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
		throw pybind11::type_error(
			"wrong data type: " + 
			std::string(_dtype) + 
			", expect float or double"
		);
	}
}

Py_Distribution CreateDistribFromDict(
	Py_EL_Grid const& mGridEL,pybind11::dict m_dict) 
{
	pybind11::array X = m_dict["values"].cast<pybind11::array>();
	size_t padding = m_dict["padding"].cast<size_t>();

	if (X.ndim() != 1) {
		throw pybind11::value_error("except one dimentional array");
	}
	
	size_t m_size = X.size();
	auto array_vars = array_variant<DISTRIB_TYPE_LIST>(X);
	return std::visit([&]<class T>(pybind11::array_t<T> const&m_array) {
		size_t stride = m_array.strides(0)/sizeof(T);
		return Py_Distribution(
			mGridEL, (const T*)X.data(), stride, m_size
		);
	}, array_vars);
	
}
Py_Distribution Py_Distribution::from_dict(pybind11::dict const& distr_dict) {
	pybind11::handle G = distr_dict["grid"];
	if (pybind11::isinstance<pybind11::dict>(G)) {
		return CreateDistribFromDict(
			Py_EL_Grid::from_dict1(
				G.cast<pybind11::dict>()
			), distr_dict
		);
	}
	else {
		return CreateDistribFromDict(
			G.cast<Py_EL_Grid>(), distr_dict
		);
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
				pybind11::print(
					"Warning: compared distributions should have grids of the same Body"
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
				return evdm::distrib_norm(D1, D2, p_deg, evdm::measure_dEdL2{});
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
			"constructor of Distrib class\n\n"
			"Parameters:\n"
			"___________\n"
			"ELGrid : GridEL\n\tgrid, where distribution is created.\n"
			"dtype : string\n\t'float' or 'double'.\n"
			"init : lambda\n\tinitialiser fuinction (density)"
			"of 3 args: (ptype,e,l), "
			"default - None, meaning zero dirtribution.",
			py::arg("ELGrid"),
			py::arg_v("dtype", "float"),
			py::arg_v("init", py::none()))
		.def(py::init(&CreateDistribFromDict),
			"constructor of Distrib class\n\n"
			"Parameters:\n"
			"___________\n"
			"ELGrid : GridEL\n\tgrid, where distribution is created.\n"
			"object : dict representation of distrib",
			py::arg("ELGrid"),
			py::arg("object")
		)
		.def(py::init(&Py_Distribution::from_dict))
		.def("__getstate__", [](py::handle self)->py::dict {
				return py::cast<Py_Distribution&>(self).get_object(self);
			}
		)
		.def("__setstate__", [](py::handle self, py::dict state) {
				self.attr("__init__")(state);
			}
		)
		.def(py::init(&Py_Distribution::CreatePyDistribFromArray))
		.def("plot", &Py_Distribution::plot1o,
			"returns tuple: (X,Y,triangles,values)\n"
			"where X,Y - arrays of vertises x and y coords,\n"
			"triangles - triangles (N,3) shape array of indexes\n"
			"values - values corresponding to vertexes.\n\n"
			"Parameters:\n"
			"___________\n"
			"ptype : int\n\twimp type.\n"
			"mes : str\n\t measure of bins "
			"(standart is dEdL, no measure - '1' or '')\n"
			"space: str\n\t"
			"if space = 'latent', L would be from 0 to 1, "
			"otherwise - from 0 to Lmax(E)\n\n"
			"to plot using matplotlib type:\n\t"
			"import matplotlib.tri as mtri\n\t"
			"plt.tricontourf(triang,Z)\n\t"
			"plt.triplot(triang,color = 'black')",
			py::arg("ptype"), 
			py::arg_v("mes","dEdL"),
			py::arg_v("space",""))
		.def("plot2", &Py_Distribution::plot2o,
			"returns tuple: (vertexes,triangles,values)\n"
			"where vertexes - (N,2) shape array of vert coords,\n"
			"triangles - triangles (N,3) shape array of indexes\n"
			"values - values corresponding to vertexes.\n\t"
			"Parameters:\n"
			"___________\n"
			"ptype : int\n\twimp type.\n\n"
			"to plot using plotly type:\n\n"
			"fig = figf.create_trisurf(vert[:,0],vert[:,1], np.log(np.abs(-vals)),\n\t"
                "colormap=\"Portland\",\n\t"
                "simplices=trs,\n\t"
                "title=\"title\")",
			py::arg("ptype"))
		.def("count", &Py_Distribution::count,
			"calculates number of particles",
			py::arg_v("ptype", -1))
		.def("as_type", &Py_Distribution::as_type, py::arg("dtype"),
			"creating Distrib with another dtype")
		.def("copy", &Py_Distribution::copy)
		.def_property_readonly("grid", &Py_Distribution::getGrid)
		.def(
			"avarage",
			&Py_Distribution::get_avarage, 
			"calc avarage over distribution of F(e,l), "
			"where l \\in [0,1]", 
			py::arg("F"),
			py::arg_v("ptype", -1)
		).def(
			"Edistrib",
			&Py_Distribution::get_E_distrib,
			"reduce distribution to energy distribution",
			py::arg_v("ptype",-1)
		)
		.def("rdens", &Py_Distribution::get_r_dense,
			"git density function, i.e. d^N/d^3r\n\n"
			"Parameters:\n"
			"___________\n"
			"ptypes : array or int\n\t considered ptypes.\n"
			"rmin : float\n\tmin radius in r distribution.\n"
			"rmax : float\n\tmax radius in r distribution.\n"
			"Nr : int\n\tnumber of points in r grid from rmin to rmax.\n"
			"Nb : int\n\tnumber of MK iterations per bin for integration, default - 10000.\n"
			"rden : function\n\toptional density function of r points location.\n"
			"bar : progressbar\n\toptional progress bar update function.",
			py::arg_v("ptypes", -1),
			py::arg_v("rmin", 0), py::arg_v("rmax", 1), py::arg_v("Nr", 100),
			py::arg_v("Nb", 10000), py::arg_v("rden", py::none()),
			py::arg_v("bar", py::none())
		)
		.def(
			"to_numpy", 
			[](py::handle self, int ptype,bool is_raw)->py::array {
				return py::cast<Py_Distribution&>(self).get_array(
					self, ptype,is_raw
				);
			},
			"gives numpy array view to distribution\n\n"
			"Parameters:\n"
			"___________\n"
			"ptype : int\n\twimp type, default -1, meaning all.\n"
			"is_raw : bool\n\tgives distribution array with padding.",
			py::arg_v("ptype",-1), py::arg_v("is_raw", false)
		)
		.def_property_readonly(
			"grid",&Py_Distribution::getGrid,"get grid of distrib")
		.def("to_object", [](py::handle self)->py::dict {
			return py::cast<Py_Distribution&>(self).get_object(self);
		},"return numpy array, so Distrib/Capture could be restored"
			"from grid and this object:\n"
			"\tm_array = m_distrib.to_object()\n"
			"\tm_distrib1 = Distrib(m_distrib.grid(),m_array)\n");
	
	py::class_<Py_DistribMeasure>(m, "Metric")
		.def(
			py::init<std::string, double, int >(),
			"Creation Lp Metric of distributions\n"
			"use: M(D1,D2), constructor arguments:\n\n"
			"Parameters:\n"
			"___________\n"
			"Mes : Measure\n"
			"\t'dEdl' - measure differential is dE and dl, where l = L/Lmax(E)\n"
			"\t'dEdL' - measure differential is dE and dL\n"
			"\t'dEdL2' - measure differential is dE and dL^2 = 2LdL.\n"
			"p_deg : float\n\tvalue of p parameter of Lp norm.\n"
			"ptype : int\n\tindex of particles to compare, if -1, than summ of all particles.",
			py::arg_v("Mes", "dEdL"),
			py::arg_v("p_deg", 1),
			py::arg_v("ptype", -1)
		)
		.def("__call__", &Py_DistribMeasure::call,"compare two distribs");
}

size_t Py_Distribution::get_padding() const
{
	return std::visit([]<_DISTRIB_TMPL_>(
		evdm::Distribution<_DISTRIB_PARS_> const& dstrb) {
			return dstrb.get_padding();
		},
		m_distrib
	);
}

pybind11::tuple Py_Distribution::get_E_distrib(int ptype) const {

	return std::visit(
		[ptype]<_DISTRIB_TMPL_>(
			evdm::Distribution<_DISTRIB_PARS_> const& distrib
		)
	{
		auto Edistrib = evdm::get_E_distrib(distrib, ptype);
		return pybind11::make_tuple(
			make_py_array(distrib.grid().inner(0).grid().unhisto()),
			make_py_array(Edistrib)
		);
	}, m_distrib
	);
}

double Py_Distribution::get_avarage(
	pybind11::handle Functor, int ptype
) const {
	typedef std::function<double(double, double)> FuncType;
	return std::visit(
		[F = Functor.cast<FuncType>(), ptype]
		<_DISTRIB_TMPL_>(
			evdm::Distribution<_DISTRIB_PARS_> const& distrib
		) ->double
		{
			return evdm::avarage_by_grid(distrib, ptype,F);
		},m_distrib
	);
}