#include "core_python.hpp"
#include "debugdef.hpp"
#include "grob_python.hpp"
size_t Py_EL_Grid::size() const
{
	return std::visit([](auto const& ELgrid) {
			return ELgrid.Grid->size();
		}, m_grid);
}

size_t Py_EL_Grid::ptypes() const
{
	return std::visit([](auto const& ELgrid) {
			return ELgrid.Grid->ptypes();
		}, m_grid);
}

const char* Py_EL_Grid::dtype() const
{
	return std::visit([](auto const& ELgrid) {
			return type_name<decltype(ELgrid.get_grid_vtype())>();
		}, m_grid);
}

std::string Py_EL_Grid::repr() const
{
	return std::visit([](auto const& ELgrid) 
	{
		const char * _dt = type_name<decltype(ELgrid.get_grid_vtype())>();
		return "GridEL( ptypes = " +
			std::to_string(ELgrid.Grid->ptypes()) +
			", size = " +
			std::to_string(ELgrid.Grid->size()) +
			", dtype = " + _dt + ")";

	}, m_grid);
}

std::string Py_EL_Grid::printd1() const
{
	std::stringstream os;
	std::visit(
		[&os](const auto& grid_v) {
			os << *grid_v.Grid;
		},m_grid
	);
	return os.str();
}

Py_BodyModel Py_EL_Grid::getBody() const
{
	return std::visit([](auto& grid) {
			return Py_BodyModel(grid.body);
		}, m_grid);
}

pybind11::array Py_EL_Grid::getPlot(bool is_internal) const
{
	return std::visit(
		[is_internal](auto const& Grid) {
			auto lines_sizes = evdm::GridEL_printing_size(*Grid.Grid);
			int N = lines_sizes.size();
			using T = decltype(Grid.get_grid_vtype());
			
			//pybind11::array_t<T> X0(N), X1(N), Y0(N), Y1(N);

			pybind11::array_t<T> Points(pybind11::array::ShapeContainer({N,2,2 }));
			DEBUG1(Points.size());
			auto vector_view = [_data = Points.mutable_data(), N, &Points](size_t shift) {
				return grob::as_container([_data, shift,&Points](size_t i)->T& {
						if (4 * i + shift >= Points.size()) {
							std::string _error = "error, out f range 4 * i + shift >= Points.size(), where  = " +
								std::to_string(shift) + ", i = " + std::to_string(i) + ", size = " + std::to_string(Points.size());
							throw pybind11::index_error(_error.c_str());
						}
						return _data[4 * i + shift];
					}, N);
			};
			
			auto X0 = vector_view(0);
			auto X1 = vector_view(2);
			auto Y0 = vector_view(1);
			auto Y1 = vector_view(3);


			if (is_internal){
				evdm::GridEL_printing_fill(*Grid.Grid, [](auto const&_e) {return 1; },
					X0, X1,
					Y0, Y1,
					lines_sizes
				);
			} else {
				auto const& LE = Grid.LE();
				evdm::GridEL_printing_fill(*Grid.Grid, 
					[&LE](auto const& _e) {return LE(-_e); },
					X0, X1,
					Y0, Y1,
					lines_sizes
				); 
			}
			return pybind11::array(Points);
		},
		m_grid
	);
}

pybind11::tuple Py_EL_Grid::getLE() const
{
	return std::visit(
		[](auto const& EL_grid_w) {
			typedef decltype(EL_grid_w.get_bm_vtype()) T;
			size_t Ne = 2*(EL_grid_w.Grid->size_e()+1);
			DEBUG1(Ne);
			auto _grid_e = grob::GridUniform<T>(-EL_grid_w.body->Phi[0], 0, Ne);
			DEBUG1(_grid_e);
			DEBUG1(EL_grid_w.LE()(-_grid_e[0]));
			return make_python_function_1D(
				grob::make_function(
					_grid_e,
					grob::as_container([le = EL_grid_w.LE(), _grid_e](size_t i) {
						return le(-_grid_e[i]);
					}, Ne)
				)
			);
		},m_grid
	);
}



Py_EL_Grid CreateELGrid(
	Py_BodyModel const& BM,
	size_t ptypes, size_t Ne,
	pybind11::handle Nl_func_or_size,
	pybind11::handle RhoE_func_opt,
	pybind11::handle RhoL_func_opt,
	std::string const & dtype) 
{
	DEBUG1(ptypes);
	DEBUG1(Ne);
	auto ret_action = [&](auto marker)-> Py_EL_Grid {
		typedef typename decltype(marker)::type T;
		//DEBUG1(type_name<T>());
		return std::visit(
			[&](auto const& BM_p) {
				//DEBUG1(type_name<decltype(BM_p->get_vtype())>());
				auto RhoE_func = [&](T t_e)->T 
				{
					return pybind11::cast<T>(RhoE_func_opt(t_e));
				};
				auto RhoL_func = [&](T t_e,T t_l)->T 
				{
					return pybind11::cast<T>(RhoE_func_opt(t_e,t_l));
				};
				int nl_size = -1;
				try {
					nl_size = (int)pybind11::cast<size_t>(Nl_func_or_size);
				}
				catch (...) {}
				if (nl_size >= 0) 
				{
					DEBUG1("nl >=0");
					auto Nl_func = [nl_size](auto t_e) {return nl_size; };
					if (RhoE_func_opt.is_none()) 
					{
						DEBUG1("RhoE_func_opt is none");
						return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, 0, 0, marker, 
							std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>{});
					}
					else {
						if (RhoL_func_opt.is_none()) 
						{
							DEBUG1("RhoL_func_opt is none");
							auto _RhoL_func = [&](T t_e, T t_l)->T {
								return 1;
							};
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, _RhoL_func, marker,
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{});
						}
						else {
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, RhoL_func, marker, 
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{});
						}
					}
				}
				else {
					DEBUG1("nl determined by function");
					auto Nl_func = [nl_size,&Nl_func_or_size](auto t_e)->size_t {
						return pybind11::cast<double>(Nl_func_or_size(t_e)); 
					};
					if (RhoE_func_opt.is_none()) 
					{
						DEBUG1("RhoE_func_opt is none");
						return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, 0, 0, marker, 
							std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>{});
					}
					else {
						if (RhoL_func_opt.is_none()) {
							DEBUG1("RhoL_func_opt is none");
							auto _RhoL_func = [&](T t_e, T t_l)->T {
								return 1;
							};
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, _RhoL_func, marker, 
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{});
						}
						else {
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, RhoL_func, marker, 
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{});
						}
					}
				}
			}, BM.m_body);
	};
	DEBUG1(__PRETTY_FUNCTION__);
	if (dtype == "float") {
		DEBUG1("float");
		return ret_action(Py_EL_Grid::m_type_marker<float>{});
	}
	else if (dtype == "double") {
		DEBUG1("double");
		return ret_action(Py_EL_Grid::m_type_marker<double>{});
	} else {
		throw pybind11::type_error("wrong data type: " + dtype + ", expect float or double");
	}

}

void Py_EL_Grid::add_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;
	using pyobj_ref = py::object const&;
	py::class_<Py_EL_Grid>(m, "GridEL")
		.def(py::init(&CreateELGrid),
			"constructs eELGrid\n"
			"body -- instance of Body class\n"
			"ptypes -- number of particle types\n"
			"Ne -- number of bins of E axis\n"
			"Nl_func -- if number then number of bins of L axis,"
			"else a function from [0.0,1.0]->int,"
			"indicating number of bins in L grid depending on e,"
			"where 0.0 correspond to Emin = phi(0), and 1.0 -- to Emax = 0\n"
			"RhoE [optional] -- bin density of e axis RhoE : e in [0.0,1.0]->float\n"
			"RhoE [optional] -- bin density of l axis RhoL : (e,l) in [0.0,1.0]x[0.0,1.0]->float\n"
			"dtype -- float or double",
			py::arg("body"),
			py::arg("ptypes"),
			py::arg("Ne"),
			py::arg("Nl_func"),
			py::arg_v("RhoE", py::none()),
			py::arg_v("RhoL", py::none()),
			py::arg_v("dtype", "float")
		)
		.def_property_readonly("size",&Py_EL_Grid::size,"number of bins in grid")
		.def_property_readonly("dtype",&Py_EL_Grid::dtype,"type of grid")
		.def_property_readonly("ptypes",&Py_EL_Grid::ptypes,"number of particles")
		.def("__repr__",&Py_EL_Grid::repr)
		.def_property_readonly("body", &Py_EL_Grid::getBody)
		.def("plot",&Py_EL_Grid::getPlot,
			"return array of shape (N,2,2):\n"
			"[ [[x_start_i,y_start_i],[x_end_i,y_end_i]],...]\n"
			"arrays could be plottted with matplotlib:\n"
			"lc = matplotlib.collections.LineCollection(result of this function)\n"
			"fig, ax = plt.subplots()\n"	
			"ax.add_collection(lc)\n"
			"if is_internal true, then l from 0 to 1, else from 0, l(e)", py::arg_v("is_internal",false)
		).def("LE", &Py_EL_Grid::getLE,
			"return tuple (E array,L(E) array)"
		).def("print_debug",&Py_EL_Grid::printd1);
}