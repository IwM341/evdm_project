#include "core_python.hpp"
#include "debugdef.hpp"
#include "grob_python.hpp"
#include <pybind11/functional.h>
#include <set>
#include <pybind11/stl.h>

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

pybind11::dict Py_EL_Grid::_get_traj_all(double e, double l) {
	return std::visit(
		[e, l]<_GRID_EL_TMPL_>(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr)
	{
		evdm::Body<_T1> const& B = *m_gr.body;
		typename evdm::Body<_T1>::LE_func_t const & LE_s = m_gr.getBodyLE();
		auto L = l * LE_s.i_lm(-e);
		auto [r0, r1] = B.find_rmin_rmax(-e, L * L, LE_s);
		auto [theta_max, tau_max, _] = B.get_internal_traj(r0, r1, 200);
		using namespace pybind11::literals;
		return pybind11::dict("rmin"_a= r0, "rmax"_a = r1, "theta"_a = theta_max,"tau_theta"_a = tau_max);
	},
		m_grid
	);
}

pybind11::dict Py_EL_Grid::_get_traj_all_inter(double e, double l) {
	return std::visit(
		[e, l]<_GRID_EL_TMPL_>(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr)
	{
		typename evdm::Body<_T1>::LE_func_t const& LE_s = m_gr.getBodyLE();
		auto L = l * LE_s.i_lm(-e);
		size_t i = m_gr.Grid->inner(0).LinearIndex(
			m_gr.Grid->inner(0).pos(e, l)
		);
		auto const& m_traj_p = (*m_gr._TrajPools)[i];
		evdm::Body<_T1> const& B = *m_gr.body;
		auto _F2 = B._DD_F_C(1);
		auto LEf = m_gr.LE();
		auto TinFunc = make_bin_traj_pool_Tin_func(
			m_traj_p, LEf, _F2);
		auto ToutFunc = make_bin_traj_pool_Tout_func(
			m_traj_p, LEf);

		auto Lme = LEf(-e);
		auto [u0, u1, theta_max] = TinFunc.u0_u1_theta1(e, l, Lme);
		auto tau_max = TinFunc.tin_theta(e, l);
		auto r0 = (u0 > 0 ? std::sqrt(u0) : -std::sqrt(-u0));
		auto r1 = (u1 > 0 ? std::sqrt(u1) : -std::sqrt(-u1));

		using namespace pybind11::literals;
		return pybind11::dict(
			"rmin"_a = r0, "rmax"_a = r1, 
			"theta"_a = theta_max, "tau_theta"_a = tau_max
		);
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

std::variant<size_t, std::pair<size_t, size_t>> Py_EL_Grid::get_index(
	double e, double L, bool linear,bool hidden
) const{
	return std::visit([e, L, linear,hidden]<_GRID_EL_TMPL_>
	(const  evdm::EL_Grid<_GRID_EL_PARS_> & m_gr)->
		std::variant<size_t, std::pair<size_t, size_t>>
	{
		const auto& grid = m_gr.Grid->inner(0);
		auto LE_f = m_gr.LE();
		auto l = (hidden ? L : L / LE_f(-e));
		auto IJ = grid.pos(e, l);
		std::variant<size_t, std::pair<size_t, size_t>>V;
		if (linear) {
			V = (size_t)grid.LinearIndex(IJ);
		} else {
			V = std::pair<size_t, size_t>(IJ.i, IJ.m);
		}
		return V;
	},m_grid);
}
size_t Py_EL_Grid::get_e_index(double e)const {
	return std::visit([e]<_GRID_EL_TMPL_>
	(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr)->size_t
	{
		const auto& grid_e = m_gr.Grid->inner(0).grid();
		return grid_e.pos(e);
	}, m_grid);
}

pybind11::dict Py_EL_Grid::get_el(
	std::variant<size_t, std::pair<size_t, size_t>> index)const {

	return std::visit([]<_GRID_EL_TMPL_,typename T>
	(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr,T m_index)->pybind11::dict
	{
		grob::MultiIndex<size_t,size_t> real_index;
		const auto& grid = m_gr.Grid->inner(0);
		if constexpr (std::is_same_v<decltype(m_index), size_t>) {
			real_index = grid.FromLinear(m_index);
		}
		else {
			real_index = grob::MultiIndex<size_t, size_t>(m_index.first,m_index.second);
		}
		auto m_bin = grid[real_index];
		auto [e0, e1] = std::get<0>(m_bin);
		auto [l0, l1] = std::get<1>(m_bin);
		auto Lm0 = m_gr.LE()(-e0);
		auto Lm1 = m_gr.LE()(-e1);

		using namespace pybind11::literals;
		return pybind11::dict(
			"e0"_a = e0, "e1"_a = e1,
			"l0"_a = l0, "l1"_a = l1,
			"L00"_a = l0 * Lm0, "L01"_a = l1 * Lm0, "L10"_a = l0 * Lm1, "L11"_a = l1 * Lm1,
			"Lm0"_a = Lm0, "Lm1"_a = Lm1);
	}, m_grid,index);

}

pybind11::array Py_EL_Grid::get_E_array()const {
	return std::visit(
		[]<_GRID_EL_TMPL_>
		(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr)->pybind11::array
	{
		return make_py_array(m_gr.Grid->inner(0).grid().unhisto());
	},m_grid);
}
pybind11::array Py_EL_Grid::get_L_array(size_t index, bool hidden)const {
	return std::visit(
		[index, hidden]<_GRID_EL_TMPL_>
		(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr)->pybind11::array
	{
		_T2 Lm = m_gr.LE()(-m_gr.Grid->inner(0).grid()[index].center());
		auto const& m_grid = m_gr.Grid->inner(0).inner(index).unhisto();
		return make_py_array(
			grob::as_container(
				[hidden, Lm, &m_grid](size_t i)->_T2 {return hidden ? m_grid[i] : Lm * m_grid[i]; },
				m_grid.size()
			)
		);
	}, m_grid);
}

struct id_func_t {
	template <typename T>
	inline const T operator()(T x) const{
		return x;
	}
};
template <typename Grid1_t,typename Grid2_t,typename LE_Functype = id_func_t>
std::pair<pybind11::array_t<typename Grid1_t::value_type>,
		pybind11::array_t<typename Grid2_t::value_type>>
	mesh_np_grids(Grid1_t const& Grid1, Grid2_t const& Grid2, LE_Functype&& LEmax_f = id_func_t{}) {
	using T1 = typename Grid1_t::value_type;
	using T2= typename Grid2_t::value_type;
	namespace py = pybind11;
	size_t N1 = Grid1.size();
	size_t N2 = Grid2.size();
	py::array_t<T1> X(py::array::ShapeContainer({ N1,N2 }));
	py::array_t<T2> Y(py::array::ShapeContainer({ N1,N2 }));
	for (size_t i = 0; i < N1; ++i) {
		for (size_t j = 0; j < N2; ++j) {
			X.mutable_at(i, j) = Grid1[i];
			Y.mutable_at(i, j) = Grid2[j]* LEmax_f(Grid1[i]);
		}
	}
	return std::make_pair(std::move(X), std::move(Y));
}


static const std::vector<const char *> PNameFunctorNames = {
	"T0", "T1", "Tin0", "Tin1", "Tout0","Tout1", 
	"rmin", "rmax", "theta0","theta1", "Tinth"
};


std::function<double(double, double)> 
	Py_EL_Grid::getTraj_callFunc(
		pybind11::handle output_param_name
	) const {

	size_t my_functor_t;
	std::string pname = output_param_name.cast<std::string>();
	size_t pi = 0;
	for (; pi < PNameFunctorNames.size(); ++pi) {
		if (pname == PNameFunctorNames[pi]) {
			my_functor_t =  pi;
			break;
		}
	}
	if (pi == PNameFunctorNames.size()) {
		throw pybind11::value_error("unexpected param name");
	}
	
	return std::visit([my_functor_t](const auto& grid) {
		auto el_find = [&mel_grid = grid.getLE_inner_grid()](double e, double l) {
			return mel_grid.LinearIndex(mel_grid.pos(e, l));
		};
		auto const & vpools = grid.TrajPools();
		std::function<double(double, double)> _default = [](double, double)->double {return 42; };
		return evdm::make_choose_alt<std::function<double(double,double)>>(
			my_functor_t,
			//getT0
			[&vpools, el_find,LEf=grid.LE()](double e, double l) {
				auto const & mpool = vpools[el_find(e, l)];
				return evdm::bin_traj_pool_func(
					mpool, evdm::make_traj_tfull_getter(LEf)
				)(e, l);
			},
			//getT1
			[&vpools, el_find, LEf = grid.LE(),
				_FC = grid.body->_DD_F_C(1)](double e, double l)->double {
				auto const& mpool = vpools[el_find(e, l)];
				return evdm::make_bin_traj_pool_Tin_func(
						mpool, LEf, _FC
					)(e, l, LEf(-e)) + evdm::make_bin_traj_pool_Tout_func(
						mpool, LEf
					)(e, l);
			},
			//getTin0 = 
			[&vpools, el_find, LEf = grid.LE()](double e, double l) {
				auto const& mpool = vpools[el_find(e, l)];
				return evdm::bin_traj_pool_func(
					mpool, evdm::make_traj_tin_getter(LEf)
				)(e, l);
			},
			//getTin1
			[&vpools, el_find, LEf = grid.LE(),
				_FC = grid.body->_DD_F_C(1)](double e, double l)->double {
				auto const& mpool = vpools[el_find(e, l)];
				return evdm::make_bin_traj_pool_Tin_func(
					mpool, LEf, _FC
				)(e, l, LEf(-e));
			},
			//getTout0 = 
			[&vpools, el_find, LEf = grid.LE()](double e, double l) {
				auto const& mpool = vpools[el_find(e, l)];	
				return evdm::bin_traj_pool_func(
					mpool, evdm::make_traj_tout_getter(LEf)
				)(e, l);
			},
			//getTout1
			[&vpools, el_find, LEf = grid.LE()](double e, double l) {
				auto const& mpool = vpools[el_find(e, l)];
				return evdm::make_bin_traj_pool_Tout_func(
					mpool, LEf
				)(e, l);
			},
			//getrmin
			_default,
			//getrmax
			_default,
			//gettheta0
			[&vpools, el_find, LEf = grid.LE()](double e, double l) {
				auto const& mpool = vpools[el_find(e, l)];
				return evdm::bin_traj_pool_func(
					mpool, evdm::traj_info_transform_t<evdm::tit_ext_v::theta>{}
				)(e, l);
			},
			//gettheta1
			[&vpools, el_find, LEf = grid.LE(),
				_FC = grid.body->_DD_F_C(1)](double e, double l)->double {
				auto const& mpool = vpools[el_find(e, l)];
				return evdm::make_bin_traj_pool_Tin_func(
					mpool, LEf, _FC
				).theta1(e, l,LEf(-e));
			},
			_default,
			//getTinTheta
			_default
		);
	}, m_grid);

	return [](double,double) ->double {return 0; };
}

std::function<double(double)> 
	Py_EL_Grid::getLE_callFunc() const 
{
	return std::visit([](auto const& grid)->std::function<double(double)> {
		return [&grid, LEf = grid.LE()](double e)->double {
			return LEf(-e);
		};
	}, m_grid);
}

pybind11::tuple Py_EL_Grid::getRminRmax(double e, double l) const {
	return std::visit([e,l](auto const& grid) {
		const auto& B = *grid.body;
		const auto& LE = grid.LE();
		auto L = LE(-e) * l;
		auto [rmin, rmax] = B.find_rmin_rmax(-e,L*L,grid.getBodyLE());
		return pybind11::make_tuple(rmin, rmax);
	}, m_grid);
}

pybind11::list Py_EL_Grid::getTrajTFuncs(pybind11::handle output_param_name,bool LE_mult) const
{
	namespace py = pybind11;
	const std::set<std::string> TinNames = {
		"tin","TIN","t_in","T_in","Tin","T_IN"
	};
	const std::set<std::string> ToutNames = {
		"tout","TOUT","t_out","T_out","Tout","T_OUT"
	};
	const std::set<std::string> TallNames = {
		"t","T","Tall","T_all","TALL","T_ALL","t_all","tall"
	};
	const std::set<std::string> RminNames = {
		"rmin","rm"
	};
	const std::set<std::string> RmaxNames = {
		"rmax","rp"
	};
	const std::set<std::string> ThetaNames = {
		"theta","theta_max","ThetaMax","thetamax","maxtheta"
	};
	return std::visit(
		[&]<typename Bvt,typename Gvt, evdm::GridEL_type gt>
			(evdm::EL_Grid<Bvt, Gvt, gt> const& elg) 
		{
			auto TinGetter = evdm::make_traj_tin_getter(elg.LE());
			auto ToutGetter = evdm::make_traj_tout_getter(elg.LE());
			auto TallGetter = evdm::make_traj_tfull_getter(elg.LE());
			auto RminGetter = evdm::traj_info_transform_t<evdm::tit_ext_v::rmin>{};
			auto RmaxGetter = evdm::traj_info_transform_t<evdm::tit_ext_v::rmax>{};
			auto RthetaGetter = evdm::traj_info_transform_t<evdm::tit_ext_v::theta>{};

			size_t var_index;
			auto pname = output_param_name.cast<std::string>();
			if (TinNames.contains(pname)) {
				var_index = 0;
			}
			else if (ToutNames.contains(pname)) {
				var_index = 1;
			}
			else if (TallNames.contains(pname)) {
				var_index = 2;
			} 
			else if (RminNames.contains(pname)) {
				var_index = 3;
			}
			else if (RmaxNames.contains(pname)) {
				var_index = 4;
			}
			else if (ThetaNames.contains(pname)) {
				var_index = 5;
			}
			else {
				throw py::value_error("unknown par name");
			}
			auto Getter =
				evdm::make_variant_alt(var_index,
					TinGetter, ToutGetter, TallGetter,RminGetter,
					RmaxGetter, RthetaGetter);
			

			std::vector<evdm::bin_traj_pool_t<Gvt>> const &TrPl = elg.TrajPools();

			return std::visit([&](auto const& getter) {
				pybind11::list Pools;
				using namespace pybind11::literals;
				for (evdm::bin_traj_pool_t<Gvt> const& pool : TrPl) {
					auto [X, Y] = mesh_np_grids(
						pool.Grid.grid(), pool.Grid.inner(0),
						[&, LEf = elg.LE()](auto e)->decltype(e) {if (!LE_mult) return 1; else return LEf(-e); });
					py::array_t<Gvt> Z(py::array::ShapeContainer{X.shape(0),X.shape(1)});
					
					auto GetterFunc = evdm::bin_traj_pool_func(pool, getter);
					auto z_ptr = Z.mutable_data();
					for(size_t i=0;i< GetterFunc.Values.size();++i){
						z_ptr[i] = GetterFunc.Values[i];
					}
					Pools.append(py::dict("z"_a = Z, "x"_a = X, "y"_a = Y));
				}
				return Pools;
			}, Getter);
		}, m_grid
	);
}



Py_EL_Grid CreateELGrid(
	Py_BodyModel const& BM,
	size_t ptypes, size_t Ne,
	pybind11::handle Nl_func_or_size,
	pybind11::handle RhoE_func_opt,
	pybind11::handle RhoL_func_opt,
	std::string const & dtype,
	pybind11::kwargs const & kwargs) 
{
	DEBUG1(ptypes);
	DEBUG1(Ne);

	auto _kwarg = [&kwargs](const char* argname, auto default_value)->decltype(default_value) {
		if(kwargs.contains(argname))
			return kwargs[argname].cast<decltype(default_value)>();
		else
			return default_value;
	};
	auto TPPars = evdm::TrajPoolInitParams_t(
		_kwarg("static_tp", false),
		_kwarg("eps", 0.05),
		_kwarg("tbin", 100),
		_kwarg("ne_max", 4),
		_kwarg("nl_max", 4)
	);

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
					return pybind11::cast<T>(RhoL_func_opt(t_e,t_l));
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
							std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>{}, 
							TPPars);
					}
					else {
						if (RhoL_func_opt.is_none()) 
						{
							DEBUG1("RhoL_func_opt is none");
							auto _RhoL_func = [&](T t_e, T t_l)->T {
								return 1;
							};
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, _RhoL_func, marker,
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{},
								TPPars);
						}
						else {
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, RhoL_func, marker, 
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{},
								TPPars);
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
							std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>{},
							TPPars);
					}
					else {
						if (RhoL_func_opt.is_none()) {
							DEBUG1("RhoL_func_opt is none");
							auto _RhoL_func = [&](T t_e, T t_l)->T {
								return 1;
							};
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, _RhoL_func, marker, 
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{},
								TPPars);
						}
						else {
							return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, RhoE_func, RhoL_func, marker, 
								std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>{},
								TPPars);
						}
					}
				}
			}, BM.m_body);
	};
	DEBUG1(__PRETTY_FUNCTION__);
	if (dtype == "float") {
		DEBUG1("float");
#ifdef GRID_EL_USE_FLOAT
		return ret_action(Py_EL_Grid::m_type_marker<float>{});
#else
		throw pybind11::type_error("unsupported grid value type 'float'");
#endif
	}
	else if (dtype == "double") {
		DEBUG1("double");
#ifdef GRID_EL_USE_DOUBLE
		return ret_action(Py_EL_Grid::m_type_marker<double>{});
#else
		throw pybind11::type_error("unsupported grid value type 'float'");
#endif
	} else {
		throw pybind11::type_error("wrong data type: " + dtype + ", expect float or double");
	}

}

void Py_EL_Grid::add_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;
	using pyobj_ref = py::object const&;
	constexpr auto val_ref_pol = py::return_value_policy::reference_internal;
	py::class_<Py_EL_Grid>(m, "GridEL")
		.def(py::init(&CreateELGrid),
			"constructs E-L Grid\n\n"
			"Parameters:\n"
			"___________\n"
			"body : Body\n\tinstance of Body class\n"
			"ptypes : int\n\tnumber of particle types.\n"
			"Ne : int\n\t number of bins of E axis.\n"
			"Nl_func : int | function\n\tif number then number of bins of L axis,"
			"else a function from [0.0,1.0]->int,"
			"indicating number of bins in L grid depending on e,"
			"where 0.0 correspond to Emin = phi(0), and 1.0 -- to Emax = 0.\n"
			"RhoE : function\n\t"
			"[optional] bin density of e axis RhoE : e in [0.0,1.0]->float.\n"
			"RhoL : function\n\t"
			"[optional] bin density of l axis RhoL : (e,l) in [0.0,1.0]x[0.0,1.0]->float.\n"
			"dtype : string\n\tfloat or double.",
			py::arg("body"),
			py::arg("ptypes"),
			py::arg("Ne"),
			py::arg("Nl_func"),
			py::arg_v("RhoE", py::none()),
			py::arg_v("RhoL", py::none()),
			py::arg_v("dtype", "float")
		)
		.def_property_readonly("size", &Py_EL_Grid::size, "number of bins in grid")
		.def_property_readonly("dtype", &Py_EL_Grid::dtype, "type of grid")
		.def_property_readonly("ptypes", &Py_EL_Grid::ptypes, "number of particles")
		.def("__repr__", &Py_EL_Grid::repr)
		.def_property_readonly("body", &Py_EL_Grid::getBody)
		.def("plot", &Py_EL_Grid::getPlot,
			"gives array of shape (N,2,2):\n"
			"[ [[x_start_i,y_start_i],[x_end_i,y_end_i]],...]\n"
			"arrays could be plottted with matplotlib:\n"
			"lc = matplotlib.collections.LineCollection(result of this function)\n"
			"fig, ax = plt.subplots()\n"
			"ax.add_collection(lc)\n"
			"if is_internal true, then l from 0 to 1, else from 0, l(e).\n\n"
			"Parameters:\n"
			"___________\n"
			"is_internal : bool\n\tif true, plot in internal representation.",
			py::arg_v("is_internal", false)
		)
		.def("LE", &Py_EL_Grid::getLE,
			"return tuple (E array,L(E) array)"
		)
		.def("print_debug", &Py_EL_Grid::printd1)
		.def("plot_trajp", &Py_EL_Grid::getTrajTFuncs,
			"returns list L of dicts able to be plot in plotly of some"
			"parametrs of trajectory (Tin,Tout,Tfull,rmin,rmax,theta)\n"
			"pname name of plotted parameter: 'Tin', 'Tout' or 'T'\n"
			"(if ld = False, L would be plotted in hidden space)\n"
			"to plot list in plotly use go.Figure(data = [go.Surface(P) for P in RetList])",
			py::arg_v("pname", "T"),
			py::arg_v("ld", true)
		)
		.def("le_functor", &Py_EL_Grid::getLE_callFunc,
			"returns functor of lmax(e) function", val_ref_pol)
		.def("traj_functor", &Py_EL_Grid::getTraj_callFunc,
			"return functor of required param:\n"
			"T0,T1,Tin0,Tin1,Tout,rmin,rmax,theta,Tinth", val_ref_pol,
			py::arg_v("pname", "T")
		)
		.def("rmp",&Py_EL_Grid::getRminRmax,
			"returns (rmin,rmax)(e,l_undim)",
			py::arg("e"), py::arg("l"))
		.def("_get_index",&Py_EL_Grid::get_index,
			py::arg("e"), py::arg("l"), py::arg_v("linear",false), py::arg_v("hidden",false))
		.def("_get_el",&Py_EL_Grid::get_el, py::arg("index"))
		.def("_get_e_array", &Py_EL_Grid::get_E_array)
		.def("_get_l_array", &Py_EL_Grid::get_L_array, 
			py::arg("e_index"), py::arg_v("hidden",false))
		.def("_get_traj_all",&Py_EL_Grid::_get_traj_all,py::arg("e"), py::arg("l"))
		.def("_get_traj_all_inter", &Py_EL_Grid::_get_traj_all_inter, py::arg("e"), py::arg("l"));
}
