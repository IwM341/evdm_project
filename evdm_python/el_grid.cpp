#include "grid_python.hpp"
#include "debugdef.hpp"
#include "grob_python.hpp"
#include <pybind11/functional.h>
#include <set>
#include <pybind11/stl.h>
#include <grob/object_serialization.hpp>
#include <variant>
//#include <format>


size_t Py_EL_Grid::size() const
{
	return std::visit([](auto const& ELgrid) {
			return ELgrid.Grid->inner(0).size();
		}, m_grid);
}

size_t Py_EL_Grid::size_e() const
{
	return std::visit(
		[]<_GRID_EL_TMPL_>(evdm::EL_Grid<_GRID_EL_PARS_> const& G) {
		return G.getLE_inner_grid().grid().size();
	},m_grid);
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
			VERBOSE_VAR1(Points.size());
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
			VERBOSE_VAR1(Ne);
			auto _grid_e = grob::GridUniform<T>(-EL_grid_w.body->Phi[0], 0, Ne);
			VERBOSE_VAR1(_grid_e);
			VERBOSE_VAR1(EL_grid_w.LE()(-_grid_e[0]));
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

pybind11::dict Py_EL_Grid::get_el_all(
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
pybind11::tuple Py_EL_Grid::get_el(
	std::variant<size_t, std::pair<size_t, size_t>> index)const {

	return std::visit([]<_GRID_EL_TMPL_, typename T>
		(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr, T m_index)->pybind11::tuple
	{
		grob::MultiIndex<size_t, size_t> real_index;
		const auto& grid = m_gr.Grid->inner(0);
		if constexpr (std::is_same_v<decltype(m_index), size_t>) {
			real_index = grid.FromLinear(m_index);
		}
		else {
			real_index = grob::MultiIndex<size_t, size_t>(
				m_index.first, m_index.second
			);
		}
		auto m_bin = grid[real_index];
		auto [e0, e1] = std::get<0>(m_bin);
		auto [l0, l1] = std::get<1>(m_bin);

		using namespace pybind11::literals;
		return pybind11::make_tuple(
			e0, e1,l0, l1);
	}, m_grid, index);

}
pybind11::array_t<size_t> Py_EL_Grid::indexes(pybind11::handle self) const {
	return std::visit(
		[self]<_GRID_EL_TMPL_>(evdm::EL_Grid<_GRID_EL_PARS_> const& G)->
		pybind11::array_t<size_t> 
	{
		return make_py_array_view(G.getLE_inner_grid().indexes(), self);
	},m_grid);
}

pybind11::array Py_EL_Grid::get_E_array(pybind11::handle self)const {
	return std::visit(
		[self]<_GRID_EL_TMPL_>
		(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr)->pybind11::array
	{
		if constexpr (_m_grid_t == evdm::GridEL_type::GridCVV) {
			return make_py_array_view(
				m_gr.Grid->inner(0).grid().unhisto(),
				self
			);
		}
		else {
			return make_py_array(
				m_gr.Grid->inner(0).grid().unhisto()
			);
		}
		
	},m_grid);
}
pybind11::array Py_EL_Grid::get_L_array(
	pybind11::handle self,int index, bool hidden
)const {
	return std::visit(
		[self,index_ = index, hidden]<_GRID_EL_TMPL_>
		(const  evdm::EL_Grid<_GRID_EL_PARS_> &m_gr)->pybind11::array
	{
		auto const& m_grid = m_gr.Grid->inner(0).inner(index_).unhisto();
		size_t size_e = m_gr.getLE_inner_grid().grid().size();
		auto index = index_;
		if (index < 0) {
			index = index + size_e;
		}
		if (index >= size_e || index < 0) {
			throw pybind11::index_error(
				"index" + std::to_string(index) + "out of range"
			);
		}

		if (hidden) {
			if constexpr (_m_grid_t == evdm::GridEL_type::GridCVV) {
				return make_py_array_view(m_grid, self);
			}
			else {
				return make_py_array(m_grid);
			}
		}
		_T2 Lm = m_gr.LE()(-m_gr.Grid->inner(0).grid()[index].center());
		return make_py_array(
			grob::as_container(
				[hidden, Lm, &m_grid](size_t i)->_T2 {
					return hidden ? m_grid[i] : Lm * m_grid[i]; 
				},
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
	/*
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
	*/
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


pybind11::list Py_EL_Grid::getTrajTFuncs(
	pybind11::handle output_param_name,
	bool LE_mult) const
{
	return pybind11::list();
	/*
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
	);*/
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
	VERBOSE_VAR1(ptypes);
	VERBOSE_VAR1(Ne);

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
		//VERBOSE_VAR1(type_name<T>());
		return std::visit(
			[&](auto const& BM_p) {
				//VERBOSE_VAR1(type_name<decltype(BM_p->get_vtype())>());
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
					VERBOSE_VAR1("nl >=0");
					auto Nl_func = [nl_size](auto t_e) {return nl_size; };
					if (RhoE_func_opt.is_none()) 
					{
						VERBOSE_VAR1("RhoE_func_opt is none");
						return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, 0, 0, marker, 
							std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>{}, 
							TPPars);
					}
					else {
						if (RhoL_func_opt.is_none()) 
						{
							VERBOSE_VAR1("RhoL_func_opt is none");
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
					VERBOSE_VAR1("nl determined by function");
					auto Nl_func = [nl_size,&Nl_func_or_size](auto t_e)->size_t {
						return pybind11::cast<double>(Nl_func_or_size(t_e)); 
					};
					if (RhoE_func_opt.is_none()) 
					{
						VERBOSE_VAR1("RhoE_func_opt is none");
						return Py_EL_Grid(BM_p, ptypes, Ne, Nl_func, 0, 0, marker, 
							std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>{},
							TPPars);
					}
					else {
						if (RhoL_func_opt.is_none()) {
							VERBOSE_VAR1("RhoL_func_opt is none");
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
	if (dtype == "float") {
		VERBOSE_VAR1("float");
#ifdef GRID_EL_USE_FLOAT
		return ret_action(Py_EL_Grid::m_type_marker<float>{});
#else
		throw pybind11::type_error("unsupported grid value type 'float'");
#endif
	}
	else if (dtype == "double") {
		VERBOSE_VAR1("double");
#ifdef GRID_EL_USE_DOUBLE
		return ret_action(Py_EL_Grid::m_type_marker<double>{});
#else
		throw pybind11::type_error("unsupported grid value type 'float'");
#endif
	} else {
		throw pybind11::type_error("wrong data type: " + dtype + ", expect float or double");
	}

}


struct NpDictSerializator {

	template <typename T>
	struct serializable_s : std::false_type {};

	template <typename T>
	struct deserializable_s : std::false_type {};

	template <typename T,typename...Args>
	struct serializable_s<
		std::vector<T, Args...>
	> : std::is_arithmetic<T> {};

	template <typename T,typename...Args>
	struct deserializable_s<
		std::vector<T, Args...>
	> : std::is_arithmetic<T> {};


	template <typename T>
	pybind11::object MakePrimitive(T const & value) {
		return pybind11::cast(
			value,
			pybind11::return_value_policy::automatic_reference
		);
	}
	template <typename Keys_t,typename Values_t>
	pybind11::object MakeDict(Keys_t && keys, Values_t && values) {
		pybind11::dict m_dict;
		for (size_t i = 0; i < keys.size(); ++i) {
			if constexpr (
				std::is_same_v<
					std::decay_t<decltype(keys[i])>,
					std::string
				>
			) {
				m_dict[keys[i].c_str()] = values[i];
			}
			else if constexpr (
				std::is_same_v<
					std::decay_t<decltype(keys[i])>, 
					std::string_view
				>
			) {
				std::string X(keys[i].begin(), keys[i].end() );
				m_dict[X.c_str()] = values[i];
			} else {
				m_dict[keys[i]] = values[i];
			}
		}
		return m_dict;
	}

	template <typename Values_t>
	pybind11::object MakeArray(Values_t&& values) {
		typedef std::decay_t<decltype(values[0])> value_type;
		if constexpr (std::is_fundamental_v<value_type>) {
			return pybind11::cast(make_py_array(values));
		} else {
			pybind11::list L;
			size_t N = values.size();
			for (size_t i = 0; i < N;++i) {
				L.append(values[i]);
			}
			return L;
		}
	}

	template <typename...Ts>
	pybind11::object Make(std::vector<Ts...> const& value) {
		return make_py_array(value);
	}

	template <typename T, typename...Args>
	std::vector<T, Args...> get_impl(
		std::type_identity<std::vector<T,Args...>>,
		pybind11::handle Obj) 
	{
		try {
			pybind11::array_t<T> m_array = Obj.cast<pybind11::array_t<T>>();
			if (m_array.ndim() != 1) {
				throw pybind11::index_error("number of dimentions in array should be 1 to cast to std::vector");
			}
			const T* _ptr = m_array.data();
			size_t _stride = m_array.strides()[0]/sizeof(T);
			size_t _size = m_array.shape()[0];
			std::vector<T, Args...> Ret(_size);
			for (size_t i = 0; i < _size; ++i) {
				Ret[i] = *(_ptr + i * _stride);
			}
			return Ret;
		}
		catch (pybind11::cast_error&) {
			
		}

		std::vector<T, Args...> Ret;
		for (auto it = Obj.begin(); it != Obj.end();++it) {
			Ret.push_back(it->cast<T>());
		}

		return Ret;
	}
	

	template <typename T>
	T GetPrimitive(pybind11::handle Obj) {
		return pybind11::cast<T>(Obj);
	}
	template <typename T>
	T Get(pybind11::handle Obj) {
		return get_impl(std::type_identity<T>{}, Obj);
	}


	auto BeginArray(pybind11::handle Obj) {
		return Obj.begin();
	}
	auto EndArray(pybind11::handle Obj) {
		return Obj.end();
	}

	auto BeginDict(pybind11::handle Obj) {
		return pybind11::cast<pybind11::dict>(Obj).begin();
	}
	auto EndDict(pybind11::handle Obj) {
		return pybind11::cast<pybind11::dict>(Obj).end();
	}
	template <typename Iter_type>
	auto GetKey(Iter_type it) {
		return it->first.template cast<std::string_view>();
	}
	template <typename Iter_type>
	auto GetValue(Iter_type it) {
		return it->second;
	}
	template <typename Iter_type>
	auto GetItem(Iter_type it) {
		return *it;
	}
	template <typename T>
	pybind11::handle GetPrimitive(pybind11::handle Obj, T & value) {
		return value = Obj.cast<T>();
		grob::GridUniform<T>::Serialize;
	}
};

pybind11::dict Py_EL_Grid::get_object(pybind11::handle self) {
	using namespace pybind11::literals;

	return std::visit([this]<class B_vt,class Gr_vt,evdm::GridEL_type _m_grid_t>(
		evdm::EL_Grid<B_vt, Gr_vt, _m_grid_t> const& grid
		) 
	{
		NpDictSerializator  S;
		auto m_dict = pybind11::dict(
			"type"_a = "evdm.GridEL",
			"grid_t"_a = type_name<Gr_vt>()
		);
		// Saving grid values
		if constexpr (_m_grid_t == evdm::GridEL_type::GridCUU) {
			m_dict["gtype"] = "CUU";
		}
		else {
			m_dict["gtype"] = "CVV";
		}
		m_dict["body"] = getBody();
		m_dict["grid"] = grid.Grid->Serialize(S);
		m_dict["LE"] = grid._LE->Serialize(S);
		m_dict["Trajs"] = grob::Serialize(*grid._TrajPools, S);
		return m_dict;
	}, m_grid);

}

Py_EL_Grid Py_EL_Grid::from_dict(Py_BodyModel BM, pybind11::dict const& grid_dict) {
	if (grid_dict["type"].cast<std::string_view>() != "evdm.GridEL") {
		throw pybind11::type_error(
			"bad attempt to restore evdm.GridEL from not evdm.GridEL value"
		);
	}

	auto Grid_var_t = type_from_str(
		grid_dict["grid_t"].cast<std::string_view>(),
		std::type_identity<std::tuple<GRID_EL_TYPE_LIST>>{}
	);
	constexpr static auto _CUU = evdm::GridEL_type::GridCUU;
	constexpr static auto _CVV = evdm::GridEL_type::GridCVV;
	
	std::string_view _kind_str = grid_dict["gtype"].cast<std::string_view>();
	size_t _kind = (_kind_str == "CUU" ? 0 : (_kind_str == "CVV" ? 1 : 2));

	std::variant< GRID_KIND_TYPE_LIST> m_grid_type;
	if (_kind >= std::tuple_size_v<std::tuple<GRID_KIND_TYPE_LIST>>) {
		throw pybind11::type_error("unsupported gtype");
	}
	m_grid_type = evdm::make_variant<
		std::variant< GRID_KIND_TYPE_LIST>
	>(_kind, [](auto m_T) {
		return typename decltype(m_T)::type{};
	});
	
	return std::visit(
		[&]<class B_vt,class Gr_vt, evdm::GridEL_type _grid_kind>(
			evdm::BodyModel<B_vt> const &body, 
			std::type_identity< Gr_vt>, 
			evdm::GridEL_type_t<_grid_kind>)->Py_EL_Grid
	{
		typedef evdm::EL_Grid<B_vt, Gr_vt, _grid_kind> EL_Grid_full_t;
		
		typedef typename EL_Grid_full_t::GridEL_t mGridEL_t;
		typedef typename EL_Grid_full_t::LE_func_t LE_func_t;
		typedef typename EL_Grid_full_t::TrajPool_t TrajPool_t;

		NpDictSerializator S;
		
		mGridEL_t m_grid = mGridEL_t::DeSerialize(grid_dict["grid"], S);

				
		if (grid_dict.contains("LE")) {
			LE_func_t LE_Func = LE_func_t::DeSerialize(grid_dict["LE"], S);
			if (grid_dict.contains("Trajs")) {
				TrajPool_t Trajs = 
					grob::DeSerialize< TrajPool_t>(grid_dict["Trajs"], S);
				return Py_EL_Grid(
					evdm::EL_Grid<B_vt, Gr_vt, _grid_kind>(
						body, std::move(m_grid),
						std::move(LE_Func),
						std::move(Trajs)
						)
				);
			}
			else {
				return Py_EL_Grid(
					evdm::EL_Grid<B_vt, Gr_vt, _grid_kind>(
						body, std::move(m_grid),
						std::move(LE_Func)
					)
				);
			}
		}
		else {
			size_t Ne = pyget<size_t>( 2*m_grid.size_e(),grid_dict,"NLe");
			if (grid_dict.contains("Trajs")) {
				TrajPool_t Trajs =
					grob::DeSerialize< TrajPool_t>(grid_dict["Trajs"], S);
				return Py_EL_Grid(
					evdm::EL_Grid<B_vt, Gr_vt, _grid_kind>(
						body, std::move(m_grid),
						Ne,
						std::move(Trajs)
						)
				);
			}
			else {
				return Py_EL_Grid(
					evdm::EL_Grid<B_vt, Gr_vt, _grid_kind>(
						body, std::move(m_grid),Ne)
				);
			}
		}
		
	}, BM.m_body, Grid_var_t, m_grid_type);
		
	
}

Py_EL_Grid Py_EL_Grid::from_dict1(pybind11::dict const& grid_dict) {
	pybind11::handle B = grid_dict["body"];
	if (pybind11::isinstance<pybind11::dict>(B)) {
		return Py_EL_Grid::from_dict(
			Py_BodyModel::from_dict(
				B.cast<pybind11::dict>()
			), grid_dict
		);
	} else {
		return Py_EL_Grid::from_dict(
			B.cast<Py_BodyModel>(), grid_dict
		);
	}
}
pybind11::array GetBins(Py_EL_Grid const& self) {
	return std::visit(
		[]<_GRID_EL_TMPL_>(
			evdm::EL_Grid<_GRID_EL_PARS_> const& grid
		)->pybind11::array {
			typedef typename evdm::EL_Grid<_GRID_EL_PARS_>::GEL_vt T;
			auto const& grid_el = grid.getLE_inner_grid();
			pybind11::array_t<T> m_array(pybind11::array::ShapeContainer({ grid_el.size(),(size_t)4}));
			T* data = m_array.mutable_data();
			size_t i = 0;
			for (auto [e0e1, l0l1] : grid_el) {
				if (4 * i + 3 >= m_array.size()) {
					throw pybind11::index_error("in bins property c++ error!");
				}
				data[4 * i] = e0e1.left;
				data[4 * i + 1] = e0e1.right;
				data[4 * i + 2] = l0l1.left;
				data[4 * i + 3] = l0l1.right;
				++i;
			}
			return m_array;
		}, self.m_grid
	);
}
void Py_EL_Grid::add_to_python_module(pybind11::module_& m)
{
	namespace py = pybind11;
	using pyobj_ref = py::object const&;
	constexpr auto val_ref_pol = py::return_value_policy::reference_internal;

	auto m_lass = py::class_<Py_EL_Grid>(m, "GridEL")
		.def(
			py::init(&CreateELGrid),
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
		.def(
			"refine",&Py_EL_Grid::refine,
			"make refeined grid by dividing E step by Er"
			"and L step by Lr",
			py::arg("Er"), py::arg("Lr")
		)
		.def(py::init(&Py_EL_Grid::from_dict))
		.def(py::init(&Py_EL_Grid::from_dict1))
		.def(
			"__getstate__", 
			[](py::handle self)->py::dict { 
				return py::cast<Py_EL_Grid&>(self).get_object(self); 
			}
		)
		.def(
			"__setstate__", 
			[](py::handle self, py::dict state) { 
				self.attr("__init__")(state);
			}
		)
		.def(
			"to_object", 
			[](py::handle self)->py::dict { 
				return py::cast<Py_EL_Grid&>(self).get_object(self); 
			},
			"serialize to dict"
		).def_property_readonly(
			"size",
			&Py_EL_Grid::size,
			"number of bins in EL grid"
		).def_property_readonly(
			"size_e",
			&Py_EL_Grid::size_e,
			"number of bins in E subgrid"
		).def_property_readonly(
			"full_size",
			[](Py_EL_Grid const& _G) { return _G.size() * _G.ptypes(); },
			"full size of grid (considering all ptypes)"
		).def_property_readonly(
			"dtype",
			&Py_EL_Grid::dtype,
			"type of grid"
		).def_property_readonly(
			"ptypes",
			&Py_EL_Grid::ptypes,
			"number of particles"
		).def("__repr__", &Py_EL_Grid::repr)
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
		.def_property_readonly("LE", &Py_EL_Grid::getLE,
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
		.def_property_readonly("le_functor", &Py_EL_Grid::getLE_callFunc,
			"returns functor of lmax(e) function", val_ref_pol)
		.def("traj_functor", &Py_EL_Grid::getTraj_callFunc,
			"return functor of required param:\n"
			"T0,T1,Tin0,Tin1,Tout,rmin,rmax,theta,Tinth", val_ref_pol,
			py::arg_v("pname", "T")
		)
		.def("rmp", &Py_EL_Grid::getRminRmax,
			"returns (rmin,rmax)(e,l_undim)",
			py::arg("e"), py::arg("l"))
		.def("_get_index", &Py_EL_Grid::get_index,
			py::arg("e"), py::arg("l"), py::arg_v("linear", false), py::arg_v("hidden", false))
		.def("_get_el", &Py_EL_Grid::get_el_all, py::arg("index"))
		.def("__getitem__", &Py_EL_Grid::get_el)
		.def_property_readonly(
			"indexes",
			[](py::handle self) {
				return self.cast<Py_EL_Grid&>().indexes(self);
		
			},
			"indexes of bins, where l = 0"
		)
		.def_property_readonly("bins", &GetBins,"array of all bins size*4: [[e0,e1,l0,l1]_0,...]")
		.def_property_readonly(
			"Epoints", 
			[](pybind11::handle self) {
				return self.cast<Py_EL_Grid>().get_E_array(self);
			},
			"array of energy points"
		).def("Lpoints", [](pybind11::handle self, size_t index) {
				return self.cast<Py_EL_Grid>().get_L_array(
					self, index, false
				);
			},py::arg("index")
		).def(
			"lpoints", 
			[](pybind11::handle self,size_t index) {
				return self.cast<Py_EL_Grid>().get_L_array(
					self, index,true
				);
			}, py::arg("index"),
			"array of l points at enegy with index index"
		)
		.def("_get_traj_all",&Py_EL_Grid::_get_traj_all,py::arg("e"), py::arg("l"))
		.def("_get_traj_all_inter", &Py_EL_Grid::_get_traj_all_inter, py::arg("e"), py::arg("l"));
}

Py_EL_Grid Py_EL_Grid::refine(size_t Ne, size_t Nl) const {
	return std::visit(
		[Ne, Nl](auto const& grid){
			return Py_EL_Grid(grid.refine(Ne,Nl));
		},m_grid
	);
}