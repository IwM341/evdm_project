#pragma once
#include <grob/grid.hpp>
#include <grob/linear_interpolator.hpp>
#include "../body_potential.hpp"
#include "../form_factors.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>

#include "../traj_pool.hpp"





namespace grob {
	template <typename T>
	concept EigenDenseClass = requires(T && x) {
		requires std::is_base_of_v<
			Eigen::EigenBase<
				std::decay_t<decltype(x.derived())>
			>,
			std::decay_t<decltype(x.derived())>
		>;
	};

	template <EigenDenseClass T>
	inline auto make_slice(T && V, size_t shift, size_t size = 0)  {
		return V(Eigen::seqN(shift, size));
	}
};

namespace evdm {

	struct interpolator_idx {
		template <typename GridType, typename ContainerType, typename Point_t>
		inline constexpr static auto interpolate(GridType const& grid,
			ContainerType const& values,
			Point_t const& point)
		{
			values[point];
		}
	};

	struct interpolator_none {
		template <typename GridType, typename ContainerType, typename Point_t>
		inline constexpr static auto interpolate(GridType const& grid,
			ContainerType const& values,
			Point_t const& point)
		{
			throw std::runtime_error("non callable interpolator");
		}
	};

	namespace detail {
		template <typename GridFuncType, typename GridType>
		void setGrid(GridFuncType& F, GridType G) {
			F.Grid = std::move(G);
			F.Values.resize(F.Grid.size());
		}
	};

	
	
	template <typename T, typename FormFactor_t>
	struct ElementInfo_t {
		T M;
		FormFactor_t dF;
		Eigen::VectorX<T> Rd;
		size_t Nmk;
	};

	template <typename T,typename FormFactor_t>
	struct MCEvolutor {


		typedef grob::GridUniform<T> Grid_fp_t;

		typedef grob::GridRangeLight<int> Grid_idx_t;

		typedef grob::MultiGrid<
			Grid_fp_t,
			grob::ConstValueVector<Grid_fp_t>
		> ELGrid_t;

		typedef grob::MultiGrid<
			Grid_idx_t,
			grob::ConstValueVector<
			grob::MultiGrid<
			Grid_fp_t,
			grob::ConstValueVector<Grid_fp_t>
			>
			>
		> Grid_pin_El_t;
		//(ptype, E,l) => total probs for P(pin, E, l)

		typedef grob::MultiGrid<
			Grid_pin_El_t,
			grob::ConstValueVector<
			Grid_idx_t
			>
		> Grid_pin_El_pout_t;
		//+ P(pin, E, l -> pout) PCF array => P(pin, E, l)(pout)

		typedef grob::MultiGrid <
			Grid_idx_t,
				grob::ConstValueVector <
					grob::MultiGrid<
						Grid_pin_El_t,
						grob::ConstValueVector<Grid_idx_t>
				>
			>
		>Grid_pout_pin_El_element_t;
		//+ P(pin, E, l,pout -> scatterelement) PCF array => P(pout,{pin, E, l})(element)

		typedef grob::MultiGrid<
			Grid_idx_t,
			grob::ConstValueVector<grob::MultiGrid<
				Grid_idx_t,
				grob::ConstValueVector<grob::MultiGrid<
					Grid_pin_El_t,
					grob::ConstValueVector<Grid_fp_t>
				>>
			>>
		> Grid_element_pout_pin_El_taupar_t;
		// 
		//ptype, ptype, E,l => probs on elements, element sample functions

		template <typename It, typename Iu>
		using iprod = grob::interpolator_product<It, Iu>;

		typedef grob::linear_interpolator linint;
		typedef interpolator_idx idxint;

		typedef iprod<linint, linint> EL_interpol;
		typedef iprod<idxint, EL_interpol> EL_ptype_interpol;

		typedef iprod<
			idxint,
			EL_ptype_interpol
		>  EL_ptype_element_interpol;


		typedef grob::MultiGrid<
			ELGrid_t,
			grob::ConstValueVector<
			Grid_fp_t
			>
		> Grid_TrajPool_t;

		grob::GridFunction<
			interpolator_none,
			Grid_pin_El_t, Eigen::VectorX<T>
		> Func_pin_El_to_Ptotal;
		//    println(X(Eigen::seqN(0,5,2))[4]);

		grob::GridFunction<
			interpolator_none,
			Grid_pin_El_pout_t, Eigen::VectorX<T>
		> Func_pin_El_to_Ppout;

		grob::GridFunction<
			interpolator_none,
			Grid_pout_pin_El_element_t, Eigen::VectorX<T>
		> Func_pout_pin_El_to_CPFelement;

		grob::GridFunction<
			interpolator_none,
			Grid_element_pout_pin_El_taupar_t, Eigen::VectorX<T>
		> Func_element_pout_pin_El_to_CPFtau;

		evdm::Body<T> B;
		decltype (B.get_le(0)) LE;

		bool ValuesFilled;
		grob::GridFunction<
			interpolator_none,
			Grid_TrajPool_t, Eigen::VectorX<T>
		> TrajPool;


		std::vector<ElementInfo_t<T, FormFactor_t>> ElementInfos;

		MCEvolutor(
			evdm::Body<T> _B,
			size_t ptypes,size_t Ne,size_t Nl,size_t Ntraj,
			std::vector< ElementInfo_t<T, FormFactor_t> > m_elements
			
		):B(std::move(_B)), ElementInfos(std::move(m_elements)), LE(B.get_le(2*Ne)){

			ValuesFilled = false;
			auto Egrid = Grid_fp_t(B.Phi.Values[0], 0, Ne);
			auto lgrid = Grid_fp_t(0, 1, Nl);
			auto ElemetGrid = Grid_idx_t(ElementInfos.size());
			auto ptypeGrid = Grid_idx_t(ptypes);
			auto TrajPoolrGrid = Grid_fp_t(0, 1, Ntraj+1);

			Grid_pin_El_t Grid_pin_El = grob::mesh_grids(ptypeGrid, grob::mesh_grids(Egrid, lgrid));
			Grid_TrajPool_t Grid_TrajPool = grob::mesh_grids(grob::mesh_grids(Egrid, lgrid), TrajPoolrGrid);
			
			detail::setGrid(Func_pin_El_to_Ptotal, Grid_pin_El);
			
			Grid_pin_El_pout_t Grid_pin_El_pout = grob::mesh_grids(Grid_pin_El, ptypeGrid);


			detail::setGrid(Func_pin_El_to_Ppout, Grid_pin_El_pout);

			Grid_pout_pin_El_element_t
				Grid_pout_pin_El_element = 
				grob::mesh_grids(ptypeGrid, grob::mesh_grids(Grid_pin_El, ElemetGrid));

			detail::setGrid(Func_pout_pin_El_to_CPFelement, Grid_pout_pin_El_element);

			Grid_element_pout_pin_El_taupar_t 
				Grid_element_pout_pin_El_taupar = 
				grob::mesh_grids(ElemetGrid,grob::mesh_grids(ptypeGrid, grob::mesh_grids(Grid_pin_El, TrajPoolrGrid)));

			detail::setGrid(Func_element_pout_pin_El_to_CPFtau,Grid_element_pout_pin_El_taupar);

			detail::setGrid(TrajPool,Grid_TrajPool);

			auto const& gridEL = Grid_TrajPool.grid();
			int NeNl = gridEL.size();
			#pragma omp parallel for
			for (size_t k = 0; k < NeNl; ++k) {
				auto IJ = gridEL.FromLinear(k);
				auto [i, j] = IJ;
				auto [e, l] = gridEL[{i, j}];
				auto Lmax = LE.i_lm(-e);
				auto rm_e = LE.i_rm(-e);
				auto L = (Lmax * l);
				auto L2 = L * L;
				auto [rm, rp] = B.find_rmin_rmax(-e, L2, LE);
				auto [theta_max,tau_max,taus] = B.template get_internal_traj<T>(rm, rp, Ntraj);
				
				auto thetas = InverseTrajectory<interpolator_none>(taus, 30);
				auto M_traj = TrajPool.inner_slice(k);
				
				assert(thetas.size() == M_traj.size());
				std::copy(thetas.Values.begin(), thetas.Values.end(), M_traj.Values.begin());
				M_traj
			}

		}

		void FillValues() {
		
		}

		typedef Eigen::Matrix2<T> IC_t;
		std::pair<grob::MultiIndex<int, int>, IC_t> get_el_info(T e, T l) {
			auto MI = TrajPool.Grid.grid().pos(e, l);
			auto [i, j] = MI;
			T e0 = TrajPool.Grid.grid().grid()[i];
			T e1 = TrajPool.Grid.grid().grid()[i + 1];
			T l0 = TrajPool.Grid.grid().inner(0)[j];
			T l1 = TrajPool.Grid.grid().inner(0)[j + 1];

			T c0_e = downbound((e1 - e) / (e1 - e0),0);
			T c0_l = downbound((l1 - l) / (l1 - l0),0);

			IC_t M;
			M(0,0) = c0_e * c0_l;
			M(0,1) = c0_e * (1-c0_l);
			M(1, 0) = (1-c0_e)* c0_l;
			M(1, 1) = (1 - c0_e) * (1 - c0_l);
			return { MI ,M };
		}
	};

};
