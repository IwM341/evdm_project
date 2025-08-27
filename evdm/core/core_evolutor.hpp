#pragma once

#include "body_potential.hpp"
#include "grid_variants.hpp"
#include "histo_norm.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>
#include "dynamics/dynamics.hpp"
#include "r_conversion.hpp"
#include "traj_pool.hpp"
#include <grob/grid.hpp>
#include <grob/linear_interpolator.hpp>


struct interpolator_idx {
	template <typename GridType, typename ContainerType, typename Point_t>
	inline constexpr static auto interpolate(GridType const& grid,
		ContainerType const& values,
		Point_t const& point)
	{
		values[point];
	}
};

template <typename T>
struct MCEvolutor {
	typedef grob::GridUniform<T> Grid_fp_t;

	typedef grob::GridRangeLight Grid_idx_t;

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
	> El_ptype_Grid_t;

	typedef grob::MultiGrid<
		Grid_idx_t,
		grob::ConstValueVector<
			grob::MultiGrid<
				Grid_idx_t,
				grob::ConstValueVector<
					grob::MultiGrid<
						Grid_fp_t,
						grob::ConstValueVector<Grid_fp_t>
					>
				>
			>
		>
	> El_ptype_element_Grid_t;

	typedef grob::MultiGrid<
		Grid_fp_t,
		grob::ConstValueVector<
			grob::MultiGrid<
				Grid_fp_t,
				grob::ConstValueVector<
					grob::MultiGrid<
						Grid_idx_t,
						grob::ConstValueVector<
							grob::MultiGrid<
								Grid_idx_t,
								grob::ConstValueVector<Grid_fp_t>
							>
						>
					>
				>
			>
		>
	> El_ptype_element_xitau_Grid_t;

	typedef grob::linear_interpolator linint;
	typedef interpolator_idx idxint;
	
	typedef grob::interpolator_product<linint, linint> EL_interpol;
	typedef grob::interpolator_product<
		idxint,
		EL_interpol
	> EL_ptype_interpol;

	typedef grob::interpolator_product<
		idxint,
		EL_ptype_interpol
	>  EL_ptype_element_interpol;

	grob::GridFunction<
		EL_interpol,
		ELGrid_t, std::vector<T>
	> P_total_El_func;

	/*
	grob::GridFunction<
		EL_interpol,
		ELGrid_t, std::vector<T>
	> P_total_El_func;
	*/


};