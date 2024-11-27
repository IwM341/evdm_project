#pragma once
#include "capture_python.hpp"

struct Py_Matrix {
	Matrix_Variant_t m_matrix;
	//Matrix_Pair_Inst<T_i ...> m_matrix
	std::vector<scatter_event_info> events;

	template <typename T>
	Py_Matrix(Py_EL_Grid const& mGridEL,
		Py_EL_Grid::m_type_marker<T> _type/*,
		size_t padding = 1*/
	) :
		m_matrix(
			std::visit([](auto const& _grid)->Matrix_Variant_t {
				return std::make_pair(
					evdm::make_Matrix<T>(_grid, 1),
					evdm::make_Distribution<T>(_grid, false,1)
				);
				}, mGridEL.m_grid)
		)
	{
	}

	/*
	/// @brief constructor from object
	Py_Matrix(
		Py_EL_Grid const& mGridEL,
		pybind11::dict const& saved_obj);
	*/

	/// @brif constructor from GridMatrix
	template <
		typename T, typename Bt, typename Gt, evdm::GridEL_type g_type
	>
	Py_Matrix(evdm::GridMatrix<T, Bt, Gt, g_type> const& mat) :
		m_matrix(
			std::make_pair(
				mat,
				evdm::make_Distribution<T>(mat.Grid, false, 1/*mat.get_padding()*/)
			)
		) {}

	template <
		typename T, typename Bt, typename Gt, evdm::GridEL_type g_type
	>
	Py_Matrix(
		evdm::GridMatrix<T, Bt, Gt, g_type> mat,
		evdm::Distribution<T, Bt, Gt, g_type> distrib
	) : m_matrix(std::make_pair(std::move(mat), std::move(distrib))) {
		/*std::visit([](auto const& _p) {
			
			if (_p.first.get_padding() != _p.second.get_padding()) {
				throw pybind11::value_error(
					"paddings of matrix and evap vector do not match"
				);
			}
		}, m_matrix);*/
	}

	template <typename T>
	Py_Matrix(Py_EL_Grid const& mGridEL,
		const T* m_data, size_t N,
		size_t _OuterStride, size_t _InnerStride,
		const T* m_data_evap, size_t Nevap, size_t _strideEvap/*,
		int padding = -1*/) :
		m_matrix(
			std::visit([&](auto const& _grid)->Matrix_Variant_t {
				size_t grid_size = _grid.Grid->size();

	size_t mat_padding = (N >= grid_size) ?
		evdm::min_deg_2(N - grid_size) : 1;


	size_t vec_padding = (Nevap >= grid_size) ?
		evdm::min_deg_2(Nevap - grid_size) : 1;
	/*if (padding > 0) {
		mat_padding = padding;
		vec_padding = padding;
	}*/

	if (Nevap != 0 && m_data_evap != nullptr) {
		return std::make_pair(
			evdm::make_Matrix<T>(_grid, m_data, grid_size,
				_OuterStride, _InnerStride, 1/*mat_padding*/),
			evdm::make_Distribution_data<T>(
				_grid, m_data_evap, _strideEvap, grid_size, 1/*vec_padding*/
				)
		);
	}
	else {
		return std::make_pair(
			evdm::make_Matrix<T>(_grid, m_data, N,
				_OuterStride, _InnerStride, 1/*mat_padding*/),
			evdm::make_Distribution<T>(
				_grid, false, 1/*vec_padding*/
				)
		);
	}
				}, mGridEL.m_grid)
		) {}

	static Py_Matrix fromArray_impl(
		Py_EL_Grid const& mGridEL,
		pybind11::array const& Mat,
		pybind11::array const& Evap/*,
		int padding*/
	);
	static Py_Matrix fromArray1(
		Py_EL_Grid const& mGridEL,
		pybind11::array const& Mat,
		pybind11::array const& Evap
	);
	static Py_Matrix fromArray2
	(Py_EL_Grid const& mGridEL,
		pybind11::array const& Mat);
	static Py_Matrix from_dict(pybind11::dict const&);

	size_t get_padding() const;


	Py_EL_Grid getGrid()const;

	template <typename T>
	Py_Matrix as_type_t() const {
		return std::visit([&](auto const& mat) {
			return Py_Matrix(mat.first.as_type<T>(), mat.second.as_type<T>());
			}, m_matrix);
	}

	Py_Matrix copy() const;
	void calc_diag();
	Py_Matrix evolve_matrix() const;
	/// @brief adds two matrix if events not intersects
	void combine(Py_Matrix const& _another);

	/// @brief same as combiine
	Py_Matrix& add(Py_Matrix const& _another);
	/// @brief same as add but with copy
	Py_Matrix operator_plus(Py_Matrix const& _another) const;

	std::vector<scatter_event_info> const& get_events()const;
	std::string repr()const;
	Py_Matrix as_type(const char* type_n)const;

	//Py_Distribution total_probs()const;
	Py_Distribution evap_distrib();


	pybind11::array get_matrix(
		pybind11::handle self, int p_in, int p_out, bool raw
	);
	pybind11::dict get_object(pybind11::handle self);
	static Py_Matrix from_object(Py_EL_Grid const& grid, pybind11::dict const&);

	pybind11::array get_diag(
		pybind11::handle self, int p_in, int p_out, bool raw
	);
	pybind11::array get_evap(
		pybind11::handle self, int p_type, bool raw
	);

	Py_Distribution get_diag_distrib();


	static void add_to_python_module(pybind11::module_& m);

};


