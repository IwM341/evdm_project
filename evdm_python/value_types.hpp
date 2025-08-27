#pragma once

#include <evdm/core/core_body.hpp>
#include <evdm/core/core_grid.hpp>
#include <evdm/core/core_distrib.hpp>
#include <evdm/core/core_matrix.hpp>
#include <evdm/core/core_annihilation.hpp>
#include <evdm/dynamics/dynamics.hpp>
#include <variant>

//#define BODY_MODEL_USE_FLOAT
//#define BODY_MODEL_USE_DOUBLE
//#define GRID_EL_USE_FLOAT
//#define GRID_EL_USE_DOUBLE
//#define GRID_EL_USE_CUU
//#define GRID_EL_USE_CVV
//#define DISTRIB_USE_FLOAT
//#define DISTRIB_USE_DOUBLE


#ifdef BODY_MODEL_USE_FLOAT
	#ifdef BODY_MODEL_USE_DOUBLE
		#define BODY_TYPE_LIST float, double
	#else
		#define BODY_TYPE_LIST float
	#endif
#else
	#ifdef BODY_MODEL_USE_DOUBLE
		#define BODY_TYPE_LIST double
	#else
		#error "BODY_MODEL_USE_FLOAT or BODY_MODEL_USE_DOUBLE should be defined"
	#endif
#endif

#ifdef GRID_EL_USE_FLOAT
	#ifdef GRID_EL_USE_DOUBLE
		#define GRID_EL_TYPE_LIST float, double
	#else
		#define GRID_EL_TYPE_LIST float
	#endif
#else
	#ifdef GRID_EL_USE_DOUBLE
		#define GRID_EL_TYPE_LIST double
	#else
		#error "GRID_EL_USE_FLOAT or GRID_EL_USE_DOUBLE should be defined"
	#endif
#endif
#ifdef DISTRIB_USE_FLOAT
	#ifdef DISTRIB_USE_DOUBLE
		#define DISTRIB_TYPE_LIST float, double
	#else
		#define DISTRIB_TYPE_LIST float
	#endif
#else
	#ifdef DISTRIB_USE_DOUBLE
		#define DISTRIB_TYPE_LIST double
	#else
		#error "DISTRIB_USE_FLOAT or DISTRIB_USE_DOUBLE should be defined"
	#endif
#endif

using GridCUU_t = evdm::GridEL_type_t< evdm::GridEL_type::GridCUU>;
using GridCVV_t = evdm::GridEL_type_t< evdm::GridEL_type::GridCVV>;


#ifdef GRID_EL_USE_CUU
	#ifdef GRID_EL_USE_CVV
		#define GRID_KIND_TYPE_LIST GridCUU_t, GridCVV_t
	#else
		#define GRID_KIND_TYPE_LIST  GridCUU_t
	#endif
#else
	#ifdef GRID_EL_USE_CVV
		#define GRID_KIND_TYPE_LIST  GridCVV_t
	#else
		#error "GRID_EL_USE_CUU  or GRID_EL_USE_CVV should be defined"
	#endif
#endif

using body_types = std::tuple<BODY_TYPE_LIST>;
using grid_types = std::tuple<GRID_EL_TYPE_LIST>;
using distrib_types = std::tuple <DISTRIB_TYPE_LIST>;
using grid_pos_types =
	std::tuple<GRID_KIND_TYPE_LIST>;



template <
	typename Mat_vt, typename Body_vt,
	typename Grid_vt, typename mGrid_type
> using Matrix_Pair =
std::pair<
	evdm::GridMatrix<Mat_vt, Body_vt, Grid_vt, mGrid_type::value>,
	evdm::Distribution<Mat_vt, Body_vt, Grid_vt, mGrid_type::value>
>;



template<typename ... input_t>
using tuple_cat_t =
decltype(std::tuple_cat(
	std::declval<input_t>()...
));


template <typename...Tuples_t>
struct decart_product_tuples;

template <typename...Tuples_t>
using decart_product_tuples_t = typename decart_product_tuples< Tuples_t...>::type;

template <typename _First_Tp, typename existing_product_t>
struct decart_product_tuples_helper;

template <typename _Fist_Tp, typename existing_product_instance_t>
struct decart_product_tuples_helper_one_inst;


template <typename T, typename _Tp>
struct tuple_insert_first;

template <typename T, typename...Args>
struct tuple_insert_first < T, std::tuple<Args...>> {
	typedef std::tuple<T, Args...> type;
};


template <typename...Args, typename existing_product_instance_t>
struct decart_product_tuples_helper_one_inst
	<std::tuple<Args...>, existing_product_instance_t>
{
	typedef std::tuple<
		typename tuple_insert_first<Args, existing_product_instance_t>::type...
	> type;
};

template <typename _First_Tp, typename...product_tuples_t>
struct decart_product_tuples_helper<_First_Tp, std::tuple<product_tuples_t...>> {
	typedef tuple_cat_t<
		typename decart_product_tuples_helper_one_inst<
		_First_Tp, product_tuples_t
		>::type...
	> type;
};

template <typename...Args>
struct decart_product_tuples<std::tuple<Args...>> {
	typedef std::tuple<std::tuple<Args>...> type;
};
template <typename Tp1, typename...Tps>
struct decart_product_tuples<Tp1, Tps...> {
	typedef typename decart_product_tuples_helper<
		Tp1, decart_product_tuples_t<Tps...>
	>::type type;
};


using BodtyModel_tp_t = decart_product_tuples_t<body_types>;
using grid_el_tp_t = decart_product_tuples_t<body_types,grid_types,grid_pos_types>;
using distrib_tp_t = decart_product_tuples_t<distrib_types,body_types, grid_types, grid_pos_types>;



template <typename _Tp>
struct convert_body_tp;

template <typename T>
struct convert_body_tp<std::tuple<T>> {
	using type = evdm::BodyModel<T>;
};

template <typename _Tp>
struct convert_grid_tp;

template <typename body_t,typename grid_vt, typename grid_type_t>
struct convert_grid_tp<std::tuple<body_t, grid_vt, grid_type_t>> {
	using type = evdm::EL_Grid<body_t, grid_vt, grid_type_t::value>;
};

template <typename _Tp>
struct convert_distrib_tp;

template <typename distrib_t, typename body_t, typename grid_vt, typename grid_type_t>
struct convert_distrib_tp<std::tuple<distrib_t,body_t, grid_vt, grid_type_t>> {
	using type = evdm::Distribution<distrib_t,body_t, grid_vt, grid_type_t::value>;
};

template <typename _Tp>
struct convert_matrix_tp;

template <typename distrib_t, typename body_t, typename grid_vt, typename grid_type_t>
struct convert_matrix_tp<std::tuple<distrib_t,body_t, grid_vt, grid_type_t>> {
	using type = evdm::GridMatrix<distrib_t, body_t, grid_vt, grid_type_t::value>;
};

template <typename _Tp>
struct convert_matrix_pair_tp;

template <typename distrib_t, typename body_t, typename grid_vt, typename grid_type_t>
struct convert_matrix_pair_tp<std::tuple<distrib_t, body_t, grid_vt, grid_type_t>> {
	using type = Matrix_Pair<distrib_t, body_t, grid_vt, grid_type_t>;
};

template <typename _Tp>
struct convert_matrix_ann_tp;

template <typename distrib_t, typename body_t, typename grid_vt, typename grid_type_t>
struct convert_matrix_ann_tp<std::tuple<distrib_t, body_t, grid_vt, grid_type_t>> {
	using type = evdm::GridAnnPreMatrix<distrib_t, body_t, grid_vt, grid_type_t::value>;
};

template <template <typename> typename Converter_tp_t,typename tuple_pf_tuples>
struct variant_type_s;

template <template <typename> typename Converter_tp_t, typename...Tuples>
struct variant_type_s<Converter_tp_t, std::tuple<Tuples...>> {
	using type = std::variant<typename Converter_tp_t<Tuples>::type...>;
};

template <template <typename> typename Converter_tp_t, typename tuple_pf_tuples>
using variant_type_t = typename variant_type_s<Converter_tp_t, tuple_pf_tuples>::type;

using BodyModel_Variant_t = 
	//std::variant< evdm::BodyModel<float>>;
variant_type_t<convert_body_tp, BodtyModel_tp_t>;


#define _GRID_EL_TMPL_ class _T1,class _T2,evdm::GridEL_type _m_grid_t
#define _GRID_EL_PARS_ _T1,_T2,_m_grid_t

#define _DISTRIB_TMPL_ class _T1,class _T2,class _T3,evdm::GridEL_type _m_grid_t
#define _DISTRIB_PARS_ _T1,_T2,_T3,_m_grid_t

#define _MATRIX_TMPL_ class _T1,class _T2,class _T3,evdm::GridEL_type _m_grid_t
#define _ANN_MATRIX_TMPL_ class _T1a,class _T2a,class _T3a,evdm::GridEL_type _ma_grid_t
#define _MATRIX_PARS_ _T1,_T2,_T3,_m_grid_t
#define _ANN_MATRIX_PARS_ _T1a,_T2a,_T3a,_ma_grid_t



using ELGrid_Variant_t =
	/*std::variant<
		evdm::EL_Grid<float, float, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<float, float, evdm::GridEL_type::GridCVV>>;*/
variant_type_t<convert_grid_tp, grid_el_tp_t>;
using Distrib_Variant_t =
/*std::variant<
	evdm::Distribution<float, float, float, evdm::GridEL_type::GridCUU>,
	evdm::Distribution<float,float, float, evdm::GridEL_type::GridCVV>>;*/
variant_type_t<convert_distrib_tp, distrib_tp_t>;

typedef evdm::Distribution<float, float, float, evdm::GridEL_type::GridCUU> Distrib_Float_t;

template <
	typename Mat_vt, typename Body_vt, 
	typename Grid_vt, evdm::GridEL_type GT 
> using Matrix_Pair_Inst =
std::pair<
	evdm::GridMatrix<Mat_vt, Body_vt, Grid_vt, GT>,
	evdm::Distribution<Mat_vt, Body_vt, Grid_vt, GT>
>;



template <typename Matrix_Pair_Inst_t>
struct MatrixTypeInfo;

template <
	typename Mat_vt, typename Body_vt,
	typename Grid_vt, evdm::GridEL_type GT 
>
struct MatrixTypeInfo<Matrix_Pair_Inst<Mat_vt, Body_vt, Grid_vt, GT>> {
	using matrix_vtype = Mat_vt;
	using body_vtype = Body_vt;
	using grid_vtype = Grid_vt;
	static constexpr evdm::GridEL_type grid_type = GT;
};


//variant_type_t<convert_distrib_tp, distrib_tp_t>;
using Matrix_Variant_t =
	/*std::variant<
	Matrix_Pair_Inst<float, float, float, evdm::GridEL_type::GridCUU>,
	Matrix_Pair_Inst<float, float, float, evdm::GridEL_type::GridCVV>
>;*/
variant_type_t<convert_matrix_pair_tp, distrib_tp_t>;

typedef Matrix_Pair_Inst<float, float, float, evdm::GridEL_type::GridCUU> 
	Matrix_Variant_FFF_t;

using PreAnn_Variant_t = 
/*std::variant <
		evdm::GridAnnPreMatrix<float, float, float, evdm::GridEL_type::GridCUU>,
		evdm::GridAnnPreMatrix<float, float, float, evdm::GridEL_type::GridCVV>
	>;*/
variant_type_t<convert_matrix_ann_tp, distrib_tp_t>;

#define SCATTER_COUNT \
	(1 + (defined SCATTER_METHOD_USE_NAIVE) + \
	 (defined SCATTER_METHOD_USE_SOFT) + \
	 (defined SCATTER_METHOD_USE_NOTHERM) + \
	 (defined SCATTER_METHOD_USE_SOFT_TRESH))

using ScatterMethodVariant_t =
	std::variant <
	evdm::ThermGaussGenerator_Full
#if (SCATTER_COUNT > 1)
	,
#endif
#ifdef SCATTER_METHOD_USE_NAIVE
	evdm::ThermGaussGenerator_Naive
#endif //SCATTER_METHOD_USE_NAIVE
#if (SCATTER_COUNT > 2)
	,
#endif 
#ifdef SCATTER_METHOD_USE_SOFT
	evdm::ThermGaussGenerator_Soft8
#endif //SCATTER_METHOD_USE_SOFT

#if (SCATTER_COUNT > 3)
	,
#endif
#ifdef SCATTER_METHOD_USE_NOTHERM
	evdm::ThermGaussGenerator_NoTherm
#endif // SCATTER_METHOD_USE_NOTHERM
#if (SCATTER_COUNT > 4)
	,
#endif
#ifdef SCATTER_METHOD_USE_SOFT_TRESH
	evdm::ThermGaussGenerator_Soft8_Treshold
#endif // SCATTER_METHOD_USE_SOFT_TRESH
>;

#define SCATTER_MEASURE_USE_COUNT \
	( (defined SCATTER_MEASURE_USE_DEDL) + \
	 (defined SCATTER_MEASURE_USE_DEDL2) )

using ScatterMeasureVariant_t =
std::variant<
#ifdef SCATTER_MEASURE_USE_DEDL
	evdm::measure_dEdL
#endif // SCATTER_MEASURE_USE_DEDL
#if (SCATTER_MEASURE_USE_COUNT > 1)
	,
#endif // (SCATTER_MEASURE_USE_COUNT > 1)
#ifdef SCATTER_MEASURE_USE_DEDL2
	evdm::measure_dEdL2
#endif
>;

