#ifndef VALUE_TYPES_HPP
#define VALUE_TYPES_HPP
#include <evdm/core.hpp>
#include <variant>

#define BODY_MODEL_USE_FLOAT
//#define BODY_MODEL_USE_DOUBLE
#define GRID_EL_USE_FLOAT
//#define GRID_EL_USE_DOUBLE
//#define GRID_EL_USE_CUU
#define GRID_EL_USE_CVV
#define DISTRIB_USE_FLOAT
//#define DISTRIB_USE_DOUBLE

#ifdef BODY_MODEL_USE_FLOAT
#ifdef BODY_MODEL_USE_DOUBLE
using body_types = std::tuple<float, double>;
#else
using body_types = std::tuple < float >;
#endif
#elif defined(BODY_MODEL_USE_DOUBLE)
using body_types = std::tuple<double>;
#else
static_assert("Body Model should have at lest 1 type: float or double");
#endif

#ifdef GRID_EL_USE_FLOAT
#ifdef GRID_EL_USE_DOUBLE
using grid_types = std::tuple<float, double>;
#else
using grid_types = std::tuple < float >;
#endif
#elif defined(GRID_EL_USE_DOUBLE)
using grid_types = std::tuple<double>;
#else
static_assert("el grid should have at lest 1 type: float or double");
#endif


#ifdef DISTRIB_USE_FLOAT
#ifdef DISTRIB_USE_DOUBLE
using distrib_types = std::tuple<float, double>;
#else
using distrib_types = std::tuple < float >;
#endif
#elif defined(DISTRIB_USE_DOUBLE)
using distrib_types = std::tuple<double>;
#else
static_assert("Body Model should have at lest 1 type: float or double");
#endif

using GridCUU_t = 
	std::integral_constant< 
		evdm::GridEL_type,
		evdm::GridEL_type::GridCUU
	>;
using GridCVV_t =
std::integral_constant<
	evdm::GridEL_type,
	evdm::GridEL_type::GridCVV
>;

#ifdef GRID_EL_USE_CUU
#ifdef GRID_EL_USE_CVV
using grid_pos_types =
	std::tuple<GridCUU_t,GridCVV_t>;
#else
using grid_pos_types =
	std::tuple<GridCUU_t>;
#endif
#elif defined(GRID_EL_USE_CVV)
using grid_pos_types =
	std::tuple<GridCVV_t>;
#else
static_assert("el grid should be lest of 1 type: CUU or CVV");
#endif

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

template <typename _Fist_Tp,typename existing_product_instance_t>
struct decart_product_tuples_helper_one_inst;


template <typename T,typename _Tp>
struct tuple_insert_first;

template <typename T, typename...Args>
struct tuple_insert_first < T, std::tuple<Args...>>{
	using type = std::tuple<T, Args...>;
};


template <typename...Args, typename existing_product_instance_t>
struct decart_product_tuples_helper_one_inst
	<std::tuple<Args...>, existing_product_instance_t>
{
	using type = std::tuple<
		typename tuple_insert_first<Args, existing_product_instance_t>::type...
	>;
};

template <typename _First_Tp, typename...product_tuples_t>
struct decart_product_tuples_helper<_First_Tp, std::tuple<product_tuples_t...>> {
	using type = tuple_cat_t<
		typename decart_product_tuples_helper_one_inst<
		_First_Tp, product_tuples_t
		>::type
	>;
};

template <typename...Args>
struct decart_product_tuples<std::tuple<Args...>> {
	using type = std::tuple<std::tuple<Args>...>;
};
template <typename Tp1,typename...Tps>
struct decart_product_tuples<Tp1, Tps...> {
	using type = 
		typename decart_product_tuples_helper<
			Tp1, decart_product_tuples_t<Tps...>
		>::type;
};

using BodtyModel_tp_t = body_types;
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

template <template <typename> typename Converter_tp_t,typename tuple_pf_tuples>
struct variant_type_t;

template <template <typename> typename Converter_tp_t, typename...Tuples>
struct variant_type_t<Converter_tp_t, std::tuple<Tuples...>> {
	using type = std::variant<typename Converter_tp_t<Tuples>::type...>;
};

using BodyModel_Variant_t = 
	std::variant< evdm::BodyModel<float>>;
//variant_type_t<convert_body_tp, BodtyModel_tp_t>;


#define _GRID_EL_TMPL_ <class _T1,class _T2,evdm::GridEL_type _m_grid_t>
#define _GRID_EL_PARS_ _T1,_T2,_m_grid_t

#define _DISTRIB_TMPL_ <class _T1,class _T2,class _T3,evdm::GridEL_type _m_grid_t>
#define _DISTRIB_PARS_ _T1,_T2,_T3,_m_grid_t

#define _MATRIX_TMPL_ <class _T1,class _T2,class _T3,evdm::GridEL_type _m_grid_t>
#define _MATRIX_PARS_ _T1,_T2,_T3,_m_grid_t

template <
	typename Mat_vt, typename Body_vt,
	typename Grid_vt, typename mGrid_type
> using Matrix_Pair = 
	std::pair<
		evdm::GridMatrix<Mat_vt, Body_vt, Grid_vt, mGrid_type::value>,
		evdm::Distribution<Mat_vt, Body_vt, Grid_vt, mGrid_type::value>
	>;


using ELGrid_Variant_t =
	std::variant< 
		evdm::EL_Grid<float, float, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<float, float, evdm::GridEL_type::GridCVV>>;
//variant_type_t<convert_grid_tp, grid_el_tp_t>;
using Distrib_Variant_t = 
	std::variant< 
		evdm::Distribution<float, float, float, evdm::GridEL_type::GridCUU>,
		evdm::Distribution<float,float, float, evdm::GridEL_type::GridCVV>>;

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
	std::variant< 
	Matrix_Pair_Inst<float, float, float, evdm::GridEL_type::GridCUU>,
	Matrix_Pair_Inst<float, float, float, evdm::GridEL_type::GridCVV>
>;
// variant_type_t<convert_matrix_tp, distrib_tp_t>;


/*
std::variant<
		//           B vtype, GEL gridtype
		evdm::EL_Grid<float,float, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<float, double, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<double, float, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<double, double, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<float, float, evdm::GridEL_type::GridCVV>,
		evdm::EL_Grid<float, double, evdm::GridEL_type::GridCVV>,
		evdm::EL_Grid<double, float, evdm::GridEL_type::GridCVV>,
		evdm::EL_Grid<double, double, evdm::GridEL_type::GridCVV>
	>
*/

#define DeclareClassVariants(AliasName,evdm_class_name) \
using AliasName = std::variant<\
	evdm::evdm_class_name<float, float, float, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<float, float, double, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<float, double, float, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<float, double, double, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<float, float, float, evdm::GridEL_type::GridCVV>,\
	evdm::evdm_class_name<float, float, double, evdm::GridEL_type::GridCVV>,\
	evdm::evdm_class_name<float, double, float, evdm::GridEL_type::GridCVV>,\
	evdm::evdm_class_name<float, double, double, evdm::GridEL_type::GridCVV>,\
\
	evdm::evdm_class_name<double, float, float, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<double, float, double, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<double, double, float, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<double, double, double, evdm::GridEL_type::GridCUU>,\
	evdm::evdm_class_name<double, float, float, evdm::GridEL_type::GridCVV>,\
	evdm::evdm_class_name<double, float, double, evdm::GridEL_type::GridCVV>,\
	evdm::evdm_class_name<double, double, float, evdm::GridEL_type::GridCVV>,\
	evdm::evdm_class_name<double, double, double, evdm::GridEL_type::GridCVV>\
>;

DeclareClassVariants(DistribVariants, Distribution)
DeclareClassVariants(MatrixVariants, GridMatrix)
#endif//VALUE_TYPES_HPP