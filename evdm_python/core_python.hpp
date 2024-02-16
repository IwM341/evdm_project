#ifndef CORE_PYTHON_HPP
#define CORE_PYTHON_HPP
#include <evdm/core.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <tuple>
#include <set>
template <typename T>
const char* type_name();

template <typename T, typename Tuple_t>
struct tuple_contains;

template <typename T, typename... Args>
struct tuple_contains<T, std::tuple<Args...>> :
	std::disjunction<std::is_same<T, Args>...> {};

template <typename IndexSequenceType,typename IndexSequenceType::value_type I>
struct index_sequence_contain;

template <typename _T,_T...Is, _T I>
struct index_sequence_contain<std::integer_sequence<_T, Is...>, I> :
	std::disjunction<std::integral_constant<bool,(Is==I)>...> {};


#ifdef BODY_MODEL_USE_FLOAT
	#ifdef BODY_MODEL_USE_DOUBLE
		using body_types = std::tuple<float, double>;
	#else
		using body_types = std::tuple < float > ;
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

#ifdef GRID_EL_USE_CUU
#ifdef GRID_EL_USE_CVV
using grid_pos_types = 
	std::integer_sequence<
		evdm::GridEL_type, 
		evdm::GridEL_type::GridCUU,
		evdm::GridEL_type::GridCVV
	>;
#else
using grid_pos_types =
	std::integer_sequence<
		evdm::GridEL_type,
		evdm::GridEL_type::GridCUU
	>;
#endif
#elif defined(GRID_EL_USE_CVV)
using grid_pos_types =
	std::integer_sequence<
		evdm::GridEL_type,
		evdm::GridEL_type::GridCVV
	>;
#else
static_assert("el grid should be lest of 1 type: CUU or CVV");
#endif

struct Py_BodyModel {
	std::variant<
		evdm::BodyModel<float>,
		evdm::BodyModel<double>
	> m_body;

	template <typename T>
	Py_BodyModel(evdm::BodyModel<T> m_body) :m_body(std::move(m_body)) {}
	Py_BodyModel(pybind11::array_t<double> const & RhoValues, double Velocity,
		pybind11::handle Temp = pybind11::none());
	Py_BodyModel(pybind11::array_t<float> const& RhoValues,float Velocity, 
		pybind11::handle Temp = pybind11::none());
	Py_BodyModel(pybind11::handle  const& RhoObject,
				std::string_view dtype,
				std::optional<size_t> _size,
				double Velocity, pybind11::handle Temp = pybind11::none());

	void setTemp(pybind11::handle mTemp);


	const char* dtype()const;
	std::string repr() const;
	size_t size()const;
	pybind11::tuple getRho()const;
	pybind11::tuple getPhi()const;
	pybind11::tuple getM()const;
	pybind11::tuple getQ()const;

	static void add_to_python_module(pybind11::module_& m);
};

struct Py_EL_Grid {
	std::variant<
		//           B vtype, GEL vtype
		evdm::EL_Grid<float,float, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<float, double, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<double, float, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<double, double, evdm::GridEL_type::GridCUU>,
		evdm::EL_Grid<float, float, evdm::GridEL_type::GridCVV>,
		evdm::EL_Grid<float, double, evdm::GridEL_type::GridCVV>,
		evdm::EL_Grid<double, float, evdm::GridEL_type::GridCVV>,
		evdm::EL_Grid<double, double, evdm::GridEL_type::GridCVV>
	> m_grid;
	
	size_t size()const;
	size_t ptypes()const;
	const char* dtype()const;
	std::string repr() const;
	std::string printd1() const;
	
	template <typename T> struct m_type_marker 
	{
		typedef T type;
	};

	template <typename Bt,typename Gt, evdm::GridEL_type _m_type>
	Py_EL_Grid(evdm::EL_Grid<Bt, Gt, _m_type> const& _m_grid) :m_grid(_m_grid) {}
	template <typename T,typename FuncType_NL,typename GEL_t, typename RhoE_ft, typename RhoL_ft>
	inline Py_EL_Grid(evdm::BodyModel<T> const& BM,
		size_t ptypes, size_t Ne,
		FuncType_NL&& Nl_func,
		RhoE_ft && RhoE_func_or_false,
		RhoL_ft && RhoL_func_or_false,
		m_type_marker<GEL_t>,
		std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>):
			m_grid(
				evdm::EL_Grid <T, GEL_t, evdm::GridEL_type::GridCUU>(BM,
					evdm::Grid_types<GEL_t, evdm::GridEL_type::GridCUU>::construct(
						ptypes, (T) - BM->Phi[0], Ne, Nl_func
					),
					2 * Ne
				)
			){}
	template <typename T, typename FuncType_NL, typename GEL_t, 
				typename RhoE_ft, typename RhoL_ft>
	inline Py_EL_Grid(evdm::BodyModel<T> const& BM,
		size_t ptypes, size_t Ne,
		FuncType_NL&& Nl_func,
		RhoE_ft&& RhoE_func_or_false,
		RhoL_ft&& RhoL_func_or_false,
		m_type_marker<GEL_t>,
		std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>) :
			m_grid(
				evdm::EL_Grid <T, GEL_t, evdm::GridEL_type::GridCVV>(BM,
					evdm::Grid_types <GEL_t, evdm::GridEL_type::GridCVV>::construct(
						ptypes, (T) - BM->Phi[0], Ne, RhoE_func_or_false,
						Nl_func, RhoL_func_or_false
					),2*Ne
				)
			) {}

	Py_BodyModel getBody()const;
	pybind11::array getPlot(bool is_internal) const;
	pybind11::tuple getLE() const;
	static void add_to_python_module(pybind11::module_& m);
};

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


struct Py_Distribution {
	DistribVariants m_distrib;

	template <typename T,typename Initializer_t>
	Py_Distribution(Py_EL_Grid const& mGridEL,
		Py_EL_Grid::m_type_marker<T> _type,
		Initializer_t&& init) :
		m_distrib(
			std::visit([&init](auto const& _grid)->DistribVariants {
				return evdm::make_Distribution<T>(_grid, init);
			}, mGridEL.m_grid)
		)
	{}

	template <typename T,typename Bt,typename Gt, evdm::GridEL_type g_type>
	Py_Distribution(
		evdm::Distribution<T, Bt, Gt, g_type> const& distrib) :
		m_distrib(distrib) {}

	Py_EL_Grid getGrid()const;
	
	template <typename T>
	Py_Distribution as_type_t() const{
		return std::visit([](auto const& distrib) {
			return Py_Distribution(distrib.as_type<T>());
		}, m_distrib);
	}
	Py_Distribution copy() const;
	
	std::string repr()const;
	Py_Distribution as_type(const char* type_n)const;
	double count(int ptype = -1)const;
	pybind11::tuple plot(size_t ptype)const;
	static void add_to_python_module(pybind11::module_& m);
};
Py_Distribution CreatePyDistrib(
	Py_EL_Grid const& mGridEL,
	const char* dtype,
	pybind11::handle Init
);
double compare_distribs(
	Py_Distribution const& D1,
	Py_Distribution const& D2,
	std::string const& measure,
	double p_deg,
	int ptype);

struct Py_DistribMeasure {
	double p_deg;
	int ptype;
	std::string measure;
	Py_DistribMeasure(std::string measure, double p_deg, int ptype = -1);
	double call(Py_Distribution const& D1, Py_Distribution const& D2);
	std::string repr()const;
};

struct scatter_event_info {
	std::string name;
	size_t ptype;
	double amount;
	double Nmk;
	bool unique;
	inline scatter_event_info():name(""), ptype(0),amount(0), Nmk(0), unique(false) {}
	inline scatter_event_info(std::string name, size_t ptype,double amount, double Nmk):
		name(name), ptype(ptype), amount(amount), Nmk(Nmk){}
	inline bool operator == (const scatter_event_info& event) const {
		return name == event.name && ptype == event.ptype;
	}
	inline bool operator < (const scatter_event_info& event) const {
		return (name < event.name) || (ptype < event.ptype);
	}
	inline bool operator > (const scatter_event_info& event) const {
		return event < *this;
	}
	inline bool operator <= (const scatter_event_info& event) const {
		return (event < *this) || (event == *this);
	}
	inline bool operator >= (const scatter_event_info& event) const {
		return event <= *this;
	}
	friend inline scatter_event_info linear_combination(
		const scatter_event_info & E1, 
		const scatter_event_info & E2,
		double a1,double a2) {
		return scatter_event_info(
			E1.name,
			E1.ptype,
			a1 * E1.amount + a2 * E2.amount,
			1 / (a1 * a1 / E1.Nmk + a2 * a2 / E2.Nmk)
		);
	}
	std::string repr()const;
	static void add_to_python_module(pybind11::module_& m);
};
struct Py_Capture : public Py_Distribution {
	using Py_Distribution::Py_Distribution;
	inline Py_Capture(Py_Distribution const& _distrib) :
		Py_Distribution(_distrib) {}

	std::vector<scatter_event_info> events;
	static void add_to_python_module(pybind11::module_& m);
	inline std::vector<scatter_event_info>const& get_events()const { return events; }
	


	Py_Capture copy() const; 

	template <typename T>
	inline Py_Capture as_type_t() const {
		return Py_Capture(Py_Distribution::as_type_t<T>());
	}

	inline Py_Capture as_type(const char* dtype) const;
	
	/// <summary>
	/// summ two captures WITH DIFFERENT EVENTS
	/// </summary>
	/// <param name="_another"></param>
	Py_Capture & add(Py_Capture const& _another);
	Py_Capture operator_plus(Py_Capture const& _another) const;
};
struct Py_Matrix{
	MatrixVariants m_matrix;
	std::vector<scatter_event_info> events;

	template <typename T>
	Py_Matrix(Py_EL_Grid const& mGridEL,
		Py_EL_Grid::m_type_marker<T> _type):
		m_matrix(
			std::visit([](auto const& _grid)->MatrixVariants {
				return evdm::make_Matrix<T>(_grid);
			}, mGridEL.m_grid)
		)
	{}

	template <typename T, typename Bt, typename Gt, evdm::GridEL_type g_type>
	Py_Matrix(evdm::GridMatrix<T, Bt, Gt, g_type> const& mat) :
		m_matrix(mat) {}

	Py_EL_Grid getGrid()const;

	template <typename T>
	Py_Matrix as_type_t() const {
		return std::visit([](auto const& mat) {
				return Py_Matrix(mat.as_type<T>());
			}, m_matrix);
	}

	Py_Matrix make_copy() const { static_assert(""); }

	void combine(Py_Matrix const& _another) { static_assert(""); }

	std::vector<scatter_event_info> const& get_events()const;
	std::string repr()const;
	Py_Matrix as_type(const char* type_n)const;

	Py_Distribution total_probs()const;

	
	
	static void add_to_python_module(pybind11::module_& m);

};
#endif//CORE_PYTHON_HPP