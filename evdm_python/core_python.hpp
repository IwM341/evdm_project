#ifndef CORE_PYTHON_HPP
#define CORE_PYTHON_HPP
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <tuple>
#include <set>
#include "value_types.hpp"

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

template <typename T>
bool check_np_type(pybind11::array const& np_array) {
	if constexpr (std::is_floating_point_v<T>) {
		if (np_array.dtype().kind() == 'f') {
			return np_array.dtype().itemsize() == sizeof(T);
		} else {
			return false;
		}
	} if constexpr (std::is_integral_v<T>) {
		if (np_array.dtype().kind() == 'i') {
			return np_array.dtype().itemsize() == sizeof(T);
		}
		else {
			return false;
		}
	}
}

struct Py_BodyModel {
	BodyModel_Variant_t	m_body;
#ifdef BODY_MODEL_USE_DOUBLE
	using max_vtype = double;
#else
	using max_vtype = float;
#endif // DISTRIB_USE_DOUBLE
#ifdef BODY_MODEL_USE_FLOAT
	using min_vtype = float;
#else
	using min_vtype = double;
#endif // DISTRIB_USE_DOUBLE

	template <typename T>
	Py_BodyModel(evdm::BodyModel<T> m_body) :m_body(std::move(m_body)) {}

#ifdef BODY_MODEL_USE_DOUBLE
	Py_BodyModel(pybind11::array_t<double> const & RhoValues, double Velocity,
		pybind11::handle Temp = pybind11::none());
#endif
#ifdef BODY_MODEL_USE_FLOAT
	Py_BodyModel(pybind11::array_t<float> const& RhoValues,float Velocity, 
		pybind11::handle Temp = pybind11::none());
#endif
	Py_BodyModel(pybind11::handle  const& RhoObject, double Velocity,
				std::string_view dtype,
				pybind11::handle Temp = pybind11::none());

	void setTemp(pybind11::handle mTemp);


	const char* dtype()const;
	std::string repr() const;
	size_t size()const;
	pybind11::tuple getRho()const;
	pybind11::tuple getPhi()const;
	pybind11::tuple getM()const;
	pybind11::tuple getQ()const;


	pybind11::dict get_object(pybind11::handle self) const;
	static Py_BodyModel from_dict(pybind11::dict const& Object);

	static void add_to_python_module(pybind11::module_& m);

};

struct Py_EL_Grid {
	ELGrid_Variant_t m_grid;
	
	/// @brief size of inner part of grid
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
		std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCUU>,
		evdm::TrajPoolInitParams_t TrajPoolInit = evdm::TrajPoolInitParams_t()):
			m_grid(
				evdm::EL_Grid <T, GEL_t, evdm::GridEL_type::GridCUU>(BM,
					evdm::Grid_types<GEL_t, evdm::GridEL_type::GridCUU>::construct(
						ptypes, (T) - BM->Phi[0], Ne, Nl_func
					),
					2 * Ne, TrajPoolInit
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
		std::integral_constant<evdm::GridEL_type, evdm::GridEL_type::GridCVV>,
		evdm::TrajPoolInitParams_t TrajPoolInit = evdm::TrajPoolInitParams_t()) :
			m_grid(
				evdm::EL_Grid <T, GEL_t, evdm::GridEL_type::GridCVV>(BM,
					evdm::Grid_types <GEL_t, evdm::GridEL_type::GridCVV>::construct(
						ptypes, (T) - BM->Phi[0], Ne, RhoE_func_or_false,
						Nl_func, RhoL_func_or_false
					),4*Ne, TrajPoolInit
				)
			) {}
	Py_EL_Grid refine(size_t Ne, size_t Nl) const;
	
	static Py_EL_Grid from_dict(Py_BodyModel BM, pybind11::dict const& grid_dict);
	static Py_EL_Grid from_dict1(pybind11::dict const& grid_dict);
	
	pybind11::dict get_object(pybind11::handle self);

	Py_BodyModel getBody()const;
	pybind11::array getPlot(bool is_internal) const;
	pybind11::tuple getLE() const;

	pybind11::tuple getRminRmax(double e, double l_undim) const;

	pybind11::list getTrajTFuncs(pybind11::handle output_param_name, bool LE_mult = false) const;
	std::function<double(double,double)> getTraj_callFunc(pybind11::handle output_param_name) const;
	std::function<double(double)> getLE_callFunc() const;

	static void add_to_python_module(pybind11::module_& m);
	std::variant<size_t,std::pair<size_t,size_t>> get_index(
		double e, double L,bool linear,bool hidden)const;
	size_t get_e_index(double e)const;
	pybind11::dict get_el(
		std::variant<size_t, std::pair<size_t, size_t>> index)const;


	pybind11::array get_E_array(pybind11::handle self)const;
	pybind11::array get_L_array(
		pybind11::handle,size_t index, bool hidden
	)const;



	pybind11::dict _get_traj_all(double e, double L);

	pybind11::dict _get_traj_all_inter(double e, double L);
};



struct Py_Distribution {
	Distrib_Variant_t m_distrib;

	template <typename T,typename Initializer_t>
	Py_Distribution(Py_EL_Grid const& mGridEL,
		Py_EL_Grid::m_type_marker<T> _type,
		Initializer_t&& init) :
		m_distrib(
			std::visit([&init](auto const& _grid)->Distrib_Variant_t {
				return evdm::make_Distribution<T>(_grid, init);
			}, mGridEL.m_grid)
		)
	{}

	template <typename ArrayType>
	Py_Distribution(Py_EL_Grid const& mGridEL,
		ArrayType && _array,size_t padding) :
		m_distrib(
			std::visit([&_array, padding](auto const& _grid)->Distrib_Variant_t {
				return evdm::make_Distribution_array(
					_grid, std::forward<ArrayType>(_array), padding
				);
			}, mGridEL.m_grid)
		)
	{}


	template <typename T>
	Py_Distribution(Py_EL_Grid const& mGridEL,
		const T * _data, size_t _stride,size_t _size,size_t padding) :
		m_distrib(
			std::visit([&](auto const& _grid)->Distrib_Variant_t {
				return evdm::make_Distribution_data<T>(_grid, _data, _stride,_size,padding);
				}, mGridEL.m_grid)
		)
	{}

	template <typename T,typename Bt,typename Gt, evdm::GridEL_type g_type>
	Py_Distribution(
		evdm::Distribution<T, Bt, Gt, g_type> const& distrib) :
		m_distrib(distrib) {}

	static Py_Distribution CreatePyDistribFromArray(
		Py_EL_Grid const& mGridEL,
		pybind11::array values
	);
	static Py_Distribution from_dict(pybind11::dict const&);

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
	
	pybind11::dict get_object(pybind11::handle self);

	pybind11::array get_array(
		pybind11::handle self, int ptype,bool raw = false);

	pybind11::tuple get_E_distrib(int ptype) const;
	double get_avarage(pybind11::handle Functor, int ptype) const;

	double count(int ptype = -1)const;
	pybind11::tuple get_r_dense(
		pybind11::handle ptypes,
		double r_min, double  r_max, size_t Nr,
		pybind11::handle Nperbin,
		pybind11::handle opt_dense_funct,
		pybind11::handle update_function) const;

	pybind11::tuple plot2o(size_t ptype)const;
	pybind11::tuple plot1o(
		size_t ptype , 
		std::string_view m_measure,
		std::string_view LE_Space)const;

	static void add_to_python_module(pybind11::module_& m);
	size_t get_padding() const;
};
Py_Distribution CreatePyDistrib(
	Py_EL_Grid const& mGridEL,
	const char* dtype,
	pybind11::handle Init
);
Py_Distribution CreateDistribFromDict(
	Py_EL_Grid const& mGridEL, pybind11::dict X
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
	size_t ptype_in;
	size_t ptype_out;
	double amount;
	double Nmk;
	bool unique;

	scatter_event_info(pybind11::dict const & _object);
	inline scatter_event_info() :
		name(""), ptype_in(0), ptype_out(0), amount(0), Nmk(0), unique(false) {}
	pybind11::dict to_object() const;
	inline scatter_event_info(
		std::string name, size_t ptype_in, size_t ptype_out,
		double amount, double Nmk,
		bool unique = true):
		name(name), ptype_in(ptype_in), ptype_out(ptype_out),
		amount(amount), Nmk(Nmk),unique(unique){}

	inline bool operator == (const scatter_event_info& event) const {
		return name == event.name && ptype_in == event.ptype_in 
			&& ptype_out == event.ptype_out;
	}
	inline bool operator < (const scatter_event_info& event) const {
		return (name < event.name) || (ptype_in < event.ptype_in)
			|| ((ptype_out < event.ptype_out));
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
			E1.ptype_in,
			E1.ptype_out,
			a1 * E1.amount + a2 * E2.amount,
			1 / (a1 * a1 / E1.Nmk + a2 * a2 / E2.Nmk),
			E1.unique || E2.unique
		);
	}
	std::string repr()const;
	static void add_to_python_module(pybind11::module_& m);
};
struct Py_Capture : public Py_Distribution {
	using Py_Distribution::Py_Distribution;
	inline Py_Capture(Py_Distribution _distrib) :
		Py_Distribution(std::move(_distrib)) {}

	static Py_Capture from_dict(pybind11::dict const&);

	std::vector<scatter_event_info> events;
	static void add_to_python_module(pybind11::module_& m);
	inline std::vector<scatter_event_info>const& get_events()const { return events; }
	
	Py_Capture(Py_EL_Grid const&m_grid, 
		pybind11::dict m_py_object);
	pybind11::dict get_object(pybind11::handle self);
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
	Matrix_Variant_t m_matrix;
	//Matrix_Pair_Inst<T_i ...> m_matrix
	std::vector<scatter_event_info> events;

	template <typename T>
	Py_Matrix(Py_EL_Grid const& mGridEL,
		Py_EL_Grid::m_type_marker<T> _type,
		size_t padding = 128
	):
		m_matrix(
			std::visit([padding](auto const& _grid)->Matrix_Variant_t {
				return std::make_pair(
					evdm::make_Matrix<T>(_grid,padding),
					evdm::make_Distribution<T>(_grid,false, padding)
				);
			}, mGridEL.m_grid)
		)
	{}
	
	/// @brief constructor from object
	Py_Matrix(
		Py_EL_Grid const& mGridEL,
		pybind11::dict const& saved_obj);

	/// @brif constructor from GridMatrix
	template <
		typename T, typename Bt, typename Gt, evdm::GridEL_type g_type
	>
	Py_Matrix(evdm::GridMatrix<T, Bt, Gt, g_type> const& mat) :
		m_matrix(
			std::make_pair(
				mat, 
				evdm::make_Distribution<T>(mat.Grid, false,mat.get_padding())
			)
		){}

	template <
		typename T, typename Bt, typename Gt, evdm::GridEL_type g_type
	>
	Py_Matrix(
		evdm::GridMatrix<T, Bt, Gt, g_type> mat,
		evdm::Distribution<T, Bt, Gt, g_type> distrib
	) : m_matrix(std::make_pair(std::move(mat), std::move(distrib))) {
		std::visit([](auto const& _p) {
			if (_p.first.get_padding() != _p.second.get_padding()) {
				throw pybind11::value_error(
					"paddings of matrix and evap vector do not match"
				);
			}
		}, m_matrix);
	}

	template <typename T>
	Py_Matrix(Py_EL_Grid const& mGridEL,
		const T* m_data, size_t N,
		size_t _OuterStride, size_t _InnerStride,
		const T* m_data_evap, size_t Nevap, size_t _strideEvap,
		int padding = -1) :
		m_matrix(
			std::visit([&](auto const& _grid)->Matrix_Variant_t {
				size_t grid_size = _grid.Grid->size();

				size_t mat_padding = (N >= grid_size) ?
					evdm::min_deg_2(N - grid_size) : 1;
				

				size_t vec_padding = (Nevap >= grid_size) ?
					evdm::min_deg_2(Nevap - grid_size) : 1;
				if (padding > 0) {
					mat_padding = padding;
					vec_padding = padding;
				}

				if (Nevap != 0 && m_data_evap != nullptr) {
					return std::make_pair(
						evdm::make_Matrix<T>(_grid, m_data, grid_size,
							_OuterStride, _InnerStride, mat_padding),
						evdm::make_Distribution_data<T>(
							_grid, m_data_evap, _strideEvap, grid_size, vec_padding
						)
					);
				}
				else {
					return std::make_pair(
						evdm::make_Matrix<T>(_grid, m_data, N,
							_OuterStride, _InnerStride, mat_padding),
						evdm::make_Distribution<T>(
							_grid, false, mat_padding
						)
					);
				}
			}, mGridEL.m_grid)
		) {}

	static Py_Matrix fromArray_impl(
		Py_EL_Grid const& mGridEL,
		pybind11::array const& Mat,
		pybind11::array const& Evap,
		int padding
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
				return Py_Matrix(mat.first.as_type<T>(),mat.second.as_type<T>());
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
		pybind11::handle self,int p_in,int p_out,bool raw
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



template <typename Grid_t>
std::vector<size_t> get_N_distrib_from_handle(Grid_t const& Grid, pybind11::handle h) {
	if (pybind11::isinstance<pybind11::array>(h)) {
		auto A = h.cast< pybind11::array_t<size_t>>();
		if (Grid.size() != A.size()) {
			pybind11::index_error("Array of values not match grid size");
		}
		return std::vector<size_t>((size_t*)A.data(), (size_t*)A.data() + A.size());
	}
	if (pybind11::isinstance<pybind11::function>(h)) {
		typedef decltype(Grid.grid()[0].left) T;
		try {
			auto F_el = h.cast < std::function<size_t(T, T)>>();
			(void)(F_el(0, 0));
			std::vector<size_t> NRet(Grid.size());
			auto it = NRet.begin();
			for (auto [de, dl] : Grid) {
				*it = F_el(de.center(), dl.center());
				++it;
			}
			return NRet;
		}
		catch (std::exception& e) {
			auto F_el = h.cast < std::function<size_t(T, T, T, T)>>();
			(void)(F_el(0, 0, 0, 0));
			std::vector<size_t> NRet(Grid.size());
			auto it = NRet.begin();
			for (auto [de, dl] : Grid) {
				auto [e0, e1] = de;
				auto [l0, l1] = dl;
				*it = F_el(e0, e1, l0, l1);
				++it;
			}
			return NRet;
		}
	}
	else {
		pybind11::value_error("expect function or array to make distrib of N");
		return {};
	}
}


#endif//CORE_PYTHON_HPP