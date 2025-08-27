#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <evdm/utils/discret_generator.hpp>
#include <evdm/utils/variant_tools.hpp>
#include "grob_python.hpp"
#include <evdm/utils/prng.hpp>

typedef std::variant<
	evdm::MarkovProcess<float>,
	evdm::MarkovProcess<double>
> MarkovVariant;

struct PyMarkovChain {
	std::variant<
		evdm::MarkovProcess<float>,
		evdm::MarkovProcess<double>
	> m_process;
	static PyMarkovChain Consrtuct(pybind11::handle m_matrix);

	pybind11::array evolute(
		pybind11::array initials,
		double maxT, int Ntraj, int maxStep
	);

};

void PyMarkov_add_to_python_module(pybind11::module& m) {
	namespace py = pybind11;
	auto markov = py::class_<PyMarkovChain>(m, "MarkovChain",
		"class for model markov chain");
	markov.def(
		py::init(&PyMarkovChain::Consrtuct),
		"constructor\n\n"
		"Parameters:\n"
		"___________\n"
		"matrix : array or csc sparse matrix\n\t"
		"probability matrix a[j,i] from i to j",
		py::arg("matrix")
	);
	markov.def(
		"evolute",
		&PyMarkovChain::evolute,
		"evolution function. Returns final distribution \n\n"
		"Parameters:\n"
		"___________\n"
		"init : array or csc sparse matrix\n\t"
		"initial distribution\n",
		"t_max : float\n\t"
		"final time of evolution\n",
		"traj_num : int\n\t"
		"number of trajectories\n",
		"max_step : int\n\t"
		"maximum iteration value\n",

		py::arg("init"), py::arg("t_max"),
		py::arg("traj_num"), py::arg("max_step")
	);

}



template <typename T>
using MatrixMap = Eigen::Map <
	Eigen::MatrixX<T>, 0,
	Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
>;

typedef std::variant<
	Eigen::Map<Eigen::SparseMatrix<float>>,
	Eigen::Map<Eigen::SparseMatrix<double>>,
	MatrixMap<float>,
	MatrixMap<double>
> MatrixTypes;

std::variant<
	pybind11::array,
	pybind11::handle
> array_detect(pybind11::handle m_handle) {

	namespace py = pybind11;

	if (py::isinstance<py::array>(m_handle))// array
	{
		return m_handle.cast<py::array>();
	}
	try {
		// Импортируем scipy.sparse
		py::module scipy_sparse = py::module::import("scipy.sparse");
		// Получаем класс csc_matrix
		py::object csc_matrix_class = scipy_sparse.attr("csc_matrix");

		// Проверяем принадлежность
		if (py::isinstance(m_handle, csc_matrix_class)) {
			return m_handle;
		}
	}
	catch (const std::exception& e) {}
	throw py::type_error(
		"expect numpy array or scipy sparse matrix"
	);
}

template <typename T>
Eigen::Map<Eigen::SparseMatrix<T>>
ConvertSpMatrix(pybind11::handle csc_matrix) {
	auto data = csc_matrix.attr("data").cast<pybind11::array_t<T>>();
	auto indices = csc_matrix.attr("indices").cast<pybind11::array_t<int>>();
	auto indptr = csc_matrix.attr("indptr").cast<pybind11::array_t<int>>();
	auto shape = csc_matrix.attr("shape").cast<std::tuple<int, int>>();

	auto data_v = csc_matrix.attr("data").cast<pybind11::array>();
	auto indices_v = csc_matrix.attr("indices").cast<pybind11::array>();
	auto indptr_v = csc_matrix.attr("indptr").cast<pybind11::array>();

	//pybind11::print("data array ", data_v.data(), ", vs ", (void*)data.data());
	//pybind11::print("indices array ", indices_v.data(), ", vs ", (void*)indices.data());
	//pybind11::print("indptr array ", indptr_v.data(), ", vs ", (void*)indptr.data());

	return Eigen::Map<Eigen::SparseMatrix<T>>(
		std::get<0>(shape), std::get<1>(shape), data.size(),
		const_cast<int*>(indptr.data()),
		const_cast<int*>(indices.data()),
		const_cast<T*>(data.data())
		);
}

template <typename T>
MatrixMap<T>
ConvertMatrix(pybind11::array_t<T> arr) {
	auto buf = arr.request();
	if (buf.ndim != 2) {
		throw std::runtime_error("Number of dimensions must be 2");
	}
	if (!buf.ptr) {
		throw std::runtime_error("Buffer pointer is null");
	}
	Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> stride(
		buf.strides[1] / sizeof(T),   // stride между столбцами
		buf.strides[0] / sizeof(T)  // stride между строками
	);

	return MatrixMap<T>(
		static_cast<T*>(buf.ptr), buf.shape[0], buf.shape[1], stride
	);
}

struct matrix_caster {
	template <typename T>
	using idt = std::type_identity<T>;

	MarkovVariant operator()(pybind11::array m_arr) {
		
		auto m_type = getType([&m_arr]<class T>(idt<T>) {
			//bool b = m_arr.dtype().equal(pybind11::dtype::of<T>());
			return m_arr.dtype().equal(pybind11::dtype::of<T>());
		},  idt<std::variant<float, double>>{}, 
			idt<void>{}
		);
		return std::visit(
			[&m_arr]<class T>(idt<T>)->MarkovVariant
		{
			auto m_mat = ConvertMatrix(pybind11::array_t<T>(m_arr));
			// m_arr.cast<Eigen::Map<Eigen::MatrixX<T>>>();
			return evdm::MarkovProcess<T>(m_mat);
		},m_type);
	}
	MarkovVariant operator()(pybind11::handle m_arr) {
		pybind11::handle dtype = m_arr.attr("dtype");

		auto m_type = getType([dtype]<class T>(idt<T>) 
		{
			return dtype.equal(pybind11::dtype::of<T>());
		}, 
			idt<std::variant<float, double>>{}, 
			idt<void>{}
		);
		return std::visit(
			[m_arr]<class T>(idt<T>)->MarkovVariant
		{
			auto m_mat = ConvertSpMatrix<T>(m_arr);
			return evdm::MarkovProcess<T>(m_mat);
		}, m_type);
	}
};

PyMarkovChain PyMarkovChain::Consrtuct(pybind11::handle m_matrix) {
	matrix_caster m;
	return PyMarkovChain{ std::visit(
		m, array_detect(m_matrix)
	)};
}

pybind11::array PyMarkovChain::evolute(
	pybind11::array initials,
	double maxT, int Ntraj, int maxStep
) {

	return std::visit([&]<class T>(
		evdm::MarkovProcess<T> const & _this
	)->pybind11::array {
		auto Init = initials.cast<
			Eigen::VectorX<T>
		>();

		evdm::VectorGenerator<T> VG(Init.begin(), Init.end());
		evdm::xorshift<T, evdm::genpol::E0I1> G;
		
		Eigen::VectorX<T> FinalDistribs;
		{
			pybind11::gil_scoped_release m_unlock;
			FinalDistribs =
				evdm::MarkovEvolute(_this, VG, (T)maxT, Ntraj, maxStep, G);
		}	
		return make_py_array(FinalDistribs);
		
		
	},m_process);
}