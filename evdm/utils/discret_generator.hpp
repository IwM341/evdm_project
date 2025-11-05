#pragma once

#include <vector>
#include <algorithm>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <stdexcept>


namespace evdm {

    struct IndexGenBigHelper {


        template <typename Storage_iterator_t>
        static void  PrepareValues(
            Storage_iterator_t values_begin,
            Storage_iterator_t values_end,
            bool do_sort = true)
        {
            if (do_sort)
                std::sort(values_begin, values_end);
            auto it = values_begin;
            auto jt = it + 1;
            for (; jt != values_end; ++jt) {
                *jt += *it;
                it = jt;
            }
            if (it != values_end) {
                auto max_value = *it;
                auto m_div = 1 / (max_value > 0 ? max_value : decltype(max_value)(1));
                for (auto kt = values_begin; kt != values_end; ++kt) {
                    *kt *= m_div;
                }
            }
        }

        template <typename Storage_iterator_t, typename Input_iterator_t>
        static Storage_iterator_t  MakeValues(
            Storage_iterator_t values_begin,
            Input_iterator_t m_begin,
            Input_iterator_t m_end)
        {
            Storage_iterator_t values_end = std::copy(m_begin, m_end, values_begin);
            PrepareValues(values_begin, values_end);
            return values_end;
        }


        template <typename iter_t, typename Vt>
        inline static size_t gen(iter_t m_begin, iter_t m_end, Vt xi) {
            auto it = std::lower_bound(m_begin, m_end, xi);
            return std::min((size_t)(it - m_begin), (size_t)(m_end - m_begin) - 1);
        }
    };

    template <typename Vt>
    struct VectorGenerator {
        std::vector<size_t> indicies;
        std::vector<Vt> xi_values;
        size_t m_size;

        inline size_t size() const {
            return m_size;
        }
        template <typename iter_t>
        VectorGenerator(iter_t m_begin, iter_t m_end) {
            m_size = (size_t)(m_end - m_begin);
            if(m_size == 0){
                throw std::runtime_error("VectorGenerator: size == 0");
            }
            size_t count_nonzeros = std::count_if(
                m_begin, m_end, [](auto i) { return i > 0; }
            );
            indicies.reserve(count_nonzeros);
            xi_values.reserve(count_nonzeros);


            std::vector<std::pair<Vt, size_t>> data_inds;
            data_inds.reserve(count_nonzeros);
            size_t i = 0;
            for (auto it = m_begin; it != m_end; ++it) {
                if (*it > 0) {
                    data_inds.push_back({ *it,i });
                }
                ++i;
            }
            std::sort(data_inds.begin(), data_inds.end(),
                [](auto const& p, auto const& q) {
                    return p.first < q.first;
                });
            for (auto p : data_inds) {
                indicies.push_back(p.second);
                xi_values.push_back(p.first);
            }
                IndexGenBigHelper::PrepareValues(xi_values.begin(), xi_values.end(), false);
        }
        inline size_t gen(Vt xi)const {
            return indicies[IndexGenBigHelper::gen(xi_values.begin(), xi_values.end(), xi)];
        }
    };

    namespace __utils {
        
        template <typename Mat_t>
        size_t NonzeroColCount(Eigen::DenseBase<Mat_t> const& M, size_t i) {
            size_t sch = 0;
            for (size_t j = 0; j < M.rows(); ++j) {
                if (M(j, i) > 0 && j != i) {
                    sch++;
                }
            }
            return sch;
        }

        template <typename Mat_t>
        size_t NonzeroColCount(Eigen::SparseMatrixBase<Mat_t> const& M, size_t i) {
            if constexpr (Mat_t::IsRowMajor) {
                throw std::runtime_error(
                    "Support only col major matrix in markov process"
                );
            }
            typedef Eigen::SparseMatrixBase<Mat_t> SpM_t;
            typedef typename Mat_t::InnerIterator InIt;
            size_t sch = 0;
            Mat_t const& _M = M.derived();
            for (InIt it(_M, i); it; ++it) {
                size_t j = it.row();
                auto val = it.value();
                if (val > 0 && j != i) {
                    sch++;
                }
            }
            return sch;
        }

        template <typename Mat_t>
        auto NonzeroColElements(
            Eigen::DenseBase<Mat_t> const& M,
            size_t i,
            size_t m_noneros)
        {
            typedef typename Eigen::DenseBase<Mat_t>::Scalar T;
            std::vector<std::pair<T, size_t>> Values;
            Values.reserve(m_noneros);
            for (size_t j = 0; j < M.rows(); ++j) {
                auto val = M(j, i);
                if (val > 0 && j != i) {
                    Values.push_back({ val,j });
                }
            }
            return Values;
        }

        template <typename Mat_t>
        auto NonzeroColElements(
            Eigen::SparseMatrixBase<Mat_t> const& M,
            size_t i,
            size_t m_noneros)
        {
            typedef Eigen::SparseMatrixBase<Mat_t> SpM_t;
            typedef typename Mat_t::InnerIterator InIt;
            typedef typename SpM_t::Scalar T;
            std::vector<std::pair<T, size_t>> Values;
            Mat_t const& _M = M.derived();
            Values.reserve(m_noneros);
            for (InIt it(_M, i); it; ++it) {
                size_t j = it.row();
                auto val = it.value();
                if (val > 0 && j != i) {
                    Values.push_back({ val,j });
                }
            }
            return Values;
        }
    };

    template <typename Vt>
    struct MarkovProcess {

        std::vector<size_t> indicies_data; // (...), (...) -N times
        std::vector<Vt> xi_values_data;

        std::vector<size_t> shift_data;
        std::vector<Vt> probabilities; // size = N
        std::vector<Vt> EvapProbs; // size = N
        size_t m_size;

        inline size_t size() const {
            return m_size;
        }

        template <typename Matrix_t,typename Evap_t = bool>
        MarkovProcess(Matrix_t const& m_mat, Evap_t const & Evap = false) {
            size_t N = m_mat.cols();
            m_size = N;
            if (m_mat.cols() != m_mat.rows()) {
                throw std::runtime_error("MarkovProcessDens: cols != rows");
            }
            probabilities.resize(N);
            EvapProbs.resize(N,0);
            if constexpr (!std::is_same_v< Evap_t, bool>) {
                for (size_t i = 0; i < Evap.size();++i) {
                    EvapProbs[i] = Evap[i];
                }
            }
            
            shift_data.resize(N + 1);
            shift_data[0] = 0;
            for (size_t i = 0; i < N; ++i) {
                size_t NonZero =
                    __utils::NonzeroColCount(m_mat, i);
                shift_data[i + 1] = shift_data[i] + NonZero;
            }
            xi_values_data.resize(shift_data.back());
            indicies_data.resize(shift_data.back());
            for (size_t i = 0; i < N; ++i) {
                size_t NonZero = shift_data[i + 1] - shift_data[i];
                auto data_inds = __utils::NonzeroColElements(m_mat, i, NonZero);
                size_t shift = shift_data[i];
                probabilities[i] = 0;
                if (NonZero) {
                    std::sort(data_inds.begin(), data_inds.end(),
                        [](auto const& p, auto const& q) {
                            return p.first < q.first;
                        });
                    for (size_t nz = 0; nz < NonZero; ++nz) {
                        xi_values_data[shift + nz] = data_inds[nz].first;
                        indicies_data[shift + nz] = data_inds[nz].second;
                    }
                    probabilities[i] = std::accumulate(
                        xi_values_data.begin() + shift,
                        xi_values_data.begin() + shift + NonZero, Vt(0)
                    );
                    IndexGenBigHelper::PrepareValues(
                        xi_values_data.begin() + shift,
                        xi_values_data.begin() + shift + NonZero, false
                    );
                }
            }
        }

        std::pair<size_t, Vt> make_step(
                size_t i0, Vt maxT, Vt xi1, Vt xi2,Vt xi3 = 0
            )const {
            if (i0 >= probabilities.size()) {
                return { i0,maxT };
            }
            Vt full_prob = (probabilities[i0] + EvapProbs[i0]);
            Vt dT = -std::log(xi2) / (probabilities[i0]+EvapProbs[i0]);
            if (dT < maxT) {
                if (xi3 * full_prob >= probabilities[i0]) {
                    return { probabilities.size(),maxT };
                }
                size_t m_index_begin = shift_data[i0];
                size_t m_index_end = shift_data[i0 + 1];
                if (m_index_end == m_index_begin) {
                    return { i0,maxT };
                }

                size_t k = IndexGenBigHelper::gen(
                    xi_values_data.begin() + m_index_begin,
                    xi_values_data.begin() + m_index_end, xi1
                );

                return { indicies_data[m_index_begin + k],dT };
            }
            else {
                return { i0,maxT };
            }
        }

        template <typename Gen_t>
        size_t evolute(size_t i0, Vt maxT, size_t max_steps, Gen_t& G)const {
            while (max_steps) {
                if (maxT > 0) {
                    auto [i1, dT] = make_step(i0, maxT, G(), G(),G());
                    maxT -= dT;
                    //dbg(i1,i0) += 1;
                    //dbg1(i1,i0) += dT;
                    i0 = i1;
                    //println("i': ",i0,", maxT: ", maxT,"\n----------->");
                }
                else {
                    return i0;
                }
                --max_steps;
            }
            return i0;
        }
    };

    template <typename Vt, typename Gen_t>
    Eigen::VectorX<Vt> MarkovEvolute(
        MarkovProcess<Vt> const& MP,
        VectorGenerator<Vt> const& Initials,
        Vt maxT,
        int Ntraj,
        size_t maxStep,
        Gen_t _G
    ) {
        Eigen::VectorX<int> Ret(MP.size());
        Ret.setZero();
        size_t seed = _G.state;
        int* ret_data = Ret.data();
        #pragma omp parallel for
        for (int i = 0; i < Ntraj; ++i) {
            auto G = _G;
            size_t m_seed = seed ^ (size_t)i + 1;
            G.set_seed(m_seed);

            size_t i0 = Initials.gen(G());
            int i1 = MP.evolute(i0, maxT, maxStep, G);
            if (i1 < MP.size()) {
                #pragma omp atomic
                ret_data[i1] += 1;
            }
            
        }
        return Ret.cast<Vt>() * (1 / (Vt)Ntraj);
    }
};