#pragma once
#include <grob/grid.hpp>
#include <grob/linear_interpolator.hpp>
#include "../body_potential.hpp"
#include "../form_factors.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>

#include "../traj_pool.hpp"
#include "../utils/quadtree.hpp"
#include "../utils/progress_bar.hpp"
#include "../dynamics/dynamics.hpp"
#include <list>
#include <algorithm>
#include "../dynamics/dynamics_evolutor.hpp"
#include <ranges>
#include "../trajectory.hpp"
#include "../dynamics/dynamic_evolutor_funcs.hpp"
namespace evdm {



	struct FFCImplStub {

		template <typename T, typename U>
		static FFCImplStub DeSerialize(T&& x, U&& y) {
			return {};
		}

		template <typename T>
		auto Serialize(T&& S) const {
			return grob::Serialize(false, S);
		}

		std::vector< ElementInfo_t> construct(size_t, size_t) const {
			return {};
		}
	};

	

	/// @brief type with quadtree (e,l) for each process ptype_in->ptype_out
	template <typename T, typename FormFactorConstructor_impl>
	struct ScatterInfoEL_in_t {
		typedef ScatterInfoElements<T> ScInfoEls_t;
		using QT_t = Quadtree<T, ScatterInfoElements<T>, ScatterSampleElements<T>>;
		using Node_t = typename QT_t::Node;

		typedef std::vector<QT_t>  QTrees_t;
		typedef std::vector<ElementsRVQuadTree<T>>  RV_Info_t;

		size_t ptypes;
		FormFactorConstructor_impl FFC_impl;
		std::vector<int> ValueIndexes;
		//-1 if no value
		std::vector<std::vector<ElementInfo_t>> FormFactors;
		//index of vectors - is ptype_in*ptypes + ptype_out
		QTrees_t QuadTrees;
		RV_Info_t RV_Info;
		ScatterInfoEL_in_t(
			size_t ptypes,
			FormFactorConstructor_impl _FFC_impl,
			std::vector<int> _ValueIndexes = {},
			QTrees_t _QuadTrees = {},
			RV_Info_t _RV_Info = {}
		) :
			ptypes(ptypes),
			FFC_impl(std::move(_FFC_impl)),
			ValueIndexes(std::move(_ValueIndexes)),
			QuadTrees(std::move(_QuadTrees)),
			RV_Info(std::move(_RV_Info))
		{
			if (ValueIndexes.size() < ptypes * ptypes) {
				ValueIndexes.resize(ptypes * ptypes, -1);
			}
			FormFactors.resize(QuadTrees.size());
			size_t m_size = QuadTrees.size();
			for (size_t i = 0; i < ptypes; ++i) {
				for (size_t j = 0; j < ptypes; ++j) {
					int index = ValueIndexes[i * ptypes + j];
					if (index >= (int)m_size) {
						throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t:"
							" invalid QuadTrees entry");
					}
					if (index >= 0) {
						FormFactors[index] = make_element_info(
							FFC_impl, i, j
						);
						//sort: first with bigger mN
						/*
						std::sort(
							FormFactors[index].begin(),
							FormFactors[index].end(),
							[](const ElementInfo_t & a,const ElementInfo_t & b){
								return a.mN > b.mN;
							}
						);*/
					}
				}
			}
		}
		struct gen_next_errors {
			constexpr static size_t NO_ERROR = 0;
			constexpr static size_t ELEMENT_FULL_PROB = 1;
			constexpr static size_t MAJORANT_ERROR = 2;
			constexpr static size_t ATTEMPT_THETA_ERROR = 4;
			constexpr static size_t ATTEMPT_ERROR = 8;
		};

		/// @brief generate next state and time delta after scatter. If no scatter - return current state and 0
		/// @return pair of next state and time delta
		template <typename Gen_t>
		std::pair<StateEL<T>, T> gen_next(
			StateEL<T> current, T max_time,
			T* p_in_p_out_ProbMatrix, const Body<T>& B, Gen_t& G,
			std::span<T> ElementBuffer,
			T MaxProbDiff_traj,
			std::pair < std::span<grob::Point<T,T>>, std::span<grob::Point<T, T>>> theta_buffers,
			bool FastDecay = false,
			size_t max_attempts = 1000,
			size_t errors = gen_next_errors::NO_ERROR,
			std::ostream * logs = nullptr
		) const {
			
			current.ptype = FastDecay ? 0 : current.ptype;
			//calculate scatter probability
			//generate ptype_out
			std::array<T, 4> ptype_out_probs;
			std::array<size_t, 4> ptype_out_indexes;
			if (ptypes > ptype_out_probs.size()) {
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::gen_next:"
					" too many ptypes for fixed size arrays");
			}
			size_t array_index = 0;
			for (size_t ptype_out = 0; ptype_out < ptypes; ++ptype_out) {
				if (hasValue(current.ptype, ptype_out)) {
					const QT_t& qt = getQT(current.ptype, ptype_out);
					typename QT_t::Node const& Node = qt.find_node(current.e, current.l);
					auto [alpha, betha] = Node.coeffs(current.e, current.l);
					auto [p00, p01, p10, p11] = Node.values_view(
						[](const ScInfoEls_t& siot) {
							return siot.full_prob;
						}
					);
					T p_elem = 0;
					if(Node.x1() != qt.X1)
					{
						p_elem = interpolate_log(p00, p01, p10, p11, alpha, betha);
					}
					else {
						p_elem = interpolate_lin(p00, p01, p10, p11, alpha, betha);
					}
					ptype_out_probs[array_index] =
						p_elem * p_in_p_out_ProbMatrix[current.ptype * ptypes + ptype_out];
					ptype_out_indexes[array_index] = ptype_out;
					++array_index;
				}
			}
			T scat_prob = std::accumulate(
				ptype_out_probs.begin(),
				ptype_out_probs.begin() + array_index,
				(T)0
			);

			if (scat_prob <= 0) {
				return { current,max_time };
			}
			//generate scatter time
			auto xi = E0I1_G(G);
			T dt_scat = -std::log(xi) / scat_prob;
			if (dt_scat >= max_time) {
				return { current,max_time };
			}
			//scatter
			//generate ptype_out
			for (size_t i = 1; i < array_index; ++i) {
				ptype_out_probs[i] += ptype_out_probs[i - 1];
			}

			T xi_pt = E0I1_G(G) * ptype_out_probs[array_index - 1];
			size_t sel_index = evdm::IndexGenBigHelper::gen(
				ptype_out_probs.begin(),
				ptype_out_probs.begin() + array_index,
				xi_pt
			);
			size_t ptype_out = ptype_out_indexes[sel_index];

			//generate e,l after scatter
			//generate element
			const QT_t& qt = getQT(current.ptype, ptype_out);
			QuadtreeNode<T, ScatterInfoElements<T>, ScatterSampleElements<T>> const& Node =
				qt.find_node(current.e, current.l);
			auto [alpha_el, beta_el] = Node.coeffs(current.e, current.l);
			//to generate element we need to:
			//1) make interpolated range from corners
			T xi_elem = G();

			auto ElementException = [current, errors]() {
				if (errors & gen_next_errors::ELEMENT_FULL_PROB) {
					std::ostringstream S;
					S << "MCEvolutor::ScatterInfoEL_in_t::gen_next:"
						" can't generate element index after scatter. cureent = ";
					S << "(" << current.ptype << ", " << current.e << ", " << current.l << ")";
					throw std::runtime_error(S.str());
				}
				};

			ScatterSampleElements<T> const& m_scatter_sampler = Node.uvalue();
			size_t sel_elem_index = m_scatter_sampler.gen_index(
				xi_elem, ElementBuffer, Node.values(), alpha_el, beta_el, ElementException);

			//3) get theta distribution for selected element
			// we will use min from corners
			auto xi_th = G();


			const QT_RV_t<T>& m_RV_Qt = RV_Info[
				getIndex(current.ptype, ptype_out)
			][sel_elem_index];
			/*
			const ScatterSampleTheta<T> & theta_sampler = m_scatter_sampler.Samples[sel_elem_index];
			auto ThetaGenException = [current,errors](){
				if (errors & gen_next_errors::ATTEMPT_THETA_ERROR) {
					std::ostringstream S;
					S << "MCEvolutor::ScatterInfoEL_in_t::gen_next:"
						" can't generate theta after scatter. cureent = ";
					S << "(" << current.ptype << ", " << current.e << ", " << current.l << ")";
					throw std::runtime_error(S.str());
				}
			};*/
			using namespace grob::literals;

			//find minimum of tuple with boost
			size_t il0_hint = m_scatter_sampler.i0_lmax();
			size_t il1_hint = m_scatter_sampler.i1_lmax();
			auto [L2max, r2_max] = B.maxL2(-current.e, il0_hint, il1_hint);

			T L2 = L2max * current.l * current.l;
			auto [rmin, rmax] = B.find_rmin_rmax(-current.e, L2, r2_max);

			T theta_max = B.get_theta_max(rmin, rmax);
			T Tout = eT(-current.e, L2);
			auto [T00, T01, T10, T11] = Node.values_view([sel_elem_index](
				const evdm::ScatterInfoElements<T>& cel) {
					return cel.Tin;
				});
			T Tin = interpolate_lin(T00, T01, T10, T11, alpha_el, beta_el);

			auto [m_points, m_buffer] = theta_buffers;
			
			T theta_undim = GenTrajProbs(
				B, current.e, current.l, rmin, rmax, theta_max, Tin, Tout, m_RV_Qt,
				m_points, m_buffer, 20, MaxProbDiff_traj, G()
			);
			/*
			T theta_undim = theta_sampler.gen_theta(
				dP_dth_functor<T>(
					B,getRVInfo(current.ptype,ptype_out)[sel_elem_index],
					current.e, rmin,rmax,theta_max,Tin+Tout),
				G,200,ElementBuffer
			);*/

			T theta_tmp = theta_max * theta_undim;
			auto cth2 = rmin * std::cos(theta_tmp / 2);
			auto sth2 = rmax * std::sin(theta_tmp / 2);
			T r = std::sqrt(upbound(cth2 * cth2 + sth2 * sth2, 1));
			T phi_r = B.Phi(r);
			T v2 = downbound(phi_r + current.e, 0);
			auto v = std::sqrt(v2);
			T vt = upbound(std::sqrt(L2) / r, v);
			//vec3<T> V_in = V_in_nd * VescMin;
			T vr = ssqrt(v2 - vt * vt);
			T vesc = B.Vesc;
			auto const& ff = FormFactors[getIndex(current.ptype, ptype_out)][sel_elem_index];
			T mu_chi = ff.mX * (1 / (ff.mN + ff.mX));
			T mu_i = ff.mN * (1 / (ff.mN + ff.mX));
			T Mu_chi_i = ff.mN * mu_chi;
			T Delta_M_chi = ff.dmX;
			T v2_delta = 2 * ff.dmX / Mu_chi_i / (vesc * vesc);
			T v1_temp_sq = downbound(B.Temp(r) / ff.mN / (vesc * vesc), 1e-10);
			T v1_tres = v2_delta > v2 ? std::sqrt(v2_delta) - v : T(0);
			T u_tresh = v2_delta > v2 ? (v1_tres * v1_tres) / (2 * v1_temp_sq) : T(0);



			typename QT_RV_t<T>::Node const& NodeRV =
				m_RV_Qt.find_node(r, v / std::sqrt(B.Phi(r)));


			for (size_t attempt = 0; attempt < max_attempts; ++attempt) {
				auto xi_rv = G();
				ScatterRVExp<T> const& m_scatter_rv = NodeRV.uvalue();

				size_t rv_index = m_scatter_rv.gen_index(G());
				// biggest majorant
				T majorant = m_scatter_rv.Majorants[rv_index];

				auto [xi0, xi1] = m_scatter_rv.box(rv_index);
				T xi = m_scatter_rv.min_value + ((xi0)+G() * (xi1 - xi0));
				T s = -std::log(xi);
				T u_in = (s + u_tresh);
				T v1 = std::sqrt(2 * u_in * v1_temp_sq);

				scatter_pre_result<T> m_s =
					gen_scatter(G(), G(), v, v1, u_in, -v2_delta);
				T Vu_in2 = v * v + v1 * v1 - 2 * v * v1 * m_s.cos_theta_in;
				T Vu_in = std::sqrt(Vu_in2);
				T Vu_out = std::sqrt(downbound(Vu_in2 - v2_delta, 0));

				T m_factor = get_ff_out(
					Vu_in, Vu_out, Mu_chi_i, v2_delta, vesc, m_s.cos_theta_out, ff.ff.sf);

				T from_factor = m_factor * m_s.gen_pre_form_factor;

				size_t p_type_out_eff = FastDecay ? 0 : ptype_out;
				if (from_factor > majorant) {
					if (errors & gen_next_errors::MAJORANT_ERROR) {
						std::ostringstream err;
						err << "evolute error: form factor "
							" is more than majorant ";
						err << from_factor / majorant;
						err << ", state = (" <<
							current.ptype << ", " << current.e << ", "
							<< current.l << ")";
						throw std::runtime_error(err.str());
					}
				}
				if (from_factor > majorant * G()) {
					auto [newE, newL] = get_el_out(
						mu_chi, mu_i, r, phi_r, vr, vt, v, v1, Vu_in, Vu_out,
						m_s.cos_theta_in, 2 * std::numbers::pi_v<T>*G(),
						m_s.cos_theta_out, 2 * std::numbers::pi_v<T>*G()
					);
					if (newE < 0) {

						T Lmax = std::sqrt(B.maxL2(-newE).first);

						T l = Lmax > 0 ? newL / Lmax : T(0);

						return { {p_type_out_eff,newE,l},  dt_scat };

					}
					else {
						return { {p_type_out_eff,newE,0}, max_time };
					}
				}
			}
			if (errors & gen_next_errors::ATTEMPT_ERROR) {
				std::ostringstream err;
				err << "evolute error: too many attempts to";
				err << ", state = (" <<
					current.ptype << ", " << current.e << ", "
					<< current.l << ")" << ", v = "  << v << ", vd^2 = " << v2_delta;
				throw std::runtime_error(err.str());
			}
			return { {ptype_out,1e10,1e10}, max_time };

		}

		T ScatterProb(T e, T l, size_t ptype_in, T* ptype_out_mlt_factors) const {
			T prob = 0;
			for (size_t ptype_out = 0; ptype_out < ptypes; ++ptype_out) {
				if (hasValue(ptype_in, ptype_out)) {
					const QT_t& qt = getQT(ptype_in, ptype_out);
					typename QT_t::Node const& Node = qt.find_node(e, l);
					auto [alpha, betha] = Node.coeffs(e, l);
					auto [p00, p01, p10, p11] = Node.values_view(
						[](const ScInfoEls_t& siot) {
							return siot.full_prob;
						}
					);
					T p_elem = interpolate_log(p00, p01, p10, p11, alpha, betha);
					prob += p_elem * ptype_out_mlt_factors[ptype_out];
				}
			}

			return prob;
		}

		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("ptypes", "FFC_impl", "ValueIndexes", "QuadTrees", "RV_Info"),
			PROPERTIES(ptypes, FFC_impl, ValueIndexes, QuadTrees, RV_Info)
		);
		DESERIALIZATOR_FUNCTION(
			ScatterInfoEL_in_t,
			PROPERTY_NAMES("ptypes", "FFC_impl", "ValueIndexes", "QuadTrees", "RV_Info"),
			PROPERTY_TYPES(ptypes, FFC_impl, ValueIndexes, QuadTrees, RV_Info)
		);


		void addQT(
			size_t ptype_in,
			size_t ptype_out,
			QT_t qt,
			ElementsRVQuadTree<T> rv_info
		) {
			int index = ValueIndexes[ptype_in * ptypes + ptype_out];
			if (index < 0) {
				ValueIndexes[ptype_in * ptypes + ptype_out] = QuadTrees.size();
				QuadTrees.push_back(std::move(qt));

				FormFactors.push_back(
					make_element_info(
						FFC_impl, ptype_in, ptype_out
					)
				);
				RV_Info.push_back(
					std::move(rv_info)
				);
			}
			else {
				QuadTrees[index] = std::move(qt);
				FormFactors[index] = make_element_info(
					FFC_impl, ptype_in, ptype_out
				);
				RV_Info[index] = std::move(rv_info);
			}
		}

		bool hasValue(size_t ptype_in, size_t ptype_out) const {
			return ValueIndexes[ptype_in * ptypes + ptype_out] >= 0;
		}
		int getIndex(size_t ptype_in, size_t ptype_out) const {
			return ValueIndexes[ptype_in * ptypes + ptype_out];
		}
		const ElementsRVQuadTree<T>& getRVInfo(size_t ptype_in, size_t ptype_out) const {
			if (!hasValue(ptype_in, ptype_out)) {
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::getRVInfo:"
					" no quadtree for given ptype_in,ptype_out");
			}
			return RV_Info[getIndex(ptype_in, ptype_out)];
		}
		ElementsRVQuadTree<T>& getRVInfo(size_t ptype_in, size_t ptype_out) {
			return const_cast<ElementsRVQuadTree<T> &>(
				static_cast<const ScatterInfoEL_in_t*>(this)->getRVInfo(ptype_in, ptype_out)
				);
		}
		const QT_t& getQT(size_t ptype_in, size_t ptype_out) const {
			if (!hasValue(ptype_in, ptype_out)) {
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::getRVInfo:"
					" no quadtree for given ptype_in,ptype_out");
			}
			return QuadTrees[getIndex(ptype_in, ptype_out)];
		}
		QT_t& getQT(size_t ptype_in, size_t ptype_out) {
			return const_cast<QT_t&>(
				static_cast<const ScatterInfoEL_in_t*>(this)->getQT(ptype_in, ptype_out)
				);
		}
		const std::vector<ElementInfo_t>& getFF(size_t ptype_in, size_t ptype_out) const {
			if (!hasValue(ptype_in, ptype_out)) {
				throw std::runtime_error("MCEvolutor::ScatterInfoEL_in_t::getRVInfo:"
					" no quadtree for given ptype_in,ptype_out");
			}
			return FormFactors[getIndex(ptype_in, ptype_out)];
		}
	};

	struct EvolutionExtraParams {
		std::list<std::string> m_errors = {};
		size_t max_theta_steps = 100;
		double prob_theta_accept = 0.05;
		bool fast_decay = false;
		std::ostream *logoutput = nullptr;
		size_t max_steps_eq = 0;
		size_t max_attempts = 1000;
	};
	template <typename T,typename FFC_impl_t>
	struct MCEvolutor {
		//content
		evdm::Body<T> B;
		//typedef decltype (B.get_le(0)) EL_Func_t;
		//EL_Func_t LE;
		ScatterInfoEL_in_t<T,FFC_impl_t> ScatterInfoEL;

		SERIALIZATOR_FUNCTION(
			PROPERTY_NAMES("body","ScatterInfoEL"), 
			PROPERTIES(B,ScatterInfoEL)
		);
		DESERIALIZATOR_FUNCTION(
			MCEvolutor, 
			PROPERTY_NAMES("body","ScatterInfoEL"), 
			PROPERTY_TYPES(B,ScatterInfoEL)
		);
		
		MCEvolutor(
			evdm::Body<T> B,
			//EL_Func_t LE,
			ScatterInfoEL_in_t<T,FFC_impl_t> ScatterInfoEL
		):
			B(std::move(B)),
			//LE(std::move(LE)),
			ScatterInfoEL(std::move(ScatterInfoEL))
		{}

		MCEvolutor(
			evdm::Body<T> B,
			size_t ptypes,
			FFC_impl_t FFC_impl
		):B(std::move(B)),ScatterInfoEL(ptypes,std::move(FFC_impl))
		{}
		size_t ptypes()const {
			return ScatterInfoEL.ptypes;
		}

		template <typename Gen_t>
		void AddProcess(
			size_t ptype_in,size_t ptype_out,
			size_t initial_rv_level, size_t max_rv_level,size_t max_rv_size, T rv_cmp_accept,
			size_t initial_el_level, size_t max_el_level,size_t max_el_size,T el_cmp_accept,
			size_t ThetaMaxSteps, T ThetaAccept, T zeroProb, Gen_t G, size_t Nmk,
			std::ostream*logs = nullptr){

			std::vector<ElementInfo_t> m_elements = make_element_info(
				ScatterInfoEL.FFC_impl,ptype_in,ptype_out
			);
			ElementsRVQuadTree<T> m_qt;
			m_qt.reserve(m_elements.size());
			

			for(size_t i=0;i<m_elements.size();++i){
				const auto & m_el = m_elements[i];
				auto rv_qt_filler = [&](T r,T v_undim)->ScatterRVExpInfo_t<T> {
					auto vmax = std::sqrt(B.Phi(r));
					return genRVinfo(ptype_in,ptype_out,r,vmax*v_undim,m_el.ff.n_e(r),m_el,G,Nmk);
				};
				m_qt.push_back( 
					make_rv_tree (
						rv_qt_filler,initial_rv_level,max_rv_level,max_rv_size,
						zeroProb/m_elements.size(),rv_cmp_accept
					)
				);
				debug_args(logs, "", "make rv tree for element: (",
					m_elements[i].A, ", ", m_elements[i].Z, ")");

			}

			// 2) using RV info fill elinfo quadtree
			//	critrium to split: 
			//		bad interpolation of probs for any element
			using EL_QT_t = Quadtree<T,ScatterInfoElements<T>,ScatterSampleElements<T>>;
			auto el_filler = [&](T e,T l)->ScatterInfoElements<T> {
				auto const & m_r_grid = B.Rho.Grid;
				auto [L2,rme] = B.maxL2(-e);
				auto get_poses = [&m_r_grid](T r)->size_t {
					size_t i = m_r_grid.pos(r);
					return i;
				};
				auto il_min= get_poses(rme);
				auto [rmin,rmax] = B.find_rmin_rmax(-e,L2*l*l,rme);
				auto i0_rmin = get_poses(rmin);
				auto i0_rmax = get_poses(rmax);
				std::vector<T> m_traj_prob_for_elements(m_elements.size());
				std::vector<grob::Point<T, T>> buffers(ThetaMaxSteps*2);
				std::span<grob::Point<T, T>> buffer1(buffers.begin(), ThetaMaxSteps);
				std::span<grob::Point<T, T>> buffer2(buffers.begin()+ ThetaMaxSteps, ThetaMaxSteps);
				T theta_max = B.get_theta_max(rmin, rmax);
				T Tin = B.get_internal_period(rmin, rmax, 200);
				T Tout = eT(-e, L2 * l * l);

				for (size_t i=0;i<m_qt.size();++i){
					m_traj_prob_for_elements[i] = GenTrajProbs(
						B, e, l, rmin, rmax, theta_max, Tin, Tout, 
						m_qt[i], buffer1, buffer2, 20, ThetaAccept
					);
				}

				return ScatterInfoElements<T>(
					Tin,
					std::move(m_traj_prob_for_elements),
					std::make_tuple(il_min,i0_rmin,i0_rmax)
				);
			};
			auto m_el_tree = make_el_tree(
				el_filler, -B.Phi[0],
				initial_el_level, max_el_level, max_el_size, zeroProb,
				el_cmp_accept
			);
			debug_args(logs, "", "make el tree");
			ScatterInfoEL.addQT(
				ptype_in,ptype_out,
				std::move(m_el_tree), std::move(m_qt)
			);
		}
		/*
		void AddSameProcess(size_t ptype_in_dsc,size_t ptype_out_dsc,
		size_t ptype_in_src,size_t ptype_out_src){
			ScatterInfoEL.addQT(ptype_in_dsc,ptype_out_dsc,
					ScatterInfoEL.getQT(ptype_in_src,ptype_out_src),
					ScatterInfoEL.getRVInfo(ptype_in_src,ptype_out_src)
			);
		}*/

		template <typename Gen_t>
		void evolute(
			std::span<StateEL<T>> init,
			Gen_t&& _G,
			std::array<T, 16> prob_matrix,
			const T evolve_time,
			const size_t max_scatter,
			EvolutionExtraParams m_params = {}
		) const{

			debug_args(m_params.logoutput, "", "evolution: ", " is FD: ", 
				m_params.fast_decay, ", time: ", evolve_time, ", max_steps_eq: ", m_params.max_steps_eq);
			if (m_params.max_steps_eq == 0) {
				m_params.max_steps_eq = max_scatter;
			}
			typedef typename ScatterInfoEL_in_t<T, FFC_impl_t>::gen_next_errors err_t;
			size_t m_err = 0;
			auto m_errors = m_params.m_errors;
			for (std::string_view err : m_errors) {
				if (err == "element") {
					m_err |= err_t::ELEMENT_FULL_PROB;
				}
				if (err == "majorant") {
					m_err |= err_t::MAJORANT_ERROR;
				}
				if (err == "theta") {
					m_err |= err_t::ATTEMPT_THETA_ERROR;
				}
				if (err == "attempt") {
					m_err |= err_t::ATTEMPT_ERROR;
				}
			}
			int N = init.size();
			const auto _seed = _G.state;
			std::array<T,80> ElementBuffer;
			std::vector<grob::Point<T,T>> prob_buffers(m_params.max_theta_steps * 2);
			
			std::atomic<bool> error_occurred(false);
			std::exception_ptr global_exception_ptr = nullptr;

			std::vector< StateTimeEL<T>> TrajectoryStory;
			#pragma omp parallel for firstprivate(ElementBuffer) firstprivate(prob_buffers) firstprivate(TrajectoryStory)
			for(int i=0;i<N;++i){
				if(error_occurred.load()){
					continue;
				}
				auto G = _G;
				G.set_seed(
					hash_combine(
						hash_float(init[i].e), hash_float(init[i].l),
						hash_seed(init[i].ptype), hash_seed(i), hash_seed(_G.state)
					)
				);

				std::span<grob::Point<T,T>> m_buf1{ prob_buffers.data(),m_params.max_theta_steps };
				std::span<grob::Point<T,T>> m_buf2{ prob_buffers.data() + m_params.max_theta_steps,m_params.max_theta_steps };

				T max_scatter1 = max_scatter+std::floor( G()*ScatterInfoEL.ptypes);
				StateEL<T> tmp_state = init[i];
				T time_remain = evolve_time;
				T time_delta = 0;
				size_t j = 0;
				for(;j< max_scatter1 && !error_occurred.load() ;++j){
					
					try {
						std::tie(tmp_state, time_delta) = ScatterInfoEL.gen_next(
							tmp_state, time_remain, prob_matrix.data(), B, G, ElementBuffer,
							m_params.prob_theta_accept, { m_buf1,m_buf2 },m_params.max_attempts,
							m_params.fast_decay, m_err
						);
						time_remain -= time_delta;
						if (time_remain <= 0) {
							break;
						}
					}
					catch (std::runtime_error &e) {
						std::ostringstream S;
						S << e.what() << "; at i = " << i << ", j = " << j;

						if (!error_occurred.exchange(true)) {
							global_exception_ptr = std::make_exception_ptr(std::runtime_error(S.str()));
						}
						break;
					}
				}
				if (j == max_scatter1 && time_remain > 0) {
					if (TrajectoryStory.size() < m_params.max_steps_eq)
						TrajectoryStory.resize(m_params.max_steps_eq);
					
					TrajectoryStory[0] = { tmp_state, 0};
					T last_time = 0;
					for(size_t k=0;k<TrajectoryStory.size()  && !error_occurred.load(); ++k){
						try {
							T delta_time;
							std::tie(tmp_state, delta_time) = ScatterInfoEL.gen_next(
								TrajectoryStory[k], time_remain, prob_matrix.data(), B, G, ElementBuffer,
								m_params.prob_theta_accept, { m_buf1,m_buf2 },
								m_params.fast_decay, m_err
							);
							last_time += delta_time;
							TrajectoryStory[k].t = last_time;
							if(k+1 < TrajectoryStory.size()){
								TrajectoryStory[k+1] = { tmp_state, 0 };
							}
							
							time_remain -= delta_time;
							if (time_remain <= 0) {
								break;
							}
							
						}
						catch (std::runtime_error& e) {
							std::ostringstream S;
							S << e.what() << "; at i = " << i << ", j = " << j;

							if (!error_occurred.exchange(true)) {
								global_exception_ptr = std::make_exception_ptr(std::runtime_error(S.str()));
							}
							break;
						}
					}
					if(time_remain > 0 && !error_occurred.load()){
						// we should generate state according to trajectory story
						auto time_view = TrajectoryStory | std::views::transform([](const StateTimeEL<T>& st) {
							return st.t;
						});
						size_t index = IndexGenBigHelper::gen(time_view.begin(), time_view.end(),  time_view.back() * G());
						tmp_state = TrajectoryStory[index];
					}
				} 
				init[i] = tmp_state;
				
			}
			if (global_exception_ptr) {
				std::rethrow_exception(global_exception_ptr);
			}
		}

		private:

		/// @brief get for particular r and v majorants and prob scatter
		template <typename Gen_t>
		ScatterRVExpInfo_t<T> genRVinfo(
			size_t ptype_in, size_t ptype_out,
			T r,T v,T n_r, ElementInfo_t const & El,Gen_t &G,size_t Nmk)
		{
			T mu_chi = El.mX*(1/(El.mN+El.mX));
			T mu_i = El.mN*(1/(El.mN+El.mX));
			T Mu_chi_i = El.mN*mu_chi;
			T Delta_M_chi = El.dmX;
			T vesc = B.Vesc;
			T v2_delta = 2*El.dmX/Mu_chi_i/(vesc*vesc);
			

			T v2 = v*v;
			T v1_temp_sq = (B.Temp(r) + 1e-12)/El.mN/(vesc*vesc);
			T v1_tres = v2_delta > v2 ? std::sqrt(v2_delta) - v : T(0);
			T u_tresh = v1_tres ? (v1_tres*v1_tres)/(2*v1_temp_sq) : T(0);


			using RVPi = ScatterRVExpInfo_t<T>;
			T sum = 0;
			std::array<T, RVPi::size> max_ffs;
			for (int i=0;i<RVPi::size;++i){
				auto [xiu_0,xiu_1] = RVPi::box(i);
				T max_ff = 0;
				if (u_tresh < 40) {
					for (int _k = 0; _k < Nmk; ++_k) {
						T xi = RVPi::min_value + (xiu_0 + (xiu_1 - xiu_0) * G());
						T s = -std::log(xi);
						T u_in = (s + u_tresh);
						T v1 = std::sqrt(2 * u_in * v1_temp_sq);
						scatter_pre_result<T> m_s =
							gen_scatter(G(), G(), v, v1, u_in, -v2_delta);

						T Vu_in2 = v * v + v1 * v1 - 2 * v * v1 * m_s.cos_theta_in;
						T Vu_in = std::sqrt(Vu_in2);
						T Vu_out = std::sqrt(downbound(Vu_in2 - v2_delta, 0));

						T m_factor = get_ff_out(
							Vu_in, Vu_out, Mu_chi_i, v2_delta, vesc, m_s.cos_theta_out, El.ff.sf);
						T from_factor = m_factor * m_s.gen_pre_form_factor;
						max_ff = std::max(max_ff, from_factor);
						sum += (xiu_1 - xiu_0)*from_factor * n_r * std::exp(-u_tresh);
					}
				}
				max_ffs[i] = max_ff*T(1.05);
			}
			return ScatterRVExpInfo_t<T>(std::move(max_ffs),sum/Nmk,u_tresh);
		}

		public:

		template <typename ValueGetter_t>
		QT_RV_t<T> make_rv_tree(
			ValueGetter_t && value_getter,
			size_t initial_level, size_t max_level, size_t max_bins,
			T zeroProb, T cmp_accept
		) const{
			// TODO 1) for each element fill RV info quadtree
			// critrium to split: 
			//		bad interpolation of probs
			//		big difference in treshold
			//		big difference in majorants
			auto GetProb = [](ScatterRVExpInfo_t<T> const & x){
				return x.full_prob;
			};
			auto comparator = [zeroProb](
				ScatterRVExpInfo_t<T> const & R1,
				ScatterRVExpInfo_t<T> const & R2) ->T
			{
				auto get_delta = [](T x,T y){
					if(x>y){
						std::swap(x,y);
					}
					return ( x > 0 ? y/x : T(1e10) ) - 1;
				};
				T _x = R1.full_prob;
				T _y = R2.full_prob;
				T delta_fp = (_x > zeroProb || _y > zeroProb) ? get_delta(_x, _y) : 0;
				return delta_fp;
			};
			QT_RV_t<T> m_tree =  make_quadtree<ScatterRVExpInfo_t<T>,ScatterRVExp<T>>(
				value_getter,0,1,0,1,initial_level,max_level,
				max_bins,zeroProb,GetProb,comparator,cmp_accept
			);
			typedef  QuadtreeNode<T,ScatterRVExpInfo_t<T>,ScatterRVExp<T>> Node_t;
			for(Node_t & nd : m_tree){
				auto const & [I1,I2,I3,I4] = nd.values();
				auto majorants_view = grob::as_container(
					[&I1,&I2,&I3,&I4](size_t i){
						return std::max(
							std::max(I1.Majorants[i],I2.Majorants[i]),
							std::max(I2.Majorants[i],I3.Majorants[i])
						);
					},I1.Majorants.size()
				);
				nd.uvalue() = ScatterRVExp<T>(majorants_view.begin(),majorants_view.end());
			}
			return m_tree;
		}

		template <typename ValueGetter>
		static Quadtree<T,ScatterInfoElements<T>,ScatterSampleElements<T>> make_el_tree(
			ValueGetter &&Fel,T Emin, 
			size_t initial_level, size_t max_level, size_t max_bins,
			T zeroProb, T cmp_accept
		){
			auto FuncGetProb = [](ScatterInfoElements<T> const & x){
				return x.full_prob;
			};
			auto cmp = [zeroProb](ScatterInfoElements<T> const & x,
				ScatterInfoElements<T> const & y)->T 
			{
				auto get_delta = [](T x,T y){	
					if(x>y){
						std::swap(x,y);
					}
					return ( x > 0 ? y/x : T(1e10) ) - 1;
				};
				T _x = x.full_prob;
				T _y = y.full_prob;
				T delta_fp = (_x > zeroProb || _y > zeroProb) ? get_delta(_x,_y) : 0;
				T max_delta = delta_fp;
				/*
				for(size_t i=0;i<std::min(x.ProbsElements.size(),y.ProbsElements.size());++i){
					T _x = x.ProbsElements[i];
					T _y = y.ProbsElements[i];
					if(_x > zeroProb ||  _y> zeroProb){
						max_delta = std::max(max_delta,get_delta(_x,_y));
					} 
				}*/

				return max_delta;
			};
			Quadtree<T,ScatterInfoElements<T>,ScatterSampleElements<T>> m_tree =  
			make_quadtree<ScatterInfoElements<T>, ScatterSampleElements<T>>(
				Fel, Emin,0,0,1,initial_level,max_level,max_bins,zeroProb,
				FuncGetProb,cmp,cmp_accept
			);

			for (QuadtreeNode<T,ScatterInfoElements<T>,ScatterSampleElements<T>> & N : m_tree){
				N.uvalue() = ScatterSampleElements<T>(N(0,0),N(0,1),N(1,0),N(1,1));
			}

			return m_tree;
		}
		
		template <typename Val_t,typename UV_t,
			typename ValueGetter_t, //e,l -> Val_t
			typename ProbGetter, //Val_t & -> probability
			typename RelComparator > // (Val_t, Val_t) -> some relative difference
		static Quadtree<T,Val_t,UV_t> make_quadtree(
			ValueGetter_t && F,T x0,T x1,T y0,T y1,
			size_t initial_level, size_t max_lvl,size_t max_size,
			T zeroProb, ProbGetter && FuncGetProb, RelComparator && cmp, 
			T cmp_accept
		){
			
			using Node_t = QuadtreeNode<T,Val_t,UV_t>;
			Quadtree<T,Val_t,UV_t> m_tree(F,x0,x1,y0,y1);
			auto ValueGetterVectorized =[&F](auto const & coords){
				int N = coords.size();
				std::vector<std::shared_ptr<Val_t>> ret(N);
				#pragma omp parallel for
				for(int i=0;i<N;++i){
					auto [x,y] = coords[i];
					ret[i] = std::make_shared<Val_t>(F(x,y));
				}
				return ret;
			};
			auto always_true_cond = [zeroProb,&FuncGetProb](auto const& nd){
				return true;
			};
			auto CMPaccept_condition = [&](T accept_error){ 
				return [accept_error,zeroProb,max_lvl,&cmp,&FuncGetProb](Node_t const& nd){
					if(nd.level() > max_lvl){
						return false;
					}
					//auto const & [P0,P1,P2,P3]
					auto const & [P0,P1,P2,P3] = nd.values();
					auto Ptuple = nd.values_view(FuncGetProb);
					auto Pmax = tuple_max(Ptuple);
					if(Pmax <= zeroProb)
						return false;
					auto cond = [&cmp,accept_error](auto const & x,auto const & y){
						return cmp(x,y) > accept_error;
					};
					return cond(P0,P1) || cond(P0,P2) || cond(P0,P3) ||
						   cond(P1,P2) || cond(P1,P3) || cond(P2,P3);
				};
			};

			
			// refine to initial level
			size_t remain_size = max_size;
			remain_size -= 1;
			for(size_t _i=0; _i<initial_level;++_i){
				remain_size -= m_tree.refine_vectorized(
					always_true_cond,ValueGetterVectorized,1,remain_size 
				);
			}

			// some funtions
			auto P_cmp_coarse_condition = CMPaccept_condition(1);
			auto P_cmp_strict_condition = CMPaccept_condition(cmp_accept);
			
			auto apply_condition = [&](auto && m_cond,size_t m_rem_size){
				size_t new_nodes = 0;
				size_t delta=1;
				while(delta){
					delta = m_tree.refine_vectorized(
						m_cond,ValueGetterVectorized,1,m_rem_size
					);
					new_nodes += delta;
					m_rem_size-=delta;
				}
				return new_nodes;
			};

			/// Step 1: we should refine most strict condition
			remain_size -= apply_condition(P_cmp_strict_condition,remain_size/2);
			/// Step 2:  we should refine less strict condition for prob
			remain_size -= apply_condition(P_cmp_coarse_condition,remain_size);
			/// Step 3: if we still have remain_size, we coult try to refine strict
			remain_size -= apply_condition(P_cmp_strict_condition,remain_size);
			
			return m_tree;
		}
		PlotInfo<T> plotinfo_rv(
		size_t ptypein,size_t ptypeout,
		size_t elementA,size_t elementZ,
		std::string_view what_to_plot = "")
		{
			std::vector<ElementInfo_t> const & qtffs = 
				ScatterInfoEL.getFF(ptypein,ptypeout);
			auto it = std::find_if(qtffs.begin(),qtffs.end(),
				[elementA,elementZ](ElementInfo_t const & info){
					return info.A == elementA && info.Z == elementZ;
				}
			);
			if(it != qtffs.end()){
				size_t index = it-qtffs.begin();
				QT_RV_t<T> const & m_qt = ScatterInfoEL.getRVInfo(ptypein,ptypeout)[index];
				
				if(what_to_plot == "prob" || what_to_plot.empty()){
					return quadtree_plotinfo(m_qt,
						[](ScatterRVExpInfo_t<T> const & V){
							return V.full_prob; 
						}
					);
				} else {
					throw std::runtime_error("can plot only full prob of rv tree");
				}
			} else {
				throw std::runtime_error("element not found");
			}
			
		}

		PlotInfo<T> plotinfo_el(
			size_t ptypein,size_t ptypeout,
			std::string_view what_to_plot = ""
		){
			auto const & qtel = 
				ScatterInfoEL.getQT(ptypein,ptypeout);
			constexpr size_t PLOT_RMIN = 0;
			constexpr size_t PLOT_RMAX = 1;
			constexpr size_t PLOT_L2MAX = 2;
			constexpr size_t PLOT_Tin = 3;
			constexpr size_t PLOT_Tout = 4;
			constexpr size_t PLOT_T = 5;
			constexpr size_t PLOT_ThetaMax = 6;
			constexpr size_t PLOT_Prob = 7;

			size_t wtp = 7;
			std::array<const char *,7> strs = {
				"rmin","rmax","L2","Tin","Tout","T","theta_max"
			};
			wtp = std::find(strs.begin(), strs.end(), what_to_plot) - strs.begin();

			return quadtree_plotinfo_xyz(qtel, 
				[&,wtp](auto e, auto l, auto const& el) -> T
				{
					if (wtp == PLOT_Prob) {
						return el.full_prob;
					}
					auto [L2m,rm] = B.maxL2(-e);
					if (wtp == PLOT_L2MAX) {
						return L2m;
					}
					auto L2 = L2m * l * l;
					auto [rmin,rmax] = B.find_rmin_rmax(-e, L2, rm,0);
					if (wtp == PLOT_RMIN) {
						return rmin;
					}
					if (wtp == PLOT_RMAX) {
						return rmax;
					}
					if (wtp == PLOT_ThetaMax) {
						return B.get_theta_max(rmin, rmax);
					}
					auto Tin = B.get_internal_period(rmin, rmax, 100);
					if (wtp == PLOT_Tin) {
						return Tin;
					}

					auto Tout = eT(-e, L2);
					if (wtp == PLOT_Tout) {
						return Tout;
					}

					if (wtp == PLOT_T) {
						return Tin+Tout;
					}

					return 0;
				});
		}

	};
	

};
