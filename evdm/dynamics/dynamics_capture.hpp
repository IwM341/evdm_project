#pragma once 

#include "../form_factors.hpp"
#include "../core.hpp"
#include <grob/grid.hpp>
#include <grob/grid_objects.hpp>
#include <numbers>
#include "../utils/progress_bar.hpp"
#include "../measure.hpp"
namespace evdm {


	using FormFactor_t = QexpFactors<std::index_sequence<2, 4, 8, 12>>;
	
	template <typename T>
	using vec3 = Eigen::Vector3<T>;
	

	template <typename T>
	using vec2 = Eigen::Vector2<T>;

	struct ScatterEvent {
		grob::GridFunction <
			grob::linear_interpolator,
			grob::GridUniform<float>,
			std::vector<float>
		> n_e;
		FormFactor_t sf;
		
		template <typename Array_t>
		ScatterEvent(size_t arr_size, Array_t && N_values, FormFactor_t _ffact):
			n_e(
				grob::GridUniform<float>(0,1, arr_size),std::vector<float>(arr_size)
			), sf(_ffact) {
			for (size_t i = 0; i < arr_size; ++i) {
				n_e.Values[i] = N_values[i];
			}
		}
		inline std::string repr()const {
			std::ostringstream S;
			S << "ScatterEvent( size = " << n_e.size() << ", sf = " << sf.repr() << ")";
			return S.str();
		}

		template <typename scale_t>
		void rescale(scale_t x) {
			sf.rescale(x);
		}
	};

	/// <summary>
	/// return gauss vector 3 with weight
	/// </summary>
	/// <param name="G">random gen</param>
	/// <param name="Vdisp">dispersion gauss</param>
	/// <param name="p">probability of using standart generator, 
	/// 1-p - using uniform generator</param>
	/// <param name="max_xi"></param>
	/// <returns></returns>
	template <class Gen_t>
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> Gauss3_BeyondD(
		Gen_t&& G, Gen_vt<Gen_t> Vdisp, 
		Gen_vt<Gen_t> p = (Gen_vt < Gen_t>)0.8,
		Gen_vt<Gen_t> max_xi = (Gen_vt < Gen_t>)8)
	{
		using T = Gen_vt<Gen_t>;
		using mvec3 = vec3<T>;
		if (G() < p) {
			auto V1 = Vdisp * sqrt(-2 * log(1 - G()));
			auto V2 = Vdisp * sqrt(-2 * log(1 - G()));
			auto phi1 = RandomPhi(G);
			auto phi2 = RandomPhi(G);
			return { 
				mvec3(V1 * cos(phi1), V1 * sin(phi1), V2 * cos(phi2)),
				(T)1 
			};
		}
		else {
			T r = std::cbrt(G()) * max_xi;
			T phi = RandomPhi(G);
			T cos_Th = RandomCos(G);
			T sin_Th = std::sqrt(1 - cos_Th* cos_Th);
			T r_xy = r * sin_Th;
			T fac = exp(-r*r / 2) / std::pow(2 * (T) std::numbers::pi, (T)1.5);
			return { 
				Vdisp * mvec3(r_xy * std::sin(phi),r_xy*std::sin(phi),r*cos_Th),
				fac
			};
		}
	}

	template <typename Gen_t>
	vec3<Gen_vt<Gen_t>> RandomNvec(Gen_t&& G) {
		auto cosT = RandomCos(G);
		auto sinT = std::sqrt(1 - cosT * cosT);
		auto phi = RandomPhi(G);

		return {sinT * cos(phi),sinT * sin(phi),cosT};
	}

	template <class Gen_t>
	/*MK generator of input velocity*/
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> HaloVelocity(
		Gen_t&& G,Gen_vt<Gen_t> VescTmp,
		Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> mU0) 
	{
		using T = Gen_vt<Gen_t>;
		auto ksi = sqrt(-2 * std::log( E0I1_G(G)) );
		auto phi = RandomPhi(G)/2;


		auto sinPhi = sin(phi);
		auto cosPhi = cos(phi);

		auto u0 = mU0 / Vdisp;
		auto ve = VescTmp / Vdisp;

		auto u = sqrt(u0 * u0 + ksi * ksi + 2 * u0 * ksi * cosPhi);


		auto sinTheta = (u != 0.0 ? ksi * sinPhi / u : 0);
		auto v = sqrt(u * u + ve * ve);
		
		auto n = RandomNvec(G);
		
		return MCResult<vec3<T>,T>(n * (v * Vdisp), 
			sinTheta * v * sqrt((T)std::numbers::pi / 2));
	}
	template <class Gen_t>
	/*MK generator of input velocity*/
	inline MCResult<vec3< Gen_vt<Gen_t>>, Gen_vt<Gen_t>> 
		HaloVelocityConstrained(
		Gen_t&& G, Gen_vt<Gen_t> VescTmp, Gen_vt<Gen_t> Vmin,
		Gen_vt<Gen_t> Vmax,Gen_vt<Gen_t>  Vdisp, Gen_vt<Gen_t>  mU0)
	{
		typedef Gen_vt<Gen_t> number_t;
		number_t umin = (Vmin > VescTmp ? sqrt(Vmin * Vmin - VescTmp * VescTmp) / Vdisp : 0.0);
		number_t umax = ((Vmax > Vmin && Vmax > VescTmp) ? sqrt(Vmax * Vmax - VescTmp * VescTmp) / Vdisp : umin);

		number_t u0 = mU0 / Vdisp;

		number_t ksi_max = umax + u0;
		number_t ksi_min = (mU0 - umax > 0 ? u0 - umax : (umin - u0 > 0 ? umin - u0 : 0));

		number_t ksi_rd = exp(-ksi_min * ksi_min / 2) - exp(-ksi_max * ksi_max / 2);

		number_t ksi = sqrt(-2 * log(exp(-ksi_max * ksi_max / 2) + G() * ksi_rd));

		number_t cosThMin = std::min(1.0, std::max(-1.0, (umin * umin - ksi * ksi - u0 * u0) / (2 * ksi * u0)));
		number_t cosThMax = std::max(-1.0, std::min(1.0, (umax * umax - ksi * ksi - u0 * u0) / (2 * ksi * u0)));

		number_t max_theta = acos(cosThMin);
		number_t min_theta = acos(cosThMax);

		number_t rd_th = (max_theta - min_theta) / (number_t)std::numbers::pi;

		number_t theta = min_theta + (max_theta - min_theta) * G();

		number_t ve = VescTmp / Vdisp;

		number_t u = sqrt(u0 * u0 + ksi * ksi + 2 * u0 * ksi * cos(theta));

		number_t sinTheta = (u != 0.0 ? ksi * sin(theta) / u : 0);

		auto v = sqrt(u * u + ve * ve);
		auto n = RandomNvec(G);

		return { n * (v * Vdisp), 
			ksi_rd * rd_th * sinTheta * v * sqrt((number_t)std::numbers::pi/ 2) 
		};
	}

	template <class Gen_t,typename VectorT>
	/*MK generator of output nu'*/
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> VuOut(
		Gen_t && G, const  VectorT& Vcm, const VectorT& Vu,
		Gen_vt<Gen_t> Vesc, Gen_vt<Gen_t> VescMin, 
		Gen_vt<Gen_t>  mp, Gen_vt<Gen_t>  mk, 
		Gen_vt<Gen_t>  mp_frac, Gen_vt<Gen_t>  mk_frac,
		Gen_vt<Gen_t>  m_cm, Gen_vt<Gen_t> deltaE = 0) {

		typedef Gen_vt<Gen_t> number_t;
		number_t VcmN = Vcm.norm();
		typedef vec3<number_t> vec3_t;
		vec3_t  n_v = Vcm / VcmN;

		number_t cosThetaVcm = n_v[2];
		number_t sinThetaVcm = sqrt(n_v[0] * n_v[0] + n_v[1] * n_v[1]);

		number_t cosPhiVcm = 1.0;
		number_t sinPhiVcm = 0.0;

		if (sinThetaVcm > 1e-10) {
			cosPhiVcm = n_v[0] / sinThetaVcm;
			sinPhiVcm = n_v[1] / sinThetaVcm;
		}

		vec3_t n_1(cosThetaVcm * cosPhiVcm, cosThetaVcm * sinPhiVcm, -sinThetaVcm);
		vec3_t n_2(-sinPhiVcm, cosPhiVcm, 0);

		/*
		std::cout << n_1*n_1 << "\t" << n_1*n_2 << "\t" << n_1*n_v << std::endl;
		std::cout << n_2*n_1 << "\t" << n_2*n_2 << "\t" << n_2*n_v << std::endl;
		std::cout << n_v*n_1 << "\t" << n_v*n_2 << "\t" << n_v*n_v << std::endl<< std::endl;
		*/

		number_t Vu1_squared = Vu.squaredNorm() - deltaE * 2/ m_cm;
		if (Vu1_squared <= 0.0)
			return MCResult<vec3_t, number_t>(vec3_t(0, 0, 0), 0);

		number_t Vu1 = sqrt(Vu1_squared);

		number_t cosTh1max = (
			Vesc * Vesc - Vu1_squared * mp_frac * mp_frac - VcmN * VcmN
			) / (2 * VcmN * Vu1* mp_frac);

		if (!(cosTh1max > -1))
			return MCResult<vec3_t, number_t>(vec3_t(0, 0, 0), 0);
		else if (cosTh1max >= 1) {
			cosTh1max = 1;
		}

		number_t cosTh1 = (1 + cosTh1max) * G() - 1;
		number_t sinTh1 = sqrt(1.0 - cosTh1 * cosTh1);
		number_t phi1 = RandomPhi(G);

		const vec3_t Vu1_vec = Vu1 * (n_v * cosTh1 + n_1 * sinTh1 * cos(phi1) + n_2 * sinTh1 * sin(phi1));
		return MCResult<vec3_t, number_t>(Vu1_vec, (1 + cosTh1max)/2 * Vu1 / VescMin);
	}

	/*MK generator of output vu', considering evaporation */
	template <class Gen_t, typename VectorT>
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> VuOut_evap(
		Gen_t&& G, const  VectorT& Vcm, const VectorT& Vu,
		Gen_vt<Gen_t> Vesc, Gen_vt<Gen_t> VescMin,
		Gen_vt<Gen_t>  mp, Gen_vt<Gen_t>  mk,
		Gen_vt<Gen_t>  mp_frac, Gen_vt<Gen_t>  mk_frac,
		Gen_vt<Gen_t>  m_cm, Gen_vt<Gen_t> deltaE = 0) {

		typedef Gen_vt<Gen_t> number_t;
		number_t VcmN = Vcm.norm();
		typedef vec3<number_t> vec3_t;
		vec3_t  n_v = Vcm / VcmN;

		number_t cosThetaVcm = n_v[2];
		number_t sinThetaVcm = sqrt(n_v[0] * n_v[0] + n_v[1] * n_v[1]);

		number_t cosPhiVcm = 1.0;
		number_t sinPhiVcm = 0.0;

		if (sinThetaVcm > 1e-10) {
			cosPhiVcm = n_v[0] / sinThetaVcm;
			sinPhiVcm = n_v[1] / sinThetaVcm;
		}

		vec3_t n_1(cosThetaVcm * cosPhiVcm, cosThetaVcm * sinPhiVcm, -sinThetaVcm);
		vec3_t n_2(-sinPhiVcm, cosPhiVcm, 0);

		/*
		std::cout << n_1*n_1 << "\t" << n_1*n_2 << "\t" << n_1*n_v << std::endl;
		std::cout << n_2*n_1 << "\t" << n_2*n_2 << "\t" << n_2*n_v << std::endl;
		std::cout << n_v*n_1 << "\t" << n_v*n_2 << "\t" << n_v*n_v << std::endl<< std::endl;
		*/

		number_t Vu1_squared = Vu.squaredNorm() - deltaE * 2 / m_cm;
		if (Vu1_squared <= 0.0)
			return MCResult<vec3_t, number_t>(vec3_t(0, 0, 0), 0);

		number_t Vu1 = sqrt(Vu1_squared);

		number_t cosTh1 = 2* G() - 1;
		number_t sinTh1 = sqrt(1.0 - cosTh1 * cosTh1);
		number_t phi1 = RandomPhi(G);

		const vec3_t Vu1_vec = Vu1 * 
			(n_v * cosTh1 + n_1 * sinTh1 * cos(phi1) + n_2 * sinTh1 * sin(phi1));
		return MCResult<vec3_t, number_t>(Vu1_vec, Vu1 / VescMin);
	}

	/// <summary>
	/// orthogonal complement
	/// </summary>
	/// <param name="V">initial vector</param>
	/// <returns>pair of vectors 3D</returns>
	template <typename T>
	std::pair<vec3<T>, vec3<T>> Perps(const vec3<T>& V) {
		T ax = std::abs(V[0]);
		T ay = std::abs(V[1]);
		T az = std::abs(V[2]);

		vec3<T> v1, v2;

		if (ax >= ay && ax >= az) {
			T n1 = std::sqrt(ay * ay + az * az);
			if (n1 == 0) {
				return { vec3(0, 1, 0), vec3<T>(0, 0, 1) };
			}
			v1[1] = V[2] / n1;
			v1[2] = -V[1] / n1;

			v2 = vec3<T>(-(ay * ay + az * az) / V[0], V[1], V[2]);
			v2 /= v2.norm();
		
			return { v1, v2 };
		}
		else if (ay >= ax && ay >= az) {
			T n1 = std::sqrt(ax * ax + az * az);
			if (n1 == 0) {
				return { vec3<T>(1, 0, 0), vec3<T>(0, 0, 1) };
			}
			v1[0] = V[2] / n1;
			v1[2] = -V[0] / n1;

			v2 = vec3<T>(V[0], -(ax * ax + az * az) / V[1], V[2]);
			v2 /= v2.norm();
			return { v1, v2 };
		}
		else {
			T n1 = std::sqrt(ax * ax + ay * ay);
			if (n1 == 0) {
				return { vec3<T>(1, 0, 0), vec3<T>(0, 1, 0) };
			}
			v1[0] = V[1] / n1;
			v1[1] = -V[0] / n1;

			v2 = vec3<T>(V[0], V[1], -(ax * ax + ay * ay) / V[2]);
			v2 /= v2.norm();
			return { v1, v2 };
		}
	}

	/// <summary>
	/// generate out velocity (with weight) in capture.
	/// </summary>
	/// <param name="mp">nuclei mass</param>
	/// <param name="mk">DM mass</param>
	/// <param name="delta_mk">DM delta mass </param>
	/// <param name="dF">form factor function</param>
	/// <param name="VescR"></param>
	/// <param name="VescMin"></param>
	/// <param name="nR"></param>
	/// <param name="TempR"></param>
	/// <param name="G"></param>
	/// <param name="Vdisp"></param>
	/// <param name="mU0"></param>
	/// <param name="pow_r"></param>
	/// <returns>tuple(\vec{V}_{out},r - location, Vdisp - disp velocity at r)</returns>
	template <
		typename Gen_t,
		typename dF_Type,
		typename VescRFuncType,
		typename N_FuncType,
		typename TempRFuncType
		>
	inline MCResult<
		std::tuple<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>, Gen_vt<Gen_t>>,
		Gen_vt<Gen_t>
	> Vout1(
		Gen_t&& G,
		Gen_vt<Gen_t> mi, Gen_vt<Gen_t> mk, 
		Gen_vt<Gen_t> mi_frac, Gen_vt<Gen_t> mk_frac,
		Gen_vt<Gen_t> m_cm, Gen_vt<Gen_t> delta_mk,
		dF_Type dF,
		VescRFuncType const& VescR, 
		Gen_vt<Gen_t> VescMin,
		N_FuncType const& nR, 
		TempRFuncType const& TempR,
		Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> mU0, 
		Gen_vt<Gen_t> pow_r = 1) 
	{
		using T = Gen_vt<Gen_t>;
		T factor = 1;

		//generate radius
		T r_nd = pow(G(), pow_r);//pow(G(),1.0/3.0);
		factor *= (3 * pow_r * pow(r_nd, (3 * pow_r - 1) / pow_r));
		//gain escape velocity from redius
		//if(r_nd < 0.3){
		//    r_nd+=0.0;
		//}
		auto Vesc = VescR(r_nd)*VescMin;

		//random input velocity
		auto VelocityMk = HaloVelocity(G, Vesc, Vdisp, mU0);
		auto V_wimp = VelocityMk.Result;//vec3::PolarCos(sqrt(VelocityMk.Result*VelocityMk.Result+Vesc*Vesc),
		//RandomCos(G),RandomPhi(G));
		factor *= VelocityMk.RemainDensity;


		auto n_nd = nR(r_nd);//TODO n_nd as a function of radius
		factor *= n_nd;
		
		//thermal generation
		auto Input_Vel_gen = Gauss3_BeyondD(G, std::sqrt(TempR(r_nd) / mi),0.8,8);
		
		factor *= Input_Vel_gen.RemainDensity;
		vec3<T> V1 = Input_Vel_gen.Result;

		//Vcm - is a vrlocity of momentum center
		vec3<T> Vcm = (V_wimp * mk_frac + V1 * mi_frac);
		//Nu is input velocity of WIMP in cm coordinatesd
		vec3<T> Vdelta =  (V_wimp - V1);

		//Ecm - is kinetic energy in cm
		auto E_cm = m_cm * Vdelta.squaredNorm() / 2;


		//Energyloss is used in inelastic processes
		//auto EnLoss_mc = dF.EnergyLoss(E_cm - delta_mk);
		T inel_enloss = 0;
		//factor *= EnLoss_mc.RemainDensity;

		// Generating out velocity
		auto Vumk = VuOut(G, Vcm, Vdelta, Vesc, VescMin, 
			mi, mk,mi_frac,mk_frac,m_cm, inel_enloss + delta_mk
		);
		vec3<T> Vu1 = Vumk.Result;
		factor *= Vumk.RemainDensity;

		// q - exchange momentum
		vec3<T> vec_Q = m_cm * (Vdelta - Vu1);
		auto q_2 = vec_Q.squaredNorm();
		
		//
		vec3<T> vec_V_T_inel = (Vdelta + Vu1) * (mi_frac-mk_frac) / 2;
		if(q_2 != 0) 
			vec_V_T_inel += delta_mk * vec_Q / q_2;
		T v_2 = vec_V_T_inel.squaredNorm();
		//////////
		// v_2 = \vec_{v}_{inel T}^{\perp 2}
		// \vec_{v}_{inel T}^{\perp 2} = 1/2( v_x_1 + v_x_2 - v_N_1 - v_N_2)+delta/q^2*q
		// v_x_1 + v_x_2 - v_N_1 - v_N_2 = (nu+nu1)*(1-mk/mi)
		////////// 

		factor *= dF.ScatterFactor(q_2, v_2,inel_enloss);

		using Tuple_t = std::tuple<vec3<T>, T, T>;
		return MCResult<Tuple_t,T>(Tuple_t(Vu1*mi_frac + Vcm, r_nd, Vesc), factor);
	}


	/// <summary>
	/// calculates capture.
	/// </summary>
	/// <param name="G"></param>
	/// <param name="se"></param>
	/// <param name="TempR"></param>
	/// <param name="H_le_sp_type"></param>
	/// <param name="mk"></param>
	/// <param name="dm"></param>
	/// <param name="mi"></param>
	/// <param name="_mp"></param>
	/// <param name="_3p_r"></param>
	/// <param name="mU0"></param>
	/// <param name="Vdisp"></param>
	/// <param name="VescMin"></param>
	/// <param name="Nmk"></param>
	template <
		typename Histo_t,
		typename LEFunc_t,
		typename Rm_type_t,
		typename Gen_t,
		typename TempFunc_t,
		typename VescFunc_t>
	inline std::pair<double,double> CaptureImpl(
		Gen_t&& G,
		ScatterEvent const & se,
		TempFunc_t const& TempR,
		Histo_t &H_le_sp_type,
		LEFunc_t const & LEf,
		Gen_vt<Gen_t> mk, Gen_vt<Gen_t> dm,
		Gen_vt<Gen_t> mi,
		Rm_type_t _3p_r,
		Gen_vt<Gen_t> mU0,
		Gen_vt<Gen_t> Vdisp,
		VescFunc_t const &VescR,
		Gen_vt<Gen_t> VescMin,
		size_t Nmk,
		Gen_vt<Gen_t> weight = 1)
	{
		using T = Gen_vt<Gen_t>;
		T sum = 0;
		T sum2 = 0;
		auto const_fact_rd = weight / Nmk;
		auto m_cm = (mk * mi) / (mk + mi);
		auto mi_frac = (mi) / (mk + mi);
		auto mk_frac = (mk) / (mk + mi);
		auto action = [&](auto&& dF) {
			for (size_t i = 0; i < Nmk; ++i) {
				auto mk_res = Vout1(
					G,mi, mk, mi_frac, mk_frac,m_cm, dm, 
					dF, VescR, VescMin, se.n_e, TempR, Vdisp, mU0, _3p_r);
				vec3<T> v_nd = std::get<0>(mk_res.Result) / VescMin;
				

				T r_nd = std::get<1>(mk_res.Result);
				T  v_esc_nd = std::get<2>(mk_res.Result) / VescMin;

				T E_nd = (v_nd.squaredNorm() - v_esc_nd * v_esc_nd);
				
				if (E_nd < 0) {
					T L_nd = r_nd * sqrt(v_nd.x() * v_nd.x() + v_nd.y() * v_nd.y());
					T l = L_nd / LEf(-E_nd);
					T dens = mk_res.RemainDensity * const_fact_rd;;
					H_le_sp_type.put_force(dens, E_nd, l);
					sum += dens;
					sum2 += dens * mk_res.RemainDensity;
				}/*
				else{
					auto l = L_nd/H.LE_func(E_nd);
				}*/
				/*
			   if(Ewas != H.values[1].values[0]){
				   PVAR(dens);
				   PVAR(v_nd);
				   PVAR(r_nd);
				   PVAR(v_esc_nd);
				   PVAR(E_nd);
				   PVAR(L_nd);
				   PVAR(H.values[1].values[1]);
				   print();
			   }*/
			   /*
		   if(std::isnan(E_nd) or std::isnan(L_nd)){
			   PVAR(r_nd);
			   PVAR(v_nd);
			   PVAR(v_esc_nd);
			   PVAR(v_esc_nd);
		   }*/
			}
		};
		std::visit(action,se.sf);

		return {sum,std::sqrt(sum2* weight - sum * sum )};
	}


	namespace _detail {
		struct _any_ctype {
			template <typename T>
			_any_ctype(T&&) {
				static_assert("metatype");
			}
		};
		template <typename _histo_type>
		auto IsGridObject(_histo_type _inst) ->
			std::tuple<decltype(_inst.Grid, _inst.Values)>;

		auto IsGridObject(_any_ctype _inst) ->
			std::false_type;

		template <typename T>
		struct is_grid_object {
			constexpr static bool value = !std::is_same_v<
				decltype(IsGridObject(std::declval<T>())),
				std::false_type
			>;
		};

	};


	template <
		typename Gen_t,
		typename dF_Type,
		typename VescR_t,
		typename NR_t,
		typename Temp_t
	>
	inline MCResult<
		vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>
	> Vout1_Scatter(
		vec3<Gen_vt<Gen_t>> V_wimp,
		Gen_t&& G,
		Gen_vt<Gen_t> mi, Gen_vt<Gen_t> mk, 
		Gen_vt<Gen_t> mi_frac, Gen_vt<Gen_t> mk_frac,
		Gen_vt<Gen_t> m_cm,Gen_vt<Gen_t> delta_mk,
		dF_Type dF,
		VescR_t const VescR,
		Gen_vt<Gen_t> VescMin,
		NR_t const n_nd,
		Temp_t const TempR
	){
		using T = Gen_vt<Gen_t>;
		T factor = 1;
		factor *= n_nd;

		//thermal generation
		auto Input_Vel_gen = Gauss3_BeyondD(G, std::sqrt(TempR / mi), 0.8, 8);
		factor *= Input_Vel_gen.RemainDensity;
		vec3<T> V1 = Input_Vel_gen.Result;
		//Vcm - is a vrlocity of momentum center
		vec3<T> Vcm = (V_wimp * mk_frac + V1 * mi_frac);
		//Nu is input velocity of WIMP in cm coordinatesd
		vec3<T> Vu =  (V_wimp - V1);

		//Ecm - is kinetic energy in cm
		auto E_cm =m_cm* Vu.squaredNorm() / 2;

		//Energyloss is used in inelastic processes
		//auto EnLoss_mc = dF.EnergyLoss(E_cm - delta_mk);
		T inel_enloss = 0;
		//factor *= EnLoss_mc.RemainDensity;

		// Generating out velocity
		auto Vumk = VuOut_evap(
			G, Vcm, Vu, VescR, VescMin, mi, mk,mi_frac,mk_frac, m_cm, 
			inel_enloss + delta_mk);
		vec3<T> Vu1 = Vumk.Result;
		factor *= Vumk.RemainDensity;

		// q - exchange momentum
		vec3<T> vec_Q = m_cm * (Vu - Vu1);
		auto q_2 = vec_Q.squaredNorm();

		//
		vec3<T> vec_V_T_inel = (Vu + Vu1) * (mi_frac-mk_frac) / 2;
		if (q_2 != 0)
			vec_V_T_inel += delta_mk * vec_Q / q_2;
		T v_2 = vec_V_T_inel.squaredNorm();
		//////////
		// v_2 = \vec_{v}_{inel T}^{\perp 2}
		// \vec_{v}_{inel T}^{\perp 2} = 1/2( v_x_1 + v_x_2 - v_N_1 - v_N_2)+delta/q^2*q
		// v_x_1 + v_x_2 - v_N_1 - v_N_2 = (nu+nu1)*(1-mk/mi)
		////////// 

		factor *= dF.ScatterFactor(q_2, v_2, inel_enloss);

		/**/
		return { vec3<T>(Vu1*mi_frac + Vcm), factor };
	}


	template <
		typename Gen_t,
		typename HistoBank_t,
		typename HistoEvap_t,
		typename LEFunc_t,
		typename BinTrajPool_t,

		typename BinMes_t,

		typename PhiFunc_t,
		typename SFunc_t,
		typename TempFunc_t,
		typename VescFunc_t
	>
	inline void ScatterImpl(
		Gen_t&& G,
		Gen_vt<Gen_t> mk, Gen_vt<Gen_t> dm,
		Gen_vt<Gen_t> mi,
		ScatterEvent const& se,
		HistoBank_t& ScatMat_bank,
		HistoEvap_t&& EvapDistrib,
		bool count_evap,
		LEFunc_t const& LEf,
		Gen_vt<Gen_t> _F2_value,
		const BinTrajPool_t&TrajPoolVec,
		BinMes_t m_mes,
		PhiFunc_t const & Phi,
		SFunc_t const & S_func,

		TempFunc_t const & TempR,
		VescFunc_t const & VescR,
		Gen_vt<Gen_t> VescMin,
		size_t Nmk,
		size_t Nmk_per_traj,
		Gen_vt<Gen_t> weight,
		progress_omp_function<> &m_progress_func
	) {
		typedef 
			bin_traj_p_vt_t<std::decay_t<decltype(TrajPoolVec[0])>>
			Traj_t;
		auto FromFactorVar_Action = [&,&Nir = se.n_e](auto&& dF) {
			const size_t N_in = ScatMat_bank.Grid.size();
			progress_omp_bar m_bar(
				m_progress_func, N_in, std::max((int)(N_in / 1000), 1)
			);
			m_bar.next();
			#pragma	omp parallel for
			auto m_cm = mk * mi / (mk + mi);
			auto mi_frac = mi / (mk + mi);
			auto mk_frac = mk / (mk + mi);
			for (size_t i = 0; i < N_in; ++i) {
				auto IJ = ScatMat_bank.Grid.FromLinear(i);

				auto el_bin = ScatMat_bank.Grid[IJ];

				auto OutHisto = ScatMat_bank.Values[i];
				
				auto * _dbg_ptr = &OutHisto.Values[0];
				(void)_dbg_ptr;
				
				const bin_traj_pool_t<Traj_t> & el_trajpool = TrajPoolVec[i];

				auto TinFunc = make_bin_traj_pool_Tin_func(
					el_trajpool, LEf, _F2_value);
				auto ToutFunc = make_bin_traj_pool_Tout_func(
					el_trajpool, LEf);

				auto m_bin_el_gen = gen_EL(m_mes, el_bin, G, LEf);


				for (size_t nm = 0; nm < Nmk; ++nm) {
					auto [e, l] = m_bin_el_gen();
					auto Lmax = LEf(-e);
					auto Ltmp = l * Lmax;

					auto [u0, u1, theta_max] = 
						TinFunc.template theta1<true>(e,l,Lmax);
					u0 = std::clamp(u0, (decltype(u0))0, (decltype(u0))1);
					u1 = std::clamp(u1, (decltype(u1))0, (decltype(u1))1);
					auto r0 = std::sqrt(u0);
					auto r1 = std::sqrt(u1);
					auto [tp_i, tp_j] = el_trajpool.Grid.pos(e, l);

					auto const& th00 = el_trajpool[{tp_i, tp_j}].theta_tau;
					//prelimenary trajectory
					
					auto Tin_Teheta = TinFunc.tin_theta(e,l);

					auto Tin = Tin_Teheta * theta_max;
					auto Tout = ToutFunc(e, l, Ltmp);

					auto mk_factor = 
						weight * Tin / ((Tin + Tout) * Nmk * Nmk_per_traj);

					for (size_t nt = 0; nt < Nmk_per_traj; ++nt) {
						auto tau = G();
						auto [theta_undim, d_theta_undim] = th00(tau);
						auto theta = theta_undim * theta_max;
						auto d_theta = d_theta_undim * theta_max;

						auto cth2 = std::cos(theta / 2);
						auto sth2 = std::sin(theta / 2);
						Traj_t r = std::sqrt(u0 * cth2 * cth2 + u1 * sth2 * sth2);

						auto renorm_factor =
							d_theta_undim / 
							(2*Tin_Teheta*std::sqrt(S_func(r,r0,r1)));

						Traj_t v2 = downbound(Phi(r) + e, 0);
						auto v = std::sqrt(v2);
						Traj_t vt = upbound(Ltmp / r, v);
						Traj_t vr = ssqrt(v2 - vt * vt);

						vec3<Traj_t > V_in(vt* VescMin, 0, vr* VescMin);
						Traj_t v_esc_nd = VescR(r);
						auto [Vout, factor] =
							Vout1_Scatter(
								V_in, G, mi, mk,mi_frac,mk_frac,m_cm, dm, dF, 
								v_esc_nd* VescMin, VescMin, 
								Nir(r), TempR(r)
							);

						vec3<Traj_t > v_nd = Vout / VescMin;//OK
						Traj_t e_out = (v_nd.squaredNorm() - v_esc_nd * v_esc_nd);
						

						auto final_factor = factor * mk_factor * renorm_factor;
						if (e_out < 0) {
							Traj_t L_nd = r * std::sqrt(v_nd.x() * v_nd.x() + v_nd.y() * v_nd.y());
							Traj_t l_out = (Lmax > 0 ? L_nd / Lmax : 0);
							OutHisto.put_force(final_factor,e_out, l_out);
						}
						else if (count_evap) {
							EvapDistrib.Values[i] += final_factor;
						}
					}
				}
				m_bar.next();
			}
		};
		std::visit(FromFactorVar_Action, se.sf);
	}


};

