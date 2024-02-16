#ifndef DYNAMICS_HPP
#define DYNAMICS_HPP
#include "../form_factors.hpp"
#include "../core.hpp"
#include <grob/grid.hpp>
#include <grob/grid_objects.hpp>



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
	inline MCResult<vec3<Gen_vt<Gen_t>>> Gauss3_BeyondD(
		Gen_t&& G, Gen_vt<Gen_t> Vdisp, 
		Gen_vt<Gen_t> p = 0.8,
		Gen_vt<Gen_t> max_xi = 8)
	{
		using mvec3 = vec3<Gen_vt<Gen_t>>;
		if (G() < p) {
			auto V1 = Vdisp * sqrt(-2 * log(1 - G()));
			auto V2 = Vdisp * sqrt(-2 * log(1 - G()));
			auto phi1 = RandomPhi(G);
			auto phi2 = RandomPhi(G);
			return { mvec3(V1 * cos(phi1), V1 * sin(phi1), V2 * cos(phi2)),1 };
		}
		else {
			auto r = std::cbrt(G()) * max_xi;
			auto phi = RandomPhi(G);
			auto cos_Th = RandomCos(G);
			auto sin_Th = std::sqrt(1 - cos_Th* cos_Th);
			auto r_xy = r * sin_Th;
			auto fac = exp(-r*r / 2) / std::pow(2 * M_PI, 1.5);
			return { 
				Vdisp * mvec3(r_xy * std::sin(phi),r_xy*std::sin(phi),r*cos_Th),
				fac};
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
	inline MCResult<vec3<Gen_vt<Gen_t>>> HaloVelocity(Gen_t&& G,Gen_vt<Gen_t> VescTmp,
		Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> mU0) {

		auto ksi = sqrt(-2 * log(E0I1_G(G())));
		auto phi = RandomPhi(G)/2;


		auto sinPhi = sin(phi);
		auto cosPhi = cos(phi);

		auto u0 = mU0 / Vdisp;
		auto ve = VescTmp / Vdisp;

		auto u = sqrt(u0 * u0 + ksi * ksi + 2 * u0 * ksi * cosPhi);


		auto sinTheta = (u != 0.0 ? ksi * sinPhi / u : 0);
		auto v = sqrt(u * u + ve * ve);
		
		auto n = RandomNvec(G);
		
		return MCResult<vec3<Gen_vt<Gen_t>>>(n * (v * Vdisp), sinTheta * v * sqrt(M_PI / 2));
	}
	template <class Gen_t>
	/*MK generator of input velocity*/
	inline MCResult<vec3< Gen_vt<Gen_t>>> HaloVelocityConstrained(
		Gen_t&& G, Gen_vt<Gen_t> VescTmp, Gen_vt<Gen_t> Vmin, Gen_vt<Gen_t> Vmax,
		Gen_vt<Gen_t>  Vdisp, Gen_vt<Gen_t>  mU0)
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

		number_t rd_th = (max_theta - min_theta) / M_PI;

		number_t theta = min_theta + (max_theta - min_theta) * G();

		number_t ve = VescTmp / Vdisp;

		number_t u = sqrt(u0 * u0 + ksi * ksi + 2 * u0 * ksi * cos(theta));

		number_t sinTheta = (u != 0.0 ? ksi * sin(theta) / u : 0);

		auto v = sqrt(u * u + ve * ve);
		auto n = RandomNvec(G);

		return { n * (v * Vdisp), ksi_rd * rd_th * sinTheta * v * sqrt(M_PI / 2) };
	}

	template <class Gen_t,typename VectorT>
	/*MK generator of output nu'*/
	inline MCResult<vec3<Gen_vt<Gen_t>>> NuOut(Gen_t && G, const  VectorT& Vcm, const VectorT& Nu,
		Gen_vt<Gen_t> Vesc, Gen_vt<Gen_t> VescMin, 
		Gen_vt<Gen_t>  mp, Gen_vt<Gen_t>  mk, Gen_vt<Gen_t> deltaE = 0) {

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

		number_t Nu1_squared = Nu.quad() - deltaE * 2 * mp / (mk * (mp + mk));
		if (Nu1_squared <= 0.0)
			return MC::MCResult<vec3_t>(vec3_t(0, 0, 0), 0);

		number_t Nu1 = sqrt(Nu1_squared);

		number_t cosTh1max = (Vesc * Vesc - Nu1_squared - VcmN * VcmN) / (2 * VcmN * Nu1);

		if (!(cosTh1max > -1))
			return MC::MCResult<vec3_t>(vec3_t(0, 0, 0), 0);
		else if (cosTh1max >= 1) {
			cosTh1max = 1;
		}

		number_t cosTh1 = (1 + cosTh1max) * G() - 1;
		number_t sinTh1 = sqrt(1.0 - cosTh1 * cosTh1);
		number_t phi1 = RandomPhi(G);

		const vec3_t vNu1 = Nu1 * (n_v * cosTh1 + n_1 * sinTh1 * cos(phi1) + n_2 * sinTh1 * sin(phi1));
		return MC::MCResult<vec3_t>(vNu1, (1 + cosTh1max)/2 * Nu1 / VescMin);
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
		typename ScatterFuncType,
		typename VescRFuncType,
		typename N_FuncType,
		typename TempRFuncType
		>
	inline MCResult<std::tuple<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>, Gen_vt<Gen_t>>, Gen_vt<Gen_t>> Vout1(
		Gen_t&& G,
		Gen_vt<Gen_t> mp, Gen_vt<Gen_t> mk, Gen_vt<Gen_t> delta_mk, 
		dF_Type dF,
		VescRFuncType const& VescR, 
		Gen_vt<Gen_t> VescMin,
		N_FuncType const& nR, 
		TempRFuncType const& TempR,
		Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> mU0, 
		Gen_vt<Gen_t> pow_r = 1) {

		Gen_vt<Gen_t> factor = 1;

		//generate radius
		Gen_vt<Gen_t> r_nd = pow(G(), pow_r);//pow(G(),1.0/3.0);
		factor *= (3 * pow_r * pow(r_nd, (3 * pow_r - 1) / pow_r));
		//gain escape velocity from redius
		//if(r_nd < 0.3){
		//    r_nd+=0.0;
		//}
		auto Vesc = VescR(r_nd);

		//random input velocity
		auto VelocityMk = HaloVelocity(G, Vesc, Vdisp, mU0);
		auto V_wimp = VelocityMk.Result;//vec3::PolarCos(sqrt(VelocityMk.Result*VelocityMk.Result+Vesc*Vesc),
		//RandomCos(G),RandomPhi(G));
		factor *= VelocityMk.RemainDensity;


		auto n_nd = nR(r_nd);//TODO n_nd as a function of radius
		factor *= n_nd;
		auto [V1, f_mlt] = Gauss3_BeyondD(G, std::sqrt(TempR(r_nd) / mp));//TODO: add thermal distribution of nuclei velocity
		factor *= f_mlt;

		//Vcm - is a vrlocity of momentum center
		auto Vcm = (V_wimp * mk + V1 * mp) / (mp + mk);
		//Nu is input velocity of WIMP in cm coordinatesd
		auto Nu = mp / (mp + mk) * (V_wimp - V1);

		//Ecm - is kinetic energy in cm
		auto E_cm = mk * (mp + mk) / mp * Nu.squaredNorm() / 2;



		auto EnLoss = dF.EnergyLoss(E_cm + delta_mk);
		factor *= EnLoss.RemainDensity;

		// Generating out velocity
		auto Numk = NuOut(G, Vcm, Nu, Vesc, VescMin, mp, mk, EnLoss.Result - delta_mk);
		vec3<Gen_vt<Gen_t>> Nu1 = Numk.Result;
		factor *= Numk.RemainDensity;

		// q - exchange momentum
		auto q = mk * (Nu - Nu1).norm();

		factor *= dF.ScatterFactor(q, EnLoss.Result);

		//factor from matrix element
		factor *= PhiFactor(q);


#ifdef VOUT_DEBUG
		if ((Nu1 + Vcm).norm() > Vesc * (1. + 1e-8) && factor != 0) {
			std::cout << "V out of bund:\n" +
				SVAR((Nu1 + Vcm).norm()) + "\n" +
				SVAR(factor) + "\n" +
				SVAR(Vesc) + "\n" +
				SVAR(Nu1 + Vcm) + "\n" +
				SVAR(V1.norm()) + "\n" +
				SVAR(V_wimp.norm()) + "\n\n";
		}
		/**/
		if (factor != 0) {
			PVAR(mp);
			PVAR(mk);
			PVAR(Vesc);
			PVAR(V_wimp.norm());
			PVAR(Vcm.norm());
			PVAR(Nu.norm());
			PVAR((Nu1 + Vcm).norm());
			PVAR((((Nu1 + Vcm).quad()) - Vesc * Vesc) / 4.227e-6);
			std::cin.get();
		}
#endif
		/**/

		return MC::MCResult<std::tuple<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>, Gen_vt<Gen_t>>>(
			std::tuple<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>, Gen_vt<Gen_t>>(Nu1 + Vcm, r_nd, Vesc),
			factor);
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
		typename Rm_type_t,
		typename Gen_t,
		typename TempFunc_t>
	inline std::pair<double,double> CaptureImpl(
		Gen_t&& G,
		ScatterEvent const & se,
		TempFunc_t TempR,
		Histo_t &H_le_sp_type,
		float mk,float dm,
		float mi,
		Rm_type_t _3p_r,
		Gen_vt<Gen_t> mU0,
		Gen_vt<Gen_t> Vdisp,
		Gen_vt<Gen_t> VescMin,
		size_t Nmk,
		Gen_vt<Gen_t> weight = 1)
	{
		double sum = 0;
		double sum2 = 0;
		auto const_fact_rd = weight / Nmk;

		auto action = [&](auto&& dF) {
			for (size_t i = 0; i < Nmk; ++i) {
				auto mk_res = Vout1(G,mi, mk, dm, dF, VescR, VescMin, se.n_e, TempR, Vdisp, mU0, _3p_r);
				auto v_nd = std::get<0>(mk_res.Result) / VescMin;
				auto r_nd = std::get<1>(mk_res.Result);
				auto v_esc_nd = std::get<2>(mk_res.Result) / VescMin;

				auto E_nd = (v_nd * v_nd - v_esc_nd * v_esc_nd);
				auto L_nd = r_nd * sqrt(v_nd.x() * v_nd.x() + v_nd.y() * v_nd.y());
				auto dens = mk_res.RemainDensity * const_fact_rd;
				//double Ewas = H.values[1].values[0];
				/*
				auto l = L_nd/H.LE_func(E_nd);
				auto [b,MI] = H.Histo.Grid.spos(E_nd,l);
				auto i_h = H.Histo.Grid.LinearIndex(MI);
				if(b){
					auto Elem = H.Histo.Grid[MI];
					auto H_i = H.Histo.Values[i_h];
				}*/
				/*if(E_nd < -4.5){
					dens = dens +0.0;
					print("E<-4.5");
				}*/
				if (H_le_sp_type.putValue(dens, E_nd, L_nd)) {
					/*
					auto H_i_1 = H.Histo.Values[i_h];
					*/
					sum += dens;
					sum2 += dens * dens;
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

		return {sum,std::sqrt(sum2- sum * sum )};
	}
};



#endif//DYNAMICS_HPP