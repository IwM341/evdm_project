#pragma once
#include "../form_factors.hpp"
#include "../measure.hpp"
#include <grob/grid.hpp>
#include <grob/grid_objects.hpp>
#include <numbers>
#include "../utils/mc.hpp"

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
		ScatterEvent(size_t arr_size, Array_t&& N_values, FormFactor_t _ffact) :
			n_e(
				grob::GridUniform<float>(0, 1, arr_size), std::vector<float>(arr_size)
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

	template <typename Gen_t>
	vec3<Gen_vt<Gen_t>> RandomNvec(Gen_t& G) {
		auto cosT = RandomCos(G);
		auto sinT = std::sqrt(1 - cosT * cosT);
		auto phi = RandomPhi(G);

		return { sinT * cos(phi),sinT * sin(phi),cosT };
	}



	template <typename T>
	T _cube_f(T x) {
		return x * x * x;
	}

	
	template <class Gen_t>
	inline MCResult < Gen_vt<Gen_t>, Gen_vt<Gen_t> > Gauss3_Soft8_abs(
		Gen_t& G, Gen_vt<Gen_t> Vdisp
	) {
		constexpr Gen_vt<Gen_t> a = 0.16;
		auto xi = E0I1_G(G);
		auto y = std::cbrt(xi);
		auto x = (a * y) / (1 + (a - 1) * y);

		auto Zx = a + (1 - a) * x;
		auto deriv_xi_inv = Zx * Zx / (3 * a * y * y);

		auto x8 = 8 * x;

		constexpr Gen_vt < Gen_t> _bf_exp =
			8 * std::numbers::inv_sqrtpi * std::numbers::sqrt2;
		auto target = _bf_exp * (x8 * x8) * std::exp(-x8 * x8 / 2);
		return { (x8 * Vdisp), target * deriv_xi_inv };
	}

	template <class Gen_t>
	inline MCResult < vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t> > Gauss3_Soft8(
		Gen_t& G, Gen_vt<Gen_t> Vdisp
	) {
		constexpr Gen_vt<Gen_t> a = 0.16;
		auto xi = E0I1_G(G);
		auto y = std::cbrt(xi);
		auto x = (a * y) / (1 + (a - 1) * y);

		auto Zx = a + (1 - a) * x;
		auto deriv_xi_inv = Zx * Zx / (3 * a * y * y);

		auto x8 = 8 * x;

		constexpr Gen_vt < Gen_t> _bf_exp = 
			8 * std::numbers::inv_sqrtpi* std::numbers::sqrt2;
		auto target = _bf_exp * (x8 * x8) * std::exp(-x8 * x8 / 2);
		return { RandomNvec(G) * (x8* Vdisp), target * deriv_xi_inv };
	}

	template <class Gen_t>
	inline MCResult<Gen_vt<Gen_t>, Gen_vt<Gen_t>> Gauss3_Soft8_abs(
		Gen_t& G, Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> Vmin) {
		constexpr Gen_vt<Gen_t> a = 0.16;

		auto xi_min = _cube_f(Vmin / (a * Vdisp + (1 - a) * Vmin));

		auto factor = (1 - xi_min);

		auto xi = xi_min + factor * G();
		auto y = std::cbrt(xi);
		auto x = (a * y) / (1 + (a - 1) * y);

		auto Zx = a + (1 - a) * x;
		auto deriv_xi_inv = Zx * Zx / (3 * a * y * y);

		auto x8 = 8 * x;

		constexpr Gen_vt < Gen_t> _bf_exp =
			8 * std::numbers::inv_sqrtpi * std::numbers::sqrt2;

		auto target = _bf_exp * (x8 * x8) * std::exp(-x8 * x8 / 2);
		return { (x8 * Vdisp), target * factor * deriv_xi_inv };
	}

	template <class Gen_t>
	inline MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> Gauss3_Soft8(
		Gen_t& G, Gen_vt<Gen_t> Vdisp, Gen_vt<Gen_t> Vmin) {
		constexpr Gen_vt<Gen_t> a = 0.16;

		auto xi_min = _cube_f(Vmin / (a* Vdisp + (1 - a) * Vmin));

		auto factor = (1 - xi_min);

		auto xi = xi_min + factor*G();
		auto y = std::cbrt(xi);
		auto x = (a * y) / (1 + (a - 1) * y);

		auto Zx = a + (1 - a) * x;
		auto deriv_xi_inv = Zx * Zx / (3 * a * y * y);

		auto x8 = 8 * x;

		constexpr Gen_vt < Gen_t> _bf_exp =
			8 * std::numbers::inv_sqrtpi * std::numbers::sqrt2;

		auto target = _bf_exp * (x8 * x8) * std::exp(-x8 * x8 / 2);
		return { RandomNvec(G) * (x8* Vdisp), target * factor*deriv_xi_inv };
	}


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
		Gen_t& G, Gen_vt<Gen_t> Vdisp,
		Gen_vt<Gen_t> p = (Gen_vt < Gen_t>)0.8,
		Gen_vt<Gen_t> max_xi = (Gen_vt < Gen_t>)8)
	{
		using T = Gen_vt<Gen_t>;
		using mvec3 = vec3<T>;
		if (G() < p) {
			auto V1 = Gauss2Norm(G, Vdisp);
			auto V2 = Gauss2Norm(G, Vdisp);
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
			T sin_Th = std::sqrt(1 - cos_Th * cos_Th);
			T r_xy = r * sin_Th;
			constexpr T e_fac = (std::numbers::sqrt2*std::numbers::inv_sqrtpi/ 3);
			T fac = exp(-r * r / 2)*(max_xi*max_xi*max_xi* e_fac);
			return {
				Vdisp * mvec3(r_xy * std::sin(phi),r_xy * std::sin(phi),r * cos_Th),
				fac
			};
		}
	}
	

	struct ThermGaussGenerator_NoTherm {
		inline static bool detect(std::string_view pot_name) {
			return pot_name == "notherm";
		}

		template <class Gen_t>
		inline static MCResult<Gen_vt<Gen_t>, Gen_vt<Gen_t>> gen_abs(
			Gen_t & G, Gen_vt<Gen_t> Delta_x_2_div_mu, Gen_vt<Gen_t> Vmk, Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			return { 0, 1 };
		}

		template <class Gen_t>
		inline static MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> gen(
			Gen_t& G, Gen_vt<Gen_t> Delta_x_2_div_mu, Gen_vt<Gen_t> Vmk, Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			return {vec3<Gen_vt<Gen_t >>(0,0,0), 1};
		}
		
	};

	struct ThermGaussGenerator_Naive{
		inline static bool detect(std::string_view pot_name) {
			return pot_name == "naive";
		}
		
		template <class Gen_t>
		inline static MCResult<Gen_vt<Gen_t>, Gen_vt<Gen_t>> gen_abs(
			Gen_t& G, Gen_vt<Gen_t> Delta_x_2_div_mu, Gen_vt<Gen_t> Vmk,
			Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			typedef Gen_vt<Gen_t> T;
			auto Vdisp = std::sqrt(Therm / Mtarget);

			auto V12 = Gauss2Norm2(G, Vdisp);
			auto V22 = Gauss2Norm2(G, Vdisp);
			auto phi2 = RandomPhi(G);
			auto cos2 = std::cos(phi2);
			return {
				std::sqrt(V12+ V22*cos2 * cos2),
				(T)1
			};
		}
		template <class Gen_t>
		inline static MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> gen(
			Gen_t& G, Gen_vt<Gen_t> Delta_x_2_div_mu, Gen_vt<Gen_t> Vmk,
			Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			typedef Gen_vt<Gen_t> T;
			auto Vdisp = std::sqrt(Therm / Mtarget);

			auto V1 = Gauss2Norm(G, Vdisp);
			auto V2 = Gauss2Norm(G, Vdisp);
			auto phi1 = RandomPhi(G);
			auto phi2 = RandomPhi(G);
			return {
				vec3<T>(V1 * std::cos(phi1), V1 * std::sin(phi1), V2 * std::cos(phi2)),
				(T)1
			};
		}
	};

	struct ThermGaussGenerator_Soft8 {
		inline static bool detect(std::string_view pot_name) {
			return pot_name == "soft";
		}

		template <class Gen_t>
		inline static MCResult<Gen_vt<Gen_t>, Gen_vt<Gen_t>> gen_abs(
			Gen_t & G, Gen_vt<Gen_t> Delta_x_2_div_mu,
			Gen_vt<Gen_t> Vmk, Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			typedef Gen_vt<Gen_t> T;
			auto Vdisp = std::sqrt(Therm / Mtarget);
			return Gauss3_Soft8_abs(G, Vdisp);
		}

		template <class Gen_t>
		inline static MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> gen(
			Gen_t& G, Gen_vt<Gen_t> Delta_x_2_div_mu, 
			Gen_vt<Gen_t> Vmk, Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			typedef Gen_vt<Gen_t> T;
			auto Vdisp = std::sqrt(Therm / Mtarget);
			return Gauss3_Soft8(G, Vdisp);
		}
	};
	struct ThermGaussGenerator_Soft8_Treshold {

		inline static bool detect(std::string_view pot_name) {
			return pot_name == "soft_tresh";
		}

		template <class Gen_t>
		inline static MCResult<Gen_vt<Gen_t>, Gen_vt<Gen_t>> gen_abs(
			Gen_t & G, Gen_vt<Gen_t> Delta_x_2_div_mu,
			Gen_vt<Gen_t> Vmk, Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			typedef Gen_vt<Gen_t> T;
			auto Vdisp = std::sqrt(Therm / Mtarget);
			auto Vtresh = std::clamp(
				std::sqrt(Delta_x_2_div_mu) - Vmk,
				(T)0, (T)(8 * Vdisp)
			);
			return Gauss3_Soft8_abs(G, Vdisp, Vtresh);
		}

		template <class Gen_t>
		inline static MCResult<vec3<Gen_vt<Gen_t>>, Gen_vt<Gen_t>> gen(
			Gen_t& G, Gen_vt<Gen_t> Delta_x_2_div_mu, 
			Gen_vt<Gen_t> Vmk, Gen_vt<Gen_t> Therm, Gen_vt<Gen_t> Mtarget)
		{
			typedef Gen_vt<Gen_t> T;
			auto Vdisp = std::sqrt(Therm / Mtarget);
			auto Vtresh = std::clamp(
				std::sqrt(Delta_x_2_div_mu) - Vmk, 
				(T)0, (T)(8 * Vdisp)
			);
			return Gauss3_Soft8(G, Vdisp, Vtresh);
		}
	};


	typedef std::variant<
		ThermGaussGenerator_NoTherm,
		ThermGaussGenerator_Naive,
		ThermGaussGenerator_Soft8,
		ThermGaussGenerator_Soft8_Treshold
	> ThermGaussGenerator_Vasriant_t;



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
				return { vec3<T>(0, 1, 0), vec3<T>(0, 0, 1) };
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
	template <typename T>
	std::tuple<vec3<T>, vec3<T>, vec3<T>> Basis3(const vec3<T>& ZVec) {
		T x2 = ZVec[0] * ZVec[0];
		T y2 = ZVec[1] * ZVec[1];
		T z2 = ZVec[2] * ZVec[2];
		T norm = std::sqrt(x2 + y2 + z2);
		if (norm == 0) {
			return { {1,0,0},{0,1,0},{0,0,1} };
		}
		T norm_inv = 1 / norm;
		vec3<T> Nv = ZVec * norm_inv;
		if (x2 >= y2 && x2 >= z2) {
			T cosTheta = Nv[0];
			T sinTheta = std::sqrt(Nv[1] * Nv[1] + Nv[2] * Nv[2]);
			if (sinTheta == 0) {
				return { {0,1,0},{0,0,1},{cosTheta,0,0} };
			}
			T norm12_inv = 1 / sinTheta;
			T sphi = Nv[2] * norm12_inv;
			T cphi = Nv[1] * norm12_inv;
			return { {-sinTheta,cosTheta * cphi,cosTheta * sphi},
					{0,-sphi,cphi},
					Nv
			};
		}
		else if (y2 >= z2 && y2 >= x2) {
			T cosTheta = Nv[1];
			T sinTheta = std::sqrt(Nv[0] * Nv[0] + Nv[2] * Nv[2]);
			if (sinTheta == 0) {
				return { {1,0,0},{0,0,1},{0,cosTheta,0} };
			}
			T norm12_inv = 1 / sinTheta;
			T sphi = Nv[2] * norm12_inv;
			T cphi = Nv[0] * norm12_inv;
			return { {cosTheta * cphi,-sinTheta,cosTheta * sphi},
					{-sphi,0,cphi},
					Nv
			};
		}
		else {
			T cosTheta = Nv[2];
			T sinTheta = std::sqrt(Nv[0] * Nv[0] + Nv[1] * Nv[1]);
			if (sinTheta == 0) {
				return { {1,0,0},{0,1,0},{0,0,cosTheta} };
			}
			T norm12_inv = 1 / sinTheta;
			T sphi = Nv[1] * norm12_inv;
			T cphi = Nv[0] * norm12_inv;
			return { {cosTheta * cphi,cosTheta * sphi,-sinTheta},
					{-sphi,cphi,0},
					Nv
			};
		}
	}
	template <typename T,typename GenType>
	vec3<T> GenVecCos(GenType& G, const vec3<T>& Vref, same_t<T> cosThetaMin = -1, same_t<T> cosThetaMax = 1) {
		auto [e1, e2, ez] = Basis3(Vref);
		T cosTheta = cosThetaMin + (cosThetaMax - cosThetaMin) * G();
		T sinTheta = std::sqrt(1 - cosTheta * cosTheta);
		T phi = RandomPhi(G);
		return e1 * (std::cos(phi) * sinTheta) +
			e2 * (std::sin(phi) * sinTheta) +
			ez * cosTheta;
	}
};
