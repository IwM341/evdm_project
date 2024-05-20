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
	vec3<Gen_vt<Gen_t>> RandomNvec(Gen_t&& G) {
		auto cosT = RandomCos(G);
		auto sinT = std::sqrt(1 - cosT * cosT);
		auto phi = RandomPhi(G);

		return { sinT * cos(phi),sinT * sin(phi),cosT };
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
		Gen_t&& G, Gen_vt<Gen_t> Vdisp,
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

};
