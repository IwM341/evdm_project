#ifndef MATH_EXTERNAL_HPP
#define MATH_EXTERNAL_HPP
#include <cmath>
#include <numbers>

namespace evdm {

	template <typename T>
	inline constexpr T ssqrt(T x){
		return std::sqrt((x > 0 ? x : (T)0));
	}

	template <typename T>
	inline constexpr T unnan(T x) {
		return !std::isnan(x) ? x : 0;
	}

	template <typename T,typename U>
	inline constexpr T upbound(T  x,U up_x ) {
		return x< up_x ? x : up_x;
	}
	template <typename T, typename U>
	inline constexpr T downbound(T  x, U down_x) {
		return x > down_x ? x : down_x;
	}

	template <typename T>
	constexpr T pi = (T)std::numbers::pi;
};

#endif