#pragma once
namespace evdm {

	namespace _detail {
		template <typename T>
		struct same_type_base {
			typedef T type;
		};
	};
	

	template <typename T>
	using same_t = typename _detail::same_type_base<T>::type;
};