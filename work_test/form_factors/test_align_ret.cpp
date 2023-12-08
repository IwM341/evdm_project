#include <iostream>
#include <variant>
#include <tuple>
#include <array>
template <typename T>
struct self_type{
	typedef T type;
};

template <typename variant_type,typename InitializerLambda,size_t...I>
constexpr inline auto make_variant_impl(size_t i,InitializerLambda && Init,std::index_sequence<I...>){
	constexpr size_t variantSize = sizeof...(I);
	std::array<variant_type (*)(InitializerLambda &&),variantSize> Funcs {
		[](InitializerLambda && _Init){
			typedef typename std::variant_alternative<I,variant_type>::type type;
			return variant_type (_Init(self_type<type>{}));
		}...
	};
	return Funcs[i](std::forward<InitializerLambda>(Init));
}
template <typename variant_type,typename InitializerLambda>
auto make_variant(size_t index,InitializerLambda && Init){
	return make_variant_impl<variant_type>(index,Init,
	std::make_index_sequence<std::variant_size<variant_type>::value>{});
}

struct alignas(4*alignof(float))  S16{
	float arr[4];
};


struct alignas(16*alignof(float))  S32{
	float arr[16];
};

struct alignas(32*alignof(float))  S64{
	float arr[32];
};

template <typename T>
T f(){
	return T{};
}

int main(void){
	
	std::variant<S64,S32,S16> x = 
		make_variant<std::variant<S64,S32,S16>>
			(1,[](auto t){return typename decltype(t)::type {};});

	return 0;
}