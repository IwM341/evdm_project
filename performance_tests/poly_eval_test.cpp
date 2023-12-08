#include <iostream>
#include <string>

#include <time.h>
#include <chrono>
#include <cstdlib>
#include <experimental/simd>
#include <array>
#include <memory>

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;



template <size_t Deg>
struct coeff_evaulator1{
	template <typename T,typename V>
	static inline constexpr T eval(const T * data,V const &x){
		V sum = 0;
		for(int i=Deg;i>=0;--i){
			sum *= x;
			sum += data[i];
		}
		return sum;
	}
};

namespace stdx = std::experimental;
template <typename T,typename Abi>
inline constexpr T eval2deg_poly(const stdx::simd<T,Abi> & even_coeffs,const stdx::simd<T,Abi> & odd_coeffs,T x){
    constexpr size_t N =  stdx::simd<T,Abi>::size();
    if constexpr(N > 1){
        auto result_tp = stdx::split<N/2,N/2>(even_coeffs + x*odd_coeffs);
        return eval2deg_poly(std::get<0>(result_tp),std::get<1>(result_tp),x*x);
    } else {
        return (even_coeffs + x*odd_coeffs)[0];
    }
}
template <typename T,typename Abi>
inline constexpr T eval2deg_poly_v2(const stdx::simd<T,Abi> & even_coeffs,const stdx::simd<T,Abi> & odd_coeffs,T x){
    constexpr size_t N =  stdx::simd<T,Abi>::size();
    if constexpr(N > 1){
        auto result = (even_coeffs + x*odd_coeffs);
		constexpr size_t N = decltype(result)::size();
		if constexpr(N <= 4){
			T sum = result[N-1];
			for(int i=N-2;i>=0;--i){
				sum *=x;
				sum+=result[i];
			}
			return sum;
		} else{
			auto result_tp = stdx::split<N/2,N/2>(result);
        	return eval2deg_poly(std::get<0>(result_tp),std::get<1>(result_tp),x*x);
		}
    } else {
        return (even_coeffs + x*odd_coeffs)[0];
    }
}


template <size_t size>
float __attribute__ ((noinline)) eval1(float * __poly,float x){
	float * poly = std::assume_aligned<4*sizeof(float)>(__poly);
	return coeff_evaulator1<size-1>::eval(poly,x);
}

template <size_t size>
float __attribute__ ((noinline)) eval(float *__poly,float x){
	float * poly = std::assume_aligned<4*sizeof(float)>(__poly);
	constexpr size_t partial_size = size/2; 
	using CF = stdx::fixed_size_simd<float, partial_size>;
	return eval2deg_poly(CF(poly,stdx::vector_aligned),CF(poly+partial_size,stdx::vector_aligned),x);
}


template <size_t size>
inline float eval2_impl(float * __poly,float x){
	float * poly = std::assume_aligned<4*sizeof(float)>(__poly);
	constexpr size_t partial_size = size/2; 
	if constexpr(size > 2){
		constexpr size_t partial_size = size/2; 
		alignas(partial_size*sizeof(float)) float res[partial_size];
		for(size_t i=0;i<partial_size;++i){
			res[i] = poly[i] + x*poly[partial_size+i];
			x = 0.0f + x*x;
		}
		return eval2_impl<partial_size>(res,x);
	} else{
		if constexpr (size==2)
			return poly[0] + x*poly[1];
		else if (size == 1)
			return poly[0];
		else return 0;
	}
}

template <size_t size>
float __attribute__ ((noinline)) eval2(float * __poly,float x){
	float * poly = std::assume_aligned<4*sizeof(float)>(__poly);
	constexpr size_t partial_size = size/2; 
	using CF = stdx::fixed_size_simd<float, partial_size>;
	return eval2deg_poly_v2(CF(poly,stdx::vector_aligned),CF(poly+partial_size,stdx::vector_aligned),x);
}

float next_f(float x){
	float y = x + 3.14;
	return y - floor(y);
}

static float rand_init = 12.23;
float rnd(){
	rand_init = next_f(rand_init);
	return rand_init;
}
#define RAND_COEFFS
#undef RAND_COEFFS
int main(int argc,char ** argv){
	//int cnt = std::stoi(argv[1]);
	
	constexpr static size_t Deg = 7;
	size_t Ntest = 100000000;

	alignas(8*alignof(float)) float coeffs[Deg+1] = {1,2,3,4,5,6,7,8}; //needs 3 operation op(a,b,x) = a + b*x
	alignas(8*alignof(float)) std::array<float,8> cf{1,5,3,7,2,6,4,8}; //needs 1 parallel and 1 single operation op(a,b,x)
    using CF = stdx::fixed_size_simd<float, 2>;
    //cout << eval2deg_poly(CF(cf.data(),stdx::vector_aligned),CF(cf.data()+2,stdx::vector_aligned),-1.0f) <<endl;

	float sum = 0;
	srand(time(0));

	rand_init = rand()/(RAND_MAX+0.0f);


	for(size_t j = 0;j<8;++j){
		coeffs[j] = rnd();
		cf[j] = rnd();
	}

	auto t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		#ifdef RAND_COEFFS
		for(size_t j = 0;j<8;++j){
			coeffs[j] = next_f(coeffs[j]);
		}
		#endif
		sum += rnd();
	}
    auto t1 = Time::now();
    fsec fs = t1 - t0;
    ms d = std::chrono::duration_cast<ms>(fs);
    std::cout << "default time = " << d.count() << "ms\n";
	std::cout << "sum = " << sum/Ntest <<std::endl << std::endl;


	

	auto def_time = d.count();

	sum = 0;
	t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		#ifdef RAND_COEFFS
		for(size_t j = 0;j<8;++j){
			cf[j] = next_f(cf[j]);
		}
		#endif
		sum += rnd();
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "eval zero time = " << d.count() << "ms\n";
	std::cout << "sum = " << sum/Ntest <<std::endl << std::endl;

	auto def_time1 = d.count();

	sum = 0;
	t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		#ifdef RAND_COEFFS
		for(size_t j = 0;j<8;++j){
			cf[j] = next_f(cf[j]);
		}
		#endif
		sum += eval<8>(cf.data(),rnd());
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "eval parallel time = " << d.count() << "ms\n";
	std::cout << "delta time = " << d.count() - def_time1  << "ms\n";
	std::cout << "sum = " << sum/Ntest <<std::endl << std::endl;

	


	sum = 0;
	t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		#ifdef RAND_COEFFS
		for(size_t j = 0;j<8;++j){
			coeffs[j] = next_f(coeffs[j]);
		}
		#endif
		sum += eval1<8>(coeffs,rnd());;
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "eval not parallel time = " << d.count() << "ms\n";
	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	std::cout << "sum = " << sum/Ntest << std::endl << std::endl;

	sum = 0;
	t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		#ifdef RAND_COEFFS
		for(size_t j = 0;j<8;++j){
			coeffs[j] = next_f(coeffs[j]);
		}
		#endif
		sum += eval2<8>(cf.data(),rnd());;
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "not simd parallel time = " << d.count() << "ms\n";
	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	std::cout << "sum = " << sum/Ntest << std::endl << std::endl;
	
	return 0;
}