#include <iostream>
#include <string>

#include <time.h>
#include <chrono>
#include <cstdlib>
#include <experimental/simd>
#include <array>
#include <memory>
#include <malloc.h>

typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;


float __attribute__ ((noinline)) sum_simd(const float * array,const size_t N){
	const float * poly = std::assume_aligned<16*sizeof(float)>(array);
	alignas(16*sizeof(float)) float tmp[16];
	std::memcpy(tmp,poly,16*sizeof(float));

	for(size_t i=1;i<N/16;++i){
		for(size_t j=0;j< 16;++j){
			tmp[j] += poly[i*16+j];
		}
	}
	float sum = 0;
	for(size_t j=0;j< 16;++j){
		sum += tmp[j];
	}
	return sum;
}

float __attribute__ ((noinline)) sum(const float *array,const size_t N){
	float sum = 0;
	for(size_t i=0;i<N;++i){
		sum += array[i];
	}
	return sum;
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
	
	const size_t N = 1000000000;
	float * array = (float * ) _aligned_malloc(N*sizeof(float),16*sizeof(float));
	const size_t Ntest = 1; 
	auto t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		for(size_t i=0;i<N;++i){
			array[i] = rnd();
		}
	}
    auto t1 = Time::now();
    fsec fs = t1 - t0;
    ms d = std::chrono::duration_cast<ms>(fs);
    std::cout << "default time = " << d.count() << "ms\n";

	

	auto def_time = d.count();

	float result = 0;
	t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		result = sum(array,N);
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "eval simple time = " << d.count() << "ms\n";
	std::cout << "result = " << result/N <<std::endl << std::endl;

	auto def_time1 = d.count();

	t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		result = sum_simd(array,N);
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "eval parallel time = " << d.count() << "ms\n";
	std::cout << "result = " << result/N <<std::endl << std::endl;



	return 0;
}