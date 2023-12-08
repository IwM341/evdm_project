#include <iostream>
#include <string>

#include <time.h>
#include <chrono>
#include <cstdlib>
#include <array>
#include <memory>
#include "../src/utils/polynom.hpp"

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

float next_f(float x){
	float y = x + 3.14;
	return y - floor(y);
}

static float rand_init = 12.23;
float rnd(){
	rand_init = next_f(rand_init);
	return rand_init;
}

template <size_t num>
struct ptype{
	typedef evdm::PolynomHorner<float,num> type;	
};
template<>
struct ptype<8>{
	typedef evdm::PolynomEstrin8<float> type;	
};
template<>
struct ptype<12>{
	typedef evdm::PolynomEstrin12<float> type;	
};

template <size_t N>
constexpr size_t get_num(){
	if constexpr (N <= 6){
		return N;
	} else if (N <= 8){
		return 8;
	} else {
		return 12;
	}
}

#ifdef __GNUC__
#define NO_INLINE __attribute__ ((noinline)) 
#else
#define NO_INLINE __declspec(noinline)
#endif

constexpr static size_t Deg = 11;
using poly_t = typename ptype<get_num<Deg + 1>()>::type;

NO_INLINE float get_rand_number(poly_t * Poly,float x){
	return x;
}
NO_INLINE float eval_horner(poly_t * Poly,float x){
	return Poly->__horner<Deg+1>(x);
}
NO_INLINE float eval_default(poly_t * Poly,float x){
	return (*Poly)(x);
}

#define RAND_COEFFS
#undef RAND_COEFFS
int main(int argc,char ** argv){
	//int cnt = std::stoi(argv[1]);

	
	size_t Ntest = 100000000;

	

	alignas(8*alignof(float)) float coeffs[16];

	float sum = 0;
	srand(time(0));

	rand_init = rand()/(RAND_MAX+0.0f);

	using namespace std;
	cout << "Coeffs: ";
	for(size_t j = 0;j<16;++j){
		coeffs[j] = rnd();
		cout << (j ? ", " : "") << coeffs[j];
	} cout << endl;
	poly_t P(coeffs,Deg+1);
	cout << "Polynom: ";
	for(size_t i=0;i<P.size();++i){
		cout << (i ? ", " : "") << P[i];
	} cout << endl;
	cout << "P(0) = " << P(0) << " vs " << coeff_evaulator1<Deg>::eval(coeffs,0.0f) << endl;
	cout << "P(1) = " << P(1) << " vs " << coeff_evaulator1<Deg>::eval(coeffs,1.0f) << endl;
	cout << "P(2) = " << P(2) << " vs " << coeff_evaulator1<Deg>::eval(coeffs,2.0f) << endl;

	auto t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		#ifdef RAND_COEFFS
		for(size_t j = 0;j<8;++j){
			coeffs[j] = next_f(coeffs[j]);
		}
		#endif
		sum += get_rand_number(&P,rnd());
	}
    auto t1 = Time::now();
    fsec fs = t1 - t0;
    ms d = std::chrono::duration_cast<ms>(fs);
    std::cout << "default time = " << d.count() << "ms\n";
	std::cout << "sum = " << sum/Ntest <<std::endl << std::endl;
	auto def_time  = d.count();

	
	sum = 0;
	t0 = Time::now();
	for(size_t i=0;i<Ntest;++i){
		#ifdef RAND_COEFFS
		for(size_t j = 0;j<8;++j){
			coeffs[j] = next_f(coeffs[j]);
		}
		#endif
		sum += eval_horner(&P,rnd());
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "Horner Scheme = " << d.count() << "ms\n";
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
		sum += eval_default(&P,rnd());
	}
	t1 = Time::now();
	fs = t1 - t0;
	d = std::chrono::duration_cast<ms>(fs);
	std::cout << "Estring Scheme = " << d.count() << "ms\n";
	std::cout << "delta time = " << d.count() - def_time  << "ms\n";
	std::cout << "sum = " << sum/Ntest << std::endl << std::endl;


	return 0;
}