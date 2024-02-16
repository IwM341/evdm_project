#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP
#include <cmath>

namespace evdm{
template <typename FuncType,typename T>			
auto integrateAB1(FuncType && f,T const & a,T const & b,size_t N){
	
	T sum = 0.0;
	T h = (b-a)/N;
	for (size_t i = 0; i < N; i++) {
		sum+= f(a + (i+(T)0.5)*h)*h;
	}
	return sum;
}

template <typename FuncType,typename T>			
auto integrateAB2(FuncType && f,T const & a,T const & b,size_t N){
	
	T sum = (T)0.0;
	T h = (b-a)/N;
    T c1 = (T)0.21132486540518711774;
    T c2 = (T)0.78867513459481288226;
	for (size_t i = 0; i < N; i++) {
		sum+= (static_cast<T>(f(a + (i+c1)*h)) + 
				static_cast<T>(f(a + (i+c2)*h)))*h/2;
	}
	return sum;
}


#define __INTTEGR_X_COEFF0 0.5
#define __INTTEGR_X_COEFF1 0.76923465505284154552
#define __INTTEGR_X_COEFF2 0.23076534494715845448
#define __INTTEGR_X_COEFF3 ((T)0.95308992296933199641)
#define __INTTEGR_X_COEFF4 ((T)0.04691007703066800359)

#define __INTTEGR_W_COEFF0 ((T)0.28444444444444444444)
#define __INTTEGR_W_COEFF1 ((T)0.23931433524968323402)
#define __INTTEGR_W_COEFF2 ((T)0.23931433524968323402)
#define __INTTEGR_W_COEFF3 ((T)0.11846344252809454376)
#define __INTTEGR_W_COEFF4 ((T)0.11846344252809454376)


template <typename FuncType,typename T>	
auto integrateAB5(FuncType && f,T const & a,T const & b,size_t N){

	double sum = 0.0;
	double xi;
	double h = (b-a)/N;
	double h0 = h*__INTTEGR_X_COEFF0,
		   h1 = h*__INTTEGR_X_COEFF1,
		   h2 = h*__INTTEGR_X_COEFF2, 
		   h3 = h*__INTTEGR_X_COEFF3, 
		   h4 = h* __INTTEGR_X_COEFF4;
	for (int i = 0; i < N; i++) {
		xi = a + i*h;
		sum += h*(__INTTEGR_W_COEFF0*f(xi + h0)+
					__INTTEGR_W_COEFF1*f(xi + h1)+
					__INTTEGR_W_COEFF2*f(xi + h2)+
					__INTTEGR_W_COEFF3*f(xi + h3)+
					__INTTEGR_W_COEFF4*f(xi + h4));//sum5*h;
	}
	return sum;
}

#undef __INTTEGR_X_COEFF0
#undef __INTTEGR_X_COEFF1
#undef __INTTEGR_X_COEFF2
#undef __INTTEGR_X_COEFF3
#undef __INTTEGR_X_COEFF4 

#undef __INTTEGR_W_COEFF0
#undef __INTTEGR_W_COEFF1
#undef __INTTEGR_W_COEFF2 
#undef __INTTEGR_W_COEFF3 
#undef __INTTEGR_W_COEFF4 
}
#endif//INTEGRATION_HPP