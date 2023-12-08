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
	
	T sum = 0.0;
	T h = (b-a)/N;
    T c1 = 0.21132486540518711774;
    T c2 = 0.78867513459481288226;
	for (size_t i = 0; i < N; i++) {
		sum+= (f(a + (i+c1)*h) + f(a + (i+c2)*h))*h/2;
	}
	return sum;
}


#define X0 0.5
#define X1 0.76923465505284154552
#define X2 0.23076534494715845448
#define X3 ((T)0.95308992296933199641)
#define X4 ((T)0.04691007703066800359)

#define W0 ((T)0.28444444444444444444)
#define W1 ((T)0.23931433524968323402)
#define W2 ((T)0.23931433524968323402)
#define W3 ((T)0.11846344252809454376)
#define W4 ((T)0.11846344252809454376)


template <typename FuncType,typename T>	
auto integrateAB5(FuncType && f,T const & a,T const & b,size_t N){

	double sum = 0.0;
	double xi;
	double h = (b-a)/N;
	double h0 = h*X0,h1 = h*X1,h2 = h*X2, h3 = h*X3, h4 = h*h4;
	for (int i = 0; i < N; i++) {
		xi = a + i*h;
		sum += h*(W0*f(xi + h0)+
					W1*f(xi + h1)+
					W2*f(xi + h2)+
					W3*f(xi + h3)+
					W4*f(xi + h4));//sum5*h;
	}
	return sum;
}
}
#endif//INTEGRATION_HPP