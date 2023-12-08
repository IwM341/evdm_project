#include <iostream>
#include "../debugdef.hpp"
#include "../../src/utils/polynom.hpp"
#include "../../src/form_factors.hpp"
#include <variant>
int main(void){
	evdm::PolynomHorner<float,6> P = {1,2,3,1};
	print(11,'\n',"hrllo world!");
	PVAR(P(1));
	PVAR(P(2));
	PVAR(P(3));

	typedef evdm::QexpFactors<std::index_sequence<1,2,3,4>> ff_t;
	
	std::array<float,4> coeffs = {1,2,3,4};
	ff_t F(true,1,coeffs);
	PVAR(F.index());

	auto V = evdm::make_variant<std::variant<int,float>>(
		1,[](auto t){return typename decltype(t)::type (3.14);}
	);
	std::visit([](auto const & x){print(x);},V);

	auto true_factor = [b_2 = 1,coeffs](float q_2,float v_2,float){
		float y = q_2*b_2/4;
		return exp(-2*y)*(coeffs[0] + coeffs[1]*y+coeffs[2]*y*y+coeffs[3]*y*y*y)/y;
	};
	print("factor(1) = ",std::visit([](auto const & Factor){
		return Factor.ScatterFactor(2,1,0);
	},F.as_variant())," vs ", true_factor(2,1,0) );

	return 0;
}