#include <iostream>
#include <grob/multigrid.hpp>
#include <evdm/grid_variants.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <grob/serialization.hpp>
int main(void){
	using S = stools::PtreeSerializator<boost::property_tree::ptree>;

	auto ELG = evdm::GridEL<float>::grid_var_CUU(2,-5,10,[](auto e)->size_t{return (e+5)*2+2;});
	std::cout << ELG << std::endl;

	

	auto ELV = evdm::GridEL<float>::grid_var_CVV(2,-5,10,[](float e){return (6+e);},
											[](float e){return 5;},[](float e,float l){return 1+l;});

	std::cout <<ELV << std::endl;
	for(auto p : ELV){
		(void)p;
		//std::cout << p << std::endl;
	}

	evdm::GridEL<float> G1 = ELV;
	evdm::GridEL<float> G2 = ELG;
	
	auto p1 = stools::Serialize(G1, S{});
	auto p2 = stools::Serialize(G2, S{});

	boost::property_tree::write_json(std::cout,p1);
	boost::property_tree::write_json(std::cout,p2);

	return 0;
}