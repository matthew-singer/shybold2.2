#include "network.h"

//void network::update(std::array<double, geneNN> val1, std::array<double, geneNN> val2, std::vector<double> inputs) {
void network::update(std::shared_ptr<genome> g, std::vector<double> inputs, std::vector<double> inputs_food,std::vector<double>outputfood_values) {
    for (int o = 0; o < output_values.size(); ++o) {
        output_values[o] = 0;
	//	std::cout<<"lo"<<'\n';
	/*	for (int k = 0;  k < outputfood_values.size(); ++k){
		outputfood_values[k] = 0;}*/
        for (int i = 0; i < inputs.size(); ++i) {
            //std::cout << "Input " << i << " " << o << " " << (o*(inputs.size()))+i << '\n';
            output_values[o] += (inputs[i] * g->getWeight( (o*inputs.size()) + i )  );
        }
            output_values[o] = -1.0 + 2.0 * (1.0/(1.0+pow(2.7183,-1*output_values[o])));
    }    
    // std::cout<<"lodown"<<'\n';
   
            
    for (int k = 0; k < outputfood_values.size(); ++k) {
      // outputfood_values[k] = 0;
      std::cout<<"is this read?"<<'\n';
      //std::cout<<"lodown234"<<'\n';
   
   //      outputfood_values[k] = 0;
      for (int i = 0; i < inputs_food.size(); ++i) {
	//std::cout << "Input " << i << " " << o << " " << (o*(inputs.size()))+i << '\n';
	outputfood_values[k] += (inputs_food[i] * g->getWeight( (k*inputs_food.size()) + i )  );
      }
      outputfood_values[k] = -1.0 + 2.0 * (1.0/(1.0+pow(2.7183,-1*output_values[k])));
    }
    
    

}
