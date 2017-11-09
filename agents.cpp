#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "agents.h"

#define PI 3.14159

void agent::updatePrey() {
    std::vector<double> inputs;
    std::vector<double> inputs_food;
    std::vector<double> outputfood_values;
    // std::vector<double> inputs_food; 
    if (input_agent.size() != 0 && input_agent[0]) { //change this (loop througth input_agents and push_back onto inputs
       saw_last = true;
       //	inputs_food ={  this->distance(input_agent[0]),  atan2(input_agent[0]->y - y , input_agent[0]->x - x), 0.0, 1.0 };
        inputs = { this->distance(input_agent[0]),  atan2(input_agent[0]->y - y, input_agent[0]->x - x), 0.0, 1.0 };
    } else {
        //make appropariate number of inputs
        saw_last = false;
	//inputs_food = { 0, 0, getRandom(), 1.0 };

        inputs = { 0, 0, getRandom(), 1.0 };
    }
    //you hand it chrome 1 values and chrome 2 values from genome
    n->update(g, inputs, inputs_food, outputfood_values);
    //move based on ouptputs
    //    move_food(v->outputfood_values[0], v - > output_foodvalues[1], v->outputfood_values[2], v->outputfood_values[3])
    move_x_y(n->output_values[0] * pred_capture,  n->output_values[1] * pred_capture);
    //move_mag_theta(n->output_values[0], n->output_values[1], n->output_values[2]);
    //move_mag_theta(n->output_values[0], n->output_values[1], 0.0, pred_capture);
    


}

void agent::updatePred(int time) {
    std::vector<double> inputs;
    std::vector<double> inputs_food;
    std::vector<double> inputs_next;
    std::vector<double> outputfood_values;
    double tmp_x = x;
    double tmp_y = y;
    double tmp_x2 = x2;
    double tmp_y2 = x2;

    bubbleSort(inputs_food);
    
    if (input_agent.size() != 0&& input_agent[0]) { //change this (loop througth input_agents and push_back onto inputs
    // if (input_agent.size() !=0){

       // std::cout<<"ho"<<'\n';
        //inputs = { this->distance(input_agent[0]),  atan2(input_agent[0]->y - y, input_agent[0]->x - x), 0.0, 1.0 };
	    inputs_food = { 10.0*(input_agent[0]->x-x)/sensing_range_pred,  10.0*(input_agent[0]->y-y)/sensing_range_pred, 10.0*(input_agent[0]->x2-x2)/sensing_range_pred,  10.0*(input_agent[0]->y2-y2), 0.0, 1.0 };
	    ////outputfood_values= {inputs_food[0]*+};
	    //   std::cout<<"howdy" <<'\n';
	    // n->update(g, inputs, inputs_food, outputfood_values);
	    
      // Original implementation -- keeping in case my fix doesn't work --
        // Inputs:  0 = x1; 1 = y 1; 2 = x2; 3 = y2
        if (inputs_food[0]* g->getWeight(16)+ inputs_food[1]*g->getWeight(17)> inputs_food[2]*g->getWeight(16)+inputs_food[3]*g->getWeight(17)){
	      inputs = {inputs_food[0], inputs_food[1],getRandom(),1.0};

	    }

        if (inputs_food[0]* g->getWeight(16)+ inputs_food[1]*g->getWeight(17)< inputs_food[2]*g->getWeight(16)+inputs_food[3]*g->getWeight(17)){
            inputs = {inputs_food[2], inputs_food[3],getRandom(),1.0};

            }





	    //	     inputs ={  outputfood_values[0]+outputfood_values[1], outputfood_values[2] + outputfood_values[3],  0.0, 1.0 };
	    // inputsnext = {};
	    //	    inputs ={  inputs_food[1]* g->getWeight(1)+ inputs_food[3]* g->getWeight(3)+ inputs_food[5]* g->getWeight(5)+ inputs_food[6]* g->getWeight(6),inputs_food[2]* g->getWeight(2)+ inputs_food[4]* g->getWeight(4)+ inputs_food[5]* g->getWeight(5)+ inputs_food[6]* g->getWeight(6), 0,1}; 



	    // std::cout<<"cowboy" <<'\n';
	    //  inputs_food = {10.0*(input_agent[0]->x-x)/sensing_range_pred,  10.0 			   *(input_agent[0]->y-y)/sensing_range_pred, 0.0, 1.0}; 
	    //std::cout<<"befreodown"<<'\n';
	    
	    //	inputs = { 10.0*(input_agent[0]->x-x)/sensing_range_pred,  10.0*(input_agent[0]->y-y)/sensing_range_pred, 0.0, 1.0 };
    } else {
      inputs = { 0, 0, getRandom(), 1.0 }; // With random movement
    //   inputs = { 0, 0, 0, 1.0 }; // Without random movement
      inputs_food = { 0, 0,0,0, getRandom(), 1.0 }; 
    }
     //	std::cout<<"lodown"<<'\n';

    //you hand it chrome 1 values and chrome 2 values from genome
	n->update(g, inputs, inputs_food, outputfood_values);
   // move_mag_theta(n->output_values[0], n->output_values[1], 0.0, pred_capture);
    move_x_y(n->output_values[0] * pred_capture,  n->output_values[1] * pred_capture);
    
    double x2dif = x2- tmp_x2;
    double y2dif = y2-tmp_y2;
    double xdif = x - tmp_x;
    double ydif = y - tmp_y;
    double d = sqrt((xdif * xdif) + (ydif * ydif));
    double d2 = sqrt ((x2dif * x2dif) +(y2dif * y2dif));
    sum_x_step += xdif;
    sum_y_step += ydif;
    sum_d += d;
    sq_sum_x_step += xdif * xdif;
    sq_sum_y_step += ydif * ydif;
    sq_sum_d += d * d;
    
    calc_stuff();
    consume(time);
    chance_of_death(time);

}

//Need to change to get sorted vector of agents inputs

void bubbleSort(std::vector<double>& inputs_food)
		{
		  std::cout<<"hi"<<'\n';
		  		  bool swapp = true;
		  while(swapp){
		    swapp = false;
		    for (size_t i = 0; i < inputs_food.size()-1; i++) {
		       
                     if (inputs_food[i]>inputs_food[i+1] ){
		       inputs_food[i]>inputs_food[i+1] ;{
			inputs_food[i] +=   inputs_food[i+1];
			inputs_food[i+1] = inputs_food[i] - inputs_food[i+1];
			inputs_food[i] -=inputs_food[i+1];
			swapp = true;
		      }
		    }
		  }
		}
		  
		}

void agent::getNearestAgentPrey(const std::shared_ptr<agent> &a) {

   if (a->alive && valid_agent(a)) {
       input_agent[0] = a;
    } 
    /*
    if (input_agent.size() == 0) {
        if (this->distance(a) < sensing_range_prey ) {
            input_agent.push_back(a);
        }
    } else if (!input_agent[0] || this->distance(a) < this->distance(input_agent[0])) {
        //change it here, sort the input_agents and push them appropriately on
        if (this->distance(a) < sensing_range_prey) {
            input_agent[0] = a;
        }
    }*/
}

void agent::getNearestAgentPred(const std::shared_ptr<agent> &a) {
   if (valid_agent(a)) {
        input_agent[0] = a;
    }  
}
   void agent::getSecondAgentPred(const std::shared_ptr<agent> &b){
     if(valid_agent(b)) {

	 input_agent[1] = b;
     }
   }
   /*if (input_agent.size() == 0) {
        if (this->distance(a) < sensing_range_pred) {
            input_agent.push_back(a);
        }
    } else if ( !input_agent[0] || this->distance(a) < this->distance(input_agent[0]) || !input_agent[0]->alive ) {
        //same thing here
        if (this->distance(a) < sensing_range_pred) {
            input_agent[0] = a;
        }
	}*/


//figure out what you want to do here
/*void agent::getSecondAgentPred(const std::shared_ptr<agent> &b) {
  if (valid_agent(b)) {
    input_agent[1] = b;
  }
  }*/
 void agent::consume(int time) {
    if (input_agent[0]) {
        if (this->distance(input_agent[0]) < pred_capture) {
            input_agent[0]->alive = false;
            input_agent[0]->lastTime = time;
            ++fitness;
        }
        /*if (this->distance(input_agent[0]) < pred_capture && input_agent[0]->alive) {
            input_agent[0]->alive = false;
            input_agent[0]->lastTime = time;
            ++fitness;
        }*/
    }
}

void agent::move_x_y(double dx, double dy) {
    x = std::max(0.0, std::min(x + dx, sizeX));
    y = std::max(0.0, std::min(y + dy, sizeY));
}

void agent::move_mag_theta(double mag, double theta, double direction_facing, double move) {
    /*std::cout << "X " << sin(theta) * mag << '\n';
    std::cout << "Y " << cos(theta) * mag << '\n';
    std::cout << "Theta " << theta << '\n';
    std::cout << "Mag " << mag << '\n';*/
    x = std::max(0.0, std::min(cos(theta) * mag * move + x, sizeX));
    y = std::max(0.0, std::min(sin(theta) * mag * move + y, sizeY));
    //angle_facing = std::fmod(angle_facing + direction_facing, 2 * PI);
    
}

void agent::output_data(std::fstream &o, bool prey, int gen, int index) {
    o << gen << '\t';
    o << idTrue << '\t';
    if (prey) {
        o << lastTime << '\t';
    } else {
        o << fitness << '\t';
    }
    
    for (int i = 0; i < geneNN; ++i) {
        o << g->getWeight(i) << '\t';
    }
    
    o << ((g->c1)->metabolic + (g->c2)->metabolic)/2.0 << '\t';
    o << ((g->c1)->radius + (g->c2)->radius) /2.0 << '\t';
    o << g->getStddev() << '\t'; 
    if (prey) {
        o << calcFitnessPrey() << '\t';
    } else {
        o << calcFitnessPred() << '\t';
    }
    o << 
        sum_x / deathtime << "\t" << sum_y / deathtime << "\t" << 
        sq_sum_x / deathtime - (sum_x /deathtime) * (sum_x/deathtime) << '\t' << 
        sq_sum_y / deathtime - (sum_y /deathtime) * (sum_y/deathtime) << '\t' << 
        sum_x_step / deathtime << "\t" << sum_y_step / deathtime << "\t" << 
        sq_sum_x_step / deathtime - (sum_x_step /deathtime) * (sum_x_step/deathtime) << '\t' << 
        sq_sum_y_step / deathtime - (sum_y_step /deathtime) * (sum_y_step/deathtime) << '\t' <<
        sum_d / deathtime << '\t' <<
        sq_sum_d / deathtime - (sum_d / deathtime) * (sum_d / deathtime) << '\t' << deathtime << '\t' << prey_pop_count << '\t' << i_s << std::endl;
}
