#include <iostream> 
#include <vector>
#include "particle.h"

void Particle::print_loc2(){
	for(int i: this->location){
		std::cout << i << " | ";
	}
	std::cout<< std::endl;
	return ;
};