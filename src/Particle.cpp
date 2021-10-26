#include <vector> 
#include "Particle.h"
#include "misc.h"
#include <iostream> 

void Particle::print_loc(){
	print(this->loc);
	std::cout << "print out the location of the particle" << std::endl;
	return ; 
}