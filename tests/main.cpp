#include <iostream> 
#include <vector> 
#include "particle.h"

int main() {

	std::vector <int> l {0,0,0}; 
	Particle p;
	p.location = l; 
	p.print_loc(); 
	p.print_loc2();

	return 0; 
}