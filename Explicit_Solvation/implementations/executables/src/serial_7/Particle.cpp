#include <iostream>
#include "Particle.h"

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                    PARTICLE METHODS. 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//============================================================
//============================================================
// 
// NAME OF FUNCTION: printCoords
//
// PARAMETERS: none, method on a particle
// 
// WHAT THE FUNCTION DOES: Given a particle, it will print out coordinates in a nice way 
//
// DEPENDENCIES: print
//
// THE CODE: 

void Particle::print_coords(){
    print(this->coords); 
    return; 
}


void print ( std::vector <Particle*>* LATTICE ){

	for ( Particle*& p: (*LATTICE) ){
		print( p->coords );
	}
}

void print(std::vector <Particle> pvec){
	for (Particle p: pvec){
		p.print_coords();
		// std::cout << p.ptype << std::endl;
	}
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of printCoords. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#