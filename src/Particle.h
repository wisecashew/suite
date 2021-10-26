#ifndef _PARTICLE_H_
#define _PARTICLE_H_ 

#include <vector>
#include <iostream> 
#include "misc.h"

class Particle{
public: 
	std::vector <int> loc; 
	// energetic interaction 
	// isotropy vs anisotropy

	// constructor 

	Particle(){

	};

	Particle(std::vector <int> &loc){
		this->loc = loc; 
		//std::cout << "Particle has been generated." << std::endl; 
	} 

	// destructor 
	~Particle(){
		//std::cout << "Particle has been erased from system memory." << std::endl;
	}

	void print_loc(); 

}; 



#endif // _PARTICLE_H_ 