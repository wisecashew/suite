#ifndef _POLYMER_H_
#define _POLYMER_H_ 

#include <vector> 
#include <iostream>
#include "Particle.h"

class Polymer{
public:
	std::vector <Particle> chain;

	// constructor 
	Polymer(std::vector <Particle> chain_)
	: chain {chain_} {	
	};

	// destructor 
	~Polymer(){
		//std::cout << "Polymer object has been destroyed from machine."<<std::endl;
	}


};


#endif
