#ifndef _PARTICLE_SETS_H_
#define _PARTICLE_SETS_H_ 
#include "classes.h"


class ParticleSets{
public:
	std::vector <Polymer>   Polymers;
	std::vector <Particle*> Lattice;
	std::vector <Particle*> Cosolvent;

	// constructor
	ParticleSets(){}; 

	// destructor
	~ParticleSets(){};

}

#endif