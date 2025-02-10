#ifndef _PARTICLE_SETS_H_
#define _PARTICLE_SETS_H_ 
#include "Polymer.hpp"

class ParticleSets{
public:
	std::vector <Particle*> Lattice;
	std::vector <Polymer*>  Polymers;
	std::vector <Particle*> Solvent;
	std::vector <Particle*> Cosolvent;

	// constructor
	ParticleSets(){}; 

	// destructor
	~ParticleSets(){};

};

#endif
