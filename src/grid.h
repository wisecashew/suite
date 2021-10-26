#ifndef _GRID_H_
#define _GRID_H_ 

#include <vector>
#include "Particle.h"
#include "Polymer.h"

#include <iostream> 

class Grid{
public:
	// attributes
	std::vector <Particle> GRID; 
	int x_len; // length of grid along x-axis 
	int y_len; // length of grid along y-axis
	int z_len; // length of grid along z-axis 

	// constructor 
	Grid(std::vector <Particle> pvec, int x_len, int y_len, int z_len){

	};

	~Grid(){

	};


	// insert a particle in the grid  
	void particle_insertion(Particle p);  

	// given a particle, will it be accepted by the GRID 
	bool particle_acceptance(Particle p); 

	// insert polymer in the grid 
	void polymer_insertion(Polymer pmer); 

	

	// methods 
};

#endif