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

	std::vector <std::vector <int>> loc_list; 

	int x_len; // length of grid along x-axis 
	int y_len; // length of grid along y-axis
	int z_len; // length of grid along z-axis 

	// int total_cap {x_len*y_len*z_len}; 

	// static int current_number; 

	// constructor 
	Grid(std::vector <Particle> g, int xl, int yl, int zl){
		this->GRID = g; 
		this->x_len = xl; 
		this->y_len = yl; 
		this->z_len = zl; 
	};

	~Grid(){
	};

	// methods 

	// insert a particle in the grid  
	void particle_insertion();  

	// given a particle, will it be accepted by the GRID 
	bool particle_acceptance(Particle p); 

	// insert polymer in the grid 
	void polymer_insertion(int dop, std::vector <int> seed); 

	void get_loclist(); 

	void print_loclist();

	// methods 
};

#endif