#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_

#include <iostream> // need this for couts and stuff 
#include <vector> // my primary data type 
#include "misc.h" // my collection of useful functions 


// #########################################################
// This is my Particle class. 
// it has the following attributes: 
// 1. loc - location. This is a 3D vector holding the location of the particle. 
//
// it has the following methods:
// 1. print_loc - just prints out location of the molecule. Useful for debugging. 

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

// END OF CLASS PARTICLE 
// #############################################################




// #############################################################
// This is my Polymer class 
// it has the following attributes:
// 1. chain - this is a std::vector of Particles. 
//
// it has the following methods:
// (nothing yet)...  


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

// END OF CLASS POLYMER
// ##############################################################




// ##############################################################
// This is my Grid class. Probably the most important class. 
// it has the following attributes: 
// 1. GRID - this is a vector of particles. Quite like Polymer, but Grid has more 
// more functionality 
// 2. loc_list - questionable if I want this. This holds the locations of 
// particles in Grid. 
// 3. x_len - length of box in x direction 
// 4. y_len - length of box in y direction 
// 5. z_len - length of box in z direction 
//
// it has the following methods: 
// 1. polymer_insertion - inserts a polymer into the grid, with seed = (0,0,0). 

class Grid{
public:
	// attributes
	std::vector <Particle> GRID;				// all the particles in the grid 
	std::vector <Polymer> pmers; 				// all the polymers in the grid 
	std::vector <std::vector <int>> loc_list;	// all the occupied position in the grid 

	int x_len; 									// length of grid along x-axis 
	int y_len; 									// length of grid along y-axis
	int z_len; 									// length of grid along z-axis 

	// constructor 
	Grid(int xl, int yl, int zl){
		std::vector <Particle> g; 
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

	std::vector <std::vector <int>> solvate(); 

	// methods 
};
// END OF CLASS POLYMER
// ##############################################################


# endif // _MC_CLASSES_H_