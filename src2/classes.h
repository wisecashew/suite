#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_

#include <iostream> // need this for couts and stuff 
#include <vector> // my primary data type 
#include <string> 
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
	std::string ptype; 
	// energetic interaction 
	// isotropy vs anisotropy

	// constructor 

	Particle(){

	};

	Particle(std::vector <int> &loc){
		this->loc = loc; 
		//std::cout << "Particle has been generated." << std::endl; 
	} 

	Particle(std::vector <int> &loc, std::string s){
		this->loc=loc;
		this->ptype=s;
	}

	// destructor 
	~Particle(){
		//std::cout << "Particle has been erased from system memory." << std::endl;
	}

	void print_loc(); 

	std::vector <Particle> loc2part(std::vector <std::vector <int>> loc_list); 

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
	Polymer(){

	}; 

	Polymer(std::vector <Particle> chain_)
	: chain {chain_} {	
		for (Particle p: chain_){
			p.ptype = "polymer";
		}
	};

	// destructor 
	~Polymer(){
		//std::cout << "Polymer object has been destroyed from machine."<<std::endl;
	}

	void print_loc(); 

};

// END OF CLASS POLYMER
// ##############################################################


class Solvent{
public:
	std::vector <Particle> solvation; 

	// constructor

	Solvent(){

	}; 

	Solvent(std::vector <Particle> collection)
	: solvation { collection }{
		for (Particle p: collection){
			p.ptype = "solvent";
		}

	}; 

	//destructor 
	~Solvent(){

	};

	void print_loc(); 
};




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
	std::vector <Particle> occupied;				// all the particles in the grid 
	Polymer polymer; 				// all the polymers in the grid 
	Solvent solvent; 
	

	int x_len; 									// length of grid along x-axis 
	int y_len; 									// length of grid along y-axis
	int z_len; 									// length of grid along z-axis 

	// constructor 
	Grid(int xl, int yl, int zl){
		std::vector <Particle> g; 
		this->occupied = g; 
		this->x_len = xl; 
		this->y_len = yl; 
		this->z_len = zl; 
	
	};

	~Grid(){
	};

	// methods 

	// print out the coordinates of the polymer
	void print_polymer(); 

	// print out the coordinates of the solvent 
	void print_solvent(); 

	// update occupied as system is solvated or polymer is inserted
	void print_occupied();  


	// insert a polymer in the grid 
	void polymer_insertion(int dop, std::vector <int> seed); 

	// solvate the system 
	void solvate();


};
// END OF CLASS GRID
// ##############################################################


# endif // _MC_CLASSES_H_