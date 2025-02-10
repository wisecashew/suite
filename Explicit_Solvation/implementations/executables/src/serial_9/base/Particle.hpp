#ifndef _PARTICLE_H_
#define _PARTICLE_H_ 
#include "misc.hpp"


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// DEFINITIONS FOR CLASS PARTICLE 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


class Particle{
public: 
	int orientation;           // orientation of the particle
	std::string ptype;         // type of the particle
	std::array <int,3> coords; // the coordinates of the particle

	bool operator<(const Particle& rhs)const{
		return (this->coords) < (rhs.coords); 
	}

	bool operator==(const Particle& rhs){
		return std::tie(coords, ptype, orientation) == std::tie(rhs.coords, rhs.ptype, rhs.orientation);
	} 


	// constructor 1 (default)
	Particle(){};  

	// constructor 2 (assigns everything)
	Particle (int orientation_, std::string type_, std::array <int, 3> crds): orientation(orientation_), ptype (type_), coords (crds) { };

	// destructor 
	~Particle(){}

	// print location of the particle 
	void print_coords(); 

};

// void print (std::vector <Particle> pvec);
// void print (std::map<std::array<int,3>,Particle*> LATTICE);

#endif
