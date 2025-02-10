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
	std::array <int,3> coords;                                        // the coordinates of the particles
	short int orientation;
	std::string ptype; 

	bool operator<(const Particle& rhs)const{
		return (this->coords) < (rhs.coords); 
	}

	bool operator==(const Particle& rhs){
		return std::tie(coords, ptype, orientation) == std::tie(rhs.coords, rhs.ptype, rhs.orientation );
	} 


	// constructor 
	Particle(){};  // default constructor

	Particle (std::array <int, 3> crds, std::string type_, int orientation_): coords (crds), orientation (orientation_), ptype (type_) { };

	// destructor 
	~Particle(){}

	// print location of the particle 
	void print_coords(); 

};

void print (std::vector <Particle> pvec);
void print (std::map<std::array<int,3>,Particle*> LATTICE);

#endif
