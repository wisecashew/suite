#ifndef _POLYMER_H_
#define _POLYMER_H_ 
#include "Particle.h"

class History{
public:
	std::map <Particle*,std::vector<int>> SpinH;
	std::map <int,int> LocH;

	// constructor
	History(){};

	// destructor
	~History(){};

}
