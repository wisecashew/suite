#ifndef _HISTORY_H_
#define _HISTORY_H_ 
#include "Particle.hpp"

class History {
public:

	std::vector <std::array<int,3>> old_cut;
	std::vector <std::array<int,3>> new_cut;
	std::map <Particle*,std::vector<int>> SpinH;
	std::map <Particle*,std::vector<int>> LocH;
	std::map <Particle*,std::vector<Particle*>>  SolventIdentities;

	// constructor
	History(){};

	// destructor
	~History(){};

};

#endif
