#ifndef _PARTICLE_H_
#define _PARTICLE_H_
#include <map>

class Particle{
public:
	std::vector <int> location;
	void print_loc(){

		for (int i: location){
			std::cout << i << " | ";
		}
		std::cout<< std::endl;
	}

	Particle(){

	};


	void print_loc2();
};

class Polymer{
	std::map <Particle, Particle> conn; 
};


#endif