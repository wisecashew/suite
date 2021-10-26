#ifndef _PARTICLE_H_
#define _PARTICLE_H_

class Particle{
public:
	std::vector <int> location;
	void print_loc(){

		for (int i: location){
			std::cout << i << " | ";
		}
		std::cout<< std::endl;
	}

	void print_loc2();
};


#endif