#include <iostream> 
#include <vector> 
#include <string> 
#include <map>
#include "classes.h"
#include "misc.h"

int main(){

	int x_len {4}, y_len{4}, z_len{4}; 
	std::vector <int> seed{0,0,0};  // this is where the first polymer particle will fall 

	// std::vector <Particle> gp; 
	int dop = 15;
	Grid Gp( x_len, y_len, z_len); // instantiate Grid 
	Gp.polymer_insertion(dop, seed); // drop the polymer in 
	Gp.solvate();

	std::cout << "coordinates of monomers: " << std::endl;
	Gp.polymer.print_loc(); // print out coordinates
	std::cout << "=====================" << std::endl;

	Gp.polymer.obtain_connectivity(); // get the map

	std::cout << "Checking efficacy of conn..." << std::endl;
	std::cout << "Particle at index 5 is connected to " << std::endl;

	std::map <Particle, std::vector<Particle>> adj = Gp.polymer.conn;

	std::vector <Particle> pvec1 = adj[Gp.polymer.chain.at(5)]; 
	
	for (Particle p: pvec1){
		p.print_loc();
	}

	std::cout << "===================" << std::endl;
	

	return 0;
}