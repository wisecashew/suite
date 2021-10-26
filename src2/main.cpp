#include <iostream> 
#include <vector> 
#include <string> 

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

	// Gp.print_polymer(); 
	// Gp.print_solvent();
	// Gp.print_occupied();

	print(Gp.polymer.chain.at(0).loc);
	std::cout << Gp.polymer.chain.at(0).ptype << std::endl; 


	return 0;
}