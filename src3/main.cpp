#include <iostream> 
#include <vector> 
#include <string> 
#include <map>
#include "classes.h"
#include "misc.h"

int main(){

	int x_len {5}, y_len{5}, z_len{5}; 
	std::vector <int> seed{0,0,0}; 		 		// this is where the first polymer particle will fall 

	int dop = 7;								// degree of polymerization 

	Grid Gp( x_len, y_len, z_len); 				// instantiate Grid 
	Gp.polymer_insertion(dop, seed); 			// drop the polymer in 
	// Gp.solvate();							// solvate the system 
	// Gp.print_occupied();

	Gp.polymer.obtain_connectivity(); 			// get the map

	Gp.print_polymer();
	Gp.FinalToZero();
	std::cout << "============= Performed reptation ============" << std::endl;
	Gp.print_polymer();

	Gp.FinalToZero();
	std::cout << "============= Performed reptation ============" << std::endl;
	Gp.print_polymer();

	Gp.ZeroToFinal();
	std::cout << "============= Performed reptation ============" << std::endl;
	Gp.print_polymer();
	

	std::cout << "Checking efficacy of conn..." << std::endl;
	std::cout << "Particle at index 5 is connected to " << std::endl;

	std::map <Particle, std::vector<Particle>> adj = Gp.polymer.conn;

	std::vector <Particle> pvec1 = adj[Gp.polymer.chain.at(5)]; 
	
	for (Particle p: pvec1){
		std::cout << "Am I here" << std::endl;
		p.print_loc();
	}
	
	Gp.polymer.find_cranks();
	Gp.print_polymer();
	std::cout << "===================" << std::endl;

	Gp.end_rotation(); 
	std::cout << "================== Post rotation 1 ================" << std::endl;
	Gp.print_polymer();

	std::cout << "===================" << std::endl;
	Gp.end_rotation();
	std::cout << "=================== Post rotation 2 ==============" << std::endl;
	Gp.print_polymer(); 
	std::cout << "===================" << std::endl;

	Gp.end_rotation(); 
	std::cout << "================== Post rotation 3 ================" << std::endl;
	Gp.print_polymer();

	std::cout << "===================" << std::endl;
	Gp.end_rotation();
	std::cout << "=================== Post rotation 4 ==============" << std::endl;
	Gp.print_polymer(); 
	std::cout << "===================" << std::endl;
	print(Gp.polymer.find_kinks());
	
	Gp.kink_jump();
	std::cout << "=================== Post kink jump 5 =================" << std::endl;
	Gp.print_polymer();

	Gp.end_rotation(); 
	Gp.print_polymer(); 

	
	print(Gp.polymer.find_kinks());
	Gp.kink_jump();
	std::cout << "=================== Post kink jump 6 =================" << std::endl;
	Gp.print_polymer();
	

	return 0;
}