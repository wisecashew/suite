#include <iostream> 
#include <vector> 
#include <string> 
#include <map>
#include "classes.h"
#include "misc.h"

int main(){

	int x_len {16}, y_len{16}, z_len{16}; 
	std::vector <int> seed{0,0,0}; 		 		// this is where the first polymer particle will fall 

	int dop = 12;								// degree of polymerization 

	Grid Gp( x_len, y_len, z_len); 				// instantiate Grid 
	Gp.polymer_insertion(dop, seed); 			// drop the polymer in 
	// Gp.solvate();							// solvate the system 
	// Gp.print_occupied();

	Gp.polymer.obtain_connectivity(); 			// get the map
	Gp.print_polymer();

	std::cout << "=======================" << std::endl;
	
	// std::cout << PolymerEnergySolvation(Gp.polymer.chain, x_len, y_len, z_len, -1) << std::endl;

	int energy = Gp.CalcEn();

	std::cout <<"Energy of the system is " << energy << std::endl; 

	// acceptance(1, 1); 

	Gp.ZeroIndexRotation_MC(); 
	Gp.print_polymer();

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;

	std::cout << "=======================" << std::endl;
	
	Gp.FinalIndexRotation_MC(); 
	Gp.print_polymer();

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;

	Gp.end_rotation_MC(); 
	Gp.print_polymer(); 

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "======================= YEET" << std::endl;


	Gp.FinalToZero_MC(); 
	Gp.print_polymer(); 

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;

	Gp.reptation_MC(); 
	Gp.print_polymer(); 

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;
	
	Gp.kink_jump_MC();
	Gp.print_polymer(); 

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;

	Gp.reptation_MC();
	Gp.print_polymer(); 

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;	

	Gp.reptation_MC();
	Gp.print_polymer(); 

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;	

	Gp.crank_shaft_MC();
	Gp.print_polymer();
	
	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;	

	/*
	Gp.reptation_MC(); 
	Gp.print_polymer();

	std::cout << "Energy of the system is " << Gp.CalcEn() << std::endl;
	std::cout << "=======================" << std::endl;
	*/
	

	return 0;
}