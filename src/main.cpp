#include <iostream> 
#include <vector> 
#include <string> 

#include "Particle.h"
#include "Polymer.h"
#include "misc.h"
#include "grid.h"

int main(){

	// instantiate a grid; 
	std::vector <Particle> g; 
	const int x_len {4}, y_len{4}, z_len{4}; 

	Grid (g, x_len, y_len, z_len); 
	
	std::vector <int> seed{0,0,0};  // this is where the first particle will fall 
	std::vector<std::vector <int>> loc_list = {seed}; // this is the grid, basically. for now, anyway. 
	std::vector<std::vector <int>> loc_list2 = {seed}; // this is the grid, basically. for now, anyway. 

	int dop = 10; // this is the degree of polymerization 
	sarw(&loc_list, dop-1); // dop-1 because seed is already a monomer in place, so will need dop-1 more monomers  

	std::cout << "\nobtained loc_list, time for loc_list2..." << std::endl;

	sarw_pbc(&loc_list2, dop-1, 3,3,3);

	for (auto v:loc_list2){
		print(v);
	}

	std::cout << "\nprinting out loc_list2..."<<std::endl;


	std::vector <Particle> chain; 
	for (int i{0}; i< dop; i++){
		print(loc_list2.at(i));
		Particle p (loc_list2.at(i));
		// p.loc = loc_list2.at(i); 
		print(p.loc); 
		chain.push_back(p);
	}

	Polymer polymer {chain}; // created polymer chain 


	return 0;
}