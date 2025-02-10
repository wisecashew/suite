#include "Simulation.hpp"
#include "lattice_directions.hpp"
#include "History.hpp"

std::set <int> Simulation::get_solvation_shell(){
	std::set   <int> solvation_shell_set; 
	std::array <std::array<int,3>,26> ne_list; 
	// int dop = static_cast<int>((*Polymers)[0].chain.size() ); 

	// get the first solvation shell 
	// auto start = std::chrono::high_resolution_clock::now(); 
	for (Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
			ne_list = obtain_ne_list ( p->coords, this->x, this->y, this->z); 
			for (std::array <int,3>& loc: ne_list){
				if (this->Lattice[lattice_index (loc, this->y, this->z)]->ptype[0] == 's'){
					solvation_shell_set.insert (lattice_index (loc, this->y, this->z)); 
				}
			}
		}
	}
	return solvation_shell_set;
}
