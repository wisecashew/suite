#include "Potts.hpp"

void Potts::simulate(){
	
	int start = this->step_number + 1;
	int stop  = this->step_number+this->max_iter+1;

	for (int i{start}; i < stop; ++i){
		this->step_number += 1;
		this->perturb_lattice_flip();
		if (this->acceptance){
			this->n_acceptances += 1;
		}

		// start dumping out information
		this->dump();

	}

	return;

}