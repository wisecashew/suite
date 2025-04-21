#include "Simulation.hpp"

void Simulation::run_isotropic(){

	std::cout << "Running (isotropic)..." << std::endl;
	int start = this->step_number+1;
	int stop  = this->step_number+this->max_iter+1;
	for (int i {start}; i < stop; ++i){
		// std::cout << "@ step " << i << "." << std::endl;
		this->step_number += 1;
		this->perturb_system_isotropic();
		if (this->IMP_BOOL){
			this->acceptances[this->move_number] += 1;
		}
		// start dumping out information
		this->dump_local();
		this->dump_lattice();
	}

	return;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Simulation::run_straight(){
	std::cout << "Running..." << std::endl;
	int start = this->step_number+1;
	int stop  = this->step_number+this->max_iter+1;
	for (int i {start}; i < stop; ++i){
		// std::cout << "@ step " << i << "." << std::endl;
		this->step_number += 1;
		this->perturb_system_straight();
		if (this->IMP_BOOL){
			this->acceptances[this->move_number] += 1;
		}
		// start dumping out information
		this->dump_local();
		this->dump_lattice();
	}

	// this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	// this->check_structures();

	return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Simulation::run_debug(){
	std::cout << "Inside the debugging run!" << std::endl;
	int start = this->step_number+1;
	int stop  = this->step_number+this->max_iter+1;
	for (int i {start}; i < stop; ++i){
		this->step_number += 1;
		std::cout << "Initial config for step is:" << std::endl;
		for (int j{0}; j< int((Polymers)[0].chain.size()); ++j){
			print((Polymers)[0].chain[j]->coords, ", "); std::cout << "o = " << (Polymers)[0].chain[j]->orientation << std::endl;
		}
		std::cout << "Contacts = "; print (contacts);
		std::cout << "-------------------------- " << std::endl;
		std::cout << "Step number: " << i << "." << std::endl; 
		std::cout << "Executing..." << std::endl << std::endl;
		
		this->perturb_system_debug();
		if (this->IMP_BOOL){
			this->acceptances[this->move_number] += 1;
		}

		this->check_structures();
		this->dump_local();
		this->dump_lattice();
	}
	return;
}
