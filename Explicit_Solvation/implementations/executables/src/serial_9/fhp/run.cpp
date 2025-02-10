#include "FHP.hpp"

void FHP::run_simple(){

	std::cout << "Running simple polymer simulation..." << std::endl;
	int start = this->step_number+1;
	int stop  = this->step_number+this->max_iter+1;
	for (int i {start}; i < stop; ++i){
		// std::cout << "@ step " << i << "." << std::endl;
		this->step_number += 1;
		this->perturb_system_simple();

		// start dumping out information
		this->dump();
		// std::cout << "attempts:    "; print(this->attempts);
		// std::cout << "acceptances: "; print(this->acceptances);
	}

	this->dump_lattice();
	this->dump_stats();

	// this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	// this->check_structures();

	return;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FHP::run_isotropic(){

	std::cout << "Running isotropic simulation (no orientation flip moves)..." << std::endl;
	int start = this->step_number+1;
	int stop  = this->step_number+this->max_iter+1;
	for (int i {start}; i < stop; ++i){
		// std::cout << "@ step " << i << "." << std::endl;
		this->step_number += 1;
		this->perturb_system_isotropic();

		// start dumping out information
		this->dump();
	}

	this->dump_lattice();
	this->dump_stats();

	// this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	// this->check_structures();

	return;
}


void FHP::run_dry(){

	std::cout << "Running dry (no solvent moves)..." << std::endl;
	int start = this->step_number+1;
	int stop  = this->step_number+this->max_iter+1;
	for (int i {start}; i < stop; ++i){
		// std::cout << "@ step " << i << "." << std::endl;
		this->step_number += 1;
		this->perturb_system_dry();

		// start dumping out information
		this->dump();
	}

	this->dump_lattice();
	this->dump_stats();

	// this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	// this->check_structures();

	return;
}

void FHP::run_straight(){

	std::cout << "Running straight..." << std::endl;
	int start = this->step_number+1;
	int stop  = this->step_number+this->max_iter+1;
	for (int i {start}; i < stop; ++i){
		// std::cout << "@ step " << i << "." << std::endl;
		this->step_number += 1;
		this->perturb_system_straight();

		// start dumping out information
		this->dump();
	}

	this->dump_lattice();
	this->dump_stats();

	// this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	// this->check_structures();

	return;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void FHP::run_debug(){
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

		std::cout << "checking structures" << std::endl;
		this->check_structures();
		std::cout << "done." << std::endl;
		std::cout << "dumping stuff." << std::endl;
		this->dump();
		std::cout << "done dumping. " << std::endl;
	}
	return;
}

void FHP::run(){
	std::cout << "Running the setup..." << std::endl;
	this->setup();
	// ===
	this->dump_opening_tiles(); 
	// === 
	std::cout << "Running the simulation..." << std::endl;
	(this->*run_ptr)();
	return;
}
