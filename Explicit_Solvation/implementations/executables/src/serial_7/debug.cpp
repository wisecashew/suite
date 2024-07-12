#include "Simulation.h"
#include "lattice_directions.h"
#include "misc.h"

// all methods will have the prefix `debug_`
// all properties and variable unique to these methods will have the prefix `db_`

void Simulation::debug_calculate_energy(double* db_energy, std::array<double,8>* db_contacts){

	*db_energy   = 0.0;
	*db_contacts = {0,0,0,0,0,0,0,0}; 
	std::array <std::array <int,3>, 26> ne_list; 

	// run energy computations for every monomer bead 
	// std::cout << "Running a debugging energy calculation..." << std::endl;
	for (Polymer& pmer: this->Polymers) {
		for (Particle*& p: pmer.chain) {
			ne_list = obtain_ne_list(p->coords, this->x, this->y, this->z);
			for (std::array <int,3>& loc: ne_list) {
				this->debug_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], db_contacts, db_energy);
				// std::cout << "\t\t@neighboring particle, coords = "; print(loc, ", "); std::cout << "ptype = " << this->Lattice[lattice_index(loc, this->y, this->z)]->ptype << ", energy = " << *db_energy << std::endl;
			}
		}
	}

	for (Particle*& p: this->Cosolvent){
		// std::cout << "Inside cosolvent loop vector..." << std::endl;
		ne_list = obtain_ne_list (p->coords, this->x, this->y, this->z);
		for (std::array <int,3>& loc: ne_list){
			if (this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "m1" || 
			this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "s2"){
				continue;
			}
			this->debug_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], db_contacts, db_energy);
		}
	}
	return;
}

void Simulation::debug_pair_interaction(Particle* p1, Particle* p2, std::array<double,8>* db_contacts, double* db_energy){

	double dot_product = 0;
	double theta_1     = 0;
	double theta_2     = 0;
	double magnitude   = 0;

	std::array <int,3> connvec = {0, 0, 0};

	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 

	std::string interaction_type = std::get<0>((this->InteractionMap)[particle_pair]);
	switch (interaction_type[0]){

	case 'i':
		if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
			(*db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])+1] += 0.5;
			*db_energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5;
		}
		else {
			(*db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])+1] += 1;
			*db_energy += std::get<1>((this->InteractionMap)[particle_pair]);
		}
		break;

	case 'p':
		dot_product = take_dot_product ( p1->orientation, p2->orientation );
		if (dot_product > 0.54){
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				(*db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])] += 0.5; 
				*db_energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 

			}
			else {
				(*db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])] += 1;
				*db_energy += std::get<1>((this->InteractionMap)[particle_pair]); 
			}	
		}
		else {
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				(*db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 0.5; 
				*db_energy += std::get<2>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(*db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 1; 
				*db_energy += std::get<2>((this->InteractionMap)[particle_pair]); 
			}
		}
		break;

	case 'a':
		connvec = subtract_arrays ( &(p2->coords), &(p1->coords) );
		modified_direction (&connvec, x, y, z); 
		magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
		theta_1   = std::acos (take_dot_product (scale_arrays (1/magnitude , &connvec), Or2Dir[p1->orientation] ) );
		theta_2   = std::acos (take_dot_product (scale_arrays (-1/magnitude, &connvec), Or2Dir[p2->orientation] ) );

		if ( theta_1 + theta_2 > M_PI/2 ){
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ) {
				// std::cout << "misaligned antiparallel m1-m1..." << std::endl;
				(*db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 0.5; 
				*db_energy += std::get<2>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(*db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 1; 	
				*db_energy += std::get<2>((this->InteractionMap)[particle_pair]);
			}
		}
		else {
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				// std::cout << "aligned antiparallel m1-m1..." << std::endl;
				(*db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])] += 0.5; 
				*db_energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(*db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])] += 1; 
				*db_energy += std::get<1>((this->InteractionMap)[particle_pair]); 
			}
		}
		break;
	}

	return;
}

void Simulation::debug_isolated_pair_particle_interaction(Particle* p1, Particle* p2, std::array<double,8>* db_contacts, double* db_energy){

	std::array<int,3>     connvec = subtract_arrays(&(p2->coords), &(p1->coords));
	modified_direction  (&connvec, this->x, this->y, this->z);

	if (std::find(adrns.begin(), adrns.end(), connvec) == adrns.end()) {
		*db_energy = 0;
		return;
	}

	int    c_idx     = -1;
	double pE        = 0;
	double dot_prod  = 0;
	double magnitude = 0;
	double theta_1   = 0;
	double theta_2   = 0;
	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype);

	switch (std::get<0>( (this->InteractionMap)[particle_pair])[0] ) {
		case 'i':
			pE     = std::get<1>((this->InteractionMap)[particle_pair]); 
			c_idx  = std::get<3>((this->InteractionMap)[particle_pair]);
		break;

		case 'p':
			dot_prod  = take_dot_product (p1->orientation, p2->orientation); 
			if (dot_prod > 0.54) {
				pE    = std::get<1>((this->InteractionMap)[particle_pair]);
				c_idx = std::get<3>((this->InteractionMap)[particle_pair]);
			}
			else {
				pE    = std::get<2>((this->InteractionMap)[particle_pair]);
				c_idx = std::get<4>((this->InteractionMap)[particle_pair]);	
			}
		break;

		case 'a':
			magnitude = std::sqrt ( connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2] );
			theta_1   = std::acos (take_dot_product (scale_arrays (1/magnitude, &connvec),  Or2Dir[p1->orientation]) ); 
			theta_2   = std::acos (take_dot_product (scale_arrays (-1/magnitude, &connvec), Or2Dir[p2->orientation]) );

			if (theta_1 + theta_2 > M_PI/2) {
				pE = std::get<2>((this->InteractionMap)[particle_pair]);
				c_idx = std::get<4>((this->InteractionMap)[particle_pair]);  
			}
			else {
				pE = std::get<1>((this->InteractionMap)[particle_pair]);
				c_idx = std::get<3>((this->InteractionMap)[particle_pair]); 
			}
		break; 
	}

	(*db_contacts)[c_idx] += 1;
	*db_energy             = pE;

	return;

}

void Simulation::debug_checks_energy_contacts(double E_final, std::array<double,8> final_contacts){

	double E_debug = 0;
	std::array<double,8> debug_contacts = {0, 0, 0, 0, 0, 0, 0, 0}; 

	// run the debugger
	this->debug_calculate_energy(&E_debug, &debug_contacts); 
	if (E_debug != E_final){
		std::cout << "E_real = " << E_debug << ", E_final = " << E_final << "." << std::endl;
		exit(EXIT_FAILURE);
	} 
	else {
		std::cout << "energies match @ " << E_debug << "!" << std::endl;
	}
	if (debug_contacts != final_contacts){
		std::cout << "real_contacts = "; print(debug_contacts); 
		std::cout << "final_contacts = "; print(final_contacts);
		exit(EXIT_FAILURE);
	}
	else {
		std::cout << "contacts match @ "; print(final_contacts, ""); std::cout << "!" << std::endl;
	}

	return; 
}

/////////////////////////////////////////////////////////////////

void Simulation::debug_monomer_orientation_sampler_forwards(EnhancedOrientationFlipper* eof, std::array<double,8>* contacts_sys, double E_sys, int iteration_idx, int polymer_idx, int monomer_idx){

	int monomer_lattice_idx = lattice_index(this->Polymers[polymer_idx].chain[monomer_idx]->coords, this->y, this->z);
	eof->initial_orientations[iteration_idx] = (this->Polymers)[polymer_idx].chain[monomer_idx]->orientation;
	this->neighbor_energetics(monomer_lattice_idx, &(eof->initial_contacts), &(eof->initial_E));

	for (int j{0}; j < eof->ntest; ++j){
		this->Polymers[polymer_idx].chain[ monomer_idx ]->orientation = rng_uniform (0, 25); 
		eof->orientations[j]   = this->Polymers[polymer_idx].chain[monomer_idx]->orientation;
		this->neighbor_energetics(monomer_lattice_idx, &(eof->perturbed_contacts), &(eof->perturbed_E));
		eof->energies[j]        = E_sys - eof->initial_E + eof->perturbed_E;
		eof->contacts_store[j]  = subtract_arrays(*contacts_sys, eof->initial_contacts);
		eof->contacts_store[j]  = add_arrays     (eof->contacts_store[j], eof->perturbed_contacts);
		eof->perturbed_contacts = {0, 0, 0, 0, 0, 0, 0, 0};
		eof->perturbed_E        = 0;
		this->debug_checks_energy_contacts(eof->energies[j], eof->contacts_store[j]);
	}

	return;

}

void Simulation::debug_choose_monomer_state_forward(EnhancedOrientationFlipper* eof, int iterator_idx, int polymer_idx, int monomer_idx){

	// get the total boltzmann sum
	for ( int k{0}; k < eof->ntest; ++k ){
		eof->boltzmann[k] = std::exp (-1/this->T*(eof->energies[k] - eof->Emin));
		eof->rboltzmann  += eof->boltzmann [k]; 
	}

	// get the rng to sample the right orientation
	eof->sampler_rng   = rng_uniform (0.0, 1.0);

	// get the appropriate state
	for ( int j{0}; j<5; ++j){
		eof->sampler_rsum += eof->boltzmann[j]/eof->rboltzmann;
		if (eof->sampler_rng < eof->sampler_rsum){
			eof->sampler_idx = j; 
			break; 
		}
	}

	// make the jump to the new state 
	std::cout << "orientation to pick = "; print(eof->orientations);
	eof->final_orientations[iterator_idx]                       = eof->orientations[eof->sampler_idx];
	this->Polymers[polymer_idx].chain[monomer_idx]->orientation = eof->orientations[eof->sampler_idx];
	eof->prob_o_to_n *= eof->boltzmann[eof->sampler_idx]/eof->rboltzmann;

	return;
}

void Simulation::debug_monomer_orientation_sampler_backwards_0(EnhancedOrientationFlipper* eof, std::array<double,8>* contacts_sys, double E_sys, int iteration_idx, int polymer_idx, int monomer_idx){

	int m_lattice_idx = lattice_index((this->Polymers)[polymer_idx].chain[monomer_idx]->coords, this->y, this->z);
	this->neighbor_energetics(m_lattice_idx, &(eof->initial_contacts), &(eof->initial_E));

	(this->Polymers)[polymer_idx].chain[monomer_idx]->orientation = eof->initial_orientations[iteration_idx];
	this->neighbor_energetics(m_lattice_idx, &(eof->perturbed_contacts), &(eof->perturbed_E));
	eof->energies[0] = E_sys - eof->initial_E + eof->perturbed_E;
	eof->contacts_store[0]  = subtract_arrays(*contacts_sys, eof->initial_contacts);
	eof->contacts_store[0]  = add_arrays     (eof->contacts_store[0], eof->perturbed_contacts);
	std::cout << "Inside fourth set of flip test..." << std::endl;
	this->debug_checks_energy_contacts(eof->energies[0], eof->contacts_store[0]);

	eof->perturbed_contacts = {0,0,0,0,0,0,0,0};
	eof->perturbed_E        = 0;

	return;

}

void Simulation::debug_monomer_orientation_sampler_backwards(EnhancedOrientationFlipper* eof, std::array<double,8>* contacts_sys, double E_sys, int iteration_idx, int polymer_idx, int monomer_idx){

	this->debug_orientation_sampler_backwards_0(eof, contacts_sys, E_sys, iteration_idx, polymer_idx, monomer_idx);

	int m_lattice_idx = lattice_index((this->Polymers)[polymer_idx].chain[monomer_idx]->coords, this->y, this->z);
	std::cout << "Inside fifth set of flip test." << std::endl;
	for (int j{1}; j<eof->ntest; ++j){
		this->Polymers[polymer_idx].chain[monomer_idx]->orientation = rng_uniform (0, 25); 
		this->neighbor_energetics(m_lattice_idx, &(eof->perturbed_contacts), &(eof->perturbed_E));
		eof->energies[j] = E_sys - eof->initial_E + eof->perturbed_E; //  Ei + Epert;
		eof->contacts_store [j]  = subtract_arrays (*contacts_sys, eof->initial_contacts);
		eof->contacts_store [j]  = add_arrays      (eof->contacts_store[j], eof->perturbed_contacts);
		eof->perturbed_contacts = {0, 0, 0, 0, 0, 0, 0, 0};
		eof->perturbed_E        = 0;
		this->debug_checks_energy_contacts(eof->energies[j], eof->contacts_store[j]);
	}

	return;

}

//////////////////////////////////////////////////////////////////

void Simulation:: debug_orientation_sampler_forwards(EnhancedOrientationFlipper* eof, std::array<double,8>* contacts_sys, double E_sys, int iterator_idx, int lat_idx){

	eof->initial_orientations[iterator_idx] = this->Lattice[lat_idx]->orientation;
	this->neighbor_energetics(lat_idx, eof->initial_contacts, eof->initial_E);

	for (int j{0}; j<eof->ntest; ++j){
		(this->Lattice)[lat_idx]->orientation = rng_uniform(0, 25);
		eof->orientations[j] = this->Lattice[lat_idx]->orientation;
		this->neighbor_energetics(lat_idx, eof->perturbed_contacts, this->perturbed_E);
		eof->energies[j] = E_sys - eof->initial_E + eof->perturbed_E;
		eof->contacts_store[j]  = subtract_arrays(contacts_sys, &(eof->initial_contacts));
		eof->contacts_store[j]  = add_arrays     (&(eof->contacts_store[j]), &(eof->perturbed_contacts));
		eof->perturbed_contacts = {0, 0, 0, 0, 0, 0, 0, 0};
		eof->perturbed_E        = 0;
		this->debug_checks_energy_contacts(eof->energies[j], eof->contacts_store[j]);
	}

	return;

}



//////////////////////////////////////////////////////////////////

void Simulation::debug_tail_rotation(int p_idx){

	// get the neighborlist of the particle at index 1
	std::array<int,3> loc_0 = this->Polymers[p_idx].chain[0]->coords;
	std::array<int,3> loc_1 = this->Polymers[p_idx].chain[1]->coords;

	// define some container for energy
	double                final_E        = 0; // energy of the final configuration of the system
	std::array <double,8> final_contacts = {0, 0, 0, 0, 0, 0, 0, 0};

	// generate a neighbor list 
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

	// get the choice for a neighbor
	int choice = rng_uniform(0, 25);

	// instantiate some containers -- 
	// this object implicitly saves us a lot of lines of code of just instantiating variables
	Container container;

	if ( this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->ptype[0] == 's' ){

		this->particle_swap_with_monomer_update(&container, loc_0, p_idx, 0, lattice_index(ne_list[choice], this->y, this->z));
		this->check_structures();
		this->Polymers[p_idx].print_chain();

		final_E = this->sysEnergy - (container.E1_initial + container.E2_initial - container.Epair_initial); 
		final_E = final_E + (container.E1_final   + container.E2_final  - container.Epair_final);
		final_contacts = add_arrays(subtract_arrays(this->contacts, add_arrays (container.c1_initial, container.c2_initial)), add_arrays (container.c1_final, container.c2_final));
		final_contacts = add_arrays(final_contacts, container.cpair_initial);
		final_contacts = subtract_arrays(final_contacts, container.cpair_final);

		// run the debugger
		std::cout << "Running the debugger..." << std::endl; 
		this->debug_checks_energy_contacts(final_E, final_contacts);
		std::cout << "Ran the debugger!" << std::endl;

		if (metropolis_acceptance(this->sysEnergy, final_E, this->T)){
			this->sysEnergy = final_E;
			this->contacts  = final_contacts;
			this->IMP_BOOL  = true;
		}
		else {
			// revert the polymer 
			this->particle_swap_with_monomer(ne_list[choice], p_idx, 0, lattice_index(loc_0, this->y, this->z));
			this->IMP_BOOL = false;
			this->check_structures();
			this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
		}
		std::cout << "Final polymer config:" << std::endl;
		this->Polymers[p_idx].print_chain();
	}
	else {
		this->IMP_BOOL = false; 
	}

	std::cout << "Done with tail rotation!" << std::endl;
	return;

}

void Simulation::debug_head_rotation(int p_idx){

	// get the neighborlist of the particle at index -2
	int deg_poly = this->Polymers[p_idx].deg_poly;
	std::array<int,3> loc_0 = this->Polymers[p_idx].chain[deg_poly-1]->coords;
	std::array<int,3> loc_1 = this->Polymers[p_idx].chain[deg_poly-2]->coords;

	double                final_E        = 0; // energy of the final configuration of the system
	std::array <double,8> final_contacts = {0, 0, 0, 0, 0, 0, 0, 0};

	// generate a neighbor list 
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

	// get the choice for a neighbor
	int choice = rng_uniform(0, 25);

	// instantiate some containers
	Container container;

	if ( this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->ptype[0] == 's' ){

		this->particle_swap_with_monomer_update(&container, loc_0, p_idx, deg_poly-1, lattice_index(ne_list[choice], this->y, this->z));
		this->check_structures();
		this->Polymers[p_idx].print_chain();
		
		// doing the quick manipulations to get the final energy and contacts
		final_E = this->sysEnergy - (container.E1_initial + container.E2_initial - container.Epair_initial); 
		final_E = final_E + (container.E1_final   + container.E2_final  - container.Epair_final);
		final_contacts = add_arrays(subtract_arrays(this->contacts, add_arrays (container.c1_initial, container.c2_initial)), add_arrays (container.c1_final, container.c2_final));
		final_contacts = add_arrays(final_contacts, container.cpair_initial);
		final_contacts = subtract_arrays(final_contacts, container.cpair_final);

		// run the debugger
		std::cout << "Running the debugger..." << std::endl; 
		this->debug_checks_energy_contacts(final_E, final_contacts);
		std::cout << "Ran the debugger!" << std::endl;

		if (metropolis_acceptance(this->sysEnergy, final_E, this->T)){
			this->sysEnergy = final_E;
			this->contacts  = final_contacts;
			this->IMP_BOOL  = true;
		}
		else {
			// revert the polymer 
			this->particle_swap_with_monomer(ne_list[choice], p_idx, deg_poly-1, lattice_index(loc_0, this->y, this->z));
			this->IMP_BOOL = false;
			this->check_structures();
			this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
		}
		std::cout << "Final polymer config:" << std::endl;
		this->Polymers[p_idx].print_chain();
	}
	else {
		this->IMP_BOOL = false; 
	}
	return;

}

void Simulation::debug_reptation_forward(int p_idx){

	int                  deg_poly   = this->Polymers[p_idx].deg_poly;
	std::array<int,3>    loc_i      = this->Polymers[p_idx].chain[0]->coords;
	std::array<int,3>    loc_f      = this->Polymers[p_idx].chain[deg_poly-1]->coords;
	std::array<double,8> contacts_i = this->contacts; 
	double               E_i        = this->sysEnergy;

	// get the possibilities of how the polymer can slither forward
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_f, this->x, this->y, this->z);
	int choice = rng_uniform(0, 25);

	// generate the position to slither to, and its copy!
	std::array <int,3> to_slither      = ne_list[choice];
	std::array <int,3> to_slither_copy = to_slither;

	if (to_slither == loc_i){
		std::cout << "With tailbiting." << std::endl;
		this->forward_reptation_with_tail_biting(&contacts_i, &E_i, p_idx);
	}
	else if (this->Lattice[lattice_index(to_slither, this->y, this->z)]->ptype[0] == 's'){
		std::cout << "No tailbiting." << std::endl;
		this->forward_reptation_without_tail_biting(&contacts_i, &to_slither, &E_i, p_idx);
	}
	else {
		this->IMP_BOOL = false;
		return;
	}
	this->check_structures(); 

	// printing out polymer
	this->Polymers[p_idx].print_chain();

	// run the debugger
	std::cout << "Running the debugger..." << std::endl; 
	this->debug_checks_energy_contacts(E_i, contacts_i);
	std::cout << "Ran the debugger!" << std::endl;

	if (metropolis_acceptance(this->sysEnergy, E_i, this->T)){
		this->sysEnergy = E_i;
		this->contacts  = contacts_i;
		this->IMP_BOOL  = true;
	}
	else {
		if (to_slither_copy == loc_i){
			this->revert_with_head_butting(&to_slither_copy, p_idx);
		}
		else{
			// instantiate a nullptr...
			Particle* tmp {nullptr};
			tmp = this->Lattice[lattice_index(loc_i, y, z)];
			this->revert_without_head_butting(tmp, &loc_i, &to_slither_copy, p_idx);
		}
		this->IMP_BOOL = false;
	}

	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	this->check_structures();

	return;

}

void Simulation::debug_reptation_backward(int p_idx){

	int                   deg_poly         = this->Polymers[p_idx].deg_poly;
	std::array <int,3>    loc_i            = this->Polymers[p_idx].chain[0]->coords;
	std::array <int,3>    loc_f            = this->Polymers[p_idx].chain[deg_poly-1]->coords;
	std::array <double,8> contacts_i       = this->contacts;
	double                E_i              = this->sysEnergy;

	// first check if tail rotation can be performed at all 
	std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( loc_i, x, y, z );
	int choice = rng_uniform (0, 25);

	std::array <int,3> to_slither      = ne_list[choice]; 
	std::array <int,3> to_slither_copy = to_slither; 


	if ( to_slither == loc_f ){
		std::cout << "With head butting..." << std::endl;
		this->backward_reptation_with_head_butting (&contacts_i, &E_i, p_idx); 
	}
	else if ( (this->Lattice)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		std::cout << "Without head butting... " << std::endl;
		this->backward_reptation_without_head_butting (&contacts_i, &to_slither, &E_i, p_idx);
	}

	else {
		this->IMP_BOOL = false; 
		return; 
	}

	this->check_structures();

	// printing out polymer
	this->Polymers[p_idx].print_chain();

	// run the debugger
	std::cout << "Running the debugger..." << std::endl;
	this->debug_checks_energy_contacts(E_i, contacts_i);
	std::cout << "Ran the debugger!" << std::endl;

	if ( metropolis_acceptance (this->sysEnergy, E_i, this->T) ){
		this->sysEnergy = E_i;
		this->contacts  = contacts_i;
	}
	else {
		// revert back to old state 
		if ( to_slither_copy == loc_f ) {
			this->revert_with_tail_biting (&to_slither_copy, p_idx);
		}
		else {
			Particle* tmp {nullptr};
			tmp = this->Lattice[ lattice_index (loc_f, y, z) ]; 
			this->revert_without_tail_biting (tmp, &loc_f, &to_slither_copy, p_idx); 
		}
		this->IMP_BOOL = false; 
	}

	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	this->check_structures();

	return;
}

void Simulation::debug_polymer_orientation_flip(int p_idx){

	// instantiate some useful variables
	int deg_poly = this->Polymers[p_idx].deg_poly;
	std::vector<int> polymer_indices (deg_poly);
	std::iota(polymer_indices.begin(), polymer_indices.end(), 0);

	// set up the random number generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(polymer_indices.begin(), polymer_indices.end(), std::default_random_engine(seed));

	// instantiate variables for boltzmann sampling
	int ntest    = 5;
	int nflip    = (deg_poly==1) ? 1 : rng_uniform (1, deg_poly-1);

	EnhancedOrientationFlipper enhanced_flipper(nflip, ntest); 

	double               E_sys        = this->sysEnergy;
	std::array<double,8> contacts_sys = this->contacts;

	double               E_post        = 0; // the energy of the system after final perturbation
	std::array<double,8> contacts_post = {0,0,0,0,0,0,0,0}; // the contacts of the system after final perturbation
	double               rng_acc       = 0; // rng for acceptance at the very end

	// relevant indices variable store
	int m_lattice_idx = 0; 

	// loop over the different monomer indices
	for ( int i{0}; i < nflip; ++i ){

		// sample different orientations and get the boltzmann factors 
		this->debug_orientation_sampler_forwards(&enhanced_flipper, &contacts_sys, E_sys, i, p_idx, polymer_indices[i]);

		// get the energy minima
		enhanced_flipper.Emin = *std::min_element(enhanced_flipper.energies.begin(), enhanced_flipper.energies.end());

		// go to the new state 
		this->debug_choose_state_forward(&enhanced_flipper, i, p_idx, polymer_indices[i]);

		// set the system energy at the final point
		E_sys        = enhanced_flipper.energies[enhanced_flipper.sampler_idx];
		contacts_sys = enhanced_flipper.contacts_store[enhanced_flipper.sampler_idx];

		// run the check
		std::cout << "Second flip test..." << std::endl;
		this->debug_checks_energy_contacts(E_sys, contacts_sys);

		// reset
		enhanced_flipper.initial_contacts = {0, 0, 0, 0, 0, 0, 0, 0};
		enhanced_flipper.initial_E        = 0;
		enhanced_flipper.rboltzmann       = 0;
		enhanced_flipper.sampler_rsum     = 0;

	}

	E_post        = E_sys;
	contacts_post = contacts_sys;

	this->debug_checks_energy_contacts(E_sys, contacts_sys);
	enhanced_flipper.perturbed_contacts = {0, 0, 0, 0, 0, 0, 0, 0};
	enhanced_flipper.perturbed_E        = 0;


	for ( int i{0}; i < nflip; ++i ){

		// sample different orientations and get the boltzmann factors for the backward flux
		this->debug_orientation_sampler_backwards(&enhanced_flipper, &contacts_sys, E_sys, i, p_idx, polymer_indices[i]);

		// get the minimum energy
		enhanced_flipper.Emin = *std::min_element(enhanced_flipper.energies.begin(), enhanced_flipper.energies.end());

		// get the backwards probability flux
		for (int k{0}; k < ntest; ++k){
			enhanced_flipper.boltzmann [k] = std::exp (-1/this->T*(enhanced_flipper.energies[k] - enhanced_flipper.Emin));
			enhanced_flipper.rboltzmann   += enhanced_flipper.boltzmann[k];
		}
		enhanced_flipper.prob_n_to_o      *= enhanced_flipper.boltzmann[0]/enhanced_flipper.rboltzmann;

		// make the jump to the old state 
		this->Polymers[p_idx].chain[polymer_indices[i]]->orientation = enhanced_flipper.initial_orientations[i];

		// set the energetics and contacts right
		E_sys        = enhanced_flipper.energies[0];
		contacts_sys = enhanced_flipper.contacts_store[0];
		std::cout << "Sixth flip test..." << std::endl;
		this->debug_checks_energy_contacts(E_sys, contacts_sys);

		// reset
		enhanced_flipper.initial_contacts = {0, 0, 0, 0, 0, 0, 0 ,0};
		enhanced_flipper.initial_E        = 0;
		enhanced_flipper.rboltzmann       = 0;

	}

	std::cout << "Seventh flip test..." << std::endl;
	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);

	// get a random number
	rng_acc = rng_uniform(0.0, 1.0);

	// run the criterion
	if ( rng_acc < std::exp (-1/this->T * (E_post - this->sysEnergy ) * enhanced_flipper.prob_n_to_o/enhanced_flipper.prob_o_to_n)){
		std::cout << "Orientation to be assigned = "; print(enhanced_flipper.final_orientations);
		for (int j{0}; j < nflip; ++j){
			(this->Polymers)[p_idx].chain[polymer_indices[j]]->orientation = enhanced_flipper.final_orientations[j];
		}
		this->sysEnergy = E_post;
		this->contacts  = contacts_post;
		this->IMP_BOOL  = true;
	}
	else {
		this->IMP_BOOL = false; 
	}

	std::cout << "Eighth flip test..." << std::endl;
	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);

	return;
}

void Simulation::debug_solvation_shell_flip(){

	// get solvation shell of the polymers 
	std::set   <int> solvation_shell_set = this->get_solvation_shell(); 

	// get the solvation shell indices in a vector 
	std::vector <int> solvation_shell_indices (solvation_shell_set.begin(), solvation_shell_set.end()); 

	// shuffle the solvation shell 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));

	// choose how many particle orientations to flip
	int ntest = 5;
	int nflip = rng_uniform(1, static_cast<int>(solvation_shell_indices.size()/2 ) );     // number of sites to be flipped 

	EnhancedOrientationFlipper enhanced_flipper (nflip, ntest);

	double               E_sys        = this->sysEnergy;
	std::array<double,8> contacts_sys = this->contacts;

	double               E_post        = 0; // the energy of the system after final perturbation
	std::array<double,8> contacts_post = {0,0,0,0,0,0,0,0}; // the contacts of the system after final perturbation
	double               rng_acc       = 0; // rng for acceptance at the very end

	for (int i{0}; i < nflip; ++i){

		// sample different orientations and get the boltzmann factors 
		this->debug_orientation_sampler_forwards(&enhanced_flipper, &contacts_sys, E_sys, i, solvation_shell_indices[i]);

		/*
		// get the flip index 
		rboltzmann = 0; 
		old_ori.push_back((this->Lattice)[solvation_shell_indices[i]]->orientation);

		// find the neighboring interaction energies 
		this->neighbor_energetics(solvation_shell_indices[i], &contacts_i, &Ei); 
		std::cout << "Inside first debug test..." << std::endl;
		for (int j{0}; j < ntest; ++j){
			(this->Lattice)[solvation_shell_indices[i]]->orientation = rng_uniform (0, 25);
			orientations [j]    = (this->Lattice)[solvation_shell_indices[i]]->orientation;
			this->neighbor_energetics(solvation_shell_indices[i], &contacts_pert, &Epert);
			energies [j]        = Esys - Ei + Epert;
			contacts_store [j]  = subtract_arrays(&contacts_sys,      &contacts_i);
			contacts_store [j]  = add_arrays     (&contacts_store[j], &contacts_pert);
			contacts_pert       = {0, 0, 0, 0, 0, 0, 0, 0};
			Epert               = 0;
			this->debug_checks_energy_contacts(energies[j], contacts_store[j]);
		}
		*/
		// std::cout << "Energies are "; print(energies);
		enhanced_flipper.Emin = *std::min_element ( enhanced_flipper.energies.begin(), enhanced_flipper.energies.end() );

		for (int k{0}; k < ntest; ++k){
			enhanced_flipper.boltzmann[k] = std::exp (-1/this->T*(enhanced_flipper.energies[k] - enhanced_flipper.Emin));
			enhanced_flipper.rboltzmann  += enhanced_flipper.boltzmann[k];
		}

		rng     = rng_uniform (0.0, 1.0);
		rsum    = 0;
		e_idx   = 0;

		for (int j{0}; j < ntest; ++j){
			rsum += boltzmann[j]/rboltzmann; 
			if (rng < rsum){
				e_idx = j; 
				break; 
			}
		}

		// make the jump to the new state 
		std::cout << "orientation to pick = "; print(orientations);
		new_ori.push_back (orientations[e_idx]);
		this->Lattice[solvation_shell_indices[i]]->orientation = orientations[e_idx];
		prob_o_to_n *= boltzmann[e_idx]/rboltzmann;
		Esys         = energies[e_idx];
		contacts_sys = contacts_store[e_idx];
		std::cout << "Second flip test..." << std::endl;
		this->debug_checks_energy_contacts(energies[e_idx], contacts_store[e_idx]);

		// reset
		contacts_i = {0, 0, 0, 0, 0, 0, 0, 0};
		Ei         = 0;
	}

	frontflow_energy   = energies[e_idx];
	frontflow_contacts = contacts_store[e_idx];
	std::cout << "Third flip test..." << std::endl;
	contacts_i    = {0, 0, 0, 0, 0, 0, 0, 0};
	contacts_pert = {0, 0, 0, 0, 0, 0, 0, 0};
	Ei            = 0;
	Epert         = 0;

	for (int i{0}; i<nflip; ++i){

		rboltzmann = 0;
		this->neighbor_energetics(solvation_shell_indices[i], &contacts_i, &Ei);

		// first iteration
		this->Lattice[solvation_shell_indices[i]]->orientation = old_ori[i];
		this->neighbor_energetics(solvation_shell_indices[i], &contacts_pert, &Epert);
		energies[0] = Esys - Ei + Epert;
		contacts_store[0] = subtract_arrays(&contacts_sys, &contacts_i);
		contacts_store[0] = add_arrays(&contacts_store[0], &contacts_pert);

		std::cout << "Inside fourth set of flip test..." << std::endl;
		this->debug_checks_energy_contacts(energies[0], contacts_store[0]);

		contacts_pert = {0, 0, 0, 0, 0, 0, 0, 0};
		Epert         = 0;
		std::cout << "Inside fifth set of flip test..." << std::endl;

		for (int j{1}; j <ntest; ++j){
			this->Lattice[solvation_shell_indices[i]]->orientation = rng_uniform(0, 25);
			this->neighbor_energetics(solvation_shell_indices[i], &contacts_pert, &Epert);
			energies[j] = Esys - Ei + Epert;
			contacts_store[j] = subtract_arrays(&contacts_sys, &contacts_i);
			contacts_store[j] = add_arrays(&contacts_store[j], &contacts_pert);
			contacts_pert     = {0,0,0,0,0,0,0,0};
			Epert             = 0;
			this->debug_checks_energy_contacts(energies[j], contacts_store[j]);
		}

		Emin = *std::min_element(energies.begin(), energies.end()); 

		for (int k{0}; k<ntest; ++k){
			boltzmann[k] = std::exp(-1/this->T*(energies[k]-Emin));
			rboltzmann  += boltzmann[k]; 
		}

		prob_n_to_o += boltzmann[0]/rboltzmann; 
		this->Lattice[solvation_shell_indices[i]]->orientation = old_ori[i]; 
		Esys         = energies[0];
		contacts_sys = contacts_store[0];
		std::cout << "Sixth flip test..." << std::endl;
		contacts_i   = {0, 0, 0, 0, 0, 0, 0, 0};
		Ei           = 0;

	}
	this->debug_checks_energy_contacts(energies[j], contacts_store[j]);

	rng_acc = rng_uniform(0.0, 1.0);
	if ( rng_acc < std::exp (-1/this->T * (frontflow_energy - this->sysEnergy)) * prob_n_to_o/prob_o_to_n){

		// if accepted, return to the new orientations 
		for (int j{0}; j < nflip; ++j){
			this->Lattice[solvation_shell_indices[j]]->orientation = new_ori[j]; 
		}

		this->sysEnergy = frontflow_energy;
		this->contacts  = frontflow_contacts;
		this->IMP_BOOL  = true;

	}
	else {
		this->IMP_BOOL = false;
	}

	std::cout << "Seventh flip test..." << std::endl;
	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);

	return;
}

void Simulation::debug_lattice_flip(){

	int nflips = rng_uniform(1, static_cast<int>(this->x * this->y * this->z / 8));
	std::array <double,8> contacts = this->contacts; 

	std::vector <int> orientations; 

	std::array <double,8> cs_i = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cs_f = {0, 0, 0, 0, 0, 0, 0, 0};

	double Es_i = 0;
	double Es_f = 0;
	double Ef   = this->sysEnergy; 

	std::vector<int> samples(this->x * this->y * this->z - 1); 
	std::iota(samples.begin(), samples.end(), 0);

	std::random_device rd;
	std::mt19937 g(rd()); 
	std::shuffle(samples.begin(), samples.end(), g); 

	for (int i{0}; i<nflips; ++i){
		
		this->neighbor_energetics(samples[i], &cs_i, &Es_i); 
		orientations.push_back(this->Lattice[samples[i]]->orientation); 
		this->Lattice[samples[i]]->orientation = rng_uniform(0, 25);
		this->neighbor_energetics(samples[i], &cs_f, &Es_f); 
		contacts = add_arrays(subtract_arrays(contacts, cs_i), cs_f);
		Ef += Es_f - Es_i;
		this->debug_checks_energy_contacts(Ef, contacts);

		// resetting...
		cs_i = {0, 0, 0, 0, 0, 0, 0, 0};
		cs_f = {0, 0, 0, 0, 0, 0, 0, 0};
		Es_i = 0;
		Es_f = 0;

	}

	if (metropolis_acceptance(this->sysEnergy, Ef, this->T)){
		this->sysEnergy = Ef; 
		this->contacts  = contacts; 
		this->IMP_BOOL  = true;
	}
	else {
		for (int i{0}; i<nflips; ++i){
			this->Lattice[samples[i]]->orientation = orientations[i];
		}
		this->IMP_BOOL = false;
	}

	return;
}

void Simulation::debug_solvent_exchange_from_shell(){

	std::set   <int>                   solvation_shell_set;
	std::array <std::array<int,3>, 26> ne_list, ne_list_  ;

	// get the first solvation shell
	for (Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain) {
			ne_list = obtain_ne_list (p->coords, this->x, this->y, this->z);
			for ( std::array <int,3>& loc: ne_list ){
				if ( (this->Lattice) [ lattice_index (loc, y, z) ]->ptype[0] == 's' ) {
					solvation_shell_set.insert (lattice_index (loc, y, z));
					ne_list_ = obtain_ne_list (loc, this->x, this->y, this->z);
					for ( std::array <int,3>& loc_: ne_list_){
						if ( (this->Lattice)[ lattice_index (loc_, y, z) ]->ptype[0] == 's' ){
							solvation_shell_set.insert ( lattice_index (loc_, y, z) );
						}
					}
				}
			}
		}
	}

	std::vector <int> solvation_shell_indices (solvation_shell_set.begin(), solvation_shell_set.end());

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));

	Container container;
	/*
	std::array <double,8> cs1_i              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_i              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i            = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs1_f              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_f              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_f            = {0,0,0,0,0,0,0,0};
	*/
	std::array <double,8> frontflow_contacts = {0,0,0,0,0,0,0,0};

	/*
	double Es1_i            = 0;
	double Es1_f            = 0;
	double Es2_i            = 0;
	double Es2_f            = 0;
	double Epair_i          = 0;
	double Epair_f          = 0;
	*/
	double frontflow_energy = 0;
	double rng_acc          = 0;

	Particle* tmp_par_ptr {nullptr};

	int exc_idx = 0; 
	int nswitch = rng_uniform(1, static_cast<int>(solvation_shell_indices.size()/4));

	for (int n{0}; n<nswitch; ++n){
		
		exc_idx = rng_uniform(0, this->x * this->y * this->z - 1);
		if (exc_idx == solvation_shell_indices[n] || this->Lattice[exc_idx]->ptype[0] == 'm'){
			continue;
		}
		
		else {

			this->particle_swap_with_update(&container, tmp_par_ptr, solvation_shell_indices[n], exc_idx);

			/*
			this->neighbor_energetics(solvation_shell_indices[n], &cs1_i, &Es1_i);
			this->neighbor_energetics(exc_idx, &cs2_i, &Es2_i);
			this->selected_pair_interaction(this->Lattice[solvation_shell_indices[n]], this->Lattice[exc_idx], &cpair_i, &Epair_i); 

			tmp_par_ptr = (this->Lattice)[ solvation_shell_indices[n]];
			
			(this->Lattice)[(solvation_shell_indices)[n]] = (this->Lattice)[exc_idx];
			(this->Lattice)[(solvation_shell_indices)[n]]->coords = location ((solvation_shell_indices)[n], this->x, this->y, this->z);

			(this->Lattice)[exc_idx] = tmp_par_ptr;
			(this->Lattice)[exc_idx]->coords = location (exc_idx, this->x, this->y, this->z);

			// flip the particle newly added to the solvation shell
			this->neighbor_energetics (solvation_shell_indices[n], &cs1_f, &Es1_f);
			this->neighbor_energetics (exc_idx, &cs2_f, &Es2_f); 
			this->selected_pair_interaction(this->Lattice[solvation_shell_indices[n]], this->Lattice[exc_idx], &cpair_f, &Epair_f); 
			*/

			// update energy and contacts
			frontflow_energy   = this->sysEnergy  - (container.E1_initial + container.E2_initial - container.Epair_initial); 
			frontflow_energy   = frontflow_energy + (container.E1_final   + container.E2_final  - container.Epair_final);
			frontflow_contacts = add_arrays(subtract_arrays(this->contacts, add_arrays (container.c1_initial, container.c2_initial)), add_arrays (container.c1_final, container.c2_final));
			frontflow_contacts = add_arrays(frontflow_contacts, container.cpair_initial);
			frontflow_contacts = subtract_arrays(frontflow_contacts, container.cpair_final);
			
			this->debug_checks_energy_contacts(frontflow_energy, frontflow_contacts);

			rng_acc = rng_uniform (0.0, 1.0);

			if (rng_acc < std::exp (-1/this->T * (frontflow_energy - this->sysEnergy))){
				this->sysEnergy = frontflow_energy;
				this->contacts  = frontflow_contacts;
				this->IMP_BOOL  = true; 
			}
			else {
				this->particle_swap(tmp_par_ptr, solvation_shell_indices[n], exc_idx);
				/*
				tmp_par_ptr = (this->Lattice)[solvation_shell_indices[n]];
				(this->Lattice)[solvation_shell_indices [n]]         = (this->Lattice)[exc_idx]; 
				(this->Lattice)[solvation_shell_indices [n]]->coords = location (solvation_shell_indices[n], this->x, this->y, this->z);

				(this->Lattice)[exc_idx]         = tmp_par_ptr;
				(this->Lattice)[exc_idx]->coords = location (exc_idx, this->x, this->y, this->z);
				*/
				this->IMP_BOOL = false; 
			}

			this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);

			// resetting 
			container.reset();

		}

	}
	return;
}

void Simulation::debug_solvent_exchange(){

	int idx1 = -1;
	int idx2 = -1; 
	int nswitches = 50; 
	double frontflow_energy = 0;
	double rng_acc          = 0;
	std::array <double,8> frontflow_contacts = this->contacts;

	Container container;
	Particle* tmp_par_ptr {nullptr};

	for (int j{0}; j<nswitches; ++j){

		idx1 = rng_uniform(0, this->x * this->y * this->z - 1); 
		idx2 = rng_uniform(0, this->x * this->y * this->z - 1); 

		if (idx1 == idx2 || this->Lattice[idx1]->ptype == "m1" || this->Lattice[idx2]->ptype == "m1"){
			continue; 
		}
		else {
			this->particle_swap_with_update(&container, tmp_par_ptr, idx1, idx2);
			/*
			this->neighbor_energetics(idx1, &cs1_i, &Es1_i);
			this->neighbor_energetics(idx2, &cs2_i, &Es2_i);
			this->selected_pair_interaction((this->Lattice)[idx1], (this->Lattice)[idx2], &cpair_i, &Epair_i); 

			tmp_par_ptr = (this->Lattice)[idx1];
			
			(this->Lattice)[idx1] = (this->Lattice)[idx2];
			(this->Lattice)[idx1]->coords = location (idx1, this->x, this->y, this->z);

			(this->Lattice)[idx2] = tmp_par_ptr;
			(this->Lattice)[idx2]->coords = location(idx2, this->x, this->y, this->z);

			// flip the particle newly added to the solvation shell
			this->neighbor_energetics (idx1, &cs1_f, &Es1_f);
			this->neighbor_energetics (idx2, &cs2_f, &Es2_f); 
			this->selected_pair_interaction(this->Lattice[idx1], this->Lattice[idx2], &cpair_f, &Epair_f); 
			*/ 
			frontflow_energy   = this->sysEnergy - (container.E1_initial + container.E2_initial - container.Epair_initial); 
			frontflow_energy   = frontflow_energy + (container.E1_final   + container.E2_final  - container.Epair_final);
			frontflow_contacts = add_arrays(subtract_arrays(this->contacts, add_arrays (container.c1_initial, container.c2_initial)), add_arrays (container.c1_final, container.c2_final));
			frontflow_contacts = add_arrays(frontflow_contacts, container.cpair_initial);
			frontflow_contacts = subtract_arrays(frontflow_contacts, container.cpair_final);

			this->debug_checks_energy_contacts(frontflow_energy, frontflow_contacts);

			rng_acc = rng_uniform (0.0, 1.0); 
			if (rng_acc < std::exp (-1/this->T * (frontflow_energy - this->sysEnergy))){
				this->sysEnergy = frontflow_energy;
				this->contacts  = frontflow_contacts;
		
			}
			else {
				this->particle_swap(tmp_par_ptr, idx1, idx2);
				/*
				tmp_par_ptr = (this->Lattice)[idx1];
				(this->Lattice)[idx1]         = (this->Lattice)[idx2]; 
				(this->Lattice)[idx1]->coords = location(idx1, this->x, this->y, this->z);

				(this->Lattice)[idx2]         = tmp_par_ptr;
				(this->Lattice)[idx2]->coords = location (idx2, this->x, this->y, this->z);
				*/
			}

			container.reset();

		}

	}

	return;
}

void Simulation::debug_regrowth(int p_idx){

	int deg_poly        = this->Polymers[p_idx].deg_poly; 
	int m_idx           = rng_uniform(1, deg_poly-2);
	int recursion_depth = 0;
	int growth          = -1;

	double rng_acc     {0};
	double prob_o_to_n {1};
	double prob_n_to_o {1};
	double back_energy {0};
	double forw_energy {this->sysEnergy};

	std::array <double,8> forw_contacts = this->contacts;
	std::array <double,8> back_contacts = {0,0,0,0,0,0,0,0};

	std::vector <std::array<int,3>> old_cut;
	std::vector <std::array<int,3>> new_cut;
	old_cut.reserve(deg_poly);
	new_cut.reserve(deg_poly);

	if (deg_poly%2 == 0){
		growth = (0.5 >= (m_idx+1)/static_cast<double>(deg_poly)) ? 0 : 1;
	}
	else {
		if ( 0.5 == (m_idx+1)/static_cast<double>(deg_poly+1) ){
			growth = rng_uniform (0, 1);
		}
		else {
			growth = (0.5 > (m_idx+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
		}
	}

	if (growth){

		std::cout << "Performing head regrowth..." << std::endl;
		for (int i{m_idx+1}; i<deg_poly; ++i){
			old_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		std::cout << "\tPerforming the forward pass..." << std::endl;
		this->debug_forward_head_regrowth(&forw_contacts, &prob_o_to_n, &forw_energy, p_idx, m_idx);
		std::cout << "\tPerforming the forward pass check..." << std::endl;
		this->Polymers[p_idx].print_chain();
		this->debug_checks_energy_contacts(forw_energy, forw_contacts);
		this->check_structures();
		std::cout << "done!" << std::endl;

		for (int i{m_idx+1}; i<deg_poly; ++i){
			new_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		if (old_cut == new_cut){
			// do nothing
			;
		}

		else if (!(this->IMP_BOOL)){
			this->debug_accept_after_head_regrowth(&new_cut, &old_cut);
			this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
			this->check_structures();
		}

		else {
			back_contacts = forw_contacts;
			back_energy   = forw_energy;
			std::cout << "\tPerforming the backward pass..." << std::endl;
			this->debug_backward_head_regrowth(&old_cut, &back_contacts, &prob_n_to_o, &back_energy, p_idx, m_idx, recursion_depth);

			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * prob_n_to_o / prob_o_to_n){

				this->debug_accept_after_head_regrowth(&old_cut, &new_cut);
				this->sysEnergy = forw_energy; 
				this->contacts = forw_contacts;
				this->IMP_BOOL = true;
			}
			else {
				this->IMP_BOOL = false;
			}
		}
		std::cout << "\tPerforming the post back flow check..." << std::endl;
		this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
		this->check_structures();
		std::cout << "done!" << std::endl;
	}

	else {

		std::cout << "\tPerforming forward tail regrowth..." << std::endl;
		for (int i{0}; i<m_idx; ++i){
			old_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		this->debug_forward_tail_regrowth(&forw_contacts, &prob_o_to_n, &forw_energy, p_idx, m_idx);
		this->debug_checks_energy_contacts(forw_energy, forw_contacts);
		this->check_structures();

		for (int i{0}; i<m_idx; ++i){
			new_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		if (old_cut == new_cut){
			// do nothing
			;
		}

		else if (!(this->IMP_BOOL)){
			this->debug_accept_after_tail_regrowth(&new_cut, &old_cut);
		}

		else {

			back_contacts = forw_contacts; 
			back_energy = forw_energy;
			std::cout << "\tPerforming backward tail regrowth..." << std::endl;
			debug_backward_tail_regrowth(&old_cut, &back_contacts, &prob_n_to_o, &back_energy, p_idx, m_idx, recursion_depth);

			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * prob_n_to_o/prob_o_to_n) {

				this->debug_accept_after_tail_regrowth(&old_cut, &new_cut);
				this->sysEnergy = forw_energy;
				this->contacts  = forw_contacts;
			}

			else {
				this->IMP_BOOL = false;
			}
		}
		this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
		this->check_structures();
	}

	return;
}

void Simulation::debug_forward_head_regrowth(std::array <double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx){

	int deg_poly   = this->Polymers[p_idx].deg_poly; 
	if (m_idx == deg_poly-1) {
		this->IMP_BOOL = true;
		return; 
	}

	std::array<double,8> current_contacts = *forw_contacts;

	std::array<double,5>               energies;
	std::array<std::array<double,8>,5> contacts_store;
	std::array<double,5>               boltzmann;
	std::array<int,3>                  loc_m = this->Polymers.at(p_idx).chain.at(m_idx+1)->coords;

	std::array<std::array<int,3>,26> ne_list = obtain_ne_list(this->Polymers.at(p_idx).chain.at(m_idx)->coords, this->x, this->y, this->z);
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
	std::shuffle(ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// set up the initial energies
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em_i       = 0;
	double Es_i       = 0;
	double Em_n       = 0;
	double Es_n       = 0;
	double Epair_i    = 0;
	double Epair_n    = 0;
	double Esys       = *forw_energy;

	// set the contacts ready
	// doubles for energy transfer 
	std::array <double,8> cm_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps
	int block_counter = 0 ;
	int idx_counter   = 0 ;
	int self_swap_idx = -1;

	while (idx_counter<5) {
		std::cout << "While index counter is " << idx_counter << "." << std::endl;
		std::cout << "Initial config:" << std::endl;
		this->Polymers[p_idx].print_chain();
		if (ne_list[idx_counter] == loc_m){
			energies[idx_counter]       = Esys;
			contacts_store[idx_counter] = current_contacts;
			block_counter              += 1;
		}
		else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm') {

			for (int u{0}; u<deg_poly; ++u){
				if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
					self_swap_idx = u;
					break;
				}
			}

			std::cout << "self swap monomer index = " << self_swap_idx << "." << std::endl;
			if (self_swap_idx <= m_idx){
				energies[idx_counter]       = 1e+8;
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1};
				block_counter              += 1;
			}

			else {
				std::cout << "Switching up monomer-monomer particles..." << std::endl;

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

				// swap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(contacts_store[idx_counter]);

				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(energies[idx_counter], contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// resetting...
				cm_i    = {0,0,0,0,0,0,0,0};
				cs_i    = {0,0,0,0,0,0,0,0};
				cpair_i = {0,0,0,0,0,0,0,0};
				cm_n    = {0,0,0,0,0,0,0,0};
				cs_n    = {0,0,0,0,0,0,0,0};
				cpair_n = {0,0,0,0,0,0,0,0};
				Em_i    = 0;
				Es_i    = 0;
				Em_n    = 0;
				Es_n    = 0;
				Epair_i = 0;
				Epair_n = 0;
			}
		}
		else {
			std::cout << "Switching up monomer-solvent particles..." << std::endl;
			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

			// swap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);
			std::cout << "contacts_store[" << idx_counter << "] = "; print(contacts_store[idx_counter]);

			this->Polymers[p_idx].print_chain();
			this->debug_checks_energy_contacts(energies[idx_counter], contacts_store[idx_counter]);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// resetting...
			cm_i    = {0,0,0,0,0,0,0,0};
			cs_i    = {0,0,0,0,0,0,0,0};
			cpair_i = {0,0,0,0,0,0,0,0};
			cm_n    = {0,0,0,0,0,0,0,0};
			cs_n    = {0,0,0,0,0,0,0,0};
			cpair_n = {0,0,0,0,0,0,0,0};
			Em_i    = 0;
			Es_i    = 0;
			Em_n    = 0;
			Es_n    = 0;
			Epair_i = 0;
			Epair_n = 0;
		}
		idx_counter += 1;
	}

	if (block_counter == 5){
		this->IMP_BOOL = false;
	}
	else {

		double Emin = *std::min_element(energies.begin(), energies.end());
		for (int i{0}; i<5; ++i){
			boltzmann[i] = std::exp(-1/this->T * (energies[i]-Emin));
			rboltzmann  += boltzmann[i];
		}

		double rng_acc = rng_uniform(0.0, 1.0);
		double rsum    = 0;
		int    e_idx   = 0;

		for (int j{0}; j<5; ++j){
			rsum += boltzmann[j]/rboltzmann;
			if (rng_acc < rsum){
				e_idx = j;
				break;
			}
		}

		*prob_o_to_n  *= boltzmann[e_idx]/rboltzmann;
		*forw_energy   = energies[e_idx];
		*forw_contacts = contacts_store[e_idx];

		// do the swap again and go to the new configuration
		this->perturb_particle_swap(lattice_index(ne_list[e_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->Polymers[p_idx].print_chain();
		this->debug_checks_energy_contacts(energies[e_idx], contacts_store[e_idx]);
		this->debug_forward_head_regrowth(forw_contacts, prob_o_to_n, forw_energy, p_idx, m_idx+1);
	}

	return;
}

void Simulation::debug_accept_after_head_regrowth(std::vector<std::array<int,3>>* old_cut, std::vector<std::array<int,3>>* new_cut){

	// define some stuff
	int L {0};

	std::vector<std::vector<std::array<int,3>>> master_linked_list;
	std::vector<std::array<int,3>> link;

	create_linked_list(*old_cut, *new_cut, link, &master_linked_list, 1);

	Particle* tmp_par_1 {nullptr};
	Particle* tmp_par_2 {nullptr};

	for (std::vector<std::array<int,3>>& linked_list: master_linked_list){
		L = static_cast<int>(linked_list.size()); 
		if (this->Lattice[lattice_index(linked_list[L-1], this->y, this->z)]->ptype[0] =='s'){
			tmp_par_1 = this->Lattice[lattice_index(linked_list.back(), this->y, this->z)];
			// go backward
			for (int i{0}; i<L; i=i+2){
				this->Lattice[lattice_index(linked_list[L-2-i], this->y, this->z)]->coords = linked_list[L-1-i];
				this->Lattice[lattice_index(linked_list[L-1-i], this->y, this->z)] = this->Lattice[lattice_index(linked_list[L-2-i], this->y, this->z)];
			}
			tmp_par_1->coords = linked_list[0];
			this->Lattice[lattice_index(linked_list[0], this->y, this->z)] = tmp_par_1;
		}
		else {
			for (int i{0}; i<L; i=i+2){

				if (i==0){
					tmp_par_1 = this->Lattice[lattice_index(linked_list[i+1], this->y, this->z)];
					this->Lattice[lattice_index(linked_list[i],   this->y, this->z)]->coords = linked_list[i+1];
					this->Lattice[lattice_index(linked_list[i+1], this->y, this->z)] = this->Lattice[lattice_index(linked_list[i], this->y, this->z)];
				}
				else {
					tmp_par_1->coords = linked_list[i+1];
					tmp_par_2 = this->Lattice[lattice_index(linked_list[i+1], this->y, this->z)];
					this->Lattice[lattice_index(linked_list[i+1], this->y, this->z)] = tmp_par_1;
					tmp_par_1 = tmp_par_2;
				}
			}
		}
	}

	return;
}

void Simulation::debug_backward_head_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth){

	int deg_poly = this->Polymers[p_idx].deg_poly;
	if (m_idx == deg_poly-1){
		this->IMP_BOOL = true; 
		return; 
	}

	// these are doubles for energies 
	double Em_i       = 0;
	double Es_i       = 0;
	double Em_n       = 0;
	double Es_n       = 0;
	double Epair_i    = 0;
	double Epair_n    = 0;
	double Esys       = *back_energy;
	double rboltzmann = 0; // running sum for boltzmann weights

	std::array <double,5>               energies;
	std::array <double,5>               boltzmann;
	std::array <double,8>               current_contacts = *back_contacts;
	std::array <std::array<double,8>,5> contacts_store; 

	std::array <int,3> loc_m = (this->Polymers)[p_idx].chain[m_idx+1]->coords; // this is the key item of interest

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((this->Polymers)[p_idx].chain[m_idx]->coords, this->x, this->y, this->z);

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[recursion_depth]) - ne_list.begin();

	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	std::array <double,8> cm_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int idx_counter         = 0;
	int self_swap_idx       = -1;

	while (idx_counter < 5){
		std::cout << "While index counter is " << idx_counter << "." << std::endl;
		std::cout << "Initial config:" << std::endl;
		this->Polymers[p_idx].print_chain();

		if (ne_list[idx_counter] == loc_m){
			energies[idx_counter]       = Esys;
			contacts_store[idx_counter] = current_contacts;
		}
		else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm'){

			for (int u{0}; u<deg_poly; ++u){
				if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
					self_swap_idx = u;
					break;
				}
			}

			if (self_swap_idx <= m_idx){
				energies[idx_counter] = 1e+8; 
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1};
			}
			else {
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

				// swap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(contacts_store[idx_counter]);

				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(energies[idx_counter], contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// resetting...
				cm_i    = {0,0,0,0,0,0,0,0};
				cs_i    = {0,0,0,0,0,0,0,0};
				cpair_i = {0,0,0,0,0,0,0,0};
				cm_n    = {0,0,0,0,0,0,0,0};
				cs_n    = {0,0,0,0,0,0,0,0};
				cpair_n = {0,0,0,0,0,0,0,0};
				Em_i    = 0;
				Es_i    = 0;
				Em_n    = 0;
				Es_n    = 0;
				Epair_i = 0;
				Epair_n = 0;
			}

		}
		else {
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

			// swap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));	

			// resetting...
			cm_i    = {0,0,0,0,0,0,0,0};
			cs_i    = {0,0,0,0,0,0,0,0};
			cpair_i = {0,0,0,0,0,0,0,0};
			cm_n    = {0,0,0,0,0,0,0,0};
			cs_n    = {0,0,0,0,0,0,0,0};
			cpair_n = {0,0,0,0,0,0,0,0};
			Em_i    = 0;
			Es_i    = 0;
			Em_n    = 0;
			Es_n    = 0;
			Epair_i = 0;
			Epair_n = 0;
		}
		idx_counter += 1;
	}

	double Emin = *std::min_element(energies.begin(), energies.end());
	for (int i{0}; i<5; ++i){
		boltzmann[i] = std::exp(-1/this->T * (energies[i]-Emin));
		rboltzmann  += boltzmann[i];
	}

	*prob_n_to_o   = (*prob_n_to_o) * boltzmann[0]/rboltzmann;
	*back_energy   = energies[0];
	*back_contacts = contacts_store[0];

	// perform the swap
	this->perturb_particle_swap(lattice_index(ne_list[0], this->y, this->z), lattice_index(loc_m, this->y, this->z));
	std::cout << "Check at end of regrowth @ recursion depth = " << recursion_depth << "." << std::endl;
	this->Polymers[p_idx].print_chain();
	this->debug_checks_energy_contacts(energies[0], contacts_store[0]);
	this->debug_backward_head_regrowth(old_cut, back_contacts, prob_n_to_o, back_energy, p_idx, m_idx+1, recursion_depth+1);

	return; 
}

void Simulation::debug_forward_tail_regrowth(std::array <double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx){

	int deg_poly = this->Polymers[p_idx].deg_poly;
	if (m_idx == 0){
		this->IMP_BOOL = true; 
		return;
	}

	std::array<double,8> current_contacts = *forw_contacts;

	std::array<double,5>               energies;
	std::array<std::array<double,8>,5> contacts_store;
	std::array<double,5>               boltzmann;
	std::array<int,3>                  loc_m = this->Polymers.at(p_idx).chain.at(m_idx-1)->coords;

	std::array<std::array<int,3>,26> ne_list = obtain_ne_list(this->Polymers.at(p_idx).chain.at(m_idx)->coords, this->x, this->y, this->z);
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
	std::shuffle(ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// set up the initial energies
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em_i       = 0;
	double Es_i       = 0;
	double Em_n       = 0;
	double Es_n       = 0;
	double Epair_i    = 0;
	double Epair_n    = 0;
	double Esys       = *forw_energy;

	// set the contacts ready
	// doubles for energy transfer 
	std::array <double,8> cm_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps
	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 

	while (idx_counter<5) {
		std::cout << "While index counter is " << idx_counter << "." << std::endl;
		std::cout << "Initial config:" << std::endl;
		this->Polymers[p_idx].print_chain();
		if (ne_list[idx_counter] == loc_m){
			energies[idx_counter]       = Esys;
			contacts_store[idx_counter] = current_contacts;
			block_counter              += 1;
		}
		else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm'){

			for (int u{0}; u<deg_poly; ++u){
				if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
					self_swap_idx = u;
					break;
				}
			}

			if (self_swap_idx >= m_idx){
				energies[idx_counter]       = 1e+8;
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1};
				block_counter              += 1;
			}

			else {
				std::cout << "Switching up monomer-monomer particles..." << std::endl;

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

				// swap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(contacts_store[idx_counter]);

				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(energies[idx_counter], contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// resetting
				cm_i    = {0,0,0,0,0,0,0,0};
				cs_i    = {0,0,0,0,0,0,0,0};
				cpair_i = {0,0,0,0,0,0,0,0};
				cm_n    = {0,0,0,0,0,0,0,0};
				cs_n    = {0,0,0,0,0,0,0,0};
				cpair_n = {0,0,0,0,0,0,0,0};
				Em_i    = 0;
				Es_i    = 0;
				Em_n    = 0;
				Es_n    = 0;
				Epair_i = 0;
				Epair_n = 0;

			}
		}
		else {
			std::cout << "Switching up monomer-solvent particles..." << std::endl;
			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

			// swap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);
			std::cout << "contacts_store[" << idx_counter << "] = "; print(contacts_store[idx_counter]);

			this->Polymers[p_idx].print_chain();
			this->debug_checks_energy_contacts(energies[idx_counter], contacts_store[idx_counter]);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// resetting
			cm_i    = {0,0,0,0,0,0,0,0};
			cs_i    = {0,0,0,0,0,0,0,0};
			cpair_i = {0,0,0,0,0,0,0,0};
			cm_n    = {0,0,0,0,0,0,0,0};
			cs_n    = {0,0,0,0,0,0,0,0};
			cpair_n = {0,0,0,0,0,0,0,0};
			Em_i    = 0;
			Es_i    = 0;
			Em_n    = 0;
			Es_n    = 0;
			Epair_i = 0;
			Epair_n = 0;

		}
		idx_counter += 1;
	}

	if (block_counter == 5){
		this->IMP_BOOL = false;
	}
	else {

		double Emin = *std::min_element(energies.begin(), energies.end());
		for (int i{0}; i<5; ++i){
			boltzmann[i] = std::exp(-1/this->T * (energies[i]-Emin));
			rboltzmann  += boltzmann[i];
		}

		double rng_acc = rng_uniform(0.0, 1.0);
		double rsum    = 0;
		int    e_idx   = 0;

		for (int j{0}; j<5; ++j){
			rsum += boltzmann[j]/rboltzmann;
			if (rng_acc < rsum){
				e_idx = j;
				break;
			}
		}

		*prob_o_to_n  *= boltzmann[e_idx]/rboltzmann;
		*forw_energy   = energies[e_idx];
		*forw_contacts = contacts_store[e_idx];

		// do the swap again
		this->perturb_particle_swap(lattice_index(ne_list[e_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->debug_forward_tail_regrowth(forw_contacts, prob_o_to_n, forw_energy, p_idx, m_idx-1);
	}

	return;
}

void Simulation::debug_accept_after_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::vector <std::array<int,3>>* new_cut){

	std::vector< std::vector<std::array<int,3>>> master_linked_list;
	std::vector <std::array<int,3>> link; 

	create_linked_list (*old_cut, *new_cut, link, &master_linked_list, 1);

	int L {0}; // to store length of linked_list 

	Particle* tmp_par_ptr_1 {nullptr}; // to help with swaps of pointers
	Particle* tmp_par_ptr_2 {nullptr}; // this guy as well 

	for ( std::vector <std::array<int,3>>& linked_list: master_linked_list ) {

		L = static_cast<int> (linked_list.size()) ;

		if ( (this->Lattice)[ lattice_index(linked_list[L-1], y, z) ]->ptype[0] == 's' ){

			tmp_par_ptr_1 = (this->Lattice) [ lattice_index (linked_list.back(), y, z) ];

			// go backwards 
			for (int i{0}; i < L ; i=i+2 ){

				(this->Lattice)[lattice_index (linked_list[L-2-i], y, z)]->coords = linked_list[L-1-i]; 
				(this->Lattice)[lattice_index (linked_list[L-1-i], y, z)] = (this->Lattice)[lattice_index( linked_list[L-2-i], y, z )]; 
				
			}
			tmp_par_ptr_1->coords = linked_list[0]; 
			(this->Lattice)[lattice_index (linked_list[0], y, z)] = tmp_par_ptr_1; 

		}

		else { // when it is a circulation of monomers 

			for ( int i{0}; i < L; i=i+2){
				
				if ( i == 0 ) {
					// store info about linked_list[1]... and consequently also linked_list[2]
					tmp_par_ptr_1 = (this->Lattice)[lattice_index ( linked_list[i+1], y, z) ];

					(this->Lattice)[ lattice_index ( linked_list[i],   y, z) ]->coords = linked_list[i+1]; 
					(this->Lattice)[ lattice_index ( linked_list[i+1], y, z) ] = (this->Lattice)[lattice_index ( linked_list[i], y, z) ]; 
				
				}
				else {

					tmp_par_ptr_1->coords = linked_list[i+1]; 
					tmp_par_ptr_2 = (this->Lattice)[ lattice_index (linked_list[i+1], y, z) ]; 
					(this->Lattice)[ lattice_index (linked_list[i+1], y, z) ] = tmp_par_ptr_1;
					tmp_par_ptr_1 = tmp_par_ptr_2; 
					
				}
			}
		}	
	}

	return;
}

void Simulation::debug_backward_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth){

	int deg_poly = this->Polymers[p_idx].deg_poly;
	if (m_idx == 0) {
		this->IMP_BOOL = true; 
		return; 
	}

	// these are doubles for energies 
	double Em_i       = 0;
	double Es_i       = 0;
	double Em_n       = 0;
	double Es_n       = 0;
	double Epair_i    = 0;
	double Epair_n    = 0;
	double Esys       = *back_energy;
	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <double,5>               energies;
	std::array <double,5>               boltzmann;
	std::array <double,8>               current_contacts = *back_contacts;
	std::array <std::array<double,8>,5> contacts_store;

	std::array <int,3> loc_m = (this->Polymers)[p_idx].chain[m_idx-1]->coords; // this is the key item of interest

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((this->Polymers)[p_idx].chain[m_idx]->coords, this->x, this->y, this->z); 
	std::cout << "ne_list = " << std::endl;
	for (std::array<int,3> ne: ne_list){
		print(ne);
	}

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[m_idx-1]) - ne_list.begin() ; 

	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	std::array <double,8> cm_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int idx_counter         = 0;
	int self_swap_idx       = -1;

	while (idx_counter < 5){
		std::cout << "While index counter is " << idx_counter << "." << std::endl;
		std::cout << "Initial config:" << std::endl;
		this->Polymers[p_idx].print_chain();

		std::cout << "Going into the if-statements..." << std::endl;
		std::cout << "ne_list[" << idx_counter << "] = "; print(ne_list[idx_counter]);
		std::cout << "loc_m = "; print(loc_m);
		if (ne_list[idx_counter] == loc_m){
			std::cout << "is this the segfault locations?" << std::endl;
			energies[idx_counter]       = *back_energy;
			contacts_store[idx_counter] = current_contacts;
		}
		else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm'){
			
			std::cout << "monomer-monomer switch." << std::endl;
			for (int u{0}; u<deg_poly; ++u){
				if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
					self_swap_idx = u;
					break;
				}
			}

			if (self_swap_idx >= m_idx){
				energies[idx_counter]       = 1e+8;
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1};
			}
			else {
				std::cout << "Running the monomer-monomer switch." << std::endl;
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

				// swap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(contacts_store[idx_counter]);

				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(energies[idx_counter], contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// resetting...
				cm_i    = {0,0,0,0,0,0,0,0};
				cs_i    = {0,0,0,0,0,0,0,0};
				cpair_i = {0,0,0,0,0,0,0,0};
				cm_n    = {0,0,0,0,0,0,0,0};
				cs_n    = {0,0,0,0,0,0,0,0};
				cpair_n = {0,0,0,0,0,0,0,0};
				Em_i    = 0;
				Es_i    = 0;
				Em_n    = 0;
				Es_n    = 0;
				Epair_i = 0;
				Epair_n = 0;

			}

		}
		else {
			std::cout << "Running the monomer-solvent switch." << std::endl;
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

			// swap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);
			std::cout << "contacts_store[" << idx_counter << "] = "; print(contacts_store[idx_counter]);

			this->Polymers[p_idx].print_chain();
			this->debug_checks_energy_contacts(energies[idx_counter], contacts_store[idx_counter]);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// resetting...
			cm_i    = {0,0,0,0,0,0,0,0};
			cs_i    = {0,0,0,0,0,0,0,0};
			cpair_i = {0,0,0,0,0,0,0,0};
			cm_n    = {0,0,0,0,0,0,0,0};
			cs_n    = {0,0,0,0,0,0,0,0};
			cpair_n = {0,0,0,0,0,0,0,0};
			Em_i    = 0;
			Es_i    = 0;
			Em_n    = 0;
			Es_n    = 0;
			Epair_i = 0;
			Epair_n = 0;
		}
		idx_counter += 1;
	}

	double Emin = *std::min_element(energies.begin(), energies.end());
	for (int i{0}; i<5; ++i){
		boltzmann[i] = std::exp(-1/this->T * (energies[i]-Emin));
		rboltzmann  += boltzmann[i];
	}

	*prob_n_to_o   = (*prob_n_to_o) * boltzmann[0]/rboltzmann;
	*back_energy   = energies[0];
	*back_contacts = contacts_store[0];

	// perform the swap
	this->perturb_particle_swap(lattice_index(ne_list[0], this->y, this->z), lattice_index(loc_m, this->y, this->z));
	std::cout << "Check at end of regrowth @ recursion depth = " << recursion_depth << "." << std::endl;
	this->Polymers[p_idx].print_chain();
	this->debug_checks_energy_contacts(energies[0], contacts_store[0]);
	this->debug_backward_tail_regrowth(old_cut, back_contacts, prob_n_to_o, back_energy, p_idx, m_idx-1, recursion_depth+1);

	return;
}

