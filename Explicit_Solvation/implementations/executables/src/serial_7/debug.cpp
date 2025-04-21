#include "Simulation.hpp"
#include "lattice_directions.hpp"
#include "misc.hpp"

// all methods will have the prefix `debug_`
// all properties and variable unique to these methods will have the prefix `db_`

void Simulation::debug_calculate_energy(double* db_energy, std::array<double,CONTACT_SIZE>* db_contacts){

	*db_energy   = 0.0;
	*db_contacts = {0,0,0, 0,0,0, 0,0}; 
	std::array <std::array <int,3>, 26> ne_list; 

	// more debug shit
	// int idx = 0;
	// run energy computations for every monomer bead 
	// std::cout << "Running a debugging energy calculation..." << std::endl;
	for (Polymer& pmer: this->Polymers) {
		for (Particle*& p: pmer.chain) {
			ne_list = obtain_ne_list(p->coords, this->x, this->y, this->z);
			for (std::array <int,3>& loc: ne_list) {
				this->debug_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], db_contacts, db_energy);
			}
		}
	}

	for (Particle*& p: this->Cosolvent){
		ne_list = obtain_ne_list (p->coords, this->x, this->y, this->z);
		for (std::array <int,3>& loc: ne_list){
			if (this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "m1" || 
			this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "s2" || p->ptype != "s2"){
				// std::cout << "This better not be anything other than s2: " << p->ptype << std::endl;
				continue;
			}
			this->debug_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], db_contacts, db_energy);
		}
	}
	return;
}

void Simulation::debug_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* db_contacts, double* db_energy){

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
		connvec = subtract_containers ( (p2->coords), (p1->coords) );
		modified_direction (&connvec, x, y, z); 
		magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
		theta_1   = branchless_acos(take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
		theta_2   = branchless_acos(take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );

		// std::cout << "theta_1 = " << theta_1 << ", theta_2 = " << theta_2 << std::endl;
		if ( (theta_1 + theta_2) > M_PI/2 ){
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ) {
				// std::cout << "Misaligned." << std::endl;
				(*db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 0.5; 
				*db_energy += std::get<2>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else if ((particle_pair.first == "s1" && particle_pair.second == "s2") || (particle_pair.first == "s2" && particle_pair.second == "s1")){
				// std::cout << "\t@ p1 = "; print(p1->coords, ", "); std::cout << "p2 = "; print(p2->coords, ", "); std::cout << "Misaligned." << std::endl;
				(*db_contacts)[std::get<4>((this->InteractionMap)[particle_pair])] += 1; 
				*db_energy += std::get<2>((this->InteractionMap)[particle_pair]); 				
			}
			else {
				(*db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 1; 	
				*db_energy += std::get<2>((this->InteractionMap)[particle_pair]);
			}
		}
		else {
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				// std::cout << "Aligned." << std::endl;
				(*db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])] += 0.5; 
				*db_energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else if ((particle_pair.first == "s1" && particle_pair.second == "s2") || (particle_pair.first == "s2" && particle_pair.second == "s1")){
				// std::cout << "\t@ p1 = "; print(p1->coords, ", "); std::cout << "p2 = "; print(p2->coords, ", "); std::cout << "Aligned." << std::endl;
				(*db_contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1; 
				*db_energy += std::get<1>((this->InteractionMap)[particle_pair]); 				
			}
			else {
				(*db_contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1; 
				*db_energy += std::get<1>((this->InteractionMap)[particle_pair]); 
			}
		}
		break;
	}

	return;
}

void Simulation::debug_checks_energy_contacts(double E_final, std::array<double,CONTACT_SIZE> final_contacts){

	double E_debug = 0;
	std::array<double,CONTACT_SIZE> debug_contacts = {0,0,0, 0,0,0, 0,0}; 

	// run the debugger
	this->debug_calculate_energy(&E_debug, &debug_contacts); 
	if ( std::abs(E_debug - E_final) > 1e-8){
		std::cout << "E_debug - E_final = " << E_debug-E_final << std::endl;
		std::cout << "Energy mismatch." << std::endl;
		std::cout << "E_real = " << E_debug << ", E_final = " << E_final << "." << std::endl;
		std::cout << "real_contacts  = "; print(debug_contacts); 
		std::cout << "final_contacts = "; print(final_contacts);
		exit(EXIT_FAILURE);
	} 
	else {
		// std::cout << "energies match @ " << E_debug << "!" << std::endl;
	}

	for (int i{0}; i<CONTACT_SIZE; ++i){
		if (debug_contacts[i] != final_contacts[i]){
			std::cout << "Contacts mismatch." << std::endl;
			std::cout << "real_contacts  = "; print(debug_contacts); 
			std::cout << "final_contacts = "; print(final_contacts);
			exit(EXIT_FAILURE);
		}
		else {
			// std::cout << "contacts match @ "; print(final_contacts, ""); std::cout << "!" << std::endl;
		}
	}

	return; 
}

/////////////////////////////////////////////////////////////////

void Simulation::debug_orientation_sampler_forwards(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iterator_idx, int lat_idx){

	this->enhanced_flipper.initial_orientations[iterator_idx] = this->Lattice[lat_idx]->orientation;
	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.initial_contacts), &(this->enhanced_flipper.initial_E));

	for (int j{0}; j<this->enhanced_flipper.ntest; ++j){
		(this->Lattice)[lat_idx]->orientation = rng_uniform(0, 25);
		this->enhanced_flipper.orientations[j] = this->Lattice[lat_idx]->orientation;
		this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
		this->enhanced_flipper.energies[j] = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E;
		this->enhanced_flipper.contacts_store[j]  = subtract_containers(*contacts_sys, (this->enhanced_flipper.initial_contacts));
		this->enhanced_flipper.contacts_store[j]  = add_containers     ((this->enhanced_flipper.contacts_store[j]), (this->enhanced_flipper.perturbed_contacts));
		this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.perturbed_E        = 0;
		std::cout << "Forward sampler (" << j << ")" << std::endl;
		this->debug_checks_energy_contacts(this->enhanced_flipper.energies[j], this->enhanced_flipper.contacts_store[j]);
	}

	return;

}

void Simulation::debug_orientation_sampler_backwards_0(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx){

	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.initial_contacts), &(this->enhanced_flipper.initial_E));

	this->Lattice[lat_idx]->orientation = this->enhanced_flipper.initial_orientations[iteration_idx];
	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
	this->enhanced_flipper.energies[0] = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E;
	this->enhanced_flipper.contacts_store[0]  = subtract_containers(*contacts_sys, this->enhanced_flipper.initial_contacts);
	this->enhanced_flipper.contacts_store[0]  = add_containers     (this->enhanced_flipper.contacts_store[0], this->enhanced_flipper.perturbed_contacts);
	std::cout << "Inside fourth set of flip test..." << std::endl;
	this->debug_checks_energy_contacts(this->enhanced_flipper.energies[0], this->enhanced_flipper.contacts_store[0]);

	this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
	this->enhanced_flipper.perturbed_E        = 0;

	return;

}

void Simulation::debug_orientation_sampler_backwards(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx){

	this->debug_orientation_sampler_backwards_0(contacts_sys, E_sys, iteration_idx, lat_idx);

	std::cout << "Inside fifth set of flip test." << std::endl;
	for (int j{1}; j<this->enhanced_flipper.ntest; ++j){
		this->Lattice[lat_idx]->orientation = rng_uniform (0, 25); 
		this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
		this->enhanced_flipper.energies[j] = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E;
		this->enhanced_flipper.contacts_store [j]  = subtract_containers (*contacts_sys, this->enhanced_flipper.initial_contacts);
		this->enhanced_flipper.contacts_store [j]  = add_containers      (this->enhanced_flipper.contacts_store[j], this->enhanced_flipper.perturbed_contacts);
		this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.perturbed_E        = 0;
		this->debug_checks_energy_contacts(this->enhanced_flipper.energies[j], this->enhanced_flipper.contacts_store[j]);
	}

	return;

}

void Simulation::debug_choose_state_forward(int iterator_idx, int lat_idx){

	// get the total boltzmann sum
	for ( int k{0}; k < this->enhanced_flipper.ntest; ++k ){
		this->enhanced_flipper.boltzmann[k] = std::exp (-1/this->T*(this->enhanced_flipper.energies[k] - this->enhanced_flipper.Emin));
		this->enhanced_flipper.rboltzmann  += this->enhanced_flipper.boltzmann [k]; 
	}

	// get the rng to sample the right orientation
	this->enhanced_flipper.sampler_rng   = rng_uniform (0.0, 1.0);

	// get the appropriate state
	for ( int j{0}; j<5; ++j){
		this->enhanced_flipper.sampler_rsum += this->enhanced_flipper.boltzmann[j]/this->enhanced_flipper.rboltzmann;
		if (this->enhanced_flipper.sampler_rng < this->enhanced_flipper.sampler_rsum){
			this->enhanced_flipper.sampler_idx = j; 
			break; 
		}
	}

	// make the jump to the new state 
	std::cout << "orientation to pick = "; print(this->enhanced_flipper.orientations);
	this->enhanced_flipper.final_orientations[iterator_idx]                       = this->enhanced_flipper.orientations[this->enhanced_flipper.sampler_idx];
	this->Lattice[lat_idx]->orientation = this->enhanced_flipper.orientations[this->enhanced_flipper.sampler_idx];
	this->enhanced_flipper.prob_o_to_n *= this->enhanced_flipper.boltzmann[this->enhanced_flipper.sampler_idx]/this->enhanced_flipper.rboltzmann;

	return;
}

//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////

void Simulation::debug_tail_rotation(int p_idx){

	// get the neighborlist of the particle at index 1
	std::array<int,3> loc_0 = this->Polymers[p_idx].chain[0]->coords;
	std::array<int,3> loc_1 = this->Polymers[p_idx].chain[1]->coords;

	// define some container for energy and contacts
	double                final_E        = 0; // energy of the final configuration of the system
	std::array <double,CONTACT_SIZE> final_contacts = {0,0,0, 0,0,0, 0,0};

	// generate a neighbor list 
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

	// get the choice for a neighbor
	int choice = rng_uniform(0, 25);

	// instantiate the pointer
	Particle* tmp_par_ptr {nullptr};

	// instantiate the switch indices
	int lat_idx_1 = lattice_index(ne_list[choice], this->y, this->z);
	int lat_idx_2 = lattice_index(loc_0,           this->y, this->z);

	// reset the container
	this->rotation_container.reset();

	if ( this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->ptype[0] == 's' ){

		this->perturb_particle_swap_with_update(tmp_par_ptr, lat_idx_1, lat_idx_2);
		this->check_structures();
		std::cout << "Post-swing structure:" << std::endl;
		this->Polymers[p_idx].print_chain();

		final_E = this->sysEnergy - (this->rotation_container.E1_initial + this->rotation_container.E2_initial - this->rotation_container.Epair_initial);
		final_E = final_E + (this->rotation_container.E1_final   + this->rotation_container.E2_final  - this->rotation_container.Epair_final);
		final_contacts = add_containers(subtract_containers(this->contacts, add_containers (this->rotation_container.c1_initial, this->rotation_container.c2_initial)), add_containers (this->rotation_container.c1_final, this->rotation_container.c2_final));
		final_contacts = add_containers(final_contacts, this->rotation_container.cpair_initial);
		final_contacts = subtract_containers(final_contacts, this->rotation_container.cpair_final);

		// run the debugger
		std::cout << "Running the debugger..." << std::endl; 
		this->debug_checks_energy_contacts(final_E, final_contacts);
		std::cout << "Ran the debugger!" << std::endl;

		if (metropolis_acceptance(this->sysEnergy, final_E, this->T)){
			std::cout << "Accepted!" << std::endl;
			this->sysEnergy = final_E;
			this->contacts  = final_contacts;
			this->IMP_BOOL  = true;
		}
		else {
			// revert the polymer 
			std::cout << "Reverting..." << std::endl;
			this->perturb_particle_swap(tmp_par_ptr, lat_idx_1, lat_idx_2);
			this->IMP_BOOL = false;
			this->check_structures();
			this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
		}
		std::cout << "Final polymer config:" << std::endl;
		this->Polymers[p_idx].print_chain();
	}
	else {
		std::cout << "No rotation." << std::endl;
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

	// define some container for energy and contacts
	double                final_E        = 0; // energy of the final configuration of the system
	std::array <double,CONTACT_SIZE> final_contacts = {0,0,0, 0,0,0, 0,0};

	// generate a neighbor list 
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

	// get the choice for a neighbor
	int choice = rng_uniform(0, 25);

	// instantiate the pointer
	Particle* tmp_par_ptr {nullptr};

	// instantiate the switch indices
	int lat_idx_1 = lattice_index(ne_list[choice], this->y, this->z);
	int lat_idx_2 = lattice_index(loc_0,           this->y, this->z);

	// reset the container containers
	this->rotation_container.reset();

	if ( this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->ptype[0] == 's' ){

		this->perturb_particle_swap_with_update(tmp_par_ptr, lat_idx_1, lat_idx_2);
		this->check_structures();
		this->Polymers[p_idx].print_chain();
		
		// doing the quick manipulations to get the final energy and contacts
		final_E = this->sysEnergy - (this->rotation_container.E1_initial + this->rotation_container.E2_initial - this->rotation_container.Epair_initial);
		final_E = final_E + (this->rotation_container.E1_final   + this->rotation_container.E2_final  - this->rotation_container.Epair_final);
		final_contacts = add_containers(subtract_containers(this->contacts, add_containers (this->rotation_container.c1_initial, this->rotation_container.c2_initial)), add_containers (this->rotation_container.c1_final, this->rotation_container.c2_final));
		final_contacts = add_containers(final_contacts, this->rotation_container.cpair_initial);
		final_contacts = subtract_containers(final_contacts, this->rotation_container.cpair_final);

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
			this->perturb_particle_swap(tmp_par_ptr, lat_idx_1, lat_idx_2);
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

	int                             deg_poly   = this->Polymers[p_idx].deg_poly;
	std::array<int,3>               loc_i      = this->Polymers[p_idx].chain[0]->coords;
	std::array<int,3>               loc_f      = this->Polymers[p_idx].chain[deg_poly-1]->coords;
	std::array<double,CONTACT_SIZE> contacts_i = this->contacts; 
	double                          E_i        = this->sysEnergy;

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

	int                              deg_poly         = this->Polymers[p_idx].deg_poly;
	std::array <int,3>               loc_i            = this->Polymers[p_idx].chain[0]->coords;
	std::array <int,3>               loc_f            = this->Polymers[p_idx].chain[deg_poly-1]->coords;
	std::array <double,CONTACT_SIZE> contacts_i       = this->contacts;
	double                           E_i              = this->sysEnergy;

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

//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////

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

	// set up the flipper object! 
	this->enhanced_flipper.reset(nflip, ntest); 

	// define some holders 
	double                          E_sys        = this->sysEnergy;
	std::array<double,CONTACT_SIZE> contacts_sys = this->contacts;

	// some more holders
	double                          E_post        = 0;                 // the energy of the system after final perturbation
	double                          rng_acc       = 0;                 // rng for acceptance at the very end
	std::array<double,CONTACT_SIZE> contacts_post;
	contacts_post.fill(0); // the contacts of the system after final perturbation

	// relevant indices variable store
	int m_lattice_idx = 0; 

	// loop over the different monomer indices
	for ( int i{0}; i < nflip; ++i ){

		// sample different orientations and get the boltzmann factors 
		m_lattice_idx = lattice_index(this->Polymers[p_idx].chain[polymer_indices[i]]->coords, this->y, this->z);
		this->debug_orientation_sampler_forwards(&contacts_sys, E_sys, i, m_lattice_idx);

		// get the energy minima
		this->enhanced_flipper.Emin = *std::min_element(this->enhanced_flipper.energies.begin(), this->enhanced_flipper.energies.end());

		// go to the new state 
		this->debug_choose_state_forward(i, m_lattice_idx);

		// set the system energy at the final point
		E_sys        = this->enhanced_flipper.energies[enhanced_flipper.sampler_idx];
		contacts_sys = this->enhanced_flipper.contacts_store[enhanced_flipper.sampler_idx];

		// run the check
		std::cout << "Second flip test..." << std::endl;
		this->debug_checks_energy_contacts(E_sys, contacts_sys);

		// reset
		this->enhanced_flipper.initial_contacts.fill(0);
		this->enhanced_flipper.initial_E        = 0;
		this->enhanced_flipper.rboltzmann       = 0;
		this->enhanced_flipper.sampler_rsum     = 0;

	}

	E_post        = E_sys;
	contacts_post = contacts_sys;

	this->debug_checks_energy_contacts(E_sys, contacts_sys);
	this->enhanced_flipper.perturbed_contacts.fill(0);
	this->enhanced_flipper.perturbed_E        = 0;


	for ( int i{0}; i < nflip; ++i ){

		// sample different orientations and get the boltzmann factors for the backward flux
		m_lattice_idx = lattice_index(this->Polymers[p_idx].chain[polymer_indices[i]]->coords, this->y, this->z);
		this->debug_orientation_sampler_backwards(&contacts_sys, E_sys, i, m_lattice_idx);

		// get the minimum energy
		this->enhanced_flipper.Emin = *std::min_element(this->enhanced_flipper.energies.begin(), this->enhanced_flipper.energies.end());

		// get the backwards probability flux
		for (int k{0}; k < ntest; ++k){
			this->enhanced_flipper.boltzmann [k] = std::exp (-1/this->T*(this->enhanced_flipper.energies[k] - this->enhanced_flipper.Emin));
			this->enhanced_flipper.rboltzmann   += this->enhanced_flipper.boltzmann[k];
		}
		this->enhanced_flipper.prob_n_to_o      *= this->enhanced_flipper.boltzmann[0]/this->enhanced_flipper.rboltzmann;

		// make the jump to the old state 
		this->Polymers[p_idx].chain[polymer_indices[i]]->orientation = this->enhanced_flipper.initial_orientations[i];

		// set the energetics and contacts right
		E_sys        = this->enhanced_flipper.energies[0];
		contacts_sys = this->enhanced_flipper.contacts_store[0];
		std::cout << "Sixth flip test..." << std::endl;
		this->debug_checks_energy_contacts(E_sys, contacts_sys);

		// reset
		this->enhanced_flipper.initial_contacts.fill(0);
		this->enhanced_flipper.initial_E        = 0;
		this->enhanced_flipper.rboltzmann       = 0;

	}

	std::cout << "Seventh flip test..." << std::endl;
	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);

	// get a random number
	rng_acc = rng_uniform(0.0, 1.0);

	// run the criterion
	if ( rng_acc < std::exp (-1/this->T * (E_post - this->sysEnergy )) * this->enhanced_flipper.prob_n_to_o/this->enhanced_flipper.prob_o_to_n){
		std::cout << "Orientation to be assigned = "; print(this->enhanced_flipper.final_orientations);
		for (int j{0}; j < nflip; ++j){
			(this->Polymers)[p_idx].chain[polymer_indices[j]]->orientation = this->enhanced_flipper.final_orientations[j];
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

	// set up the flipper object
	this->enhanced_flipper.reset(nflip, ntest);

	// define some holders 
	double               E_sys                   = this->sysEnergy;
	std::array<double,CONTACT_SIZE> contacts_sys = this->contacts;

	double               E_post                   = 0; // the energy of the system after final perturbation
	std::array<double,CONTACT_SIZE> contacts_post;
	contacts_post.fill(0); // the contacts of the system after final perturbation
	double               rng_acc                  = 0; // rng for acceptance at the very end

	for (int i{0}; i < nflip; ++i){

		// sample different orientations and get the boltzmann factors 
		this->debug_orientation_sampler_forwards(&contacts_sys, E_sys, i, solvation_shell_indices[i]);

		// get the energy minima
		this->enhanced_flipper.Emin = *std::min_element (this->enhanced_flipper.energies.begin(), this->enhanced_flipper.energies.end() );

		// go to the new state
		this->debug_choose_state_forward(i, solvation_shell_indices[i]);

		// set the system energy at the final point
		E_sys         = this->enhanced_flipper.energies[enhanced_flipper.sampler_idx];
		contacts_sys  = this->enhanced_flipper.contacts_store[this->enhanced_flipper.sampler_idx];

		// run the check
		this->debug_checks_energy_contacts(E_sys, contacts_sys);

		// reset
		this->enhanced_flipper.initial_contacts.fill(0);
		this->enhanced_flipper.initial_E        = 0;
		this->enhanced_flipper.rboltzmann       = 0;
		this->enhanced_flipper.sampler_rsum     = 0;

	}

	E_post = E_sys;
	contacts_post = contacts_sys;

	this->debug_checks_energy_contacts(E_sys, contacts_sys);
	this->enhanced_flipper.perturbed_contacts.fill(0);
	this->enhanced_flipper.perturbed_E        = 0;


	for (int i{0}; i<nflip; ++i){

		// sample different orientations and get the boltzmann factors for the backward flux
		this->debug_orientation_sampler_backwards(&contacts_sys, E_sys, i, solvation_shell_indices[i]);

		// get the minimum energy
		this->enhanced_flipper.Emin = *std::min_element(this->enhanced_flipper.energies.begin(), this->enhanced_flipper.energies.end());
		
		// get the backwards probability flux
		for (int k{0}; k < ntest; ++k){
			this->enhanced_flipper.boltzmann [k] = std::exp (-1/this->T*(this->enhanced_flipper.energies[k] - this->enhanced_flipper.Emin));
			this->enhanced_flipper.rboltzmann   += this->enhanced_flipper.boltzmann[k];
		}
		this->enhanced_flipper.prob_n_to_o      *= this->enhanced_flipper.boltzmann[0]/this->enhanced_flipper.rboltzmann;

		// make the jump to the old state 
		this->Lattice[solvation_shell_indices[i]]->orientation = this->enhanced_flipper.initial_orientations[i];

		// set the energetics and contacts right
		E_sys        = this->enhanced_flipper.energies[0];
		contacts_sys = this->enhanced_flipper.contacts_store[0];
		std::cout << "Sixth flip test..." << std::endl;
		this->debug_checks_energy_contacts(E_sys, contacts_sys);

		// reset
		this->enhanced_flipper.initial_contacts.fill(0);
		this->enhanced_flipper.initial_E        = 0;
		this->enhanced_flipper.rboltzmann       = 0;

	}

	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
	if (this->sysEnergy != E_sys || this->contacts != contacts_sys){
		std::cerr << "Something's off." << std::endl;
		std::cout << "this->sysEnergy = " << this->sysEnergy << ", E_sys = " << E_sys << "." << std::endl;
		std::cout << "this->contacts  = "; print(this->contacts, ", ");
		std::cout << "contacts_sys = "; print(contacts_sys);
		exit(EXIT_FAILURE);
	}

	rng_acc = rng_uniform(0.0, 1.0);
	if ( rng_acc < std::exp (-1/this->T * (E_post - this->sysEnergy)) * this->enhanced_flipper.prob_n_to_o/this->enhanced_flipper.prob_o_to_n){

		// if accepted, return to the new orientations 
		for (int j{0}; j < nflip; ++j){
			this->Lattice[solvation_shell_indices[j]]->orientation = this->enhanced_flipper.final_orientations[j]; 
		}

		this->sysEnergy = E_post;
		this->contacts  = contacts_post;
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

	// number of particles to flip
	int nflips = rng_uniform(1, static_cast<int>(this->x * this->y * this->z / 8));
	std::array <double,CONTACT_SIZE> contacts = this->contacts; 

	// set up the orientations
	std::vector <int> orientations (nflips, 0); 

	// set up the contacts store
	std::array <double,CONTACT_SIZE> cs_i = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cs_f = {0,0,0, 0,0,0, 0,0};

	// set up the energies
	double Es_i = 0;
	double Es_f = 0;
	double Ef   = this->sysEnergy; 

	// set up the samples
	std::vector<int> samples(this->x * this->y * this->z - 1); 
	std::iota(samples.begin(), samples.end(), 0);

	// randomize
	std::random_device rd;
	std::mt19937 g(rd()); 
	std::shuffle(samples.begin(), samples.end(), g); 

	std::cout << "orientations.size() = " << orientations.size() << std::endl;

	for (int i{0}; i<nflips; ++i){

		this->neighbor_energetics(samples[i], &cs_i, &Es_i); 
		orientations[i] =  (this->Lattice[samples[i]])->orientation; // .push_back((this->Lattice[samples[i]])->orientation); 
		(this->Lattice[samples[i]])->orientation = rng_uniform(0, 25);
		this->neighbor_energetics(samples[i], &cs_f, &Es_f); 
		contacts = add_containers(subtract_containers(contacts, cs_i), cs_f);
		Ef += Es_f - Es_i;
		std::cout << "i = " << i << ", " << this->Lattice[samples[i]]->ptype << ", "; print(this->Lattice[samples[i]]->coords);
		this->debug_checks_energy_contacts(Ef, contacts);

		// resetting...
		cs_i = {0,0,0, 0,0,0, 0,0};
		cs_f = {0,0,0, 0,0,0, 0,0};
		Es_i = 0;
		Es_f = 0;

	}

	double rng = rng_uniform(0.0, 1.0);

	if ( rng < std::exp(-1/this->T * (Ef - this->sysEnergy)) ) {
		this->sysEnergy = Ef;
		this->contacts  = contacts;
		this->IMP_BOOL  = true;
	}
	else {
		for (int i{0}; i<nflips; ++i){
			(this->Lattice[samples[i]])->orientation = orientations[i];
		}
		this->IMP_BOOL = false;
	}

	return;
}

// this method entails a series of exchanges
void Simulation::debug_solvent_exchange_from_shell(){

	std::set    <int> solvation_shell_set = this->get_solvation_shell();
	std::vector <int> solvation_shell_indices(solvation_shell_set.begin(), solvation_shell_set.end());

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));

	// reset the container
	this->rotation_container.reset();

	// set up some holders
	double post_energy = 0;
	double rng_acc     = 0;
	std::array <double,CONTACT_SIZE> post_contacts = {0,0,0, 0,0,0, 0,0};

	// define the holder to facilitate the swap
	Particle* tmp_par_ptr {nullptr};

	int exc_idx = 0; 
	int nswitch = rng_uniform(1, static_cast<int>(solvation_shell_indices.size()/4));

	for (int n{0}; n<nswitch; ++n){
		
		exc_idx = rng_uniform(0, this->x * this->y * this->z - 1);
		if (exc_idx == solvation_shell_indices[n] || this->Lattice[exc_idx]->ptype[0] == 'm'){
			continue;
		}
		
		else {

			this->perturb_particle_swap_with_update(tmp_par_ptr, solvation_shell_indices[n], exc_idx);

			// update energy and contacts
			post_energy   = this->sysEnergy  - (this->rotation_container.E1_initial + this->rotation_container.E2_initial - this->rotation_container.Epair_initial); 
			post_energy   = post_energy + (this->rotation_container.E1_final + this->rotation_container.E2_final  - this->rotation_container.Epair_final);
			post_contacts = add_containers(subtract_containers(this->contacts, add_containers (this->rotation_container.c1_initial, this->rotation_container.c2_initial)), add_containers (this->rotation_container.c1_final, this->rotation_container.c2_final));
			post_contacts = add_containers(post_contacts, this->rotation_container.cpair_initial);
			post_contacts = subtract_containers(post_contacts, this->rotation_container.cpair_final);
			
			this->debug_checks_energy_contacts(post_energy, post_contacts);

			rng_acc = rng_uniform (0.0, 1.0);

			if (rng_acc < std::exp (-1/this->T * (post_energy - this->sysEnergy))){
				this->sysEnergy = post_energy;
				this->contacts  = post_contacts;
				this->IMP_BOOL  = true; 
			}
			else {
				this->perturb_particle_swap(tmp_par_ptr, solvation_shell_indices[n], exc_idx);
				this->IMP_BOOL = false;
			}

			// run the debugger
			this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);

			// resetting 
			this->rotation_container.reset();
		}
	}
	return;
	
}

void Simulation::debug_solvent_exchange(){

	// these are the indices on the lattice for the particles to swap
	int idx1 = -1;
	int idx2 = -1; 

	// this is the number of particles to switch up
	int nswitches = 50; 

	// reset the container
	this->rotation_container.reset();

	// these are the variables 
	double post_energy = 0;
	double rng_acc          = 0;
	std::array <double,CONTACT_SIZE> post_contacts = this->contacts;

	// set up the pointer to facilitate the swap
	Particle* tmp_par_ptr {nullptr};

	for (int j{0}; j<nswitches; ++j){

		idx1 = rng_uniform(0, this->x * this->y * this->z - 1); 
		idx2 = rng_uniform(0, this->x * this->y * this->z - 1); 

		if (idx1 == idx2 || this->Lattice[idx1]->ptype == "m1" || this->Lattice[idx2]->ptype == "m1"){
			continue; 
		}
		else {
			this->perturb_particle_swap_with_update(tmp_par_ptr, idx1, idx2);
			
			post_energy   = this->sysEnergy  - (this->rotation_container.E1_initial + this->rotation_container.E2_initial - this->rotation_container.Epair_initial); 
			post_energy   = post_energy + (this->rotation_container.E1_final   + this->rotation_container.E2_final   - this->rotation_container.Epair_final);
			post_contacts = add_containers(subtract_containers(this->contacts, add_containers (this->rotation_container.c1_initial, this->rotation_container.c2_initial)), add_containers (this->rotation_container.c1_final, this->rotation_container.c2_final));
			post_contacts = add_containers(post_contacts, this->rotation_container.cpair_initial);
			post_contacts = subtract_containers(post_contacts, this->rotation_container.cpair_final);

			this->debug_checks_energy_contacts(post_energy, post_contacts);

			rng_acc = rng_uniform (0.0, 1.0); 
			if (rng_acc < std::exp (-1/this->T * (post_energy - this->sysEnergy))){
				this->sysEnergy = post_energy;
				this->contacts  = post_contacts;
		
			}
			else {
				this->perturb_particle_swap(tmp_par_ptr, idx1, idx2);
			}
			this->rotation_container.reset();
		}

	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////

void Simulation::debug_regrowth(int p_idx){

	// set up some containers 
	int deg_poly        = this->Polymers[p_idx].deg_poly; 
	int m_idx           =  deg_poly-4;
	int growth          = -1;
	double forw_energy  =  0;
	double rng_acc      =  0;
	std::array <double,CONTACT_SIZE> forw_contacts = {0,0,0, 0,0,0, 0,0};

	// get that reset going
	this->enhanced_swing.reset(deg_poly);
	this->enhanced_swing.current_energy   = this->sysEnergy;
	this->enhanced_swing.current_contacts = this->contacts;

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
			this->enhanced_swing.initial_locations.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		std::cout << "\tPerforming the forward pass..." << std::endl;
		this->debug_forward_head_regrowth(p_idx, m_idx);

		for (int i{m_idx+1}; i<deg_poly; ++i) {
			this->enhanced_swing.final_locations.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		std::cout << "initial size = " << this->enhanced_swing.initial_locations.size() << std::endl;
		std::cout << "final size   = " << this->enhanced_swing.initial_locations.size() << std::endl;

		if (this->enhanced_swing.initial_locations == this->enhanced_swing.final_locations){
			// do nothing
			;
		}

		else if (!(this->IMP_BOOL)){
			this->debug_accept_after_head_regrowth(this->IMP_BOOL);
			this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
			this->check_structures();
		}

		else {
			forw_energy   = this->enhanced_swing.current_energy;
			forw_contacts = this->enhanced_swing.current_contacts;
			std::cout << "\tPerforming the backward pass..." << std::endl;
			this->debug_backward_head_regrowth(p_idx, m_idx, 0);

			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * this->enhanced_swing.prob_n_to_o / this->enhanced_swing.prob_o_to_n){
				this->debug_accept_after_head_regrowth(this->IMP_BOOL);
				this->sysEnergy = forw_energy; 
				this->contacts  = forw_contacts;
			}
			else {
				this->IMP_BOOL = false;
			}
		}
		this->Polymers[0].print_chain();
		this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
		this->check_structures();
	}

	else {

		std::cout << "\tPerforming tail regrowth..." << std::endl;
		for (int i{0}; i<m_idx; ++i){
			// old_cut.push_back ((this->Polymers).at(p_idx).chain.at(i)->coords) ;
			this->enhanced_swing.initial_locations.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		this->debug_forward_tail_regrowth(p_idx, m_idx);

		for (int i{0}; i<m_idx; ++i){
			// new_cut.push_back ((this->Polymers).at(p_idx).chain.at(i)->coords) ;
			this->enhanced_swing.final_locations.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		if (!(this->IMP_BOOL)){
			std::cout << "Rejecting due to trapping." << std::endl;
			this->debug_accept_after_tail_regrowth(this->IMP_BOOL);
			this->Polymers[0].print_chain();
		}

		else if (this->enhanced_swing.initial_locations == this->enhanced_swing.final_locations){
			// do nothing
			;
		}

		else {

			this->debug_checks_energy_contacts(this->enhanced_swing.current_energy, this->enhanced_swing.current_contacts);
			this->check_structures();

			forw_energy   = this->enhanced_swing.current_energy;
			forw_contacts = this->enhanced_swing.current_contacts;
			std::cout << "\tPerforming backward tail regrowth..." << std::endl;
			this->debug_backward_tail_regrowth(p_idx, m_idx, 0);

			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * this->enhanced_swing.prob_n_to_o/this->enhanced_swing.prob_o_to_n) {
				this->debug_accept_after_tail_regrowth(this->IMP_BOOL);
				this->sysEnergy = forw_energy;
				this->contacts  = forw_contacts;
			}

			else {
				this->IMP_BOOL = false;
			}
			
		}
		this->Polymers[0].print_chain();
		this->debug_checks_energy_contacts(this->sysEnergy, this->contacts);
		this->check_structures();
	}
	

	return;

}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::debug_forward_head_regrowth(int p_idx, int m_idx){

	std::cout << "@m_idx = " << m_idx << std::endl;
	if (m_idx == this->Polymers[p_idx].deg_poly-1){
		this->IMP_BOOL = true;
	}

	else { 

		// reset the running boltzmann sum
		this->enhanced_swing.rboltzmann = 0;

		// reset locally
		this->enhanced_swing.reset_local();

		// set up some contact counters
		std::array <int,3>                loc_m   = this->Polymers[p_idx].chain[m_idx+1]->coords;
		std::array <std::array<int,3>,26> ne_list = obtain_ne_list(this->Polymers[p_idx].chain[m_idx]->coords, this->x, this->y, this->z);
		
		// set up the shuffle
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
		std::shuffle(ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

		// set up some variables
		int block_counter =  0;
		int idx_counter   =  0;
		int self_swap_idx = -1;

		while(idx_counter < this->enhanced_swing.ntest){

			if (ne_list[idx_counter] == loc_m){
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy;
				this->enhanced_swing.contacts_store[idx_counter] = this->enhanced_swing.current_contacts;
				// block_counter                                   += 1;
			}
			else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm') {
				for (int u{0}; u < this->Polymers[p_idx].deg_poly; ++u){
					if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
						self_swap_idx = u; 
						break;
					}
				}

				std::cout << "self swap monomer index = " << self_swap_idx << "." << std::endl;
				if (self_swap_idx < m_idx){
					this->enhanced_swing.energies[idx_counter]       = 1e+8;
					this->enhanced_swing.contacts_store[idx_counter] = {-1,-1,-1, -1,-1,-1, -1,-1};
					block_counter                                   += 1;
				}
				else{

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

					// swap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

					// set up the final energies
					this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
					this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
					this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
					std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

					std::cout << "Running check..." << std::endl;
					this->Polymers[p_idx].print_chain();
					this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

					// revert back to original structure, reswap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// reset some local variables
					this->enhanced_swing.reset_local();

				}
			}
			else {

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

				// swap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

				// set up the final energies
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
				this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
				this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

				std::cout << "Running check..." << std::endl;
				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// reset some local variables
				this->enhanced_swing.reset_local();

			}
			idx_counter += 1;
		}

		if (block_counter == this->enhanced_swing.ntest){
			this->IMP_BOOL = false;
		}
		else {
			this->enhanced_swing.Emin = *std::min_element(this->enhanced_swing.energies.begin(), this->enhanced_swing.energies.end());
			for (int i{0}; i<this->enhanced_swing.ntest; ++i){
				this->enhanced_swing.boltzmann[i] = std::exp(-1/this->T * (this->enhanced_swing.energies[i]-this->enhanced_swing.Emin));
				this->enhanced_swing.rboltzmann  += this->enhanced_swing.boltzmann[i];
			}
			this->enhanced_swing.sampler_rsum = 0;
			this->enhanced_swing.sampler_rng  = rng_uniform(0.0, 1.0);
			for (int j{0}; j<this->enhanced_swing.ntest; ++j){
				this->enhanced_swing.sampler_rsum += this->enhanced_swing.boltzmann[j]/this->enhanced_swing.rboltzmann;
				if (this->enhanced_swing.sampler_rng < this->enhanced_swing.sampler_rsum){
					this->enhanced_swing.sampler_idx = j;
					break;
				}
			}
			this->enhanced_swing.prob_o_to_n     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
			this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
			this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

			// do the swap to get to the new configuration
			this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
			this->Polymers[p_idx].print_chain();
			this->debug_checks_energy_contacts(this->enhanced_swing.energies[this->enhanced_swing.sampler_idx], this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx]);
			this->debug_forward_head_regrowth(p_idx, m_idx+1);

		}

	}

	return;
}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::debug_accept_after_head_regrowth(bool not_trap_bool){

	// define some stuff
	int L{0}; 

	// get the linked list going
	std::vector<std::vector<std::array<int,3>>> master_linked_list;
	std::vector<std::array<int,3>> link;

	// get the linked list 
	if (not_trap_bool){
		create_linked_list (this->enhanced_swing.initial_locations, this->enhanced_swing.final_locations, link, &master_linked_list, 1);
	}
	else {
		create_linked_list (this->enhanced_swing.final_locations,   this->enhanced_swing.initial_locations, link, &master_linked_list, 1);
	}

	// set up some particle pointers 
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

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::debug_backward_head_regrowth(int p_idx, int m_idx, int recursion_depth){

	std::cout << "Check for completion..." << std::endl;
	if (m_idx == this->Polymers[p_idx].deg_poly - 1){
		std::cout << "Completed!" << std::endl;
		this->IMP_BOOL = true;
	}
	else{
		std::cout << "Running the move..." << std::endl;
		// reset the running boltzmann sum
		this->enhanced_swing.rboltzmann = 0;

		// reset locally
		this->enhanced_swing.reset_local();

		// set up some contact counters
		std::array <int,3>                loc_m   = this->Polymers[p_idx].chain[m_idx+1]->coords;
		std::array <std::array<int,3>,26> ne_list = obtain_ne_list(this->Polymers[p_idx].chain[m_idx]->coords, this->x, this->y, this->z);
		
		// set up the shuffle
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
		std::shuffle(ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

		// get the old position
		int ne_idx = std::find (ne_list.begin(), ne_list.end(), this->enhanced_swing.initial_locations[recursion_depth]) - ne_list.begin();
		std::array <int,3> tmp = ne_list[0]; 
		ne_list[0]             = ne_list[ne_idx];
		ne_list[ne_idx]        = tmp; 

		// set up some variables
		int idx_counter   =  0;
		int self_swap_idx = -1;

		std::cout << "Outside the while loop..." << std::endl;
		std::cout << "ntest = " << this->enhanced_swing.ntest << std::endl;
		while(idx_counter < this->enhanced_swing.ntest){
			std::cout << "Inside the while loop..." << std::endl;
			std::cout << "ne_list["<< idx_counter << "] = "; print(ne_list[idx_counter]);
			std::cout << "loc_m = "; print(loc_m);
			if (ne_list[idx_counter] == loc_m){
				std::cout << "Hit the same monomer" << std::endl;
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy;
				this->enhanced_swing.contacts_store[idx_counter] = this->enhanced_swing.current_contacts;
			}
			else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm') {
				std::cout << "Switching with monomer..." << std::endl;
				for (int u{0}; u < this->Polymers[p_idx].deg_poly; ++u){
					if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
						self_swap_idx = u; 
						break;
					}
				}

				std::cout << "self swap monomer index = " << self_swap_idx << "." << std::endl;
				if (self_swap_idx <= m_idx){
					this->enhanced_swing.energies[idx_counter]       = 1e+8;
					this->enhanced_swing.contacts_store[idx_counter] = {-1,-1, -1,-1,-1, -1,-1};
				}
				else{

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

					// swap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

					// set up the final energies
					this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
					this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
					this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
					std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

					std::cout << "Checking the structures..." << std::endl;
					this->Polymers[p_idx].print_chain();
					this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

					// revert back to original structure, reswap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// reset some local variables
					this->enhanced_swing.reset_local();

				}
			}
			else {
				std::cout << "Switching with solvent." << std::endl;

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

				std::cout << "Doing the swap." << std::endl;
				// swap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				std::cout << "Get the energies (post)." << std::endl;
				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

				std::cout << "Store the energies." << std::endl;
				// set up the final energies
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
				this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
				this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

				std::cout << "Checking the structures..." << std::endl;
				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// reset some local variables
				this->enhanced_swing.reset_local();

			}
			idx_counter += 1;
		}
		this->enhanced_swing.Emin = *std::min_element(this->enhanced_swing.energies.begin(), this->enhanced_swing.energies.end());
		for (int i{0}; i<this->enhanced_swing.ntest; ++i){
			this->enhanced_swing.boltzmann[i] = std::exp(-1/this->T * (this->enhanced_swing.energies[i]-this->enhanced_swing.Emin));
			this->enhanced_swing.rboltzmann += this->enhanced_swing.boltzmann[i];
		}

		this->enhanced_swing.sampler_idx      = 0;
		this->enhanced_swing.prob_o_to_n     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
		this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
		this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

		// do the swap to get to the new configuration
		this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->Polymers[p_idx].print_chain();
		this->debug_checks_energy_contacts(this->enhanced_swing.energies[this->enhanced_swing.sampler_idx], this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx]);
		this->debug_backward_head_regrowth(p_idx, m_idx+1, recursion_depth+1);

	}
	return;
}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::debug_forward_tail_regrowth(int p_idx, int m_idx){
	std::cout << "@m_idx = " << m_idx << std::endl;
	if (m_idx == 0){
		this->IMP_BOOL = true;
	}

	else { 

		// reset the running boltzmann sum
		this->enhanced_swing.rboltzmann = 0;

		// reset locally
		this->enhanced_swing.reset_local();

		// set up some contact counters
		std::array <int,3>                loc_m   = this->Polymers[p_idx].chain[m_idx-1]->coords;
		std::array <std::array<int,3>,26> ne_list = obtain_ne_list(this->Polymers[p_idx].chain[m_idx]->coords, this->x, this->y, this->z);
		
		// set up the shuffle
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
		std::shuffle(ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

		// set up some variables
		int block_counter =  0;
		int idx_counter   =  0;
		int self_swap_idx = -1;

		while(idx_counter < this->enhanced_swing.ntest){

			if (ne_list[idx_counter] == loc_m) {
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy;
				this->enhanced_swing.contacts_store[idx_counter] = this->enhanced_swing.current_contacts;
				// block_counter                                   += 1;
			}
			else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm') {
				for (int u{0}; u < this->Polymers[p_idx].deg_poly; ++u){
					if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
						self_swap_idx = u; 
						break;
					}
				}

				std::cout << "self swap monomer index = " << self_swap_idx << "." << std::endl;
				if (self_swap_idx > m_idx){
					this->enhanced_swing.energies[idx_counter]       = 1e+8;
					this->enhanced_swing.contacts_store[idx_counter] = {-1,-1,-1, -1,-1,-1, -1,-1};
					block_counter                                   += 1;
				}
				else{

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

					// swap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

					// set up the final energies
					this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
					this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
					this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
					std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

					this->Polymers[p_idx].print_chain();
					this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

					// revert back to original structure, reswap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// reset some local variables
					this->enhanced_swing.reset_local();

				}
			}
			else {

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

				// swap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

				// set up the final energies
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
				this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
				this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// reset some local variables
				this->enhanced_swing.reset_local();

			}
			idx_counter += 1;
		}

		if (block_counter == this->enhanced_swing.ntest){
			std::cout << "block_counter = " << block_counter << ", ntest = " << this->enhanced_swing.ntest << std::endl;
			this->IMP_BOOL = false;
		}
		
		else {
			this->enhanced_swing.Emin = *std::min_element(this->enhanced_swing.energies.begin(), this->enhanced_swing.energies.end());
			for (int i{0}; i<this->enhanced_swing.ntest; ++i){
				this->enhanced_swing.boltzmann[i] = std::exp(-1/this->T * (this->enhanced_swing.energies[i]-this->enhanced_swing.Emin));
				this->enhanced_swing.rboltzmann += this->enhanced_swing.boltzmann[i];
			}
			this->enhanced_swing.sampler_rsum = 0;
			this->enhanced_swing.sampler_rng  = rng_uniform(0.0, 1.0);
			for (int j{0}; j<this->enhanced_swing.ntest; ++j){
				this->enhanced_swing.sampler_rsum += this->enhanced_swing.boltzmann[j]/this->enhanced_swing.rboltzmann;
				if (this->enhanced_swing.sampler_rng < this->enhanced_swing.sampler_rsum){
					this->enhanced_swing.sampler_idx = j;
					break;
				}
			}
			this->enhanced_swing.prob_o_to_n     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
			this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
			this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

			// do the swap to get to the new configuration
			this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
			this->Polymers[p_idx].print_chain();
			this->debug_checks_energy_contacts(this->enhanced_swing.energies[this->enhanced_swing.sampler_idx], this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx]);
			this->debug_forward_tail_regrowth(p_idx, m_idx-1);

		}

	}

	return;
}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::debug_accept_after_tail_regrowth(bool not_trap_bool){

	// set up some vectors
	std::vector< std::vector<std::array<int,3>>> master_linked_list;
	std::vector <std::array<int,3>> link; 

	if (not_trap_bool){
		create_linked_list (this->enhanced_swing.initial_locations, this->enhanced_swing.final_locations, link, &master_linked_list, 1);
	}
	else {
		create_linked_list (this->enhanced_swing.final_locations, this->enhanced_swing.initial_locations, link, &master_linked_list, 1);
	}

	int L {0}; // to store length of linked_list 

	Particle* tmp_par_ptr_1 {nullptr}; // to help with swaps of pointers
	Particle* tmp_par_ptr_2 {nullptr}; // this guy as well 

	for ( std::vector <std::array<int,3>>& linked_list: master_linked_list ) {

		L = static_cast<int> (linked_list.size()) ;

		if ( (this->Lattice)[ lattice_index(linked_list[L-1], y, z) ]->ptype[0] == 's' ){

			tmp_par_ptr_1 = (this->Lattice) [ lattice_index (linked_list.back(), this->y, this->z) ];

			// go backwards 
			for (int i{0}; i < L ; i=i+2 ){

				(this->Lattice)[lattice_index (linked_list[L-2-i], this->y, this->z)]->coords = linked_list[L-1-i]; 
				(this->Lattice)[lattice_index (linked_list[L-1-i], this->y, this->z)] = (this->Lattice)[lattice_index( linked_list[L-2-i], this->y, this->z)]; 
				
			}
			tmp_par_ptr_1->coords = linked_list[0]; 
			(this->Lattice)[lattice_index (linked_list[0], this->y, this->z)] = tmp_par_ptr_1; 

		}

		else { // when it is a circulation of monomers 

			for ( int i{0}; i < L; i=i+2){
				
				if ( i == 0 ) {
					// store info about linked_list[1]... and consequently also linked_list[2]
					tmp_par_ptr_1 = (this->Lattice)[lattice_index ( linked_list[i+1], this->y, this->z) ];

					(this->Lattice)[ lattice_index ( linked_list[i],   this->y, this->z) ]->coords = linked_list[i+1]; 
					(this->Lattice)[ lattice_index ( linked_list[i+1], this->y, this->z) ] = (this->Lattice)[lattice_index ( linked_list[i], this->y, this->z) ]; 
				
				}
				else {

					tmp_par_ptr_1->coords = linked_list[i+1]; 
					tmp_par_ptr_2 = (this->Lattice)[ lattice_index (linked_list[i+1], this->y, this->z) ]; 
					(this->Lattice)[ lattice_index (linked_list[i+1], this->y, this->z) ] = tmp_par_ptr_1;
					tmp_par_ptr_1 = tmp_par_ptr_2; 
					
				}
			}
		}	
	}

	return;

}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::debug_backward_tail_regrowth(int p_idx, int m_idx, int recursion_depth){

	if (m_idx == 0){
		this->IMP_BOOL = true;
	}
	else{
		// reset the running boltzmann sum
		this->enhanced_swing.rboltzmann = 0;

		// reset locally
		this->enhanced_swing.reset_local();

		// set up some contact counters
		std::array <int,3>                loc_m   = this->Polymers[p_idx].chain[m_idx-1]->coords;
		std::array <std::array<int,3>,26> ne_list = obtain_ne_list(this->Polymers[p_idx].chain[m_idx]->coords, this->x, this->y, this->z);
		
		// set up the shuffle
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); 
		std::shuffle(ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

		// get the old position
		int ne_idx = std::find (ne_list.begin(), ne_list.end(), this->enhanced_swing.initial_locations[m_idx-1]) - ne_list.begin();
		std::array <int,3> tmp = ne_list[0]; 
		ne_list[0]             = ne_list[ne_idx];
		ne_list[ne_idx]        = tmp; 

		// set up some variables
		int idx_counter   =  0;
		int self_swap_idx = -1;

		while(idx_counter < this->enhanced_swing.ntest){
			if (ne_list[idx_counter] == loc_m){
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy;
				this->enhanced_swing.contacts_store[idx_counter] = this->enhanced_swing.current_contacts;
			}
			else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm') {
				for (int u{0}; u < this->Polymers[p_idx].deg_poly; ++u){
					if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
						self_swap_idx = u; 
						break;
					}
				}

				std::cout << "self swap monomer index = " << self_swap_idx << "." << std::endl;
				if (self_swap_idx > m_idx){
					this->enhanced_swing.energies[idx_counter]       = 1e+8;
					this->enhanced_swing.contacts_store[idx_counter] = {-1,-1,-1, -1,-1,-1, -1,-1};
				}
				else{

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

					// swap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// get the energies
					this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
					this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
					this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

					// set up the final energies
					this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
					this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
					this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
					std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

					this->Polymers[p_idx].print_chain();
					this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

					// revert back to original structure, reswap particles
					this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

					// reset some local variables
					this->enhanced_swing.reset_local();

				}
			}
			else {

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.initial_cont_part2switch), &(this->enhanced_swing.initial_E_part2switch));
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.initial_cont_monomer),     &(this->enhanced_swing.initial_E_monomer));
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.initial_cont_pair), &(this->enhanced_swing.initial_E_pair));

				// swap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &(this->enhanced_swing.final_cont_part2switch), &(this->enhanced_swing.final_E_part2switch)); 
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &(this->enhanced_swing.final_cont_monomer),     &(this->enhanced_swing.final_E_monomer)); 
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &(this->enhanced_swing.final_cont_pair), &(this->enhanced_swing.final_E_pair));

				// set up the final energies
				this->enhanced_swing.energies[idx_counter]       = this->enhanced_swing.current_energy - (this->enhanced_swing.initial_E_part2switch+this->enhanced_swing.initial_E_monomer-this->enhanced_swing.initial_E_pair) + (this->enhanced_swing.final_E_part2switch+this->enhanced_swing.final_E_monomer-this->enhanced_swing.final_E_pair);
				this->enhanced_swing.contacts_store[idx_counter] = add_containers(subtract_containers(this->enhanced_swing.current_contacts, add_containers (this->enhanced_swing.initial_cont_monomer, this->enhanced_swing.initial_cont_part2switch)), add_containers (this->enhanced_swing.final_cont_monomer, this->enhanced_swing.final_cont_part2switch));
				this->enhanced_swing.contacts_store[idx_counter] = subtract_containers(add_containers(this->enhanced_swing.contacts_store[idx_counter], this->enhanced_swing.initial_cont_pair), this->enhanced_swing.final_cont_pair);
				std::cout << "contacts_store[" << idx_counter << "] = "; print(this->enhanced_swing.contacts_store[idx_counter]);

				this->Polymers[p_idx].print_chain();
				this->debug_checks_energy_contacts(this->enhanced_swing.energies[idx_counter], this->enhanced_swing.contacts_store[idx_counter]);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// reset some local variables
				this->enhanced_swing.reset_local();

			}
			idx_counter += 1;
		}

		this->enhanced_swing.Emin = *std::min_element(this->enhanced_swing.energies.begin(), this->enhanced_swing.energies.end());
		for (int i{0}; i<this->enhanced_swing.ntest; ++i){
			this->enhanced_swing.boltzmann[i] = std::exp(-1/this->T * (this->enhanced_swing.energies[i]-this->enhanced_swing.Emin));
			this->enhanced_swing.rboltzmann += this->enhanced_swing.boltzmann[i];
		}

		this->enhanced_swing.sampler_idx      = 0;
		this->enhanced_swing.prob_o_to_n     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
		this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
		this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

		// do the swap to get to the new configuration
		this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->Polymers[p_idx].print_chain();
		this->debug_checks_energy_contacts(this->enhanced_swing.energies[this->enhanced_swing.sampler_idx], this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx]);
		this->debug_backward_tail_regrowth(p_idx, m_idx-1, recursion_depth+1);

	}
	return;

}

/////////////////////////////////////////////////////////////////////////////////////
