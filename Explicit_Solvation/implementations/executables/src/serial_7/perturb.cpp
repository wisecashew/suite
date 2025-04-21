#include "Simulation.hpp"
#include "lattice_directions.hpp"
#include "History.hpp"

//////////////////////////////////////////////////////////
//
//                Perturb the system
//
//////////////////////////////////////////////////////////

void Simulation::perturb_particle_swap_with_monomer(std::array<int,3> monomer_location, int polymer_idx, int monomer_idx, int other_particle_idx){

	// setting up the switch 
	this->Polymers[polymer_idx].chain[monomer_idx]->coords = location(other_particle_idx, this->x, this->y, this->z);
	this->Lattice[other_particle_idx]->coords = monomer_location;
	
	// do the switch 
	// take the pointer of the solvent, and put it where the monomer was on the lattice 
	this->Lattice[lattice_index(monomer_location, this->y, this->z)] = this->Lattice[other_particle_idx];
	this->Lattice[other_particle_idx]                                = this->Polymers[polymer_idx].chain[monomer_idx];

	return;
}

void Simulation::perturb_particle_swap_with_monomer_update(std::array<int,3> monomer_location, int polymer_idx, int monomer_idx, int other_particle_idx){

	// get the initial energetics
	this->neighbor_energetics(other_particle_idx,                                &(this->rotation_container.c1_initial), &(this->rotation_container.E1_initial));
	this->neighbor_energetics(lattice_index(monomer_location, this->y, this->z), &(this->rotation_container.c2_initial), &(this->rotation_container.E2_initial));
	this->selected_pair_interaction(this->Polymers[polymer_idx].chain[monomer_idx], this->Lattice[other_particle_idx], &(this->rotation_container.cpair_initial), &(this->rotation_container.Epair_initial));

	// perform the swap
	this->perturb_particle_swap_with_monomer(monomer_location, polymer_idx, monomer_idx, other_particle_idx);

	// get the energies in the new neighborhood
	this->neighbor_energetics(other_particle_idx,                                &(this->rotation_container.c1_final), &(this->rotation_container.E1_final));
	this->neighbor_energetics(lattice_index(monomer_location, this->y, this->z), &(this->rotation_container.c2_final), &(this->rotation_container.E2_final));
	this->selected_pair_interaction(this->Polymers[polymer_idx].chain[monomer_idx], this->Lattice[other_particle_idx], &(this->rotation_container.cpair_final), &(this->rotation_container.Epair_final));

	return;
}

void Simulation::perturb_particle_swap(Particle* tmp_par_ptr, int lat_idx_1, int lat_idx_2){

	tmp_par_ptr = (this->Lattice)[ lat_idx_1];
	
	(this->Lattice)[lat_idx_1] = (this->Lattice)[lat_idx_2];
	(this->Lattice)[lat_idx_1]->coords = location (lat_idx_1, this->x, this->y, this->z);

	(this->Lattice)[lat_idx_2] = tmp_par_ptr;
	(this->Lattice)[lat_idx_2]->coords = location (lat_idx_2, this->x, this->y, this->z);

	return;
}

void Simulation::perturb_particle_swap_with_update(Particle* tmp_par_ptr, int lat_idx_1, int lat_idx_2){

	// get the energies in the current neighborhood
	// std::cout << "Initial computation." << std::endl;
	this->neighbor_energetics(lat_idx_1, &(this->rotation_container.c1_initial), &(this->rotation_container.E1_initial));
	this->neighbor_energetics(lat_idx_2, &(this->rotation_container.c2_initial), &(this->rotation_container.E2_initial));
	this->selected_pair_interaction(this->Lattice[lat_idx_1], this->Lattice[lat_idx_2], &(this->rotation_container.cpair_initial), &(this->rotation_container.Epair_initial)); 

	// perform the swap
	this->perturb_particle_swap(tmp_par_ptr, lat_idx_1, lat_idx_2);

	// get the energies in the new neighborhood
	// std::cout << "Post swap computation." << std::endl;
	this->neighbor_energetics (lat_idx_1, &(this->rotation_container.c1_final), &(this->rotation_container.E1_final));
	this->neighbor_energetics (lat_idx_2, &(this->rotation_container.c2_final), &(this->rotation_container.E2_final)); 
	this->selected_pair_interaction(this->Lattice[lat_idx_1], this->Lattice[lat_idx_2], &(this->rotation_container.cpair_final), &(this->rotation_container.Epair_final));

	return;
}

//////////////////////////////////////////////////////////
// modular functions for sampling based on orientations
//////////////////////////////////////////////////////////

void Simulation::perturb_orientation_sampler_forwards(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iterator_idx, int lat_idx){

	this->enhanced_flipper.initial_orientations[iterator_idx] = this->Lattice[lat_idx]->orientation;
	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.initial_contacts), &(this->enhanced_flipper.initial_E));
	std::vector <int> test_orientations (26, 0);
	std::iota(test_orientations.begin(), test_orientations.end(), 0);

	// set up the random number generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(test_orientations.begin(), test_orientations.end(), std::default_random_engine(seed));

	// auto it = std::find(test_orientations.begin(), test_orientations.end(), this->enhanced_flipper.initial_orientations[iterator_idx]);
	// test_orientations.erase(it);

	for (int j{0}; j<this->enhanced_flipper.ntest; ++j){
		(this->Lattice)[lat_idx]->orientation = test_orientations[j];
		this->enhanced_flipper.orientations[j] = this->Lattice[lat_idx]->orientation;
		this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
		this->enhanced_flipper.energies[j] = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E;
		this->enhanced_flipper.contacts_store[j]  = subtract_containers(*contacts_sys, (this->enhanced_flipper.initial_contacts));
		this->enhanced_flipper.contacts_store[j]  = add_containers     ((this->enhanced_flipper.contacts_store[j]), (this->enhanced_flipper.perturbed_contacts));
		this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.perturbed_E        = 0;
	}

	return;
}

void Simulation::perturb_choose_state_forward(int iterator_idx, int lat_idx){

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
	this->enhanced_flipper.final_orientations[iterator_idx] = this->enhanced_flipper.orientations[this->enhanced_flipper.sampler_idx];
	this->Lattice[lat_idx]->orientation  = this->enhanced_flipper.orientations[this->enhanced_flipper.sampler_idx];
	this->enhanced_flipper.prob_o_to_n  *= this->enhanced_flipper.boltzmann[this->enhanced_flipper.sampler_idx]/this->enhanced_flipper.rboltzmann;

	return;
}

void Simulation::perturb_orientation_sampler_backwards_0(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx){

	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.initial_contacts), &(this->enhanced_flipper.initial_E));

	this->Lattice[lat_idx]->orientation = this->enhanced_flipper.initial_orientations[iteration_idx];
	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
	this->enhanced_flipper.energies[0] = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E;
	this->enhanced_flipper.contacts_store[0]  = subtract_containers(*contacts_sys, this->enhanced_flipper.initial_contacts);
	this->enhanced_flipper.contacts_store[0]  = add_containers     (this->enhanced_flipper.contacts_store[0], this->enhanced_flipper.perturbed_contacts);

	this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
	this->enhanced_flipper.perturbed_E        = 0;

	return;
}

void Simulation::perturb_orientation_sampler_backwards(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx){

	this->perturb_orientation_sampler_backwards_0(contacts_sys, E_sys, iteration_idx, lat_idx);
	std::vector <int> test_orientations (26, 0);
	std::iota(test_orientations.begin(), test_orientations.end(), 0);

	// set up the random number generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(test_orientations.begin(), test_orientations.end(), std::default_random_engine(seed));

	auto it = std::find(test_orientations.begin(), test_orientations.end(), this->enhanced_flipper.initial_orientations[iteration_idx]);
	test_orientations.erase(it);

	for (int j{1}; j<this->enhanced_flipper.ntest; ++j){
		this->Lattice[lat_idx]->orientation = test_orientations[j];
		this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
		this->enhanced_flipper.energies[j] = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E;
		this->enhanced_flipper.contacts_store [j]  = subtract_containers (*contacts_sys, this->enhanced_flipper.initial_contacts);
		this->enhanced_flipper.contacts_store [j]  = add_containers      (this->enhanced_flipper.contacts_store[j], this->enhanced_flipper.perturbed_contacts);
		this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.perturbed_E        = 0;
	}

	return;
}

//////////////////////////////////////////////////////////
// rotation an end of a polymer
//////////////////////////////////////////////////////////

void Simulation::perturb_tail_rotation(int p_idx){

	// add attempt
	this->attempts[0] += 1;

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

		final_E = this->sysEnergy - (this->rotation_container.E1_initial + this->rotation_container.E2_initial - this->rotation_container.Epair_initial);
		final_E = final_E + (this->rotation_container.E1_final   + this->rotation_container.E2_final  - this->rotation_container.Epair_final);
		final_contacts = add_containers(subtract_containers(this->contacts, add_containers (this->rotation_container.c1_initial, this->rotation_container.c2_initial)), add_containers (this->rotation_container.c1_final, this->rotation_container.c2_final));
		final_contacts = add_containers(final_contacts, this->rotation_container.cpair_initial);
		final_contacts = subtract_containers(final_contacts, this->rotation_container.cpair_final);

		if (metropolis_acceptance(this->sysEnergy, final_E, this->T)){
			this->sysEnergy = final_E;
			this->contacts  = final_contacts;
			this->IMP_BOOL  = true;
			this->acceptances[0] += 1;
		}
		else {
			// revert the polymer 
			this->perturb_particle_swap(tmp_par_ptr, lat_idx_1, lat_idx_2);
			this->IMP_BOOL = false;
		}
	}
	else {
		this->IMP_BOOL = false; 
	}

	return;
}

void Simulation::perturb_head_rotation(int p_idx){

	// add attempt
	this->attempts[0] += 1;

	// get the neighborlist of the particle at index-2
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

	// make the move
	if ( this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->ptype[0] == 's' ){

		this->perturb_particle_swap_with_update(tmp_par_ptr, lat_idx_1, lat_idx_2);
		
		// doing the quick manipulations to get the final energy and contacts
		final_E = this->sysEnergy - (this->rotation_container.E1_initial + this->rotation_container.E2_initial - this->rotation_container.Epair_initial);
		final_E = final_E + (this->rotation_container.E1_final   + this->rotation_container.E2_final  - this->rotation_container.Epair_final);
		final_contacts = add_containers(subtract_containers(this->contacts, add_containers (this->rotation_container.c1_initial, this->rotation_container.c2_initial)), add_containers (this->rotation_container.c1_final, this->rotation_container.c2_final));
		final_contacts = add_containers(final_contacts, this->rotation_container.cpair_initial);
		final_contacts = subtract_containers(final_contacts, this->rotation_container.cpair_final);

		if (metropolis_acceptance(this->sysEnergy, final_E, this->T)){
			this->sysEnergy = final_E;
			this->contacts  = final_contacts;
			this->IMP_BOOL  = true;
			this->acceptances[0] += 1;
		}
		else {
			// revert the polymer 
			this->perturb_particle_swap(tmp_par_ptr, lat_idx_1, lat_idx_2);
			this->IMP_BOOL = false;
		}
	}

	else {
		this->IMP_BOOL = false; 
	}

	return;
}

//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// perform reptation
//////////////////////////////////////////////////////////

void Simulation::revert_with_tail_biting(std::array<int,3>* to_slither, int p_idx){

	int deg_poly = this->Polymers[p_idx].deg_poly;
	for (int i{0}; i<deg_poly; ++i){
		// push everything forward
		// starting from the head
		if (i == deg_poly-1){
			this->Polymers[p_idx].chain[deg_poly-1]->coords = *to_slither;
			this->Lattice[lattice_index(*to_slither, this->y, this->z)] = this->Polymers[p_idx].chain[deg_poly-1];
		}
		else {
			this->Polymers[p_idx].chain[i]->coords = this->Polymers[p_idx].chain[i+1]->coords;
			this->Lattice[lattice_index(this->Polymers[p_idx].chain[i]->coords, this->y, this->z)] = this->Polymers[p_idx].chain[i];
		}
	}
	return; 
}

void Simulation::revert_without_tail_biting(Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* loc_0, int p_idx){

	int deg_poly = this->Polymers[p_idx].deg_poly;
	for (int i{0}; i<deg_poly; ++i){
		if (i == deg_poly-1){
			this->Polymers[p_idx].chain[deg_poly-1]->coords = *to_slither; 
			this->Lattice[ lattice_index (*to_slither, this->y, this->z) ] = this->Polymers[p_idx].chain[deg_poly-1];
		}
		else {
			this->Polymers[p_idx].chain[i]->coords = this->Polymers[p_idx].chain[i+1]->coords; 
			this->Lattice[ lattice_index (this->Polymers[p_idx].chain[i]->coords, this->y, this->z) ] = this->Polymers[p_idx].chain[i]; 
		} 
	}

	// put the solvent molecule back 
	this->Lattice[lattice_index (*loc_0, this->y, this->z)] = tmp; 
	this->Lattice[lattice_index (*loc_0, this->y, this->z)]->coords = *loc_0;

	return; 
}

void Simulation::revert_with_head_butting(std::array <int,3>* to_slither, int p_idx){

	int deg_poly = this->Polymers[p_idx].deg_poly;
	for (int i{0}; i<deg_poly; ++i){
		// push everything forward 
		// starting from the head 
		if ( deg_poly-1-i == 0 ){
			this->Polymers[p_idx].chain[deg_poly-1-i]->coords   = *to_slither; // to_slither = locf 
			this->Lattice[ lattice_index (*to_slither, this->y, this->z) ] = this->Polymers[p_idx].chain[deg_poly-1-i]; 
		}
		else {
			this->Polymers[p_idx].chain[deg_poly-1-i]->coords = this->Polymers[p_idx].chain[deg_poly-1-i-1]->coords; 
			this->Lattice[lattice_index(this->Polymers[p_idx].chain[deg_poly-1-i]->coords, this->y, this->z) ] = this->Polymers[p_idx].chain[deg_poly-1-i]; 
		}
	}
	return;
}

void Simulation::revert_without_head_butting(Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* locf, int p_idx){

	int deg_poly = this->Polymers[p_idx].deg_poly;
	for (int i{0}; i<deg_poly; ++i){
		// push everything forward 
		// starting from the head 
		if ( deg_poly-1-i == 0 ){
			this->Polymers[p_idx].chain[deg_poly-1-i]->coords   = *to_slither; // to_slither = locf 
			this->Lattice[ lattice_index (*to_slither, this->y, this->z) ] = this->Polymers[p_idx].chain[deg_poly-1-i]; 
		}
		else {
			this->Polymers[p_idx].chain[deg_poly-1-i]->coords = this->Polymers[p_idx].chain[deg_poly-1-i-1]->coords; 
			this->Lattice[ lattice_index (this->Polymers[p_idx].chain[deg_poly-1-i]->coords, this->y, this->z) ] = this->Polymers[p_idx].chain[deg_poly-1-i]; 
		}
	}

	this->Lattice[lattice_index(*locf, y, z)] = tmp; 
	this->Lattice[lattice_index(*locf, y, z)]->coords = *locf; 

	return;
}

void Simulation::forward_reptation_with_tail_biting(std::array<double,CONTACT_SIZE>* contacts_i, double* E_i, int p_idx){

	int deg_poly   = this->Polymers[p_idx].deg_poly;
	double Em1_i   = 0;
	double Em2_i   = 0;
	double Em1_f   = 0;
	double Em2_f   = 0;
	double Epair_i = 0;
	double Epair_f = 0;
	double Esys    = *E_i;

	std::array <double,CONTACT_SIZE> cm1_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm2_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_i   = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm1_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm2_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_f   = {0,0,0, 0,0,0, 0,0};
	std::array <int,3>    loc_m1    = {0,0,0};
	std::array <int,3>    loc_m2    = {0,0,0};
	std::array <double,CONTACT_SIZE> contacts  = *contacts_i;

	Particle* tmp_par_ptr {nullptr};

	for (int idx{0}; idx<deg_poly-1; ++idx) {

		// get locations of the two starting points
		loc_m1 = this->Polymers[p_idx].chain[idx]->coords;
		loc_m2 = this->Polymers[p_idx].chain[idx+1]->coords;

		// find energies
		this->neighbor_energetics(lattice_index(loc_m1, this->y, this->z), &cm1_i, &Em1_i);
		this->neighbor_energetics(lattice_index(loc_m2, this->y, this->z), &cm2_i, &Em2_i);
		this->selected_pair_interaction(this->Polymers[p_idx].chain[idx], this->Polymers[p_idx].chain[idx+1], &cpair_i, &Epair_i);

		// do the swap 
		this->perturb_particle_swap(tmp_par_ptr, lattice_index(loc_m1, this->y, this->z), lattice_index(loc_m2, this->y, this->z));

		// reevaluate energies
		this->neighbor_energetics(lattice_index(loc_m1, this->y, this->z), &cm1_f, &Em1_f);
		this->neighbor_energetics(lattice_index(loc_m2, this->y, this->z), &cm2_f, &Em2_f);
		this->selected_pair_interaction(this->Polymers[p_idx].chain[idx], this->Polymers[p_idx].chain[idx+1], &cpair_f, &Epair_f);

		// update energies and contacts
		contacts = add_containers(subtract_containers(contacts, add_containers(cm1_i, cm2_i)), add_containers(cm1_f, cm2_f));
		contacts = add_containers(contacts, cpair_i);
		contacts = subtract_containers(contacts, cpair_f);
		Esys = Esys - (Em1_i+Em2_i-Epair_i) + (Em1_f+Em2_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		// resetting everything
		cpair_i  = {0,0,0, 0,0,0, 0,0};
		cpair_f  = {0,0,0, 0,0,0, 0,0};
		cm1_i    = {0,0,0, 0,0,0, 0,0};
		cm1_f    = {0,0,0, 0,0,0, 0,0};
		cm2_i    = {0,0,0, 0,0,0, 0,0};
		cm2_f    = {0,0,0, 0,0,0, 0,0};
		Em1_i    = 0;
		Em2_i    = 0;
		Em1_f    = 0;
		Em2_f    = 0;
		Epair_i  = 0;
		Epair_f  = 0;
	}

	*E_i        = Esys;
	*contacts_i = contacts;
	return; 
}

void Simulation::forward_reptation_without_tail_biting(std::array<double,CONTACT_SIZE>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx){

	int    deg_poly = this->Polymers[p_idx].deg_poly;
	double Em_i     = 0;
	double Es_i     = 0;
	double Em_f     = 0;
	double Es_f     = 0;
	double Epair_i  = 0;
	double Epair_f  = 0;
	double Esys     = *E_i;

	std::array <double,CONTACT_SIZE> cm_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cs_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_i  = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cs_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_f  = {0,0,0, 0,0,0, 0,0};
	std::array <int,3>    loc_m    = {0,0,0};
	std::array <int,3>    loc_s    = *to_slither;
	std::array <double,CONTACT_SIZE> contacts = *contacts_i;

	for (int idx{0}; idx<deg_poly; ++idx){

		// std::cout << "idx = " << idx << ". " << std::endl;
		// get locations of the two starting points 
		loc_s = *to_slither; 
		loc_m = this->Polymers[p_idx].chain[deg_poly-1-idx]->coords;

		// std::cout << "Location of solvent = "; print(loc_s);
		// std::cout << "Location of monomer = "; print(loc_m);

		// find energies
		this->neighbor_energetics(lattice_index(loc_m, this->y, this->z), &cm_i, &Em_i);
		this->neighbor_energetics(lattice_index(loc_s, this->y, this->z), &cs_i, &Es_i);
		this->selected_pair_interaction(this->Lattice[lattice_index(loc_m, this->y, this->z)], this->Lattice[lattice_index(loc_s, this->y, this->z)], &cpair_i, &Epair_i);

		// do the swap 
		this->Lattice[lattice_index(loc_m, this->y, this->z)] = this->Lattice[lattice_index(loc_s, this->y, this->z)];
		this->Lattice[lattice_index(loc_s, this->y, this->z)] = this->Polymers[p_idx].chain[deg_poly-1-idx];
		this->Lattice[lattice_index(loc_m, this->y, this->z)]->coords = loc_m;
		this->Lattice[lattice_index(loc_s, this->y, this->z)]->coords = loc_s;

		// reevaluate energies
		this->neighbor_energetics(lattice_index(loc_s, this->y, this->z), &cs_f, &Es_f);
		this->neighbor_energetics(lattice_index(loc_m, this->y, this->z), &cm_f, &Em_f);
		this->selected_pair_interaction(this->Lattice[lattice_index(loc_s, this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_f, &Epair_f);

		// update energies and contacts
		contacts = add_containers(subtract_containers(contacts, add_containers(cm_i, cs_i)), add_containers(cm_f, cs_f));
		contacts = add_containers(contacts, cpair_i);
		contacts = subtract_containers(contacts, cpair_f);

		Esys = Esys - (Em_i+Es_i-Epair_i) + (Em_f+Es_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		// the particle is now at the old position of the monomer bead
		*to_slither = loc_m; 

		// resetting everything
		cpair_i = {0,0,0, 0,0,0, 0,0};
		cpair_f = {0,0,0, 0,0,0, 0,0};
		cs_i    = {0,0,0, 0,0,0, 0,0};
		cs_f    = {0,0,0, 0,0,0, 0,0};
		cm_i    = {0,0,0, 0,0,0, 0,0};
		cm_f    = {0,0,0, 0,0,0, 0,0};
		Em_i    = 0;
		Es_i    = 0;
		Em_f    = 0;
		Es_f    = 0;
		Epair_i = 0;
		Epair_f = 0;
	}

	*E_i        = Esys;
	*contacts_i = contacts;

	return; 
}

void Simulation::backward_reptation_with_head_butting(std::array<double,CONTACT_SIZE>* contacts_i, double* E_i, int p_idx){

	int    deg_poly    = this->Polymers[p_idx].deg_poly;
	double Em1_i       = 0;
	double Em2_i       = 0;
	double Em1_f       = 0;
	double Em2_f       = 0;
	double Epair_i     = 0;
	double Epair_f     = 0;
	double Esys        = *E_i;

	std::array <double,CONTACT_SIZE> cm1_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm2_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_i   = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm1_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm2_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_f   = {0,0,0, 0,0,0, 0,0};
	std::array <int,3>    loc_m1    = {0,0,0}; 
	std::array <int,3>    loc_m2    = {0,0,0}; 
	std::array <double,CONTACT_SIZE> contacts  = *contacts_i;

	for (int i{0}; i<deg_poly-1; ++i){

		// get locations of the the two starting points
		loc_m1 = this->Polymers[p_idx].chain[deg_poly-(i+1)]->coords;
		loc_m2 = this->Polymers[p_idx].chain[deg_poly-(i+2)]->coords;

		// find energies
		this->neighbor_energetics(lattice_index(loc_m1, this->y, this->z), &cm1_i, &Em1_i);
		this->neighbor_energetics(lattice_index(loc_m2, this->y, this->z), &cm1_i, &Em2_i);
		this->selected_pair_interaction(this->Lattice[lattice_index(loc_m1, this->y, this->z)], this->Lattice[lattice_index(loc_m2, this->y, this->z)], &cpair_i, &Epair_i);

		// make the swap
		this->Lattice[lattice_index(loc_m1, this->y, this->z)] = this->Lattice[lattice_index(loc_m2, this->y, this->z)];
		this->Lattice[lattice_index(loc_m2, this->y, this->z)] = this->Polymers[p_idx].chain[deg_poly-(i+1)];
		this->Lattice[lattice_index(loc_m1, this->y, this->z)]->coords = loc_m1;
		this->Lattice[lattice_index(loc_m2, this->y, this->z)]->coords = loc_m2;

		// find new energies
		this->neighbor_energetics(lattice_index(loc_m1, this->y, this->z), &cm1_f, &Em1_f);
		this->neighbor_energetics(lattice_index(loc_m2, this->y, this->z), &cm2_f, &Em2_f);
		this->selected_pair_interaction(this->Lattice[lattice_index(loc_m1, this->y, this->z)], this->Lattice[lattice_index(loc_m2, this->y, this->z)], &cpair_f, &Epair_f);

		contacts = add_containers(subtract_containers(contacts, add_containers(cm1_i, cm2_i)), add_containers(cm1_f, cm2_f));
		contacts = add_containers(contacts, cpair_i);
		contacts = subtract_containers(contacts, cpair_f);
		Esys = Esys-(Em1_i+Em2_i-Epair_i)+(Em1_f+Em2_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		// resetting everything
		cpair_i  = {0,0,0, 0,0,0, 0,0};
		cpair_f  = {0,0,0, 0,0,0, 0,0};
		cm1_i    = {0,0,0, 0,0,0, 0,0};
		cm1_f    = {0,0,0, 0,0,0, 0,0};
		cm2_i    = {0,0,0, 0,0,0, 0,0};
		cm2_f    = {0,0,0, 0,0,0, 0,0};
		Em1_i    = 0;
		Em2_i    = 0;
		Em1_f    = 0;
		Em2_f    = 0;
		Epair_i  = 0;
		Epair_f  = 0;

	}

	*contacts_i = contacts;
	*E_i        = Esys;
	return; 
}

void Simulation::backward_reptation_without_head_butting(std::array<double,CONTACT_SIZE>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx){

	int    deg_poly    = this->Polymers[p_idx].deg_poly;
	double Em_i        = 0;
	double Es_i        = 0; 
	double Em_f        = 0; 
	double Es_f        = 0; 
	double Epair_i     = 0;
	double Epair_f     = 0;
	double Esys        = *E_i; 

	std::array <double,CONTACT_SIZE> cm_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cs_i     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_i  = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cm_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cs_f     = {0,0,0, 0,0,0, 0,0};
	std::array <double,CONTACT_SIZE> cpair_f  = {0,0,0, 0,0,0, 0,0};
	std::array <int,3>    loc_m    = {0,0,0};
	std::array <int,3>    loc_s    = *to_slither;
	std::array <double,CONTACT_SIZE> contacts = *contacts_i;

	for (int i{0}; i<deg_poly; ++i){

		loc_s = *to_slither; 
		loc_m = this->Polymers[p_idx].chain[i]->coords;

		// find energetics
		this->neighbor_energetics(lattice_index(loc_m, this->y, this->z), &cm_i, &Em_i);
		this->neighbor_energetics(lattice_index(loc_s, this->y, this->z), &cs_i, &Es_i);
		this->selected_pair_interaction(this->Lattice[lattice_index(loc_m, this->y, this->z)], this->Lattice[lattice_index(loc_s, this->y, this->z)], &cpair_i, &Epair_i);

		// make the swap
		this->Lattice[lattice_index(loc_m, this->y, this->z)] = this->Lattice[lattice_index(loc_s, this->y, this->z)];
		this->Lattice[lattice_index(loc_s, this->y, this->z)] = this->Polymers[p_idx].chain[i];
		this->Lattice[lattice_index(loc_m, this->y, this->z)]->coords = loc_m;
		this->Lattice[lattice_index(loc_s, this->y, this->z)]->coords = loc_s;

		// find new energies
		this->neighbor_energetics(lattice_index(loc_m, this->y, this->z), &cm_f, &Em_f);
		this->neighbor_energetics(lattice_index(loc_s, this->y, this->z), &cs_f, &Es_f);
		this->selected_pair_interaction(this->Lattice[lattice_index(loc_m, this->y, this->z)], this->Lattice[lattice_index(loc_s, this->y, this->z)], &cpair_f, &Epair_f);

		contacts = add_containers(subtract_containers(contacts, add_containers(cm_i, cs_i)), add_containers(cm_f, cs_f));
		contacts = add_containers(contacts, cpair_i);
		contacts = subtract_containers(contacts, cpair_f);		
		Esys = Esys-(Em_i+Es_i-Epair_i)+(Em_f+Es_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		*to_slither = loc_m; 

		// resetting everything
		cpair_i = {0,0,0, 0,0,0, 0,0};
		cpair_f = {0,0,0, 0,0,0, 0,0};
		cs_i    = {0,0,0, 0,0,0, 0,0};
		cs_f    = {0,0,0, 0,0,0, 0,0};
		cm_i    = {0,0,0, 0,0,0, 0,0};
		cm_f    = {0,0,0, 0,0,0, 0,0};
		Em_i    = 0;
		Es_i    = 0;
		Em_f    = 0;
		Es_f    = 0;
		Epair_i = 0;
		Epair_f = 0;

	}
	
	*E_i        = Esys;
	*contacts_i = contacts;

	return; 
}

void Simulation::perturb_reptation_forward(int p_idx){

	// update the number of attempts
	this->attempts[1] += 1;

	// instantiate deg_poly because this will be used later
	int                  deg_poly   = this->Polymers[p_idx].deg_poly;
	
	// instantiate a copy of current energy
	double               E_i        = this->sysEnergy;

	// instantiate the locations of the first and final monomer
	std::array<int,3>    loc_i      = this->Polymers[p_idx].chain[0]->coords;
	std::array<int,3>    loc_f      = this->Polymers[p_idx].chain[deg_poly-1]->coords;

	// instantiate a copy of current contacts
	std::array<double,CONTACT_SIZE> contacts_i = this->contacts; 

	// get the possibilities of how the polymer can slither forward
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_f, this->x, this->y, this->z);
	int choice = rng_uniform(0, 25);

	// generate the position to slither to, and its copy!
	std::array <int,3> to_slither      = ne_list[choice];
	std::array <int,3> to_slither_copy = to_slither;

	// this will lead to the head biting the tail unless careful
	if (to_slither == loc_i){
		this->forward_reptation_with_tail_biting(&contacts_i, &E_i, p_idx);
	}
	// you have to displace a solvent molecule to move forward
	else if (this->Lattice[lattice_index(to_slither, this->y, this->z)]->ptype[0] == 's'){
		this->forward_reptation_without_tail_biting(&contacts_i, &to_slither, &E_i, p_idx);
	}

	// if you bump into some central monomer, then the move is rejected
	else {
		// nothing to do
		this->IMP_BOOL = false;
		return;
	}

	if (metropolis_acceptance(this->sysEnergy, E_i, this->T)){
		this->sysEnergy = E_i;
		this->contacts  = contacts_i;
		this->IMP_BOOL  = true;
		this->acceptances[1] += 1;
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

	return;
}

void Simulation::perturb_reptation_backward(int p_idx){

	// update the number of attempts
	this->attempts[1] += 1;

	// instantiate deg_poly because this will be used later
	int                   deg_poly         = this->Polymers[p_idx].deg_poly;

	// instantiate a copy of current energy
	double                E_i              = this->sysEnergy;

	// instantiate the locations of the first and final monomer
	std::array <int,3>    loc_i            = this->Polymers[p_idx].chain[0]->coords;
	std::array <int,3>    loc_f            = this->Polymers[p_idx].chain[deg_poly-1]->coords;

	// instantiate a copy of current contacts
	std::array <double,CONTACT_SIZE> contacts_i       = this->contacts;


	// first check if tail rotation can be performed at all 
	std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( loc_i, x, y, z );
	int choice = rng_uniform (0, 25);

	std::array <int,3> to_slither      = ne_list[choice]; 
	std::array <int,3> to_slither_copy = to_slither; 

	if (to_slither == loc_f){
		this->backward_reptation_with_head_butting (&contacts_i, &E_i, p_idx);
	}
	else if ( (this->Lattice)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		this->backward_reptation_without_head_butting (&contacts_i, &to_slither, &E_i, p_idx);
	}

	else {
		// failed reptation
		this->IMP_BOOL = false; 
		return; 
	}

	if ( metropolis_acceptance (this->sysEnergy, E_i, this->T) ){
		this->sysEnergy = E_i;
		this->contacts  = contacts_i;
		this->IMP_BOOL  = true;
		this->acceptances[1] += 1;
	}
	else {
		if (to_slither_copy == loc_f){
			this->revert_with_tail_biting(&to_slither_copy, p_idx);
		}
		else {
			Particle* tmp {nullptr};
			tmp = this->Lattice[lattice_index(loc_f, y, z)];
			this->revert_without_tail_biting(tmp, &loc_f, &to_slither_copy, p_idx);
		}
		this->IMP_BOOL = false;
	}

	return;
}

//////////////////////////////////////////////////////////
// perform orientation flip
//////////////////////////////////////////////////////////

void Simulation::perturb_polymer_orientation_flip(int p_idx){

	// update the number of attempts
	this->attempts[2] += 1;

	// instantiate some useful variables
	int deg_poly = this->Polymers[p_idx].deg_poly;
	std::vector<int> polymer_indices (deg_poly);
	std::iota(polymer_indices.begin(), polymer_indices.end(), 0);

	// std::cout << "polymer indices = "; print(polymer_indices);

	// set up the random number generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(polymer_indices.begin(), polymer_indices.end(), std::default_random_engine(seed));

	// instantiate variables for boltzmann sampling
	int ntest    = 5;
	int nflip    = (deg_poly==1) ? 1 : rng_uniform(1, deg_poly-1);

	// set up the flipper object
	this->enhanced_flipper.reset(nflip, ntest);

	// define some holders 
	double               E_sys        = this->sysEnergy;
	std::array<double,CONTACT_SIZE> contacts_sys = this->contacts;

	// some more holders
	double               E_post        = 0;                                // the energy of the system after final perturbation
	double               rng_acc       = 0;                                // rng for acceptance at the very end
	std::array<double,CONTACT_SIZE> contacts_post = {0,0,0, 0,0,0, 0,0}; // the contacts of the system after final perturbation

	// relevant indices variable store
	int m_lattice_idx = 0; 

	// loop over the different monomer indices
	for (int i{0}; i < nflip; ++i){

		// sample different orientations and get the boltzmann factors 
		m_lattice_idx = lattice_index(this->Polymers[p_idx].chain[polymer_indices[i]]->coords, this->y, this->z);

		// make this without debug
		this->perturb_orientation_sampler_forwards(&contacts_sys, E_sys, i, m_lattice_idx);

		// get the energy minima
		this->enhanced_flipper.Emin = *std::min_element(this->enhanced_flipper.energies.begin(), this->enhanced_flipper.energies.end());

		// make this without debug
		this->perturb_choose_state_forward(i, m_lattice_idx);

		// set the system energy at the final point
		E_sys        = this->enhanced_flipper.energies[enhanced_flipper.sampler_idx];
		contacts_sys = this->enhanced_flipper.contacts_store[enhanced_flipper.sampler_idx];

		// reset
		this->enhanced_flipper.initial_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.initial_E        =  0;
		this->enhanced_flipper.rboltzmann       =  0;
		this->enhanced_flipper.sampler_rsum     =  0;

	}
	
	// store up the energies and contacts
	E_post        = E_sys;
	contacts_post = contacts_sys;

	// reset
	this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
	this->enhanced_flipper.perturbed_E        = 0;


	for (int i{0}; i < nflip; ++i){

		// sample different orientations and get the boltzmann factors for the backward flux
		m_lattice_idx = lattice_index(this->Polymers[p_idx].chain[polymer_indices[i]]->coords, this->y, this->z);

		// make this without debug
		this->perturb_orientation_sampler_backwards(&contacts_sys, E_sys, i, m_lattice_idx);

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

		// reset
		this->enhanced_flipper.initial_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.initial_E        = 0;
		this->enhanced_flipper.rboltzmann       = 0;

	}

	// get a random number
	rng_acc = rng_uniform(0.0, 1.0);

	// run the criterion
	if ( rng_acc < std::exp (-1/this->T * (E_post - this->sysEnergy)) * this->enhanced_flipper.prob_n_to_o/this->enhanced_flipper.prob_o_to_n){

		// if (E_post > this->sysEnergy){
		//	std::cout << "prob_n_to_o = " << this->enhanced_flipper.prob_n_to_o << std::endl;
		//	std::cout << "prob_o_to_n = " << this->enhanced_flipper.prob_o_to_n << std::endl;
		//	std::cout << "E_post = " << E_post << ", sysEnergy = " << this->sysEnergy << std::endl;
		// }

		for (int j{0}; j < nflip; ++j){
			(this->Polymers)[p_idx].chain[polymer_indices[j]]->orientation = this->enhanced_flipper.final_orientations[j];
		}
		this->sysEnergy       = E_post;
		this->contacts        = contacts_post;
		this->IMP_BOOL        = true;
		this->acceptances[2] += 1;
	}
	else {
		this->IMP_BOOL = false; 
	}

	return;
}

void Simulation::perturb_solvation_shell_flip(){

	// update the number of attempts
	this->attempts[3] += 1;

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
	double               E_sys        = this->sysEnergy;
	std::array<double,CONTACT_SIZE> contacts_sys = this->contacts;

	double               E_post        = 0; // the energy of the system after final perturbation
	std::array<double,CONTACT_SIZE> contacts_post = {0,0,0, 0,0,0, 0,0}; // the contacts of the system after final perturbation
	double               rng_acc       = 0; // rng for acceptance at the very end

	for (int i{0}; i < nflip; ++i){

		// sample different orientations and get the boltzmann factors 
		// make without debug
		this->perturb_orientation_sampler_forwards(&contacts_sys, E_sys, i, solvation_shell_indices[i]);

		// get the energy minima
		this->enhanced_flipper.Emin = *std::min_element (this->enhanced_flipper.energies.begin(), this->enhanced_flipper.energies.end() );

		// go to the new state
		// make without debug
		this->perturb_choose_state_forward(i, solvation_shell_indices[i]);

		// set the system energy at the final point
		E_sys         = this->enhanced_flipper.energies[enhanced_flipper.sampler_idx];
		contacts_sys  = this->enhanced_flipper.contacts_store[this->enhanced_flipper.sampler_idx];

		// reset
		this->enhanced_flipper.initial_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.initial_E        = 0;
		this->enhanced_flipper.rboltzmann       = 0;
		this->enhanced_flipper.sampler_rsum     = 0;

	}

	E_post = E_sys;
	contacts_post = contacts_sys;

	this->enhanced_flipper.perturbed_contacts = {0,0,0, 0,0,0, 0,0};
	this->enhanced_flipper.perturbed_E        = 0;


	for (int i{0}; i<nflip; ++i){

		// sample different orientations and get the boltzmann factors for the backward flux
		// make without debug
		this->perturb_orientation_sampler_backwards(&contacts_sys, E_sys, i, solvation_shell_indices[i]);

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

		// reset
		this->enhanced_flipper.initial_contacts = {0,0,0, 0,0,0, 0,0};
		this->enhanced_flipper.initial_E        = 0;
		this->enhanced_flipper.rboltzmann       = 0;

	}

	rng_acc = rng_uniform(0.0, 1.0);
	if ( rng_acc < std::exp (-1/this->T * (E_post - this->sysEnergy)) * this->enhanced_flipper.prob_n_to_o/this->enhanced_flipper.prob_o_to_n){

		// if accepted, return to the new orientations 
		for (int j{0}; j < nflip; ++j){
			this->Lattice[solvation_shell_indices[j]]->orientation = this->enhanced_flipper.final_orientations[j]; 
		}

		this->sysEnergy       = E_post;
		this->contacts        = contacts_post;
		this->IMP_BOOL        = true;
		this->acceptances[3] += 1;

	}
	else {
		this->IMP_BOOL = false;
	}

	return;
}

void Simulation::perturb_lattice_flip(){

	// update the number of attempts
	this->attempts[4] += 1;

	// number of particles to flip
	int nflips = rng_uniform(1, static_cast<int>(this->x * this->y * this->z / 50));
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

	for (int i{0}; i<nflips; ++i){
		
		// get the neighbor energy 
		this->neighbor_energetics(samples[i], &cs_i, &Es_i);

		// store the old orientations 
		orientations[i] = (this->Lattice[samples[i]])->orientation; // .push_back(this->Lattice[samples[i]]->orientation); 

		// perturb the orientation
		this->Lattice[samples[i]]->orientation = rng_uniform(0, 25);

		// get the new neighbor energy
		this->neighbor_energetics(samples[i], &cs_f, &Es_f); 

		// update the contacts and energy
		contacts = add_containers(subtract_containers(contacts, cs_i), cs_f);
		Ef += Es_f - Es_i;

		// resetting...
		cs_i = {0,0,0, 0,0,0, 0,0};
		cs_f = {0,0,0, 0,0,0, 0,0};
		Es_i = 0;
		Es_f = 0;

	}

	double rng = rng_uniform(0.0, 1.0);

	if (rng < std::exp(-1/this->T * (Ef - this->sysEnergy))){
		this->sysEnergy       = Ef;
		this->contacts        = contacts;
		this->IMP_BOOL        = true;
		this->acceptances[4] += 1;
	}
	else {
		// go back to old orientations
		for (int i{0}; i<nflips; ++i){
			this->Lattice[samples[i]]->orientation = orientations[i];
		}
		this->IMP_BOOL = false;
	}

	return;
}

void Simulation::perturb_solvent_exchange_from_shell(){

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

		// update the number of attempts
		this->attempts[5] += 1;

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

			rng_acc = rng_uniform (0.0, 1.0);

			if (rng_acc < std::exp (-1/this->T * (post_energy - this->sysEnergy))){
				this->sysEnergy = post_energy;
				this->contacts  = post_contacts;
				this->IMP_BOOL  = true;
				this->acceptances[5] += 1;
			}
			else {
				this->perturb_particle_swap(tmp_par_ptr, solvation_shell_indices[n], exc_idx);
				this->IMP_BOOL = false;
			}

			// resetting 
			this->rotation_container.reset();
		}
	}
	return;
}

void Simulation::perturb_solvent_exchange(){

	// these are the indices on the lattice for the particles to swap
	int idx1 = -1;
	int idx2 = -1; 

	// this is the number of particles to switch up
	int nswitches = rng_uniform(1, static_cast<int>(10 * std::log10(this->x*this->y*this->z)));

	// reset the container
	this->rotation_container.reset();

	// these are the variables 
	double post_energy = 0;
	double rng_acc     = 0;
	std::array <double,CONTACT_SIZE> post_contacts = this->contacts;

	// set up the pointer to facilitate the swap
	Particle* tmp_par_ptr {nullptr};

	for (int j{0}; j<nswitches; ++j){
		// update the number of attempts
		this->attempts[6] += 1;

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

			rng_acc = rng_uniform (0.0, 1.0); 
			if (rng_acc < std::exp (-1/this->T * (post_energy - this->sysEnergy))){
				this->sysEnergy = post_energy;
				this->contacts  = post_contacts;
				this->acceptances[6] += 1;
		
			}
			else {
				this->perturb_particle_swap(tmp_par_ptr, idx1, idx2);
			}
			this->rotation_container.reset();
		}

	}

	return; 
}

//////////////////////////////////////////////////////////
// perform regrowth
//////////////////////////////////////////////////////////

void Simulation::perturb_regrowth(int p_idx){

	// update the number of attempts
	this->attempts[7] += 1;

	// set up some containers 
	int deg_poly        = this->Polymers[p_idx].deg_poly; 
	int m_idx           =  rng_uniform(1, deg_poly-2);
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
		// std::cout << "Head regrowth..." << std::endl;
		for (int i{m_idx+1}; i<deg_poly; ++i){
			this->enhanced_swing.initial_locations.push_back(this->Polymers[p_idx].chain[i]->coords);
		}

		// perform the forward head regrowth
		this->perturb_forward_head_regrowth(p_idx, m_idx);

		for (int i{m_idx+1}; i<deg_poly; ++i){
			this->enhanced_swing.final_locations.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		if (this->enhanced_swing.initial_locations == this->enhanced_swing.final_locations){
			// do nothing
			;
		}

		else if (!(this->IMP_BOOL)){
			this->perturb_accept_after_head_regrowth(this->IMP_BOOL);
		}

		else {
			forw_energy   = this->enhanced_swing.current_energy;
			forw_contacts = this->enhanced_swing.current_contacts;
			this->perturb_backward_head_regrowth(p_idx, m_idx, 0);

			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * this->enhanced_swing.prob_n_to_o / this->enhanced_swing.prob_o_to_n){
				this->perturb_accept_after_head_regrowth(this->IMP_BOOL);
				this->sysEnergy = forw_energy; 
				this->contacts  = forw_contacts;
				this->acceptances[7] += 1;
			}
			else {
				this->IMP_BOOL = false;
			}
		}
	}

	else {
		// std::cout << "Tail regrowth..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		for (int i{0}; i<m_idx; ++i){
			this->enhanced_swing.initial_locations.push_back(this->Polymers[p_idx].chain[i]->coords);
		}

		// perform the tail regrowth
		this->perturb_forward_tail_regrowth(p_idx, m_idx);

		// std::cout << "Filling up final locations..." << std::endl;
		for (int i{0}; i<m_idx; ++i){
			this->enhanced_swing.final_locations.push_back(this->Polymers[p_idx].chain[i]->coords);
		}

		if (this->enhanced_swing.initial_locations == this->enhanced_swing.final_locations){
			// std::cout << "No change in final and initial." << std::endl;
			// do nothing
			;
		}

		else if (!(this->IMP_BOOL)){
			// std::cout << "Reverting because of premature termination." << std::endl;
			this->perturb_accept_after_tail_regrowth(this->IMP_BOOL);
		}

		else {
			
			forw_energy   = this->enhanced_swing.current_energy;
			forw_contacts = this->enhanced_swing.current_contacts;
			this->perturb_backward_tail_regrowth(p_idx, m_idx, 0);
			
			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * this->enhanced_swing.prob_n_to_o/this->enhanced_swing.prob_o_to_n) {
				this->perturb_accept_after_tail_regrowth(this->IMP_BOOL);
				this->sysEnergy = forw_energy;
				this->contacts  = forw_contacts;
				this->acceptances[7] += 1;
			}

			else {
				this->IMP_BOOL = false;
			}
		}
	}
	return;
}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::perturb_forward_head_regrowth(int p_idx, int m_idx){

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

				if (self_swap_idx < m_idx){
					this->enhanced_swing.energies[idx_counter]       = 1e+8;
					this->enhanced_swing.contacts_store[idx_counter] = {-1,-1,-1, -1,-1,-1, -1,-1};
					block_counter                                   += 1;
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
			for (int i{0}; i<this->enhanced_swing.ntest; ++i) {
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
			
			// update probability fluxes
			this->enhanced_swing.prob_o_to_n     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
			this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
			this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

			// do the swap to get to the new configuration
			this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
			this->perturb_forward_head_regrowth(p_idx, m_idx+1);

		}

	}

	return;
}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::perturb_accept_after_head_regrowth(bool not_trap_bool){

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

void Simulation::perturb_backward_head_regrowth(int p_idx, int m_idx, int recursion_depth){

	if (m_idx == this->Polymers[p_idx].deg_poly - 1){
		this->IMP_BOOL = true;
	}
	else{
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

				if (self_swap_idx <= m_idx){
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
		this->enhanced_swing.prob_n_to_o     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
		this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
		this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

		// do the swap to get to the new configuration
		this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->perturb_backward_head_regrowth(p_idx, m_idx+1, recursion_depth+1);

	}
	return;
}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::perturb_forward_tail_regrowth(int p_idx, int m_idx){

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

			// update probability flux and energy
			this->enhanced_swing.prob_o_to_n     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
			this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
			this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

			// do the swap to get to the new configuration
			this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
			this->perturb_forward_tail_regrowth(p_idx, m_idx-1);
		}
	}

	return;
}

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::perturb_accept_after_tail_regrowth(bool not_trap_bool){

	// set up some vectors
	std::vector< std::vector<std::array<int,3>>> master_linked_list;
	std::vector <std::array<int,3>> link; 

	// get the linked list
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

//#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void Simulation::perturb_backward_tail_regrowth(int p_idx, int m_idx, int recursion_depth){

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

				if (self_swap_idx > m_idx){
					this->enhanced_swing.energies[idx_counter]       = 1e+8;
					this->enhanced_swing.contacts_store[idx_counter] = {-1,-1,-1, -1,-1,-1, -1,-1};
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
		this->enhanced_swing.prob_n_to_o     *= this->enhanced_swing.boltzmann[this->enhanced_swing.sampler_idx]/this->enhanced_swing.rboltzmann;
		this->enhanced_swing.current_energy   = this->enhanced_swing.energies[this->enhanced_swing.sampler_idx];
		this->enhanced_swing.current_contacts = this->enhanced_swing.contacts_store[this->enhanced_swing.sampler_idx];

		// do the swap to get to the new configuration
		this->perturb_particle_swap(this->enhanced_swing.tmp_par_ptr, lattice_index(ne_list[this->enhanced_swing.sampler_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->perturb_backward_tail_regrowth(p_idx, m_idx-1, recursion_depth+1);

	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//
//                Run the simulation
//
//////////////////////////////////////////////////////////

void Simulation::perturb_system_isotropic(){

	int r = rng_uniform(0, 4);
	switch (r) {
		case 0:
			this->perturb_tail_rotation(0);
			break;
		case 1:
			this->perturb_head_rotation(0);
			break;
		case 2:
			this->perturb_reptation_forward(0);
			break;
		case 3:
			this->perturb_reptation_backward(0);
			break;
		case 4:
			this->perturb_regrowth(0);
			break;
		default:
			std::cout << "Bad random number generated. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
			break;
	}
	return; 
}

// ~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%

void Simulation::perturb_system_straight(){

	int r = rng_uniform(0, 9);
	switch (r) {
		case 0:
			this->perturb_tail_rotation(0);
			break;
		case 1:
			this->perturb_head_rotation(0);
			break;
		case 2:
			this->perturb_reptation_forward(0);
			break;
		case 3:
			this->perturb_reptation_backward(0);
			break;
		case 4:
			this->perturb_polymer_orientation_flip(0);
			break;
		case 5:
			this->perturb_solvation_shell_flip();
			break;
		case 6:
			this->perturb_lattice_flip();
			break;
		case 7:
			this->perturb_solvent_exchange_from_shell();
			break;
		case 8:
			this->perturb_solvent_exchange();
			break;
		case 9:
			this->perturb_regrowth(0);
			break;
		default:
			std::cout << "Bad random number generated. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
			break;
	}
	return;
}

// ~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%

void Simulation::perturb_system_debug(){
	std::cout << "Running the debugging-based perturbation." << std::endl;
	int r = rng_uniform(0,8);
	
	switch (r) {
		case 0:
			std::cout << "Doing the tail rotation..." << std::endl;
			this->debug_tail_rotation(0);
			break;
		case 1:
			std::cout << "Doing the head rotation..." << std::endl;
			this->debug_head_rotation(0);
			break;
		
		case 2:
			std::cout << "Doing the forward reptation..." << std::endl;
			this->debug_reptation_forward(0);
			break;
		case 3:
			std::cout << "Doing the backward reptation..." << std::endl;
			this->debug_reptation_backward(0);
			break;
		case 4:
			std::cout << "Flip polymer orientations..." << std::endl;
			this->debug_polymer_orientation_flip(0);
			break;
		
		case 5:
			std::cout << "Flip the solvation shell..." << std::endl;
			this->debug_solvation_shell_flip();
			break;
		
		case 6:
			std::cout << "Perform a lattice flip..." << std::endl;
			this->debug_lattice_flip();
			break;
		
		case 7:
			std::cout << "Perform a solvent exchange..." << std::endl;
			this->debug_solvent_exchange_from_shell();
			break;
		
		case 8:
			std::cout << "Perform a solvent exchange..." << std::endl;
			this->debug_solvent_exchange();
			break;
		
		case 9:
			std::cout << "Perform a regrowth..." << std::endl;
			this->debug_regrowth(0);
			break;
		
		default:
			std::cout << "Bad random number generated. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
			break;
	}
	return;
}

// ~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%

