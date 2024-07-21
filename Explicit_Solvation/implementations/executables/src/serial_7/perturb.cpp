#include "Simulation.h"
#include "lattice_directions.h"
#include "History.h"

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
	std::cout << "Initial computation." << std::endl;
	this->neighbor_energetics(lat_idx_1, &(this->rotation_container.c1_initial), &(this->rotation_container.E1_initial));
	this->neighbor_energetics(lat_idx_2, &(this->rotation_container.c2_initial), &(this->rotation_container.E2_initial));
	this->selected_pair_interaction(this->Lattice[lat_idx_1], this->Lattice[lat_idx_2], &(this->rotation_container.cpair_initial), &(this->rotation_container.Epair_initial)); 

	// perform the swap
	this->perturb_particle_swap(tmp_par_ptr, lat_idx_1, lat_idx_2);

	// get the energies in the new neighborhood
	std::cout << "Post swap computation." << std::endl;
	this->neighbor_energetics (lat_idx_1, &(this->rotation_container.c1_final), &(this->rotation_container.E1_final));
	this->neighbor_energetics (lat_idx_2, &(this->rotation_container.c2_final), &(this->rotation_container.E2_final)); 
	this->selected_pair_interaction(this->Lattice[lat_idx_1], this->Lattice[lat_idx_2], &(this->rotation_container.cpair_final), &(this->rotation_container.Epair_final));

	return;
}

//////////////////////////////////////////////////////////
// modular functions for sampling based on orientations
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// rotation an end of a polymer
//////////////////////////////////////////////////////////

void Simulation::perturb_tail_rotation(int p_idx){

	// get the neighborlist of the particle at index 1
	std::array<int,3> loc_0 = this->Polymers[p_idx].chain[0]->coords;
	std::array<int,3> loc_1 = this->Polymers[p_idx].chain[1]->coords;

	// generate a neighbor list 
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

	// instantiate some containers
	std::array <double,8> cs_i           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cm_i           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cpair_i        = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cs_f           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cm_f           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cpair_f        = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> final_contacts = {0, 0, 0, 0, 0, 0, 0, 0};

	// get the choice for a neighbor
	int choice = rng_uniform(0, 25);

	// define all the local variables you need to perform valid tail rotation
	double E_solvent_i = 0; // energy of the solvent particle at the initial position
	double E_solvent_f = 0; // energy of the solvent particle at the final position
	double E_monomer_i = 0; // energy of the monomer particle at the initial position
	double E_monomer_f = 0; // energy of the monomer particle at the final position
	double E_pair_i    = 0; // initial pair energy
	double E_pair_f    = 0; // final pair energy
	double E_final     = 0; // energy of the final configuration of the system

	if (this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->ptype[0] == 's'){

		// find the energetic interaction for the solvent molecule 
		this->neighbor_energetics(lattice_index(ne_list[choice], this->y, this->z), &cs_i, &E_solvent_i);
		this->neighbor_energetics(lattice_index(loc_0,           this->y, this->z), &cm_i, &E_monomer_i);
		this->selected_pair_interaction(this->Polymers[p_idx].chain[0], this->Lattice[lattice_index(ne_list[choice], this->y, this->z)], &cpair_i, &E_pair_i);

		// setting up the switch 
		this->Polymers[p_idx].chain[0]->coords                                   = ne_list[choice];
		this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->coords = loc_0;
		
		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		this->Lattice[lattice_index(loc_0, this->y, this->z) ]           = this->Lattice[lattice_index(ne_list[choice], this->y, this->z)];
		this->Lattice[lattice_index(ne_list[choice], this->y, this->z)]  = this->Polymers[p_idx].chain[0];

		// get new neighbor energetics 
		this->neighbor_energetics(lattice_index(loc_0,                                  this->y, this->z), &cs_f, &E_solvent_f);
		this->neighbor_energetics(lattice_index(this->Polymers[p_idx].chain[0]->coords, this->y, this->z), &cm_f, &E_monomer_f);
		this->selected_pair_interaction(this->Polymers[p_idx].chain[0], this->Lattice[lattice_index(loc_0, this->y, this->z)], &cpair_f, &E_pair_f);

		E_final = this->sysEnergy-(E_solvent_i+E_monomer_i-E_pair_i)+(E_solvent_f+E_monomer_f-E_pair_f);
		final_contacts = add_arrays(subtract_arrays(this->contacts, add_arrays (cs_i, cm_i)), add_arrays(cs_f, cm_f));
		final_contacts = add_arrays(final_contacts, cpair_i);
		final_contacts = subtract_arrays(final_contacts, cpair_f);

		if (metropolis_acceptance(this->sysEnergy, E_final, this->T)) {
			this->sysEnergy = E_final;
			this->contacts  = final_contacts;
			this->IMP_BOOL  = true;
		}
		else {
			// revert the polymer 
			this->Polymers[p_idx].chain[0]->coords = loc_0;
			this->Lattice[lattice_index(loc_0, this->y, this->z)]->coords = ne_list[choice];
			
			// do the switch 
			this->Lattice[lattice_index(ne_list[choice], this->y, this->z) ] = this->Lattice[lattice_index(loc_0, this->y, this->z)];
			this->Lattice[lattice_index(loc_0, this->y, this->z) ]           = this->Polymers[p_idx].chain[0];
			this->IMP_BOOL = false;
		}
	}
	else {
		this->IMP_BOOL = false;
	}

	return;
}

void Simulation::perturb_head_rotation(int p_idx){

	// get the neighborlist of the particle at index-2
	int deg_poly = this->Polymers[p_idx].deg_poly;
	std::array<int,3> loc_0 = this->Polymers[p_idx].chain[deg_poly-1]->coords;
	std::array<int,3> loc_1 = this->Polymers[p_idx].chain[deg_poly-2]->coords;

	// generate a neighbor list 
	std::array<std::array<int,3>, 26> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

	// 
	std::array <double,8> cs_i           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cm_i           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cpair_i        = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cs_f           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cm_f           = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> cpair_f        = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> final_contacts = {0, 0, 0, 0, 0, 0, 0, 0};

	// get the choice for a neighbor
	int choice = rng_uniform(0, 25);

	// define all the local variables you need to perform valid tail rotation
	double E_solvent_i = 0; // energy of the solvent particle at the initial position
	double E_solvent_f = 0; // energy of the solvent particle at the final position
	double E_monomer_i = 0; // energy of the monomer particle at the initial position
	double E_monomer_f = 0; // energy of the monomer particle at the final position
	double E_pair_i    = 0; // initial pair energy
	double E_pair_f    = 0; // final pair energy
	double E_final     = 0; // energy of the final configuration of the system

	if (this->Lattice[lattice_index (ne_list[choice], this->y, this->z)]->ptype[0] == 's'){

		// find the energetic interaction for the solvent molecule 
		this->neighbor_energetics(lattice_index(ne_list[choice], this->y, this->z), &cs_i, &E_solvent_i);
		this->neighbor_energetics(lattice_index(loc_0,           this->y, this->z), &cm_i, &E_monomer_i);
		this->selected_pair_interaction(this->Polymers[p_idx].chain[deg_poly-1], this->Lattice[lattice_index(ne_list[choice], this->y, this->z)], &cpair_i, &E_pair_i);

		// setting up the switch 
		this->Polymers[p_idx].chain[deg_poly-1]->coords                           = ne_list[choice];
		this->Lattice[ lattice_index (ne_list[choice], this->y, this->z)]->coords = loc_0; 
		
		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		this->Lattice[ lattice_index (loc_0, this->y, this->z) ]           = this->Lattice[ lattice_index (ne_list[choice], this->y, this->z)];
		this->Lattice[ lattice_index (ne_list[choice], this->y, this->z)]  = this->Polymers[p_idx].chain[deg_poly-1];

		// get new neighbor energetics 
		this->neighbor_energetics(lattice_index(loc_0,                                           this->y, this->z), &cs_f, &E_solvent_f);
		this->neighbor_energetics(lattice_index(this->Polymers[p_idx].chain[deg_poly-1]->coords, this->y, this->z), &cm_f, &E_monomer_f);
		this->selected_pair_interaction(this->Polymers[p_idx].chain[deg_poly-1], this->Lattice[lattice_index(loc_0, this->y, this->z)], &cpair_f, &E_pair_f);
	
		// get the new energies
		E_final = this->sysEnergy-(E_solvent_i+E_monomer_i-E_pair_i)+(E_solvent_f+E_monomer_f-E_pair_f);
		final_contacts = add_arrays     (subtract_arrays(this->contacts, add_arrays(cs_i, cm_i)), add_arrays(cs_f, cm_f));
		final_contacts = add_arrays     (final_contacts, cpair_i);
		final_contacts = subtract_arrays(final_contacts, cpair_f);

		if (metropolis_acceptance(this->sysEnergy, E_final, this->T)) {
			this->sysEnergy = E_final;
			this->contacts  = final_contacts;
			this->IMP_BOOL  = true;
		}
		else {
			// revert the polymer 
			this->Polymers[p_idx].chain[deg_poly-1]->coords = loc_0;
			this->Lattice[lattice_index(loc_0, this->y, this->z)]->coords = ne_list[choice];
			
			// do the switch 
			this->Lattice[lattice_index(ne_list[choice], this->y, this->z)] = this->Lattice[lattice_index(loc_0, this->y, this->z)];
			this->Lattice[lattice_index(loc_0, this->y, this->z)]           = this->Polymers[p_idx].chain[deg_poly-1];
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

void Simulation::forward_reptation_with_tail_biting(std::array<double,8>* contacts_i, double* E_i, int p_idx){

	int deg_poly   = this->Polymers[p_idx].deg_poly;
	double Em1_i   = 0;
	double Em2_i   = 0;
	double Em1_f   = 0;
	double Em2_f   = 0;
	double Epair_i = 0;
	double Epair_f = 0;
	double Esys    = *E_i;

	std::array <double,8> cm1_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm1_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_f   = {0,0,0,0,0,0,0,0};
	std::array <int,3>    loc_m1    = {0,0,0};
	std::array <int,3>    loc_m2    = {0,0,0};
	std::array <double,8> contacts  = *contacts_i;

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
		// this->Lattice[lattice_index(loc_m1, this->y, this->z)] = this->Lattice[lattice_index(loc_m2, this->y, this->z)];
		// this->Lattice[lattice_index(loc_m2, this->y, this->z)] = this->Polymers[p_idx].chain[idx];
		// this->Lattice[lattice_index(loc_m1, this->y, this->z)]->coords = loc_m1;
		// this->Lattice[lattice_index(loc_m2, this->y, this->z)]->coords = loc_m2;

		// reevaluate energies
		this->neighbor_energetics(lattice_index(loc_m1, this->y, this->z), &cm1_f, &Em1_f);
		this->neighbor_energetics(lattice_index(loc_m2, this->y, this->z), &cm2_f, &Em2_f);
		this->selected_pair_interaction(this->Polymers[p_idx].chain[idx], this->Polymers[p_idx].chain[idx+1], &cpair_f, &Epair_f);

		// update energies and contacts
		contacts = add_arrays(subtract_arrays(contacts, add_arrays(cm1_i, cm2_i)), add_arrays(cm1_f, cm2_f));
		contacts = add_arrays(contacts, cpair_i);
		contacts = subtract_arrays(contacts, cpair_f);
		Esys = Esys - (Em1_i+Em2_i-Epair_i) + (Em1_f+Em2_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		// resetting everything
		cpair_i  = {0, 0, 0, 0, 0, 0, 0, 0};
		cpair_f  = {0, 0, 0, 0, 0, 0, 0, 0};
		cm1_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm1_f    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm2_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm2_f    = {0, 0, 0, 0, 0, 0, 0, 0};
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

void Simulation::forward_reptation_without_tail_biting(std::array<double,8>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx){

	int    deg_poly = this->Polymers[p_idx].deg_poly;
	double Em_i     = 0;
	double Es_i     = 0;
	double Em_f     = 0;
	double Es_f     = 0;
	double Epair_i  = 0;
	double Epair_f  = 0;
	double Esys     = *E_i;

	std::array <double,8> cm_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i  = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_f  = {0,0,0,0,0,0,0,0};
	std::array <int,3>    loc_m    = {0,0,0};
	std::array <int,3>    loc_s    = *to_slither;
	std::array <double,8> contacts = *contacts_i;

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
		contacts = add_arrays(subtract_arrays(contacts, add_arrays(cm_i, cs_i)), add_arrays(cm_f, cs_f));
		contacts = add_arrays(contacts, cpair_i);
		contacts = subtract_arrays(contacts, cpair_f);

		Esys = Esys - (Em_i+Es_i-Epair_i) + (Em_f+Es_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		// the particle is now at the old position of the monomer bead
		*to_slither = loc_m; 

		// resetting everything
		cpair_i = {0, 0, 0, 0, 0, 0, 0, 0};
		cpair_f = {0, 0, 0, 0, 0, 0, 0, 0};
		cs_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cs_f    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm_f    = {0, 0, 0, 0, 0, 0, 0, 0};
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

void Simulation::backward_reptation_with_head_butting(std::array<double,8>* contacts_i, double* E_i, int p_idx){

	int    deg_poly    = this->Polymers[p_idx].deg_poly;
	double Em1_i       = 0;
	double Em2_i       = 0;
	double Em1_f       = 0;
	double Em2_f       = 0;
	double Epair_i     = 0;
	double Epair_f     = 0;
	double Esys        = *E_i;

	std::array <double,8> cm1_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm1_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_f   = {0,0,0,0,0,0,0,0};
	std::array <int,3>    loc_m1    = {0,0,0}; 
	std::array <int,3>    loc_m2    = {0,0,0}; 
	std::array <double,8> contacts  = *contacts_i;

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

		contacts = add_arrays(subtract_arrays(contacts, add_arrays(cm1_i, cm2_i)), add_arrays(cm1_f, cm2_f));
		contacts = add_arrays(contacts, cpair_i);
		contacts = subtract_arrays(contacts, cpair_f);
		Esys = Esys-(Em1_i+Em2_i-Epair_i)+(Em1_f+Em2_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		// resetting everything
		cpair_i  = {0, 0, 0, 0, 0, 0, 0, 0};
		cpair_f  = {0, 0, 0, 0, 0, 0, 0, 0};
		cm1_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm1_f    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm2_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm2_f    = {0, 0, 0, 0, 0, 0, 0, 0};
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

void Simulation::backward_reptation_without_head_butting(std::array<double,8>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx){

	int    deg_poly    = this->Polymers[p_idx].deg_poly;
	double Em_i        = 0;
	double Es_i        = 0; 
	double Em_f        = 0; 
	double Es_f        = 0; 
	double Epair_i     = 0;
	double Epair_f     = 0;
	double Esys        = *E_i; 

	std::array <double,8> cm_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i  = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_f     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_f  = {0,0,0,0,0,0,0,0};
	std::array <int,3>    loc_m    = {0,0,0};
	std::array <int,3>    loc_s    = *to_slither;
	std::array <double,8> contacts = *contacts_i;

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

		contacts = add_arrays(subtract_arrays(contacts, add_arrays(cm_i, cs_i)), add_arrays(cm_f, cs_f));
		contacts = add_arrays(contacts, cpair_i);
		contacts = subtract_arrays(contacts, cpair_f);		
		Esys = Esys-(Em_i+Es_i-Epair_i)+(Em_f+Es_f-Epair_f);

		// std::cout << "Run the debugger..." << std::endl;
		// this->Polymers[p_idx].print_chain();
		// this->debug_checks_energy_contacts(Esys, contacts);

		*to_slither = loc_m; 

		// resetting everything
		cpair_i = {0, 0, 0, 0, 0, 0, 0, 0};
		cpair_f = {0, 0, 0, 0, 0, 0, 0, 0};
		cs_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cs_f    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm_i    = {0, 0, 0, 0, 0, 0, 0, 0};
		cm_f    = {0, 0, 0, 0, 0, 0, 0, 0};
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

	// instantiate deg_poly because this will be used later
	int                  deg_poly   = this->Polymers[p_idx].deg_poly;
	
	// instantiate a copy of current energy
	double               E_i        = this->sysEnergy;

	// instantiate the locations of the first and final monomer
	std::array<int,3>    loc_i      = this->Polymers[p_idx].chain[0]->coords;
	std::array<int,3>    loc_f      = this->Polymers[p_idx].chain[deg_poly-1]->coords;

	// instantiate a copy of current contacts
	std::array<double,8> contacts_i = this->contacts; 

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

	// instantiate deg_poly because this will be used later
	int                   deg_poly         = this->Polymers[p_idx].deg_poly;

	// instantiate a copy of current energy
	double                E_i              = this->sysEnergy;

	// instantiate the locations of the first and final monomer
	std::array <int,3>    loc_i            = this->Polymers[p_idx].chain[0]->coords;
	std::array <int,3>    loc_f            = this->Polymers[p_idx].chain[deg_poly-1]->coords;

	// instantiate a copy of current contacts
	std::array <double,8> contacts_i       = this->contacts;


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

	// instantiate 
	std::vector<int> old_ori; 
	std::vector<int> new_ori; 

	std::array<double,5>               energies           = {0,0,0,0,0};
	std::array<double,5>               boltzmann          = {0,0,0,0,0};
	std::array<int,5>                  orientations       = {0,0,0,0,0};

	// define contact stores
	std::array<std::array<double,8>,5> contacts_store     = {this->contacts, this->contacts, this->contacts, this->contacts, this->contacts};
	std::array<double,8>               contacts_sys		  = this->contacts;
	std::array<double,8>               contacts_i         = {0,0,0,0,0,0,0,0};
	std::array<double,8>               contacts_pert      = {0,0,0,0,0,0,0,0};
	std::array<double,8>               frontflow_contacts = {0,0,0,0,0,0,0,0};

	// these are running sums and cumulant probabilities to keep track of
	double                             rboltzmann         = 0;
	double                             frontflow_energy   = 0;
	double                             prob_o_to_n        = 1;
	double                             prob_n_to_o        = 1;
	double                             Emin               = 0;

	// instantiate some newer variables
	double rng     = 0;
	double rng_acc = 0;
	double rsum    = 0;
	double Ei      = 0;
	double Epert   = 0;
	double Esys    = this->sysEnergy;

	// 
	int e_idx         = 0;
	int m_lattice_idx = 0; 

	// loop over all solvent_indices 
	for ( int i{0}; i < nflip; ++i ){

		rboltzmann        = 0; 
		old_ori.push_back ((this->Polymers)[p_idx].chain[ polymer_indices[i] ]->orientation); 
		m_lattice_idx     = lattice_index(this->Polymers[p_idx].chain[polymer_indices[i]]->coords, y, z);
		this->neighbor_energetics (m_lattice_idx, &contacts_i, &Ei); 

		for (int j{0}; j < ntest; ++j){
			this->Polymers[p_idx].chain[polymer_indices[i]]->orientation = rng_uniform (0, 25); 
			orientations[j]   = this->Polymers[p_idx].chain[polymer_indices[i]]->orientation;
			this->neighbor_energetics(m_lattice_idx, &contacts_pert, &Epert);
			energies[j]       = Esys - Ei + Epert; 
			contacts_store[j] = subtract_arrays (&contacts_sys, &contacts_i); 
			contacts_store[j] = add_arrays      (&contacts_store[j], &contacts_pert);
			contacts_pert     = {0, 0, 0, 0, 0, 0, 0, 0};
			Epert             = 0;
		}

		// get the energy minima
		Emin = *std::min_element(energies.begin(), energies.end());

		for ( int k{0}; k < ntest; ++k ){
			boltzmann[k] = std::exp (-1/this->T*(energies[k] - Emin));
			rboltzmann  += boltzmann [k]; 
		}

		rng   = rng_uniform (0.0, 1.0); 
		rsum  = 0; 
		e_idx = 0;

		for ( int j{0}; j<5; ++j){
			rsum += boltzmann[j]/rboltzmann;
			if ( rng < rsum ){
				e_idx = j; 
				break; 
			}
		}

		// make the jump to the new state 
		new_ori.push_back (orientations[e_idx]);
		this->Polymers[p_idx].chain[polymer_indices[i]]->orientation = orientations[e_idx];
		prob_o_to_n *= boltzmann[e_idx]/rboltzmann;
		Esys         = energies[e_idx];
		contacts_sys = contacts_store [e_idx];
		
		// resetting... 
		contacts_i = {0, 0, 0, 0, 0, 0, 0, 0};
		Ei         = 0;
	}

	frontflow_energy   = energies[e_idx];
	frontflow_contacts = contacts_store[e_idx];

	contacts_i    = {0, 0, 0, 0, 0, 0, 0, 0};
	contacts_pert = {0, 0, 0, 0, 0, 0, 0, 0};
	Ei            = 0;
	Epert         = 0;

	for (int i{0}; i < nflip; ++i){

		rboltzmann = 0; 
		m_lattice_idx = lattice_index((this->Polymers)[p_idx].chain[polymer_indices[i]]->coords, y, z);
		this->neighbor_energetics(m_lattice_idx, &contacts_i, &Ei);

		(this->Polymers)[p_idx].chain[polymer_indices[i]]->orientation = old_ori[i];
		this->neighbor_energetics(m_lattice_idx, &contacts_pert, &Epert);
		energies[0] = Esys - Ei + Epert;
		contacts_store[0]  = subtract_arrays(&contacts_sys, &contacts_i);
		contacts_store[0]  = add_arrays     (&contacts_store[0], &contacts_pert);

		contacts_pert = {0, 0, 0, 0, 0, 0, 0, 0};
		Epert         = 0;

		for (int j{1}; j<ntest; ++j){
			this->Polymers[p_idx].chain[polymer_indices[i]]->orientation = rng_uniform (0, 25); 
			this->neighbor_energetics(m_lattice_idx, &contacts_pert, &Epert);
			energies[j] = Esys - Ei + Epert;
			contacts_store [j]  = subtract_arrays (&contacts_sys, &contacts_i);
			contacts_store [j]  = add_arrays      (&contacts_store[j], &contacts_pert);
			contacts_pert = {0, 0, 0, 0, 0, 0, 0, 0};
			Epert         = 0;
		}

		Emin = *std::min_element(energies.begin(), energies.end());

		for (int k{0}; k < ntest; ++k){
			boltzmann [k] = std::exp (-1/this->T*(energies[k] - Emin));
			rboltzmann   += boltzmann[k];
		}
		prob_n_to_o      *= boltzmann[0]/rboltzmann;

		// make the jump to the old state 
		this->Polymers[p_idx].chain[polymer_indices[i]]->orientation = old_ori[i];
		Esys = energies[0];
		contacts_sys = contacts_store[0];

		contacts_i = {0, 0, 0, 0, 0, 0, 0, 0};
		Ei         = 0;
	}

	rng_acc = rng_uniform(0.0, 1.0);

	if ( rng_acc < std::exp(-1/this->T * (frontflow_energy-this->sysEnergy)) * prob_n_to_o/prob_o_to_n) {

		// if accepted, return to the new orientations
		for (int j{0}; j < nflip; ++j){
			(this->Polymers)[p_idx].chain[polymer_indices[j]]->orientation = new_ori[j];
		}

		this->sysEnergy = frontflow_energy;
		this->contacts  = frontflow_contacts;
		this->IMP_BOOL  = true;

	}
	else {
		this->IMP_BOOL = false; 
	}

	return;

}

void Simulation::perturb_solvation_shell_flip(){

	std::set   <int> solvation_shell_set;
	std::array <std::array<int,3>,26> ne_list;

	// get the first solvation shell 
	for ( Polymer& pmer: this->Polymers){
		for ( Particle*& p: pmer.chain ){
		ne_list = obtain_ne_list ( p->coords, this->x, this->y, this->z); 
			for ( std::array <int,3>& loc: ne_list ){
				if ( this->Lattice[lattice_index (loc, y, z)]->ptype[0] == 's' ){
					solvation_shell_set.insert (lattice_index (loc, y, z)); 
				}
			}
		}
	}

	std::vector <int> solvation_shell_indices (solvation_shell_set.begin(), solvation_shell_set.end()); 
	int nflip = rng_uniform(1, static_cast<int>(solvation_shell_indices.size()/2 ) );     // number of sites to be flipped 

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));
    
	std::vector <int> old_ori; // vector to hold old orientations
	std::vector <int> new_ori; // vector to hold new orientations

	// energy store for boltzmann sampling 
	// instantiating a bunch of variables for the process 
	int                                 ntest              = 5; 
	std::array <double,5>               energies           = {0,0,0,0,0}; 
	std::array <double,5>               boltzmann          = {0,0,0,0,0};
	std::array <int,5>                  orientations       = {0,0,0,0,0}; 

	double                              Esys    = this->sysEnergy; 
	double                              Ei      = 0; 
	double                              Epert   = 0; 
	double                              Emin    = 0; 

	double                              rboltzmann         = 0;  
	double                              frontflow_energy   = 0; 
	double                              prob_n_to_o        = 1; 
	double                              prob_o_to_n        = 1; 

	std::array<std::array<double,8>,5>  contacts_store     = {this->contacts, this->contacts, this->contacts, this->contacts, this->contacts}; 
	std::array <double,8>               contacts_sys       = this->contacts; 

	std::array <double,8>               contacts_i         = {0,0,0,0,0,0,0,0}; 
	std::array <double,8>               contacts_pert      = {0,0,0,0,0,0,0,0}; 
	std::array <double,8>               frontflow_contacts = {0,0,0,0,0,0,0,0};

	double rng     = 0;
	double rng_acc = 0;
	double rsum    = 0;
	int    e_idx   = 0;

	for (int i{0}; i < nflip; ++i){
		// get the flip index 
		rboltzmann = 0; 
		old_ori.push_back((this->Lattice)[solvation_shell_indices[i]]->orientation);

		// find the neighboring interaction energies 
		this->neighbor_energetics(solvation_shell_indices[i], &contacts_i, &Ei); 

		for (int j{0}; j < ntest; ++j){
			(this->Lattice)[solvation_shell_indices[i]]->orientation = rng_uniform (0, 25); 
			orientations [j]    = (this->Lattice)[solvation_shell_indices[i]]->orientation; 
			this->neighbor_energetics(solvation_shell_indices[i], &contacts_pert, &Epert);
			energies [j]        = Esys - Ei + Epert; 
			contacts_store [j]  = subtract_arrays(&contacts_sys, &contacts_i); 
			contacts_store [j]  = add_arrays(&contacts_store[j], &contacts_pert); 
			contacts_pert       = {0, 0, 0, 0, 0, 0, 0, 0};
			Epert               = 0;
		}

		// get the minimal energy
		Emin = *std::min_element ( energies.begin(), energies.end() ); 

		// find the running boltzmann sum
		for (int k{0}; k < ntest; ++k){
			boltzmann[k] = std::exp (-1/this->T*(energies[k] - Emin)); 
			rboltzmann  += boltzmann[k]; 
		}

		// instantiate some variables
		rng     = rng_uniform (0.0, 1.0); 
		rsum    = 0; 
		e_idx   = 0; 

		for (int j{0}; j < ntest; ++j){
			rsum += boltzmann[j]/rboltzmann; 
			if ( rng < rsum ) {
				e_idx = j; 
				break; 
			}
		}

		// make the jump to the new state 
		new_ori.push_back (orientations[e_idx]); 
		this->Lattice[solvation_shell_indices[i]]->orientation = orientations[e_idx]; 
		prob_o_to_n *= boltzmann[e_idx]/rboltzmann;
		Esys         = energies[e_idx];
		contacts_sys = contacts_store[e_idx];

		// reset 
		contacts_i = {0, 0, 0, 0, 0, 0, 0, 0};
		Ei         = 0;
	}

	frontflow_energy   = energies[e_idx];
	frontflow_contacts = contacts_store[e_idx];

	// reset...
	contacts_i    = {0, 0, 0, 0, 0, 0, 0, 0};
	contacts_pert = {0, 0, 0, 0, 0, 0, 0, 0};
	Ei            = 0;
	Epert         = 0;

	for (int i{0}; i<nflip; ++i){

		this->neighbor_energetics(solvation_shell_indices[i], &contacts_i, &Ei);
		rboltzmann = 0;

		// first iteration
		this->Lattice[solvation_shell_indices[i]]->orientation = old_ori[i];
		this->neighbor_energetics(solvation_shell_indices[i], &contacts_pert, &Epert);
		energies[0] = Esys - Ei + Epert;
		contacts_store[0] = subtract_arrays(&contacts_sys, &contacts_i);
		contacts_store[0] = add_arrays(&contacts_store[0], &contacts_pert);

		// reset...
		contacts_pert = {0, 0, 0, 0, 0, 0, 0, 0}; 
		Epert         = 0;

		for (int j{1}; j <ntest; ++j){

			// solvation shell flips
			this->Lattice[solvation_shell_indices[i]]->orientation = rng_uniform(0, 25);
			this->neighbor_energetics(solvation_shell_indices[i], &contacts_pert, &Epert); 
			energies      [j] = Esys - Ei + Epert;
			contacts_store[j] = subtract_arrays(&contacts_sys, &contacts_i);
			contacts_store[j] = add_arrays(&contacts_store[j], &contacts_pert);

			// reset...
			contacts_pert     = {0, 0, 0, 0, 0, 0, 0, 0};
			Epert             = 0;
		}

		// get the minimal energy of the system
		Emin = *std::min_element(energies.begin(), energies.end()); 

		for (int k{0}; k<ntest; ++k){
			boltzmann[k] = std::exp(-1/this->T*(energies[k]-Emin));
			rboltzmann  += boltzmann[k]; 
		}

		prob_n_to_o += boltzmann[0]/rboltzmann; 
		this->Lattice[solvation_shell_indices[i]]->orientation = old_ori[i]; 
		Esys         = energies[0];
		contacts_sys = contacts_store[0]; 

		// reset...
		contacts_i   = {0,0,0,0,0,0,0,0};
		Ei           = 0;

	}

	rng_acc = rng_uniform(0.0, 1.0);
	if (rng_acc < std::exp (-1/this->T * (frontflow_energy - this->sysEnergy)) * prob_n_to_o/prob_o_to_n){

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

	return;
}

void Simulation::perturb_lattice_flip(){
	
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

void Simulation::perturb_solvent_exchange_from_shell(){

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

	std::array <double,8> cs1_i              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_i              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i            = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs1_f              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_f              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_f            = {0,0,0,0,0,0,0,0};
	std::array <double,8> frontflow_contacts = this->contacts;

	double Es1_i            = 0;
	double Es1_f            = 0;
	double Es2_i            = 0;
	double Es2_f            = 0;
	double Epair_i          = 0;
	double Epair_f          = 0;
	double frontflow_energy = 0;
	double Esys             = this->sysEnergy;

	Particle* tmp_par_ptr {nullptr};

	int exc_idx    = 0; 
	int nswitch    = rng_uniform(1, static_cast<int>(solvation_shell_indices.size()/4));
	double rng_acc = 0.0;

	for (int n{0}; n<nswitch; ++n){
		
		exc_idx = rng_uniform(0, this->x * this->y * this->z - 1);
		if (exc_idx == solvation_shell_indices[n] || this->Lattice[exc_idx]->ptype[0] == 'm'){
			continue;
		}
		
		// get the energies at the initial stage
		this->neighbor_energetics(solvation_shell_indices[n], &cs1_i, &Es1_i);
		this->neighbor_energetics(exc_idx, &cs2_i, &Es2_i);
		this->selected_pair_interaction(this->Lattice[solvation_shell_indices[n]], this->Lattice[exc_idx], &cpair_i, &Epair_i); 

		// swap particles
		tmp_par_ptr = (this->Lattice)[ solvation_shell_indices[n]];
		
		(this->Lattice)[(solvation_shell_indices)[n]] = (this->Lattice)[exc_idx];
		(this->Lattice)[(solvation_shell_indices)[n]]->coords = location ((solvation_shell_indices)[n], this->x, this->y, this->z);

		(this->Lattice)[exc_idx] = tmp_par_ptr;
		(this->Lattice)[exc_idx]->coords = location (exc_idx, this->x, this->y, this->z);

		// get the new energetics
		this->neighbor_energetics (solvation_shell_indices[n], &cs1_f, &Es1_f);
		this->neighbor_energetics (exc_idx, &cs2_f, &Es2_f); 
		this->selected_pair_interaction(this->Lattice[solvation_shell_indices[n]], this->Lattice[exc_idx], &cpair_f, &Epair_f); 

		// update contacts and energetics
		frontflow_energy   = Esys - (Es1_i + Es2_i -Epair_i) + (Es1_f + Es2_f - Epair_f);  
		frontflow_contacts = add_arrays(subtract_arrays(this->contacts, add_arrays (cs1_i, cs2_i)), add_arrays (cs1_f, cs2_f));
		frontflow_contacts = add_arrays(frontflow_contacts, cpair_i);
		frontflow_contacts = subtract_arrays(frontflow_contacts, cpair_f);

		// get some random number for metropolis
		rng_acc = rng_uniform (0.0, 1.0);

		if (rng_acc < std::exp (-1/this->T * (frontflow_energy - this->sysEnergy))){
			this->sysEnergy = frontflow_energy;
			this->contacts  = frontflow_contacts;
			this->IMP_BOOL  = true; 
		}
		else {
			tmp_par_ptr = (this->Lattice)[solvation_shell_indices[n]];
			(this->Lattice)[solvation_shell_indices [n]]         = (this->Lattice)[exc_idx]; 
			(this->Lattice)[solvation_shell_indices [n]]->coords = location (solvation_shell_indices[n], this->x, this->y, this->z);

			(this->Lattice)[exc_idx]         = tmp_par_ptr;
			(this->Lattice)[exc_idx]->coords = location (exc_idx, this->x, this->y, this->z);
			this->IMP_BOOL = false; 
		}

		// resetting 
		cs1_i   = {0, 0, 0, 0, 0, 0, 0, 0};
		cs1_f   = {0, 0, 0, 0, 0, 0, 0, 0};
		cs2_i   = {0, 0, 0, 0, 0, 0, 0, 0};
		cs2_f   = {0, 0, 0, 0, 0, 0, 0, 0};
		cpair_i = {0, 0, 0, 0, 0, 0, 0, 0};
		cpair_f = {0, 0, 0, 0, 0, 0, 0, 0};
		Es1_i   = 0;
		Es1_f   = 0; 
		Es2_i   = 0;
		Es2_f   = 0;
		Epair_i = 0; 
		Epair_f = 0;

	}

	return;
}

void Simulation::perturb_solvent_exchange(){

	int idx1 = -1;
	int idx2 = -1; 

	// instantiate some containers 
	std::array <double,8> frontflow_contacts = this->contacts;
	std::array <double,8> sys_contacts       = this->contacts;
	std::array <double,8> cs1_i              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_i              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i            = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs1_f              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_f              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_f            = {0,0,0,0,0,0,0,0};

	// instantiate some doubles
	double Esys             = this->sysEnergy;
	double frontflow_energy = 0;
	double Es1_i            = 0;
	double Es1_f            = 0;
	double Es2_i            = 0;
	double Es2_f            = 0;
	double Epair_i          = 0;
	double Epair_f          = 0;
	double rng_acc          = 0;

	Particle* tmp_par_ptr {nullptr};

	// arbitrary
	// TO DO: tune aggression
	int nswitches = 50; 

	for (int j{0}; j<nswitches; ++j){

		idx1 = rng_uniform(0, this->x * this->y * this->z - 1);
		idx2 = rng_uniform(0, this->x * this->y * this->z - 1);

		if (this->Lattice[idx1]->ptype[0] == 'm' || this->Lattice[idx2]->ptype[0] == 'm' || idx1 == idx2){
			continue; 
		}

		// get the initial energetics
		this->neighbor_energetics(idx1, &cs1_i, &Es1_i);
		this->neighbor_energetics(idx2, &cs2_i, &Es2_i);
		this->selected_pair_interaction((this->Lattice)[idx1], (this->Lattice)[idx2], &cpair_i, &Epair_i); 

		// swap the particles
		tmp_par_ptr = (this->Lattice)[idx1];
		
		(this->Lattice)[idx1] = (this->Lattice)[idx2];
		(this->Lattice)[idx1]->coords = location (idx1, this->x, this->y, this->z);

		(this->Lattice)[idx2] = tmp_par_ptr;
		(this->Lattice)[idx2]->coords = location(idx2, this->x, this->y, this->z);

		// flip the particle newly added to the solvation shell
		this->neighbor_energetics (idx1, &cs1_f, &Es1_f);
		this->neighbor_energetics (idx2, &cs2_f, &Es2_f); 
		this->selected_pair_interaction(this->Lattice[idx1], this->Lattice[idx2], &cpair_f, &Epair_f); 

		// get the new energies and contacts
		frontflow_energy   = Esys - (Es1_i + Es2_i -Epair_i) + (Es1_f + Es2_f - Epair_f);
		frontflow_contacts = add_arrays(subtract_arrays(sys_contacts, add_arrays (cs1_i, cs2_i)), add_arrays (cs1_f, cs2_f));
		frontflow_contacts = add_arrays(frontflow_contacts, cpair_i);
		frontflow_contacts = subtract_arrays(frontflow_contacts, cpair_f);

		rng_acc = rng_uniform (0.0, 1.0); 
		if (rng_acc < std::exp (-1/this->T * (frontflow_energy - Esys))){
			Esys         = frontflow_energy;
			sys_contacts = frontflow_contacts;
	
		}
		else {
			tmp_par_ptr = (this->Lattice)[idx1];
			(this->Lattice)[idx1]         = (this->Lattice)[idx2]; 
			(this->Lattice)[idx1]->coords = location(idx1, this->x, this->y, this->z);

			(this->Lattice)[idx2]         = tmp_par_ptr;
			(this->Lattice)[idx2]->coords = location (idx2, this->x, this->y, this->z);
	
		}

		cpair_i = {0, 0, 0, 0, 0, 0, 0, 0};
		cpair_f = {0, 0, 0, 0, 0, 0, 0, 0};
		cs1_i   = {0, 0, 0, 0, 0, 0, 0, 0};
		cs1_f   = {0, 0, 0, 0, 0, 0, 0, 0};
		cs2_i   = {0, 0, 0, 0, 0, 0, 0, 0};
		cs2_f   = {0, 0, 0, 0, 0, 0, 0, 0};
		Es1_i   = 0;
		Es2_i   = 0;
		Es1_f   = 0;
		Es2_f   = 0;
		Epair_i = 0;
		Epair_f = 0;
	}

	this->sysEnergy = Esys; 
	this->contacts  = sys_contacts;

	return; 
}

//////////////////////////////////////////////////////////
// perform regrowth
//////////////////////////////////////////////////////////

void Simulation::perturb_regrowth(int p_idx){

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

		for (int i{m_idx+1}; i<deg_poly; ++i){
			old_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		this->forward_head_regrowth(&forw_contacts, &prob_o_to_n, &forw_energy, p_idx, m_idx);

		for (int i{m_idx+1}; i<deg_poly; ++i){
			new_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		if (old_cut == new_cut){
			// do nothing
			;
		}

		else if (!(this->IMP_BOOL)){
			this->accept_after_head_regrowth(&new_cut, &old_cut);
		}

		else {
			back_contacts = forw_contacts;
			back_energy   = forw_energy;
			this->backward_head_regrowth(&old_cut, &back_contacts, &prob_n_to_o, &back_energy, p_idx, m_idx, recursion_depth);

			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * prob_n_to_o / prob_o_to_n) {
				this->accept_after_head_regrowth(&old_cut, &new_cut);
				this->sysEnergy = forw_energy;
				this->contacts  = forw_contacts;
				this->IMP_BOOL  = true;
			}
			else {
				this->IMP_BOOL = false;
			}
		}
	}

	else {

		for (int i{0}; i<m_idx; ++i){
			old_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		this->forward_tail_regrowth(&forw_contacts, &prob_o_to_n, &forw_energy, p_idx, m_idx);

		for (int i{0}; i<m_idx; ++i){
			new_cut.push_back(this->Polymers.at(p_idx).chain.at(i)->coords);
		}

		if (old_cut == new_cut){
			// do nothing
			;
		}

		else if (!(this->IMP_BOOL)){
			this->accept_after_tail_regrowth(&new_cut, &old_cut);
		}

		else {

			back_contacts = forw_contacts; 
			back_energy = forw_energy;
			backward_tail_regrowth(&old_cut, &back_contacts, &prob_n_to_o, &back_energy, p_idx, m_idx, recursion_depth);

			rng_acc = rng_uniform(0.0, 1.0);
			if (rng_acc < std::exp(-1/this->T * (forw_energy - this->sysEnergy)) * prob_n_to_o/prob_o_to_n) {
				this->accept_after_tail_regrowth(&old_cut, &new_cut);
				this->sysEnergy = forw_energy;
				this->contacts  = forw_contacts;
			}

			else {
				this->IMP_BOOL = false;
			}
		}
	}

	return;
}

void Simulation::forward_head_regrowth(std::array <double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx){

	int deg_poly   = this->Polymers[p_idx].deg_poly; 
	if (m_idx == deg_poly-1) {
		this->IMP_BOOL = true;
		return; 
	}

	// set up containers for contacts
	std::array<double,8> current_contacts = *forw_contacts;
	std::array<std::array<double,8>,5> contacts_store;

	// set up contacts for energies and boltzmann factors
	std::array<double,5>               energies;
	std::array<double,5>               boltzmann;

	// set up a varable for the relevant monomer location
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

	Particle* tmp_par_ptr {nullptr};

	while (idx_counter<5) {
		if (ne_list[idx_counter] == loc_m){
			energies[idx_counter]        = Esys;
			contacts_store[idx_counter]  = current_contacts;
			block_counter               += 1;
		}
		else if (this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm') {

			for (int u{0}; u<deg_poly; ++u){
				if (this->Polymers[p_idx].chain[u]->coords == ne_list[idx_counter]){
					self_swap_idx = u;
					break;
				}
			}

			if (self_swap_idx <= m_idx){
				energies[idx_counter]       = 1e+8;
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1};
				block_counter              += 1;
			}

			else {

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

				// swap particles
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

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
			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

			// swap particles
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

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

	if (block_counter == 5) {
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
		this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[e_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->forward_head_regrowth(forw_contacts, prob_o_to_n, forw_energy, p_idx, m_idx+1);
	}

	return;
}

void Simulation::accept_after_head_regrowth(std::vector<std::array<int,3>>* old_cut, std::vector<std::array<int,3>>* new_cut){

	// define some stuff
	int L {0};

	// instantiate the linked list
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

void Simulation::backward_head_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth){

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

	// define containers for contacts
	std::array <double,8>               current_contacts = *back_contacts;
	std::array <std::array<double,8>,5> contacts_store; 

	// define contacts for energies and boltzmann factors
	std::array <double,5>               energies;
	std::array <double,5>               boltzmann;

	// define a variable for location of relevant monomer index
	std::array <int,3> loc_m = (this->Polymers)[p_idx].chain[m_idx+1]->coords; // this is the key item of interest

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((this->Polymers)[p_idx].chain[m_idx]->coords, this->x, this->y, this->z);

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 
	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[recursion_depth]) - ne_list.begin();

	// do a quick switch
	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	// containers for contacts
	std::array <double,8> cm_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int idx_counter         = 0;
	int self_swap_idx       = -1;

	// get the empty pointer
	Particle* tmp_par_ptr {nullptr};

	while (idx_counter < 5){

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
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

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
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));	

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
	this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[0], this->y, this->z), lattice_index(loc_m, this->y, this->z));
	this->backward_head_regrowth(old_cut, back_contacts, prob_n_to_o, back_energy, p_idx, m_idx+1, recursion_depth+1);

	return;
}

void Simulation::forward_tail_regrowth(std::array <double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx){

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

	Particle* tmp_par_ptr {nullptr};

	while(idx_counter < 5) {

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

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

				// swap particles
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

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

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

			// swap particles
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

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
		this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[e_idx],this->y, this->z), lattice_index(loc_m, this->y, this->z));
		this->forward_tail_regrowth(forw_contacts, prob_o_to_n, forw_energy, p_idx, m_idx-1);
	}


	return;
}

void Simulation::accept_after_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::vector <std::array<int,3>>* new_cut){

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

void Simulation::backward_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth){

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

	// containers for contacts
	std::array <double,8>               current_contacts = *back_contacts;
	std::array <std::array<double,8>,5> contacts_store;

	// stores for energy
	std::array <double,5>               energies;
	std::array <double,5>               boltzmann;

	// store for location
	std::array <int,3> loc_m = (this->Polymers)[p_idx].chain[m_idx-1]->coords; // this is the key item of interest

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((this->Polymers)[p_idx].chain[m_idx]->coords, this->x, this->y, this->z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[m_idx-1]) - ne_list.begin() ; 

	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	// define containers for contacts
	std::array <double,8> cm_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_i    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_i = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n    = {0,0,0,0,0,0,0,0};
	std::array <double,8> cpair_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int idx_counter         = 0;
	int self_swap_idx       = -1;

	Particle* tmp_par_ptr {nullptr};

	while (idx_counter < 5){

		if (ne_list[idx_counter] == loc_m){
			energies[idx_counter]       = *back_energy;
			contacts_store[idx_counter] = current_contacts;
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
			}
			else {
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_i, &Es_i);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_i, &Em_i);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_i, &Epair_i);

				// swap particles
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

				// get the energies
				this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
				this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
				this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

				// set up the final energies
				energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
				contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
				contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

				// revert back to original structure, reswap particles
				this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

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
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

			// get the energies
			this->neighbor_energetics(lattice_index(ne_list[idx_counter], this->y, this->z), &cs_n, &Es_n);
			this->neighbor_energetics(lattice_index(loc_m,                this->y, this->z), &cm_n, &Em_n);
			this->selected_pair_interaction(this->Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)], this->Lattice[lattice_index(loc_m, this->y, this->z)], &cpair_n, &Epair_n);

			// set up the final energies
			energies[idx_counter]       = Esys - (Es_i+Em_i-Epair_i) + (Es_n+Em_n-Epair_n);
			contacts_store[idx_counter] = add_arrays(subtract_arrays(current_contacts, add_arrays (cm_i, cs_i)), add_arrays (cm_n, cs_n));
			contacts_store[idx_counter] = subtract_arrays(add_arrays(contacts_store[idx_counter], cpair_i), cpair_n);

			// revert back to original structure, reswap particles
			this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[idx_counter],this->y, this->z), lattice_index(loc_m, this->y, this->z));

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
	this->perturb_particle_swap(tmp_par_ptr, lattice_index(ne_list[0], this->y, this->z), lattice_index(loc_m, this->y, this->z));
	this->backward_tail_regrowth(old_cut, back_contacts, prob_n_to_o, back_energy, p_idx, m_idx-1, recursion_depth+1);

	return; 
}

//////////////////////////////////////////////////////////
//   
//                Run the simulation
//
//////////////////////////////////////////////////////////

void Simulation::perturb_system_straight(){
	// std::cout << "Running it straight." << std::endl;
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

void Simulation::perturb_system_debug(){
	std::cout << "Running the debugging-based perturbation." << std::endl;
	int r = 8; // rng_uniform(0, 1);
	
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