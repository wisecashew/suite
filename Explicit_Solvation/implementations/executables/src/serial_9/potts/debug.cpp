#include "Potts.hpp"

// this is a very simplistic, readable pair interaction function
void Potts::debug_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& db_contacts, double& db_energy){

	double dot_product = 0;
	double theta_1     = 0;
	double theta_2     = 0;
	double magnitude   = 0;

	std::array <int,3> connvec = {0, 0, 0};

	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 

	std::string interaction_type = std::get<0>((this->InteractionMap)[particle_pair]);
	switch (interaction_type[0]){

		case 'i':
			(db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])+1] += 0.5;
			db_energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5;
			break;

		case 'p':
			dot_product = take_dot_product ( p1->orientation, p2->orientation );
			if (dot_product > 0.54){
				(db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])] += 0.5; 
				db_energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 0.5; 
				db_energy += std::get<2>((this->InteractionMap)[particle_pair])*0.5; 
			}
			break;

		case 'a':
			connvec = subtract_containers ( p2->coords, p1->coords );
			modified_direction (&connvec, x, y, z); 
			magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
			theta_1   = branchless_acos(take_dot_product (scale_containers (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
			theta_2   = branchless_acos(take_dot_product (scale_containers (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );

			if ( (theta_1 + theta_2) > M_PI/2 ){
				(db_contacts) [std::get<4>((this->InteractionMap)[particle_pair])] += 0.5; 
				db_energy += std::get<2>((this->InteractionMap)[particle_pair])*0.5; 

			}
			else {
				(db_contacts) [std::get<3>((this->InteractionMap)[particle_pair])] += 0.5; 
				db_energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 
			}
			break;
	}

	return;
}

// calculate energy of the system
void Potts::debug_calculate_energy(double& db_energy, std::array<double,CONTACT_SIZE_POTTS>& db_contacts){

	db_energy   = 0;
	db_contacts = {0, 0};
	std::array <std::array<int,3>,26> ne_list; 

	for (Particle*& p: this->Lattice){
		ne_list = obtain_ne_list(p->coords, this->x, this->y, this->z);
		for (std::array <int,3>& loc: ne_list){
			this->debug_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], db_contacts, db_energy);
		}
	}

	return;

}

// this is essentially a checker function
void Potts::debug_checks_energy_contacts(double& E_final, std::array<double,CONTACT_SIZE_POTTS>& final_contacts){

	double E_debug = 0;
	std::array <double,CONTACT_SIZE_POTTS> debug_contacts = {0, 0};

	// run the debugger
	this->debug_calculate_energy(E_debug, debug_contacts);

	if (std::abs(E_debug - E_final) > 1e-8){
		std::cout << "E_debug - E_final = " << E_debug-E_final << std::endl;
		std::cout << "Energy mismatch." << std::endl;
		std::cout << "E_real = " << E_debug << ", E_final = " << E_final << "." << std::endl;
		std::cout << "real_contacts  = "; print(debug_contacts); 
		std::cout << "final_contacts = "; print(final_contacts);
		exit(EXIT_FAILURE);	
	}
	

	for (int i{0}; i<CONTACT_SIZE_POTTS; ++i){
		if (debug_contacts[i] != final_contacts[i]){
			std::cout << "Contacts mismatch." << std::endl;
			std::cout << "real_contacts  = "; print(debug_contacts); 
			std::cout << "final_contacts = "; print(final_contacts);
			exit(EXIT_FAILURE);
		}
		else {
			;
		}
	}

	std::cout << "energies match. debug = " << E_debug <<", computed = " << E_final << "." << std::endl;
	std::cout << "contacts match @ "; print(final_contacts, ""); std::cout << "!" << std::endl;

	return;

}

// start flipping the lattice
void Potts::debug_lattice_flip(){

	// update the number of attempts
	this->n_attempts += 1;

	// number of particles to flip
	int nflips = rng_uniform(1, static_cast<int>(this->x * this->y * this->z / 8));
	std::array <double,CONTACT_SIZE_POTTS> contacts = this->contacts; 

	// set up the orientations
	std::vector <int> orientations (nflips, 0); 

	// set up the contacts store
	std::array <double,CONTACT_SIZE_POTTS> cs_i;
    cs_i.fill(0);
	std::array <double,CONTACT_SIZE_POTTS> cs_f;
    cs_f.fill(0);

	// set up the energies
	double Es_i = 0;
	double Es_f = 0;
	double Ef   = this->Energy; 

	// set up the samples
	std::vector<int> samples(this->x * this->y * this->z - 1); 
	std::iota(samples.begin(), samples.end(), 0);

	// randomize
	std::random_device rd;
	std::mt19937 g(rd()); 
	std::shuffle(samples.begin(), samples.end(), g); 

	for (int i{0}; i<nflips; ++i){
		
		// get the neighbor energy 
		this->neighbor_energetics(samples[i], cs_i, Es_i);

		// store the old orientations 
		orientations[i] = (this->Lattice[samples[i]])->orientation;

		// perturb the orientation
		this->Lattice[samples[i]]->orientation = rng_uniform(0, 25);

		// get the new neighbor energy
		this->neighbor_energetics(samples[i], cs_f, Es_f); 

		// update the contacts and energy
		contacts = add_containers(subtract_containers(contacts, cs_i), cs_f);
		Ef += Es_f - Es_i;

		// resetting...
		cs_i.fill(0);
		cs_f.fill(0);
		Es_i = 0;
		Es_f = 0;

	}

	// run the debugging check
	this->debug_checks_energy_contacts(Ef, contacts);

	// set up the rng
	double rng = rng_uniform(0.0, 1.0);
	if (rng < std::exp(-1/this->T * (Ef - this->Energy))){
		this->Energy          = Ef;
		this->contacts        = contacts;
		this->acceptance      = true;
	}
	else {
		// go back to old orientations
		for (int i{0}; i<nflips; ++i){
			this->Lattice[samples[i]]->orientation = orientations[i];
		}
		this->acceptance = false;
	}
	return;
}

// start running the debugging simulation
void Potts::debug_simulation() {

	int start = this->step_number + 1;
	int stop  = this->step_number+this->max_iter+1;

	// run the primary loop
	for (int i{start}; i < stop; ++i) {
		this->step_number += 1;
		this->debug_lattice_flip();
		if (this->acceptance){
			this->n_acceptances += 1;
		}

		// start dumping out information
		this->dump();

	}

	return;
}

// start running the debugging
void Potts::debug(){

	// set up the simulation
	this->setup();

	// print out the opening tiles
	this->print_opening_tiles();

	// run the simulation!
	this->debug_simulation();

	// final dumps of lattice and move statistics
	this->dump_lattice();
	this->dump_stats();

	return;
}
