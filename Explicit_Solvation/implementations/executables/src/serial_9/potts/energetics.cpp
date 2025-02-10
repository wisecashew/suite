#include "Potts.hpp"

void initial_i_sp_sp(Potts* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	(contacts)[0] += 0.5;
	(energy)      += 0.5 * sim->energy_surface[0];
	return; 
}

void initial_p_sp_sp(Potts* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = (dot_product > 0.54); 
	int  c_idx          = 1 - condition;
	(contacts)[c_idx] += 0.5; 
	(energy)          += 0.5 * ((sim->energy_surface[0]) * condition + (sim->energy_surface[1]) * !condition);
	return; 
}

void initial_a_sp_sp(Potts* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	std::array <int,3> connvec = subtract_containers (p2->coords, p1->coords);
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude    = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1      = branchless_acos (take_dot_product (scale_containers ( 1/magnitude, connvec), Or2Dir[p1->orientation] ) );
	double theta_2      = branchless_acos (take_dot_product (scale_containers (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition      = (theta_1 + theta_2) > M_PI/2; 
	int c_idx           = 0 + condition;
	(contacts)[c_idx] += 0.5;
	energy            += 0.5 * ((sim->energy_surface[0]) * !condition + (sim->energy_surface[1]) * condition);
	return;
}

void initial_symm_sp_sp(Potts* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = (std::fabs(dot_product) > 0.54); 
	int  c_idx          = 1 - condition;
	(contacts)[c_idx] += 0.5; 
	(energy)          += 0.5 * ((sim->energy_surface[0]) * condition + (sim->energy_surface[1]) * !condition);
	return; 
}

// NEED TO INITIALIZE THIS
void Potts::initial_energetics_map(){

	std::pair  <std::string, std::string> spsp_pair  = std::make_pair("sp", "sp");
	if (std::get<0>((this->InteractionMap)[spsp_pair]) == "isotropic"){
		// std::cout << "pairwise s1-s2 is isotropic." << std::endl;
		this->InitialInteractionFunctionMap[spsp_pair] = &initial_i_sp_sp;
	}
	else if (std::get<0>((this->InteractionMap)[spsp_pair]) == "parallel"){
		// std::cout << "pairwise s1-s2 is parallel." << std::endl;
		this->InitialInteractionFunctionMap[spsp_pair] = &initial_p_sp_sp;
	}
	else if (std::get<0>((this->InteractionMap)[spsp_pair]) == "antiparallel"){
		// std::cout << "pairwise s1-s2 is antiparallel." << std::endl;
		this->InitialInteractionFunctionMap[spsp_pair] = &initial_a_sp_sp;
	}

	else if (std::get<0>((this->InteractionMap)[spsp_pair]) == "symmetric"){
		// std::cout << "pairwise sp-sp is symmetric." << std::endl;
		this->InitialInteractionFunctionMap[spsp_pair] = &initial_symm_sp_sp;
	}
	else {
		std::cout << "No interactions for sp-sp." << std::endl;
		exit(EXIT_FAILURE);
	}

	return;
}

void Potts::initial_energetics(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	std::pair<std::string, std::string> particle_pair = std::make_pair(p1->ptype, p2->ptype);
	auto it = this->InitialInteractionFunctionMap.find(particle_pair);
	it->second(this, p1, p2, contacts, energy);
	return;
}

// THE ABOVE IS ONLY FOR THE INITIAL INTERACTIONS

void neighbor_i_sp_sp(Potts* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	(contacts)[0] += 1;
	(energy)      += 1 * sim->energy_surface[0];
	return; 
}

void neighbor_p_sp_sp(Potts* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = (dot_product > 0.54); 
	int  c_idx          = 1 - condition;
	(contacts)[c_idx]  += 1; 
	(energy)           += 1 * ((sim->energy_surface[0]) * condition + (sim->energy_surface[1]) * !condition);
	return; 
}

void neighbor_a_sp_sp(Potts* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = (std::fabs(dot_product) > 0.54);
	int  c_idx          = 1 - condition;
	(contacts)[c_idx]  += 1;
	(energy)           += 1 * ((sim->energy_surface[0]) * condition + (sim->energy_surface[1]) * !condition);
	return; 
}

void neighbor_symm_sp_sp(Potts* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){
	std::array <int,3> connvec = subtract_containers (p2->coords, p1->coords);
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude    = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1      = branchless_acos (take_dot_product (scale_containers (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2      = branchless_acos (take_dot_product (scale_containers (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition      = (theta_1 + theta_2) > M_PI/2; 
	int c_idx           = 0 + condition;
	(contacts)[c_idx]  += 1;
	(energy)           += (sim->energy_surface[0]) * !condition + (sim->energy_surface[1]) * condition;
	return;
}

// NEED TO INITIALIZE THIS
void Potts::neighbor_energetic_map(){

	std::pair  <std::string, std::string> spsp_pair  = std::make_pair("sp", "sp");
	if (std::get<0>((this->InteractionMap)[spsp_pair]) == "isotropic"){
		// std::cout << "neighbor s1-s2 is isotropic." << std::endl;
		this->NeighborhoodInteractionFunctionMap[spsp_pair] = &neighbor_i_sp_sp;
	}
	else if (std::get<0>((this->InteractionMap)[spsp_pair]) == "parallel"){
		// std::cout << "neighbor s1-s2 is parallel." << std::endl;
		this->NeighborhoodInteractionFunctionMap[spsp_pair] = &neighbor_p_sp_sp;
	}
	else if (std::get<0>((this->InteractionMap)[spsp_pair]) == "antiparallel"){
		// std::cout << "neighbor s1-s2 is antiparallel." << std::endl;
		this->NeighborhoodInteractionFunctionMap[spsp_pair] = &neighbor_a_sp_sp;
	}
	else if (std::get<0>((this->InteractionMap)[spsp_pair]) == "symmetric"){
		// std::cout << "neighbor s1-s2 is antiparallel." << std::endl;
		this->NeighborhoodInteractionFunctionMap[spsp_pair] = &neighbor_symm_sp_sp;
	}

	else {
		std::cout << "No interactions for sp-sp." << std::endl; 
		exit(EXIT_FAILURE);
	}

	return;
}

void Potts::selected_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& contacts, double& energy){

	std::array <int,3> connvec = subtract_containers(p2->coords, p1->coords);
	modified_direction (&connvec, x, y, z);

	// check if given particles are neighbors
	if (std:: find (adrns.begin(), adrns.end(), connvec) != adrns.end()){
		std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 
		auto it = this->NeighborhoodInteractionFunctionMap.find(particle_pair);
		it->second(this, p1, p2, contacts, energy);
	}
	return; 

}

void Potts::neighbor_energetics(int lat_idx, std::array<double,CONTACT_SIZE_POTTS>& contacts_store, double& neighbor_energy){

	std::array <std::array<int,3>,26> ne_list = obtain_ne_list(this->Lattice[lat_idx]->coords, this->x, this->y, this->z);
	for (std::array<int,3>& loc: ne_list){
		this->selected_pair_interaction(this->Lattice[lat_idx], this->Lattice[lattice_index(loc, this->y, this->z)], contacts_store, neighbor_energy);
	}
	return;
}

// GET THE ENERGY OF THE SIMULATION
void Potts::energy_compute(){

	double energy {0.0};
	this->contacts.fill(0.0);
	std::array <std::array<int,3>,26> ne_list;

	for (Particle*& p: this->Lattice){
		ne_list = obtain_ne_list(p->coords, this->x, this->y, this->z);
		for (std::array <int,3>& loc: ne_list){
			this->initial_energetics(p, this->Lattice[lattice_index(loc, this->y, this->z)], this->contacts, energy);
		}
	}

	this->Energy = energy;
	return;

}