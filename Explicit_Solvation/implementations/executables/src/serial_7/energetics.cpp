#include "Simulation.hpp"
#include "lattice_directions.hpp"
#include "misc.hpp"

//////////////////////////////////////////////////////////
//
//                Run the energetic calculation
//
//////////////////////////////////////////////////////////

void interaction_s1_s1(Simulation*, Particle*, Particle*,   std::array<double,CONTACT_SIZE>*, double*) {
	return; 
}

void interaction_s2_s2(Simulation*, Particle*, Particle*,   std::array<double,CONTACT_SIZE>*, double*) {
	return; 
}

void interaction_i_m1_m1(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	(*contacts)[1] += 0.5;
	(*energy)      += 0.5 * sim->energy_surface[0];
	return; 
}

void interaction_i_m1_s1(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	(*contacts)[3] += 1;
	(*energy)      += sim->energy_surface[2];
	return;
}

void interaction_i_m1_s2(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	(*contacts)[5] += 1;
	(*energy)      += sim->energy_surface[4];
	return; 
}

void interaction_i_s1_s2(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	(*contacts)[7] += 1;
	(*energy)      += sim->energy_surface[6];
	return; 
}

void interaction_p_m1_m1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = (dot_product > 0.54); 
	int  c_idx          = 1 - condition;
	(*contacts)[c_idx] += 0.5; 
	(*energy)          += 0.5 * ((sim->energy_surface[0]) * condition + (sim->energy_surface[1]) * !condition);
	return; 
}

void interaction_p_m1_s1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = dot_product > 0.54;
	int  c_idx          = 3 - condition;
	(*contacts)[c_idx] += 1;
	(*energy)          += (sim->energy_surface[2]) * condition + (sim->energy_surface[3]) * !condition;
	return; 
}

void interaction_p_m1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = dot_product > 0.54; 
	int  c_idx          = 5 - condition;
	(*contacts)[c_idx] += 1;
	(*energy)          += (sim->energy_surface[4]) * condition + (sim->energy_surface[5]) * !condition;
	return; 
}

void interaction_p_s1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = dot_product > 0.54;
	int  c_idx          = 7 - condition;
	(*contacts)[c_idx] += 1;
	(*energy)          += (sim->energy_surface[6]) * condition + (sim->energy_surface[7]) * !condition;
	return; 
}

void interaction_a_m1_m1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	// std::cout << "Entering function." << std::endl;
	std::array <int,3> connvec = subtract_containers (p2->coords, p1->coords);
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	// std::cout << "Evaluating the if condition." << std::endl;
	double magnitude    = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1      = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2      = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool   condition    = (theta_1 + theta_2) > M_PI/2;
	int    c_idx        = 0 + condition;
	(*contacts)[c_idx] += 0.5;
	// std::cout << "Doing energy compute..." << std::endl;
	(*energy)          += 0.5 * ((sim->energy_surface[0]) * !condition + (sim->energy_surface[1]) * condition);
	return; 
}

void interaction_a_m1_s1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	std::array <int,3> connvec = subtract_containers ((p2->coords), (p1->coords));
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude    = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1      = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2      = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition      = (theta_1 + theta_2) > M_PI/2; 
	int c_idx           = 2 + condition;
	(*contacts)[c_idx] += 1;
	*energy            += (sim->energy_surface[2]) * !condition + (sim->energy_surface[3]) * condition;
	return; 
}

void interaction_a_m1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	std::array <int,3> connvec = subtract_containers ((p2->coords), (p1->coords));
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude    = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1      = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2      = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition      = (theta_1 + theta_2) > M_PI/2; 
	int  c_idx          = 4 + condition;
	(*contacts)[c_idx] += 1;
	(*energy)          += (sim->energy_surface[4]) * !condition + (sim->energy_surface[5]) * condition;
	return; 
}

void interaction_a_s1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	std::array <int,3> connvec = subtract_containers ((p2->coords), (p1->coords));
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude    = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1      = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2      = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition      = (theta_1 + theta_2) > M_PI/2; 
	int c_idx           = 6 + condition;
	// std::cout << "\t@ p1 = "; print(p1->coords, ", "); std::cout << "p2 = "; print(p2->coords, ", "); std::cout << "Alignment: " << !condition << "." << std::endl;
	(*contacts)[c_idx] += 1;
	(*energy)          += (sim->energy_surface[6]) * !condition + (sim->energy_surface[7]) * condition;
	return; 
}

void Simulation::initialize_pairwise_function_map(){

	std::pair  <std::string, std::string> mm_pair    = std::make_pair("m1", "m1");
	std::pair  <std::string, std::string> ms1_pair   = std::make_pair("m1", "s1");
	std::pair  <std::string, std::string> s1m_pair   = std::make_pair("s1", "m1");
	std::pair  <std::string, std::string> ms2_pair   = std::make_pair("m1", "s2");
	std::pair  <std::string, std::string> s2m_pair   = std::make_pair("s2", "m1");
	std::pair  <std::string, std::string> s1s1_pair  = std::make_pair("s1", "s1");
	std::pair  <std::string, std::string> s2s2_pair  = std::make_pair("s2", "s2");
	std::pair  <std::string, std::string> s1s2_pair  = std::make_pair("s1", "s2");
	std::pair  <std::string, std::string> s2s1_pair  = std::make_pair("s2", "s1");

	this->PairwiseFunctionMap[s1s1_pair] = &interaction_s1_s1;
	this->PairwiseFunctionMap[s2s2_pair] = &interaction_s2_s2;

	if (std::get<0>((this->InteractionMap)[mm_pair]) == "isotropic"){
		// std::cout << "pairwise m1-m1 is isotropic." << std::endl;
		this->PairwiseFunctionMap[mm_pair] = &interaction_i_m1_m1;
	}
	else if (std::get<0>((this->InteractionMap)[mm_pair]) == "parallel"){
		// std::cout << "pairwise m1-m1 is parallel." << std::endl;
		this->PairwiseFunctionMap[mm_pair] = &interaction_p_m1_m1;
	}
	else if (std::get<0>((this->InteractionMap)[mm_pair]) == "antiparallel"){
		// std::cout << "pairwise m1-m1 is antiparallel." << std::endl;
		this->PairwiseFunctionMap[mm_pair] = &interaction_a_m1_m1;
	}
	else {
		std::cout << "No interaction provided for m-m." << std::endl; 
		exit(EXIT_FAILURE); 
	}

	if (std::get<0>((this->InteractionMap)[ms1_pair]) == "isotropic"){
		// std::cout << "pairwise m1-s1 is isotropic." << std::endl;
		this->PairwiseFunctionMap[ms1_pair] = &interaction_i_m1_s1;
		this->PairwiseFunctionMap[s1m_pair] = &interaction_i_m1_s1;
	}
	else if (std::get<0>((this->InteractionMap)[ms1_pair]) == "parallel"){
		// std::cout << "pairwise m1-s1 is parallel." << std::endl;
		this->PairwiseFunctionMap[ms1_pair] = &interaction_p_m1_s1;
		this->PairwiseFunctionMap[s1m_pair] = &interaction_p_m1_s1;
	}
	else if (std::get<0>((this->InteractionMap)[ms1_pair]) == "antiparallel"){
		// std::cout << "pairwise m1-s1 is antiparallel." << std::endl;
		this->PairwiseFunctionMap[ms1_pair] = &interaction_a_m1_s1;
		this->PairwiseFunctionMap[s1m_pair] = &interaction_a_m1_s1;
	}
	else {
		std::cout << "No interaction provided for m-s1." << std::endl; 
		exit(EXIT_FAILURE);
	}

	if (std::get<0>((this->InteractionMap)[ms2_pair]) == "isotropic"){
		// std::cout << "pairwise m1-s2 is isotropic." << std::endl;
		this->PairwiseFunctionMap[ms2_pair] = &interaction_i_m1_s2;
		this->PairwiseFunctionMap[ms2_pair] = &interaction_i_m1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[ms2_pair]) == "parallel"){
		// std::cout << "pairwise m1-s2 is parallel." << std::endl;
		this->PairwiseFunctionMap[ms2_pair] = &interaction_p_m1_s2;
		this->PairwiseFunctionMap[s2m_pair] = &interaction_p_m1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[ms2_pair]) == "antiparallel"){
		// std::cout << "pairwise m1-s2 is antiparallel." << std::endl;
		this->PairwiseFunctionMap[ms2_pair] = &interaction_a_m1_s2;
		this->PairwiseFunctionMap[s2m_pair] = &interaction_a_m1_s2;
	}
	else {
		std::cout << "No interactions provided for m-s2." << std::endl; 
		exit(EXIT_FAILURE);
	}

	if (std::get<0>((this->InteractionMap)[s1s2_pair]) == "isotropic"){
		// std::cout << "pairwise s1-s2 is isotropic." << std::endl;
		this->PairwiseFunctionMap[s1s2_pair] = &interaction_i_s1_s2;
		this->PairwiseFunctionMap[s2s1_pair] = &interaction_i_s1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[s1s2_pair]) == "parallel"){
		// std::cout << "pairwise s1-s2 is parallel." << std::endl;
		this->PairwiseFunctionMap[s1s2_pair] = &interaction_p_s1_s2;
		this->PairwiseFunctionMap[s2s1_pair] = &interaction_p_s1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[s1s2_pair]) == "antiparallel"){
		// std::cout << "pairwise s1-s2 is antiparallel." << std::endl;
		this->PairwiseFunctionMap[s1s2_pair] = &interaction_a_s1_s2;
		this->PairwiseFunctionMap[s2s1_pair] = &interaction_a_s1_s2;
	}
	else {
		std::cout << "No interactions for s1-s2." << std::endl;
		exit(EXIT_FAILURE);
	}
	

	return;
}

void Simulation::accelerate_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy){
	std::pair<std::string, std::string> particle_pair = std::make_pair(p1->ptype, p2->ptype);

	auto it = this->PairwiseFunctionMap.find(particle_pair);
	it->second(this, p1, p2, contacts, energy);
	return;
}

// the following function has been replaced by accelerate_pair_interaction

// %~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~
// %~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~
// %~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~
// %~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~

void neighbor_s1_s1(Simulation*, Particle*, Particle*, std::array<double,CONTACT_SIZE>*, double*) {
	return;
}

void neighbor_s2_s2(Simulation*, Particle*, Particle*, std::array<double,CONTACT_SIZE>*, double*) {
	return;
}

void neighbor_i_m1_m1(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy) {
	(*contacts)[1] += 1;
	(*energy)      += sim->energy_surface[1];
	return; 
}

void neighbor_i_m1_s1(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
    (*contacts)[3] += 1;
    *energy_incr   += sim->energy_surface[3];
	return; 
}

void neighbor_i_m1_s2(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	(*contacts)[5] += 1;
	*energy_incr   += sim->energy_surface[5];
	return; 
}

void neighbor_i_s1_s2(Simulation* sim, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	(*contacts)[7]  += 1;
	*energy_incr    += sim->energy_surface[7];
	return; 
}

void neighbor_p_m1_m1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	double dot_product = take_dot_product(p1->orientation, p2->orientation);
	bool condition     = (dot_product > 0.54); 
	int  c_idx         = 1 - condition;
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[0]) * condition + (sim->energy_surface[1]) * !condition;
	return; 
}

void neighbor_p_m1_s1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = dot_product > 0.54;
	int  c_idx          = 3 - condition; 
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[2]) * condition + (sim->energy_surface[3]) * !condition;
	return; 
}

void neighbor_p_m1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = dot_product > 0.54; 
	int c_idx           = 5 - condition;
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[4]) * condition + (sim->energy_surface[5]) * !condition;
	return; 
}

void neighbor_p_s1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	double dot_product  = take_dot_product(p1->orientation, p2->orientation);
	bool condition      = dot_product > 0.54;
	int c_idx           = 7 - condition;
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[6]) * condition + (sim->energy_surface[7]) * !condition;
	return; 
}

void neighbor_a_m1_m1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	// std::cout << "\t(neighbor): Particles @ "; print(p1->coords, ", "); print(p2->coords);
	std::array <int,3> connvec = subtract_containers ((p2->coords), (p1->coords));
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1   = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );// std::acos (take_dot_product (scale_containers (1/magnitude , &connvec), Or2Dir[p1->orientation] ) );
	double theta_2   = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );// std::acos (take_dot_product (scale_containers (-1/magnitude, &connvec), Or2Dir[p2->orientation] ) );
	bool condition   = (theta_1 + theta_2) > M_PI/2;
	int       c_idx  = 0 + condition;
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[0]) * !(condition) + (sim->energy_surface[1]) * condition;
	return; 
}

void neighbor_a_m1_s1(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	std::array <int,3> connvec = subtract_containers ((p2->coords), (p1->coords));
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1   = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2   = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition   = (theta_1 + theta_2) > M_PI/2; 
	int  c_idx       = 2 + condition;
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[2]) * !condition + (sim->energy_surface[3]) * condition;
	return; 
}

void neighbor_a_m1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	std::array <int,3> connvec = subtract_containers ((p2->coords), (p1->coords));
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1   = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2   = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition   = (theta_1 + theta_2) > M_PI/2; 
	int  c_idx       = 4 + condition;
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[4]) * !condition + (sim->energy_surface[5]) * condition;
	return; 
}

void neighbor_a_s1_s2(Simulation* sim, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr) {
	std::array <int,3> connvec = subtract_containers ((p2->coords), (p1->coords));
	modified_direction (&connvec, sim->x, sim->y, sim->z); 
	double magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
	double theta_1   = branchless_acos (take_dot_product (scale_arrays (1/magnitude , connvec), Or2Dir[p1->orientation] ) );
	double theta_2   = branchless_acos (take_dot_product (scale_arrays (-1/magnitude, connvec), Or2Dir[p2->orientation] ) );
	bool condition   = (theta_1 + theta_2) > M_PI/2; 
	int  c_idx       = 6 + condition;
	(*contacts)[c_idx] += 1;
	*energy_incr       += (sim->energy_surface[6]) * !condition + (sim->energy_surface[7]) * condition;
	return; 
}

void Simulation::initialize_neighbor_function_map(){

	std::pair  <std::string, std::string> mm_pair    = std::make_pair("m1", "m1");
	std::pair  <std::string, std::string> ms1_pair   = std::make_pair("m1", "s1");
	std::pair  <std::string, std::string> s1m_pair   = std::make_pair("s1", "m1");
	std::pair  <std::string, std::string> ms2_pair   = std::make_pair("m1", "s2");
	std::pair  <std::string, std::string> s2m_pair   = std::make_pair("s2", "m1");
	std::pair  <std::string, std::string> s1s1_pair  = std::make_pair("s1", "s1");
	std::pair  <std::string, std::string> s2s2_pair  = std::make_pair("s2", "s2");
	std::pair  <std::string, std::string> s1s2_pair  = std::make_pair("s1", "s2");
	std::pair  <std::string, std::string> s2s1_pair  = std::make_pair("s2", "s1");

	this->NeighborFunctionMap[s1s1_pair] = &neighbor_s1_s1;
	this->NeighborFunctionMap[s2s2_pair] = &neighbor_s2_s2;

	if (std::get<0>((this->InteractionMap)[mm_pair]) == "isotropic"){
		// std::cout << "neighbor m1-m1 is isotropic" << std::endl;
		this->NeighborFunctionMap[mm_pair] = &neighbor_i_m1_m1;
	}
	else if (std::get<0>((this->InteractionMap)[mm_pair]) == "parallel"){
		// std::cout << "neighbor m1-m1 is parallel." << std::endl;
		this->NeighborFunctionMap[mm_pair] = &neighbor_p_m1_m1;
	}
	else if (std::get<0>((this->InteractionMap)[mm_pair]) == "antiparallel"){
		// std::cout << "neighbor m1-m1 is antiparallel." << std::endl;
		this->NeighborFunctionMap[mm_pair] = &neighbor_a_m1_m1;
	}
	else {
		std::cout << "No interaction for mm." << std::endl; 
		exit(EXIT_FAILURE); 
	}

	if (std::get<0>((this->InteractionMap)[ms1_pair]) == "isotropic"){
		// std::cout << "neighbor m1-s1 is isotropic." << std::endl;
		this->NeighborFunctionMap[ms1_pair] = &neighbor_i_m1_s1;
		this->NeighborFunctionMap[s1m_pair] = &neighbor_i_m1_s1;
	}
	else if (std::get<0>((this->InteractionMap)[ms1_pair]) == "parallel"){
		// std::cout << "neighbor m1-s1 is parallel." << std::endl;
		this->NeighborFunctionMap[ms1_pair] = &neighbor_p_m1_s1;
		this->NeighborFunctionMap[s1m_pair] = &neighbor_p_m1_s1;
	}
	else if (std::get<0>((this->InteractionMap)[ms1_pair]) == "antiparallel"){
		// std::cout << "neighbor m1-s1 is antiparallel." << std::endl;
		this->NeighborFunctionMap[ms1_pair] = &neighbor_a_m1_s1;
		this->NeighborFunctionMap[s1m_pair] = &neighbor_a_m1_s1;
	}
	else {
		std::cout << "No interactions for m-s1." << std::endl; 
		exit(EXIT_FAILURE);
	}

	if (std::get<0>((this->InteractionMap)[ms2_pair]) == "isotropic"){
		// std::cout << "neighbor m1-s2 is isotropic." << std::endl;
		this->NeighborFunctionMap[ms2_pair] = &neighbor_i_m1_s2;
		this->NeighborFunctionMap[s2m_pair] = &neighbor_i_m1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[ms2_pair]) == "parallel"){
		// std::cout << "neighbor m1-s2 is parallel." << std::endl;
		this->NeighborFunctionMap[ms2_pair] = &neighbor_p_m1_s2;
		this->NeighborFunctionMap[s2m_pair] = &neighbor_p_m1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[ms2_pair]) == "antiparallel"){
		// std::cout << "neighbor m1-s2 is antiparallel." << std::endl;
		this->NeighborFunctionMap[ms2_pair] = &neighbor_a_m1_s2;
		this->NeighborFunctionMap[s2m_pair] = &neighbor_a_m1_s2;
	}
	else {
		std::cout << "No interactions for m-s2." << std::endl; 
		exit(EXIT_FAILURE);
	}

	if (std::get<0>((this->InteractionMap)[s1s2_pair]) == "isotropic"){
		// std::cout << "neighbor s1-s2 is isotropic." << std::endl;
		this->NeighborFunctionMap[s1s2_pair] = &neighbor_i_s1_s2;
		this->NeighborFunctionMap[s2s1_pair] = &neighbor_i_s1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[s1s2_pair]) == "parallel"){
		// std::cout << "neighbor s1-s2 is parallel." << std::endl;
		this->NeighborFunctionMap[s1s2_pair] = &neighbor_p_s1_s2;
		this->NeighborFunctionMap[s2s1_pair] = &neighbor_p_s1_s2;
	}
	else if (std::get<0>((this->InteractionMap)[s1s2_pair]) == "antiparallel"){
		// std::cout << "neighbor s1-s2 is antiparallel." << std::endl;
		this->NeighborFunctionMap[s1s2_pair] = &neighbor_a_s1_s2;
		this->NeighborFunctionMap[s2s1_pair] = &neighbor_a_s1_s2;
	}
	else {
		std::cout << "No interactions for s1-s2." << std::endl; 
		exit(EXIT_FAILURE);
	}

	return;
}

// pair interaction from ParticlePairEnergyContribution
// need to figure out when and where this function is being used... 
void Simulation::selected_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy){

	std::array <int,3> connvec = subtract_containers((p2->coords), (p1->coords));
	modified_direction (&connvec, x, y, z);

	// check if given particles are neighbors
	if (std:: find (adrns.begin(), adrns.end(), connvec) != adrns.end()){
		std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 
		auto it = this->NeighborFunctionMap.find(particle_pair);
		if (it != this->NeighborFunctionMap.end()) {
			// Function pointer found, call the function
			it->second(this, p1, p2, contacts, energy);
		}
		else {
			// Function pointer not found, handle the error
			std::cerr << "Error: Neighbor function not found for particle pair (" << p1->ptype << ", " << p2->ptype << ")" << std::endl;
			// Optionally, you can exit or throw an exception
			exit(EXIT_FAILURE);
		}
	}
	return; 
}

// %~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~
// %~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~
// this function goes to a particle and measures the energy of interactions between that 
// particle and all it's neighbors

void Simulation::neighbor_energetics(int lat_idx, std::array<double,CONTACT_SIZE>* contacts_store, double* neighbor_energy){
	std::array <std::array<int,3>,26> ne_list = obtain_ne_list(this->Lattice[lat_idx]->coords, this->x, this->y, this->z);
	for (std::array<int,3>& loc: ne_list){
		this->selected_pair_interaction(this->Lattice[lat_idx], this->Lattice[lattice_index(loc, this->y, this->z)], contacts_store, neighbor_energy);
	}
	return;
}

// %~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~%~

void Simulation::accelerate_calculate_energy_cosolvent(){

	double energy       {0.0};
	this->contacts = {0,0,0, 0,0,0, 0,0};
	std::array <std::array <int,3>, 26> ne_list; 

	// run energy computations for every monomer bead 
	// auto start = std::chrono::high_resolution_clock::now();
	std::cout << "Running through the polymers for an energy computation. Number of polymers = " << this->Polymers.size() << std::endl;
	for (Polymer& pmer: this->Polymers) {
		for (Particle*& p: pmer.chain){
			ne_list = obtain_ne_list(p->coords, this->x, this->y, this->z); // get neighbor list 
			for (std::array <int,3>& loc: ne_list) {
				this->accelerate_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], &this->contacts, &energy);
			}
		}
	}

	std::cout << "Running through the cosolvent for an energy computation..." << std::endl;
	std::cout << "Length of Cosolvent vector is " << this->Cosolvent.size() << "." << std::endl;
	for (Particle*& p: this->Cosolvent){
		ne_list = obtain_ne_list ( p->coords, this->x, this->y, this->z );
		for (std::array <int,3>& loc: ne_list){
			if (this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "m1" || this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "s2"){
				continue; 
			}
			else {
				this->accelerate_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], &this->contacts, &energy);
			}
		}
	}

	this->sysEnergy += energy;
	std::cout << "Completed computation." << std::endl;
	return; 
}

void Simulation::accelerate_calculate_energy_solvent(){

	std::cout << "Running through the solvents for an energy computation." << std::endl;
	double energy {0.0};
	this->contacts = {0,0,0, 0,0,0, 0,0};

	std::array <std::array <int,3>, 26> ne_list; 

	// run energy computations for every monomer bead 
	// auto start = std::chrono::high_resolution_clock::now();

	for (Polymer& pmer: this->Polymers) {
		for (Particle*& p: pmer.chain){
			ne_list = obtain_ne_list(p->coords, this->x, this->y, this->z); // get neighbor list 
			for (std::array <int,3>& loc: ne_list) {
				this->accelerate_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], &this->contacts, &energy);
			}
		}
	}

	std::cout << "Size of solvent vector is " << this->Solvent.size() << "." << std::endl;

	for (Particle*& p: this->Solvent){
		ne_list = obtain_ne_list ( p->coords, this->x, this->y, this->z );
		for (std::array <int,3>& loc: ne_list){
			if (this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "m1" || this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "s1"){
				continue; 
			}
			else {
				// std::cout << "Solvent-Cosolvent contact." << std::endl;
				this->accelerate_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], &this->contacts, &energy);
			}
		}
	}

	this->sysEnergy += energy;
	return; 
}

