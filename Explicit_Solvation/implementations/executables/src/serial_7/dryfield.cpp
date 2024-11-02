#include "Simulation.h"
#include "lattice_directions.h"
#include "misc.h"

//////////////////////////////////////////////////////////
// modular functions for sampling based on orientations
//////////////////////////////////////////////////////////

void Simulation::perturb_orientation_sampler_forwards_dryfield(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iterator_idx, int lat_idx){

	
	this->enhanced_flipper.initial_orientations[iterator_idx] = this->Lattice[lat_idx]->orientation;
	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.initial_contacts), &(this->enhanced_flipper.initial_E));
	std::vector <int> test_orientations (26, 0);
	std::iota(test_orientations.begin(), test_orientations.end(), 0);

	// set up the random number generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(test_orientations.begin(), test_orientations.end(), std::default_random_engine(seed));

	// magnetization 
	double d_magnitude = norm(this->polymer_magnetization);
	double c_magnitude {0};
	std::array <double,3> d_magnetization   = subtract_arrays(&this->polymer_magnetization - &Or2Dir[this->enhanced_flipper.initial_orientations[iterator_idx]]);
	std::array <double,3> c_magnetization   = {0, 0, 0};

	for (int j{0}; j<this->enhanced_flipper.ntest; ++j){
		(this->Lattice)[lat_idx]->orientation  = test_orientations[j];
		c_magnetization                        = add_arrays(&(d_magnetization), &(Or2Dir[test_orientations[j]]));
		this->enhanced_flipper.orientations[j] = this->Lattice[lat_idx]->orientation;
		this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
		this->enhanced_flipper.energies[j]        = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E - this->Bfield*this->polymer_magnetization * this->Bfield * norm(c_magnetization);
		this->enhanced_flipper.contacts_store[j]  = subtract_arrays(contacts_sys, &(this->enhanced_flipper.initial_contacts));
		this->enhanced_flipper.contacts_store[j]  = add_arrays     (&(this->enhanced_flipper.contacts_store[j]), &(this->enhanced_flipper.perturbed_contacts));
		this->enhanced_flipper.perturbed_contacts = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		this->enhanced_flipper.perturbed_E        = 0;
	}

	return;
}

void Simulation::perturb_choose_state_forward_dryfield(int iterator_idx, int lat_idx){

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

void Simulation::perturb_orientation_sampler_backwards_0_dryfield(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx){

	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.initial_contacts), &(this->enhanced_flipper.initial_E));

	this->Lattice[lat_idx]->orientation = this->enhanced_flipper.initial_orientations[iteration_idx];
	this->neighbor_energetics(lat_idx, &(this->enhanced_flipper.perturbed_contacts), &(this->enhanced_flipper.perturbed_E));
	this->enhanced_flipper.energies[0] = E_sys - this->enhanced_flipper.initial_E + this->enhanced_flipper.perturbed_E;
	this->enhanced_flipper.contacts_store[0]  = subtract_arrays(*contacts_sys, this->enhanced_flipper.initial_contacts);
	this->enhanced_flipper.contacts_store[0]  = add_arrays     (this->enhanced_flipper.contacts_store[0], this->enhanced_flipper.perturbed_contacts);

	this->enhanced_flipper.perturbed_contacts = {0,0,0,0,0,0,0,0,0,0};
	this->enhanced_flipper.perturbed_E        = 0;

	return;
}

void Simulation::perturb_orientation_sampler_backwards_dryfield(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx){

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
		this->enhanced_flipper.contacts_store [j]  = subtract_arrays (*contacts_sys, this->enhanced_flipper.initial_contacts);
		this->enhanced_flipper.contacts_store [j]  = add_arrays      (this->enhanced_flipper.contacts_store[j], this->enhanced_flipper.perturbed_contacts);
		this->enhanced_flipper.perturbed_contacts = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		this->enhanced_flipper.perturbed_E        = 0;
	}

	return;
}

//////////////////////////////////////////////////////////
// end of modular functions for sampling based on orientations
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// get the dry field energy
//////////////////////////////////////////////////////////

void Simulation::accelerate_calculate_energy_dry_field(){
	
	std::cout << "Running through the polymers for an energy computation." << std::endl;
	double net_align {0.0};
	double energy    {0.0};
	this->contacts = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	
	std::array <double,3> magnetization = {0.0, 0.0, 0.0};
	std::array <std::array <int,3>, 26> ne_list; 

	// run energy computations for every monomer bead 
	for (Polymer& pmer: this->Polymers) {
		for (Particle*& p: pmer.chain){
			magnetization = add_arrays(&alignment, &(Or2Dir[p->orientation]))
			ne_list       = obtain_ne_list(p->coords, this->x, this->y, this->z); // get neighbor list 
			for (std::array <int,3>& loc: ne_list) {
				this->accelerate_pair_interaction(p, this->Lattice[lattice_index(loc, this->y, this->z)], &this->contacts, &energy);
			}
		}
	}

	this->polymer_magnetization  = magnetization;
	this->sysEnergy             += energy;
	this->sysEnergy             += this->Bfield * std::sqrt(this->polymer_magnetization[0]*this->polymer_magnetization[0] + this->polymer_magnetization[1]*this->polymer_magnetization[1] + this->polymer_magnetization[2]*this->polymer_magnetization[2]);
	return; 
}

