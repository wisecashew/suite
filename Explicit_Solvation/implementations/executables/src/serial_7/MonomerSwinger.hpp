#ifndef _MONOMERSWINGER_H_
#define _MONOMERSWINGER_H_
#include "misc.hpp"

class EnhancedMonomerSwinger {
public:

	// properties
	int ntest;
	std::array <int,3>    zero_3array;
	std::array <double,CONTACT_SIZE> zero_Narray; // a zero'd out array of contacts
	std::vector <std::array<int,3>> initial_locations;
	std::vector <std::array<int,3>> final_locations;
	std::vector <double> boltzmann; 

	// holders for energies as the orientation of each particle is perturbed
	double initial_E_monomer;              // energy of the system at the start
	double initial_E_part2switch;
	double final_E_monomer;
	double final_E_part2switch;
	double initial_E_pair;
	double final_E_pair; 
	std::vector <double> energies; // this is the container holding energies of each perturbation

	// containers for contacts as the orientation of each particle is perturbed
	std::array <double,CONTACT_SIZE>              initial_cont_monomer;     // initial contacts, prior to perturbation
	std::array <double,CONTACT_SIZE>              initial_cont_part2switch; // initial contacts, prior to perturbation
	std::array <double,CONTACT_SIZE>              final_cont_monomer;       // 
	std::array <double,CONTACT_SIZE>              final_cont_part2switch;   // 
	std::array <double,CONTACT_SIZE>              initial_cont_pair;        //
	std::array <double,CONTACT_SIZE>              final_cont_pair;          // 
	std::vector <std::array<double,CONTACT_SIZE>> contacts_store;     // this is the container holding all the net contacts when running through the trialed states

	// these are running sums and cumulant probabilities to keep track of
	double rboltzmann;
	double prob_o_to_n;
	double prob_n_to_o;
	double Emin;

	// create a holder for the current energy and contacts
	double               current_energy;
	std::array<double,CONTACT_SIZE> current_contacts;

	// rng holder for an internal decision
	int    sampler_idx;
	double sampler_rng;
	double sampler_rsum;

	// define the Particle pointer
	Particle* tmp_par_ptr;

	// constructor 
	EnhancedMonomerSwinger(){};

	// destructor
	~EnhancedMonomerSwinger(){};

	// set up the object 
	void setup(int loc_size, int ntest){
		this->ntest           = ntest;
		this->zero_3array     = {0,0,0};
		this->zero_Narray     = {0,0,0, 0,0,0, 0,0};
		this->initial_locations.reserve(loc_size);
		this->final_locations.reserve(loc_size);
		this->initial_cont_monomer     = this->zero_Narray;     // initial contacts, prior to perturbation
		this->initial_cont_part2switch = this->zero_Narray;     // initial contacts, prior to perturbation
		this->final_cont_monomer       = this->zero_Narray;      
		this->final_cont_part2switch   = this->zero_Narray;
		this->initial_cont_pair        = this->zero_Narray;
		this->final_cont_pair          = this->zero_Narray; 
		this->current_contacts         = this->zero_Narray;
		this->contacts_store.resize(this->ntest, this->zero_Narray);
		this->boltzmann.resize(this->ntest,0.0);
		this->energies.resize(this->ntest,0);
		this->rboltzmann       = 0;
		this->prob_o_to_n      = 1;
		this->prob_n_to_o      = 1;
		this->Emin             = 0;
		this->current_energy   = 0;
		this->sampler_idx      = 0;
		this->sampler_rng      = 0;
		this->sampler_rsum     = 0;
		this->tmp_par_ptr      = nullptr;
		return; 
	}

	// reset the object 
	void reset(int loc_size){ // fixed ntest
		this->initial_locations.clear();
		this->initial_locations.reserve(loc_size);
		this->final_locations.clear();
		this->final_locations.reserve  (loc_size);
		this->initial_cont_monomer     = this->zero_Narray;     // initial contacts, prior to perturbation
		this->initial_cont_part2switch = this->zero_Narray;     // initial contacts, prior to perturbation
		this->final_cont_monomer       = this->zero_Narray;     
		this->final_cont_part2switch   = this->zero_Narray;
		this->initial_cont_pair        = this->zero_Narray;
		this->final_cont_pair          = this->zero_Narray; 
		this->current_contacts         = this->zero_Narray;
		std::fill(this->contacts_store.begin(), this->contacts_store.end(), this->zero_Narray);
		std::fill(this->boltzmann.begin(), this->boltzmann.end(), 0);
		std::fill(this->energies.begin(), this->energies.end(), 0);
		this->rboltzmann       = 0;
		this->prob_o_to_n      = 1;
		this->prob_n_to_o      = 1;
		this->Emin             = 0;
		this->current_energy   = 0;
		this->sampler_idx      = 0;
		this->sampler_rng      = 0;
		this->sampler_rsum     = 0;
		this->tmp_par_ptr      = nullptr;
		return;
	}

	void reset_local(){
		// reset everything
		this->initial_E_monomer     = 0;
		this->initial_E_part2switch = 0;
		this->initial_E_pair        = 0;
		this->final_E_monomer       = 0;
		this->final_E_part2switch   = 0;
		this->final_E_pair          = 0;

		// reset the contacts
		this->initial_cont_monomer     = {0,0,0, 0,0,0, 0,0};
		this->initial_cont_part2switch = {0,0,0, 0,0,0, 0,0};
		this->initial_cont_pair        = {0,0,0, 0,0,0, 0,0};
		this->final_cont_monomer       = {0,0,0, 0,0,0, 0,0};
		this->final_cont_part2switch   = {0,0,0, 0,0,0, 0,0};
		this->final_cont_pair          = {0,0,0, 0,0,0, 0,0};
		return;
	}

};

#endif
