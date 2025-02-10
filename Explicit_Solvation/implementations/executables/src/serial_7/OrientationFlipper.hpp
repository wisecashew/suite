#ifndef _ORIENTATIONFLIPPER_H_
#define _ORIENTATIONFLIPPER_H_
#include "misc.hpp"

class EnhancedOrientationFlipper {
public:

	// properties
	// defining how particles were and how they were changed, and how many of them were
	int nparticles;                              // number of particles to flip
	std::vector <int> initial_orientations;      // these are the initial orientations of each particle
	std::vector <int> final_orientations;        // these are the final orientations of each particle

	// defining how new states will be sampled per particles
	int ntest;                                         // number of configurations to test
	std::vector <int>    orientations;                 // this is the container holding all the orientations of the trialed states
	std::vector <double> boltzmann;                    // this is the container holding the boltzmann factors of each perturbation
	
	// holders for energies as the orientation of each particle is perturbed
	double initial_E;              // energy of the system at the start
	double perturbed_E;            // energy of the system after perturbation
	std::vector <double> energies; // this is the container holding energies of each perturbation

	// containers for magnetization as the orientation of each particle is perturbed
	std::array <double,3>              initial_magnetization;   // initial contacts, prior to perturbation
	std::array <double,3>              perturbed_magnetization; // local change in contacts on perturbation
	std::vector <std::array<double,3>> magnetization_store;     // this is the container holding all the net contacts when running through the trialed states

	// containers for contacts as the orientation of each particle is perturbed
	std::array <double,CONTACT_SIZE>              initial_contacts;   // initial contacts, prior to perturbation
	std::array <double,CONTACT_SIZE>              perturbed_contacts; // local change in contacts on perturbation
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

	// constructor 
	EnhancedOrientationFlipper(){};

	// destructor
	~EnhancedOrientationFlipper(){};

	// set up the object
	void setup(int nparticles, int ntest){
		this->nparticles = nparticles;
		this->ntest      = ntest;
		this->current_contacts.fill(0);
		this->initial_contacts.fill(0);
		this->perturbed_contacts.fill(0);
		this->initial_orientations.resize(nparticles,0);
		this->final_orientations.resize(nparticles,0);
		this->boltzmann.resize(ntest,0.0);
		this->orientations.resize(ntest,0);
		this->energies.resize(ntest,0);
		this->contacts_store.resize(ntest, this->current_contacts);
		this->rboltzmann       = 0;
		this->prob_o_to_n      = 1;
		this->prob_n_to_o      = 1;
		this->Emin             = 0;
		this->current_energy   = 0;
		this->sampler_idx      = 0;
		this->sampler_rng      = 0;
		this->sampler_rsum     = 0;
		return;
	}

	// reset the object
	void reset(){
		std::fill(this->initial_orientations.begin(), this->initial_orientations.end(), 0);
		std::fill(this->final_orientations.begin(), this->final_orientations.end(), 0); 
		std::fill(this->boltzmann.begin(), this->boltzmann.end(), 0);
		std::fill(this->orientations.begin(), this->orientations.end(), 0);
		std::fill(this->orientations.begin(), this->orientations.end(), 0); 
		std::fill(this->energies.begin(), this->energies.end(), 0);
		this->rboltzmann       = 0;
		this->prob_o_to_n      = 1;
		this->prob_n_to_o      = 1;
		this->Emin             = 0;
		this->initial_E          = 0;
		this->perturbed_E        = 0;
		this->current_energy   = 0;
		this->sampler_idx      = 0;
		this->sampler_rng      = 0;
		this->sampler_rsum     = 0;
		this->current_contacts.fill(0);
		this->initial_contacts.fill(0);
		this->perturbed_contacts.fill(0);
		this->initial_magnetization.fill(0);
		this->perturbed_magnetization.fill(0);
		std::fill(this->contacts_store.begin(), this->contacts_store.end(), this->current_contacts);
		std::fill(this->magnetization_store.begin(), this->magnetization_store.end(), std::array<double,3>{0, 0, 0});
		return;
	}

	// more useful
	void reset(int nparticles, int ntest){
		this->nparticles = nparticles;
		this->ntest      = ntest;
		this->initial_orientations.resize(nparticles,0);
		this->final_orientations.resize(nparticles,0);
		this->current_contacts.fill(0);
		this->initial_contacts.fill(0);
		this->perturbed_contacts.fill(0);
		this->initial_magnetization   = {0,0,0};
		this->perturbed_magnetization = {0,0,0};
		this->boltzmann.resize(ntest,0.0);
		this->orientations.resize(ntest,0);
		this->energies.resize(ntest,0);
		this->contacts_store.resize(ntest, this->initial_contacts);
		this->magnetization_store.resize(ntest, std::array<double,3>{0,0,0});
		this->rboltzmann         = 0;
		this->prob_o_to_n        = 1;
		this->prob_n_to_o        = 1;
		this->Emin               = 0;
		this->initial_E          = 0;
		this->perturbed_E        = 0;
		this->current_energy     = 0;
		this->sampler_idx        = 0;
		this->sampler_rng        = 0;
		this->sampler_rsum       = 0;

		return;
	};

};

#endif
