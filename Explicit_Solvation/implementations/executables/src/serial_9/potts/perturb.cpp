#include "Potts.hpp"

void Potts::perturb_lattice_flip(){

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

	// set up the rng
	double rng = rng_uniform(0.0, 1.0);
	if (rng < std::exp(-1/this->T * (Ef - this->Energy))){
		this->Energy       = Ef;
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

