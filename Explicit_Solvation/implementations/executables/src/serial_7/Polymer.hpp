#ifndef _POLYMER_H_
#define _POLYMER_H_
#include "Particle.hpp"


class Polymer{
public:
	short deg_poly; 
	std::vector <Particle*> chain;                                                // all the particles in the polymer chain 
	std::map <Particle*, std::vector <Particle*> > ConnectivityMap;               // the set of beads that one particular polymer bead is connected to 


	// constructor 
	Polymer (short deg_poly_, std::vector <Particle*> particleChain):  deg_poly (deg_poly_), chain (particleChain) {
		// this->ChainToConnectivityMap(); 
		this->chain.reserve(deg_poly_); 
		this->chain_to_connectivity_map();
	}

	Polymer (std::vector<std::array <int,3>> locations, std::string type_m="m1"){
		std::vector <int> pmer_spins; 
		short size_ = locations.size(); 
		for (short i=0; i<size_; i++){
			// unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
			std::random_device rd;
			std::mt19937 generator(rd () ); 
			// std::cout << "check..." << std::endl;
			std::uniform_int_distribution<int> distribution (0, 25); 
			pmer_spins.push_back(distribution(generator));
			
		}

		std::vector <Particle*> ptc_vec; 

		for (int i=0;i< size_ ; i++ ){
			Particle* p = new Particle (locations.at(i), type_m, pmer_spins.at(i)); 
			ptc_vec.push_back(p); 
		}
		this->chain.reserve(ptc_vec.size());
		this->chain = ptc_vec;
		this->deg_poly = ptc_vec.size();
		this->chain_to_connectivity_map();
	}

	Polymer (std::vector<std::array<int,3>> locations, std::vector <int> spins, std::string type_m="m1"){
		std::vector <Particle*> ptc_vec; 
		short size_ = locations.size(); 

		for (int i=0; i < size_ ; i++ ){
			Particle* p = new Particle (locations.at(i), type_m, spins.at(i)); 
			ptc_vec.push_back(p); 
		}
		this->deg_poly = size_;
		this->chain.reserve(size_);
		this->chain    = ptc_vec; 
		this->chain_to_connectivity_map();
	}

	// destructor 
	~Polymer (){

	};

	// print positions of monomer present in the polymer
	void print_chain(); 

	// print orientation of monomer present in the polymer
	void print_orientation(); 

	// obtain connectivity map given the chain 
	void chain_to_connectivity_map(); 

};

#endif
