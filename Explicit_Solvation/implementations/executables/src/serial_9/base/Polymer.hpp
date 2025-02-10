#ifndef _POLYMER_H_
#define _POLYMER_H_
#include "Particle.hpp"


class Polymer{
public:
	int deg_poly; 
	std::vector <Particle*> chain;                                                // all the particles in the polymer chain 
	std::map <Particle*, std::vector <Particle*>> ConnectivityMap;               // the set of beads that one particular polymer bead is connected to 


	// constructor 1
	Polymer (int deg_poly_, std::vector <Particle*> particleChain):  deg_poly (deg_poly_){
		this->chain.reserve(deg_poly_); 
		this->chain = particleChain;
		this->chain_to_connectivity_map();
	}

	// constructor 2
	Polymer (std::vector<std::array <int,3>> locations, std::string type_m="m1"){
		std::vector <int> pmer_spins; 
		int size_ = static_cast<int>(locations.size()); 
		for (int i=0; i<size_; i++){
			// unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
			std::random_device rd;
			std::mt19937 generator(rd()); 

			// std::cout << "check..." << std::endl;
			std::uniform_int_distribution<int> distribution (0, 25); 
			pmer_spins.push_back(distribution(generator));
		}

		// define space in the chain
		this->chain.reserve(size_);

		// start throwing in particles in the chain
		for (int i=0;i < size_ ; i++ ){
			Particle* p = new Particle (pmer_spins.at(i), type_m, locations.at(i)); 
			this->chain.push_back(p); 
		}
		
		// start defining other things
		this->deg_poly = this->chain.size();
		this->chain_to_connectivity_map();
	}

	// constructor 3
	Polymer (std::vector<std::array<int,3>> locations, std::vector <int> spins, std::string type_m="m1"){

		int size_ = static_cast<int>(locations.size()); 
		this->chain.reserve(size_);

		for (int i=0; i < size_ ; i++ ){
			Particle* p = new Particle (spins.at(i), type_m, locations.at(i)); 
			this->chain.push_back(p); 
		}

		this->deg_poly = size_;
		this->chain_to_connectivity_map();
	}

	// destructor 
	~Polymer (){

	};

	// print functinos
	void print_chain();        // print coordinates of the chain 
	void print_orientation();  // print orientations of the chain

	// obtain connectivity map given the chain 
	void chain_to_connectivity_map(); 

};

#endif
