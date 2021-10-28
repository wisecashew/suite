#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_


// #########################################################
// This is my Particle class. 
// it has the following attributes: 
// 1. loc - location. This is a 3D vector holding the location of the particle. 
//
// it has the following methods:
// 1. print_loc - just prints out location of the molecule. Useful for debugging. 

class Particle{
public: 
	std::vector <int> loc; 
	std::string ptype; 
	// energetic interaction 
	// isotropy vs anisotropy

	bool operator<(const Particle& rhs )const{
		return loc < rhs.loc; 
	};

	// constructor 

	Particle(std::vector <int> loc, std::string s) {
		this->loc = loc; 
		this->ptype = s;
	}

	// destructor 
	~Particle(){
		//std::cout << "Particle has been erased from system memory." << std::endl;
	}

	// print current location of the particle 
	void print_loc(); 

	// given a particle, obtain its neighboring site locations 
	// std::vector <std::vector <int>> nlist (int x_len, int y_len, int z_len);


}; 

// END OF CLASS PARTICLE 
// #############################################################




// #############################################################
// This is my Polymer class 
// it has the following attributes:
// 1. chain - this is a std::vector of Particles. 
//
// it has the following methods:
  


class Polymer{
public:
	std::vector <Particle> chain;
	std::map <Particle, std::vector <Particle> > conn; 

	// constructor 

	// default constructor
	Polymer(){

	};

	Polymer(std::vector <Particle> chain){
		this->chain=chain;
	};
	

	// destructor 
	~Polymer(){
		//std::cout << "Polymer object has been destroyed from machine."<<std::endl;
	}

	// print position of monomers in grid 
	void print_loc(); 

	// get the list of all the points the polymer is occupying 
	// void get_plocs(); 

	// obtain the connectivity mapping 
	void obtain_connectivity(); 

	// find if there are kinks in the polymer structure 
	std::vector <int> find_kinks(); 

	// find if there are any cranks in the polymer structure 
	std::vector <int> find_cranks(); 

};

// END OF CLASS POLYMER
// ##############################################################


class Solvent{
public:
	std::vector <Particle> solvation; 

	// constructor

	Solvent(){

	}; 

	Solvent(std::vector <Particle> collection)
	: solvation { collection }{
		for (Particle p: collection){
			p.ptype = "solvent";
		}

	}; 

	//destructor 
	~Solvent(){

	};

	void print_loc(); 
};




// ##############################################################
// This is my Grid class. Probably the most important class. 
// it has the following attributes: 
// 1. GRID - this is a vector of particles. Quite like Polymer, but Grid has more 
// more functionality 
// 2. loc_list - questionable if I want this. This holds the locations of 
// particles in Grid. 
// 3. x_len - length of box in x direction 
// 4. y_len - length of box in y direction 
// 5. z_len - length of box in z direction 
//
// it has the following methods: 
// 1. polymer_insertion - inserts a polymer into the grid, with seed = (0,0,0). 

class Grid{
public:
	// attributes
	std::vector <Particle> occupied;				// all the particles in the grid 
	Polymer polymer; 								// all the polymers in the grid 
	Solvent solvent; 
	

	int x_len; 									// length of grid along x-axis 
	int y_len; 									// length of grid along y-axis
	int z_len; 									// length of grid along z-axis 

	// constructor 
	Grid(int xl, int yl, int zl){
		std::vector <Particle> g; 
		this->occupied = g; 
		this->x_len = xl; 
		this->y_len = yl; 
		this->z_len = zl; 
	
	};

	~Grid(){
	};

	// methods 

	// print out the coordinates of the polymer
	void print_polymer(); 

	// print out the coordinates of the solvent 
	void print_solvent(); 

	// update occupied as system is solvated or polymer is inserted
	void print_occupied();  


	// insert a polymer in the grid 
	void polymer_insertion(int dop, std::vector <int> seed); 

	// solvate the system 
	void solvate();

	// perform Monte Carlo moves on polymer 
	// ====================================
	// this code has a master function:
	void end_rotation(); 

	// and two child functions 
	void ZeroIndexRotation(); 
	void FinalIndexRotation(); 
	// ===================================

	// requires the use of find_kinks - part of class polymer 
	void kink_jump(); 
	// ===================================

	// requires the use of find_crank - part of class polymer 
	void crank_shaft(); 


	// perform Reptation moves on polymer 
	// ===================================
	// this code has a master function 
	void reptation(); 

	// requires two child functions
	void FinalToZero(); 
	void ZeroToFinal(); 
	// ====================================

	// calculate energy of the grid 
	int CalculateEnergy(); 





};
// END OF CLASS GRID
// ##############################################################


# endif // _MC_CLASSES_H_