#ifndef _MISC_H_
#define _MISC_H_
#include "classes.h"


// this is a bunch of miscellaneous functions I use 
// to manipulate std::vectors. Furthermore, this will be used 
// as a skeleton for building functions for grids 

// add two vectors 
std::vector <int> add_vectors(std::vector <int>* v1, std::vector <int>* v2); 
std::array <int,3> add_arrays(std::array <int,3>* a1, std::array <int,3>* a2);

// subtract two vectors 
std::vector <int> subtract_vectors(std::vector <int>* v1, std::vector <int>* v2); 
std::array <int,3> subtract_arrays(std::array <int,3>* a1, std::array <int,3>* a2);

// check if the new location added to random walk has been visited before 
bool check_avoidance(std::vector <int> to_check, std::vector<std::vector <int>> loc_list); 

// print out a vector 
void print(std::vector <int> v); 
void print(std::vector <double> v);

// print out an array 
void print(std::array <int,3> a); 
void print(std::array <double,3> a);

// print out a list of vectors
void print(std::vector <std::vector <int>> v); 
void print(std::vector <std::vector <double>> v); 
void print(std::vector <Particle> pvec);



// impose periodic boundary conditions on a vector 
void impose_pbc(std::vector <int>* vect, int x_len, int y_len, int z_len); 
void impose_pbc(std::array <int,3>* arr, int x_len, int y_len, int z_len);

// run a sarw without checking for pbc 
void sarw(std::vector<std::vector<int>>* loc_list, int dop); 

// run a sarw with periodic boundary conditions 
void sarw_pbc(std::vector<std::vector<int>>* loc_list, int dop, int x_len, int y_len, int z_len); 


// create a lattice 
std::vector <std::array <int,3>> create_lattice_pts(int x_len, int y_len, int z_len); 


// obtain the neighbor vector (list) 
// implicitly assumed in 3d 

std::vector <std::vector <int>> obtain_ne_list(std::vector <int> loc, int x_len, int y_len, int z_len); 
std::array <std::array <int,3>, 6> obtain_ne_list(std::array <int,3> loc, int x_len, int y_len, int z_len);

// given a list of locations for a particular type of particle, create a vector of particles 
std::vector <Particle> loc2part (std::vector <std::vector <int>> loc_list, std::string s); 

// given a vector of particles, obtain a vector of locations 
std::vector <std::vector <int> > part2loc (std::vector <Particle> pVec); 

// run an acceptance criterion for two polymers 
bool acceptance(int dE, double kT); 

// calculate energy of polymer-solvent interaction only 
int PolymerEnergySolvent(std::vector <Particle> polymer, int x_len, int y_len, int z_len, int intr_energy);


int PolymerEnergySolvation(std::vector <Particle> polymer, int x_len, int y_len, int z_len, int intr_energy, int intr_energymm);


// extract polymer information from file
int ExtractNumberOfPolymers(std::string filename);
std::vector <Polymer> ExtractPolymersFromFile(std::string filename);

std::array <std::array <int,3>,3> HingeSwingDirections(std::array <int,3>* HingeToHinge, std::array <int,3>* HingeToKink, int x, int y, int z); 

void modified_direction(std::array<int,3>* a, int x, int y, int z);
int modified_modulo(int divident, int divisor);


int rng_uniform(int start, int end);
double rng_uniform(double start, double end);

// extracting info from files 
std::array <double, 8> ExtractTopologyFromFile(std::string filename);
double NumberExtractor(std::string s);

// checking if info is accurate 
bool isSymmetric(std::vector <std::vector <double>> mat);

// making a polymer from a bunch of locations 
Polymer makePolymer(std::vector <std::array <int,3>> locations, std::string type_m="monomer");

// flipping the orientation of a bunch of particles 
void ClusterFlip(std::vector <Particle>*);

// metropolis acceptance criterion
bool MetropolisAcceptance(double E1, double E2, double kT); 

// sending a string to a file
void StringToFile(std::string filename, std::string to_send);

// extract number of polymers from the topology file 
int ExtractNumberOfPolymers(std::string filename);

// check input of main driver code 
void InputParser(bool a, bool r, int Nacc, int dfreq, int max_iter, 
    std::string positions, std::string topology, std::string dfile, 
    std::string efile, std::string restart_traj); 


// create the solvent vector 
std::vector <Particle> CreateSolventVector(int x, int y, int z, std::vector <Polymer>* PolymerVector); 

bool checkValidityOfCoords(std::array <int,3> v, int x, int y, int z); 
bool checkForOverlaps(std::vector <Polymer> PolymerVector);
bool checkConnectivity(std::vector <Polymer> PolymerVector, int x, int y, int z); 
bool checkForSolventMonomerOverlap(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolventVector); 
std::vector<Polymer> ExtractPolymersFromFile(std::string filename, int x, int y, int z);
std::vector <Polymer> ExtractPolymersFromTraj(std::string trajectory, std::string position);
double CalculateEnergy(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int x, int y, int z, double Emm_a, double Emm_n, double Ems_a, double Ems_n, double* m_neighbor, int* a_contacts, int* n_contacts);
void dumpPositionsOfPolymers (std::vector <Polymer>* PolymersInGrid, int step, std::string filename);
void dumpEnergy (double sysEnergy, int step, double m_contacts, int a_contacts, int n_contacts, std::string filename, bool first_call); 
std::vector <Polymer> TailRotation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL);
std::vector <Polymer> HeadRotation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> EndRotation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> KinkJump(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> CrankShaft(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL);
std::vector <Polymer> ForwardReptation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL);
std::vector <Polymer> BackwardReptation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> Reptation(std::vector<Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> Translation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> MoveChooser(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int x, int y, int z, bool v, bool* IMP_BOOL);
int ExtractIndexOfFinalMove(std::string trajectory);
double ExtractEnergyOfFinalMove(std::string energy_file); 
void dumpPositionOfSolvent(std::vector <Particle>* SolventVector, int step, std::string filename);

Particle ParticleReporter (std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVect, std::array <int,3> to_check);

void ChainRegrowth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVect, int index_of_polymer, int x, int y, int z, bool* IMP_BOOL );
void TailSpin(std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* b, bool* IMP_BOOL, bool* first_entry_bool );
void HeadSpin(std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer, int deg_poly,int x, int y, int z, bool* b, bool* IMP_BOOL, bool* first_entry_bool);
bool checkOccupancy(std::array <int,3>* loc, std::vector <Polymer>* PVec);
bool checkOccupancyTail(std::array <int,3>* loc, std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer);
bool checkOccupancyHead(std::array <int,3>* loc, std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer);

std::vector <std::array <int,3>> extract_positions_tail (std::vector <Particle>* chain, int pivot_idx);
std::vector <std::array <int,3>> extract_positions_head (std::vector <Particle>* chain, int pivot_idx);
void OrientationFlip (std::vector <Particle>* SolvVect, int x, int y, int z, int size_of_region);

#endif 
