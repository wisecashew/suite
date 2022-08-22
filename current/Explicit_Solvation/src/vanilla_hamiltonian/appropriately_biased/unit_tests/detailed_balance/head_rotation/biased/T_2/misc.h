#ifndef _MISC_H_
#define _MISC_H_
#include "classes.h"


// this is a bunch of miscellaneous functions I use 
// to manipulate std::vectors. Furthermore, this will be used 
// as a skeleton for building functions for grids 


// check if the new location added to random walk has been visited before 
bool check_avoidance             (std::vector <int> to_check, std::vector<std::vector <int>> loc_list); 

// run a sarw without checking for pbc 
void sarw                        (std::vector<std::vector<int>>* loc_list, int dop); 

// run a sarw with periodic boundary conditions 
void sarw_pbc                    (std::vector<std::vector<int>>* loc_list, int dop, int x_len, int y_len, int z_len); 


// given a list of locations for a particular type of particle, create a vector of particles 
std::vector <Particle> loc2part           (std::vector <std::vector <int>> loc_list, std::string s); 

// given a vector of particles, obtain a vector of locations 
std::vector <std::vector <int> > part2loc (std::vector <Particle> pVec); 

// run an acceptance criterion for two polymers 
bool acceptance(int dE, double kT); 

// calculate energy of polymer-solvent interaction only 
int PolymerEnergySolvent   (std::vector <Particle> polymer, int x_len, int y_len, int z_len, int intr_energy);
int PolymerEnergySolvation (std::vector <Particle> polymer, int x_len, int y_len, int z_len, int intr_energy, int intr_energymm);

// checking if info is accurate 
bool isSymmetric(std::vector <std::vector <double>> mat);


// flipping the orientation of a bunch of particles 
void ClusterFlip(std::vector <Particle>*);

// sending a string to a file
void StringToFile(std::string filename, std::string to_send);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// imposing modified modulo methods 
void   modified_direction (std::array<int,3>* a, int x, int y, int z);
int    modified_modulo    (int divident, int divisor);
double modified_modulo    (double divident, int divisor);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// imposing periodic boundary conditions 
void impose_pbc (std::vector <int>* vect , int x_len, int y_len, int z_len); 
void impose_pbc (std::array  <int,3>* arr, int x_len, int y_len, int z_len);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// obtain neighbor list 
std::vector <std::vector <int>>      obtain_ne_list (std::vector <int>   loc, int x_len, int y_len, int z_len); 
std::array  <std::array <int,3>, 26> obtain_ne_list (std::array  <int,3> loc, int x_len, int y_len, int z_len);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// random number generation 
int    rng_uniform (int start, int end);
double rng_uniform (double start, double end);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
int                  lattice_index ( std::array<int,3> location, int y, int z );
std::array <int,3>   location      ( int lattice_index, int x, int y, int z);


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// std::vector and std::array arithmetic
// adding the two 
std::vector <int>      add_vectors (std::vector <int>* v1, std::vector <int>* v2); 
std::array  <int,3>    add_arrays  (std::array <int,3>* a1, std::array <int,3>* a2);
std::array  <double,3> add_arrays  (std::array <int,3>* a1, std::array <double,3>* a2);
std::array  <double,3> add_arrays  (std::array <double,3>* a1, std::array <double,3>* a2);

// subtracting the two 
std::vector <int>      subtract_vectors (std::vector <int>* v1, std::vector <int>* v2); 
std::array  <int,3>    subtract_arrays  (std::array <int,3>* a1, std::array <int,3>* a2);
std::array  <double,3> subtract_arrays  (std::array <double,3>* a1, std::array <double,3>* a2);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// getting distance between two points 
double              distance_between_points (std::array <int,3>*    a1, std::array <int,3>* a2   , int xlen, int ylen, int zlen);
double              distance_between_points (std::array <double,3>* a1, std::array <double,3>* a2, int xlen, int ylen, int zlen);
double              take_dot_product        (int o1, int o2); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// check input of main driver code 
void InputParser (int dfreq, int max_iter, bool r, std::string positions, \
    std::string topology, std::string dfile, std::string efile, std::string mfile, \
    std::string stats_file, std::string lattice_file_read); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// print methods  
// print out a vector 

template <typename T>
void print (T vec, std::string end="\n"){

    for (int i{0}; i < static_cast<int>(vec.size()); ++i){
        std::cout << vec[i] << " | "; 
    }
    std::cout << end; 
    return;
}

void print (std::vector <int> v, std::string c="\n"); 
void print (std::vector <double> v, std::string c="\n");

// print out an array 
void print (std::array <int,3> a, std::string c="\n"); 
void print (std::array <int,6> a, std::string c="\n"); 
void print (std::array <double,3> a, std::string c="\n");
void print (std::array <double,4> a, std::string c="\n"); 
void print (std::array <double, 6> v, std::string c="\n");
void print (std::array <std::array<int,3>,6> aa );

// print out a list of vectors
void print (std::vector <std::vector <int>> v); 
void print (std::vector <std::vector <double>> v); 
void print (std::vector <std::vector <std::array<int,3>>> v);
void print (std::vector <Particle> pvec);
void print (std::map<std::array<int,3>,Particle*> LATTICE);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// moves to check locations in *Polymers and *Solvent 
bool MonomerReporter         (std::vector <Polymer>* Polymers  , std::array <int,3>* to_check);
bool MonomerReporter         (std::vector <Particle*>* LATTICE , std::array <int,3>* to_check  , int y, int z);
bool MonomerReporter         (std::vector <Particle*>* Polymers, std::array <int,3>* to_check_1, std::array <int,3>* to_check_2, int y, int z);
bool MonomerReporter         (std::vector <Particle*>* LATTICE , std::array <int,3>* to_check_1, std::array <int,3>* to_check_2, int y, int z);
bool MonomerNeighborReporter (std::vector <Polymer>* Polymers  , std::array <int,3>* to_check  , int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// creation methods 
std::vector <Particle>           CreateSolventVector     (int x, int y, int z, std::vector <Polymer>* Polymers); 
void                             AddSolvent              (int x, int y, int z, std::vector <Particle*>* LATTICE);
void                             SetUpLatticeFromScratch (int x, int y, int z, std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, std::string positions);
void                             SetUpLatticeFromRestart (int x, int y, int z, std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, int step_number, std::string lattice_file_read, std::string dfile, std::string positions); 
Polymer                          makePolymer             (std::vector <std::array <int,3>> locations, std::string type_m="m1"); 
Polymer                          makePolymer             (std::vector <std::array <int,3>> locations, std::vector<int> pmer_spins, std::string type_m="m1"); 
// create a lattice 
std::vector <std::array <int,3>> create_lattice_pts  (int x_len, int y_len, int z_len); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// validity checks 
void CheckStructures                      (int x, int y, int z, std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE);
bool checkValidityOfCoords                (std::array <int,3> v, int x, int y, int z); 
bool checkForOverlaps                     (std::vector <Polymer> Polymers);
bool checkForOverlaps                     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE);
bool checkConnectivity                    (std::vector <Polymer> Polymers, int x, int y, int z); 
bool checkForSolventMonomerOverlap        (std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, int y, int z); 
bool checkPointersOnLattice               (std::vector <Particle*>* LATTICE, int x, int y, int z);
int  IsSolvent                            (std::vector <Polymer>* Polymers, std::array <int,3>* to_check);
int  SolventIndexReporter                 (std::vector <Particle>* Solvent, std::array <int,3>* to_check);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// extraction methods 
std::vector <Polymer>    ExtractPolymersFromFile    (std::string filename, int x, int y, int z);
std::vector <Polymer>    ExtractPolymersFromFile    (std::string filename);
std::vector <Polymer>    ExtractPolymersFromTraj    (std::string trajectory, std::string position, int step_number, int x, int y, int z);
int                      ExtractIndexOfFinalMove    (std::string trajectory);
double                   ExtractEnergyOfFinalMove   (std::string energy_file); 
int                      ExtractNumberOfPolymers    (std::string filename);
std::array <double, 8>   ExtractTopologyFromFile    (std::string filename);
std::vector <Particle*>  ExtractLatticeFromRestart  (std::string rfile, int* step_num, int x, int y, int z);
double                   NumberExtractor            (std::string s);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// energy calculator and metropolis 
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double CalculateEnergy      (std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, std::array<double,4>* E, std::array<double,4>* contacts, int x, int y, int z);
bool   MetropolisAcceptance (double E1, double E2, double kT); 
bool   MetropolisAcceptance (double E1, double E2, double kT, double rweight); 
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// dump methods 
void dumpPositionsOfPolymers (std::vector <Polymer> * PolymersInGrid, int step, std::string filename);
void dumpEnergy              (double sysEnergy, int step, std::array<double,4>* contacts, std::string filename);
void dumpPositionOfSolvent   (std::vector <Particle*>* LATTICE, int step, std::string filename);
void dumpOrientation         (std::vector <Polymer> * Polymers, std::vector<Particle*>* LATTICE, int step, std::string filename, int x, int y, int z);
void dumpMoveStatistics      (std::array  <int,9>   * attempts, std::array <int,9>* acceptances, int step, std::string stats_file); 
void dumpLATTICE             (std::vector <Particle*>* LATTICE, int step, int y, int z, std::string filename);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// Polymer moves

std::vector <Polymer> Translation          (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL); 
void                  SolventFlip          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z, double* rweight, int Nsurr, std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory ); 
void                  SolventFlipSingular  (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z, double* rweight, int Nsurr, std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory );
void                  PolymerFlip          (std::vector <Polymer>* Polymers, double* rweight, int Nsurr, std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory ); 
void                  PolymerFlipSingular  (std::vector <Polymer>* Polymers, double* rweight, int Nsurr, std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory ); 
void                  SiteFlipSingular     (std::vector <Particle*>* LATTICE, int x, int y, int z, std::pair <std::vector <std::array<int,2>>, std::vector <std::array<int,2>>>* memory ) ;

// methods relevant to chain regrowth 
bool checkOccupancy                                     (std::array <int,3>* loc, std::vector <Polymer>* Polymers);
bool checkOccupancyTail                                 (std::array <int,3>* loc, std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer);
bool checkOccupancyHead                                 (std::array <int,3>* loc, std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer);
std::vector <std::array <int,3>> extract_positions_tail (std::vector <Particle*>* chain, int pivot_idx);
std::vector <std::array <int,3>> extract_positions_head (std::vector <Particle*>* chain, int pivot_idx);
void                             create_linked_list     (std::vector<std::array<int,3>> v1, std::vector<std::array<int,3>> v2, std::vector <std::array<int,3>> link, std::vector <std::vector <std::array<int,3>>>* master_linked_list, int beginning);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// Important auxiliary function 
Particle                          ParticleReporter     (std::vector <Polymer>* Polymers, std::vector <Particle>* SolvVect, std::array <int,3> to_check);
void                              ParticleReporter     (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, std::pair <std::string, int>* properties, std::array <int,3>* to_check);
std::array <std::array <int,3>,3> HingeSwingDirections (std::array <int,3>* HingeToHinge, std::array <int,3>* HingeToKink, int x, int y, int z); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// Polymer moves with Rosenbluth sampling methods 
void                  TailRotation 				   (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  HeadRotation                 (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  EndRotation                  (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  KinkJump                     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  BondVibration                (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  CrankShaft                   (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  ForwardReptation 			   (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  BackwardReptation 		   (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  Reptation  			       (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory);
void                  ChainRegrowth			       (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory, int* monomer_index, int* back_or_front); 
void                  TailSpin			           (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight); 
void                  HeadSpin			           (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int deg_poly,int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight);
void                  PerturbSystem                (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z, bool v, bool* IMP_BOOL, double* rweight, std::array <int,9>* attempts, int* move_number, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory3, std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory2, int* monomer_index, int* back_or_front, int Nsurr);
void                  ReversePerturbation          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int y, int z, bool v, int move_number, std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory3, std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory2, int monomer_index, int back_or_front);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  TailRotation_SIMPLE          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  HeadRotation_SIMPLE          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  EndRotation_SIMPLE           (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  ForwardReptation_SIMPLE      (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  BackwardReptation_SIMPLE     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  Reptation_SIMPLE             (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  ChainRegrowth                (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double* prob_o_to_n, double temperature, int index, int x, int y, int z); 
void                  HeadRegrowth                 (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  TailRegrowth                 (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  BackFlowFromHeadRegrowth     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::vector <std::array <int,3>>* old_cut, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z); 
void                  BackFlowFromTailRegrowth     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::vector <std::array <int,3>>* old_cut, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z); 
void                  PerturbSystem                (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, std::array <int,9>* attempts, bool* IMP_BOOL, bool v, double* sysEnergy, double* prob_o_to_n, double temperature, int* move_number, int x, int y, int z);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  TailRegrowth_UNBIASED        (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, bool* IMP_BOOL, int deg_poly, int p_index, int m_index, int x, int y, int z);
void                  HeadRegrowth_UNBIASED        (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, bool* IMP_BOOL, int deg_poly, int p_index, int m_index, int x, int y, int z);
void                  ChainRegrowth_UNBIASED       (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  PerturbSystem_UNBIASED       (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, std::array <int,9>* attempts, bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, int* move_number, int x, int y, int z);

void                  TailRotation_BIASED          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  HeadRotation_BIASED          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  EndRotation_BIASED           (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  PerturbSystem_BIASED         (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,4>* E, std::array <double,4>* contacts, std::array<int,9>* attempts, bool* IMP_BOOL, bool v, double* sysEnergy, double temperature,  int* move_number, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  forward_reptation_with_tail_biting        (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <int,3>* to_slither, int deg_poly, int index, int y, int z);
void                  backward_reptation_with_head_butting      (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <int,3>* to_slither, int deg_poly, int index, int y, int z);
void                  forward_reptation_without_tail_biting     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* loc0, int deg_poly, int index, int y, int z);
void                  backward_reptation_without_head_butting   (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* locf, int deg_poly, int index, int y, int z);


void        forward_reptation_without_tail_biting  (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::vector <double>* energies, std::vector <std::array<int,4>>* contacts_store, std::array <double,4>* E, std::array <double,4>* c_contacts, std::array<int,3> to_slither, std::array <int,3> loc0, std::array <int,3> locf, int deg_poly, int index, int x, int y, int z);
void        backward_reptation_without_head_biting (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::vector <double>* energies, std::vector <std::array<int,4>>* contacts_store, std::array <double,4>* E, std::array <double,4>* c_contacts, std::array<int,3> to_slither, std::array <int,3> loc0, std::array <int,3> locf, int deg_poly, int index, int x, int y, int z);
void        forward_reptation_with_tail_biting     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::vector <double>* energies, std::vector <std::array<int,4>>* contacts_store, std::array <double,4>* E, std::array <double,4>* c_contacts, std::array <int,3> to_slither, int deg_poly, int index, int x, int y, int z); 
void        backward_reptation_with_head_biting    (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::vector <double>* energies, std::vector <std::array<int,4>>* contacts_store, std::array <double,4>* E, std::array <double,4>* c_contacts, std::array <int,3> to_slither, int deg_poly, int index, int x, int y, int z);


template <typename T>
void reset(T &x){
    x = T();
}


std::array <double,3>  scale_arrays ( double scalar, std::array <double,3>* array );
std::array <double,3>  scale_arrays ( double scalar, std::array <int,3>*    array ); 


// !~!~!~!!~!!~!~!~!!~~!!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
// METHOD GRAVEYARD
// !~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~

/*
std::vector <Polymer> TailRotation         (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL);
std::vector <Polymer> HeadRotation         (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> EndRotation          (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> KinkJump             (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> CrankShaft           (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL);
std::vector <Polymer> ForwardReptation     (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL);
std::vector <Polymer> Reptation            (std::vector<Polymer>*  Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL); 
std::vector <Polymer> BackwardReptation    (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL); 
void                  ChainRegrowth        (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index_of_polymer, int x, int y, int z, bool* IMP_BOOL );
void                  TailSpin             (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* b, bool* IMP_BOOL, bool* first_entry_bool );
void                  HeadSpin             (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int deg_poly,int x, int y, int z, bool* b, bool* IMP_BOOL, bool* first_entry_bool);
std::vector <Polymer> MoveChooser          (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, bool v, bool* IMP_BOOL);
*/ 

#endif 
