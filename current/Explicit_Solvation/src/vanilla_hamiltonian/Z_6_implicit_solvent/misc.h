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
void modified_direction (std::array<int,3>* a, int x, int y, int z);
int  modified_modulo    (int divident, int divisor);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// imposing periodic boundary conditions 
void impose_pbc (std::vector <int>* vect , int x_len, int y_len, int z_len); 
void impose_pbc (std::array  <int,3>* arr, int x_len, int y_len, int z_len);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// obtain neighbor list 
std::vector <std::vector <int>>     obtain_ne_list (std::vector <int>   loc, int x_len, int y_len, int z_len); 
std::array  <std::array <int,3>, 6> obtain_ne_list (std::array  <int,3> loc, int x_len, int y_len, int z_len);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// random number generation 
int    rng_uniform (int start, int end);
double rng_uniform (double start, double end);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// std::vector and std::array arithmetic
// adding the two 
std::vector <int>   add_vectors (std::vector <int>* v1, std::vector <int>* v2); 
std::array  <int,3> add_arrays  (std::array <int,3>* a1, std::array <int,3>* a2);

// subtracting the two 
std::vector <int>   subtract_vectors (std::vector <int>* v1, std::vector <int>* v2); 
std::array  <int,3> subtract_arrays  (std::array <int,3>* a1, std::array <int,3>* a2);


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// check input of main driver code 
void InputParser (int dfreq, int max_iter, std::string positions, std::string topology, std::string dfile, std::string efile, std::string mfile, std::string stats_file); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// print methods  
// print out a vector 
void print (std::vector <int> v); 
void print (std::vector <double> v);

// print out an array 
void print (std::array <int,3> a); 
void print (std::array <double,3> a);

// print out a list of vectors
void print (std::vector <std::vector <int>> v); 
void print (std::vector <std::vector <double>> v); 
void print (std::vector <Particle> pvec);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// creation methods 
std::vector <Particle>           CreateSolventVector (int x, int y, int z, std::vector <Polymer>* Polymers); 
Polymer                          makePolymer         (std::vector <std::array <int,3>> locations, std::string type_m="monomer");
// create a lattice 
std::vector <std::array <int,3>> create_lattice_pts  (int x_len, int y_len, int z_len); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// validity checks 
bool checkValidityOfCoords                (std::array <int,3> v, int x, int y, int z); 
bool checkForOverlaps                     (std::vector <Polymer> Polymers);
bool checkConnectivity                    (std::vector <Polymer> Polymers, int x, int y, int z); 
bool checkForSolventMonomerOverlap        (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent); 
int  IsSolvent                            (std::vector <Polymer>* Polymers, std::array <int,3>* to_check);
int  SolventIndexReporter                 (std::vector <Particle>* Solvent, std::array <int,3>* to_check);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// extraction methods 
std::vector <Polymer>  ExtractPolymersFromFile    (std::string filename, int x, int y, int z);
std::vector <Polymer>  ExtractPolymersFromFile    (std::string filename);
std::vector <Polymer>  ExtractPolymersFromTraj    (std::string trajectory, std::string position);
int                    ExtractIndexOfFinalMove    (std::string trajectory);
double                 ExtractEnergyOfFinalMove   (std::string energy_file); 
int                    ExtractNumberOfPolymers    (std::string filename);
std::array <double, 8> ExtractTopologyFromFile    (std::string filename);
double                 NumberExtractor            (std::string s);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// energy calculator and metropolis 
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double CalculateEnergy      (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, double Emm_a, double Emm_n, double Ems_a, double Ems_n, double* mm_aligned, double* mm_naligned, int* ms_aligned, int* ms_naligned);
bool   MetropolisAcceptance (double E1, double E2, double kT); 
bool   MetropolisAcceptance (double E1, double E2, double kT, double rweight); 
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// dump methods 
void dumpPositionsOfPolymers (std::vector <Polymer> * PolymersInGrid, int step, std::string filename);
void dumpEnergy              (double sysEnergy, int step, double mm_aligned, double mm_naligned, int ms_aligned, int ms_naligned, std::string filename);
void dumpPositionOfSolvent   (std::vector <Particle>* Solvent, int step, std::string filename);
void dumpOrientation         (std::vector <Polymer> * Polymers, std::vector <Particle>* Solvent, int step, std::string filename, int x, int y, int z);
void dumpMoveStatistics      (std::array  <int,9>   * attempts, std::array <int,9>* acceptances, int step, std::string stats_file); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// Polymer moves

std::vector <Polymer> Translation          (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL); 
void                  OrientationFlip      (std::vector <Particle>* SolvVect, int x, int y, int z, int size_of_region);
void                  SolventFlip          (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, double* rweight);
void                  PolymerFlip          (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, double* rweight);
void                  PolymerFlipLocal     (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, double* rweight); 

// methods relevant to chain regrowth 
bool checkOccupancy                                     (std::array <int,3>* loc, std::vector <Polymer>* Polymers);
bool checkOccupancyTail                                 (std::array <int,3>* loc, std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer);
bool checkOccupancyHead                                 (std::array <int,3>* loc, std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer);
std::vector <std::array <int,3>> extract_positions_tail (std::vector <Particle>* chain, int pivot_idx);
std::vector <std::array <int,3>> extract_positions_head (std::vector <Particle>* chain, int pivot_idx);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// Important auxiliary function 
Particle                          ParticleReporter     (std::vector <Polymer>* Polymers, std::vector <Particle>* SolvVect, std::array <int,3> to_check);
std::array <std::array <int,3>,3> HingeSwingDirections (std::array <int,3>* HingeToHinge, std::array <int,3>* HingeToKink, int x, int y, int z); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// Polymer moves with Rosenbluth sampling methods 

void                  TailRotation_Rosenbluth      (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  HeadRotation_Rosenbluth      (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  EndRotation_Rosenbluth       (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  KinkJump_Rosenbluth          (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  CrankShaft_Rosenbluth        (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  ForwardReptation_Rosenbluth  (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  BackwardReptation_Rosenbluth (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  Reptation_Rosenbluth         (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight);
void                  ChainRegrowth_Rosenbluth     (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index_of_polymer, int x, int y, int z, bool* IMP_BOOL); 
void                  TailSpin_Rosenbluth          (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight); 
void                  HeadSpin_Rosenbluth          (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int deg_poly,int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight);
void                  MoveChooser_Rosenbluth       (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, bool v, bool* IMP_BOOL, double* rweight, std::array <int,9>* attempts, int* move_number, int Nsurr);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


// !~!~!~!!~!!~!~!~!!~~!!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~
// METHOD GRAVEYARD
// !~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~!~


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

#endif 
