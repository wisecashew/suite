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


// given a list of locations for a particular type of Particle, create a vector of Particles 
std::vector <Particle> loc2part           (std::vector <std::vector <int>> loc_list, std::string s); 

// given a vector of Particles, obtain a vector of locations 
std::vector <std::vector <int> > part2loc (std::vector <Particle> pVec); 

// run an acceptance criterion for two polymers 
bool acceptance(int dE, double kT); 

// sending a string to a file
void StringToFile(std::string filename, std::string to_send);

// splitting up a string into components separated by a delimiter
std::vector<std::string> split (const std::string &s, char delim);

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
std::array  <double,8> add_arrays  (std::array <double,8>* a1, std::array <double,8>* a2);
std::array <double,8>  add_arrays (std::array<double,8> a1, std::array <double,8> a2); 


// subtracting the two 
std::vector <int>      subtract_vectors (std::vector <int>* v1, std::vector <int>* v2); 
std::array  <int,3>    subtract_arrays  (std::array <int,3>* a1, std::array <int,3>* a2);
std::array  <int,8>    subtract_arrays  (std::array <int,8>* a1, std::array <int,8>* a2);
std::array  <double,3> subtract_arrays  (std::array <double,3>* a1, std::array <double,3>* a2);
std::array  <double,8> subtract_arrays  (std::array <double,8>* a1, std::array <double,8>* a2);
std::array  <double,8> subtract_arrays  (std::array<double,8> a1, std::array <double,8> a2); 

// scaling arrays 
std::array <double,3>  scale_arrays ( double scalar, std::array <double,3>* array );
std::array <double,3>  scale_arrays ( double scalar, std::array <int,3>*    array ); 


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// getting distance between two points 
double              distance_between_points (std::array <int,3>*    a1, std::array <int,3>* a2   , int xlen, int ylen, int zlen);
double              distance_between_points (std::array <double,3>* a1, std::array <double,3>* a2, int xlen, int ylen, int zlen);
double              take_dot_product        (int o1, int o2); 
double              take_dot_product        (std::array <double,3> o1, std::array <double,3> o2); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// check input of main driver code 
void InputParser (int dfreq, int lfreq, int max_iter, bool r, std::string positions, \
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
void                             AddSolvent              (std::vector <Particle*>* LATTICE, int x, int y, int z); 
void                             AddCosolvent            (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, bool solvation, double frac, int dop, int x, int y, int z);
void                             SetUpLatticeFromScratch (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector<Particle*>* LATTICE, std::string positions, bool solvation, double frac, int x, int y, int z); 
void                             SetUpLatticeFromRestart (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector<Particle*>* LATTICE, int* step_number, std::string lattice_file_read, std::string dfile, std::string positions, int x, int y, int z); 
Polymer                          makePolymer             (std::vector <std::array <int,3>> locations, std::string type_m="m1"); 
Polymer                          makePolymer             (std::vector <std::array <int,3>> locations, std::vector<int> pmer_spins, std::string type_m="m1"); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// validity checks 
void CheckStructures                      (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector<Particle*>* LATTICE, int x, int y, int z); 
bool checkValidityOfCoords                (std::array <int,3> v, int x, int y, int z); 
bool checkForOverlaps                     (std::vector <Polymer> Polymers);
bool checkForOverlaps                     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE);
bool checkConnectivity                    (std::vector <Polymer> Polymers, int x, int y, int z); 
bool checkForSolventMonomerOverlap        (std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, int y, int z); 
bool checkPointersOnLattice               (std::vector <Particle*>* LATTICE, int x, int y, int z);
bool checkSolvationShells                 (std::vector <Polymer> Polymers, int y, int z); 
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
std::array <double, 13>  ExtractTopologyFromFile    (std::string filename, std::map <std::pair <std::string, std::string>, std::tuple<std::string, double, double, int, int>>* InteractionMap); 
std::vector <Particle*>  ExtractLatticeFromRestart  (std::string rfile, int* step_num, int x, int y, int z);
double                   NumberExtractor            (std::string s);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// energy calculator and metropolis 
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double CalculateEnergy               (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array<double,8>* E, std::array<double,8>* contacts, int x, int y, int z);
double CalculateEnergy_parallel      (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array<double,8>* E, std::array<double,8>* contacts, int x, int y, int z);
double NeighborEnergy                (std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, int ss_index, int x, int y, int z); 
double PairEnergy                    (Particle* p1, Particle* p2, std::array <double,8>* E, int* c_idx, int x, int y, int z);

// energy calculation 
double CalculateEnergyRevamped                 (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array<double,8>* contacts, int x, int y, int z); 
double NeighborEnergetics                      (std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, int ss_index, std::array <double,8>* contacts, int x, int y, int z); 
double IsolatedPairParticleInteraction         (Particle* p1, Particle* p2, std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, int* c_idx, int x, int y, int z); 
void   ParticlePairEnergyContribution          (Particle* p1, Particle* p2, std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, double* pair_energy, std::array <double,8>* contacts, int x, int y, int z); 
void   PairInteractionForRevampedEnergyCalc    (Particle* p1, Particle* p2, std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, double* pair_energy, std::array <double,8>* contacts, int x, int y, int z);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bool   MetropolisAcceptance          (double E1, double E2, double kT); 
bool   MetropolisAcceptance          (double E1, double E2, double kT, double rweight); 

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// dump methods 
void dumpPositionsOfPolymers (std::vector <Polymer> * PolymersInGrid, int step, std::string filename);
void dumpEnergy              (double sysEnergy, int step, std::array<double,8>* contacts, std::string filename);
void dumpPositionOfSolvent   (std::vector <Particle*>* LATTICE, int step, std::string filename);
void dumpOrientation         (std::vector <Polymer> * Polymers, std::vector<Particle*>* LATTICE, int step, std::string filename, int x, int y, int z);
void dumpMoveStatistics      (std::array  <int,9>   * attempts, std::array <int,9>* acceptances, int step, std::string stats_file); 
void dumpLATTICE             (std::vector <Particle*>* LATTICE, int step, int y, int z, std::string filename);
void dumpSolvation           (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int step, std::string filename, int x, int y, int z );
void BiasTheStart            (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z);
void AlignTheSolvationShell  (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z);
void AlignTheLattice         (std::vector <Particle*>* LATTICE); 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
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
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// Orientation flip moves 

void SolventFlip_UNBIASED_debug   (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void SolventFlip_UNBIASED         (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void PolymerFlip_UNBIASED            (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array<double,4>* E, std::array<double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void SolventExchange_BIASED          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z);
void SolventExchange_BIASED_debug    (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z);
void SolventExchange_UNBIASED_debug  (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z);
//
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// End rotation moves

void                  TailRotation_UNBIASED_debug        (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  HeadRotation_UNBIASED_debug        (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  EndRotation_UNBIASED_debug         (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  TailRotation_UNBIASED        (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  HeadRotation_UNBIASED        (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  EndRotation_UNBIASED         (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  TailRotation_BIASED          (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  HeadRotation_BIASED          (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  EndRotation_BIASED           (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

// necessary functions for reptation 
void                  forward_reptation_with_tail_biting            (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <int,3>* to_slither, int deg_poly, int index, int y, int z);
void                  forward_reptation_with_tail_biting_new        (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, double* frontflow_energy, int deg_poly, int index, int x, int y, int z);

void                  backward_reptation_with_head_butting          (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <int,3>* to_slither, int deg_poly, int index, int y, int z);
void                  backward_reptation_with_head_butting_new      (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, double* frontflow_energy, int deg_poly, int index, int x, int y, int z);

void                  forward_reptation_without_tail_biting         (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* loc0, int deg_poly, int index, int y, int z);
void                  forward_reptation_without_tail_biting_new     (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, std::array <int,3>* to_slither, double* frontflow_energy, int deg_poly, int index, int x, int y, int z);

void                  backward_reptation_without_head_butting       (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* locf, int deg_poly, int index, int y, int z);
void                  backward_reptation_without_head_butting_new   (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, std::array <int,3>* to_slither, double* frontflow_energy, int deg_poly, int index, int x, int y, int z);

// necessary perform reptation 
void                  ForwardReptation_UNBIASED_debug           (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  BackwardReptation_UNBIASED_debug          (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  Reptation_UNBIASED_debug                  (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  ForwardReptation_UNBIASED                 (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  BackwardReptation_UNBIASED                (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  Reptation_UNBIASED                        (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  acceptance_after_head_regrowth      (std::vector <Particle*>* LATTICE, std::vector <std::array<int,3>>* old_cut, std::vector <std::array<int,3>>* new_cut, int y, int z); 
void                  acceptance_after_tail_regrowth      (std::vector <Particle*>* LATTICE, std::vector <std::array<int,3>>* old_cut, std::vector <std::array<int,3>>* new_cut, int y, int z); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  TailRegrowth_UNBIASED              (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, bool* IMP_BOOL, int deg_poly, int p_index, int m_index, int x, int y, int z);
void                  HeadRegrowth_UNBIASED              (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, bool* IMP_BOOL, int deg_poly, int p_index, int m_index, int x, int y, int z);
void                  ChainRegrowth_UNBIASED             (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int p_index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  TailRegrowth_BIASED_debug                 (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  HeadRegrowth_BIASED_debug                 (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  BackFlowFromTailRegrowth_BIASED_debug     (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::vector <std::array <int,3>>* old_cut, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z); 
void                  BackFlowFromHeadRegrowth_BIASED_debug     (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::vector <std::array <int,3>>* old_cut, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z); 
void                  ChainRegrowth_BIASED_debug                (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  TailRegrowth_BIASED                       (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  HeadRegrowth_BIASED                       (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  BackFlowFromTailRegrowth_BIASED           (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::vector <std::array <int,3>>* old_cut, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z); 
void                  BackFlowFromHeadRegrowth_BIASED           (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::vector <std::array <int,3>>* old_cut, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z); 
void                  ChainRegrowth_BIASED                      (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  SolvationShellFlip_BIASED_remake2         (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z );
void                  SolvationShellFlip_BIASED_remake2_debug   (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z );
void                  PolymerFlip_BIASED                        (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  PolymerFlip_BIASED                        (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);
void                  PolymerFlip_BIASED_debug                  (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  ChainRegrowthPlusOrientationFlip_BIASED             (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int p_index, int x, int y, int z); 
void                  HeadRegrowthPlusOrientationFlip_BIASED              (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  TailRegrowthPlusOrientationFlip_BIASED              (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z); 
void                  BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED  (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::vector <std::array<int,3>>* old_cut, std::vector <int>* old_ori, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z);
void                  BackFlowFromTailRegrowthPlusOrientationFlip_BIASED  (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::vector <std::array<int,3>>* old_cut, std::vector <int>* old_ori, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void                  PerturbSystem_UNBIASED       (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array <double,8>* E, std::array <double,8>* contacts, std::array <int,9>* attempts, bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, int* move_number, int x, int y, int z);
void                  PerturbSystem_BIASED         (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* contacts, std::array <int,9>* attempts, bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, int* move_number, int x, int y, int z);
void                  PerturbSystem_BIASED_debug   (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, std::array <double,8>* E, std::array <double,8>* contacts, std::array <int,9>* attempts, bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, int* move_number, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

template <typename T>
void reset(T &x){
    x = T();
}

#endif 
