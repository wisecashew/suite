#ifndef _SIMULATION_H_
#define _SIMULATION_H_
#include "Polymer.h"
#include "Containers.h"
#include "OrientationFlipper.h"

// Note: the hardest part of this simulation, whenever you feel like making some changes, is going 
// to be getting the regrowth right. The biggest issue is the swapping of monomers.
// this is where the linked lists come into play. Be wary while making changes! 
// Usually, making other changes to the code is harmless, but this is always the part where you will be hit by
// segfaults and bad internal structures.

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/
//
//			Definitions for `class Simulation`
//			THIS IS THE GOD CLASS
//
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/

// instantiating so i can use it for defining InteractionFunction and NeighborFunction
class Simulation;

// defining a function class to avoid if-statements and switch-cases when running energy calculations
// ~~~ there is probably a way of doing this without taking as many lines of code as I have. This is a job for future developers. ~~~
// ~~~ the solution likely is somewhere in the land of meta-programming and template definitions. ~~~
typedef void (*InteractionFunction)(Simulation*, Particle*, Particle*, std::array<double,8>*, double*);
void interaction_i_m1_m1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_i_m1_s1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_i_m1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_i_s1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_m1_m1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_m1_s1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_m1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_s1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_m1_m1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_m1_s1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_m1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_s1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);

typedef void (*NeighborFunction)(Simulation*, Particle*, Particle*, std::array<double,8>* contacts, double* energy);
void interaction_i_m1_m1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_i_m1_s1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_i_m1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_i_s1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_m1_m1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_m1_s1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_m1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_p_s1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_m1_m1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_m1_s1(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_m1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);
void interaction_a_s1_s2(Simulation* S, Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy_incr);


class Simulation {
private:
	void (Simulation::*calculate_energy_ptr)();
	void (Simulation::*dump_local_ptr)();
	void (Simulation::*run_ptr)();
	void (Simulation::*system_ptr)();

public:

	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Properties defintions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
	// to first order, let's define certain properties 
	// class properties are organized according to their data type
	// as in, all the `int`s come first, then `doubles`, then `bools`, and so on...
	int x;                // length along x-axis of box. 
	int y;                // length along y-axis of box. 
	int z;                // length along z-axis of box. 
	int max_iter;         // number of moves to perform in the simulation.
	int dfreq;            // frequency of dumps to energydump and coords file. 
	int lfreq;            // frequency of dumping out the entire lattice.
	int Npoly;            // number of polymers in the box.
	int step_number {0};  // current step number of the simulation.
	int move_number {-1}; // number of last selected move.
	// end of `int` properties

	double T         {0}; // temperature of the system
	double frac_c    {0}; // fraction of cosolvent in the system
	double sysEnergy {0}; // energy of the system
	// end of `double` properties

	bool r;         // if true, the simulation is being restarted
	bool s;         // if true, cosolvent is first placed right around polymer. 
	bool v;         // if true, verbose output, if false, non-verbose output. Useful for debugging. 
	bool A;         // if true, the entire lattice has orientation 0.
	bool S;         // if true, the solvation shell including the polymer has orientation 0. 
	bool IMP_BOOL;  // if true, the suggested perturbation has been accepted.  
	// end of `bool` properties

	std::array <int,9>    attempts       = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <int,9>    acceptances    = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> contacts       = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,8> energy_surface = {0, 0, 0, 0, 0, 0, 0, 0};
	// end of `std::array` properties

	std::string positions;          // name of file with initial coords of polymer 
	std::string topology;           // name of file with topology of system 
	std::string dfile;              // name of coordinate dump file 
	std::string efile;              // name of energy dump file 
	std::string mfile;              // name of orientation dump file 
	std::string stats_file;         // name of file with move statisitics 
	std::string lattice_file_write; // name of file where lattice will be dumped to 
	std::string lattice_file_read;  // name of file from which lattice will be read 
	std::string SSfile;             // name of file to which solvation shell will be written
	// end of `std::string` properties

	// define the entire set of particles that we are interested in 
	std::vector <Particle*> Lattice;
	std::vector <Polymer>   Polymers;
	std::vector <Particle*> Solvent;
	std::vector <Particle*> Cosolvent;
	// end of `std::vector` properties, specifically, vectors of particles. 
	// define the custom objects for some of the more complicated moves
	Container rotation_container;                // this is the object important for rotation moves
	EnhancedOrientationFlipper enhanced_flipper; // this is the object important for biased orientation flips


	// define maps for interactions 
	std::map <std::pair<std::string, std::string>, std::tuple<std::string, double, double, int, int>> InteractionMap;
	std::map <std::pair<std::string, std::string>, InteractionFunction> PairwiseFunctionMap; 
	std::map <std::pair<std::string, std::string>, NeighborFunction>    NeighborFunctionMap;

	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Constructor/Destructor definitions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
	// constructor
	Simulation(int max_iter, int dfreq, int lfreq, 
		bool r, bool s, bool v, bool A, bool S, 
		const std::string& positions,
		const std::string& topology, 
		const std::string& dfile, 
		const std::string& efile,
		const std::string& mfile, 
		const std::string& stats_file, 
		const std::string& lattice_file_write,
		const std::string& lattice_file_read,
		const std::string& SSfile)
	: max_iter(max_iter), dfreq(dfreq), lfreq(lfreq), 
	r(r), s(s), v(v), A(A), S(S), 
	positions(positions), 
	topology(topology), 
	dfile(dfile), 
	efile(efile),
	mfile(mfile), 
	stats_file(stats_file), 
	lattice_file_write(lattice_file_write),
	lattice_file_read(lattice_file_read),
	SSfile(SSfile) {};

	// destructor 
	~Simulation(){};

	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Method definitions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

	// opening tiles
	void opening_tiles(){
		std::pair  <std::string, std::string> mm_pair    = std::make_pair("m1", "m1");
		std::pair  <std::string, std::string> ms1_pair   = std::make_pair("m1", "s1");
		std::pair  <std::string, std::string> ms2_pair   = std::make_pair("m1", "s2");
		std::pair  <std::string, std::string> s1s2_pair  = std::make_pair("s1", "s2");
		std::cout << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
		std::cout << "Preparing for take-off." << std::endl << std::endl;
		std::cout << "Chemical information: " << std::endl;
		std::cout << "Number of polymers in system = " << this->Npoly << "." << std::endl << std::endl;
		std::cout << "Geometric information about simulation cell: " << std::endl;
		std::cout << "x = " << this->x <<", y = " << this->y << ", z = "<< this->z << "." << std::endl << std::endl;
		std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
		std::cout << "Temperature = " << this->T << "." << std::endl; 
		std::cout << "Fraction of Solvent II = " << this->frac_c << "." << std::endl;
		std::cout << "Monomer-Monomer energetics: "    << std::get<0>(this->InteractionMap[mm_pair])  << ", Emm_a = "  << this->energy_surface[0] <<", Emm_n = "  << this->energy_surface[1] << ".\nMonomer-Solvent I energetics: "    << std::get<0>(this->InteractionMap[ms1_pair]) << ", Ems1_a = "  << this->energy_surface[2] << ", Ems1_n = "  << this->energy_surface[3] <<"." << std::endl;
		std::cout << "Monomer-Solvent II energetics: " << std::get<0>(this->InteractionMap[ms2_pair]) << ", Ems2_a = " << this->energy_surface[4] <<", Ems2_n = " << this->energy_surface[5] << ".\nSolvent I-Solvent II energetics: " << std::get<0>(this->InteractionMap[s1s2_pair])<< ", Es1s2_a = " << this->energy_surface[6] << ", Es1s2_n = " << this->energy_surface[7] <<"." << std::endl;
		std::cout << "Energy of system is " << this->sysEnergy << "." << std::endl;
		std::cout << "Off to a good start." << std::endl << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
		std::cout << "Initiation complete. We are ready to go. The engine will output information every " << this->dfreq << " configuration(s)." << std::endl; 
		std::cout << "Number of iteration to perform: " << this->max_iter << "." << std::endl;
		std::cout << "Time to fly." << std::endl << std::endl; 
		std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
		return;
	}

	// get the solvation shell
	std::set <int> Simulation::get_solvation_shell();

	// extraction methods (in extraction.cpp)
	std::vector <std::string> extract_content_from_file(std::string filename);
	void extract_number_of_polymers   ();
	void extract_polymers_from_file   ();
	void extract_polymers_from_restart();
	void extract_lattice_from_restart ();
	void extract_topology_from_file   ();

	// lattice set up methods (in setup.cpp)
	void set_up_local_dump();
	void set_up_energy_calculator();
	void set_up_lattice_from_scratch();
	void set_up_from_scratch();
	void set_up_lattice_for_restart ();
	void set_up_polymers_for_restart();
	void set_up_files_for_restart();
	void set_up_for_restart();
	void set_up_system();
	void set_up_run();
	void add_solvent();
	void add_cosolvent();
	void align_lattice();
	void align_solvation_shell();

	// energy calculation methods (in energy.cpp) 
	void selected_pair_interaction(Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy);
	void neighbor_energetics(int lat_idx, std::array<double,8>* contacts_store, double* neighbor_energy);
	double isolated_pair_particle_interaction(Particle* p1, Particle* p2, int* c_idx);
	void accelerate_pair_interaction(Particle* p1, Particle* p2, std::array<double,8>* contacts, double* energy);
	void pair_interaction(Particle* p1, Particle* p2, double* energy);
	void modified_pair_interaction(Particle* p1, Particle* p2, double* energy);
	void accelerate_calculate_energy_solvent();
	void accelerate_calculate_energy_cosolvent();
	void accelerate_calculate_energy(){
		(this->*calculate_energy_ptr)();
	}
	void calculate_energy();
	void initialize_pairwise_function_map(); 
	void initialize_neighbor_function_map();

	// verification methods (in verify.cpp)
	bool check_validity_of_coords(std::array<int,3> v);
	bool check_for_overlaps_within_polymers_raw();
	bool check_for_overlaps_within_polymers();
	bool check_for_solvent_monomer_overlap();
	bool check_for_overlaps_on_lattice();
	bool check_pointers_on_lattice();
	bool check_connectivity_raw();
	bool check_connectivity();
	void check_structures();

	// dump methods (in dump.cpp)
	void dump_energy();
	void dump_polymers();
	void dump_solvation_shell_orientations();
	void dump_solvation_shell();
	void dump_statistics();
	void dump_local_no_ss();
	void dump_local_all(); 
	void dump_local(){
		(this->*dump_local_ptr)();
	};
	void dump_lattice();

	// perturbation methods (in perturb.cpp)
	//////////////////////////////////////////////////////////
	// swap particles
	void perturb_particle_swap            (int lat_idx_1, int lat_idx_2);
	void particle_swap                    (Particle*, int lat_idx_1, int lat_idx_2);
	void particle_swap_with_monomer       (std::array<int,3> monomer_location, int polymer_idx, int monomer_idx, int other_particle_idx);
	void particle_swap_with_update        (Container*, Particle*, int lat_idx_1, int lat_idx_2);
	void particle_swap_with_monomer_update(Container* container, std::array<int,3> monomer_location, int polymer_idx, int monomer_idx, int other_particle_idx);

	//////////////////////////////////////////////////////////
	// rotation of the end
	void perturb_tail_rotation(int p_idx);
	void perturb_head_rotation(int p_idx);

	//////////////////////////////////////////////////////////
	// all about reptation
	void perturb_reptation_forward(int p_idx);
	void perturb_reptation_backward(int p_idx);
	void revert_with_tail_biting(std::array<int,3>* to_slither, int p_idx);
	void revert_without_tail_biting(Particle* tmp, std::array<int,3>* to_slither, std::array<int,3>* loc_0, int p_idx);
	void revert_with_head_butting(std::array<int,3>* to_slither, int p_idx);
	void revert_without_head_butting(Particle*tmp, std::array<int,3>* to_slither, std::array<int,3>* locf, int p_idx);
	void forward_reptation_with_tail_biting(std::array<double,8>* contacts_i, double* E_i, int p_idx);
	void forward_reptation_without_tail_biting(std::array<double,8>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx);
	void backward_reptation_with_head_butting(std::array<double,8>* contacts_i, double* E_i, int p_idx);
	void backward_reptation_without_head_butting(std::array<double,8>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx);

	//////////////////////////////////////////////////////////
	// all about regrowth
	void perturb_regrowth(int p_idx);
	void forward_head_regrowth(std::array <double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx);
	void accept_after_head_regrowth(std::vector<std::array<int,3>>* old_cut, std::vector<std::array<int,3>>* new_cut);
	void backward_head_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth);
	void forward_tail_regrowth(std::array <double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx);
	void accept_after_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::vector <std::array<int,3>>* new_cut);
	void backward_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth);

	//////////////////////////////////////////////////////////
	// perturb orientations
	void perturb_polymer_orientation_flip(int p_idx);
	void perturb_solvation_shell_flip();
	void perturb_lattice_flip();
	void perturb_solvent_exchange();
	void perturb_solvent_exchange_from_shell();
	
	
	//////////////////////////////////////////////////////////

	// final perturbation
	void perturb_system_straight();
	void perturb_system_debug();

	// debugging (found in debug.cpp)
	void debug_calculate_energy(double* db_energy, std::array<double,8>* db_contacts); 
	void debug_pair_interaction(Particle* p1, Particle* p2, std::array<double,8>* db_contacts, double* db_energy);
	void debug_isolated_pair_particle_interaction(Particle* p1, Particle* p2, std::array<double,8>* db_contacts, double* db_energy);
	void debug_checks_energy_contacts(double E_final, std::array<double,8> final_contacts);

	// debugging for perturbations
	void debug_tail_rotation(int p_idx);
	void debug_head_rotation(int p_idx);
	void debug_reptation_forward(int p_idx);
	void debug_reptation_backward(int p_idx);
	void debug_orientation_sampler_forwards   (EnhancedOrientationFlipper* eof, std::array<double,8>* contacts_sys, double E_sys, int iteration_idx, int polymer_idx, int monomer_idx);
	void debug_orientation_sampler_backwards_0(EnhancedOrientationFlipper* eof, std::array<double,8>* contacts_sys, double E_sys, int iteration_idx, int polymer_idx, int monomer_idx);
	void debug_orientation_sampler_backwards  (EnhancedOrientationFlipper* eof, std::array<double,8>* contacts_sys, double E_sys, int iteration_idx, int polymer_idx, int monomer_idx);
	void debug_choose_state_forward(EnhancedOrientationFlipper* eof, int iteration_idx, int polymer_idx, int monomer_idx);
	void debug_polymer_orientation_flip(int p_idx);
	void debug_solvation_shell_flip();
	void debug_lattice_flip();
	void debug_solvent_exchange_from_shell();
	void debug_solvent_exchange();
	void debug_regrowth(int p_idx);
	void debug_forward_head_regrowth(std::array<double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx);
	void debug_accept_after_head_regrowth(std::vector<std::array<int,3>>* old_cut, std::vector<std::array<int,3>>* new_cut);
	void debug_backward_head_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth);
	void debug_forward_tail_regrowth(std::array <double,8>* forw_contacts, double* prob_o_to_n, double* forw_energy, int p_idx, int m_idx);
	void debug_accept_after_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::vector <std::array<int,3>>* new_cut);
	void debug_backward_tail_regrowth(std::vector <std::array<int,3>>* old_cut, std::array <double,8>* back_contacts, double* prob_n_to_o, double* back_energy, int p_idx, int m_idx, int recursion_depth);

	// execution methods (run.cpp)
	void run_straight();
	void run_low_temp();
	void run_debug   ();
	void run(){
		(this->*run_ptr)();
	}

};

std::set <int> Simulation::get_solvation_shell(){

	std::set   <int> solvation_shell_set; 
	std::array <std::array<int,3>,26> ne_list; 
	// int dop = static_cast<int>((*Polymers)[0].chain.size() ); 

	// get the first solvation shell 
	// auto start = std::chrono::high_resolution_clock::now(); 
	for (Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
			ne_list = obtain_ne_list ( p->coords, this->x, this->y, this->z); 
			for (std::array <int,3>& loc: ne_list){
				if (this->Lattice[lattice_index (loc, this->y, this->z)]->ptype[0] == 's'){
					solvation_shell_set.insert (lattice_index (loc, this->y, this->z)); 
				}
			}
		}
	}

	return solvation_shell_set;
}

#endif