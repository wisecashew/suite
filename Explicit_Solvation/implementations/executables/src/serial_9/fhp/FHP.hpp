#ifndef _FHP_HPP_
#define _FHP_HPP_
#include "../base/Simulation.hpp"
#include "MonomerSwinger.hpp"
#include "OrientationFlipper.hpp"
#include "Containers.hpp"


class FHP;

typedef void (*InitialInteraction) (FHP*, Particle*, Particle*, std::array<double,CONTACT_SIZE_FHP>&, double&);
void interaction_i_m1_m1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_i_m1_s1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_i_m1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_i_s1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_p_m1_m1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_p_m1_s1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_p_m1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_p_s1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_a_m1_m1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_a_m1_s1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_a_m1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void interaction_a_s1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);

typedef void (*NeighborInteraction)(FHP*, Particle*, Particle*, std::array<double,CONTACT_SIZE_FHP>&, double&);
void neighbor_i_m1_m1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_i_m1_s1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_i_m1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_i_s1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_p_m1_m1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_p_m1_s1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_p_m1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_p_s1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_a_m1_m1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_a_m1_s1   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_a_m1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);
void neighbor_a_s1_s2   (FHP* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy_incr);

class FHP : public Simulation {
private:
	void (FHP::*calculate_energy_ptr)();
	void (FHP::*dump_local_ptr)();
	void (FHP::*dump_stats_ptr)();
	void (FHP::*run_ptr)();
	void (FHP::*system_ptr)();
	

public:

	// additional `int` properties
	int N_poly {0};
	// end of `int` properties

	// additional `double` properties
	double frac_c {0};
	double mag_field {0};
	// end of `double` properties

	// additional `boolean` properties
	bool s;         // if true, cosolvent is first placed right around polymer. 
	bool S;         // if true, the solvation shell including the polymer has orientation 0.
	bool dry;       // if true, run a dry simulation
	bool isotropic; // if true, run a isotropic simulation
	bool polymer;	// if true, run the simple run
	bool field;     // if true, simulation has a magnetic field
	bool IMP_BOOL;  // if true, the suggested perturbation has been accepted
	// end of `boolean` properties

	// additional array-like properties
	std::array <int,9> attempts    = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <int,9> acceptances = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,CONTACT_SIZE_FHP> contacts       = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::array <double,CONTACT_SIZE_FHP> energy_surface = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	// end of array-like properties

	// additional `std::string` propeties
	std::string inp_polymer_coords;
	std::string out_coord_dump;
	std::string out_orientation_dump;
	std::string out_solvation_shell_dump;
	// end of `std::string` properties

	// define the entire set of particles that we are interested in
	std::vector <Polymer>   Polymers;
	std::vector <Particle*> Solvent;
	std::vector <Particle*> Cosolvent;
	// end of `std::vector` properties, specifically, vectors of particles

	// define the custom objects for some of the more complicated moves
	Container                  rotation_container; // this is the object important for rotation moves
	EnhancedOrientationFlipper enhanced_flipper;   // this is the object important for biased orientation flips
	EnhancedMonomerSwinger     enhanced_swing;     // this is the object important for swinging monomers
	// end of definition of custom objects

	// define interaction maps
	// these maps are complicated objects so here are some footholds.
		// For interaction map, 
		// the key is the pair ("particle 1 type", "particle type 2"); 
		// the return value is ("interaction type", "aligned_energy", "misaligned_energy", "aligned contact idx", "misaligned contact idx")
		std::map <std::pair<std::string, std::string>, std::tuple<std::string, double, double, int, int>> InteractionMap;
		// For initial interaction map and neighborhood interaction function map, 
		// the key is the pair ("particle 1 type", "particle type 2"); 
		std::map <std::pair<std::string, std::string>, InitialInteraction>  InitialInteractionFunctionMap; 
		std::map <std::pair<std::string, std::string>, NeighborInteraction> NeighborhoodInteractionFunctionMap;
	
	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Constructor/Destructor definitions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
	// default
	FHP (){};

	// constructor
	FHP (int argc, char** argv);

	// destructor
	~FHP(){};

	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Method definitions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

	// extraction of topology
	std::vector <std::string> extract_content_from_file(std::string filename);
	void extract_number_of_polymers();
	void extract_polymers_from_file();
	void extract_polymers_for_restart();
	void extract_topology_from_file();

	// initialize the system
	void setup_parser(
		int dfreq,    int lfreq,         \
		int max_iter, bool restart_bool, \
		std::string topology,             std::string energy_dump_file,      \
		std::string coords_dump_file,     std::string orientation_dump_file, \
		std::string solvation_shell_file, std::string stats_file,            \
		std::string lattice_file_write,   std::string lattice_file_read);
	void setup_add_solvent();
	void setup_add_cosolvent();
	void setup_align_lattice();
	void setup_align_solvation_shell();
	void setup_from_scratch();
	void setup_lattice_from_restart();
	void setup_polymers_for_restart();
	void setup_lattice();
	void setup_energetics();
	void setup_energy_calculator();
	void setup_dump_stats();
	void setup_dump();
	void setup_run();
	void setup();

	// checker functions
	bool check_validity_of_coords(std::array<int,3> v);
	bool check_for_overlaps_within_polymers_raw();
	bool check_for_overlaps_within_polymers();
	bool check_for_solvent_monomer_overlap(); 
	bool check_for_overlaps_on_lattice(); 
	bool check_connectivity_raw();
	bool check_connectivity();
	bool check_pointers_on_lattice();
	void check_structures();

	// energetics of the simulation
	void initialize_interaction_function_map();
	void initialize_neighbor_function_map();
	void initial_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy);
	void neighbor_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& contacts, double& energy);
	void neighbor_energetics(int, std::array<double,CONTACT_SIZE_FHP>&, double&);
	void calculate_energy_cosolvent();
	void calculate_energy_solvent();
	void calculate_energy_dryfield();
	void energy_compute(){
		(this->*calculate_energy_ptr)();
		return;
	};

	// dumps for the system
	void dump_opening_tiles();
	void dump_energy();
	void dump_polymers();
	void dump_solvation_shell_orientations();
	void dump_solvation_shell();
	void dump_stats_isotropic();
	void dump_stats_all();
	void dump_stats();
	void dump_lattice();
	// void dump_opening_tiles();
	void dump_local_all();
	void dump_local_no_ss();
	void dump_lattice_end();
	void dump(){
		(this->*dump_local_ptr)();
		return;
	};

	// get the solvation shell
	std::set <int> get_solvation_shell();

	// perturbation methods (in perturb.cpp)
	//////////////////////////////////////////////////////////
	// swap particles
	void perturb_particle_swap                    (Particle*, int lat_idx_1, int lat_idx_2);
	void perturb_particle_swap_with_monomer       (std::array<int,3> monomer_location, int polymer_idx, int monomer_idx, int other_particle_idx);
	void perturb_particle_swap_with_update        (Particle*, int lat_idx_1, int lat_idx_2);
	void perturb_particle_swap_with_monomer_update(std::array<int,3> monomer_location, int polymer_idx, int monomer_idx, int other_particle_idx);

	//////////////////////////////////////////////////////////
	// rotation of the end
	void perturb_tail_rotation(int p_idx);
	void perturb_head_rotation(int p_idx);

	//////////////////////////////////////////////////////////
	// all about reptation
	void perturb_reptation_forward              (int p_idx);
	void perturb_reptation_backward             (int p_idx);
	void revert_with_tail_biting                (std::array<int,3>& to_slither, int p_idx);
	void revert_without_tail_biting             (Particle* tmp, std::array<int,3>& to_slither, std::array<int,3>& loc_0, int p_idx);
	void revert_with_head_butting               (std::array<int,3>& to_slither, int p_idx);
	void revert_without_head_butting            (Particle* tmp, std::array<int,3>& to_slither, std::array<int,3>& locf, int p_idx);
	void forward_reptation_with_tail_biting     (std::array<double,CONTACT_SIZE_FHP>& contacts_i, double& E_i, int p_idx);
	void forward_reptation_without_tail_biting  (std::array<double,CONTACT_SIZE_FHP>& contacts_i, std::array<int,3>& to_slither, double& E_i, int p_idx);
	void backward_reptation_with_head_butting   (std::array<double,CONTACT_SIZE_FHP>& contacts_i, double& E_i, int p_idx);
	void backward_reptation_without_head_butting(std::array<double,CONTACT_SIZE_FHP>& contacts_i, std::array<int,3>& to_slither, double& E_i, int p_idx);

	//////////////////////////////////////////////////////////
	// all about regrowth
	void perturb_regrowth(int p_idx);
	void perturb_forward_head_regrowth(int p_idx, int m_idx);
	void perturb_accept_after_head_regrowth(bool not_trap_bool);
	void perturb_backward_head_regrowth(int p_idx, int m_idx, int recursion_depth);
	void perturb_forward_tail_regrowth(int p_idx, int m_idx);
	void perturb_accept_after_tail_regrowth(bool not_trap_bool);
	void perturb_backward_tail_regrowth(int p_idx, int m_idx, int recursion_depth);

	//////////////////////////////////////////////////////////
	void perturb_choose_state_forward(int iteration_idx, int lat_idx);
	void perturb_orientation_sampler_backwards_0(std::array<double,CONTACT_SIZE_FHP>& contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void perturb_orientation_sampler_backwards  (std::array<double,CONTACT_SIZE_FHP>& contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	
	//////////////////////////////////////////////////////////
	// perturb orientations
	void perturb_polymer_orientation_flip(int p_idx);
	void perturb_solvation_shell_flip();
	void perturb_lattice_flip();
	void perturb_solvent_exchange();
	void perturb_solvent_exchange_from_shell();
	
	//////////////////////////////////////////////////////////
	void perturb_orientation_sampler_forwards(std::array<double,CONTACT_SIZE_FHP>& contacts_sys, double E_sys, int iterator_idx, int lat_idx);
	
	//////////////////////////////////////////////////////////
	// final perturbation
	void perturb_system_simple();
	void perturb_system_isotropic();
	void perturb_system_dry();
	void perturb_system_straight();
	void perturb_system_debug();

	// debugging scripts
	void debug_calculate_energy                  (double& db_energy,          std::array<double,CONTACT_SIZE_FHP>& db_contacts);
	void debug_pair_interaction                  (Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_FHP>& db_contacts, double& db_energy);
	void debug_checks_energy_contacts            (double E_final,             std::array<double,CONTACT_SIZE_FHP>& final_contacts);

	// debugging for perturbations
	void debug_tail_rotation      (int p_idx);
	void debug_head_rotation      (int p_idx);
	void debug_reptation_forward  (int p_idx);
	void debug_reptation_backward (int p_idx);
	void debug_orientation_sampler_forwards    (std::array<double,CONTACT_SIZE_FHP>& contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void debug_orientation_sampler_backwards_0 (std::array<double,CONTACT_SIZE_FHP>& contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void debug_orientation_sampler_backwards   (std::array<double,CONTACT_SIZE_FHP>& contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void debug_choose_state_forward            (int iterator_idx, int lat_idx);
	void debug_polymer_orientation_flip        (int p_idx);
	void debug_solvation_shell_flip();
	void debug_lattice_flip();
	void debug_solvent_exchange_from_shell();
	void debug_solvent_exchange();
	void debug_regrowth             (int p_idx);
	void debug_forward_head_regrowth(int p_idx, int m_idx);
	void debug_accept_after_head_regrowth(bool not_trap_bool);
	void debug_backward_head_regrowth    (int p_idx, int m_idx, int recursion_depth);
	void debug_forward_tail_regrowth     (int p_idx, int m_idx);
	void debug_accept_after_tail_regrowth(bool not_trap_bool);
	void debug_backward_tail_regrowth    (int p_idx, int m_idx, int recursion_depth);

	// execute the run 
	void run_isotropic();
	void run_straight();
	void run_simple();
	void run_debug();
	void run_dry();
	void run();

};


#endif
