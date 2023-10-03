#ifndef _SIMULATION_H_
#define _SIMULATION_H_
#include "ParticleSets.h"
#include "History.h"

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
//
//           DEFINITIONS FOR CLASS Simulation
//
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// Note: the hardest part of this simulation, whenever you feel like making some changes, is going 
// to be getting the regrowth right. The biggest issue is the swapping of monomers.
// this is where the linked lists come into play. Be wary while making changes! 
// Usually, making other changes to the code is harmless, but this is always the part where you will be hit by
// segfaults and bad internal structures.

class Simulation {
public:
	// to first order, let's define certain properties 
	int x;
	int y;
	int z;
	int max_iter;
	int dfreq;
	int lfreq;
	int Npoly;
	double T;
	double frac;
	bool IMP_BOOL;
	bool r;
	bool s;
	bool v;
	bool A; 
	bool S;
	std::array<double,8> E;
	std::string positions;          // name of file with initial coords of polymer 
	std::string topology;           // name of file with topology of system 
	std::string dfile;              // name of coordinate dump file 
	std::string efile;              // name of energy dump file 
	std::string mfile;              // name of orientation dump file 
	std::string stats_file;         // name of file with move statisitics 
	std::string lattice_file_write; // name of file where lattice will be dumped to 
	std::string lattice_file_read ; // name of file from which lattice will be read 

	int    step_number;
	int    move_number;
	double sysEnergy;

	ParticleSets PSETS;

	std::array <int,9>    attempts    = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <int,9>    acceptances = {0, 0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> contacts;

  std::map <std::pair<std::string, std::string>, std::tuple<std::string, double, double, int, int>> InteractionMap; 

	// constructor
	Simulation(int max_iter, int dfreq, int lfreq, 
		 bool r, bool s, bool v, bool A, bool S, const std::string& positions,
		 const std::string& topology, const std::string& dfile, const std::string& efile,
		 const std::string& mfile, const std::string& stats_file, const std::string& lattice_file_write,
		 const std::string& lattice_file_read)
	: max_iter(max_iter), dfreq(dfreq), lfreq(lfreq), 
	r(r), s(s), v(v), A(A), S(S), positions(positions), topology(topology), dfile(dfile), efile(efile),
	mfile(mfile), stats_file(stats_file), lattice_file_write(lattice_file_write),
	lattice_file_read(lattice_file_read) {

		// get number of polymers and other bits of topology
		this->extract_number_of_polymers();
		this->extract_topology_from_file(); 

		// start setting up lattice. this includes probing the positions file and actually getting the polymers
		this->set_up_lattice(); 
		this->calculate_energetics();

	};

	// destructor 
	~Simulation(){};

	// extraction methods
	std::vector <std::string> extract_content_from_file(std::string filename);
	void extract_number_of_polymers(); 
	void extract_polymers_from_file();
	void extract_polymers_from_restart();
	void extract_lattice_from_restart();
	void extract_topology_from_file();

	// lattice set up methods
	void add_solvent();
	void add_cosolvent();
	void align_lattice();
	void align_solvation_shell();
	void set_up_lattice();
	void set_up_lattice_from_restart();

	// energy calculation methods
	void pair_interaction(Particle* p1, Particle* p2, double* energy);
	void modified_pair_interaction(Particle* p1, Particle* p2, double* energy);
	double neighbor_energetics(std::array<double,8>* contacts, int lat_idx);
	double isolated_pair_particle_interaction(Particle* p1, Particle* p2, int* c_idx);
	void calculate_energetics();

	// verification methods
	bool check_validity_of_coords(std::array<int,3> v);
	bool check_for_overlaps_within_polymers_raw(std::vector<Polymer>* Polymers);
	bool check_for_overlaps_within_polymers();
	bool check_for_solvent_monomer_overlap();
	bool check_for_overlaps_on_lattice();
	bool check_pointers_on_lattice();
	bool check_connectivity_raw(std::vector<Polymer>* Polymers);
	bool check_connectivity();
	void check_structures();

	// dump methods
	void dump_energy(int step_num);
	void dump_polymers(int step_num);
	void dump_solvation_shell_orientations(int step_num);
	void dump_statistics(int step_num);
	void dump_lattice (int step_num);

	// perturbation methods
	void swing_monomer                (int m, int deg_poly, History* history_store, double* energy_forw, std::array <double,8>* contacts, double* prob_forw);
	void kick_orientation             (int m, History* history_store, double* energy_forw, std::array <double,8>* contacts, double* prob_forw);
	void kick_orientation_sequence    (int m, History* history_store, double* energy_forw, std::array <double,8>* contacts, double* prob_forw);
	void unswing_monomer              (int m, History* history_store, double* energy_back, std::array <double,8>* contacts, double* prob_back);
	void restore_orientation          (int m, History* history_store, double* energy_back, std::array <double,8>* contacts, double* prob_back);
	void restore_orientation_sequence (int m, History* history_store, double* energy_back, std::array <double,8>* contacts, double* prob_back);
	void restore_old_structure        (History* history_store);
	void adopt_new_structure          (History* history_store);
	void iced_regrowth();
	void perturb_system();
	void perturb_system_V();
	void perturb_system_D();

	// execution methods
	void run(); 
	void run_VERBOSE();
	void run_DEBUG();

};

#endif



