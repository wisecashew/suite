#ifndef _POTTS_HPP_
#define _POTTS_HPP_
#include "../base/Simulation.hpp"

constexpr int CONTACT_SIZE_POTTS = 2;

class Potts : public Simulation {
public:

	// some types for Potts
	typedef void (*InitialInteraction) (Potts*, Particle*, Particle*, std::array<double,CONTACT_SIZE_POTTS>&, double&);
	typedef void (*NeighborInteraction)(Potts*, Particle*, Particle*, std::array<double,CONTACT_SIZE_POTTS>&, double&);
	
	// number of attempts
	int n_attempts    = 0;
	int n_acceptances = 0;

	// thermodynamic contacts and energy of the system
	std::array <double,CONTACT_SIZE_POTTS> contacts       = {0.0, 0.0};
	std::array <double,CONTACT_SIZE_POTTS> energy_surface = {0.0, 0.0};

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
	// constructor
	Potts(int argc, char** argv){
		// INSTANTIATE USER INPUT VARIABLES 
		int opt;           // storage variable to hold options from getopt() [line 36] 
		int lfreq    {-1}; // lattice frequency dump 
		int dfreq    {-1}; // frequency at which coordinates will be dumped out 
		int max_iter {-1}; // number of iteration to perform
		bool restart_bool              = false; // boolean for restarts (default: no restarts) 
		bool align_lattice_bool        = false; // boolean for aligned lattice 
		std::string topology           {"__blank__"}; // name of file with topology of system 
		std::string energy_dump_file   {"__blank__"}; // name of energy dump file 
		std::string stats_file         {"__blank__"}; // name of file with move statisitics 
		std::string lattice_file_write {"__blank__"}; // name of file where lattice will be dumped to 
		std::string lattice_file_read  {"__blank__"}; // name of file from which lattice will be read 

		// define the long options
		static struct option long_options[] = {
			{"frequency-of-lattice-dump", required_argument, 0, 'l'},
			{"total-moves",               required_argument, 0, 'M'},
			{"topology",                  required_argument, 0, 't'},
			{"energetics",                required_argument, 0, 'u'},
			{"move-statistics",           required_argument, 0, 's'},
			{"latice-dump",               required_argument, 0, 'L'},
			{"read-lattice",              required_argument, 0, 'R'},
			{"help",                   no_argument, 0, 'h'},
			{"verbose",                no_argument, 0, 'v'},
			{"restart",                no_argument, 0, 'r'},
			{"align-lattice",          no_argument, 0, 'A'},
			{0, 0, 0, 0}  // End of options
		};

		// loop to obtain inputs and assign them to the appropriate variables 
		int option_index = 0;
		while ( (opt = getopt_long(argc, argv, ":l:M:t:u:s:L:R:hAr", long_options, &option_index)) != -1 ) {
			switch (opt) {
			
			case 'l':
				lfreq = atoi(optarg); // check
				break;

			case 'M':
				max_iter = atoi(optarg); // check
				break; 

			case 'h':
				std::cout << 
				"\n" << 
				"Welcome to McLattE v2.0 for molecular simulations on a cubic lattice (z=26)! \n" << 
				"You are currently running a Potts simulation. \n" << 
				"Last updated: Nov 03, 2024, 03:23 PM. \n" << 
				"Author: satyend@princeton.edu \n" <<
				"\n" << 
				"----------------------------------------------------------------------------------------------------------------------------------\n" << 
				"These are all the inputs the engine accepts for a single run, as of right now: \n\n" <<
				"help                                     [-h, --help]                      (NO ARG)                             Prints out this message. \n"<<
				"restart flag                             [-r, --restart]                   (NO ARG)       (NOT REQUIRED)        Restarts simulation from final spot of a previous simulation. \n"<<
				"lattice orientation bias flag            [-A, --align-lattice]             (NO ARG)       (NOT REQUIRED)        All particles in lattice have orientation 0. \n"<<
				"Number of maximum moves                  [-M, --total-moves]               (INTEGER ARG)  (REQUIRED)            Number of MC moves to be run on the system. \n" <<
				"Dump frequency of entire lattice         [-l, --frequency-of-lattice-dump] (INTEGER ARG)  (NOT REQUIRED)        Frequency at which lattice should be dumped out. \n"<<
				"Energy and geometry file                 [-t, --topology]                  (STRING ARG)   (REQUIRED)            Name of input file with energetic interactions and geometric bounds.\n" <<
				"Energy and contacts dump file            [-u, --energetics]                (STRING ARG)   (REQUIRED)            Name of output file with energy of system at each step in a file.\n"<<
				"Lattice file to write to                 [-L, --lattice-dump]              (STRING ARG)   (NOT REQUIRED)        Trajectory file of a previous simulation which can be used to write out current simulation.\n" <<
				"Lattice file to read from                [-R, --read-lattice]              (STRING ARG)   (NOT REQUIRED)        Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
				"Move statistics file                     [-s, --move-statistics]           (STRING ARG)   (NOT REQUIRED)        Name of output file with move statistics. \n" << std::endl;
				exit(EXIT_SUCCESS);
				break;

			case 'A':
				align_lattice_bool = true; // check
				break;

			case 't':
				topology = optarg; // check
				break;

			case 'u':
				energy_dump_file = optarg; // check
				break;

			case 's':
				stats_file = optarg; // check
				break;

			case 'r':
				std::cout << "Simulation will be restarted from the end of previous simulation.\n" ;
				restart_bool = true; // check
				break;

			case '?':
				std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
				exit(EXIT_FAILURE); 
				break;

			case ':':
				std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
				exit(EXIT_FAILURE);
				break; 

			case 'L': 
				// std::cout << "Name of file to write lattice down at end of simulation." << std::endl;
				lattice_file_write = optarg; // check
				break;

			case 'R':
				lattice_file_read = optarg; // check
				std::cout << "Name of file to write lattice down at end of simulation is " << lattice_file_read <<  "." << std::endl;
				break;

			default:
				std::cout << "A bad option has been provided. Exiting..." << std::endl;
				exit(EXIT_FAILURE);
			}

		}
		
		//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
		// Parse inputs... 
		// This command will take all of the above inputs and make sure they are valid. 
		std::cout << "\tRunning the input parser..." << std::endl;
		this->setup_parser(lfreq, max_iter, restart_bool, inp_topology, out_energy_dump, out_stats_dump, out_lattice_file_write, lattice_file_read);
		std::cout << "\tRan the parser. Setting up the properties of the Potts object." << std::endl;
		this->max_iter                 = max_iter;              // total moves to perform
		this->dump_freq                = dfreq;                 // total frequency of dumping out statistics 
		this->lattice_dfreq            = lfreq;                 // total frequency of lattice dumps
		this->r                        = restart_bool;          // boolean for restart
		this->A                        = align_lattice_bool;    // boolean to check if you want the entire lattice to be aligned
		this->inp_lattice_file_read    = lattice_file_read;     // file where lattice is read from for restart
		this->inp_topology             = topology;              // file holding energetic parameters and simulation cell conditions
		this->out_stats_dump           = stats_file;            // file holding statistics
		this->out_energy_dump          = energy_dump_file;      // file holding energy and other simulation conditions
		this->out_lattice_file_write   = lattice_file_write;    // file where the lattice is dumped 
		std::cout << "\tSet up the properties." << std::endl;
	}

	~Potts(){};

	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Method definitions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

	// extraction of topology
	void extract_topology();

	// initialize the system
	void setup_parser(int lfreq, int max_iter, \
		bool restart_bool, \
		std::string topology,   std::string energy_dump_file, \
		std::string stats_file, std::string lattice_file_write, \
		std::string lattice_file_read);
	void setup_lattice_from_scratch();
	void setup_lattice_from_restart();
	void setup_lattice();
	void setup_energetics();
	void setup();

	// energetics of the simulation
	void initial_energetics_map(); 
	void neighbor_energetic_map();
	void initial_energetics(Particle*, Particle*, std::array<double,CONTACT_SIZE_POTTS>&, double&);
	void selected_pair_interaction(Particle*, Particle*, std::array<double,CONTACT_SIZE_POTTS>&, double&);
	void neighbor_energetics(int, std::array<double,CONTACT_SIZE_POTTS>&, double&);
	void energy_compute();

	// general print methods
	void print_opening_tiles();

	// dumps for the system
	void dump_energy();
	void dump_stats();
	void dump_lattice();
	void dump();

	// perturbation
	void perturb_lattice_flip();

	// launch the fundamental loop
	void simulate();

	// execute the run 
	void run();

	// debugging scripts
	void debug_pair_interaction      (Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE_POTTS>& db_contacts, double& db_energy);
	void debug_calculate_energy      (double& db_energy, std::array<double,CONTACT_SIZE_POTTS>& db_contacts);
	void debug_checks_energy_contacts(double& E_final,   std::array<double,CONTACT_SIZE_POTTS>& final_contacts);
	void debug_lattice_flip();
	void debug_simulation();
	void debug();

};


#endif
