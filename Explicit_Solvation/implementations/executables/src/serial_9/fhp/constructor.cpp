#include "FHP.hpp"

FHP::FHP (int argc, char** argv){
	
	// INSTANTIATE USER INPUT VARIABLES 
	int opt;           // storage variable to hold options from getopt() [line 36] 
	int lfreq    {-1}; // lattice frequency dump 
	int dfreq    {-1}; // frequency at which coordinates will be dumped out 
	int max_iter {-1}; // number of iteration to perform
	bool restart_bool              = false; // boolean for restarts (default: no restarts) 
	bool polymer_bool              = false; // boolean to only perturb polymer configuraiton
	bool dry_bool                  = false; // boolean for a dry simulation
	bool isotropic_bool            = false; // boolean for a spin-less simulation
	bool align_lattice_bool        = false; // boolean for aligned lattice 
	bool solvation_align_bool      = false; // boolean for cosolvent solvation
	bool cosolvent_solvation_bool  = false; // boolean for cosolvent solvation
	bool verbose                   = false; // boolean for debugging
	std::string polymer_coords        {"__blank__"}; // name of polymer coords file
	std::string topology              {"__blank__"}; // name of file with topology of system 
	std::string coord_dump_file       {"__blank__"}; // name of coordinate dump file
	std::string energy_dump_file      {"__blank__"}; // name of energy dump file 
	std::string orientation_dump_file {"__blank__"}; // name of orientation dump file
	std::string solvation_shell_file  {"__blank__"}; // name of solvation shell dump file
	std::string stats_file            {"__blank__"}; // name of file with move statisitics 
	std::string lattice_file_write    {"__blank__"}; // name of file where lattice will be dumped to 
	std::string lattice_file_read     {"__blank__"}; // name of file from which lattice will be read 

	// define the long options
	static struct option long_options[] = {
		{"frequency-of-sim-dump",     required_argument, 0, 'f'},
		{"frequency-of-lattice-dump", required_argument, 0, 'l'},
		{"total-moves",               required_argument, 0, 'M'},
		{"topology",                  required_argument, 0, 't'},
		{"coords",                    required_argument, 0, 'o'},
		{"polymer-coords",            required_argument, 0, 'p'},
		{"energetics",                required_argument, 0, 'u'},
		{"move-statistics",           required_argument, 0, 's'},
		{"orientation",               required_argument, 0, 'e'},
		{"latice-dump",               required_argument, 0, 'L'},
		{"read-lattice",              required_argument, 0, 'R'},
		{"solvation",                 required_argument, 0, 'H'},
		{"restart",                no_argument, 0, 'r'},
		{"align-lattice",          no_argument, 0, 'A'},
		{"align-solvation",        no_argument, 0, 'S'},
		{"solvate-with-cosolvent", no_argument, 0, 'y'},
		{"dry",                    no_argument, 0, 'd'},
		{"isotropic",              no_argument, 0, 'i'},
		{"polymer",                no_argument, 0, 'P'},
		{"verbose",                no_argument, 0, 'v'},
		{"help",                   no_argument, 0, 'h'},
		{0, 0, 0, 0}  // End of options
	};

	// loop to obtain inputs and assign them to the appropriate variables 
	int option_index = 0;
	while ( (opt = getopt_long(argc, argv, ":f:l:M:t:o:p:u:s:e:L:R:H:rASydiPvh", long_options, &option_index)) != -1 ) {
		switch (opt) {
		
		case 'f':
			dfreq = atoi(optarg); // check
			break;

		case 'l':
			lfreq = atoi(optarg); // check
			break;

		case 'M':
			max_iter = atoi(optarg); // check
			break;

		case 't':
			topology = optarg; // check
			break;

		case 'o':
			coord_dump_file = optarg; // check
			break;

		case 'p':
			polymer_coords = optarg; // check
			break;

		case 'u':
			energy_dump_file = optarg; // check
			break;

		case 's':
			stats_file = optarg; // check
			break;

		case 'e': 
			orientation_dump_file = optarg; // check
			break;

		case 'L': 
			lattice_file_write = optarg; // check
			break;

		case 'R':
			lattice_file_read = optarg; // check
			// std::cout << "Name of file to write lattice down at end of simulation is " << lattice_file_read <<  "." << std::endl;
			break;

		case 'H':
			solvation_shell_file = optarg; // check
			break;

		case 'r':
			std::cout << "Simulation will be restarted from the end of previous simulation.\n" ;
			restart_bool = true; // check
			break;

		case 'A':
			align_lattice_bool = true; // check
			break;

		case 'S':
			solvation_align_bool = true; // check
			break;

		case 'y':
			cosolvent_solvation_bool = true; // check
			break;

		case 'd':
			dry_bool = true; // check
			break;

		case 'i':
			isotropic_bool = true; // check
			break;

		case 'P':
			polymer_bool = true; // check
			break;

		case 'v':
			verbose = true; // check
			break;

		case '?':
			std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
			exit(EXIT_FAILURE); 
			break;

		case ':':
			std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
			exit(EXIT_FAILURE);
			break; 

		case 'h':
			std::cout << 
			"\n" << 
			"Welcome to McLattE v2.0 for molecular simulations on a cubic lattice (z=26)! \n" << 
			"You are currently running a FHP simulation. \n" << 
			"Last updated: Nov 03, 2024, 03:23 PM. \n" << 
			"Author: satyend@princeton.edu \n" <<
			"\n" << 
			"----------------------------------------------------------------------------------------------------------------------------------\n" << 
			"This is a Flory-Huggins-Potts simulation.\n" << 
			"----------------------------------------------------------------------------------------------------------------------------------\n" << 
			"These are all the inputs the engine accepts for a single run, as of right now: \n\n" <<
			"Help                                        [-h, --help]                      (NO ARG)                             Prints out this message. \n" <<
			"Verbose                                     [-v, --verbose]                   (NO ARG)                             Runs the debugging simulation. \n" <<
			"Dry simulation                              [-d, --dry]                       (NO ARG)       (NOT REQUIRED)        Run a dry simulation.\n" <<
			"Isotropic simulation                        [-i, --isotropic]                 (NO ARG)       (NOT REQUIRED)        Run a simulation with no orientation-flip moves.\n" <<
			"Polymer configuration moves                 [-P, --polymer]                   (NO ARG)       (NOT REQUIRED)        Run a simulation with only polymer configuraitonal moves.\n"
			"Restart flag                                [-r, --restart]                   (NO ARG)       (NOT REQUIRED)        Restarts simulation from final spot of a previous simulation. \n"<<
			"Lattice orientation bias flag               [-A, --align-lattice]             (NO ARG)       (NOT REQUIRED)        All particles in lattice have orientation 0. \n"<<
			"Align only particles in the solvation shell [-S, --align-solvation]           (NO ARG)       (NOT REQUIRED)        All particles in the solvation shell will have orientation 0. \n" <<
			"Solvate polymer with cosolvent              [-y, --solvate-with-cosolvent]    (NO ARG)       (NOT REQUIRED)        The polymer will be solvated with the cosolvent.\n" <<
			"Frequency of dumping out information        [-f, frequency-of-sim-dump]       (INTEGER ARG)  (REQUIRED)            Frequency at which coordinates should be dumped out. \n" << 
			"Dump frequency of entire lattice            [-l, --frequency-of-lattice-dump] (INTEGER ARG)  (NOT REQUIRED)        Frequency at which lattice should be dumped out. \n"<<
			"Number of maximum moves                     [-M, --total-moves]               (INTEGER ARG)  (REQUIRED)            Number of MC moves to be run on the system. \n" <<
			"Energy and geometry file                    [-t, --topology]                  (STRING ARG)   (REQUIRED)            Name of input file with energetic interactions and geometric bounds.\n" <<
			"Name of output file                         [-o, --coords]                    (STRING ARG)   (REQUIRED)            Name of output file which will contain coordinates of polymer.\n"<<
			"Polymer coordinates                         [-p, --initial-coords]            (STRING ARG)   (REQUIRED)            Name of input file with coordinates of polymer.\n" <<
			"Energy and contacts dump file               [-u, --energetics]                (STRING ARG)   (REQUIRED)            Name of output file with energy of system at each step in a file.\n"<<
			"Move statistics file                        [-s, --move-statistics]           (STRING ARG)   (NOT REQUIRED)        Name of output file with move statistics. \n" << 
			"Orientation file                            [-e, --orientation]               (STRING ARG)   (NOT REQUIRED)        Name of output file which will contain orientation of ALL particles in system.\n" <<
			"Solvation shell information dump file       [-H, --solvation]                 (STRING ARG)   (NOT REQUIRED)        Name of output file which will contain information of the solvation shell.\n" <<
			"Lattice file to write to                    [-L, --lattice-dump]              (STRING ARG)   (NOT REQUIRED)        Trajectory file of a previous simulation which can be used to write out current simulation.\n" <<
			"Lattice file to read from                   [-R, --read-lattice]              (STRING ARG)   (NOT REQUIRED)        Trajectory file of a previous simulation which can be used to start current simulation.\n" << std::endl;
			exit(EXIT_SUCCESS);
			break;

		default:
			std::cout << "A bad option has been provided. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

	}
	
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	// Parse inputs... 
	// This command will take all of the above inputs and make sure they are valid. 
	// std::cout << "\tRunning the input parser..." << std::endl;
	this->setup_parser(dfreq, lfreq, max_iter, restart_bool, topology, energy_dump_file, coord_dump_file, orientation_dump_file, solvation_shell_file, stats_file, lattice_file_write, lattice_file_read);
	this->max_iter                 = max_iter;                 // total moves to perform
	this->dump_freq                = dfreq;                    // total frequency of dumping out statistics 
	this->lattice_dfreq            = lfreq;                    // total frequency of lattice dumps
	this->v                        = verbose;                  // boolean for debugging
	this->r                        = restart_bool;             // boolean for restart
	this->dry                      = dry_bool;                 // boolean for a dry, solvent-free simulation.
	this->polymer                  = polymer_bool;			   // boolean to only allow polymer configuration moves
	this->isotropic                = isotropic_bool;           // boolean for an isotropic simulation (no flip moves)
	this->A                        = align_lattice_bool;       // boolean to check if you want the entire lattice to be aligned
	this->s                        = cosolvent_solvation_bool; // puts all the cosolvent particles around polymer
	this->S                        = solvation_align_bool;     // solvation shell is entirely orientation 0
	this->inp_lattice_file_read    = lattice_file_read;        // file where lattice is read from for restart
	this->inp_topology             = topology;                 // file holding energetic parameters and simulation cell conditions
	this->inp_polymer_coords       = polymer_coords;           // file holding initial coordinates of polymer
	this->out_coord_dump           = coord_dump_file;          // file holding the dumped coordinates of the polymer
	this->out_stats_dump           = stats_file;               // file holding statistics
	this->out_energy_dump          = energy_dump_file;         // file holding energy and other simulation conditions
	this->out_lattice_file_write   = lattice_file_write;       // file where the lattice is dumped
	this->out_solvation_shell_dump = solvation_shell_file;     // file where information about the solvation shell is dumped
	this->out_orientation_dump     = orientation_dump_file;    // file where orientations of particles is dumped
}
