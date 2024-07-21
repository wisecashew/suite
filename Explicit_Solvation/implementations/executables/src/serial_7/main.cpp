#include "Simulation.h"

int main (int argc, char** argv) {

	// INSTANTIATE USER INPUT VARIABLES 
	int opt;           // storage variable to hold options from getopt() [line 36] 
	int lfreq    {-1}; // lattice frequency dump 
	int dfreq    {-1}; // frequency at which coordinates will be dumped out 
	int max_iter {-1}; // number of iteration to perform
	bool verbosity_bool           = false; // boolean for verbosity of output  (default: not verbose)
	bool restart_bool             = false; // boolean for restarts (default: no restarts) 
	bool solvation_align_bool     = false; // boolean for biased solvation shell 
	bool align_lattice_bool       = false; // boolean for aligned lattice 
	bool cosolvent_solvation_bool = false; // boolean for solvation 
	std::string positions          {"__blank__"}; // name of file with initial coords of polymer 
	std::string topology           {"__blank__"}; // name of file with topology of system 
	std::string dfile              {"__blank__"}; // name of coordinate dump file 
	std::string efile              {"__blank__"}; // name of energy dump file 
	std::string mfile              {"__blank__"}; // name of orientation dump file 
	std::string stats_file         {"__blank__"}; // name of file with move statisitics 
	std::string lattice_file_write {"__blank__"}; // name of file where lattice will be dumped to 
	std::string lattice_file_read  {"__blank__"}; // name of file from which lattice will be read 
	std::string SSfile             {"__blank__"}; // name of file to which solvaiton shell will be written

	// define the long options
	static struct option long_options[] = {
		{"frequency-of-sim-dump",     required_argument, 0, 'f'},
		{"frequency-of-lattice-dump", required_argument, 0, 'l'},
		{"total-moves",               required_argument, 0, 'M'},
		{"topology",                  required_argument, 0, 't'},
		{"coords",                    required_argument, 0, 'o'},
		{"initial-coords",            required_argument, 0, 'p'},
		{"energetics",                required_argument, 0, 'u'},
		{"move-statistics",           required_argument, 0, 's'},
		{"orientation",               required_argument, 0, 'e'},
		{"latice-dump",               required_argument, 0, 'L'},
		{"read-lattice",              required_argument, 0, 'R'},
		{"solvation-dump",            required_argument, 0, 'H'},
		{"help",                   no_argument, 0, 'h'},
		{"verbose",                no_argument, 0, 'v'},
		{"restart",                no_argument, 0, 'r'},
		{"align-solvation",        no_argument, 0, 'S'},
		{"align-lattice",          no_argument, 0, 'A'},
		{"solvate-with-cosolvent", no_argument, 0, 'y'},
		{"low-temperature",        no_argument, 0, 'T'},
		{0, 0, 0, 0}  // End of options
	};

	// loop to obtain inputs and assign them to the appropriate variables 
	int option_index = 0;
	while ( (opt = getopt_long(argc, argv, ":l:s:L:R:f:H:M:o:u:p:t:e:vhSAyrT", long_options, &option_index)) != -1 ) {
		switch (opt) {
		
		case 'f':
			dfreq = atoi(optarg); //check
			break;

		case 'l':
			lfreq = atoi(optarg); // check
			break;

		case 'M':
			max_iter = atoi(optarg); // check
			break; 

		case 'h':
			std::cout << 
			"\n" << 
			"Welcome to flogotts (FLOry-huGgins-pOTTs) v1.2.0 for molecular simulations on a cubic lattice (z=26)! \n" << 
			"Last updated: Jun 26, 2024, 11:01 PM. \n" << 
			"Author: satyend@princeton.edu \n" <<
			"\n" << 
			"----------------------------------------------------------------------------------------------------------------------------------\n" << 
			"These are all the inputs the engine accepts for a single run, as of right now: \n\n" <<
			"help                                     [-h, --help]                      (NO ARG)                             Prints out this message. \n"<<
			"verbose flag                             [-v, --verbose]                   (NO ARG)       (NOT REQUIRED)        Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
			"restart flag                             [-r, --restart]                   (NO ARG)       (NOT REQUIRED)        Restarts simulation from final spot of a previous simulation. \n"<<
			"solvation bias flag                      [-y, --solvate-with-cosolvent]    (NO ARG)       (NOT REQUIRED)        Solvated cosolvent right around polymer. \n"<<
			"solvation shell orientation bias flag    [-S, --align-solvation]           (NO ARG)       (NOT REQUIRED)        All particles around polymer have orientation 0. \n"<<
			"lattice orientation bias flag            [-A, --align-lattice]             (NO ARG)       (NOT REQUIRED)        All particles in lattice have orientation 0. \n"<<
			// "Low temperature simulation               [-T, --low-temperature]           (NO ARG)       (NOT REQUIRED)        Simulation will have moves geared for low temperature simulations. \n" <<
			"Dump frequency                           [-f, --frequency-of-sim-dump]     (INTEGER ARG)  (REQUIRED)            Frequency at which coordinates should be dumped out. \n"<<
			"Number of maximum moves                  [-M, --total-moves]               (INTEGER ARG)  (REQUIRED)            Number of MC moves to be run on the system. \n" <<
			"Dump frequency of entire lattice         [-l, --frequency-of-lattice-dump] (INTEGER ARG)  (NOT REQUIRED)        Frequency at which lattice should be dumped out. \n"<<
			"Polymer coordinates                      [-p, --initial-coords]            (STRING ARG)   (REQUIRED)            Name of input file with coordinates of polymer.\n" <<
			"Energy and geometry                      [-t, --topology]                  (STRING ARG)   (REQUIRED)            Name of input file with energetic interactions and geometric bounds.\n" <<
			"Energy of grid                           [-u, --energetics]                (STRING ARG)   (REQUIRED)            Name of output file with energy of system at each step in a file.\n"<<
			"Name of output file                      [-o, --coords]                    (STRING ARG)   (REQUIRED)            Name of output file which will contain coordinates of polymer.\n"<<
			"Lattice file to write to                 [-L, --lattice-dump]              (STRING ARG)   (NOT REQUIRED)        Trajectory file of a previous simulation which can be used to write out current simulation.\n" <<
			"Lattice file to read from                [-R, --read-lattice]              (STRING ARG)   (NOT REQUIRED)        Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
			"Orientation file                         [-e, --orientation]               (STRING ARG)   (NOT REQUIRED)        Name of output file which will contain orientation of ALL particles in system.\n" << 
			"Move statistics file                     [-s, --move-statistics]           (STRING ARG)   (NOT REQUIRED)        Name of output file with move statistics. \n" <<
			"Solvation shell file                     [-H, --solvation-dump]            (STRING ARG)   (NOT REQUIRED)        Name of output file where solvation shell is written out. \n\n";
			exit(EXIT_SUCCESS);
			break;


		case 'p':
			positions = optarg; // check
			break;

		case 'S':
			solvation_align_bool = true; // check
			break;

		case 'A':
			align_lattice_bool = true; // check
			break;

		case 'y': 
			cosolvent_solvation_bool = true; // check
			break;

		case 't':
			topology = optarg; // check
			break;

		case 'o':
			dfile = optarg; // check
			break;

		case 'u':
			efile = optarg; // check
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

		case 'v':
			std::cout << "Output to console will be verbose. " << std::endl;
			verbosity_bool = true; // check
			break;

		case 'e':
			mfile = optarg; // check
			break; 

		case ':':
			std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
			exit(EXIT_FAILURE);
			break; 

		case 'L': 
			// std::cout << "Name of file to write lattice down at end of simulation." << std::endl;
			lattice_file_write=optarg; // check
			break;

		case 'R':
			lattice_file_read=optarg; // check
			std::cout << "Name of file to write lattice down at end of simulation is " << lattice_file_read <<  "." << std::endl;
			break;

		case 'H': 
			SSfile=optarg; // check
			break;

		default:
			std::cout << "A bad option has been provided. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

	}
	
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	auto start = std::chrono::high_resolution_clock::now();

	// Parse inputs... 
	// This command will take all of the above inputs and make sure they are valid. 
	input_parser(dfreq, lfreq, max_iter, restart_bool, positions, topology, dfile, efile, mfile, stats_file, lattice_file_read, lattice_file_write, SSfile); 

	// this is the ultimate simulation object
	// constructing the Simulation object
	Simulation mySim(max_iter, // total moves to perform
	dfreq,                     // total frequency of dumping out statistics 
	lfreq,                     // total frequency of lattice dumps
	restart_bool,              // boolean for restart
	cosolvent_solvation_bool,  // boolean to make sure where the cosolvent isbeing planted
	verbosity_bool,            // boolean to check for verbose output
	align_lattice_bool,        // boolean to check if you want the entire lattice to be aligned
	solvation_align_bool,      // boolean to check if you only want the solvation shell to be aligned
	positions,                 // file holding all positions
	topology,                  // file holding energetic parameters and simulation cell conditions
	dfile,                     // file holding polymer coordinates 
	efile,                     // file holding energy and other simulation conditions
	mfile,                     // file holding all orientations
	stats_file,                // file holding statistics regarding move selection
	lattice_file_write,        // file where the lattice is dumped 
	lattice_file_read,         // file where lattice is read from for restart
	SSfile);                   // file where solvation shell information is dumped
	// end of object instantiation

	// right now, i have only instantiated certain files and certain pathways. 
	mySim.extract_topology_from_file();       // I have the geometry and energies.
	mySim.set_up_system();                    // now that i have all the info, i can set up the simulation lattice
	mySim.set_up_local_dump();                // sets up the dump function
	mySim.set_up_energy_calculator();         // sets up the energy function
	mySim.set_up_run();                       // sets up the run function
	mySim.initialize_pairwise_function_map(); // initialize the pairwise function
	mySim.initialize_neighbor_function_map(); // initialize the neighbor function
	mySim.accelerate_calculate_energy();      // get the energy of the system
	mySim.dump_local();                       // dump out the conditions at step number 0

	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	// OPENING TILES
	mySim.opening_tiles();

	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

	// run the simulation!
	mySim.run();

	// final dumps of lattice and move statistics
	// mySim.dump_lattice();
	mySim.dump_statistics(); 

	auto stop     = std::chrono::high_resolution_clock::now(); 
	auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 

	std::cout << "Time required for computation is " << duration.count() << " microseconds." << std::endl;
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

	return 0;

}
