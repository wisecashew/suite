#ifndef _SIMULATION_H_
#define _SIMULATION_H_
#include "Polymer.hpp"
#include "Containers.hpp"
#include "OrientationFlipper.hpp"
#include "MonomerSwinger.hpp"

// Note: the hardest part of this simulation engine is going 
// to be getting the regrowth right. The biggest issue is the swapping of monomers.
// this is where the linked lists come into play. Be wary while making changes!
// Usually, making other changes to the code is harmless; but this is always the part where you will be hit by
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
typedef void (*InteractionFunction)(Simulation*, Particle*, Particle*, std::array<double,CONTACT_SIZE>*, double*);
void interaction_m1_s0     (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_i_sp_sp   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_i_m1_m1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_i_m1_s1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_i_m1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_i_s1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_p_sp_sp   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_p_m1_m1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_p_m1_s1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_p_m1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_p_s1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_a_sp_sp   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_a_m1_m1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_a_m1_s1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_a_m1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_a_s1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void interaction_symm_sp_sp(Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);

typedef void (*NeighborFunction)(Simulation*, Particle*, Particle*, std::array<double,CONTACT_SIZE>* contacts, double* energy);
void neighbor_m1_s0     (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_i_sp_sp   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_i_m1_m1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_i_m1_s1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_i_m1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_i_s1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_p_sp_sp   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_p_m1_m1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_p_m1_s1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_p_m1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_p_s1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_a_sp_sp   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_a_m1_m1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_a_m1_s1   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_a_m1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_a_s1_s2   (Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);
void neighbor_symm_sp_sp(Simulation* S, Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy_incr);


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
	// as in, all the `int`s come first, then `double`s, then `bool`s, and so on...
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
	bool isotropic; // if true, the simulation will be isotropic
	bool potts;     // if true, the simulation will be a standard potts model simulation
	bool dry;       // if true, the simulation will be solvent free
	bool IMP_BOOL;  // if true, the suggested perturbation has been accepted.  
	// end of `bool` properties

	// define the polymer magnetization
	std::array <double,3> polymer_magnetization = {0, 0, 0};

	// some other properties
	std::array <int,9>    attempts        = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <int,9>    acceptances     = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,CONTACT_SIZE> contacts       = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	std::array <double,CONTACT_SIZE> energy_surface = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
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
	Container                  rotation_container; // this is the object important for rotation moves
	EnhancedOrientationFlipper enhanced_flipper;   // this is the object important for biased orientation flips
	EnhancedMonomerSwinger     enhanced_swing;     // this is the object important for swinging monomers

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
	Simulation(int argc, char** argv){
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
		bool isotropic_bool           = false; // boolean for isotropic simulation
		bool potts_bool               = false; // boolean for a potts simulation
		bool dry_bool                 = false; // boolean for a dry, coarse-grained simulation
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
			{"isotropic",              no_argument, 0, 'I'},
			{"potts",                  no_argument, 0, 'P'},
			{"dry",                    no_argument, 0, 'D'},
			{0, 0, 0, 0}  // End of options
		};

		// loop to obtain inputs and assign them to the appropriate variables 
		int option_index = 0;
		while ( (opt = getopt_long(argc, argv, ":l:s:L:R:f:H:M:o:u:p:t:e:vhSAyrTIPD", long_options, &option_index)) != -1 ) {
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
				"isotropic polymer simulation flag        [-I, --isotropic]                 (NO ARG)       (NOT REQUIRED)        Simulation will have no orientation and solvent moves. \n"<<
				"Run a Potts simulation                   [-P, --potts]                     (NO ARG)       (NOT REQUIRED)        A Potts model simulation will be run. \n" <<
				"Run a dry simulation                     [-D, --dry]                     (NO ARG)         (NOT REQUIRED)        A dry simulation will be run. \n" <<
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

			case 'I':
				isotropic_bool = true; // check
				break;

			case 'D':
				dry_bool = true; // check
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

			case 'P':
				potts_bool = true;
				break;

			default:
				std::cout << "A bad option has been provided. Exiting..." << std::endl;
				exit(EXIT_FAILURE);
			}

		}
		
		//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
		// Parse inputs... 
		// This command will take all of the above inputs and make sure they are valid. 
		input_parser(dfreq, lfreq, max_iter, restart_bool, potts_bool, dry_bool, positions, topology, dfile, efile, mfile, stats_file, lattice_file_read, lattice_file_write, SSfile); 
		this->max_iter                 = max_iter;                  // total moves to perform
		this->dfreq                    = dfreq;                     // total frequency of dumping out statistics 
		this->lfreq                    = lfreq;                     // total frequency of lattice dumps
		this->r                        = restart_bool;              // boolean for restart
		this->s                        = cosolvent_solvation_bool;  // boolean to make sure where the cosolvent isbeing planted
		this->v                        = verbosity_bool;            // boolean to check for verbose output
		this->A                        = align_lattice_bool;        // boolean to check if you want the entire lattice to be aligned
		this->S                        = solvation_align_bool;      // boolean to check if you only want the solvation shell to be aligned
		this->isotropic                = isotropic_bool;            // boolean to check if simulation is isotropic
		this->potts                    = potts_bool;                // boolean to run a potts simulation
		this->dry                      = dry_bool;                  // boolean to run a dry simulation
		this->positions                = positions;                 // file holding all positions
		this->topology                 = topology;                  // file holding energetic parameters and simulation cell conditions
		this->dfile                    = dfile;                     // file holding polymer coordinates 
		this->efile                    = efile;                     // file holding energy and other simulation conditions
		this->mfile                    = mfile;                     // file holding all orientations
		this->stats_file               = stats_file;                // file holding statistics regarding move selection
		this->lattice_file_write       = lattice_file_write;        // file where the lattice is dumped 
		this->lattice_file_read        = lattice_file_read;         // file where lattice is read from for restart
		this->SSfile                   = SSfile;                    // file where solvation shell information is dumped

	}

	// destructor 
	~Simulation(){};

	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Method definitions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

	// opening tiles
	void opening_tiles(){
		
		std::cout << std::endl;

		if (this->potts){
			std::pair  <std::string, std::string> spsp_pair  = std::make_pair("sp", "sp");
			std::cout << "--------------------------------------------------------------------" << std::endl;
			std::cout << "This is a Potts model simulation." << std::endl;
			std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
			std::cout << "Preparing for take-off." << std::endl << std::endl;
			std::cout << "Geometric information about simulation cell: " << std::endl;
			std::cout << "x = " << this->x <<", y = " << this->y << ", z = "<< this->z << "." << std::endl << std::endl;
			std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
			std::cout << "Temperature = " << this->T << "." << std::endl;
			std::cout << "Particle-particle energetics: "    << std::get<0>(this->InteractionMap[spsp_pair])  << ", Espsp_a = "  << this->energy_surface[8] <<", Espsp_n = "  << this->energy_surface[9] << ".\n";
		}
		
		else {
			std::pair  <std::string, std::string> mm_pair    = std::make_pair("m1", "m1");
			std::pair  <std::string, std::string> ms1_pair   = std::make_pair("m1", "s1");
			std::pair  <std::string, std::string> ms2_pair   = std::make_pair("m1", "s2");
			std::pair  <std::string, std::string> s1s2_pair  = std::make_pair("s1", "s2");

			std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
			if (this->isotropic){
				std::cout << "Simulation will have no spin-based moves." << std::endl;
			}
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
		}

		std::cout << "Energy of system is " << this->sysEnergy << "." << std::endl;
		std::cout << "Off to a good start." << std::endl << std::endl;
		std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
		std::cout << "Initiation complete. We are ready to go. The engine will output information every " << this->dfreq << " configuration(s)." << std::endl; 
		std::cout << "Number of iteration to perform: " << this->max_iter << "." << std::endl;
		std::cout << "Time to fly." << std::endl << std::endl; 
		std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
		
		return;
	}

	// extraction methods (in extraction.cpp)
	std::vector <std::string> extract_content_from_file(std::string filename);
	void extract_number_of_polymers   ();
	void extract_polymers_from_file   ();
	void extract_polymers_from_restart();
	void extract_lattice_from_restart ();
	void extract_topology_from_file   ();
	void extract_topology_for_potts   ();

	// lattice set up methods (in setup.cpp)
	void set_up_local_dump();
	void set_up_energy_calculator();
	void set_up_lattice_from_scratch();
	void set_up_from_scratch();
	void set_up_from_scratch_potts();
	void set_up_lattice_for_restart ();
	void set_up_polymers_for_restart();
	void set_up_files_for_restart();
	void set_up_for_restart();
	void set_up_system();
	void set_up_style_of_run();
	void set_up_FHP();
	void set_up_Potts();
	void set_up_add_solvent();
	void set_up_add_cosolvent();
	void set_up_align_lattice();
	void set_up_align_solvation_shell();

	// get statistics
	std::set <int> Simulation::get_solvation_shell();

	// energy calculation methods (in energy.cpp) 
	void selected_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy);
	void neighbor_energetics(int lat_idx, std::array<double,CONTACT_SIZE>* contacts_store, double* neighbor_energy);
	double isolated_pair_particle_interaction(Particle* p1, Particle* p2, int* c_idx);
	void accelerate_pair_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* contacts, double* energy);
	void pair_interaction(Particle* p1, Particle* p2, double* energy);
	void modified_pair_interaction(Particle* p1, Particle* p2, double* energy);
	void accelerate_calculate_energy_solvent();
	void accelerate_calculate_energy_cosolvent();
	void accelerate_calculate_energy_potts();
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
	void dump_potts();
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
	void dump_lattice_end();

	// perturbation methods (in perturb.cpp)
	//////////////////////////////////////////////////////////
	// swap particles
	void perturb_particle_swap                    (int lat_idx_1, int lat_idx_2);
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
	void revert_with_tail_biting                (std::array<int,3>* to_slither, int p_idx);
	void revert_without_tail_biting             (Particle* tmp, std::array<int,3>* to_slither, std::array<int,3>* loc_0, int p_idx);
	void revert_with_head_butting               (std::array<int,3>* to_slither, int p_idx);
	void revert_without_head_butting            (Particle*tmp, std::array<int,3>* to_slither, std::array<int,3>* locf, int p_idx);
	void forward_reptation_with_tail_biting     (std::array<double,CONTACT_SIZE>* contacts_i, double* E_i, int p_idx);
	void forward_reptation_without_tail_biting  (std::array<double,CONTACT_SIZE>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx);
	void backward_reptation_with_head_butting   (std::array<double,CONTACT_SIZE>* contacts_i, double* E_i, int p_idx);
	void backward_reptation_without_head_butting(std::array<double,CONTACT_SIZE>* contacts_i, std::array<int,3>* to_slither, double* E_i, int p_idx);

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
	void perturb_orientation_sampler_backwards_0(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void perturb_orientation_sampler_backwards  (std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx);

	//////////////////////////////////////////////////////////
	// perturb orientations
	void perturb_polymer_orientation_flip(int p_idx);
	void perturb_solvation_shell_flip();
	void perturb_lattice_flip();
	void perturb_solvent_exchange();
	void perturb_solvent_exchange_from_shell();
	
	//////////////////////////////////////////////////////////
	void perturb_orientation_sampler_forwards(std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iterator_idx, int lat_idx);
	
	//////////////////////////////////////////////////////////
	// final perturbation
	void perturb_potts();
	void perturb_system_isotropic();
	void perturb_system_straight();
	void perturb_system_debug();
	void perturb_system_dry();
	void perturb_system_dry_debug();

	// debugging (found in debug.cpp)
	void debug_calculate_energy                  (double* db_energy,          std::array<double,CONTACT_SIZE>* db_contacts); 
	void debug_pair_interaction                  (Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* db_contacts, double* db_energy);
	void debug_isolated_pair_particle_interaction(Particle* p1, Particle* p2, std::array<double,CONTACT_SIZE>* db_contacts, double* db_energy);
	void debug_checks_energy_contacts            (double E_final,             std::array<double,CONTACT_SIZE>  final_contacts);

	// debugging for perturbations
	void debug_tail_rotation      (int p_idx);
	void debug_head_rotation      (int p_idx);
	void debug_reptation_forward  (int p_idx);
	void debug_reptation_backward (int p_idx);
	void debug_orientation_sampler_forwards    (std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void debug_orientation_sampler_backwards_0 (std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void debug_orientation_sampler_backwards   (std::array<double,CONTACT_SIZE>* contacts_sys, double E_sys, int iteration_idx, int lat_idx);
	void debug_choose_state_forward            (int iterator_idx, int lat_idx);
	void debug_polymer_orientation_flip        (int p_idx);
	void debug_polymer_orientation_flip_dry    (int p_idx);
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

	// execution methods (run.cpp)
	void run_straight ();
	void run_isotropic();
	void run_debug    ();
	void run_potts    ();
	void run(){
		(this->*run_ptr)();
	}

};

#endif
