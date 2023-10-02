#include "Simulation.h"

// obtained all necessary libraries
// start obtaining inputs...

int main (int argc, char** argv) {

	// INSTANTIATE USER INPUT VARIABLES 
	int opt;          // storage variable to hold options from getopt() [line 36] 
    int lfreq   {-1}; // lattice frequency dump 
	int dfreq   {-1}; // frequency at which coordinates will be dumped out 
	int max_iter{-1}; // number of iteration to perform
	bool v = false;   // boolean for verbosity of output  (default: not verbose)
	bool r = false;   // boolean for restarts (default: no restarts) 
	bool S = false;   // boolean for biased solvation shell 
	bool A = false;   // boolean for aligned lattice 
	bool s = false;   // boolean for solvation 
	std::string positions          {"__blank__"}; // name of file with initial coords of polymer 
	std::string topology           {"__blank__"}; // name of file with topology of system 
	std::string dfile              {"__blank__"}; // name of coordinate dump file 
	std::string efile              {"__blank__"}; // name of energy dump file 
	std::string mfile              {"__blank__"}; // name of orientation dump file 
	std::string stats_file         {"__blank__"}; // name of file with move statisitics 
	std::string lattice_file_write {"__blank__"}; // name of file where lattice will be dumped to 
	std::string lattice_file_read  {"__blank__"}; // name of file from which lattice will be read 

	// loop to obtain inputs and assign them to the appropriate variables 
	while ( (opt = getopt(argc, argv, ":l:s:L:R:f:M:o:u:p:t:e:vhSAyr")) != -1 )
	{
	    switch (opt) 
	    {
		case 'f':
			dfreq = atoi(optarg); 
			break;

        case 'l':
            lfreq = atoi(optarg); 
            break;

		case 'M':
			max_iter = atoi(optarg); 
			break; 

		case 'h':
			std::cout << 
			"\n" << 
			"Welcome to my [M]onte [C]arlo [Latt]ice [E]ngine [McLattE] (v1.0.0) for polymers and solvents on a cubic lattice (Z=26). \n" << 
			"Last updated: Aug 4, 2023, 03:01 PM. \n" << 
			"Author: satyend@princeton.edu \n" <<
			"\n" << 
			"----------------------------------------------------------------------------------------------------------------------------------\n" << 
			"These are all the inputs the engine accepts for a single run, as of right now: \n\n" <<
			"help                                     [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
			"verbose flag                             [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
			"restart flag                             [-r]           (NO ARG REQUIRED)              Restarts simulation from final spot of a previous simulation. \n"<<
			"solvation bias flag                      [-y]           (NO ARG REQUIRED)              Solvated cosolvent right around polymer. \n"<<
			"solvation shell orientation bias flag    [-S]           (NO ARG REQUIRED)              All particles around polymer have orientation 0. \n"<<
            "lattice orientation bias flag            [-A]           (NO ARG REQUIRED)              All particles in lattice have orientation 0. \n"<<
			"Dump frequency                           [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
            "Dump frequency of entire lattice         [-l]           (INTEGER ARGUMENT REQUIRED)    Frequency at which lattice should be dumped out. \n"<<                
			"Number of maximum moves                  [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
			"Polymer coordinates                      [-p]           (STRING ARGUMENT REQUIRED)     Name of input file with coordinates of polymer.\n" <<
			"Energy and geometry                      [-t]           (STRING ARGUMENT REQUIRED)     Name of input file with energetic interactions and geometric bounds.\n" <<
			"Energy of grid                           [-u]           (STRING ARGUMENT REQUIRED)     Name of output file with energy of system at each step in a file.\n"<<
			"Lattice file to write to                 [-L]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to write out current simulation.\n" <<
			"Lattice file to read from                [-R]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
			"Orientation file                         [-e]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain orientation of ALL particles in system.\n" << 
			"Move statistics file                     [-s]           (STRING ARGUMENT REQUIRED)     Name of output file with move statistics. \n" <<
			"Name of output file                      [-o]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain coordinates of polymer.\n\n";  
			exit(EXIT_SUCCESS);
			break;


		case 'p':
			positions = optarg;
			break;    

		case 'S':
			S = true;
			break;

		case 'A':
			A = true;
			break;

		case 'y': 
			s = true;
			break;

		case 't':
			topology = optarg; 
			break;

		case 'o':
			dfile = optarg;
			break;

		case 'u':
			efile = optarg;
			break;

		case 's':
			stats_file = optarg;
			break;

		case 'r':
			std::cout << "Simulation will be restarted from the end of previous simulation.\n" ;
			r = true;
			break;

		case '?':
			std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
			exit(EXIT_FAILURE); 
			break;

		case 'v':
			std::cout << "Output to console will be verbose. " << std::endl;
			v = true;
			break;

		case 'e':
			mfile=optarg;
			break; 

		case ':':
			std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
			exit(EXIT_FAILURE);           
			break; 

		case 'L': 
			// std::cout << "Name of file to write lattice down at end of simulation." << std::endl;
			lattice_file_write=optarg; 
			break;

		case 'R':
			// std::cout << "Name of file to read lattice to restart simulation." << std::endl;
			lattice_file_read=optarg; 
			std::cout << "Name of file to write lattice down at end of simulation is " << lattice_file_read <<  "." << std::endl;
			break;                

		}
	}

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // Parse inputs... 
    // This command will take all of the above inputs and make sure they are valid. 
    auto start = std::chrono::high_resolution_clock::now(); 
    input_parser(dfreq, lfreq, max_iter, r, positions, topology, dfile, efile, mfile, stats_file, lattice_file_read); 

    Simulation mySim(max_iter, dfreq, lfreq, r, s, v, A, S, positions, topology,
        dfile, efile, mfile, stats_file, lattice_file_write, lattice_file_read);

    std::pair  <std::string, std::string> mm_pair    = std::make_pair ("m1", "m1");
    std::pair  <std::string, std::string> ms1_pair   = std::make_pair ("m1", "s1");
    std::pair  <std::string, std::string> ms2_pair   = std::make_pair ("m1", "s2");
    std::pair  <std::string, std::string> s1s2_pair  = std::make_pair ("s1", "s2");

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // OPENING TILES

    std::cout << std::endl;
    std::cout << "Preparing for take-off...\n\n" ; 
    std::cout << "Chemical information: " << std::endl;
    std::cout << "Number of polymers in system = " << mySim.Npoly << ".\n\n";
    std::cout << "Geometric information about simulation cell: " << std::endl;
    std::cout << "x = " << mySim.x <<", y = " << mySim.y << ", z = "<< mySim.z << "." << std::endl << std::endl;
    std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
    std::cout << "Temperature = " << mySim.T << "." << std::endl; 
    std::cout << "Fraction of Solvent II = " << mySim.frac << "." << std::endl;
    std::cout << "Monomer-Monomer energetics: "    << std::get<0>(mySim.InteractionMap[mm_pair])  << ", Emm_a = "  << mySim.E[0] <<", Emm_n = "  << mySim.E[1] << ".\nMonomer-Solvent I energetics: "    << std::get<0>(mySim.InteractionMap[ms1_pair]) << ", Ems1_a = "  << mySim.E[2] << ", Ems1_n = "  << mySim.E[3] <<".\n";
    std::cout << "Monomer-Solvent II energetics: " << std::get<0>(mySim.InteractionMap[ms2_pair]) << ", Ems2_a = " << mySim.E[4] <<", Ems2_n = " << mySim.E[5] << ".\nSolvent I-Solvent II energetics: " << std::get<0>(mySim.InteractionMap[s1s2_pair])<< ", Es1s2_a = " << mySim.E[6] << ", Es1s2_n = " << mySim.E[7] <<".\n";  
    std::cout << "Energy of system is " << mySim.sysEnergy << ".\n" << std::endl;
    std::cout << "Off to a good start. \n\n";
    std::cout << "--------------------------------------------------------------------\n" << std::endl;
    
    auto stop     = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 

    std::cout << "Time required for set up is " << duration.count() << " microseconds." << std::endl;

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    if (!mySim.r){
        mySim.dump_energy(0);
        mySim.dump_solvation_shell_orientations(0);
        mySim.dump_polymers(0);
        mySim.dump_lattice(0);
    }

    std::cout << "Initiation complete. We are ready to go. The engine will output information every " << mySim.dfreq << " configuration(s)." << std::endl; 
    std::cout << "Number of iteration to perform: " << mySim.max_iter << "." << std::endl;

    return 0;

}
