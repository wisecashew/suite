#include <iostream> 
#include <fstream>
#include <vector> 
#include <string> 
#include <array>
#include <map>
#include <utility>
#include <array>
#include <random>
#include <chrono>
#include <getopt.h> 
#include <stdlib.h> 
#include "classes.h"
#include "misc.h"

// obtained all necessary libraries. 
// TO-DO: Check later for redundant calls. 

int main (int argc, char** argv) {

    // INSTANTIATE USER INPUT VARIABLES 
    int opt;          // storage variable to hold options from getopt() [line 36] 
    int dfreq   {-1}; // frequency at which coordinates will be dumped out 
    int max_iter{-1}; // number of iteration to perform
    bool v = false;   // boolean for verbosity of output  (default: not verbose)
    bool r = false;   // boolean for restarts (default: no restarts) 
    bool b = false;   // boolean for biased starting 
    std::string positions          {"__blank__"}; // name of file with initial coords of polymer 
    std::string topology           {"__blank__"}; // name of file with topology of system 
    std::string dfile              {"__blank__"}; // name of coordinate dump file 
    std::string efile              {"__blank__"}; // name of energy dump file 
    std::string mfile              {"__blank__"}; // name of orientation dump file 
    std::string stats_file         {"__blank__"}; // name of file with move statisitics 
    std::string lattice_file_write {"__blank__"}; // name of file where lattice will be dumped to 
    std::string lattice_file_read  {"__blank__"}; // name of file from which lattice will be read 

    // loop to obtain inputs and assign them to the appropriate variables 
    while ( (opt = getopt(argc, argv, ":s:L:R:f:M:o:u:p:t:e:vhbr")) != -1 )
    {
        switch (opt) 
        {
            case 'f':
                dfreq = atoi(optarg); 
                break;

            case 'M':
                max_iter = atoi(optarg); 
                break; 

            case 'h':
                std::cout << 
                "\n" << 
                "Welcome to the Monte Carlo simulation engine (v0.9.1) for polymers and solvents on a cubic lattice (Z=26). \n" << 
		        "Last updated: Sep 8, 2022, 01:09. \n" << 
                "Author: satyend@princeton.edu \n" <<
                "\n" << 
                "----------------------------------------------------------------------------------------------------------------------------------\n" << 
                "These are all the inputs the engine accepts for a single run, as of right now: \n\n" <<
                "help                      [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose flag              [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
                "restart flag              [-r]           (NO ARG REQUIRED)              Restarts simulation from final spot of a previous simulation. \n"<<
                "Dump Frequency            [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
                "Number of maximum moves   [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
                "Polymer coordinates       [-p]           (STRING ARGUMENT REQUIRED)     Name of input file with coordinates of polymer.\n" <<
                "Energy and geometry       [-t]           (STRING ARGUMENT REQUIRED)     Name of input file with energetic interactions and geometric bounds.\n" <<
                "Energy of grid            [-u]           (STRING ARGUMENT REQUIRED)     Name of output file with energy of system at each step in a file.\n"<<
                "Lattice file to write to  [-L]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to write out current simulation.\n" <<
                "Lattice file to read from [-R]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
                "Orientation file          [-e]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain orientation of ALL particles in system.\n" << 
                "Move statistics file      [-s]           (STRING ARGUMENT REQUIRED)     Name of output file with move statistics. \n" <<
                "Name of output file       [-o]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain coordinates of polymer.\n\n";  
                exit(EXIT_SUCCESS);
                break;


            case 'p':
                // std::cout <<"Option p was called with argument " << optarg << std::endl;
                positions = optarg;
                break;    

            case 'b':
                b = true;
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
    // INSTANTIATE SIMULATIONS VARIABLES

    int    step_number    {0};
    // int    move_number    {0};  
    double sysEnergy      {0};
    // bool   IMP_BOOL       {true}; 

    // std::array <int,9>    attempts    = {0,0,0,0,0,0,0,0,0};
    // std::array <int,9>    acceptances = {0,0,0,0,0,0,0,0,0}; 
    std::array <double,8> contacts    = {0,0,0,0,0,0,0,0}; 
    
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // Parse inputs... 
    // This command will take all of the above inputs and make sure they are valid. 
    InputParser ( dfreq, max_iter, r, positions, topology, dfile, efile, mfile, stats_file, lattice_file_read ); 

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // ExtractNumberOfPolymers extracts the total number of chains in the input file. used to reserve vector [line 171] 
    const int N = ExtractNumberOfPolymers(positions); 
    
    // EXTRACT TOPOLOGY FROM FILE 
    std::array <double,13> info_vec {ExtractTopologyFromFile(topology)}; 

    // info_vec is the vector with all the information that constitutes the toplogy of the simulation
    // assign values from info vec to relevant variables 
    const int x             =  info_vec[0] ;
    const int y             =  info_vec[1] ; 
    const int z             =  info_vec[2] ; 
    const double T          =  info_vec[3] ; 
    const double frac       =  info_vec[12]; 
    std::array <double,8> E =  {info_vec[4], info_vec[5], info_vec[6], info_vec[7], info_vec[8], info_vec[9], info_vec[10], info_vec[11]}; 
    
    // initialize custom data structures 
    // this data structure will hold the coordinates of the polymer
    std::vector <Polymer> Polymers; 
    Polymers.reserve(N);

    // this data structure will hold the coordinates of the cosolvent 
    std::vector <Particle*> Cosolvent; 

    // this data structure will hold the coordinates of the solvent 
    std::vector <Particle*> LATTICE;
    LATTICE.reserve (x*y*z); 

    std::vector <int> solvation_shells; 
    solvation_shells.reserve(2*26*26*N); 

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // OPENING TILES

    std::cout << std::endl;
    std::cout << "Preparing for take-off...\n\n" ; 
    std::cout << "Chemical information: " << std::endl;
    std::cout << "Number of polymers is " << N << ".\n\n";
    std::cout << "Geometric information about simulation cell: " << std::endl;
    std::cout << "x = " << x <<", y = " << y << ", z = "<< z << "." << std::endl << std::endl;
    std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
    std::cout << "Temperature = " << T << "." << std::endl; 
    std::cout << "Emm_a = " << E[0] <<", Emm_n = " << E[1] << ", Ems1_a = "<< E[2] << ", Ems1_n = " << E[3] <<".\n";
    std::cout << "Ems2_a = " << E[4] <<", Ems2_n = " << E[5] << ", Es1s2_a = "<< E[6] << ", Es1s2_n = " << E[7] <<".\n";  
    std::cout << "Off to a good start. \n\n";
    std::cout << "--------------------------------------------------------------------\n" << std::endl;
    std::cout << "Running some more checks on input... \n\n" ; 
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    
    // set timers for simulation set-up.  
    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 

    if ( !r ){
        std::cout << "Setting up the lattice from scratch! " << std::endl;
        SetUpLatticeFromScratch (&Polymers, &Cosolvent, &LATTICE, positions, frac, x, y, z);
    
        stop = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 

        std::cout << "Solvation took " << duration.count () << " microseconds." << std::endl;
        std::cout << "Cell has been solvated! \n\n" ;

        if (b) {
            std::cout << "Simulation will have a biased start..." << std::endl;
            BiasTheStart (&Polymers, &LATTICE, x, y, z);
        }

        // CheckStructures (&Polymers, &Cosolvent, &LATTICE, x, y, z);

        dumpPositionsOfPolymers(&Polymers, step_number, dfile); 
    }

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    else {
        std::cout << "Setting up system from a restart file!" << std::endl;
        SetUpLatticeFromRestart (&Polymers, &Cosolvent, &LATTICE, &step_number, lattice_file_read, dfile, positions, x, y, z); 
        
        stop = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 

        CheckStructures (&Polymers, &Cosolvent, &LATTICE, x, y, z);
        
        std::cout << "System set-up took " << duration.count () << " microseconds." << std::endl;
        std::cout << "Simulation cell has been made! \n\n" ;

    }


    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    // THERMODYNAMICS OF SET-UP
    std::cout <<"\nCalculating energy..." << std::endl;
    
    start = std::chrono::high_resolution_clock::now(); 
    sysEnergy   = CalculateEnergy(&Polymers, &Cosolvent, &LATTICE, &E, &contacts, x, y, z); 
    stop = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 
    std::cout << "Time required for serial computation = " << duration.count() << " microseconds. " << std::endl;
    std::cout << "Energy of system is " << sysEnergy << ".\n" << std::endl;
    

    // if i am not restarting, i do not need to dump anything. All the information is already present. 
    if (!r) {
        dumpEnergy      (sysEnergy, step_number, &contacts, efile); 
        dumpOrientation (&Polymers, &LATTICE, step_number, mfile, x, y, z); 
    }
    
    else {
        dumpOrientation (&Polymers, &LATTICE, step_number, mfile, x, y, z);    
    }
    
    std::cout << "Initiation complete. We are ready to go. The engine will output information every " << dfreq << " configuration(s)." << std::endl; 
    std::cout << "Number of iteration to perform: " << max_iter << "." << std::endl;
    // std::cout << "acceptances[1] = " << acceptances[1] << std::endl;
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    //    
    //    ALL INPUTS HAVE BEEN PROCESSED AND DATA STRUCTURES HAVE BEEN SET UP. 
    //    
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    // BEGIN: main loop of simulation! lfg  
    int c_idx = 0, c_idx_n = 0;
    double sysEnergy_ = 0, Ef = 0; 
    double Em = 0, Es = 0, Em_n = 0, Es_n = 0, Epair = 0, Epair_n = 0; 
    std::array <double,8> cs, cm, cs_n, cm_n; 
    std::array <int,3> s_loc = {1,1,0}; 
    std::array <int,3> m_loc = {3,0,0}; 

    Es    = NeighborEnergy (&LATTICE, &E, &cs, lattice_index(s_loc, y, z), x, y, z); 
    Em    = NeighborEnergy (&LATTICE, &E, &cm, lattice_index(m_loc, y, z), x, y, z); 
    Epair = PairEnergy  (LATTICE[lattice_index(s_loc, y, z)], LATTICE[ lattice_index(m_loc, y, z)], &E, &c_idx, x, y, z );

    (Polymers)[0].chain[3]->orientation = 1;

    // prep the swap 
    LATTICE[lattice_index(s_loc, y, z)]->coords = m_loc; 
    Polymers[0].chain[3]->coords                = s_loc; 

    // perform the swap 
    LATTICE[ lattice_index(m_loc, y, z) ] = LATTICE[ lattice_index(s_loc, y, z) ];
    LATTICE[ lattice_index(s_loc, y, z) ] = Polymers[0].chain[3];

    Es_n = NeighborEnergy (&LATTICE, &E, &cs_n, lattice_index(s_loc, y, z), x, y, z);
    Em_n = NeighborEnergy (&LATTICE, &E, &cm_n, lattice_index(m_loc, y, z), x, y, z);     
    Epair_n = PairEnergy  (LATTICE [lattice_index(s_loc, y, z)], LATTICE[ lattice_index(m_loc, y, z)], &E, &c_idx_n, x, y, z); 
    
    Ef = sysEnergy - (Es + Em) + (Es_n + Em_n) - Epair_n + Epair; 
    std::array <double,8> contactsf = add_arrays (subtract_arrays (contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
    contactsf[c_idx]   += 1; 
    contactsf[c_idx_n] -= 1; 


    sysEnergy_ = CalculateEnergy(&Polymers, &Cosolvent, &LATTICE, &E, &contacts, x, y, z); 

    std::cout << "Ef = "         << Ef << std::endl;
    std::cout << "sysEnergy_ = " << sysEnergy_ << std::endl;

    std::cout << "contacts = "; print (contacts); 
    std::cout << "contactsf = "; print (contactsf); 

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Es    = NeighborEnergy (&LATTICE, &E, &cs, lattice_index(s_loc, y, z), x, y, z); 
    Em    = NeighborEnergy (&LATTICE, &E, &cm, lattice_index(m_loc, y, z), x, y, z); 
    Epair = PairEnergy  (LATTICE[lattice_index(s_loc, y, z)], LATTICE[ lattice_index(m_loc, y, z)], &E, &c_idx, x, y, z );

    (Polymers)[0].chain[3]->orientation = 0;
    (LATTICE)[lattice_index(m_loc, y, z)]->coords = s_loc;
    (Polymers)[0].chain[3]->coords  = m_loc;
    
    // perform the swap (since coords were changed, this swap works)
    (LATTICE)[lattice_index(s_loc, y, z)]    = (LATTICE)[lattice_index (m_loc, y, z)];
    (LATTICE)[lattice_index(m_loc, y, z)]    = (Polymers)[0].chain[3];

    Es_n = NeighborEnergy (&LATTICE, &E, &cs_n, lattice_index(s_loc, y, z), x, y, z);
    Em_n = NeighborEnergy (&LATTICE, &E, &cm_n, lattice_index(m_loc, y, z), x, y, z);     
    Epair_n = PairEnergy  (LATTICE [lattice_index(s_loc, y, z)], LATTICE[ lattice_index(m_loc, y, z)], &E, &c_idx_n, x, y, z); 
    
    Ef = sysEnergy_ - (Es + Em) + (Es_n + Em_n) - Epair_n + Epair; 
    contactsf = add_arrays (subtract_arrays (contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
    contactsf[c_idx]   += 1; 
    contactsf[c_idx_n] -= 1; 


    sysEnergy_ = CalculateEnergy(&Polymers, &Cosolvent, &LATTICE, &E, &contacts, x, y, z); 

    std::cout << "Ef = "         << Ef << std::endl;
    std::cout << "sysEnergy_ = " << sysEnergy_ << std::endl;

    std::cout << "contacts = "; print (contacts); 
    std::cout << "contactsf = "; print (contactsf); 

    stop = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 
	
    std::cout << "\n\nTime taken for simulation: " << duration.count()/1e+06 << " seconds.\n"; 
    std::cout << "That is all she wrote. Hope it worked." << std::endl;
    std::cout << "--------------------------------------------------------------------\n\n";

    if (v){Epair++; Epair_n++;};

    return 0;

}
