#include <iostream> 
#include <fstream>
#include <vector> 
#include <string> 
#include <array>
#include <map>
#include <utility>
#include <array>
#include <random>
#include <numeric>
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
    char end ='x';
    bool v = false;   // boolean for verbosity of output  (default: not verbose)
    bool r = false;   // boolean for restarts (default: no restarts) 
    std::string positions          {"__blank__"}; // name of file with initial coords of polymer 
    std::string target             {"__blank__"}; // name of file with target coords of polymer
    std::string topology           {"__blank__"}; // name of file with topology of system 
    std::string dfile              {"__blank__"}; // name of coordinate dump file 
    std::string efile              {"__blank__"}; // name of energy dump file 
    std::string mfile              {"__blank__"}; // name of orientation dump file 
    std::string stats_file         {"__blank__"}; // name of file with move statisitics 
    std::string lattice_file_write {"__blank__"}; // name of file where lattice will be dumped to 
    std::string lattice_file_read  {"__blank__"}; // name of file from which lattice will be read 

    while ( (opt = getopt(argc, argv, ":s:L:R:f:M:o:u:p:t:e:vhr")) != -1 )
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
                "Welcome to the Monte Carlo unit tester (v0.0) for polymers and solvent on a cubic lattice (Z=26). \n" << 
		        "Last updated: Aug 17, 2022, 14:56. \n" << 
                "Author: satyend@princeton.edu \n" <<
                "\n" << 
                "----------------------------------------------------------------------------------------------------------------------------------\n" << 
                "These are all the inputs the engine accepts for a single run, as of right now: \n\n" <<
                "help                      [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose flag              [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
                "restart flag              [-r]           (NO ARG REQUIRED)              Restarts simulation from final spot of a previous simulation. \n"<<
                "Dump Frequency            [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<
                "End                       [-y]           (CHARACTER ARGUMENT REQUIRED)  Pick which end the rotation is going to be tested (h or t). \n" << 
                "Number of maximum moves   [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
                "Polymer coordinates       [-p]           (STRING ARGUMENT REQUIRED)     Name of input file with coordinates of polymer.\n" <<
                "Target coordinates        [-q]           (STRING ARGUMENT REQUIRED)     Name of input file of coordinates of target polymer.\n" <<
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
            
            case 'q':
                
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
            
            case 'y':
                end = optarg;
                break;
            
            case 'q':
                target = optarg; 
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

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // INSTANTIATE SIMULATIONS VARIABLES

    int    step_number    {0};
    int    move_number    {0};  
    double sysEnergy      {0};
    bool   IMP_BOOL       {true}; 

    std::array <int,9>    attempts    = {0,0,0,0,0,0,0,0,0};
    std::array <int,9>    acceptances = {0,0,0,0,0,0,0,0,0}; 
    std::array <double,4> contacts    = {0,0,0,0}; 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // Parse inputs... 
    // This command will take all of the above inputs and make sure they are valid. 
    InputParser ( dfreq, max_iter, r, positions, \
        topology, dfile, efile, mfile, stats_file, \
        lattice_file_read ); 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // ExtractNumberOfPolymers extracts the total number of chains in the input file. used to reserve vector [line 171] 
    const int N = ExtractNumberOfPolymers(positions); 
    
    // EXTRACT TOPOLOGY FROM FILE 
    std::array <double,8> info_vec {ExtractTopologyFromFile(topology)}; 
    const int x             =  info_vec[0];
    const int y             =  info_vec[1]; 
    const int z             =  info_vec[2]; 
    const double T          =  info_vec[3]; 
    std::array <double,4> E = {info_vec[4], info_vec[5], info_vec[6], info_vec[7]}; 
    
    // initialize custom data structures 
    std::vector <Polymer> Polymers; 
    Polymers.reserve(N);

    std::vector <Particle*> LATTICE;
    LATTICE.reserve (x*y*z); 
    
    std::vector <Polymer> Polymers_target = ExtractPolymersFromFile (target, x, y, z); 

   // compare the two polymers 
   for (int i{0 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // OPENING TILES

    std::cout << std::endl;
    std::cout << "Preparing for take-off...\n\n" ; 
    std::cout << "Chemical information: " << std::endl;
    std::cout << "Number of polymers is " << N << ".\n\n";
    std::cout << "Geometric information about simulation cell: " << std::endl;
    std::cout << "x = " << x <<", y = " << y << ", z = "<< z << "." << std::endl << std::endl;
    std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
    std::cout << "Temperature = " << T << "." << std::endl; 
    std::cout << "Emm_a = " << E[0] <<", Emm_n = " << E[1] << ", Ems_a = "<< E[2] << ", Ems_n = " << E[3] <<".\n \n";  
    std::cout << "Off to a good start. \n\n";
    std::cout << "--------------------------------------------------------------------\n" << std::endl;
    std::cout << "Running some more checks on input... \n\n" ; 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    
    // set timers for simulation set-up.  
    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 


    if ( !r ){
        std::cout << "Setting up the lattice from scratch! " << std::endl;
        SetUpLatticeFromScratch (x, y, z, &Polymers, &LATTICE, positions); 
    
        stop = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

        std::cout << "Solvation took " << duration.count () << " milliseconds." << std::endl;
        std::cout << "Cell has been solvated! \n\n" ;
    
        dumpPositionsOfPolymers(&Polymers, step_number, dfile); 
    }

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    else {
        std::cout << "Setting up system from a restart file!" << std::endl;
        SetUpLatticeFromRestart (x, y, z, &Polymers, &LATTICE, step_number, lattice_file_read, dfile, positions ); 
        
        stop = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

        std::cout << "System set-up took " << duration.count () << " milliseconds." << std::endl;
        std::cout << "Simulation cell has been made! " << std::endl << std::endl;

    }    
    
    // THERMODYNAMICS OF SET-UP

    std::cout <<"\nCalculating energy..." << std::endl;

    sysEnergy = CalculateEnergy(&Polymers, &LATTICE, &E, &contacts, x, y, z); 

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

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    // ALL INPUTS HAVE BEEN PROCESSED AND DATA STRUCTURES HAVE BEEN SET UP. 
    
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    // BEGIN: main loop of simulation! yee-haw!  
    double sysEnergy_ = sysEnergy;
    std::array <double,4> contacts_  = contacts; 
    int deg_poly = static_cast<int> (Polymers[0].chain.size()); 

    int hitcheck {1}; 
    int nhits {0}; 
    std::vector <std::array<int,3>> old_cut;
    std::vector <std::array<int,3>> new_cut; 
    
    int m_index = 5; 

    for ( int i{0}; i<deg_poly-m_index-1; ++i ){
        old_cut.push_back ( (Polymers)[0].chain[i+m_index+1]->coords ); 
    }

    std::vector <std::array<int,3>> target = {{6,0,0}, {7,0,0}};

    for (int i = step_number+1; i< (step_number+max_iter+1); ++i) {

        if ( v && (i%dfreq==0) ){
            std::cout << "Initial config is:" << std::endl;
            Polymers[0].printChainCoords();
            std::cout << " -------------------------- " << std::endl;
            std::cout << "Step number: " << i << "." << std::endl; 
            std::cout << "Executing..." << std::endl << std::endl;
        }

        // perform move on the system! 
        PerturbSystem (&Polymers, &LATTICE, &E, &contacts, &attempts, &IMP_BOOL, v, &sysEnergy, T, &move_number, x, y, z); 
        new_cut.clear();
        for ( int y{0}; y<deg_poly-m_index-1; ++y ){
            new_cut.push_back ( (Polymers)[0].chain[y+m_index+1]->coords ); 
        }

        if ( v ){
            if (IMP_BOOL){
                std::cout << "IMP_BOOL = " << IMP_BOOL << std::endl;
                std::cout << "perturbed config is:" << std::endl;
                Polymers[0].printChainCoords();
                std::cout << "Accepted!" << std::endl;
            }
            else {
                std::cout << "IMP_BOOL = " << IMP_BOOL << std::endl;
                std::cout << "Rejected..." << std::endl;   
            }
            std::cout << "Checking if data structures are in good conditions..." << std::endl; 
            CheckStructures (x, y, z, &Polymers, &LATTICE);
        }

        if ( ( i % dfreq == 0 ) ){
            dumpPositionsOfPolymers (&Polymers, i, dfile); 
            if ( i % (dfreq*10) == 0 ) {
                dumpOrientation (&Polymers, &LATTICE, i, mfile, x, y, z); 
            }
            dumpEnergy (sysEnergy, i, &contacts, efile);
        }

        // std::cout << "start running checks... " << std::endl;
        // run checks 
        if (IMP_BOOL){
            (acceptances)[move_number] += 1;  
            for (int s{m_index+1}; s<deg_poly; ++s){
                if ( (Polymers)[0].chain[s]->coords != target[s-(m_index+1)] ){
                    hitcheck = 0;
                    break; 
                }
            }
            if (hitcheck){
                nhits += 1; 
            }

            acceptance_after_head_regrowth ( &LATTICE, &new_cut, &old_cut, y, z);
        }

        sysEnergy = sysEnergy_;
        contacts  = contacts_;
        IMP_BOOL = true;     
        hitcheck = 1;  
    }


    dumpMoveStatistics (&attempts, &acceptances, max_iter, stats_file);  
    std::cout << "Number of hits = " << nhits << std::endl;

    if ( lattice_file_write != "__blank__" ) {
        dumpLATTICE ( &LATTICE, step_number+max_iter, y, z, lattice_file_write ); 
    }

    stop = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 
	
    std::cout << "\n\nTime taken for simulation: " << duration.count() << " milliseconds." << std::endl; 
    std::cout << "That is all she wrote. Hope it worked." << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;

    return 0;

}
