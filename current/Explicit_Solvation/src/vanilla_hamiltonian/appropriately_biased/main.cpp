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


int main (int argc, char** argv) {

    // set up 
    int opt; 
    int dfreq {-1}, max_iter{-1};
    std::string positions {"__blank__"}, topology {"__blank__"}, dfile {"__blank__"}, efile{"__blank__"}, mfile {"__blank__"}, stats_file {"__blank__"}, \
    lattice_file_write {"__blank__"}, lattice_file_read {"__blank__"};  
    bool v = false, r = false;

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
                "\nWelcome to my Monte Carlo simulation engine (v0.6) for polymers in a simple cubic box (Z=26). \nThis is a simulation engine which incorporates directional bonding between monomer and solvent. Here, we incorporate the relative orientations of the particles affecting the energy of interaction. \nIn this implementation, I have employed reversing moves to avoid copying. Did wonders for efficiency." <<
		        "\nLast updated: Aug 1, 2022, 12:56 AM. Added restarting capabilities. \nAuthor: satyend@princeton.edu\n" <<
                "\n----------------------------------------------------------------------------------------------------------------------------------\n" << 
                "These are all the inputs the engine accepts for a single run, as of right now:\n\n" <<
                "help                      [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose flag              [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
                "restart flag              [-r]           (NO ARG REQUIRED)              Restarts simulation from final spot of a previous simulation. \n"<<
                "Dump Frequency            [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
                "Number of maximum moves   [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
                "Polymer coordinates       [-p]           (STRING ARGUMENT REQUIRED)     Name of input file with coordinates of polymer.\n" <<
                "Energy and geometry       [-t]           (STRING ARGUMENT REQUIRED)     Name of input file with energetic interactions and geometric bounds.\n" <<
                // "Solvent coordinates       [-S]           (STRING ARGUMENT REQUIRED)     Name of output file with coordinates of solvent. \n" <<
                "Energy of grid            [-u]           (STRING ARGUMENT REQUIRED)     Name of output file with energy of system at each step in a file.\n"<<
                "Lattice file to write to  [-L]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
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

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // INSTANTIATING NECESSARY DATA STRUCTURES UP FRONT.

    int    step_number  {0}; // step_number++;
    double sysEnergy    {0}; // sysEnergy++;
    std::array <int,9>    attempts    = {0,0,0,0,0,0,0,0,0};
    std::array <int,9>    acceptances = {0,0,0,0,0,0,0,0,0}; 
    std::array <double,8> contacts    = {0,0,0,0,0,0,0,0}; 

    std::vector <Polymer> Polymers; 
    Polymers.reserve(N);

    std::vector <Particle*> LATTICE;
    LATTICE.reserve (x*y*z); 


    bool IMP_BOOL = true; 
    
    double rweight =  0; 
    int move_number = 0; 
    int monomer_index = -1; 
    int back_or_front = -1; 


    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // Parse inputs... 
    // This command will take all of the above inputs and make sure they are valid. 
    InputParser ( dfreq, max_iter, r, positions, \
        topology, dfile, efile, mfile, stats_file, \
        lattice_file_read ); 


    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    
    // GO! 

    // ################################################################################
    // ExtractNumberOfPolymers extracts the total number of chains in the input file 
    const int N = ExtractNumberOfPolymers(positions); 
    
    // EXTRACT TOPOLOGY FROM FILE 
    std::array <double,8> info_vec {ExtractTopologyFromFile(topology)}; 
    const int x             = info_vec[0];
    const int y             = info_vec[1]; 
    const int z             = info_vec[2]; 
    const double T          = info_vec[3]; 
    std::array <double,4> E = {info_vec[4], info_vec[5], info_vec[6], info_vec[7]}; 
    

    // ################################################################################ 
    // OPENING TILES
    // 

    std::cout << std::endl;
    std::cout << "Preparing for take-off...\n\n" ; 
    std::cout << "Chemical information: " << std::endl;
    std::cout << "Number of polymers is " << N << ".\n\n";
    std::cout << "Geometric information about simulation cell: " << std::endl;
    std::cout << "x = " << x <<", y = " << y << ", z = "<< z << ".\n" << std::endl;
    std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
    std::cout << "Temperature = " << T << "." << std::endl; 
    std::cout << "Emm_a = " << E[0] <<", Emm_n = " << E[1] << ", Ems_a = "<< E[2] << ", Ems_n = " << E[3] <<".\n \n";  
    std::cout << "Off to a good start. \n\n";
    std::cout << "--------------------------------------------------------------------\n" << std::endl;
    std::cout << "Running some more checks on input... \n\n" ; 

    // ################################################################################ 
    
    // SET TIMERS UP! 
    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 


    if ( !r ){
        std::cout << "Setting up the lattice from scratch! " << std::endl;
        SetUpLatticeFromScratch (x, y, z, &Polymers, &LATTICE, &positions); 
    
        stop = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

        std::cout << "Solvation took " << duration.count () << " milliseconds." << std::endl;
        std::cout << "Cell has been solvated! \n\n" ;
    
        dumpPositionsOfPolymers(&Polymers, step_number, dfile); 
    }

    // ######################################################################

    else {
        std::cout << "Setting up system from a restart file!" << std::endl;
        SetUpLatticeFromRestart (x, y, z, &Polymers, &LATTICE, step_number, lattice_file_read, dfile, positions ); 
        
        stop = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

        std::cout << "System set-up took " << duration.count () << " milliseconds." << std::endl;
        std::cout << "Simulation cell has been made! \n\n" ;

    }    
    
    // THERMODYNAMICS OF SET-UP

    std::cout <<"\nCalculating energy..." << std::endl;
    sysEnergy = CalculateEnergy(&Polymers, &LATTICE, x, y, z, &E, &contacts); 
    std::cout << "Energy of system is " << sysEnergy << ".\n" << std::endl;
    
    // if i am not restarting, i do not need to dump anything. All the information is already present. 
    if (!r) {
        dumpEnergy      (sysEnergy, step_number, &contacts, efile); 
        dumpOrientation (&Polymers, &LATTICE, step_number, mfile, x, y, z); 
    }
    
    else {
        dumpOrientation (&Polymers, &LATTICE, step_number, mfile, x, y, z);    
    }
    
    printf("Initiation complete. We are ready to go. The engine will output information every %d configuration(s).\n", dfreq); 
    
    std::cout << "Number of iteration to perform: " << max_iter << "." << std::endl;

    for (int i = step_number+1; i< (step_number+max_iter+1); ++i) {

        if ( v && (i%dfreq==0) ){
            printf("Step number: %d.\n", i);
            printf("Executing...\n");
        }

        // choose a move... 
        PerturbSystem_BIASED (&Polymers, &LATTICE, index, x, y, z, &sysEnergy, temperature, &move_number, &IMP_BOOL, &attempts); 
        
        if ( IMP_BOOL ) {
            if ( MetropolisAcceptance (sysEnergy, sysEnergy_, T, rweight) ){
                metropolis = true; 
                acceptances[move_number-1]+=1;

                if ( v ){
                    CheckStructures (x, y, z, &IMP_BOOL, &Polymers, &LATTICE);
                }
            }

            else { 

                ReversePerturbation (&Polymers, &LATTICE, y, z, v, move_number, &memory3, &memory2, monomer_index, back_or_front);
                
                if ( v ){
                    CheckStructures (x, y, z, &IMP_BOOL, &Polymers, &LATTICE); 
                }
            }
        }

        else {
            if (v){
                printf ("IMP_BOOL is zero. Nothing will be done.\n\n");
                CheckStructure (x, y, z, &IMP_BOOL, &Polymers, &LATTICE); 
            }
        }

        if ( ( i % dfreq == 0) ){
            dumpPositionsOfPolymers (&Polymers, i, dfile); 
            if ( i%(dfreq*10) == 0 ) {
                dumpOrientation (&Polymers, &LATTICE, i, mfile, x, y, z); 
            }
            dumpEnergy (sysEnergy, i, mm_aligned_copy, mm_naligned_copy, ms_aligned_copy, ms_naligned_copy, efile);
        }

        IMP_BOOL = true; 
           
    }

    dumpMoveStatistics      (&attempts, &acceptances, max_iter, stats_file);  
    
    if ( lattice_file_write != "__blank__" ) {
        dumpLATTICE ( &LATTICE, step_number+max_iter, y, z, lattice_file_write ); 
    }

    stop = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 
	
    std::cout << "\n\nTime taken for simulation: " << duration.count() << " milliseconds.\n"; 
    std::cout << "That is all she wrote. Hope it worked." << std::endl;
    std::cout << "--------------------------------------------------------------------\n\n";

    return 0;

}
