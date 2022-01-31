#include <iostream> 
#include <fstream>
#include <vector> 
#include <string> 
#include <array>
#include <map>
#include <array>
#include <random>
#include <chrono>
#include <getopt.h> 
#include <stdlib.h> 
#include "classes.h"
#include "misc.h"


int main(int argc, char** argv) {

    // set up 
    int opt; 
    int Nacc {-1}, dfreq {-1}, max_iter{-1};
    std::string positions {"blank"}, topology {"blank"}, dfile {"blank"}, efile{"blank"}, restart_traj{"blank"};  
    bool v = false, a = false, r = false;

    while ( (opt = getopt(argc, argv, ":f:M:N:T:o:u:p:t:e:vhar")) != -1 )
    {
        switch (opt) 
        {
            case 'f':
                // std::cout << "Option dreq was called with argument " << optarg << std::endl; 
                dfreq = atoi(optarg); 
                break;



            case 'N':
                // std::cout << "Option Nmov was called with argument " << optarg << std::endl; 
                Nacc = atoi(optarg); 
                break; 

            case 'M':
                max_iter = atoi(optarg); 
                break; 

            case 'h':
                std::cout << 
                "This is the main driver for the monte carlo simulation of polymers in a box.\n" <<
		        "This version number 0.2.1 of the Monte Carlo Engine. Set up on Jan 26, 2022, 04:35 PM.\n" <<
                "These are all the options we have available right now: \n" <<
                "help                     [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose                  [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
                "Data only for accepts    [-a]           (NO ARG REQUIRED)              If you only want energy and coords for every accepted structure, use this option. \n"
                "Restart simulation       [-r]           (NO ARG REQUIRED)              Pick up a simulation back from some kind of a starting point.\n"
                "Dump Frequency           [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
                "Number of maximum moves  [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
                "Required accepted moves  [-N]           (INTEGER ARGUMENT REQUIRED)    Number of accepted moves for a good simulation.\n" <<  
                "Position coordinates     [-p]           (STRING ARGUMENT REQUIRED)     File with position coordinates.\n" <<
                "Energy of grid           [-u]           (STRING ARGUMENT REQUIRED)     Dump energy of grid at each step in a file.\n"<<
                "Energy and geometry      [-t]           (STRING ARGUMENT REQUIRED)     File with energetic interactions and geometric bounds ie the topology.\n" <<
                "Previous trajectory file [-T]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n"<<
                "Name of output file      [-o]           (STRING ARGUMENT REQUIRED)     Name of file which will contain coordinates of polymer.\n";  
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

            case 'T':
                restart_traj = optarg;
                break;

            case 'u':
                efile = optarg;
                break;


            case '?':
                std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
                exit(EXIT_FAILURE); 
                break;

            case 'v':
                std::cout << "Output to console will be verbose. " << std::endl;
                v = true;
                break;

            case 'a':
                std::cout <<"Only accepted structures will be outputted." << std::endl;
                a = true; 
                break; 

            case 'r':
                std::cout <<"Will attempt to restart simulation by taking coordinates from a previous trajectory file." << std::endl;
                r = true; 
                break; 

            case ':':
                std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
                exit(EXIT_FAILURE);
                
                
                break; 
        }
    }

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // parse inputs 
    
    InputParser (a, r, Nacc, dfreq, max_iter, positions, topology, dfile, efile, restart_traj); 

    // driver 

    auto start = std::chrono::high_resolution_clock::now(); 

    // Grid G;                             // setting up the grid object
    // bool call {true};                   // setting up the call for first input or not 
    // int step_number;                    // defining step_number for output reasons 

    std::array <double,7> info_vec {ExtractTopologyFromFile(topology)}; 
    const int N = ExtractNumberOfPolymers(positions); 

    const int x = info_vec[0];
    const int y = info_vec[1]; 
    const int z = info_vec[2]; 
    const double T = info_vec[3]; 
    const double Emm_a = info_vec[4]; 
    const double Emm_n = info_vec[5]; 
    const double Ems = info_vec[6]; 

    std::cout << "Number of polymers is " << N << ".\n";
    std::cout << "x = " << x <<", y = " << y << ", z = "<< z << ", T = " << T << ".\n"; 
    std::cout << "Emm_a = " << Emm_a <<", Emm_n = " << Emm_n << ", Ems = "<< Ems << ".\n";  

    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~" << std::endl;
    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~\n" << std::endl;

    int step_number = 0;
    double sysEnergy {0}; 
    std::vector <Polymer> PolymerVector; 
    PolymerVector.reserve(N); 
    if (r){
        std::cout << "Restart mode activated.\n";
        PolymerVector = ExtractPolymersFromTraj(restart_traj, positions, x, y, z); 
        step_number = ExtractIndexOfFinalMove(restart_traj);
    }
    else {

        PolymerVector = ExtractPolymersFromFile(positions, x, y, z); 
        dumpPositionsOfPolymers(&PolymerVector, step_number, dfile); 
        sysEnergy = CalculateEnergy (&PolymerVector, x, y, z, Emm_a, Emm_n, Ems); 
        dumpEnergy (sysEnergy, step_number, efile, true);
    }
    
    double sysEnergy_ {0}; 
    std::cout << "Energy of system is " << sysEnergy << std::endl;

     
	int acc_counter = 0; 
    
    

    bool IMP_BOOL = true; 

    std::vector <Polymer> PolymerVector_; 

    printf("Simulation will output information of every %d configuration.\n", dfreq); 
    for (int i = step_number+1; i< (step_number+max_iter+1); i++) {

        if ( v && (i%dfreq==0) ){
            printf("Move number %d.\n", i);
        }
        // choose a move 
        PolymerVector_ = MoveChooser(&PolymerVector, x, y, z, v, &IMP_BOOL);  
        sysEnergy_ = CalculateEnergy(&PolymerVector_, x, y, z, Emm_a, Emm_n, Ems); 
        if ( v && (i%dfreq==0) ){
            printf("Executing...\n");
        }
        

        if ( MetropolisAcceptance (sysEnergy, sysEnergy_, T) && IMP_BOOL ) {

            // replace old config with new config
            if ( v ){
                printf("Checking validity of coords...");
                printf("checkForOverlaps says: %d.\n", checkForOverlaps(PolymerVector_)); 
                if (!checkForOverlaps(PolymerVector_)){
                    printf("Something is fucked up overlaps-wise. \n");
                    exit(EXIT_FAILURE);
                }
                printf("checkConnectivity says: %d\n", checkConnectivity(PolymerVector_, x, y, z)); 
                if (! checkConnectivity(PolymerVector_, x, y, z) ){
                    printf("Something is fucked up connectivity-wise. \n");
                    exit(EXIT_FAILURE);
                }
                printf("Accepted.\n");
                printf("Energy of the system is %.2f.\n", sysEnergy_);
                printf("This should be 1 as IMP_BOOL must be true on acceptance: %d\n", IMP_BOOL);
            }

            // std::cout << "IMP_BOOL is " << IMP_BOOL << std::endl;
            
            PolymerVector = std::move(PolymerVector_);
            sysEnergy = sysEnergy_; 
			++acc_counter;
        }

        //if (!IMP_BOOL){
        //    std::cout << "IMP_BOOL is " << IMP_BOOL << std::endl;            
        //}


        else {
            if ( v && (i%dfreq==0) ){
                printf("Not accepted.\n");
                printf("Energy of the suggested system is %.2f, while energy of the initial system is %.2f.\n", sysEnergy_, sysEnergy);
            }
            
        }

        if ( ( i % dfreq == 0) ){
            
            dumpPositionsOfPolymers(&PolymerVector, i, dfile); 
            dumpEnergy (sysEnergy, i, efile, false);

            
        }
        // std::cout <<"are you hitting  this line?"<<std::endl;
        IMP_BOOL = true; 
    }
    
    auto stop = std::chrono::high_resolution_clock::now(); 
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 
	
	printf("\nNumber of moves accepted is %d.", acc_counter);
    printf("\n\nTime taken for simulation: %lld milliseconds\n", duration.count() ); 

    return 0;

}
