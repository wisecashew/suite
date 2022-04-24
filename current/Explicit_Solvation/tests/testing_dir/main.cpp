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
                "verbose                  [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. Usually meant to debug. \n"<<
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

    Grid G;                             // setting up the grid object
    bool call {true};                   // setting up the call for first input or not 
    int step_number;                    // defining step_number for output reasons 


    // define grid objects 

    if (r){
        G = CreateGridObjectRestart(positions, topology, restart_traj);
        step_number = G.ExtractIndexOfFinalMove (restart_traj) ; 
        G.CalculateEnergy();   
        call = false; 
    }
    else {
        G = CreateGridObject(positions, topology);
        std::cout << "Temperature of box is " << G.kT << "." << std::endl;
        step_number = 0; 
        G.dumpPositionsOfPolymers(step_number, dfile); 
        G.dumpEnergyOfGrid(step_number, efile, true); 
        call = false; 
        G.CalculateEnergy();
    }


    std::cout << "Energy surface check: " << std::endl; 
    std::cout << "monomer-monomer aligned interaction is " << G.Emm_a <<
    "\nmonomer-monomer misaligned interaction is " << G.Emm_n <<"\nmonomer-solvent " << G.Ems
    << std::endl; 
    
    if (v){
        std::cout << "Energy of box is: " << G.Energy << std::endl;
        std::cout << "Next..." << std::endl;
    }

    // grid objects have been validated. 
    

    Grid G_ ;

    bool IMP_BOOL  {true}; 

    if (a) {
    printf("Simulation will output only information about accepted configurations.\n"); 

    int acceptance_count {0} ; 
	int temp_var {-1};							// this variable is to make sure acceptance_count is not recounted  
    for (int i = step_number; i< (step_number+max_iter+1); i++) {


        if ( v && (i%dfreq==0) ){
            printf("Move number %d.\n", i);
        }
        // choose a move 
        G_ = MoveChooser(&G, v, &IMP_BOOL);  

        if ( v && (i%dfreq==0) ){
            printf("Executing...\n");
        }


        if ( MetropolisAcceptance (G.Energy, G_.Energy, G.kT) && IMP_BOOL ) {
            // accepted
            // replace old config with new config
            if ( v ){ 
                printf("Accepted.\n");
                printf("Energy of the system is %f.\n", G_.Energy);
                printf("%d\n", IMP_BOOL);
            }
			++acceptance_count; 
            G = std::move(G_);
        }


        else {
            if ( v && (i%dfreq==0) ){
                printf("Not accepted.\n");
                printf("Energy of the suggested system is %f, while energy of the initial system is %f.\n", G.Energy, G_.Energy);
            }
            // continue;
        }

        if ( ( acceptance_count % dfreq == 0) && ( temp_var != acceptance_count ) ){
			printf("acceptance_count is %d and dfreq is %d.\n", acceptance_count, dfreq); 
            G.dumpPositionsOfPolymers (i, dfile) ;
            G.dumpEnergyOfGrid(i, efile, call) ; 
			temp_var = acceptance_count ;
        }
		// printf("acceptance_count is ACTUALL %d ON LINE 199.\n", acceptance_count); 
        // G.PolymersInGrid.at(0).printChainCoords();

        if (acceptance_count == Nacc){
            printf("Simulation has accepted %d moves. Required to stop after accepting %d moves.\n", acceptance_count, Nacc);
			printf("Total number of suggested moves is %d.\n", i);
			break; 
        }

        IMP_BOOL = true; 
    }
    
    
    auto stop = std::chrono::high_resolution_clock::now(); 
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

    printf("\n\nTime taken for simulation: %ld milliseconds\n", static_cast<long>( duration.count() ) );

    }


    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#


    else {

        printf("Simulation will output information of every %d configuration.\n", dfreq); 
        for (int i = step_number+1; i< (step_number+max_iter+1); i++) {

            

            if ( v && (i%dfreq==0) ){
                printf("Move number %d.\n", i);
            }
            // choose a move 
            G_ = MoveChooser(&G, v, &IMP_BOOL);  

            if ( v && (i%dfreq==0) ){
                printf("Executing...\n");
            }
            

            if ( MetropolisAcceptance (G.Energy, G_.Energy, G.kT) && IMP_BOOL ) {
                // accepted
                // replace old config with new config
                if ( v ){ 
                    printf("Accepted.\n");
                    printf("Energy of the system is %f.\n", G_.Energy);
                    printf("This should be 1 as IMP_BOOL must be true on acceptance: %d\n", IMP_BOOL);
                }

                G = std::move(G_);
            }


            else {
                if ( v && (i%dfreq==0) ){
                    printf("Not accepted.\n");
                    printf("Energy of the suggested system is %f, while energy of the initial system is %f.\n", G_.Energy, G.Energy);
                }
                
            }

            if ( ( i % dfreq == 0) ){
                
                G.dumpPositionsOfPolymers (i, dfile) ;
                G.dumpEnergyOfGrid(i, efile, call) ; 
                
            }
            // std::cout <<"are you hitting  this line?"<<std::endl;
            IMP_BOOL = true;
	    std::cout << "Printing out occupancy map..." << std::endl;
	    int nkey = 0;  
	    for (auto it = G.OccupancyMap.begin(); it != G.OccupancyMap.end(); it++){
	    	std::cout << "Key is: ";
		print(it->first);
		int n{0}; 
	        for (auto p: G.PolymersInGrid[0].chain){
			
			if (p.coords == it->first){
				std::cout << "Found key." << std::endl;
				break; 
			}
			++n;
		}
		std::cout << "n is " << n << std::endl;
		if (static_cast<size_t>(n)==G.PolymersInGrid[0].chain.size()){
			std::cout << "Something is fucked." << std::endl; 
			exit (EXIT_FAILURE); 	
	    	}
		++nkey; 
		
	    } 
	    std::cout << "number of keys is " << nkey << std::endl; 
        }
    
    
    auto stop = std::chrono::high_resolution_clock::now(); 
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

    printf("\n\nTime taken for simulation: %ld milliseconds\n", static_cast<long>(duration.count() ) ); 
   
    }


    return 0;

}
