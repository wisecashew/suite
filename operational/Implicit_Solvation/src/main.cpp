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
    std::string positions {"blank"}, topology {"blank"}, dfile {"blank"}, efile{"energydump.txt"}, restart_traj{"blank"};  
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
		        "This version number 0.2.1 of the Monte Carlo Engine. Set up on Jan 18, 2022, 11:04 PM.\n" <<
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
                "Previous trajectory file [-T]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation. Can only be used with -r flag.\n"<<
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
            case ':':
                
                std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << std::endl;
                exit(EXIT_FAILURE);
                
                
                break; 
        }
    }

    
    // Check what kind of statistics do you want! 
    // With option a, you only sample accepted positions. Do this only if you are not getting decorrelated chains! 


    if (!a) {
        if (Nacc != -1){
            std::cerr << "ERROR: You do not need to provide a -N option if you are not looking for accepted configuration statistics. Use -a if you want to use -N. Safeguarding against uncontrolled behavior. Exiting..." << std::endl;
            exit (EXIT_FAILURE); 
        }

        if (dfreq == -1 || max_iter == -1){
            std::cerr << "ERROR: No value for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
            exit (EXIT_FAILURE);
        }
    }

    if (!r) {
        if (restart_traj != "blank"){
            std::cerr << "ERROR: You cannot ask for a trajectory coordinate file without the -r flag. Exiting..." << std::endl;
            exit (EXIT_FAILURE); 
        }
    }

    else {
        if (Nacc == -1 || dfreq == -1 || max_iter == -1 ){
            std::cerr << "ERROR: No value for option N (number of accepted MC moves to have) and/or for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
            exit (EXIT_FAILURE);
        }
    }


    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#   


    if (positions=="blank" || topology == "blank" || dfile=="blank"){
        std::cerr << "positions is " << positions <<", topology is " << topology <<", dfile is " << dfile << std::endl;
        std::cerr << "ERROR: No value for option p (positions file) and/or for option t (energy and geometry file) and/or for option o (name of output dump file) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);    
    }

    if (v){
        std::cout << "\nVERBOSE OUTPUT HAS BEEN TOGGLED.\n" << std::endl;
    }
     
    if (r){
        std::ofstream dump_file(dfile, std::ios::app); 
        std::ofstream energy_dump_file (efile, std::ios::app); 
    }

    if (!r){
        std::ofstream dump_file (dfile);
        std::ofstream energy_dump_file (efile);
    }

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    // driver 

    auto start = std::chrono::high_resolution_clock::now(); 

    Grid G = CreateGridObject(positions, topology);
    
    std::cout << "Temperature of box is " << G.kT << "." << std::endl;
    G.dumpPositionsOfPolymers(0, dfile);
    bool call {true}; 
    G.dumpEnergyOfGrid(0, efile, true); 
    call = false; 
    G.CalculateEnergy();


    std::cout << "Energy surface check: " << std::endl; 
    std::cout << "monomer-monomer aligned interaction is " << G.Emm_a <<
    "\nmonomer-monomer misaligned interaction is " << G.Emm_n <<"\nmonomer-solvent " << G.Ems
    << std::endl; 
    
    if (v){
        std::cout << "Energy of box is: " << G.Energy << std::endl;
        std::cout << "Next..." << std::endl;
    }

    
    

    Grid G_ ;

    bool IMP_BOOL  {true}; 

    if (a) {
    printf("Simulation will output only information about accepted configurations.\n"); 

    int acceptance_count {0} ; 
	int temp_var {-1};							// this variable is to make sure acceptance_count is not recounted  
    for (int i{1}; i< (max_iter+1); i++) {


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

    printf("\n\nTime taken for simulation: %ld milliseconds\n", duration.count() ); //  << duration.count() << " milliseconds" << std::endl;
    // printf("Number of acceptances is %d.\n", acceptance_count) ; //  << acceptance_count << std::endl;
    }



    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#



    else {
        printf("Simulation will output information of every %d configuration.\n", dfreq); 
        for (int i{1}; i< (max_iter+1); i++) {


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
                // ++acceptance_count; 
                G = std::move(G_);
            }


            else {
                if ( v && (i%dfreq==0) ){
                    printf("Not accepted.\n");
                    printf("Energy of the suggested system is %f, while energy of the initial system is %f.\n", G.Energy, G_.Energy);
                }
                // continue;
            }

            if ( ( i % dfreq == 0) ){
                // printf("acceptance_count is %d and dfreq is %d.\n", acceptance_count, dfreq); 
                G.dumpPositionsOfPolymers (i, dfile) ;
                G.dumpEnergyOfGrid(i, efile, call) ; 
                // temp_var = acceptance_count ;
            }

            IMP_BOOL = true; 
        }
    
    
    auto stop = std::chrono::high_resolution_clock::now(); 
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

    printf("\n\nTime taken for simulation: %ld milliseconds\n", duration.count() ); 
   
    }


    return 0;

}
