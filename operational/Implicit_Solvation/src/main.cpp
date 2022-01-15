#include <iostream> 
#include <fstream>
#include <vector> 
#include <string> 
#include <map>
#include <random>
#include <chrono>
#include <getopt.h> 
#include <stdlib.h> 
#include "classes.h"
#include "misc.h"


int main(int argc, char** argv) {

    // set up 
    int opt; 
    int Nmov {-1}, dfreq {-1};
    std::string positions {"blank"}, topology {"blank"}, dfile {"blank"}, efile{"energydump.txt"};  
    bool v = false;

    while ( (opt = getopt(argc, argv, ":f:N:o:u:p:t:e:vh")) != -1 )
    {
        switch (opt) 
        {
            case 'f':
                // std::cout << "Option dreq was called with argument " << optarg << std::endl; 
                dfreq = atoi(optarg); 
                break;



            case 'N':
                // std::cout << "Option Nmov was called with argument " << optarg << std::endl; 
                Nmov = atoi(optarg); 
                break; 


            case 'h':
                std::cout << 
                "This is the main driver for the monte carlo simulation of polymers in a box.\n" <<
		        "This version number 1.1.1 of the Monte Carlo Engine. Set up on Jan 14, 2022, 12:30 AM.\n" <<
                "This version of the engine has no std::swap in the copy constructor for the Grid object.\n" <<
                "This version implements pointers for Grid based functions.\n" << 
                "These are all the options we have available right now: \n" <<
                "help                     [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose                  [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. Usually meant to debug. \n"<<
                "Dump Frequency           [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
                "Number of MC moves       [-N]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" << 
                "Position coordinates     [-p]           (STRING ARGUMENT REQUIRED)     File with position coordinates.\n" <<
                "Energy of grid           [-u]           (STRING ARGUMENT REQUIRED)     Dump energy of grid at each step in a file.\n"<<
                "Energy and geometry      [-t]           (STRING ARGUMENT REQUIRED)     File with energetic interactions and geometric bounds ie the topology.\n" <<
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


            case ':':
                
                std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << std::endl;
                exit(EXIT_FAILURE);
                
                
                break; 
        }
    }

    
    if (Nmov == -1 || dfreq == -1){
        std::cerr << "ERROR: No value for option N (number of MC moves to perform) and/or for option f (frequency of dumping) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);
    }
    else if (positions=="blank" || topology == "blank" || dfile=="blank"){
        std::cerr << "positions is " << positions <<", topology is " << topology <<", dfile is " << dfile << std::endl;
        std::cerr << "ERROR: No value for option p (positions file) and/or for option t (energy and geometry file) and/or for option o (name of output dump file) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);    
    }

    if (v){
        std::cout << "\nVERBOSE OUTPUT HAS BEEN TOGGLED.\n" << std::endl;
    }

    std::ofstream dump_file (dfile);
    std::ofstream energy_dump_file (efile);
    energy_dump_file.close(); 
    dump_file.close(); 

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
    std::cout << "monomer-monomer aligned interaction is " << G.Emm_a <<"\nmonomer-monomer misaligned interaction is " << G.Emm_n <<"\nmonomer-solvent " << G.Ems
    << std::endl; 
    
    if (v){
        std::cout << "Energy of box is: " << G.Energy << std::endl;
        std::cout << "Next..." << std::endl;
    }


    

    Grid G_ ;

    bool IMP_BOOL = true; 
    int acceptance_count = 0; 
    for (int i{1}; i< (Nmov+1); i++){


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
            }
			acceptance_count++; 
            G = std::move(G_);
        }


        else {
            if ( v && (i%dfreq==0) ){
                printf("Not accepted.\n");
                printf("Energy of the suggested system is %f, while energy of the initial system is %f.\n", G.Energy, G_.Energy);
            }
            // continue;
        }

        if (i % dfreq == 0){
            G.dumpPositionsOfPolymers (i, dfile) ;
            G.dumpEnergyOfGrid(i, efile, call) ; 
        }
        // G.PolymersInGrid.at(0).printChainCoords();

        IMP_BOOL = true; 
    }
    
    
    auto stop = std::chrono::high_resolution_clock::now(); 
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

    std::cout << "\n\nTime taken for simulation: " << duration.count() << " milliseconds" << std::endl;
    std::cout << "Number of acceptances is " << acceptance_count << std::endl;
    
    return 0;

}
