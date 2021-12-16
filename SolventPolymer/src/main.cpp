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


    int opt; 
    int Nmov {-1}, dfreq {-1};
    std::string positions {"blank"}, topology {"blank"};  

    while ( (opt = getopt(argc, argv, ":f:N:p:t:h")) != -1 )
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
                "These are all the options we have available right now: \n" <<
                "help                     [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "Dump Frequency           [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<
                "Number of MC moves       [-N]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" << 
                "Position coordinates     [-p]           (STRING ARGUMENT REQUIRED)     File with position coordinates\n" <<
                "Energy and geometry      [-t]           (STRING ARGUMENT REQUIRED)     File with energetic interactions and geometric bounds ie the topology\n"; 
                exit(EXIT_SUCCESS);
                break;


            case 'p':
                // std::cout <<"Option p was called with argument " << optarg << std::endl;
                positions = optarg;
                break;    

            case 't':
                topology = optarg; 
                break;


            case '?':
                std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
                exit(EXIT_FAILURE); 
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
    else if (positions=="blank" || topology == "blank"){
        std::cerr << "ERROR: No value for option p (positions file) and/or for option t (energy and geometry file) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);    
    }

    
    // ---------------------------------------

    std::ofstream dump_file ("dumpfile.txt");
    dump_file.close(); 

    // ----------------------------------------


    Grid G = CreateGridObject(positions, topology);

    G.instantiateOccupancyMap();
    
    G.CalculateEnergy();
    std::cout << "Energy of box is: " << G.Energy << std::endl;

    return 0;

}
