#include <iostream> 
#include <vector> 
#include <string> 
#include <map>
#include <chrono>
#include <getopt.h> 
#include <stdlib.h> 
#include "classes.h"
#include "misc.h"

int main(int argc, char** argv) {

    int opt; 
    int Nmov, dfreq;
    std::string positions, energy;  

    while ( (opt = getopt(argc, argv, ":f:N:p:h")) != -1 )
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
                "Energy and geometry      [-e]           (STRING ARGUMENT REQUIRED)     File with energetic interactions and geometric bounds\n"; 
                exit(EXIT_SUCCESS);
                break;


            case 'p':
                // std::cout <<"Option p was called with argument " << optarg << std::endl;
                positions = optarg;
                break;    

            case 'e':
                energy = optarg; 
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

    // const int x{15}, y{15}, z{15}, kT{1}; 

    std::vector <double> info_vec = ExtractTopologyFromFile("energy.txt"); 

    for (auto i: info_vec){
        std::cout << i << " | ";
    }
    /*Grid G = CreateGridObject(x, y, z, kT, positions);

    G.PolymersInGrid.at(0).printChainCoords();
    auto start = std::chrono::high_resolution_clock::now(); 

    for (int i{0}; i < Nmov; i++){

        G.TheElementaryGridEvolver();
    
    }

    auto stop = std::chrono::high_resolution_clock::now(); 

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start); 

    std::cout << "Your simulation ran for " << duration.count() << " seconds." << std::endl;

*/
    return 0;

}