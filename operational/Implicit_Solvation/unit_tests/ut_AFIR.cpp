#include <iostream> 
#include <fstream>
#include <vector> 
#include <string> 
#include <map>
#include <random>
#include <chrono>
#include <getopt.h> 
#include <stdexcept>
#include <stdlib.h> 
#include <array>
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

    /*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~
    ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~*/

    Grid G = CreateGridObject(positions, topology);

    G.instantiateOccupancyMap();
    
    G.CalculateEnergy();
    std::cout << "Energy of box is: " << G.Energy << std::endl;

    bool b = true; 
    
    Grid G3 = FinalIndexRotationAgg(&G, 0, &b) ; 


    std::cout << "On the surface..." << std::endl; 
    std::cout << "\n\nCoordinates of monomer unit of Polymer 0 in G are: " << std::endl;

    G.PolymersInGrid.at(0).printChainCoords();

    std::cout << "\n\nCoordinates of monomer unit of Polymer 0 in G3 are: " << std::endl;
    G3.PolymersInGrid.at(0).printChainCoords();    

    std::cout << "\n\n" << std::endl;

    std::cout << "\n\nCoordinates of monomer unit of Polymer 1 in G are: " << std::endl;

    G.PolymersInGrid.at(1).printChainCoords();

    std::cout << "\n\nCoordinates of monomer unit of Polymer 1 in G3 are: " << std::endl;
    G3.PolymersInGrid.at(1).printChainCoords();

    std::cout << "\n\n*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~\n\n" << std::endl;

    std::cout << "This is my unit test." << std::endl;
    std::cout << "Check 1: Check number of elements in PolymersInGrid in both G and G3." << std::endl; 
    std::cout << "number of elements in PolymersInGrid in G is " << G.PolymersInGrid.size() << " and in G3 is " << G3.PolymersInGrid.size() << "." << std::endl;
    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~" << std::endl; 

    if (G.PolymersInGrid.size() != G3.PolymersInGrid.size() ){
        std::cerr << "Error in # of PolymersInGrid." << std::endl;
        exit(EXIT_FAILURE);
    }

    for (size_t i=0; i<G.PolymersInGrid.size(); i++){
        std::cout << "Check: Check number of elements in Polymers " << i << " in both G and G3." << std::endl; 
        std::cout << "Number of elements in Polymer " << i << " in G is " << G.PolymersInGrid.at(i).chain.size() 
        << " and in G3 is " << G3.PolymersInGrid.at(i).chain.size() << "." << std::endl; 
        std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~" << std::endl; 
        if (G.PolymersInGrid.at(i).chain.size() != G3.PolymersInGrid.at(i).chain.size() ){
        std::cerr << "Error in # of monomers in Polymer "<< i << "." << std::endl;
        exit(EXIT_FAILURE);
        }
    }


    std::cout << "Find discrepancy in PolymersInGrid." << std::endl; 
    for (size_t i=0; i<G.PolymersInGrid.size(); i++){
        for (size_t j=0; j<G.PolymersInGrid.at(i).chain.size(); j++){
            if (!(G.PolymersInGrid.at(i).chain.at(j) == G3.PolymersInGrid.at(i).chain.at(j)) ){
                std::cout << "---------------FOUND THE DISPLACED PARTICLE---------------" << std::endl;
                std::cout << "The orientation of the particle originally is " << G.PolymersInGrid.at(i).chain.at(j).orientation << std::endl;
                std::cout << "The orientation of the after is " << G3.PolymersInGrid.at(i).chain.at(j).orientation << std::endl;
                std::cout << "The position of the particle originally is ";
                print( G.PolymersInGrid.at(i).chain.at(j).coords ) ; 
                std::cout << "The position of the particle after is ";
                print( G3.PolymersInGrid.at(i).chain.at(j).coords ) ; 
                std::cout << "The type of the particle originally is " << G.PolymersInGrid.at(i).chain.at(j).ptype << std::endl;
                std::cout << "The type of the particle after is " << G3.PolymersInGrid.at(i).chain.at(j).ptype << std::endl;
                
                if (!(G.PolymersInGrid.at(i).chain.at(j).orientation == G3.PolymersInGrid.at(i).chain.at(j).orientation) ){
                    std::cerr << "There is a problem in orientation." << std::endl;
                    exit (EXIT_FAILURE);
                }
                else if (!(G.PolymersInGrid.at(i).chain.at(j).ptype == G3.PolymersInGrid.at(i).chain.at(j).ptype)) {
                    std::cerr << "There is a problem in type." << std::endl;
                    exit (EXIT_FAILURE);   
                }

            }
        }
    }
    std::cout << "\n\n"; 

    std::cout << "In the original case, G..." << std::endl;
    for (size_t i=0; i<G.PolymersInGrid.size(); i++){
        std::cout << "For Polymer " << i << ":" << std::endl;
        for (size_t j=0; j<G.PolymersInGrid.at(i).chain.size(); j++){

            // print out the particle and its connections 
            std::vector <Particle> pvec; 
            pvec = G.PolymersInGrid.at(i).ConnectivityMap[G.PolymersInGrid.at(i).chain.at(j)]; 
            std::cout << "The original particle is: ";
            print(G.PolymersInGrid.at(i).chain.at(j).coords ) ;
            std::cout << "The connections are: \n"; 
            for (auto v: pvec){
                print(v.coords); 
            } 

        }
    }
    std::cout << "\n\nIn the MOVED case, G3..." << std::endl;
    for (size_t i=0; i<G3.PolymersInGrid.size(); i++){
        std::cout << "\n\nFor Polymer " << i << ":" << std::endl;
        for (size_t j=0; j<G3.PolymersInGrid.at(i).chain.size(); j++){

            // print out the particle and its connections 
            std::vector <Particle> pvec; 
            pvec = G3.PolymersInGrid.at(i).ConnectivityMap[G3.PolymersInGrid.at(i).chain.at(j)]; 
            std::cout << "The original particle is: ";
            print(G3.PolymersInGrid.at(i).chain.at(j).coords ) ;
            std::cout << "The connections are: \n"; 
            for (auto v: pvec){
                print(v.coords); 
            } 

        }
    }    

    

    std::cout<<"\n\ntime for the main event: OccupancyMap." << std::endl << std::endl;
    int c {0}; 
    for (auto it = G3.OccupancyMap.begin(); it != G3.OccupancyMap.end(); it++){
        print(it->first);
        c++;
    }
    std::cout << "number of keys in G3 is " << c << std::endl << std::endl;

    c = 0;
    for (auto it = G.OccupancyMap.begin(); it != G.OccupancyMap.end(); it++){
        print(it->first);
        c++;
    }
    
    std::cout << "number of key in G is " << c << std::endl << std::endl;

    for (auto it = G.OccupancyMap.begin(), it3 = G3.OccupancyMap.begin(); it != G.OccupancyMap.end(); it++, it3++){

        std::cout << "Key:   "; 
        for (auto k: it->first){
            std::cout << k << " | ";
        }

        std::cout << ", Key3:   "; 
        for (auto k: it3->first){
            std::cout << k << " | ";
        }

        std::cout << "\nValue: "; 
        for (auto k: it->second.coords){
            std::cout << k << " | ";
        }

        std::cout << ", Value3: "; 
        for (auto k: it3->second.coords){
            std::cout << k << " | ";
        }
        std::cout << "\n";
    }
    

    return 0;


}
