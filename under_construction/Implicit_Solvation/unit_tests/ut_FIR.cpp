#include <iostream> 
#include <fstream>
#include <vector> 
#include <string> 
#include <map>
#include <random>
#include <chrono>
#include <array>
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

    /*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~
    ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~*/

    Grid G = CreateGridObject(positions, topology);

    G.instantiateOccupancyMap();
    
    G.CalculateEnergy();
    std::cout << "Energy of box is: " << G.Energy << std::endl;
    
    bool b = true; 

    Grid G3 = FinalIndexRotation(&G, 0, &b) ; 


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
                    exit(EXIT_FAILURE);
                }
                else if (!(G.PolymersInGrid.at(i).chain.at(j).ptype == G3.PolymersInGrid.at(i).chain.at(j).ptype)) {
                    std::cerr << "There is a problem in type." << std::endl;
                    exit(EXIT_FAILURE);   
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

    for (std::map <std::array <int,3>, Particle>::iterator iter=G.OccupancyMap.begin(); iter!=G.OccupancyMap.end(); ++iter){

        std::array <int,3> key = iter->first; 

        std::cout << "key is "; 
        print(key); 
        std::vector <int> v = {1, 0, 0};
        if (key.at(0) == v.at(0) && key.at(1) == v.at(1) && key.at(2) == v.at(2)){
            std::cout << "problem key is "; 
            print(v); 
            continue; 
        }

        if ( (G.OccupancyMap.at(key) == G3.OccupancyMap.at(key) ) ){
            
            std::cout << "The coordinates of the particle originally were: ";
            print(G.OccupancyMap.at(key).coords); 
            std::cout << "The coordinates of the particle after is: ";
            print(G3.OccupancyMap.at(key).coords);
            std::cout << "The orientation of the particle was " << G.OccupancyMap.at(key).orientation << std::endl;
            std::cout << "The orientation of the particle is " << G3.OccupancyMap.at(key).orientation << std::endl; 
            std::cout << "The particle type was " << G.OccupancyMap.at(key).ptype << std::endl;
            std::cout << "The particle type was " << G3.OccupancyMap.at(key).ptype << std::endl;

            if (!(G.OccupancyMap.at(key).orientation == G3.OccupancyMap.at(key).orientation) ){
                std::cerr << "There is a problem in orientation." << std::endl;
                exit(EXIT_FAILURE);
            }
            else if ( !(G.OccupancyMap.at(key).ptype == G3.OccupancyMap.at(key).ptype) ) {
                std::cerr << "There is a problem in type." << std::endl;
                exit(EXIT_FAILURE);   
            }
            else {
                continue; 
            }

        }
            
    }

    std::cout << "\n\n";
    
   

    return 0;


}
