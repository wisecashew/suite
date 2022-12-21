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
#include <set>
#include "classes.h"
#include "misc.h"


int main(int argc, char** argv) {

    // set up 
    int opt; 
    int dfreq {-1}, max_iter{-1};
    std::string positions {"blank"}, topology {"blank"}, dfile {"blank"}, efile{"blank"}, mfile {"blank"}, stats_file {"blank"}, solvent_file {"blank"};  
    bool v = false;

    while ( (opt = getopt(argc, argv, ":s:S:f:M:o:u:p:t:e:vh")) != -1 )
    {
        switch (opt) 
        {
            case 'f':
                // std::cout << "Option dreq was called with argument " << optarg << std::endl; 
                dfreq = atoi(optarg); 
                break;



            // case 'N':
                // std::cout << "Option Nmov was called with argument " << optarg << std::endl; 
                // Nacc = atoi(optarg); 
                // break; 

            case 'M':
                max_iter = atoi(optarg); 
                break; 

            case 'h':
                std::cout << 
                "\nWelcome to my Monte Carlo simulation engine (v0.1) for polymers in a simple cubic box (Z=6). \nThis is a quasi-implicit solvent simulation engine which incorporates directional bonding between monomer and solvent." <<
		        "\nLast updated: May 10, 2022, 07:42 PM. \nAuthor: satyend@princeton.edu\n" <<
                "\n----------------------------------------------------------------------------------------------------------------------------------\n" << 
                "These are all the inputs the engine accepts for a single run, as of right now:\n\n" <<
                "help                     [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose                  [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
                "Dump Frequency           [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
                "Number of maximum moves  [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
                // "Required accepted moves  [-N]           (INTEGER ARGUMENT REQUIRED)    Number of accepted moves for a good simulation.\n" <<  
                "Polymer coordinates      [-p]           (STRING ARGUMENT REQUIRED)     Name of input file with coordinates of polymer.\n" <<
                "Energy and geometry      [-t]           (STRING ARGUMENT REQUIRED)     Name of input file with energetic interactions and geometric bounds.\n" <<
                "Solvent coordinates      [-S]           (STRING ARGUMENT REQUIRED)     Name of output file with coordinates of solvent. \n" <<
                "Energy of grid           [-u]           (STRING ARGUMENT REQUIRED)     Name of output file with energy of system at each step in a file.\n"<<
                // "Previous trajectory file [-T]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
                "Orientation file         [-e]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain orientation of ALL particles in system.\n" << 
                "Move statistics file     [-s]           (STRING ARGUMENT REQUIRED)     Name of output file with move statistics. \n" <<
                "Name of output file      [-o]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain coordinates of polymer.\n\n";  
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

            //case 'T':
            //    restart_traj = optarg;
            //   break;

            case 'u':
                efile = optarg;
                break;

            case 's':
                stats_file = optarg;
                break;

            case 'S': 
                solvent_file = optarg; 
                break;

            case '?':
                std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
                exit(EXIT_FAILURE); 
                break;

            case 'v':
                std::cout << "Output to console will be verbose. " << std::endl;
                v = true;
                break;

            // case 'a':
            //    std::cout <<"Only accepted structures will be outputted." << std::endl;
            //    a = true; 
            //    break; 

            // case 'r':
            //    std::cout <<"Will attempt to restart simulation by taking coordinates from a previous trajectory file." << std::endl;
            //    r = true; 
            //    break; 
            
            case 'e':
                mfile=optarg;
                break; 

            case ':':
                std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
                exit(EXIT_FAILURE);           
                break; 
        }
    }

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // parse inputs 
    // This command will take all of the above inputs and make sure they are valid. 
    // InputParser (dfreq, max_iter, solvent_file, positions, topology, dfile, efile, mfile, stats_file); 

    // driver 

    

    /*
    // ExtractTopologyFromFile extracts all the topology from the input file 
    std::array <double,8> info_vec {ExtractTopologyFromFile(topology)}; 

    // ExtractNumberOfPolymers extracts the total number of chains in the input file 
    const int N = ExtractNumberOfPolymers(positions); 

    // assign values from info_vec to variables 
    const int x = info_vec[0];
    const int y = info_vec[1]; 
    const int z = info_vec[2]; 
    const double T = info_vec[3]; 
    const double Emm_a = info_vec[4]; 
    const double Emm_n = info_vec[5]; 
    const double Ems_a = info_vec[6];
    const double Ems_n = info_vec[7]; 
    
    // std::array <int,9> attempts    = {0,0,0,0,0,0,0,0,0};
    // std::array <int,9> acceptances = {0,0,0,0,0,0,0,0,0}; 

    std::cout << std::endl;
    std::cout << "Preparing for take-off...\n\n" ; 
    std::cout << "Chemical information: " << std::endl;
    std::cout << "Number of polymers is " << N << ".\n\n";
    std::cout << "Geometric information about simulation cell: " << std::endl;
    std::cout << "x = " << x <<", y = " << y << ", z = "<< z << ".\n" << std::endl;
    std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
    std::cout << "Temperature = " << T << "." << std::endl; 
    std::cout << "Emm_a = " << Emm_a <<", Emm_n = " << Emm_n << ", Ems_a = "<< Ems_a << ", Ems_n = " << Ems_n <<".\n \n";  

    std::cout << "Off to a good start. \n\n";
    std::cout << "--------------------------------------------------------------------\n" << std::endl;

    std::cout << "Running some more checks on input... \n\n" ; 
    */
    // THIS MIGHT NEED TO CHANGE 
    

    // int step_number = 0;
    // double sysEnergy {0}; 
    // std::vector <Polymer> Polymers; 
    
    // Polymers.reserve(N); 

    // Polymers = ExtractPolymersFromFile(positions, x, y, z); 
        
    // double sysEnergy_ {0};

    // std::vector <Particle> Solvent = CreateSolventVector(x, y, z, &Polymers);
    
    /////////////////////////////////////////////////
    
    // for (Polymer& pmer:Polymers){
    //    for (Particle& p: pmer.chain){
    //        p.orientation = 0; 
    //    }
    // }

    /////////////////////////////////////////////////

    // dumpPositionsOfPolymers(&Polymers, step_number, dfile); 

    // Polymer Pol (8, Polymers[0].chain);
    // Particle part (Polymers[0].chain[0]->coords, Polymers[0].chain[0]->ptype, Polymers[0].chain[0]->orientation); 
    // std::tuple <short,char,std::array<int,3>> particle_c (0, 's', {0,0,0});

    std::cout << "Printing out sizes..." << std::endl; 

    std::cout << "Size of char is " << sizeof (char) << std::endl;
    std::cout << "Size of short is " << sizeof (short int) << std::endl;
    std::cout << "Size of int is " << sizeof(int) << std::endl;
    std::cout << "Size of double is " << sizeof(double) << std::endl;
    std::cout << "Size of array is " << sizeof (std::array <int,3>) << std::endl;
    std::cout << "Size of vector is " << sizeof (std::vector <int>) << std::endl;
    std::cout << "Size of tuple <short,char,std::array<int,3>> is " << sizeof (std::tuple <short,char,std::array<int,3>>) << std::endl;;
    std::cout << "Size of std::string is " << sizeof(std::string) << std::endl;
    std::cout << "Size of class Polymer is " << sizeof(Polymer) << std::endl; 
    std::cout << "Size of class Particle is " << sizeof(Particle) << std::endl; 
    std::cout << "Size of std::vector<Particle*> is " << sizeof(std::vector <Particle*>) << std::endl;
    std::cout << "Size of std::vector<std::tuple> is " << sizeof(std::vector <std::tuple <short,char,std::array<int,3>>>) << std::endl;
    std::cout << "Size of std::tuple <std::array<int,3>, Particle*> is " << sizeof( std::tuple <std::array<int,3>, Particle*> ) << std::endl;
    std::cout << "Size of std::set <std::tuple <std::array<int,3>, Particle*>> is " << sizeof(std::set <std::tuple <std::array<int,3>, Particle*>>) << std::endl;
    std::cout << "Size of std::map <std::array<int,3, Particle*> is " << sizeof(std::map <std::array<int,3>, Particle*>) << std::endl;

    // std::cout << "Size of instance of tuple is " << sizeof(particle_c) << std::endl;
    // std::cout << "Size of instance of class Polymer is " << sizeof(Pol) << std::endl; 
    // std::cout << "Size of instance of class* Polymer is " << sizeof(&Pol) << std::endl; 
    // std::cout << "Size of instance of class Particle is " << sizeof(part) << std::endl;
    // std::cout << "Size of instance of class* Particle is " << sizeof(&part) << std::endl;
    
    std::cout << "Creating the solvated box. This takes time..." << std::endl;

    std::vector <std::vector <std::array<int,3>>> master; 
    std::vector <std::array<int,3>> v1 = {{1,0,0},{3,0,0},{2,0,0},{6,0,0},{8,0,0},{7,0,0},{5,0,0}}; 
    std::vector <std::array<int,3>> v2 = {{2,0,0},{4,0,0},{5,0,0},{7,0,0},{9,0,0},{11,0,0},{12,0,0}}; 
    std::vector <std::array<int,3>> link; 

    create_linked_list (v1, v2, link, &master, 1); 

    for (auto elem: master){

        for ( auto e2: elem){
            for ( auto e3: e2 ){
                std::cout << e3 << " ";
            }
            std::cout << " -> ";
        }
        std::cout <<std::endl;

    }



    std::cout << "--------------------------------------------------------------------\n\n";

    return 0;

}
