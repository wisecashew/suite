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
    int dfreq {-1}, max_iter{-1};
    std::string positions {"blank"}, topology {"blank"}, dfile {"blank"}, efile{"blank"}, mfile {"blank"}, stats_file {"blank"};  
    bool v = false;

    while ( (opt = getopt(argc, argv, ":s:f:M:o:u:p:t:e:vh")) != -1 )
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
                "This is the main driver for the monte carlo simulation of polymers in a box. This is an explicit solvent simulation.\n" <<
		        "This version number 0.2.3 of the Monte Carlo Engine. \nThis will keep the energydump with the number of monomer-monomer contacts, aligned and misaligned interactions. \nSet up on Mar 30, 2022, 01:00 AM.\n" <<
                "This code employs Rosenbluth sampling and both local and multiple solvent flips.\n" << 
                "These are all the options we have available right now: \n" <<
                "help                     [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose                  [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
                "Dump Frequency           [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
                "Number of maximum moves  [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
                // "Required accepted moves  [-N]           (INTEGER ARGUMENT REQUIRED)    Number of accepted moves for a good simulation.\n" <<  
                "Position coordinates     [-p]           (STRING ARGUMENT REQUIRED)     File with position coordinates.\n" <<
                "Energy of grid           [-u]           (STRING ARGUMENT REQUIRED)     Dump energy of grid at each step in a file.\n"<<
                "Energy and geometry      [-t]           (STRING ARGUMENT REQUIRED)     File with energetic interactions and geometric bounds ie the topology.\n" <<
                // "Previous trajectory file [-T]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
                "Orientation file         [-e]           (STRING ARGUMENT REQUIRED)     Name of file which will contain orientation of monomer and neighboring solvent particles.\n" << 
                "Move statistics file     [-s]           (STRING ARGUMENT REQUIRED)     Name of file with move statistics. \n" <<
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

            //case 'T':
            //    restart_traj = optarg;
            //   break;

            case 'u':
                efile = optarg;
                break;

            case 's':
                stats_file = optarg;
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
    InputParser (dfreq, max_iter, positions, topology, dfile, efile, mfile, stats_file); 

    // driver 

    auto start = std::chrono::high_resolution_clock::now(); 


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
    
    std::array <int,9> attempts    = {0,0,0,0,0,0,0,0,0};
    std::array <int,9> acceptances = {0,0,0,0,0,0,0,0,0}; 

    std::cout << "Number of polymers is " << N << ".\n";
    std::cout << "x = " << x <<", y = " << y << ", z = "<< z << ", T = " << T << ".\n"; 
    std::cout << "Emm_a = " << Emm_a <<", Emm_n = " << Emm_n << ", Ems_a = "<< Ems_a << ", Ems_n = " << Ems_n <<".\n";  

    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~" << std::endl;
    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~\n" << std::endl;


    // THIS MIGHT NEED TO CHANGE 
    

    int step_number = 0;
    double sysEnergy {0}; 
    std::vector <Polymer> Polymers; 
    
    Polymers.reserve(N); 

    Polymers = ExtractPolymersFromFile(positions, x, y, z); 
        
    double sysEnergy_ {0};

    std::vector <Particle> Solvent = CreateSolventVector(x, y, z, &Polymers);
    
    /////////////////////////////////////////////////
    
    for (Polymer& pmer:Polymers){
        for (Particle& p: pmer.chain){
            p.orientation = 0; 
        }
    }

    /////////////////////////////////////////////////

    dumpPositionsOfPolymers(&Polymers, step_number, dfile); 
    
    double mm_aligned  = 0,  mm_aligned_copy  = 0; 
    double mm_naligned = 0,  mm_naligned_copy = 0;
    int ms_aligned     = 0,  ms_aligned_copy  = 0;
    int ms_naligned    = 0,  ms_naligned_copy = 0; 
    
    sysEnergy = CalculateEnergy(&Polymers, &Solvent, x, y, z, Emm_a, Emm_n, Ems_a, Ems_n, &mm_aligned, &mm_naligned, &ms_aligned, &ms_naligned); 

    mm_aligned_copy   = mm_aligned ; 
    mm_naligned_copy  = mm_naligned; 
    ms_aligned_copy   = ms_aligned ;
    ms_naligned_copy  = ms_naligned;

    dumpEnergy      (sysEnergy, step_number, mm_aligned, mm_naligned, ms_aligned, ms_naligned, efile); 
    dumpOrientation (&Polymers, &Solvent, step_number, mfile, x, y, z); 
    
    // defined single orientation solvents and polymers 

    std::cout << "Energy of system is " << sysEnergy << std::endl;
     
	int acc_counter = 0; 
    
    bool IMP_BOOL = true; 
    bool metropolis = false;
 
    std::vector <Polymer> Polymers_c ; 
    std::vector <Particle> Solvent_c; 
    
    int Nsurr = ms_aligned+ms_naligned; 
    double rweight =  0; 
    int move_number = 0; 

    printf("Simulation will output information of every %d configuration.\n", dfreq); 
    
    for (int i = step_number+1; i< (step_number+max_iter+1); i++) {

        if ( v && (i%dfreq==0) ){
            printf("Move number %d.\n", i);
        }

        if ( !(metropolis) ){
            // m_neighbors = m_neicopy; a_contacts = a_contcopy; n_contacts = n_contcopy;
            mm_aligned = mm_aligned_copy; mm_naligned = mm_naligned_copy; ms_aligned = ms_aligned_copy; ms_naligned = ms_naligned_copy;
            Nsurr = ms_aligned + ms_naligned;
        }
        else {
            // m_neicopy = m_neighbors; a_contcopy = a_contacts; n_contcopy = n_contacts; 
            mm_aligned_copy = mm_aligned; mm_naligned_copy = mm_naligned; ms_aligned_copy = ms_aligned; ms_naligned_copy = ms_naligned;
            Nsurr = ms_aligned + ms_naligned;
            metropolis = false; 
        }
        
        // make a copy of the polymer+solvation shell 
        Solvent_c  = Solvent; 
        Polymers_c = Polymers;

        // choose a move... 
        // i think movechooser_rosenbluth gotta be void... 
        MoveChooser_Rosenbluth (&Polymers_c, &Solvent_c, x, y, z, v, &IMP_BOOL, &rweight, &attempts, &move_number, Nsurr); 
        
        if (IMP_BOOL){ 
            sysEnergy_ = CalculateEnergy (&Polymers_c, &Solvent_c, x, y, z, Emm_a, Emm_n, Ems_a, Ems_n, &mm_aligned, &mm_naligned, &ms_aligned, &ms_naligned); 
        }

        if ( v && (i%dfreq==0) ){
            printf("Executing...\n");
        }
        
        if ( IMP_BOOL && MetropolisAcceptance (sysEnergy, sysEnergy_, T, rweight) ) {
            metropolis = true; 
            acceptances[move_number-1]+=1;
            // std::cout << "bro..." << std::endl;
            // replace old config with new config
            if ( v ){
                printf("Checking validity of coords...");
                printf("checkForOverlaps says: %d.\n", checkForOverlaps(Polymers_c)); 
                std::cout <<"Is this the last thing?" << std::endl;
                if (!checkForOverlaps(Polymer_c)){
                    printf("Something is fucked up overlaps-wise. \n");
                    exit(EXIT_FAILURE);
                }

                if (!checkForSolventMonomerOverlap (&Polymers_c, &Solvent_c) ){
                    printf("Something is fucked up solvent-monomer overlaps-wise. \n");
                    exit(EXIT_FAILURE);
                }

                printf("checkConnectivity says: %d\n", checkConnectivity(Polymer_c, x, y, z)); 
                std::cout << "Is this another last thing?" << std::endl;
                if (! checkConnectivity(Polymer_c, x, y, z) ){
                    printf("Something is fucked up connectivity-wise. \n");
                    exit(EXIT_FAILURE);
                }
                printf("Accepted.\n");
                printf("Energy of the system is %.2f.\n", sysEnergy_);
                printf("This should be 1 as IMP_BOOL must be true on acceptance: %d\n", IMP_BOOL);


            }

            // std::cout << "IMP_BOOL is " << IMP_BOOL << std::endl;
            
            Polymers = std::move(PolymerVector_c);
            Solvent = std::move(SolventVector_c); 
            sysEnergy = sysEnergy_; 
        }


        else {
            if ( v && (i%dfreq==0) ){
                printf("Not accepted.\n");
                printf("Energy of the suggested system is %.2f, while energy of the initial system is %.2f.\n", sysEnergy_, sysEnergy);
            }
            if (v && !IMP_BOOL){
                std::cout << "There was no change in the state of the system." << std::endl;
            }
            
        }
        
        if ( ( i % dfreq == 0) ){
           
            dumpMoveStatistics     (&attempts, &acceptances, i, stats_file);  
            dumpPositionsOfPolymers(&Polymer, i, dfile); 
            dumpOrientation ( &Polymer, &Solvent, i, mfile, x, y, z); 
            if ( metropolis ){
                dumpEnergy (sysEnergy, i, mm_aligned, mm_naligned, ms_aligned, ms_naligned, efile);
            }
            else {
                dumpEnergy (sysEnergy, i, mm_aligned_copy, mm_naligned_copy, ms_aligned_copy, ms_naligned_copy, efile);
            } 
        }
        
        IMP_BOOL = true; 
           
    }
    
    dumpPositionOfSolvent(&Solvent, max_iter, "solvent_coords");
    
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 
	
    std::cout << "\n\nTime taken for simulation: " << duration.count() << " milliseconds.\n"; 

    return 0;

}
