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
#include "classes.h"
#include "misc.h"


int main(int argc, char** argv) {

    // set up 
    int opt; 
    int dfreq {-1}, max_iter{-1};
    std::string positions {"__blank__"}, topology {"__blank__"}, dfile {"__blank__"}, efile{"__blank__"}, mfile {"__blank__"}, stats_file {"__blank__"}, \
    lattice_file_write {"__blank__"}, lattice_file_read {"__blank__"};  
    bool v = false, r = false;

    while ( (opt = getopt(argc, argv, ":s:L:R:f:M:o:u:p:t:e:vhr")) != -1 )
    {
        switch (opt) 
        {
            case 'f':
                dfreq = atoi(optarg); 
                break;

            case 'M':
                max_iter = atoi(optarg); 
                break; 

            case 'h':
                std::cout << 
                "\nWelcome to my Monte Carlo simulation engine (v0.5) for polymers in a simple cubic box (Z=26). \nThis is a simulation engine which incorporates a single type of comonomer along with a singular solvent. \nIn this implementation, I have employed reversing moves to avoid copying. Did wonders for efficiency." <<
                "\nThis system has cosolvent present. \n"
                "\nLast updated: June 15, 2022, 12:39 PM. Added restarting capabilities. \nAuthor: satyend@princeton.edu\n" <<
                "\n----------------------------------------------------------------------------------------------------------------------------------\n" << 
                "These are all the inputs the engine accepts for a single run, as of right now:\n\n" <<
                "help                      [-h]           (NO ARG REQUIRED)              Prints out this message. \n"<<
                "verbose flag              [-v]           (NO ARG REQUIRED)              Prints out a lot of information in console. MEANT FOR DEBUGGING PURPOSES. \n"<<
                "restart flag              [-r]           (NO ARG REQUIRED)              Restarts simulation from final spot of a previous simulation. \n"<<
                "Dump Frequency            [-f]           (INTEGER ARGUMENT REQUIRED)    Frequency at which coordinates should be dumped out. \n"<<                
                "Number of maximum moves   [-M]           (INTEGER ARGUMENT REQUIRED)    Number of MC moves to be run on the system. \n" <<
                "Polymer coordinates       [-p]           (STRING ARGUMENT REQUIRED)     Name of input file with coordinates of polymer.\n" <<
                "Energy and geometry       [-t]           (STRING ARGUMENT REQUIRED)     Name of input file with energetic interactions and geometric bounds.\n" <<
                // "Solvent coordinates       [-S]           (STRING ARGUMENT REQUIRED)     Name of output file with coordinates of solvent. \n" <<
                "Energy of grid            [-u]           (STRING ARGUMENT REQUIRED)     Name of output file with energy of system at each step in a file.\n"<<
                "Lattice file to write to  [-L]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
                "Lattice file to read from [-R]           (STRING ARGUMENT REQUIRED)     Trajectory file of a previous simulation which can be used to start current simulation.\n" <<
                "Orientation file          [-e]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain orientation of ALL particles in system.\n" << 
                "Move statistics file      [-s]           (STRING ARGUMENT REQUIRED)     Name of output file with move statistics. \n" <<
                "Name of output file       [-o]           (STRING ARGUMENT REQUIRED)     Name of output file which will contain coordinates of polymer.\n\n";  
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

            case 's':
                stats_file = optarg;
                break;

            case 'r':
                std::cout << "Simulation will be restarted from the end of previous simulation."; 
                r = true;
                break;

            case '?':
                std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
                exit(EXIT_FAILURE); 
                break;

            case 'v':
                std::cout << "Output to console will be verbose. " << std::endl;
                v = true;
                break;

            case 'e':
                mfile=optarg;
                break; 

            case ':':
                std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
                exit(EXIT_FAILURE);           
                break; 

            case 'L': 
                // std::cout << "Name of file to write lattice down at end of simulation." << std::endl;
                lattice_file_write=optarg; 

                break;

            case 'R':
                // std::cout << "Name of file to read lattice to restart simulation." << std::endl;
                lattice_file_read=optarg; 
                std::cout << "Name of file to write lattice down at end of simulation is " << lattice_file_read <<  "." << std::endl;
                break;
        }
    }

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // parse inputs 
    // This command will take all of the above inputs and make sure they are valid. 
    InputParser (dfreq, max_iter, r, positions, \
        topology, dfile, efile, mfile, stats_file, \
        lattice_file_read ); 

    // driver 

    // ExtractTopologyFromFile extracts all the topology from the input file 
    std::array <double,9> info_vec {ExtractTopologyFromFile(topology)}; 

    // ExtractNumberOfPolymers extracts the total number of chains in the input file 
    const int N = ExtractNumberOfPolymers(positions); 

    // assign values from info_vec to variables 
    const int x = info_vec[0];
    const int y = info_vec[1]; 
    const int z = info_vec[2]; 
    const double T = info_vec[3]; 
    
    std::array <int,3> attempts    = {0,0,0}; 
    std::array <int,3> acceptances = {0,0,0}; 
    
    std::array <double,5> E = { info_vec[4], info_vec[5], info_vec[6], info_vec[7], info_vec[8]};


    // define the interaction map 
    
    std::cout << std::endl;

    std::cout << "Preparing for take-off...\n\n" ; 
    std::cout << "Chemical information: " << std::endl;
    std::cout << "Number of polymers is " << N << ".\n\n";
    std::cout << "Geometric information about simulation cell: " << std::endl;
    std::cout << "x = " << x <<", y = " << y << ", z = "<< z << ".\n" << std::endl;
    std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
    std::cout << "Temperature = " << T << "." << std::endl; 
    std::cout << "Em1m1 = " << E[0] <<", Em2m2 = " << E[1] <<", Em1m2 = " << E[2] << ", Em1s1" << E[3] << ", Em1s2 = " << E[4] << ".\n \n";  

    std::cout << "Off to a good start. \n\n";
    std::cout << "--------------------------------------------------------------------\n" << std::endl;

    std::cout << "Running some more checks on input... \n\n" ; 

    // THIS MIGHT NEED TO CHANGE 

    int step_number   = 0; // step_number++;
    double sysEnergy  {0}; 
    double sysEnergy_ {0}; 
    
    std::vector <Particle*> LATTICE;
    LATTICE.reserve (x*y*z); 
    
    std::vector <Polymer> Polymers; 
    Polymers.reserve(N); 

    auto start = std::chrono::high_resolution_clock::now(); 
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

    
    // start creating LATTICE! 
    // if i am not restarting, get stuff from 
    // positions and set up LATTICE. 

    if ( !r ){
        std::cout << "Revving up the system from scratch..." << std::endl;
        Polymers = ExtractPolymersFromFile (positions, x, y, z); 

        std::cout << "Solvating the simulation cell... This can take some time. \n";

        AddSolvent1 (x, y, z, &LATTICE);

        std::vector <int> monomer_indices; 

        // populate the lattice 
        for (Polymer& pmer: Polymers){
            for (Particle*& p: pmer.chain ){
                
                LATTICE.at(lattice_index (p->coords, y, z) ) = p; 
                monomer_indices.push_back ( lattice_index (p->coords, y, z)); 
            
            }
        }
    
        dumpPositionsOfPolymers(&Polymers, step_number, dfile); 

    }

    else {

        std::cout << "Restarting from (some) previous simulation..." << std::endl;

        LATTICE  = ExtractLatticeFromRestart ( lattice_file_read, &step_number, x, y, z );
        Polymers = ExtractPolymersFromTraj   ( dfile, positions, step_number, x, y, z );  

        for ( Polymer& pmer: Polymers) {
            for ( Particle*& p: pmer.chain){

                std::cout << "ptypes: " << LATTICE[lattice_index(p->coords, y, z) ]->ptype << ", " << p->ptype << std::endl;
                std::cout << "orientation: " << LATTICE[lattice_index(p->coords, y, z) ]->orientation << ", " << p->orientation << std::endl;
                std::cout << "locs: "; print(LATTICE[lattice_index(p->coords, y, z) ]->coords); print( p->coords );

                if ( !( LATTICE [ lattice_index(p->coords, y, z) ]->ptype == p->ptype \
                    && LATTICE [ lattice_index(p->coords, y, z) ]->coords == p->coords \
                    && LATTICE [ lattice_index(p->coords, y, z) ]->orientation == p->orientation) ) {
                    std::cerr << "There is a problem with extraction for restart..." << std::endl;
                    exit (EXIT_FAILURE); 
                } 
                LATTICE [ lattice_index(p->coords, y, z) ] = p;
            }
        }

        stop = std::chrono::high_resolution_clock::now(); 
        duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

        std::cout << "System set-up took " << duration.count () << " milliseconds." << std::endl;
        std::cout << "Simulation cell has been made! \n\n" ;

    }
    
    Polymers[0].printChainCoords(); 

    std::array <double,5> contacts = {0, 0, 0, 0, 0};

    sysEnergy = CalculateEnergy(&Polymers, &LATTICE, x, y, z, &E, &contacts); 

    std::cout <<"\nCalculating energy..." << std::endl;
    std::cout << "Energy of system is " << sysEnergy << ".\n" << std::endl;

    std::array <double,5> contacts_copy = contacts; 

    if ( !r ) {
        dumpEnergy      (sysEnergy, step_number, &contacts, efile); 
        dumpOrientation (&Polymers, &LATTICE, step_number, mfile, x, y, z); 
    }
    if (r) {
        dumpOrientation (&Polymers, &LATTICE, step_number, mfile, x, y, z); 
    }
    
    std::cout <<"\nCalculating energy..." << std::endl;
    std::cout << "Energy of system is " << sysEnergy << ".\n" << std::endl;
    
    bool IMP_BOOL = true; 
    bool metropolis = false;
    
    double rweight =  0; 
    int move_number = 0; 
    int monomer_index = -1; 
    int back_or_front = -1; 

    std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>> memory3; 
    // std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>> memory2;
    // std::vector <std::array<int,2>> s_memory; 

    printf("Initiation complete. We are ready to go. The engine will output information every %d configuration(s).\n", dfreq); 
    
    int Nsurr = contacts[2] + contacts[3] + contacts[4] + contacts[5]; 
    std::cout << "Number of surrounding solvent molecules is " << Nsurr << "." << std::endl;
    
    for (int i = step_number+1; i< (step_number+max_iter+1); i++) {

        if ( v && (i%dfreq==0) ){
            printf("Step number: %d.\n", i);
        }

        if ( !(metropolis) ){
            
            contacts = contacts_copy; 
            Nsurr = contacts[2]+contacts[3]+contacts[4]+contacts[5];
        }
        else {
            contacts_copy = contacts; 
            Nsurr = contacts[2]+contacts[3]+contacts[4]+contacts[5];
            metropolis = false; 
        }

        // choose a move... 

        PerturbSystem (&Polymers, &LATTICE, x, y, z, v, &IMP_BOOL, &rweight, &attempts, &move_number, &memory3, &monomer_index, &back_or_front); 


        if (IMP_BOOL){ 
            sysEnergy_ = CalculateEnergy (&Polymers, &LATTICE, x, y, z, &E, &contacts); 
        }

        if ( v && (i%dfreq==0) ){
            printf("Executing...\n");
        }
        
        if ( IMP_BOOL ) {
            if ( MetropolisAcceptance (sysEnergy, sysEnergy_, T, rweight) ){
                metropolis = true; 
                acceptances[move_number]+=1;

                if ( v ){
                    printf("Checking validity of coords...");
                    printf("checkForOverlaps says: %d.\n", checkForOverlaps(Polymers)); 
                    if (!checkForOverlaps(Polymers)){
                        printf("Something is fucked up overlaps-wise. \n");
                        exit(EXIT_FAILURE);
                    }

                    if (! checkForOverlaps (&Polymers, &LATTICE)){
                        printf("Random monomer floating!!");
                        exit (EXIT_FAILURE);
                    }
                    std::cout << "No random monomers floating around!!" << std::endl;

                    if (!checkForSolventMonomerOverlap (&Polymers, &LATTICE, y, z) ){
                        printf("Something is fucked up solvent-monomer overlaps-wise. \n");
                        exit(EXIT_FAILURE);
                    }

                    printf("checkConnectivity says: %d\n", checkConnectivity(Polymers, x, y, z)); 
                    if (!checkConnectivity(Polymers, x, y, z) ){
                        printf("Something is fucked up connectivity-wise. \n");
                        exit(EXIT_FAILURE);
                    }
                    printf("Accepted!!\n");
                    printf("Energy of the system is %.2f.\n", sysEnergy_);
                    printf("This should be 1 as IMP_BOOL must be true on acceptance: %d\n", IMP_BOOL);

                    if (checkPointersOnLattice (&LATTICE, x, y, z) ){
                        printf("We good. LATTICE is in good shape. \n\n\n"); 
                    }
                    else {
                        std::cerr <<"Something is fucked with pointers on LATTICE." << std::endl; 
                        exit (EXIT_FAILURE);
                    }

                }

            sysEnergy = sysEnergy_; 

            }

            else {

                ReversePerturbation (&Polymers, &LATTICE, y, z, v, move_number, &memory3, monomer_index, back_or_front);

                // std::cout << "Polymer coordinates after reversal are: " << std::endl;
                // Polymers[0].printChainCoords();

                // std::cout << "After reversing perturbation..." << std::endl;

                if ( v ){
                    std::cout << "Rejected. Reversal in progress..." << std::endl;
                    std::cout << "move_number is " << move_number << "." << std::endl;
                    printf("Checking validity of coords...");
                    printf("checkForOverlaps says: %d.\n", checkForOverlaps(Polymers)); 
                    if (!checkForOverlaps(Polymers)){
                        printf("Something is fucked up overlaps-wise. \n");
                        exit(EXIT_FAILURE);
                    }

                    if (! checkForOverlaps (&Polymers, &LATTICE)){
                        printf("Random monomer floating!!");
                        exit (EXIT_FAILURE);
                    }

                    std::cout << "No random monomers floating around!!" << std::endl;

                    if (!checkForSolventMonomerOverlap (&Polymers, &LATTICE, y, z) ){
                        printf("Something is fucked up solvent-monomer overlaps-wise. \n");
                        exit(EXIT_FAILURE);
                    }

                    printf("checkConnectivity says: %d\n", checkConnectivity(Polymers, x, y, z)); 
                    if (!checkConnectivity(Polymers, x, y, z) ){
                        printf("Something is fucked up connectivity-wise. \n");
                        exit(EXIT_FAILURE);
                    }
                    printf("Reversed successfully!!\n");
                    printf("Energy of the system is %.2f.\n", sysEnergy_);
                    printf("This should be 1 as IMP_BOOL must be true on acceptance: %d\n", IMP_BOOL);

                    if (checkPointersOnLattice (&LATTICE, x, y, z) ){
                        printf("We good. LATTICE is in good shape. \n"); 
                    }
                    else {
                        std::cerr <<"Something is fucked with pointers on LATTICE." << std::endl; 
                        exit (EXIT_FAILURE);
                    }

                }


                if ( v && (i%dfreq==0) ){
                    printf ("Not accepted.\n");
                    printf ("Energy of the suggested system is %.2f, while energy of the initial system is %.2f.\n\n\n", sysEnergy_, sysEnergy);
                }
                if (v && !IMP_BOOL){
                    std::cout << "There was no change in the state of the system." << std::endl;
                }
                
            }

        }

        else {
            if (v){
                printf("IMP_BOOL is zero. Nothing will be done.\n\n\n");

                if (!checkForOverlaps(Polymers)){
                        printf("Something is fucked up overlaps-wise. \n");
                        exit(EXIT_FAILURE);
                    }

                if (! checkForOverlaps (&Polymers, &LATTICE)){
                    printf("Random monomer floating!!");
                    exit (EXIT_FAILURE);
                }

                std::cout << "No random monomers floating around!!" << std::endl;

                if (!checkForSolventMonomerOverlap (&Polymers, &LATTICE, y, z) ){
                    printf("Something is fucked up solvent-monomer overlaps-wise. \n");
                    exit(EXIT_FAILURE);
                }

                printf("checkConnectivity says: %d\n", checkConnectivity(Polymers, x, y, z)); 
                if (!checkConnectivity(Polymers, x, y, z) ){
                    printf("Something is fucked up connectivity-wise. \n");
                    exit(EXIT_FAILURE);
                }
                printf("Reversed successfully!!\n");
                printf("Energy of the system is %.2f.\n", sysEnergy_);
                printf("This should be 1 as IMP_BOOL must be true on acceptance: %d\n", IMP_BOOL);

                if (checkPointersOnLattice (&LATTICE, x, y, z) ){
                    printf("We good. LATTICE is in good shape. \n"); 
                }
                else {
                    std::cerr <<"Something is fucked with pointers on LATTICE." << std::endl; 
                    exit (EXIT_FAILURE);
                }
            }
        }

        if ( ( i % dfreq == 0) ){
           
            dumpPositionsOfPolymers (&Polymers, i, dfile); 
            
            if ( metropolis ){
                dumpEnergy (sysEnergy, i, &contacts, efile);
            }
            else {
                dumpEnergy (sysEnergy, i, &contacts_copy, efile);
            } 
        }

        // reset the memory carrier, and IMP_BOOL
        reset (memory3);
        // reset (memory2);
        // reset (s_memory); 
        IMP_BOOL = true; 
           
    }
    dumpMoveStatistics      (&attempts, &acceptances, max_iter, stats_file);  
    
    if ( lattice_file_write != "__blank__" ){
        dumpLATTICE ( &LATTICE, step_number+max_iter, y, z, lattice_file_write ); 
    }

    stop = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 
	
    std::cout << "\n\nTime taken for simulation: " << duration.count() << " milliseconds.\n"; 
    std::cout << "That is all she wrote. Hope it worked." << std::endl;
    std::cout << "--------------------------------------------------------------------\n\n";

    return 0;

}
