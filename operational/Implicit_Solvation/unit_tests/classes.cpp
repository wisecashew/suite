#include <iostream> 
#include <vector>
#include <string> 
#include <iterator>
#include <map>
#include <algorithm> 
#include <regex>
#include <fstream>
#include <sstream>
#include <regex>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <random>
#include "classes.h"
#include "misc.h" 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                 METHODS FOR CLASS GRID 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//============================================================
//============================================================
//
// NAME OF FUNCTION: dumpPositionsOfPolymers 
//
// PARAMETERS: (int step, std::string filename), and some attributes present in the Grid Object 
// 'step' is the current time step we are at. This is an integer which is likely defined in the driver code.  
// 'filename' is the file to which I am going to print out coordinate information about the polymer. 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers in the current Grid object \, and prints them 
// out to a text file.  
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 


void Grid::dumpPositionsOfPolymers (int step, std::string filename){
    std::ofstream dump_file(filename, std::ios::app); 
    dump_file <<"Dumping coordinates at step " << step << ".\n";
    int count = 0; 
    for (Polymer pmer: this->PolymersInGrid){
        
        dump_file <<"Dumping coordinates of Polymer # " << count << ".\n";
        for (Particle p: pmer.chain){
            for (std::vector <int>::const_iterator i = p.coords.begin(); i!=p.coords.end(); ++i){
                dump_file << *i << " | "; 
            }
            dump_file << "\n"; 
        }
        count ++; 
    }
    dump_file <<"~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n";
    dump_file.close();

    

    return; 

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of dumpPositionsOfPolymer. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
//
// NAME OF FUNCTION: dumpEnergyOfGrid 
//
// PARAMETERS: (int step, std::string filename), and some attributes present in the Grid Object 
// 'step' is the current time step we are at. This is an integer which is likely defined in the driver code.  
// 'filename' is the file to which I am going to print out coordinate information about the polymer. 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers in the current Grid object \, and prints them 
// out to a text file.  
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 


void Grid::dumpEnergyOfGrid (int step, std::string filename, bool first_call){
    std::ofstream dump_file(filename, std::ios::app); 
    if (first_call){
        dump_file <<"This file contains energy of the Grid at certain points in the simulation.\n";
        dump_file <<"Energy | Step_Number" << std::endl;
    }
    
    dump_file <<this->Energy << " | " << step << std::endl;

    
    dump_file.close();

    

    return; 

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of dumpPositionsOfPolymer. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
//
// NAME OF FUNCTION: CreateGridObject 
//
// PARAMETERS: (std::string positions, std::string topology), and some attributes present in the Grid Object 
// 'positions' is the name of the file with coordinates of the polymer, so something like "coords.txt"
// 'topology' is the name of the file with the geometric and energetic information of the system, so something like "topology.txt"
// WHAT THE FUNCTION DOES: It looks at the topology file, which contains the geometric dimensions of the box 
// the polymer solution is in and the energy surface ie the interaction energies between the different particle 
// types present in the box. 
//
// DEPENDENCIES: ExtractTopologyFromFile, ExtractNumberOfPolymers, ExtractPolymersFromFile, create_lattice_pts, InstantiateOccupancyMap, CalculateEnergy. 
//
// THE CODE: 


Grid CreateGridObject(std::string positions, std::string topology){ 
    std::vector <double> info_vec; 
    std::vector <std::vector <double>> energy_mat; 
    info_vec = ExtractTopologyFromFile(topology);
    
    // Grid G(x, y, z, kT, Emm, Ems, Ess)
    Grid G (info_vec.at(0), info_vec.at(1), info_vec.at(2), info_vec.at(3), info_vec.at(4), info_vec.at(5), info_vec.at(6) ) ; 

    int N = G.ExtractNumberOfPolymers(positions);
    std::cout << "You have given us " << N << " polymers to work with."<<std::endl;
    G.ExtractPolymersFromFile(positions);

    // create solvent object... 

    std::vector <std::vector <int>> lattice = create_lattice_pts (G.x, G.y, G.z); 

    // begin defining solvent positions 
    // std::string type_s = "solvent"; 
    // get rid of all the locations that have been occupied by the polymer 

   // for (Polymer pmer: G.PolymersInGrid){
   //    for (Particle p: pmer.chain){
   //        lattice.erase(std::remove(lattice.begin(), lattice.end(), p.coords), lattice.end()); 
   //    }
   // }
    
    // now that those locations have been cleared from lattice, these will now become locations of the solvent 
    
    // std::vector <Particle> solvents; 
    // for (std::vector <int> loc: lattice){
    //    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    //    std::mt19937 generator(seed); 
    //    std::uniform_int_distribution<int> distribution (0,1); 
    //    Particle p (loc, type_s, distribution(generator));         
    //    solvents.push_back(p); 
    // } 

    // G.SolventInGrid = solvents; 

    

    std::cout << "Polymers extracted successfuly!" << std::endl;
    std::cout << "Solvent positions defined successfully!" << std::endl;
    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~"<< std::endl<<std::endl;
    // std::cout <<"Energy of the system is " << G.Energy << std::endl;

    G.instantiateOccupancyMap();
    G.CalculateEnergy();

    return G; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of CreateGridObject. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
//
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: instantiateOccupancyMap
//
// PARAMETERS: Some attributes of the Grid Object
// WHAT THE FUNCTION DOES: It looks at the PolymersInGrid object and SolventInGrid object and creates a map
// that takes in the location of a point on the lattice, and returns the particle occupying that location
//
// DEPENDENCIES: No custom functions required apart from those defined previously in object Grid. 
//
// THE CODE: 

void Grid::instantiateOccupancyMap(){
    
    for (Polymer pmer: this->PolymersInGrid){
        for (Particle p: pmer.chain){
            // Particle pdash = p; 
            // std::vector <int> ploc = p.coords; 
            this->OccupancyMap[p.coords] = p;
        }
    }
    
    return; 
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of instantiateOccupancyMap. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: checkValidityOfCoords
//
// PARAMETERS: std::vector <int> v, some attributes of the Grid Object
// WHAT THE FUNCTION DOES: It look at a singular location coordinate and checks if it is a "good" location, 
// and does not lie outside the box, or something unphysical like that. 

// DEPENDENCIES: No custom functions required apart from those defined previously in object Grid. 
//
// THE CODE: 

bool Grid::checkValidityOfCoords(std::vector <int> v){
    if (v.at(0)>this->x || v.at(0) < 0 || v.at(1)>this->y || v.at(1)<0 || v.at(2)>this->z || v.at(2)<0){
        return false;
    }
    else {
        return true;
    }
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of checkValidityOfCoords. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: checkForOverlaps
//
// PARAMETERS: std::vector <int> PolymerVector
// 
// WHAT THE FUNCTION DOES: It looks at a vector of polymers that are supposed to go into the Grid. If
// two (or more) monomer units occupy the same spot, it will freak out and exit out of compilation. 
// 
// DEPENDENCIES: No custom functions required apart from those defined previously in object Grid. 
//
// THE CODE: 


bool Grid::checkForOverlaps(std::vector <Polymer> PolymerVector){
    
    std::vector <std::vector <int>> loc_list; 

    for (Polymer pmer: PolymerVector){
        for (Particle p: pmer.chain){
            // check if element exists in vector 
                if (std::find(loc_list.begin(), loc_list.end(), p.coords) != loc_list.end() ){
                    std::cerr << "you have a repeated element." << std::endl;
                    // std::cout << "current element is: " << std::endl;
                    // print(p.coords); 
                    //print(loc_list);
                    return false; 
                    }
            
                else{
                    loc_list.push_back(p.coords);  
                }
            }
        }    
    
    std::cout << "Input file has no overlaps!" << std::endl;
    return true;

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of checkForOverlaps. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// NAME OF FUNCTION: checkConnectivity
//
// PARAMETERS: std::vector <int> PolymerVector
// 
// WHAT THE FUNCTION DOES: It looks at a vector of polymers that are supposed to go into the Grid. If
// two adjacent monomer units in the polymer vector do not have adjacent coordinates, you have bad connectivity.  
// 
// DEPENDENCIES: subtract_vectors 
//
// THE CODE: 

bool Grid::checkConnectivity(std::vector <Polymer> PolymerVector) {
    std::vector <int> d1 = {0, 0, 1}; 
    std::vector <int> dx = {0, 0, this->x-1}; 
    std::vector <int> dy = {0, 0, this->y-1}; 
    std::vector <int> dz = {0, 0, this->z-1}; 
    for (Polymer pmer: PolymerVector){
        size_t length = pmer.chain.size(); 
        std::vector <int> connection; 
        for (int i{1}; i<static_cast<int>(length); i++){
            
            connection = subtract_vectors(&(pmer.chain.at(i).coords), &(pmer.chain.at(i-1).coords));
            impose_pbc(&connection, this->x, this->y, this->z);
            std::sort(connection.begin(), connection.end()); 
            if (connection == d1 || connection == dx || connection == dy || connection == dz){
                continue; 
            }
            else {
                std::cerr << "Shit, you have bad connectivity inside one (or maybe more) polymers. Check input file." << std::endl;
                return false; 
            }

        }
    }

    std::cout << "Input polymers are well-connected!" << std::endl;
    return true;
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of checkConnectivity. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: ExtractNumberOfPolymers 
//
// PARAMETERS: std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the positions file, and extracts the number of polymers from that file. 
// 
// DEPENDENCIES: subtract_vectors 
//
// THE CODE: 


int Grid::ExtractNumberOfPolymers(std::string filename){
    std::regex start("START"); 
    std::regex end ("END"); 
    std::ifstream myfile (filename); 

    if ( !myfile ){
        std::cerr << "File named " << filename << " could not be opened!" << std::endl; 
        exit(EXIT_FAILURE);
    }

    std::string myString; 
    // std::vector <std::string> StartContents; 
    // std::vector <std::string> EndContents; 

    int numberOfStarts {0}, numberOfEnds {0}; 

    // std::cout << "This is from ExtractNumberOfPolymers..." << std::endl;
    if (myfile.is_open()) {
        while ( std::getline(myfile, myString) ) {
            
            // std::cout << myString << std::endl; 
            // std::getline(myfile, myString); // pipe file's content into stream 

            if (std::regex_search(myString, start)){
                numberOfStarts++; 
            }
            else if (std::regex_search(myString, end)){
                numberOfEnds++; 
            }
            else if ( myString.empty()){
                std::cerr << "ERROR: Empty line found. Bad positions file. " << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    if (numberOfStarts==numberOfEnds){
        return numberOfStarts;
    }
    else {
        std::cerr << "Number of starts is not the same as number of ends. Bad input file." << std::endl;
        return 0;
    }

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractNumberOfPolymers. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// NAME OF FUNCTION: ExtractContentFromFile
//
// PARAMETERS: std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the positions file, and extracts all the data from the file 
// in the form of a vector of strings.  
// 
// DEPENDENCIES: subtract_vectors 
//
// THE CODE: 

std::vector <std::string> Grid::ExtractContentFromFile(std::string filename){
    std::ifstream myfile (filename); 

    if ( !myfile ){
        std::cerr << "File named " << filename << " could not be opened!" << std::endl; 
        exit(EXIT_FAILURE);
    }

    std::string mystring; 
    std::vector <std::string> contents; 

    if (myfile.is_open() ){
        while ( std::getline(myfile, mystring) ) {
            // pipe file's content into stream 
            contents.push_back(mystring); 
        }
    }

    return contents; 

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractContentFromFile. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: ExtractContentFromFile
//
// PARAMETERS: std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the positions file, and extracts all the data from the file 
// in the form of a vector of strings.  
// 
// DEPENDENCIES: ExtractContentFromFile, makePolymer, checkValidityOfCoords 
//
// THE CODE: 


void Grid::ExtractPolymersFromFile(std::string filename){

    int NumberOfPolymers = this->ExtractNumberOfPolymers(filename); 

    std::vector <Polymer> PolymerVector; 
    PolymerVector.reserve(NumberOfPolymers);
    std::vector <std::vector <int>> locations; 

    std::vector <std::string> contents = this->ExtractContentFromFile(filename); // this extracts every line of the file

    std::regex start ("START"), end ("END"); 


    int startCount{0}, endCount {0}; 
    
    
    for (std::string s: contents){
        
        
         
        std::stringstream ss(s); 
        if (std::regex_search(s, start) ){
            startCount++;
            continue; 
        }

        else if (std::regex_search(s, end) ) {
            endCount++;
            
            Polymer pmer = makePolymer(locations);
            PolymerVector.push_back(pmer);
            
            locations.clear();
            
        }

        else{
            std::vector <int> loc; 
            for (int i=0; ss>>i; ){
                
                loc.push_back(i);

            }

            if (!this->checkValidityOfCoords(loc)){
            std::cerr << "Coordinates are out of bounds. Bad input file." << std::endl;
            exit(EXIT_FAILURE); 
            }
        
            locations.push_back(loc); 
            
        
        }
    }
    
    

    if(!(this->checkForOverlaps(PolymerVector))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Overlap detected." << std::endl;
        exit(EXIT_FAILURE);  
    }

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // throw in a check for connectivity of polymer chains 
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    if (!(this->checkConnectivity(PolymerVector))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Monomer units are not adjacent to one another on the lattice." << std::endl;
        exit(EXIT_FAILURE); 
    }
    
    this->PolymersInGrid = PolymerVector;

    return; 
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractContentFromFile. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: CalculateEnergy
//
// PARAMETERS: some attributes of the Grid Object
// WHAT THE FUNCTION DOES: Calculates energy of the current Grid. Critical to correctly evolve system. 
// includes nearest neighbor interactions with and without directional effects.  
//
// DEPENDENCIES: obtain_ne_list 
//
// OPTIMIZATION OPPORTUNITY: I am double counting monomer-monomer interactions. This can possibly be avoided. 
//
// THE CODE: 

void Grid::CalculateEnergy(){
    double Energy {0.0}; 
    
    // polymer-polymer interaction energies 
    for (Polymer pmer: this->PolymersInGrid){
        for (Particle p: pmer.chain){
            
            std::vector <Particle> part_vec = pmer.ConnectivityMap[p];
            
            std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z); // get neighbor list 

            // consider bonded interactions  
            
            // for (Particle ple: part_vec){
            //     ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), ple.coords), ne_list.end());
            // }

            // std::cout << "printing out neighboring monomer units included in energy calculation..." << std::endl; 
            // std::cout << "current monomer unit is "; print(p.coords);
            for (std::vector <int> loc: ne_list){
                // std::cout << "neighbor list is: "; print (loc); 
                if (this->OccupancyMap.find(loc) != this->OccupancyMap.end() ){ // loc is in the map 
                    // std::cout << "inside if loop is: "; print(loc); 

                    if (p.orientation == this->OccupancyMap[loc].orientation){
                        Energy += (this->Emm_a)*0.5;     
                    }
                    else {
                        Energy += (this->Emm_n)*0.5; 
                    }
                    
                }
                else {
                    
                    Energy += this->Ems; 
                }
            }

        }
    }

    // std::cout << "Energy is " << Energy << std::endl; 
    this->Energy = Energy; 
    return; 

}
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of calculateEnergy. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: ClusterParticleMaker
//
// PARAMETERS: a well-defined Grid  
// 
// WHAT THE FUNCTION DOES: it looks at the Polymers In Grid and charts out all the solvent molecules that interact 
// with the polymers. The idea is that you dont need to update the orientation of the entire box, because it doesnt 
// matter if solvent molecules not in contact with the polymer are flipped.    
// 
// DEPENDENCIES: obtain_ne_list, OccupancyMap
//
// THE CODE: 

std::vector <Particle> Grid::ClusterParticleMaker(){
    
    std::vector <Particle> Particles; 

    std::vector <std::vector <int>> solvent_locations;  

    for (Polymer pmer: this->PolymersInGrid){

        for (Particle p: pmer.chain){
 
                Particles.push_back (p); 
                std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z);

                // std::cout << "The polymer particle I added to list is: "; 
                // print(p.coords); 

                for (std::vector <int> v: ne_list){

                    if (this->OccupancyMap[v].ptype=="solvent" ){
                        
                        if (std::find(solvent_locations.begin(), solvent_locations.end(), this->OccupancyMap[v].coords) != solvent_locations.end() ){
                            
                        }

                        else {
                            solvent_locations.push_back(this->OccupancyMap[v].coords) ;
                            Particles.push_back(this->OccupancyMap[v]); 
                            // std::cout << "The solvent molecule is at: "; 
                            // print(this->OccupancyMap[v].coords); 
                        }
                    }
                }
        }
    }
    // std::cout << "Locations of particle to Ising flip are: " << std::endl;
    // for (auto part: Particles){
    //     print(part.coords); 
    // }

    return Particles; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ClusterParticleMaker. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: ClusterMaker
//
// PARAMETERS: a well-defined Grid  
// 
// WHAT THE FUNCTION DOES: Once the particles that can POTENTIALLY be part of a cluster have been identified, 
// the next thing to do is make the network. You pick a particle at random from the potential cluster, consider all
// the neighbors with in the cluster with the same orientation. You form a bond with those neighbors with a probability of 
// exp(-2\beta J). Once the network has been formed, flip orientation of every particle in network. 
// 
// DEPENDENCIES: ExtractContentFromFile, energy predictor 
//
// THE CODE: 
/*
std::vector <Particle> Grid::ClusterMaker(std::vector <Particle> Particles, std::vector <Particle> final, std::vector <Particle> to_send_, int count){

    if ( static_cast<int>(to_send_.size()) == 0 && count == 1){

        // std::cout << "There is nothing to send..." << std::endl;
        return final; 
    }

    if (count==0){
        
        // std::cout << "ENTRY POINT!" << std::endl;
        std::vector <Particle> to_send; 
        int r = rng_uniform(0, Particles.size()-1); 
        // std::cout << "Length of Particles is " << Particles.size() << std::endl;
        
        
        Particle p = Particles.at(r); 

        // std::cout << "The absolute first particle is "; 
        // print(p.coords);

        final.push_back(p); 
        std::vector <std::vector <int>> location_list; 

        for (Particle pa: Particles){
            location_list.push_back(pa.coords); 
        }

        std::vector <Particle> neighbors;

        std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z); 

        

        if (p.ptype=="solvent"){
            for (std::vector <int> v: ne_list){
                if (std::find(location_list.begin(), location_list.end(), v) != location_list.end() ){

                    if (this->OccupancyMap[v].orientation == p.orientation){
                        neighbors.push_back(this->OccupancyMap[v]); 
                    }
                }
            }
        }
        else {
            for (std::vector <int> v: ne_list){
                if (this->OccupancyMap[v].orientation == p.orientation){
                        neighbors.push_back(this->OccupancyMap[v]); 
                }
            }
        }

        // std::cout << "Number of good neighbors is " << neighbors.size() << std::endl;
        for (Particle pcle: neighbors){
            double coupling = 2*1/this->kT*EnergyPredictor(p, pcle); 
            // std::cout << "coupling is " << coupling << std::endl;
            double prob = std::exp(coupling); 
            // std::cout << "Probability of coupling is " << prob << std::endl;
            double check = rng_uniform(0.0,1.0);
            // std::cout << "check is " << check << std::endl;
            // std::cout << "for particle "; 
            // print(pcle.coords); 
            if (check < prob){
                // std::cout << "prob is " << prob << std::endl; 
                // std::cout << "check is " << check << std::endl;

                final.push_back(pcle); 
                to_send.push_back(pcle); 
                // std::cout << "neighbor approved :"; 
                // print(pcle.coords);
            }
            else {
                continue; 
            }
        }

        count = 1 ;
        
        std::cout << "This has to be hit only once!!!!!!!!!!!!!! " << std::endl; 
        std::cout << "Count is " <<count << std::endl;
        std::cout << "printing out first final"<<std::endl;
        for (auto ppp: final){
            print(ppp.coords); 

        }
        std::cout << "printing out first to_send" <<std::endl;
        for (auto pppp: to_send){
            print(pppp.coords);
        }
        
        final = this->ClusterMaker(Particles, final, to_send, count); 

    

    }

    else if (count==1){
        
        for (Particle p: to_send_){
            
            std::vector <std::vector <int>> location_list; // this will contain the original master cluster 
            std::vector <Particle> to_send; 

            for (Particle px: Particles){
                location_list.push_back(px.coords); 
            }
            std::vector <Particle> neighbors;

            std::vector <std::vector <int>> final_list; 

            for (Particle pf: final){
                final_list.push_back(pf.coords); 

            }

            std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z); // obtaining the neighbors of the main guy 

            if (p.ptype=="solvent"){
                for (std::vector <int> v: ne_list){
                    
                    if ( std::find(location_list.begin(), location_list.end(), v) != location_list.end() ){

                            if ( !(std::find ( final_list.begin(), final_list.end(), v) != final_list.end() ) ){
                        
                                if ( this->OccupancyMap[v].orientation == p.orientation){

                                    neighbors.push_back(this->OccupancyMap[v]); // added to potential neighbor list  
                                }
                        }
                    }
                }
            }
            else {
                for (std::vector <int> v: ne_list){
                    if ( !(std::find ( final_list.begin(), final_list.end(), v) != final_list.end() ) ){
                        if (this->OccupancyMap[v].orientation == p.orientation){
                            neighbors.push_back(this->OccupancyMap[v]);     // added to potential neighbor list 
                        }
                    }
                }
            }

            // std::cout << "For particle "; 
            // print(p.coords);
            // std::cout << "Number of good neighbors is " << neighbors.size() << std::endl;

            for (Particle pcle: neighbors){
                // std::cout << "Energy is " << EnergyPredictor(p, pcle) << std::endl;
                double coupling = 2*1/this->kT*EnergyPredictor(p, pcle); 
                double prob = std::exp(coupling); 
                // std::cout << "Probability of coupling is " << prob << std::endl;

                double check = rng_uniform(0.0,1.0);
                // std::cout << "check is " << check << std::endl;
                // std::cout << "for particle "; 
                // print(pcle.coords); 
                if (check < prob){
                    final.push_back(pcle); 
                    to_send.push_back(pcle); 
                }
                else {
                    continue; 
                }
            }
            count = 1;
            
            
            std::cout << "Given current particle as: ";
            print(p.coords);  
            std::cout << "printing out the next iteration of final:" << std::endl;

            for (auto ppp: final){
                std::cout << ppp.orientation << ", ";
                print(ppp.coords); 
            }
            
            std::cout << "size of final cluster is " << final.size() << std::endl;
            std::cout << "printing out new to_send" <<std::endl;
            for (auto pppp: to_send){
                print(pppp.coords);
            }
            

            final = this->ClusterMaker(Particles, final, to_send, count); 
        }



    }

    
    std::cout << "The particles in cluster are " << std::endl; 
    for (auto part: final){
        std::cout << "Particle type is " << part.ptype << ", location is: ";
        print(part.coords);
    }
    

    return final; 

}
*/
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ClusterParticleMaker. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#




/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                    POLYMER METHODS. 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//============================================================
//============================================================
// 
// NAME OF FUNCTION: printChainCoords
//
// PARAMETERS: none, method on a Polymer
// 
// WHAT THE FUNCTION DOES: Given a Polymer, it will print out coordinates of the polymers in a nice way 
//
// DEPENDENCIES: print
//
// THE CODE: 


void Polymer::printChainCoords(){
    for (Particle p: this->chain){
        p.printCoords(); 
    }
    return; 
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//                     End of printChainCoords. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: printOrientation
//
// PARAMETERS: none, method on a Polymer
// 
// WHAT THE FUNCTION DOES: Given a Polymer, it will print out Ising orientation of monomer units
//
// DEPENDENCIES: print
//
// THE CODE: 


void Polymer::printOrientation(){
    for (Particle p:this->chain){
        std::cout << p.orientation << " | ";
    }
    std::cout << std::endl;
    return;
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//                  End of printOrientation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// !!!!!!!! CRITICAL OBJECT FOR ACCURATE SIMULATION !!!!!!!
//
// NAME OF FUNCTION: ChainToConnectivityMap
//
// PARAMETERS: none, method on a Polymer
// 
// WHAT THE FUNCTION DOES: Given a Polymer, it will construct a connectivity map. Given a monomer unit, it 
// will give locations of particles it is bonded to. 
// 
// !!!!!! NOTE: THIS ONLY WORKS FOR LINEAR POLYMERS !!!!!!!!
//
// DEPENDENCIES: print
//
// THE CODE: 

void Polymer::ChainToConnectivityMap(){

    const int chainLength = this->chain.size(); 

    for (int i{0}; i<chainLength; i++){
        if (i==0){
            
            this->ConnectivityMap[this->chain.at(i)] = { (this->chain.at(i+1)),  };

        }
        else if (i==(chainLength-1)){

            this->ConnectivityMap[this->chain.at(i)] = { (this->chain.at(i-1)), }; 
        }
        else {
            this->ConnectivityMap[this->chain.at(i)] = { (this->chain.at(i-1)), (this->chain.at(i+1)) };
        }

    }

    return; 
}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//                End of ChainToConnectivityMap. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
//
// NAME OF FUNCTION: findKinks
//
// PARAMETERS: none, method on a Polymer
// 
// WHAT THE FUNCTION DOES: Given a Polymer, it will find kinks in the polymer structure.  
//
// DEPENDENCIES: subtract_vectors 
//
// THE CODE: 


std::vector <int> Polymer::findKinks(){
    const int chainLength = this->chain.size(); 
    std::vector <int> kink_indices; 

    // obtain all location where kinks exist in the polymer 
    for (int i{0}; i<chainLength-2; i++){
        std::vector <int> v1 = subtract_vectors(&(this->chain.at(i+1).coords), &(this->chain.at(i).coords) ); 
        std::vector <int> v2 = subtract_vectors(&(this->chain.at(i+2).coords), &(this->chain.at(i+1).coords) );

        if (v1==v2){
            continue;
        }
        else {
            kink_indices.push_back(i);
        }
    }

    // if (kink_indices.size() == 0){
    //    std::cout << "No kinks in polymer..." << std::endl;
    // }

    return kink_indices; 

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//                End of findKinks. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
//
// NAME OF FUNCTION: findCranks
//
// PARAMETERS: none, method on a Polymer
// 
// WHAT THE FUNCTION DOES: Given a Polymer, it will find kinks in the polymer structure.  
//
// DEPENDENCIES: subtract_vectors 
//
// THE CODE: 


std::vector <int> Polymer::findCranks(){
    const int chainLength = this->chain.size(); 
    std::vector <int> crank_indices; 

    // obtain all location where the crank exists in the polymer
    for (int i{0}; i< chainLength-3; i++){
        std::vector <int> v1 = subtract_vectors(&(this->chain.at(i+1).coords), &(this->chain.at(i).coords) ); 
        std::vector <int> v2 = subtract_vectors(&(this->chain.at(i+2).coords), &(this->chain.at(i+3).coords) );

        //

        if (v1==v2){
            crank_indices.push_back(i); 
        } 
        else {
            continue;
        }
    }

    if (crank_indices.size() == 0) {
        // std::cout << "there are no cranks in the current structure..." << std::endl;
    } 

    return crank_indices;
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//                End of findCranks. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                    PARTICLE METHODS. 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//============================================================
//============================================================
// 
// NAME OF FUNCTION: printCoords
//
// PARAMETERS: none, method on a particle
// 
// WHAT THE FUNCTION DOES: Given a particle, it will print out coordinates in a nice way 
//
// DEPENDENCIES: print
//
// THE CODE: 

void Particle::printCoords(){
    print(this->coords); 
    return; 
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of printCoords. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                    MONTE CARLO MOVES. 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//============================================================
//============================================================
// 
// NAME OF FUNCTION: IsingFlip
//
// PARAMETERS: a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform an IsingFlip on a certain cluster of particles. Clusters
// are made using ClusterParticleMaker + ClusterMaker. Clusters are flipped using ClusterFlip. 
//
// DEPENDENCIES: ClusterParticleMaker, ClusterMaker, ClusterFlip
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully.    
//
// UPDATE: Not sure how much I need this guy! Stay tuned. Jan 8, 2022. 
//
//
// THE CODE: 

/*
Grid IsingFlip(Grid InitialG){

    

    std::vector <Particle> Particles = InitialG.ClusterParticleMaker(); 
    

    std::vector <Particle> an_empty_vector_for_ClusterMaker_1;
    std::vector <Particle> an_empty_vector_for_ClusterMaker_2; 

    
    std::vector <Particle> cluster = InitialG.ClusterMaker(Particles, an_empty_vector_for_ClusterMaker_1, an_empty_vector_for_ClusterMaker_2, 0);

    

    ClusterFlip (&cluster);                       // flip the cluster! 

    Grid NewG (InitialG);                         // make a copy of the original Grid 

    for (Particle P: cluster){
        NewG.OccupancyMap[P.coords] = P;          // update the OccupancyMap of NewG

    }

    // update PolymersInGrid 
    for (Polymer& pmer: NewG.PolymersInGrid){
        for (Particle& p: pmer.chain){
            p = NewG.OccupancyMap[p.coords]; 
        }
    }

    // update SolventInGrid 
    for (Particle& p: NewG.SolventInGrid){
        p = NewG.OccupancyMap[p.coords]; 
    }

    return NewG; 


}
*/
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of IsingFlip. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: ZeroIndexRotation
//
// PARAMETERS: index of a polymer to perform ZeroIndexRotation, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a rotation of the monomer at the first index (the tail)
// of the polymer. As it stands, this function only rotates one molecule at the tail.
//
// PLANNED EXTENSION: Multiple molecules at the tail need to be rotated.    
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully.    
//
// THE CODE: 

Grid ZeroIndexRotation(Grid* InitialG, int index){

    Grid NewG (*InitialG);

    // get the neighborlist of particle at index 1 
    std::vector <int> loc_0 = NewG.PolymersInGrid.at(index).chain.at(0).coords; 
    std::vector <int> loc_1 = NewG.PolymersInGrid.at(index).chain.at(1).coords;
    std::vector <int> loc_2 = NewG.PolymersInGrid.at(index).chain.at(2).coords; 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(loc_1, (*InitialG).x, (*InitialG).y, (*InitialG).z) ; 
    // std::cout << "neighbort list is " << std::endl; 
    // print(ne_list); 
    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0 ), ne_list.end() );     // gets rid of the locations that clearly can't be swung into 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2 ), ne_list.end() );     // gets rid of the locations that clearly can't be swung into
    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine());
    // find a location that is not occupied by monomer 
    // use the occupancy map 
	size_t tries = 0; 
    for (std::vector <int> to_rot: ne_list){

        if (!( NewG.OccupancyMap.find(to_rot) != (NewG).OccupancyMap.end() ) ){

            // update OccupancyMap
             
            Particle p1 (NewG.OccupancyMap.at(loc_0));
            p1.coords = to_rot;   
            NewG.OccupancyMap[to_rot] = p1;  
            
            int p_erase = NewG.OccupancyMap.erase(loc_0); 
            p_erase++;
            // update positions in polymers in grid 
            NewG.PolymersInGrid.at(index).chain.at(0).coords = to_rot; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            break;

        }
		else {
			tries++; 
		}

    }
	
	if (tries == ne_list.size() ){
		std::cout << "no place to rotate! configuration will be accepted by default!" << std::endl;
	}
    return NewG;
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ZeroIndexRotation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: FinalIndexRotation
//
// PARAMETERS: index of polymer to perform FinalIndexRotation on, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a rotation of the monomer at the final index (the head)
// of the polymer. As it stands, this function only rotates one molecule at the head.
//
// PLANNED EXTENSION: Multiple molecules at the head need to be rotated.    
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 


Grid FinalIndexRotation(Grid* InitialG, int index){

    Grid NewG (*InitialG);

    int dop = NewG.PolymersInGrid.at(index).chain.size(); 

    // get the neighborlist of particle at index DoP-2
    std::vector <int> loc_0 = NewG.PolymersInGrid.at(index).chain.at(dop-1).coords;
    std::vector <int> loc_1 = NewG.PolymersInGrid.at(index).chain.at(dop-2).coords;
    std::vector <int> loc_2 = NewG.PolymersInGrid.at(index).chain.at(dop-3).coords;

    // obtain neighbor list 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(loc_1, NewG.x, NewG.y, NewG.z); 
    // std::cout << "neighbort list is " << std::endl; 
    // print(ne_list); 
    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0), ne_list.end() ); 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2), ne_list.end() );
    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine());

    // find a location that is unoccupied by monomer 
    // use the occupancy map  
	size_t tries = 0; 
    for (std::vector <int> to_rot: ne_list){

        if ( !(NewG.OccupancyMap.find(to_rot) != NewG.OccupancyMap.end()) ) {

            // update OccupancyMap 
            Particle p1 (NewG.OccupancyMap[loc_0]); 
            p1.coords = to_rot; 

            NewG.OccupancyMap[to_rot] = p1; 

            int p_erase = NewG.OccupancyMap.erase(loc_0); 
            p_erase++; 

            // update positions in polymers in grid 
            NewG.PolymersInGrid.at(index).chain.at(dop-1).coords = to_rot; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            break;            

        }

		else {
		tries++; 
		}
    }

	if (tries == ne_list.size() ){
		std::cout << "no place to rotate! position will be accepted by default!" << std::endl;
	}
    return NewG;
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of FinalIndexRotation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: EndRotation
//
// PARAMETERS: index of polymer to perform EndRotation on, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a rotation of the monomer at the terminal index (the tail or head)
// of the polymer. As it stands, this function only rotates one molecule at the tail or head. The choice of head 
// or tail rotation comes from a random distribution. 
//
// PLANNED EXTENSION: Multiple molecules at the termini need to be rotated.    
//
// DEPENDENCIES: ZeroIndexRotation, FinalIndexRotation
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 


Grid EndRotation(Grid* InitialG, int index){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return ZeroIndexRotation(InitialG, index); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return FinalIndexRotation(InitialG, index); 

    }
    

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of EndRotation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// NAME OF FUNCTION: KinkJump
//
// PARAMETERS: index of which polymer to perform KinkJump on, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a kink jump when it finds a kink. 
// The kinks are shuffled before any is chosen. 
//
// DEPENDENCIES: findKinks, ChainToConnectivityMap, add_vectors, subtract_vectors
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 


Grid KinkJump(Grid* InitialG, int index){

    Grid NewG (*InitialG); 

    std::vector <int> k_idx = NewG.PolymersInGrid.at(index).findKinks(); 

    if (k_idx.size() == 0 ){
        // std::cout << "No kinks found in polymer..." << std::endl;
        return NewG;
    }

    std::shuffle(std::begin(k_idx), std::end(k_idx), std::default_random_engine() ); 
	size_t tries = 0; 
    for (int idx: k_idx){

        // std::cout << "idx right before kink spot is " << idx << std::endl; 

        std::vector <int> d1 = subtract_vectors(&(NewG.PolymersInGrid.at(index).chain.at(idx+1).coords), &(NewG.PolymersInGrid.at(index).chain.at(idx).coords) );
        std::vector <int> d2 = subtract_vectors(&(NewG.PolymersInGrid.at(index).chain.at(idx+2).coords), &(NewG.PolymersInGrid.at(index).chain.at(idx+1).coords) ); 

        std::vector <int> to_check = add_vectors( &(NewG.PolymersInGrid.at(index).chain.at(idx).coords), &d2); 
        impose_pbc(&to_check, (NewG).x, (NewG).y, (NewG).z); 

        if ( !( NewG.OccupancyMap.find(to_check) != NewG.OccupancyMap.end()) ){

            //update occupancy map
            Particle p1 ( NewG.OccupancyMap.at( NewG.PolymersInGrid.at(index).chain.at(idx+1).coords) );  // this will be the monomer particle
            p1.coords = to_check; 
            
            NewG.OccupancyMap[to_check] = p1; 

            int p_erase = NewG.OccupancyMap.erase( NewG.PolymersInGrid.at(index).chain.at(idx+1).coords ); 
            p_erase++;
			// std::cout << "kink will jump to: "; 
			// print(p1.coords); 
            // update PolymersInGrid 
            NewG.PolymersInGrid.at(index).chain.at(idx+1) = p1; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            break;

        }
        else {
			tries++;
            // std::cout << "This spot is taken by a monomer... no kink jumping." << std::endl;
        }
	
    }
	
	if (tries == k_idx.size()) {
		std::cout << "No spot was available to kink jump! Position will be accepted by default!" << std::endl;
	}
    return NewG; 

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of KinkJump. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: CrankShaft
//
// PARAMETERS: index of a polymer to perform CrankShaft on, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a kink jump when it finds a kink. 
// The kinks are shuffled before any is chosen. 
//
// DEPENDENCIES: findKinks, ChainToConnectivityMap, add_vectors, subtract_vectors
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 


Grid CrankShaft(Grid* InitialG, int index){

    Grid NewG(*InitialG); 

    std::vector <int> c_idx = NewG.PolymersInGrid.at(index).findCranks(); 

    if ( c_idx.size()==0 ){
        // std::cout << "No cranks in this polymer..." << std::endl; 
        return NewG; 
    }

    std::shuffle (std::begin (c_idx), std::end(c_idx), std::default_random_engine() ); 
	size_t tries = 0; 
    for (int idx: c_idx){

        std::vector <int> HingeToKink = subtract_vectors(&(NewG.PolymersInGrid.at(index).chain.at(idx+2).coords), &(NewG.PolymersInGrid.at(index).chain.at(idx+3).coords) ); 
        std::vector <int> HingeToHinge = subtract_vectors(&(NewG.PolymersInGrid.at(index).chain.at(idx+3).coords ), &(NewG.PolymersInGrid.at(index).chain.at(idx).coords ) );

        std::vector <std::vector <int>> drns = HingeSwingDirections (& (HingeToHinge), &(HingeToKink), NewG.x, NewG.y, NewG.z); 
        int choice = rng_uniform(0,2); 

        std::vector <int> d1 = drns.at(choice);  //subtract_vectors( &(NewG.PolymersInGrid.at(index).chain.at(idx+3).coords), &(NewG.PolymersInGrid.at(index).chain.at(idx+2).coords) );

        std::vector <int> to_check_1 = add_vectors ( &(NewG.PolymersInGrid.at(index).chain.at(idx).coords), &d1 ); 
        std::vector <int> to_check_2 = add_vectors ( &(NewG.PolymersInGrid.at(index).chain.at(idx+3).coords), &d1 ); 
        impose_pbc(&to_check_1, NewG.x, NewG.y, NewG.z); 
        impose_pbc(&to_check_2, NewG.x, NewG.y, NewG.z);  

        if ( !(NewG.OccupancyMap.find(to_check_1) != NewG.OccupancyMap.end()) && !(NewG.OccupancyMap.find(to_check_2) != NewG.OccupancyMap.end() ) ){

            // update OccupancyMap
            Particle p1 ( NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+1).coords] );
            p1.coords = to_check_1;

            Particle p2 ( NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+2).coords]); 
            p2.coords = to_check_2; 
            

            NewG.OccupancyMap[to_check_1] = p1; 
            NewG.OccupancyMap[to_check_2] = p2; 

            int p_erased_1 = NewG.OccupancyMap.erase( NewG.PolymersInGrid.at(index).chain.at(idx+1).coords );
            int p_erased_2 = NewG.OccupancyMap.erase( NewG.PolymersInGrid.at(index).chain.at(idx+2).coords );
            p_erased_1++; 
            p_erased_2++;

            // update PolymersInGrid 
            NewG.PolymersInGrid.at(index).chain.at(idx+1).coords = p1.coords; 
            NewG.PolymersInGrid.at(index).chain.at(idx+2).coords = p2.coords; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 
 
            break;

        }
        else {
			tries++; 
            // std::cout << "This position is occupied by monomer... No cranking." << std::endl;
        }


    }
	if (tries == c_idx.size()){
		std::cout << "No position was available for a crankshaft! position will be accepted by default!" << std::endl;
	}
    return NewG;

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of CrankShaft. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
//
// !!!!!!! 
// NOTE FOR REPTATION CODE: The code is about all updating only the positions 
// of each particle in PolymersInGrid. One of the termini will move to a spot occupied by solvent. 
// If the monomer at index 0 is moving to a new spot, monomer at index 1 will take on coordinates of monomer 0 in the initial config, 
// monomer at index 2 will take coordinates of monomer at index 1 in the initial config, and so on. 
// Reverse the process if monomer at final index is moving to a new spot. 
// !!!!!!!
// 
// NAME OF FUNCTION: ForwardReptation 
//
// PARAMETERS: index of a polymer to reptate forward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: index of polymer on which to perform ForwardReptation, given a Grid, it will perform a forward reptation. 
// which means that the final index (size-1) is moving somewhere in its vicinity. 
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

Grid ForwardReptation(Grid* InitialG, int index){

    Grid NewG (*InitialG); 
    // get size of polymer chain 
    int size = NewG.PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list( NewG.PolymersInGrid.at(index).chain.at(size-1).coords, NewG.x, NewG.y, NewG.z ); 

    // get rid of second to last monomer position from ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), NewG.PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() ); 

	size_t tries = 0; 
    for (std::vector <int> v: ne_list){
        if (NewG.OccupancyMap[v].ptype != "monomer"){
            Particle p = NewG.OccupancyMap[v];
            p.ptype = "monomer"; 
            
            int p_erase = NewG.OccupancyMap.erase(NewG.PolymersInGrid.at(index).chain.at(0).coords) ; // there is no particle at index 0 
            p_erase++;

            for (int i{0}; i < size; i++){

                if (i < size - 1){
                    // the particles are staying the same, they are simply changing locations. 
                    NewG.PolymersInGrid.at(index).chain.at(i).orientation = (*InitialG).PolymersInGrid.at(index).chain.at(i).orientation;
                    NewG.PolymersInGrid.at(index).chain.at(i).coords = (*InitialG).PolymersInGrid.at(index).chain.at(i+1).coords;  
                    NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(i).coords] = NewG.PolymersInGrid.at(index).chain.at(i);


                }
                else if (i == size-1){
                    NewG.PolymersInGrid.at(index).chain.at(i).orientation = (*InitialG).PolymersInGrid.at(index).chain.at(i).orientation;
                    NewG.PolymersInGrid.at(index).chain.at(i).coords = v;
                    NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(i).coords] = NewG.PolymersInGrid.at(index).chain.at(i);
                }

            }

            // update polymer connectivity maps in grid 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 
            break;

        }
        else {
			tries++; 
            std::cout << "No place to slither... " << std::endl;
        }
    }
	if (tries == ne_list.size() ){
		std::cout << "There is no place to slither! Will be accepted by default!" << std::endl;
	}
    return NewG; 

} 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ForwardReptation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// !!!!!!! 
// NOTE FOR REPTATION CODE: The code is about all updating only the positions 
// of each particle in PolymersInGrid. One of the termini will move to a spot occupied by solvent. 
// If the monomer at index 0 is moving to a new spot, monomer at index 1 will take on coordinates of monomer 0 in the initial config, 
// monomer at index 2 will take coordinates of monomer at index 1 in the initial config, and so on. 
// Reverse the process if monomer at final index is moving to a new spot. 
// !!!!!!!
//
// NAME OF FUNCTION: BackwardReptation 
//
// PARAMETERS: index of a polymer to reptate forward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: index of polymer on which to perform BackwardReptation, given a Grid, it will perform backward reptation,
// which means that the particle at index 0 is moving somewhere in its vicinity.  
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

Grid BackwardReptation(Grid* InitialG, int index){
    
    Grid NewG(*InitialG); 

    // get initial size of the polymer chain 
    int size = NewG.PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list( NewG.PolymersInGrid.at(index).chain.at(0).coords, NewG.x, NewG.y, NewG.z); 

    // get rid of second monomer from the ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), NewG.PolymersInGrid.at(index).chain.at(1).coords), ne_list.end() ); 
	size_t tries = 0; 
    for (std::vector <int> v: ne_list){

        if (NewG.OccupancyMap[v].ptype != "monomer"){
            Particle p = NewG.OccupancyMap[v]; 
            p.ptype = "monomer"; 
            
            int p_erase = NewG.OccupancyMap.erase((*InitialG).PolymersInGrid.at(index).chain.at( size-1 ).coords); // 
            p_erase++;

            for (int i{0}; i < size; i++){

                if (i != 0){
                    NewG.PolymersInGrid.at(index).chain.at(i).orientation = (*InitialG).PolymersInGrid.at(index).chain.at(i).orientation;
                    NewG.PolymersInGrid.at(index).chain.at(i).coords = (*InitialG).PolymersInGrid.at(index).chain.at(i-1).coords;  
                    NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(i).coords] = NewG.PolymersInGrid.at(index).chain.at(i);
                }
                else if (i == 0){
                    NewG.PolymersInGrid.at(index).chain.at(i).orientation = (*InitialG).PolymersInGrid.at(index).chain.at(i).orientation;
                    NewG.PolymersInGrid.at(index).chain.at(i).coords = v;
                    NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(i).coords] = NewG.PolymersInGrid.at(index).chain.at(i);
                }

            }

            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 
            break;
        }

        else {
            // std::cout << "No place to slither..." << std::endl;
        }

    }

	if (tries == ne_list.size() ){
		std::cout << "No place to slither! Position will be accepted by default!" << std::endl;
	}

    return NewG; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ForwardReptation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
//
// !!!!!!! 
// NOTE FOR REPTATION CODE: The code is about all updating only the positions 
// of each particle in PolymersInGrid. One of the termini will move to a spot occupied by solvent. 
// If the monomer at index 0 is moving to a new spot, monomer at index 1 will take on coordinates of monomer 0 in the initial config, 
// monomer at index 2 will take coordinates of monomer at index 1 in the initial config, and so on. 
// Reverse the process if monomer at final index is moving to a new spot. 
// !!!!!!!
// 
// NAME OF FUNCTION: Reptation 
//
// PARAMETERS: index of a polymer to reptate forward or backward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: index of polymer on which to perform Reptation, given a Grid, it will perform a kink jump when it finds a kink. 
// The kinks are shuffled before any is chosen. 
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

Grid Reptation(Grid* InitialG, int index){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return BackwardReptation(InitialG, index); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return ForwardReptation(InitialG, index); 

    }
    
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of Reptation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// NAME OF FUNCTION: MoveChooser 
//
// PARAMETERS: a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a certain Monte Carlo move on the Grid. 
// The move could be from any of the following: 
// 1. End Rotation
// 2. Kink Jump
// 3. Crank Shaft
// 4. Reptation
// 5. Ising Flip
//
// DEPENDENCIES: rng_uniform, CalculateEnergy, EndRotation, KinkJump, CrankShaft, Reptation, IsingFlip  
//
// THE CODE: 

Grid MoveChooser(Grid* InitialG,  bool v){

    int index = rng_uniform(0, static_cast<int>((*InitialG).PolymersInGrid.size())-1); 
    // std::cout << "Index of polymer in grid to move is " << index << "." << std::endl; 
    Grid G_ ; 
    int r = rng_uniform(1, 4);
    switch (r) {
        case (1):
            if (v){
               std::cout << "Performing end rotations." << std::endl; 
            }
            // 
            G_ = EndRotation(InitialG, index);
            G_.CalculateEnergy(); 
            break;     
        
        case (2):
            if (v){
               std::cout << "Performing crank shaft." << std::endl; 
            }
            // std::cout << "Performing crank shaft." << std::endl;
            G_ = CrankShaft(InitialG, index);
            G_.CalculateEnergy();
            break; 

        case (3):
            if (v){
               std::cout << "Performing reptation." << std::endl; 
            }
            // std::cout << "Performing reptation." << std::endl;
            G_ = Reptation(InitialG, index); 
            G_.CalculateEnergy();
            break; 

        case (4):
            if (v){
               std::cout << "Performing kink jump." << std::endl; 
            }
            // std::cout << "Performing kink jump." << std::endl;
            G_ = KinkJump(InitialG, index);
            G_.CalculateEnergy ( );        
            break; 
    }

    return G_;
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of MoveChooser. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
