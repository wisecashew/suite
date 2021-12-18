#include <iostream> 
#include <vector>
#include <string> 
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

// Methods for class grid 


void Grid::dumpPositionsOfPolymers (int step){
    std::ofstream dump_file("dumpfile.txt", std::ios::app); 
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

Grid CreateGridObject(std::string positions, std::string topology){ 
    std::vector <double> info_vec; 
    std::vector <std::vector <double>> energy_mat; 
    info_vec = ExtractTopologyFromFile(topology);
    
    // Grid G(x, y, z, kT, Emm, Ems, Ess)
    Grid G (info_vec.at(0), info_vec.at(1), info_vec.at(2), info_vec.at(3), info_vec.at(4), info_vec.at(5), info_vec.at(6), info_vec.at(7)) ; 

    int N = G.ExtractNumberOfPolymers(positions);
    std::cout << "You have given us " << N << " polymers to work with."<<std::endl;
    G.ExtractPolymersFromFile(positions);

    // create solvent object... 

    std::vector <std::vector <int>> lattice = create_lattice_pts (G.x, G.y, G.z); 

   // begin defining solvent positions 
   std::string type_s = "solvent"; 
    // get rid of all the locations that have been occupied by the polymer 

   for (Polymer pmer: G.PolymersInGrid){
       for (Particle p: pmer.chain){
           lattice.erase(std::remove(lattice.begin(), lattice.end(), p.coords), lattice.end()); 
       }
   }
    
    // now that those locations have been cleared from lattice, these will now become locations of the solvent 
    
    std::vector <Particle> solvents; 
    for (std::vector <int> loc: lattice){
        unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
        std::mt19937 generator(seed); 
        std::uniform_int_distribution<int> distribution (0,1); 
        Particle p (loc, type_s, distribution(generator));         
        solvents.push_back(p); 
    } 

    G.SolventInGrid = solvents; 

    

    std::cout << "Polymers extracted successfuly!" << std::endl;
    std::cout << "Solvent positions defined successfully!" << std::endl;
    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~"<< std::endl<<std::endl;
    // std::cout <<"Energy of the system is " << G.Energy << std::endl;

    return G; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


void Grid::instantiateOccupancyMap(){
    
    for (Polymer pmer: this->PolymersInGrid){
        for (Particle p: pmer.chain){
            // Particle pdash = p; 
            // std::vector <int> ploc = p.coords; 
            this->OccupancyMap[p.coords] = p;
        }
    }

    for (Particle p: this->SolventInGrid){
        this->OccupancyMap[p.coords] = p; 
    }   
    
    return; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


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


bool Grid::checkConnectivity(std::vector <Polymer> PolymerVector) {
    std::vector <int> d1 = {0, 0, 1}; 
    std::vector <int> dx = {0, 0, this->x-1}; 
    std::vector <int> dy = {0, 0, this->y-1}; 
    std::vector <int> dz = {0, 0, this->z-1}; 
    for (Polymer pmer: PolymerVector){
        size_t length = pmer.chain.size(); 
        std::vector <int> connection; 
        for (int i{1}; i<length; i++){
            
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


int Grid::ExtractNumberOfPolymers(std::string filename){
    std::regex start("START"); 
    std::regex end ("END"); 
    std::ifstream myfile (filename); 

    std::string myString; 
    // std::vector <std::string> StartContents; 
    // std::vector <std::string> EndContents; 

    int numberOfStarts {0}, numberOfEnds {0}; 

    if (myfile.is_open()) {
        while (myfile.good()) {
            std::getline(myfile, myString); // pipe file's content into stream 

            if (std::regex_search(myString, start)){
                numberOfStarts++; 
            }
            else if (std::regex_search(myString, end)){
                numberOfEnds++; 
            }
            else if ( myString.empty()){
                std::cerr << "ERROR: Empty line found. Bad POSITIONS.TXT file. " << std::endl;
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



std::vector <std::string> Grid::ExtractContentFromFile(std::string filename){
    std::ifstream myfile (filename); 
    std::string mystring; 
    std::vector <std::string> contents; 

    if (myfile.is_open() ){
        while (myfile.good()) {
            std::getline(myfile, mystring); // pipe file's content into stream 
            contents.push_back(mystring); 
        }
    }

    return contents; 

};


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



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
        
            // print(loc);
            // std::cout << "checking coordinates of polymers..." << std::endl;
            // print(loc);

            locations.push_back(loc); 
        
        
        }
    }



    
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    
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
    // this->CalculateEnergy();
    return; 
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#




void Grid::CalculateEnergy(){
    double Energy {0.0}; 
    
    // polymer-polymer interaction energies 
    for (Polymer pmer: this->PolymersInGrid){
        for (Particle p: pmer.chain){
            
            std::vector <Particle> part_vec = pmer.ConnectivityMap[p];

            std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z); // get neighbor list 

            // curate neighbor list - get rid of neighbors you are connected to 
            
            for (Particle ple: part_vec){
                ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), ple.coords), ne_list.end());
            }


            for (std::vector <int> loc: ne_list){

                std::string type = this->OccupancyMap[loc].ptype; 
                if (type=="monomer"){
                    Energy += (this->Emm)*0.5; 
                }
                else {
                    
                    if (p.orientation == this->OccupancyMap[loc].orientation){
                        Energy += (this->Ems_a); 
                    }
                    else{
                        Energy += (this->Ems_n); 
                    }
                }
            }

        }
    }

    for (Particle p: this->SolventInGrid){
        std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z); // get neighbor list 
        for (std::vector <int> loc: ne_list){
            if (this->OccupancyMap[loc].ptype=="solvent"){
                Energy += 0.5*(this->Ess);
            }
            else {
                continue; 
            }
        }

    }

    this->Energy = Energy; 
    return; 

}


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
                            // std::cout << "this solvent location has been accounted for." << std::endl;
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


    return Particles; 
}



std::vector <Particle> Grid::ClusterMaker(std::vector <Particle> Particles, std::vector <Particle> final, int count){

    if (Particles.size() == 0 ){
        return final; 
    }

    if (count==0){

        std::vector <Particle> to_send; 
        int r = rng_uniform(0, Particles.size()); 

        Particle p = Particles.at(r); 

        std::cout << "I am just making the cluster here. The orientation of the first particle is " << p.orientation << std::endl << std::endl;

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


        for (Particle pcle: neighbors){
            double coupling = 2*1/this->kT*EnergyPredictor(p, pcle); 
            // std::cout << "coupling is " << coupling << std::endl;
            double prob = std::exp(coupling); 
            double check = rng_uniform(0.0,1.0);

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
    for (Particle p: final)
    this->ClusterMaker(to_send, final, 1); 

    }

    else if (count==1){
        // std::cout << "print at line 525, classes.cpp" << std::endl;
        for (Particle p: Particles){
            
            std::vector <std::vector <int>> location_list; 
            std::vector <Particle> to_send; 

            for (Particle px: Particles){
                location_list.push_back(px.coords); 
            }
            std::vector <Particle> neighbors;

            std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z); 

            if (p.ptype=="solvent"){
                for (std::vector <int> v: ne_list){
                    if (std::find(location_list.begin(), location_list.end(), v) != location_list.end() ){
                        // std::cout << "particle in the cluster located." << std::endl;
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


            for (Particle pcle: neighbors){
                double coupling = 2*1/this->kT*EnergyPredictor(p, pcle); 
                double prob = std::exp(coupling); 

                if (rng_uniform(0.0,1.0) < prob){
                    final.push_back(pcle); 
                    to_send.push_back(pcle); 
                }
                else {
                    continue; 
                }
            }

            this->ClusterMaker(to_send, final, 1); 
        }



    }

    return final; 

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//void Grid::ClusterMove(){

    // cluster particles 
    // flip them with a probability of half 
    // check if this configuration is accepted by Metropolis 
    // update Grid object 

//}


// this is a function that clusters particles to run the swendsen-wang algorithm on them 










//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


/*
void Grid::ZeroIndexRotation(int index){

    // get the neighborlist of particle at index 1 
    std::vector <int> loc_0 = this->PolymersInGrid.at(index).chain.at(0).coords; 
    std::vector <int> loc_1 = this->PolymersInGrid.at(index).chain.at(1).coords; 
    std::vector <int> loc_2 = this->PolymersInGrid.at(index).chain.at(2).coords; 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0 ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2 ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the occupancy map 

    for (std::vector <int> to_rot: ne_list){
        if (this->OccupancyMap[to_rot]==0){

            // update OccupancyMap 
            this->OccupancyMap[to_rot] = 1; // particle goes into new location 
            this->OccupancyMap[loc_0] = 0; 

            // update position in polymer 
            this->PolymersInGrid.at(index).chain.at(0).coords = to_rot; 
            this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
            return;
        }
        else {
            std::cout << "this location is occupied..." << std::endl;
        }

    }

    return; 
}
*/
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

/*void Grid::ZeroIndexRotation_MC(int index){
    // create a dummy grid object 

    Grid G_temp (this->x, this->y, this->z, this->kT, this->Emm, this->Ems, this->Ess);
    G_temp.PolymersInGrid = this->PolymersInGrid;
    G_temp.OccupancyMap = this->OccupancyMap; 

    // get energy of current system
    double Ei = this->Energy; 

    std::cout << "Initial energy is " << Ei << std::endl;

    // get the neighborlist of particle at index 1 
    std::vector <int> loc_0 = this->PolymersInGrid.at(index).chain.at(0).coords; 
    std::vector <int> loc_1 = this->PolymersInGrid.at(index).chain.at(1).coords; 
    std::vector <int> loc_2 = this->PolymersInGrid.at(index).chain.at(2).coords; 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z);

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0 ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2 ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the occupancy map 

    int choice = rng_uniform(0, static_cast<int>(ne_list.size()-1));
    std::vector <int> to_rot = ne_list.at(choice); 

    // for (std::vector <int> to_rot: ne_list){//
        if (this->OccupancyMap[to_rot]==0){
            // make the copy 
            G_temp.PolymersInGrid.at(index).chain.at(0).coords=to_rot;
            G_temp.OccupancyMap[to_rot] = 1; 
            G_temp.OccupancyMap[loc_0] = 0;
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 
            // get the new energy 
            G_temp.CalculateEnergy(); 
            double Ef = G_temp.Energy; 

            std::cout << "Test config energy is " << Ef << std::endl;

            if (acceptance (Ef-Ei, this->kT)){
                // update OccupancyMap 
                this->OccupancyMap[to_rot] = 1; // particle goes into new location 
                this->OccupancyMap[loc_0] = 0; 

                // update position in polymer 
                this->PolymersInGrid.at(index).chain.at(0).coords = to_rot; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
                std::cout << "monte carlo said yes!" << std::endl;
                this->CalculateEnergy(); 
                return;
            }
            else {
                    // G_temp.PolymersInGrid = this->PolymersInGrid;
                    // G_temp.OccupancyMap = this->OccupancyMap; 
                    std::cout << "monte carlo said no..." << std::endl;
                    return;
            }
        }
        else {
            // print(to_rot);
            std::cout << "configuration contains overlaps. move rejected..." << std::endl;
        }

    // }//

    return; 
}   */ 







//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

/*
void Grid::FinalIndexRotation(int index){
    int DoP = this->PolymersInGrid.at(index).chain.size(); 

    // get this particles neighbor list 
    std::vector <int> loc_0 = this->PolymersInGrid.at(index).chain.at(DoP-1).coords; 
    std::vector <int> loc_1 = this->PolymersInGrid.at(index).chain.at(DoP-2).coords; 
    std::vector <int> loc_2 = this->PolymersInGrid.at(index).chain.at(DoP-3).coords; 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0), ne_list.end() ); 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2), ne_list.end() ); 

    // find a location that is unoccupied 
    // use the occupancy map 

    for (std::vector <int> to_rot: ne_list){
        if (this->OccupancyMap[to_rot]==0){

            // update OccupancyMap 
            this->OccupancyMap[to_rot] = 1; 
            this->OccupancyMap[loc_0] = 0; 

            // update connectivity map 

            // update position in polymer and connectivity map 
            this->PolymersInGrid.at(index).chain.at(DoP-1).coords = to_rot; 
            this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
            return;
        }

        else {
            std::cout << "this location is occupied..." << std::endl;
        }
    }

    return; 

}
*/
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

/*
void Grid::FinalIndexRotation_MC(int index){
    int DoP = this->PolymersInGrid.at(index).chain.size(); 

    // create a dummy grid object 

    Grid G_temp (this->x, this->y, this->z, this->kT, this->Emm, this->Ems, this->Ess);
    G_temp.PolymersInGrid = this->PolymersInGrid; 
    G_temp.OccupancyMap = this->OccupancyMap; 

    // get energy of the current system
    double Ei = this->Energy; 

    std::cout << "Initial energy is " << Ei << std::endl;

    // get this particles neighbor list 
    std::vector <int> loc_0 = this->PolymersInGrid.at(index).chain.at(DoP-1).coords; 
    std::vector <int> loc_1 = this->PolymersInGrid.at(index).chain.at(DoP-2).coords; 
    std::vector <int> loc_2 = this->PolymersInGrid.at(index).chain.at(DoP-3).coords; 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(loc_1, this->x, this->y, this->z); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0), ne_list.end() ); 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2), ne_list.end() ); 

    // find a location that is unoccupied 
    // use the occupancy map 

    int choice = rng_uniform(0, static_cast<int>(ne_list.size()-1));
    std::vector <int> to_rot = ne_list.at(choice); 


    // for (std::vector <int> to_rot: ne_list){ //
        if (this->OccupancyMap[to_rot]==0){
            // make the copy 
            G_temp.PolymersInGrid.at(index).chain.at(DoP-1).coords = to_rot; 
            G_temp.OccupancyMap[to_rot] = 1; 
            G_temp.OccupancyMap[loc_0] = 0; 
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // get the new energy 
            G_temp.CalculateEnergy(); 
            double Ef = G_temp.Energy;

            std::cout << "Test config energy is " << Ef << std::endl;
            // std::cout << "Chain configuration is" <<std::endl;
            // G_temp.PolymersInGrid.at(index).printChainCoords();

            if (acceptance (Ef-Ei, this->kT)){
                // update OccupancyMap 
                this->OccupancyMap[to_rot] = 1; 
                this->OccupancyMap[loc_0] = 0; 

                // update position in polymer and connectivity map 
                this->PolymersInGrid.at(index).chain.at(DoP-1).coords = to_rot; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
                this->CalculateEnergy(); 
                std::cout << "monte carlo said yes!" << std::endl;
                return;
            }
            else {
                // G_temp.PolymersInGrid = this->PolymersInGrid; 
                // G_temp.OccupancyMap = this->OccupancyMap; 
                std::cout << "monte carlo said no..." << std::endl;
                return; 
            }
        }

        else {
            // print(to_rot); 
            std::cout << "configuration contains overlaps. move rejected..." << std::endl;
            // std::cout << "this location is occupied..." << std::endl;
        }
    //} //

    return; 

}

*/

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


/*
void Grid::EndRotation(int index){

    int r = (rand()%2);
    std::cout <<"random number r is " << r <<"!"<< std::endl;


    if (r==0){
        std::cout << "Zero Index is being rotated!" << std::endl;
        this->ZeroIndexRotation(index); 
    }    

    else{
        std::cout << "Final index is being rotated!" << std::endl;
        this->FinalIndexRotation(index);
    }

    return;

}
*/


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



/*
void Grid::EndRotation_MC(int index){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        this->ZeroIndexRotation_MC(index); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        this->FinalIndexRotation_MC(index); 

    }

    return; 

}
*/

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// START OF KINK JUMP
/*
void Grid::KinkJump(int index) {

    std::vector <int> k_idx = this->PolymersInGrid.at(index).findKinks(); 

    if (k_idx.size() == 0){
        std::cout << "No kinks in polymer..." << std::endl;
        return;
    }

    std::shuffle(std::begin(k_idx), std::end(k_idx), std::default_random_engine() ); 

    for (int idx: k_idx){
        std::cout << "idx right before kink spot is " << idx << std::endl;

        std::vector <int> d1 = subtract_vectors (&(this->PolymersInGrid.at(index).chain.at(idx+1).coords), &(this->PolymersInGrid.at(index).chain.at(idx).coords) );
        std::vector <int> d2 = subtract_vectors (& (this->PolymersInGrid.at(index).chain.at(idx+2).coords), &(this->PolymersInGrid.at(index).chain.at(idx+1).coords) ); 

        std::vector <int> to_check = add_vectors (& (this->PolymersInGrid.at(index).chain.at(idx).coords), &d2); 

        if (this->OccupancyMap[to_check]==0){
            // update occupancy map 
            this->OccupancyMap[to_check] = 1; 
            this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(idx+1).coords] = 0; 

            // update positions in polymer, and connectivity map 
            this->PolymersInGrid.at(index).chain.at(idx+1).coords = to_check; 
            this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
            return; 

        }

        else {
            std::cout << "this location is occupied..." << std::endl;
            return; 
        }

    }

}
*/


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

/*
void Grid::KinkJump_MC(int index) {

    std::vector <int> k_idx = this->PolymersInGrid.at(index).findKinks(); 

    if (k_idx.size() == 0 ){
        std::cout << "No kinks in the polymer..." << std::endl;
        return; 
    }

    Grid G_temp (this->x, this->y, this->z, this->kT, this->Emm, this->Ems, this->Ess); 
    G_temp.PolymersInGrid = this->PolymersInGrid;
    G_temp.OccupancyMap = this->OccupancyMap; 

    // get energy of current system 
    double Ei = this->Energy; 
    std::cout << "Initial energy is " << Ei << std::endl; 

    std::shuffle(std::begin(k_idx), std::end(k_idx), std::default_random_engine() ); 

    int choice = rng_uniform(0, static_cast<int>(k_idx.size()-1));
    int idx = k_idx.at(choice); 

    // for (int idx: k_idx){
        // std::cout << "idx right before kink spot is " << idx << std::endl;

        std::vector <int> d1 = subtract_vectors (&(this->PolymersInGrid.at(index).chain.at(idx+1).coords), &(this->PolymersInGrid.at(index).chain.at(idx).coords) );
        std::vector <int> d2 = subtract_vectors (& (this->PolymersInGrid.at(index).chain.at(idx+2).coords), &(this->PolymersInGrid.at(index).chain.at(idx+1).coords) ); 
        std::vector <int> to_check = add_vectors (& (this->PolymersInGrid.at(index).chain.at(idx).coords), &d2); 
        impose_pbc(&to_check, this->x, this->y, this->z); 
        if (this->OccupancyMap[to_check]==0){

            G_temp.OccupancyMap[to_check] = 1; 
            G_temp.OccupancyMap[G_temp.PolymersInGrid.at(index).chain.at(idx+1).coords] = 0; 
            G_temp.PolymersInGrid.at(index).chain.at(idx+1).coords = to_check; 
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // get the new energy 
            G_temp.CalculateEnergy(); 
            double Ef = G_temp.Energy; 
            // std::cout << "Coordinates of the polymer are: " << std::endl;
            // G_temp.PolymersInGrid.at(index).printChainCoords();

            std::cout << "Test config energy for kink is " << Ef << std::endl;

            if (acceptance(Ef-Ei, this->kT)){
                // update OccupancyMap
                this->OccupancyMap[to_check] = 1; 
                this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(idx+1).coords] = 0; 

                // update position in polymer 
                this->PolymersInGrid.at(index).chain.at(idx+1).coords = to_check; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
                this->CalculateEnergy(); 
                std::cout << "monte carlo said yes!" << std::endl;
                return;  
            }

            else {
                std::cout << "monte carlo said no..." << std::endl;
                return;
            }



        }

        else {
            // print(to_check);
            std::cout << "this location is occupied..." << std::endl;
        }


    //}

    return; 

}

*/ 






//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
/*
void Grid::CrankShaft(int index){

    std::vector <int> c_idx = this->PolymersInGrid.at(index).findCranks();

    std::shuffle (std::begin (c_idx), std::end(c_idx), std::default_random_engine() ); 

    for (int idx: c_idx) { 
        std::vector <int> d1 = subtract_vectors (&(this->PolymersInGrid.at(index).chain.at(idx+3).coords), &(this->PolymersInGrid.at(index).chain.at(idx+2).coords));
        // impose_pbc(&d1, this->x, this->y, this->z); 

        std::vector <int> to_check_1 = add_vectors(&(this->PolymersInGrid.at(index).chain.at(idx).coords), &d1);
        std::vector <int> to_check_2 = add_vectors(&(this->PolymersInGrid.at(index).chain.at(idx+3).coords), &d1);
        impose_pbc(&to_check_1, this->x, this->y, this->z); 
        impose_pbc(&to_check_2, this->x, this->y, this->z); 
        if ( (this->OccupancyMap[to_check_1]==0) && (this->OccupancyMap[to_check_2]==0) ){
            // update idx + 1
            this->OccupancyMap[to_check_1] = 1; 
            this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(idx+1).coords] = 0;
            this->PolymersInGrid.at(index).chain.at(idx+1).coords = to_check_1; 
            this->PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // update idx + 2 
            this->OccupancyMap[to_check_2] = 1; 
            this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(idx+2).coords] = 0;
            this->PolymersInGrid.at(index).chain.at(idx+2).coords = to_check_2; 
            this->PolymersInGrid.at(index).ChainToConnectivityMap();

            return; 
        }

        else {
            std::cout << "could not perform the crank shaft..." << std::endl;
        }

    }


    return;
}
*/


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
/*
void Grid::CrankShaft_MC(int index){

    double Ei = this->Energy; 
    std::cout << "initial energy is " << Ei << std::endl;
    
    std::vector <int> c_idx = this->PolymersInGrid.at(index).findCranks(); 
    if ( c_idx.size()==0 ){
        std::cout << "No cranks in this polymer..." << std::endl; 
        return; 
    }

    std::shuffle (std::begin (c_idx), std::end(c_idx), std::default_random_engine() ); 

    Grid G_temp (this->x, this->y, this->z, this->kT, this->Emm, this->Ems, this->Ess); 
    G_temp.PolymersInGrid = this->PolymersInGrid; 
    G_temp.OccupancyMap = this->OccupancyMap; 

    // get energy of the current system 

    int choice_idx = rng_uniform(0, static_cast<int>(c_idx.size()-1)); 
    int idx = c_idx.at(choice_idx); 

    // for (int idx: c_idx){//

        // there are three ways to crank shaft 

        // complete 180, positive 90, negative 90 

        // first I need to get the currated vector of directions 
        std::vector <int> HingeToKink = subtract_vectors (&(this->PolymersInGrid.at(index).chain.at(idx+2).coords), &(this->PolymersInGrid.at(index).chain.at(idx+3).coords));
        std::vector <int> HingeToHinge = subtract_vectors(&(this->PolymersInGrid.at(index).chain.at(idx+3).coords), &(this->PolymersInGrid.at(index).chain.at(idx).coords) );

        std::vector <std::vector <int>> drns = HingeSwingDirections(&(HingeToHinge), &(HingeToKink), this->x, this->y, this->z); 

        // pick a direction to swing in 

        int choice = rng_uniform(0, 2); 
        // std::cout << "choice is " << choice << std::endl; 
        
        std::vector <int> d1 = drns.at(choice); 
        // print(d1);

        std::vector <int> to_check_1 = add_vectors (&(this->PolymersInGrid.at(index).chain.at(idx).coords), &d1);
        std::vector <int> to_check_2 = add_vectors (&(this->PolymersInGrid.at(index).chain.at(idx+3).coords), &d1);
        impose_pbc(&to_check_1, this->x, this->y, this->z); 
        impose_pbc(&to_check_2, this->x, this->y, this->z); 

        if ( (this->OccupancyMap[to_check_1]==0) && (this->OccupancyMap[to_check_2]==0) ){
            // update idx+1 
            G_temp.OccupancyMap[to_check_1] = 1; 
            G_temp.OccupancyMap[G_temp.PolymersInGrid.at(index).chain.at(idx+1).coords] = 0; 
            G_temp.PolymersInGrid.at(index).chain.at(idx+1).coords = to_check_1; 
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // update idx+2 
            G_temp.OccupancyMap[to_check_2] = 1; 
            G_temp.OccupancyMap[G_temp.PolymersInGrid.at(index).chain.at(idx+2).coords] = 0; 
            G_temp.PolymersInGrid.at(index).chain.at(idx+2).coords = to_check_2; 
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // get the new energy 
            G_temp.CalculateEnergy(); 
            double Ef = G_temp.Energy; 
            std::cout << "Test config energy for crank is " << Ef << std::endl;

            if (acceptance(Ef-Ei, this->kT)){
                // update OccupancyMap idx+1
                this->OccupancyMap[to_check_1] = 1; 
                this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(idx+1).coords] = 0;
                this->PolymersInGrid.at(index).chain.at(idx+1).coords = to_check_1; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 

                // update OccupancyMap idx+2
                this->OccupancyMap[to_check_2] = 1; 
                this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(idx+2).coords] = 0; 
                this->PolymersInGrid.at(index).chain.at(idx+2).coords = to_check_2; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
                this->CalculateEnergy(); 
                std::cout << "monte carlo said yes!" << std::endl; 
                return; 

            }
            else {
                std::cout << "monte carlo said no..." << std::endl;
                return;
            }



        }
        else {
            std::cout << "this location is occupied..." << std::endl;
        }



    //}    //





    return; 
}
*/
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*
void Grid::ZeroToFinalReptation(int index){

    // get size of polymer chain 
    int size = this->PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->PolymersInGrid.at(index).chain.at(size-1).coords, this->x, this->y, this->z );

    // get rid of second to last monomer position from ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() );

    std::map <std::vector <int>, int> OCmap_temp = this->OccupancyMap; 
    OCmap_temp [this->PolymersInGrid.at(index).chain.at(0).coords] = 0; 

    for (std::vector <int> v: ne_list){
        if (OCmap_temp[v]==0){
            this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(0).coords] = 0;
            this->PolymersInGrid.at(index).chain.erase(this->PolymersInGrid.at(index).chain.begin()); 
            Particle p (v); 
            this->PolymersInGrid.at(index).chain.push_back(p) ; 
            this->OccupancyMap[v]=1;
            this->PolymersInGrid.at(index).ChainToConnectivityMap(); 

            return; 
        }
        else {
            std::cout << "no space to slither..." << std::endl;
        }
    }
    return;


}
*/
/*
void Grid::ZeroToFinalReptation_MC(int index){


    Grid G_temp (this->x, this->y, this->z, this->kT, this->Emm, this->Ems, this->Ess);
    G_temp.PolymersInGrid = this->PolymersInGrid; 
    G_temp.OccupancyMap = this->OccupancyMap; 

    this->CalculateEnergy(); 
    double Ei = this->Energy; 
    std::cout << "Initial energy is " << Ei << std::endl;

    // get size of polymer chain 
    int size = this->PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->PolymersInGrid.at(index).chain.at(size-1).coords, this->x, this->y, this->z );

    // get rid of second to last monomer position from ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() );

    int choice = rng_uniform(0, static_cast<int>(ne_list.size()-1));
    std::vector <int> v = ne_list.at(choice); 


    // for (std::vector <int> v: ne_list){ //
        if (this->OccupancyMap[v] == 0){

            G_temp.OccupancyMap[G_temp.PolymersInGrid.at(index).chain.at(0).coords] = 0;
            G_temp.PolymersInGrid.at(index).chain.erase(G_temp.PolymersInGrid.at(index).chain.begin()); 
            Particle p (v); 
            G_temp.PolymersInGrid.at(index).chain.push_back(p);
            G_temp.OccupancyMap[v]=1; 
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 
            G_temp.CalculateEnergy(); 

            double Ef = G_temp.Energy; 
            std::cout << "Test config energy is " << Ef << std::endl;

            if (acceptance(Ef-Ei, this->kT)){

                this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(0).coords] = 0; 
                this->PolymersInGrid.at(index).chain.erase(this->PolymersInGrid.at(index).chain.begin()); 
                this->PolymersInGrid.at(index).chain.push_back(p); 
                this->OccupancyMap[v] = 1; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
                this->CalculateEnergy(); 

                std::cout << "monte carlo said yes!" << std::endl; 
                return;


            }

            else {
                std::cout << "monte carlo said no..." << std::endl;
                return; 
            }


        }

        else {
            // print(v); 
            std::cout << "this location is occupied..." << std::endl;
        }
    //} // 

    return;

}
*/

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*
void Grid::FinalToZeroReptation(int index){

    // get size of the polymer chain 
    int size = this->PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->PolymersInGrid.at(index).chain.at(0).coords, this->x, this->y, this-> z);

    // get rid of second to last monomer position from ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->PolymersInGrid.at(index).chain.at(1).coords), ne_list.end() ); 

    std::map <std::vector <int>, int> OCmap_temp = this->OccupancyMap; 
    OCmap_temp[this->PolymersInGrid.at(index).chain.at(size-1).coords] = 0; 

    for (std::vector <int> v: ne_list){
        if (OCmap_temp[v]==0){
            this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(size-1).coords] = 0; 
            this->PolymersInGrid.at(index).chain.pop_back(); 
            Particle p (v); 
            this->PolymersInGrid.at(index).chain.insert(this->PolymersInGrid.at(index).chain.begin(), p);
            this->OccupancyMap[v] = 1; 
            this->PolymersInGrid.at(index).ChainToConnectivityMap(); 


            return;  
        }
        else {
            std::cout << "no spack to slither..." << std::endl;
        }
    }


    return; 
}

*/
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/
/*
void Grid::FinalToZeroReptation_MC(int index){

    Grid G_temp (this->x, this->y, this->z, this->kT, this->Emm, this->Ems, this->Ess); 
    G_temp.PolymersInGrid = this->PolymersInGrid; 
    G_temp.OccupancyMap = this->OccupancyMap; 

    this->CalculateEnergy(); 
    double Ei = this->Energy; 
    std::cout << "Initial energy is " << Ei << std::endl;

    // get size of the polymer chain 
    int size = this->PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->PolymersInGrid.at(index).chain.at(0).coords, this->x, this->y, this->z); 

    // get rid of second to last monomer position from ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->PolymersInGrid.at(index).chain.at(1).coords), ne_list.end() ); 

    int choice = rng_uniform(0, static_cast<int>(ne_list.size()-1));
    std::vector <int> v = ne_list.at(choice); 

    //for (std::vector <int> v: ne_list){ //
        if (this->OccupancyMap[v]==0){

            G_temp.OccupancyMap[G_temp.PolymersInGrid.at(index).chain.at(size-1).coords] = 0; 
            G_temp.PolymersInGrid.at(index).chain.pop_back(); 
            Particle p (v); 
            G_temp.PolymersInGrid.at(index).chain.insert(G_temp.PolymersInGrid.at(index).chain.begin(), p); 
            G_temp.OccupancyMap[v] = 1; 
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 
            G_temp.CalculateEnergy(); 

            double Ef = G_temp.Energy; 
            std::cout << "Test config energy is " << Ef << std::endl;

            if (acceptance(Ef-Ei, this->kT)){

                this->OccupancyMap[this->PolymersInGrid.at(index).chain.at(size-1).coords] = 0; 
                this->PolymersInGrid.at(index).chain.pop_back(); 
                this->PolymersInGrid.at(index).chain.insert(this->PolymersInGrid.at(index).chain.begin(), p); 
                this->OccupancyMap[v] = 1; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
                this->CalculateEnergy(); 

                std::cout << "monte carlo said yes!" << std::endl;

                return; 
            }

            else {
                std::cout << "monte carlo said no..." << std::endl;
                return; 
            }

        }

        else {
            // print(v); 
            std::cout << "this location is occupied... " << std::endl;
        }
    // }

    return; 
}



void Grid::Reptation_MC(int index){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        std::cout << "Zero-To-Final slither!" << std::endl;
        this->ZeroToFinalReptation_MC(index); 
    }
    else {
        std::cout << "Final-To-Zero slither!" << std::endl;
        this->FinalToZeroReptation_MC(index); 

    }

    return; 


}

*/


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

/*
void Grid::MonteCarloExecuter(int move_index, int polymer_index){

    if (move_index == 1){
        std::cout << "Running EndRotation_MC on Polymer #" << polymer_index << std::endl;
        this->EndRotation_MC(polymer_index);
        return;
    }
    else if (move_index == 2){
        std::cout << "Running Reptation_MC on Polymer #" << polymer_index << std::endl;
        this->Reptation_MC(polymer_index); 
        return; 
    }
    else if (move_index == 3){
        std::cout << "Running KinkJump_MC on Polymer #" << polymer_index << std::endl;
        this->KinkJump_MC(polymer_index); 
        return;
    }
    else if (move_index == 4){
        std::cout << "Running CrankShaft_MC on Polymer #" << polymer_index << std::endl;
        this->CrankShaft_MC(polymer_index);
        return; 
    }

    else {
        std::cerr << "ERROR: bad move_index provided. it has to be between 1 and 4." << std::endl;
        exit(EXIT_FAILURE); 
    }

    return;


}


void Grid::TheElementaryGridEvolver(){

    // pick a polymer for a monte carlo move 
    int polymer_idx = rng_uniform(0, static_cast<int>(this->PolymersInGrid.size()-1));  
    int move_idx = rng_uniform(1, 4); // look at MonteCarloExecuter for the mapping 

    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#" << std::endl << std::endl;

    this->MonteCarloExecuter(move_idx, polymer_idx); 

    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#" << std::endl<< std::endl;

    return; 

}

*/




/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// Methods for class polymer 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void Polymer::printChainCoords(){
    for (Particle p: this->chain){
        p.printCoords(); 
    }
    return; 
}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void Polymer::printOrientation(){
    for (Particle p:this->chain){
        std::cout << p.orientation << " | ";
    }
    std::cout << std::endl;
    return;
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

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
};


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


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
        std::cout << "there are no cranks in the current structure..." << std::endl;
    } 

    return crank_indices;
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// Methods for class particle

void Particle::printCoords(){
    print(this->coords); 
    return; 
}



/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 




/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// MONTE CARLO MOVES 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


Grid IsingFlip(Grid InitialG){

    std::vector <Particle> Particles = InitialG.ClusterParticleMaker(); 
    std::vector <Particle> an_empty_vector_for_ClusterMaker; 
    std::vector <Particle> cluster = InitialG.ClusterMaker(Particles, an_empty_vector_for_ClusterMaker, 0);

    ClusterFlip (&cluster);                       // flip the cluster! 

    Grid NewG (InitialG);                         // make a copy of the original Grid 

    for (Particle P: cluster){
        NewG.OccupancyMap[P.coords] = P;          // update the OccupancyMap of NewG
    }


    // update PolymersInGrid 
    for (Polymer pmer: NewG.PolymersInGrid){
        for (Particle p: pmer.chain){
            p = NewG.OccupancyMap[p.coords]; 
        }
    }

    // update SolventInGrid 
    for (Particle p: NewG.SolventInGrid){
        p = NewG.OccupancyMap[p.coords]; 
    }


    return NewG; 


}


//#~#~#~#~#~#~~#~#~###~##~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

Grid ZeroIndexRotation(Grid InitialG, int index){

    Grid NewG (InitialG);

    // get the neighborlist of particle at index 1 
    std::vector <int> loc_0 = InitialG.PolymersInGrid.at(index).chain.at(0).coords; 
    std::vector <int> loc_1 = InitialG.PolymersInGrid.at(index).chain.at(1).coords;
    std::vector <int> loc_2 = InitialG.PolymersInGrid.at(index).chain.at(2).coords; 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(loc_1, InitialG.x, InitialG.y, InitialG.z); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0 ), ne_list.end() ); 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2 ), ne_list.end() ); 

    // find a location that is not occupied by monomer 
    // use the occupancy map 

    for (std::vector <int> to_rot: ne_list){

        if (InitialG.OccupancyMap[to_rot].ptype=="solvent"){

            // update OccupancyMap 
            Particle p1 (NewG.OccupancyMap[loc_0]);
            p1.coords = to_rot;  
            Particle p2 (NewG.OccupancyMap[to_rot]); 
            p2.coords = loc_0; 

            NewG.OccupancyMap[to_rot] = p1;  
            NewG.OccupancyMap[loc_0] = p2; 

            // update positions in polymers in grid 

            NewG.PolymersInGrid.at(index).chain.at(0).coords = to_rot; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // update position of the solvent particle you displaced during this move 
            int c =0; 
            for (Particle P: NewG.SolventInGrid){
                c++; 
                if (P.coords == to_rot){ 
                     
                    break;
                }

            }
            NewG.SolventInGrid.at(c-1).coords = p2.coords; 
            std::cout << "coordinates of p2 are ";
            print(p2.coords);
            break;
            // std::cout << "Something is up, line 1700: classes.cpp" << std::endl;

        }


    }

    return NewG;
}

//#~#~#~#~#~#~~#~#~###~##~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

Grid FinalIndexRotation(Grid InitialG, int index){

    Grid NewG (InitialG);

    int dop = InitialG.PolymersInGrid.at(index).chain.size(); 

    // get the neighborlist of particle at index DoP-2
    std::vector <int> loc_0 = InitialG.PolymersInGrid.at(index).chain.at(dop-1).coords;
    std::vector <int> loc_1 = InitialG.PolymersInGrid.at(index).chain.at(dop-2).coords;
    std::vector <int> loc_2 = InitialG.PolymersInGrid.at(index).chain.at(dop-3).coords;

    // obtain neighbor list 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(loc_1, InitialG.x, InitialG.y, InitialG.z); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_0), ne_list.end() ); 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), loc_2), ne_list.end() );

    // find a location that is unoccupied by monomer 
    // use the occupancy map  

    for (std::vector <int> to_rot: ne_list){

        if (NewG.OccupancyMap[to_rot].ptype == "solvent") {

            // update OccupancyMap 
            Particle p1 (NewG.OccupancyMap[loc_0]); 
            Particle p2(NewG.OccupancyMap[to_rot]); 
            p1.coords = to_rot; 
            p2.coords = loc_0; 

            NewG.OccupancyMap[to_rot] = p1; 
            NewG.OccupancyMap[loc_0] = p2;

            // update positions in polymers in grid 
            NewG.PolymersInGrid.at(index).chain.at(dop-1).coords = to_rot; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // update position of the solvent particle you displaced during this move
            for (Particle P: NewG.SolventInGrid){
                if (P.coords == to_rot){
                    P.coords = loc_0; 
                    break;
                }

            }
            break;
            // std::cout << "Something is up, line 1700: classes.cpp" << std::endl;


        }
    }


    return NewG;
}

//#~#~#~#~#~#~~#~#~###~##~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


Grid EndRotation(Grid InitialG, int index){

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


//#~#~#~#~#~#~~#~#~###~##~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


Grid KinkJump(Grid InitialG, int index){

    Grid NewG(InitialG); 

    std::vector <int> k_idx = InitialG.PolymersInGrid.at(index).findKinks(); 

    if (k_idx.size() == 0 ){
        std::cout << "No kinks found in polymer..." << std::endl;
        return InitialG;
    }

    std::shuffle(std::begin(k_idx), std::end(k_idx), std::default_random_engine() ); 

    for (int idx: k_idx){

        std::cout << "idx right before kink spot is " << idx << std::endl; 

        std::vector <int> d1 = subtract_vectors(&(InitialG.PolymersInGrid.at(index).chain.at(idx+1).coords), &(InitialG.PolymersInGrid.at(index).chain.at(idx).coords) );
        std::vector <int> d2 = subtract_vectors(&(InitialG.PolymersInGrid.at(index).chain.at(idx+2).coords), &(InitialG.PolymersInGrid.at(index).chain.at(idx+1).coords) ); 

        std::vector <int> to_check = add_vectors( &(InitialG.PolymersInGrid.at(index).chain.at(idx).coords), &d2); 

        if (NewG.OccupancyMap[to_check].ptype=="solvent"){

            //update occupancy map
            Particle p1 ( NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+1).coords] );  // this will be the monomer particle
            Particle p2 ( NewG.OccupancyMap[to_check] ) ;                                             // this will be the solvent particle
            p1.coords = to_check; 
            p2.coords = NewG.PolymersInGrid.at(index).chain.at(idx+1).coords;

            NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+1).coords] = p2;
            NewG.OccupancyMap[to_check] = p1; 

            // update PolymersInGrid 
            NewG.PolymersInGrid.at(index).chain.at(idx+1) = p1; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            // update position of the solvent particle you displaced during this move 
            for (Particle P: NewG.SolventInGrid){
                if (P.coords == to_check){
                    P.coords = InitialG.PolymersInGrid.at(index).chain.at(idx+1).coords; 
                    break;
                }

            }
            break;



        }
        else {
            std::cout << "This spot is taken by a monomer... no kink jumping." << std::endl;
        }



    }

    return NewG; 

}

//#~#~#~#~#~#~~#~#~###~##~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


Grid CrankShaft(Grid InitialG, int index){

    Grid NewG(InitialG); 

    std::vector <int> c_idx = InitialG.PolymersInGrid.at(index).findCranks(); 

    if ( c_idx.size()==0 ){
        std::cout << "No cranks in this polymer..." << std::endl; 
        return InitialG; 
    }

    std::shuffle (std::begin (c_idx), std::end(c_idx), std::default_random_engine() ); 

    for (int idx: c_idx){

        std::vector <int> d1 = subtract_vectors( &(NewG.PolymersInGrid.at(index).chain.at(idx+3).coords), &(NewG.PolymersInGrid.at(index).chain.at(idx+2).coords) );

        std::vector <int> to_check_1 = add_vectors ( &(NewG.PolymersInGrid.at(index).chain.at(idx).coords), &d1 ); 
        std::vector <int> to_check_2 = add_vectors ( &(NewG.PolymersInGrid.at(index).chain.at(idx+3).coords), &d1 ); 
        impose_pbc(&to_check_1, InitialG.x, InitialG.y, InitialG.z); 
        impose_pbc(&to_check_2, InitialG.x, InitialG.y, InitialG.z);  

        if ( (InitialG.OccupancyMap[to_check_1].ptype=="solvent") && (InitialG.OccupancyMap[to_check_2].ptype=="solvent") ){

            // update OccupancyMap
            Particle p1 ( NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+1).coords] );
            Particle s1 ( NewG.OccupancyMap[to_check_1]); 

            p1.coords = to_check_1;
            s1.coords = NewG.PolymersInGrid.at(index).chain.at(idx+1).coords; 

            Particle p2 ( NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+2).coords]); 
            Particle s2 ( NewG.OccupancyMap[to_check_2] );

            p2.coords = to_check_2; 
            s2.coords = NewG.PolymersInGrid.at(index).chain.at(idx+2).coords; 

            NewG.OccupancyMap[to_check_1] = p1; 
            NewG.OccupancyMap[to_check_2] = p2; 

            NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+1).coords] = s1; 
            NewG.OccupancyMap[NewG.PolymersInGrid.at(index).chain.at(idx+2).coords] = s2; 

            // update PolymersInGrid 
            NewG.PolymersInGrid.at(index).chain.at(idx+1) = p1; 
            NewG.PolymersInGrid.at(index).chain.at(idx+2) = p2; 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 


            // update position of the solvent particle you displaced during this move 
            int count = 0; 
            for (Particle P: NewG.SolventInGrid){
                if (P.coords == to_check_1){
                    P.coords = s1.coords;             //InitialG.PolymersInGrid.at(index).chain.at(idx+1).coords; 
                    count ++ ;
                }
                else if (P.coords == to_check_2){
                    P.coords = s2.coords;             //InitialG.PolymersInGrid.at(index).chain.at(idx+2).coords; 
                    count++; 
                }

                else if (count==2){
                    break;
                }

            }
            break;

        }
        else {
            std::cout << "This position is occupied by monomer... No cranking." << std::endl;
        }


    }

    return NewG;

}

//#~#~#~#~#~#~~#~#~###~##~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

Grid ZeroToFinal_Reptation(Grid InitialG, int index){

    Grid NewG (InitialG); 
    // get size of polymer chain 
    int size = InitialG.PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list( InitialG.PolymersInGrid.at(index).chain.at(size-1).coords, InitialG.x, InitialG.y, InitialG.z ); 

    // get rid of second to last monomer position from ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), InitialG.PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() ); 


    for (std::vector <int> v: ne_list){
        if (InitialG.OccupancyMap[v].ptype == "solvent"){
            Particle p = InitialG.OccupancyMap[v];
            p.ptype = "monomer"; 
            
            Particle s = InitialG.OccupancyMap[InitialG.PolymersInGrid.at(index).chain.at(0).coords]; 
            s.ptype = "solvent"; 

            // update occupancy map 
            NewG.OccupancyMap[v] = p;
            NewG.OccupancyMap[InitialG.PolymersInGrid.at(index).chain.at(0).coords] = s; 

            // update polymers in grid 
            NewG.PolymersInGrid.at(index).chain.erase(NewG.PolymersInGrid.at(index).chain.begin() ); 
            NewG.PolymersInGrid.at(index).chain.push_back(p); 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); 

            //update solvents in grid 
            for (Particle P: NewG.SolventInGrid){
                if (P.coords == v){
                    P.coords = s.coords; 
                    break;
                }
            }

            break;

        }
        else {
            std::cout << "No place to slither... " << std::endl;
        }
    }

    return NewG; 

}; 

//#~#~#~#~#~#~~#~#~###~##~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


Grid FinalToZero_Reptation(Grid InitialG, int index){
    
    Grid NewG(InitialG); 

    // get initial size of the polymer chain 
    int size = InitialG.PolymersInGrid.at(index).chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list( InitialG.PolymersInGrid.at(index).chain.at(0).coords, InitialG.x, InitialG.y, InitialG.z); 

    // get rid of second monomer from the ne list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), InitialG.PolymersInGrid.at(index).chain.at(1).coords), ne_list.end() ); 

    for (std::vector <int> v: ne_list){

        if (InitialG.OccupancyMap[v].ptype == "solvent"){
            Particle p = InitialG.OccupancyMap[v]; 
            p.ptype = "monomer"; 

            Particle s = InitialG.OccupancyMap[InitialG.PolymersInGrid.at(index).chain.at(size-1).coords]; 
            s.ptype = "solvent"; 

            // update occupancy map 
            NewG.OccupancyMap[v] = p; 
            NewG.OccupancyMap[InitialG.PolymersInGrid.at(index).chain.at(size-1).coords] = s; 

            // update polymers in grid 
            NewG.PolymersInGrid.at(index).chain.pop_back(); // pop the final particle 
            NewG.PolymersInGrid.at(index).chain.insert(NewG.PolymersInGrid.at(index).chain.begin(), p); // add one particle to the head of the chain 
            NewG.PolymersInGrid.at(index).ChainToConnectivityMap(); // make sure the connectivity map is maintained 

            // update solvents 

            for (Particle P: NewG.SolventInGrid){
                if (P.coords == v){
                    P.coords = s.coords;
                    break; 
                }
            }

            break;
        }

        else {
            std::cout << "No place to slither..." << std::endl;
        }

    }

    return NewG; 
};


Grid Reptation(Grid InitialG, int index){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return FinalToZero_Reptation(InitialG, index); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return ZeroToFinal_Reptation(InitialG, index); 

    }
    
};
