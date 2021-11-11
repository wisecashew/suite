#include <iostream> 
#include <vector>
#include <string> 
#include <map>
#include <algorithm> 
#include <regex>
#include <fstream>
#include <sstream>
#include <regex>
#include "classes.h"
#include "misc.h" 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// Methods for class grid 

void Grid::instantiateOccupancyMap(){
    for (int i{0}; i<this->x; i++){
        for (int j{0}; j<this->y; j++){
            for (int k{0}; k<this->z; k++){
                std::vector <int> lat_pt = {i, j, k}; 
                this->OccupancyMap[lat_pt] = 0;
            }
        }
    }
    return; 
}



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


int Grid::ExtractNumberOfPolymers(std::string filename){
    std::regex start("START"); 
    std::regex end ("END"); 
    std::ifstream myfile (filename); 

    std::string myString; 
    std::vector <std::string> StartContents; 
    std::vector <std::string> EndContents; 

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
    std::vector <Particle> ParticleVector;

    std::vector <std::string> contents = this->ExtractContentFromFile(filename); 

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
        }

        else{
            std::vector <int> loc; 
            for (int i=0; ss>>i; ){

                loc.push_back(i);
            }

        if (!this->checkValidityOfCoords(loc)){
            std::cerr << "Coordinates are out of bounds. Bad input file." << std::endl;
            return;
        }
        Particle p (loc); 
        // print(loc);
        this->OccupancyMap[loc] = 1; 

        ParticleVector.push_back(p); 
        
        } 

        if (startCount > endCount){
            
            continue; 
        }

        else{
            Polymer pmer (ParticleVector);

            PolymerVector.push_back(pmer);
            
            ParticleVector.clear(); 
            
        }



    }
    

    this->PolymersInGrid = PolymerVector;
    return; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#




void Grid::CalculateEnergy(){
    double Energy {0.0}; 
    for (Polymer pmer: this->PolymersInGrid){
        for (Particle p: pmer.chain){

            std::vector <Particle> part_vec = pmer.ConnectivityMap[p];

            std::vector <std::vector <int>> ne_list = obtain_ne_list(p.coords, this->x, this->y, this->z); // get neighbor list 

            // curate neighbor list - get rid of neighbors you are connected to 
            
            for (Particle ple: part_vec){
                ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), ple.coords), ne_list.end());
            }


            int nNeighbors {0}; 
            for (std::vector <int> nei: ne_list){
                
                nNeighbors += this->OccupancyMap[nei]; 

            }
            Energy += nNeighbors*(this->mmintrxenergy)*0.5; 

        }
    }

    this->Energy = Energy; 
    return; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



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

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void Grid::ZeroIndexRotation_MC(int index){
    // create a dummy grid object 

    Grid G_temp (this->x, this->y, this->z, this->kT);
    G_temp.PolymersInGrid = this->PolymersInGrid;
    G_temp.OccupancyMap = this->OccupancyMap; 

    // get energy of current system
    double Ei = this->Energy; 

    
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
            // make the copy 
            G_temp.PolymersInGrid.at(index).chain.at(0).coords=to_rot;
            G_temp.OccupancyMap[to_rot] = 1; 
            G_temp.OccupancyMap[loc_0] = 0;
            G_temp.PolymersInGrid.at(index).ChainToConnectivityMap(); 
            // get the new energy 
            G_temp.CalculateEnergy(); 
            double Ef = G_temp.Energy; 

            if (acceptance (Ef-Ei, this->kT)){
                // update OccupancyMap 
                this->OccupancyMap[to_rot] = 1; // particle goes into new location 
                this->OccupancyMap[loc_0] = 0; 

                // update position in polymer 
                this->PolymersInGrid.at(index).chain.at(0).coords = to_rot; 
                this->PolymersInGrid.at(index).ChainToConnectivityMap(); 
                return;
            }
            else {
                    G_temp.PolymersInGrid = this->PolymersInGrid;
                    G_temp.OccupancyMap = this->OccupancyMap; 
            }
        }
        else {
            std::cout << "this location is occupied..." << std::endl;
        }

    }

    return; 
}    







//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


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

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

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

    if (kink_indices.size() == 0){
        std::cout << "No kinks in polymer..." << std::endl;
    }

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