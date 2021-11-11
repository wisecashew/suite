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






//
/*void Grid::plantPolymersInGrid(std::string filename){

    this->PolymersInGrid = ExtractPolymersFromFile(filename, ); 




    for (Polymer pmer: this->PolymersInGrid){
        for (Particle p: pmer.chain){
            if (checkValidityOfCoords(p.coords)){
            this->OccupancyMap[p.coords]=1;
            }
            else {
                std::cout << "Coordinates are out of bounds. Bad input file." << std::endl;
                return;
            }
        }
    }

    return; 

}*/ 





bool Grid::checkValidityOfCoords(std::vector <int> v){
    if (v.at(0)>this->x || v.at(0) < 0 || v.at(1)>this->y || v.at(1)<0 || v.at(2)>this->z || v.at(2)<0){
        return false;
    }
    else {
        return true;
    }
}



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

















/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// Methods for class polymer 

void Polymer::printChainCoords(){
    for (Particle p: this->chain){
        p.printCoords(); 
    }
    return; 
}



void Polymer::ChainToConnectivityMap(){

    const int chainLength = this->chain.size(); 

    for (int i{0}; i<chainLength; i++){
        if (i==0){

            this->ConnectivityMap[this->chain.at(i)] = { (this->chain.at(i+1)) };

        }
        else if (i==(chainLength-1)){

            this->ConnectivityMap[this->chain.at(i)] = { (this->chain.at(i-1))}; 
        }
        else {
            this->ConnectivityMap[this->chain.at(i)] = { (this->chain.at(i-1)), (this->chain.at(i+1)) };
        }

    }

    return; 
};


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












/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// Methods for class particle

void Particle::printCoords(){
    print(this->coords); 
    return; 
}