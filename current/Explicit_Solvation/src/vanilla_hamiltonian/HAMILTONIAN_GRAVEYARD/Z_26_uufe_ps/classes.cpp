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
#include <array> 
#include "classes.h"
#include "misc.h" 


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
    for (Particle*& p: this->chain){
        p->printCoords(); 
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
    for (Particle*& p:this->chain){
        std::cout << p->orientation << " | ";
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

    for (int i{0}; i<chainLength; ++i){
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
    const int chainLength = this->deg_poly; 
    // std::cout << "chain length is " << chainLength << std::endl;
    std::vector <int> kink_indices; 

    // obtain all location where kinks exist in the polymer 
    for (int i{0}; i<chainLength-2; ++i){
        std::array <int,3> a1 = subtract_arrays(&(this->chain[(i+1)]->coords), &(this->chain[(i)]->coords) ); 
        std::array <int,3> a2 = subtract_arrays(&(this->chain[(i+2)]->coords), &(this->chain[(i+1)]->coords) );

        if (a1==a2){
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
    for (int i{0}; i< chainLength-3; ++i){
        std::array <int,3> a1 = subtract_arrays(&(this->chain[i+1]->coords), &(this->chain[i]->coords) ); 
        std::array <int,3> a2 = subtract_arrays(&(this->chain[i+2]->coords), &(this->chain[i+3]->coords) );

        //

        if (a1==a2){
            crank_indices.push_back(i); 
        } 
        else {
            continue;
        }
    }

    if (crank_indices.size() == 0) {
        // std::cout << "there are no cranks in the current structure..." << std::endl;
    } 
    // printf("Crank indices are: "); 
    // print(crank_indices);
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



