#include <iostream>
#include "Polymer.hpp"


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                    POLYMER METHODS. 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//============================================================
//============================================================
// 
// NAME OF FUNCTION: print_chain
//
// PARAMETERS: none, method on a Polymer
// 
// WHAT THE FUNCTION DOES: Given a Polymer, it will print out coordinates of the polymers in a nice way 
//
// DEPENDENCIES: print
//
// THE CODE: 


void Polymer::print_chain(){
	for (Particle*& p: this->chain){
		p->print_coords(); 
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
// NAME OF FUNCTION: print_orientation
//
// PARAMETERS: none, method on a Polymer
// 
// WHAT THE FUNCTION DOES: Given a Polymer, it will print out Ising orientation of monomer units
//
// DEPENDENCIES: print
//
// THE CODE: 


void Polymer::print_orientation(){
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


void Polymer::chain_to_connectivity_map(){

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
