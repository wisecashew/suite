#ifndef _POLYMER_H_
#define _POLYMER_H_ 
#include "misc.h"


class Polymer{
public:
    short deg_poly; 
    std::vector <Particle*> chain;                                                // all the particles in the polymer chain 
    std::map <Particle*, std::vector <Particle*> > ConnectivityMap;               // the set of beads that one particular polymer bead is connected to 


    // constructor 
    Polymer (short deg_poly_, std::vector <Particle*> particleChain):  deg_poly (deg_poly_), chain (particleChain) {
        // this->ChainToConnectivityMap(); 
        this->chain.reserve(deg_poly_); 
    }


    // destructor 
    ~Polymer (){

    };

    // print positions of monomer present in the polymer
    void print_chain(); 

    // print orientation of monomer present in the polymer
    void print_orientation(); 

    // obtain connectivity map given the chain 
    void chain_to_connectivity_map(); 

};

#endif