#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_ 
#include "misc.h"


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                 DEFINITIONS FOR CLASS PARTICLE 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


class Particle{
public: 
    std::array <int,3> coords;                                        // the coordinates of the particles
    short int orientation;
    std::string ptype; 

    bool operator<(const Particle& rhs)const{
        return (this->coords) < (rhs.coords); 
    }

    bool operator==(const Particle& rhs){
        return std::tie(coords, ptype, orientation) == std::tie(rhs.coords, rhs.ptype, rhs.orientation );
    } 


    // constructor 
    Particle(){};  // default constructor

    Particle (std::array <int, 3> crds, std::string type_, int orientation_): coords (crds), orientation (orientation_), ptype (type_) { };

    // destructor 
    ~Particle(){

    }

    // print location of the particle 
    void printCoords(); 

};

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// =====       END OF CLASS PARTICLE       =====

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                 DEFINITIONS FOR CLASS POLYMER 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


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

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// =====       END OF CLASS POLYMER      =====

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


#endif 
