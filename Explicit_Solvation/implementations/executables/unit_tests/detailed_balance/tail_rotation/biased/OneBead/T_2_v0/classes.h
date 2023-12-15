#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_ 
#include <array>


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
    void printChainCoords(); 

    // print orientation of monomer present in the polymer
    void printOrientation(); 

    // obtain connectivity map given the chain 
    void ChainToConnectivityMap(); 

    // find if there are kinks in the polymer structure 
    std::vector <int> findKinks(); 

    // find if there are any cranks in the polymer structure 
    std::vector <int> findCranks(); 


};

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// =====       END OF CLASS POLYMER      =====

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

int ExtractNumberOfPolymers(std::string filename);
std::vector <std::string> ExtractContentFromFile(std::string filename); 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// =====       END OF MONTE CARLO MOVES      =====

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// extract topology from the topology file 
/*
class MonomerSwing {
public:
    
    std::array <int,3> old_loc = {-1, -1, -1};
    std::array <int,3> new_loc = {-1, -1, -1}; 

    // constructor 
    MonomerSwing () {};
    MonomerSwing (std::array <int,3> ol, std::array <int,3> nl): old_loc (ol), new_loc (nl); 

    // destructor 
    ~MonomerSwing() {};
};


class OrientationFlip {
public:

    std::array <int,26> locations; 
    std::array <int,26> old_orientations; 
    std::array <int,26> new_orientations; 

    // constructor 
    OrientationFlip() {};

    // destructor 
    ~OrientationFlip () {}; 

};
*/
#endif 
