#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_ 


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

class Particle{
public: 
    std::vector <int> coords;                                        // the coordinates of the particles

    bool operator<(const Particle& rhs)const{
        return coords < rhs.coords; 
    }


    // constructor 
    Particle (std::vector <int> crds): coords (crds){

    }

    // destructor 
    ~Particle(){

    }

    // print location of the particle 
    void printCoords(); 

};



/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

class Polymer{
public: 
    std::vector <Particle> chain;                                     // all the particles in the polymer chain
    std::map <Particle, std::vector <Particle>> ConnectivityMap;      // the particles one particular polymer bead is connected to 


    // constructor 
    Polymer (std::vector <Particle> particleChain): chain (particleChain) {
        this->ChainToConnectivityMap(); 
    }


    // destructor 
    ~Polymer (){

    };

    // print positions of monomer present in the polymer
    void printChainCoords(); 

    // obtain connectivity map given the chain 
    void ChainToConnectivityMap(); 

    // find if there are kinks in the polymer structure 
    std::vector <int> findKinks(); 

    // find if there are any cranks in the polymer structure 
    std::vector <int> findCranks(); 


};


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
This is the Grid Class. 
This is the master class. Everything cool happens to this guy over here. 
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/

class Grid{

public:
    std::vector <Polymer> PolymersInGrid;                 // all the polymers in the grid 
    const int x;                                        // length of x-edge of grid 
    const int y;                                        // length of y-edge of grid 
    const int z;                                        // length of z-edge of grid 
    const double mmintrxenergy {-1}; 
    double Energy; 
    std::map <std::vector <int>, int> OccupancyMap;     // a map that checks occupancy of a spot 
    double kT; 

    Grid(int xlen, int ylen, int zlen, double kT_): x (xlen), y (ylen), z (zlen), mmintrxenergy (-1), kT (kT_) {        // Constructor of class
        this->instantiateOccupancyMap(); 
    };



    ~Grid(){                                    // Destructor of class 

    }; 

    // get the initial completed unoccupied map
    void instantiateOccupancyMap();      

    // update OccupancyMap 
    void updateOccupancyMap();      

    // plant the polymer from input file 
    void plantPolymersInGrid(std::string filename);  

    // check validity of input coords
    bool checkValidityOfCoords(std::vector <int> v);


    int ExtractNumberOfPolymers(std::string filename);

    std::vector <std::string> ExtractContentFromFile(std::string filename); 

    void ExtractPolymersFromFile(std::string filename);
    

    void CalculateEnergy(); 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // terminal end rotation 
    void ZeroIndexRotation(int polymer_index);
    void FinalIndexRotation(int polymer_index); 
    void EndRotation(int polymer_index); 

    // include MC criterion 
    void ZeroIndexRotation_MC(int polymer_index);
    void FinalIndexRotation_MC(int polymer_index); 
    void EndRotation_MC(int polymer_index); 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // kink jump
    void KinkJump(int polymer_index); 


}; 













#endif 