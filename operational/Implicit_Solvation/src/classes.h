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
    std::string ptype; 
    int orientation;

    bool operator<(const Particle& rhs)const{
        return (this->coords) < (rhs.coords); 
    }

    bool operator==(const Particle& rhs){
        return std::tie(coords, ptype, orientation) == std::tie(rhs.coords, rhs.ptype, rhs.orientation );
    } 


    // constructor 
    Particle(){};  // default constructor

    Particle (std::array <int, 3> crds, std::string type_, int orientation_): coords (crds), ptype (type_), orientation (orientation_){

    }

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
    int deg_poly; 
    std::vector <Particle> chain;                                               // all the particles in the polymer chain
    std::map <Particle, std::vector <Particle> > ConnectivityMap;               // the particles one particular polymer bead is connected to 


    // constructor 
    Polymer (int deg_poly_, std::vector <Particle> particleChain): deg_poly (deg_poly_), chain (particleChain) {
        this->ChainToConnectivityMap(); 
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

// =====       DEFINITIONS FOR CLASS GRID       ===== 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
This is the Grid Class. 
This is the master class. Everything cool happens to this guy over here. 
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/

class Grid{

public:
    int Np;                                                  // number of polymers in the grid 
    std::vector <Polymer> PolymersInGrid;                    // all the polymers in the grid 
    int x;                                                   // length of x-edge of grid 
    int y;                                                   // length of y-edge of grid 
    int z;                                                   // length of z-edge of grid 
    double kT;                                               // energy factor 
    double Emm_n ;                                           // monomer-solvent when Not aligned 
    double Emm_a ;                                           // monomer-solvent when Aligned
    double Ems;                                              // monomer-solvent interaction
    double Energy;                                           // energy of grid 
    std::map <std::array <int,3>, Particle> OccupancyMap;    // a map that gives the particle given the location
    

    Grid() {}; 

    Grid(int Np_, int xlen, int ylen, int zlen, double kT_, double Emm_a_, double Emm_n_, double Ems_): Np(Np_), x (xlen), y (ylen), z (zlen), kT (kT_), Emm_n(Emm_n_), Emm_a (Emm_a_), Ems (Ems_) {        // Constructor of class
        this->PolymersInGrid.reserve(Np_); 
        // this->instantiateOccupancyMap(); 
    };

    /*
    // Destructor of class 
    ~Grid(){                                    

    }; 

    // assignment operator that allows for a correct transfer of properties. Important to functioning of program. 
    Grid& operator=(Grid other){
        std::swap(PolymersInGrid, other.PolymersInGrid); 
        std::swap(Energy, other.Energy); 
        std::swap(OccupancyMap, other.OccupancyMap);
        return *this; 
    } 
    */
    // get the initial completed unoccupied map
    void instantiateOccupancyMap();      

    // update OccupancyMap 
    void updateOccupancyMap();      

    // plant the polymer from input file 
    void plantPolymersInGrid(std::string filename);  

    // extract polymer coordinates from a file 
    void ExtractPolymersFromFile(std::string filename);
    
    // calculate energy of Grid 
    void CalculateEnergy(); 

    // dump coordinates of polymers in Grid and energy of grid into a text file 
    void dumpPositionsOfPolymers (int step, std::string filename="dumpfile.txt"); 
    void dumpEnergyOfGrid (int step, std::string filename, bool first_call);
    void ExtractPolymersFromTraj(std::string trajectory, std::string filename);
    int ExtractIndexOfFinalMove(std::string trajectory);

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // a function to calculate energy of interaction between two particles. 

    // double EnergyPredictor(Particle p1, Particle p2){
    //    if (p1.orientation == p2.orientation){
    //        return this->Emm_a, 
    //    }
    // };
    // 
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    // check validity of input coords
    bool checkValidityOfCoords(std::array <int,3> v);
    bool checkForOverlaps(std::vector <Polymer> PolymerVector); 
    bool checkConnectivity(std::vector <Polymer> PolymerVector); 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    // create clusters meant for Ising flips
    std::vector <Particle> ClusterParticleMaker();
    std::vector <Particle> ClusterMaker(std::vector <Particle> Particles, std::vector <Particle> final, std::vector <Particle> to_send_, int count=0); 
    

}; 


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// =====       END OF CLASS GRID      =====

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 



/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// =====       DEFINITIONS OF MONTE CARLO MOVES       ===== 

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

Grid CreateGridObject(std::string positions, std::string topology);
Grid CreateGridObjectRestart(std::string positions, std::string topology, std::string trajectory);
Grid IsingFlip(Grid InitialG);
Grid ZeroIndexRotation(Grid* InitialG, int index, bool*);
Grid FinalIndexRotation(Grid* InitialG, int index, bool*);
Grid EndRotation(Grid* InitialG, int index, bool*); 
Grid KinkJump(Grid* InitialG, int index, bool*); 
Grid CrankShaft(Grid* InitialG, int index, bool*); 
Grid ForwardReptation(Grid* InitialG, int index, bool*); 
Grid BackwardReptation(Grid* InitialG, int index, bool*);
Grid Reptation(Grid* InitialG, int index, bool*);
Grid Translation(Grid* InitialG, int index, std::vector <int> direction);  
Grid MoveChooser(Grid* InitialG, bool v, bool* IMP_BOOL); 
Grid ZeroIndexRotationAgg(Grid* InitialG, int index, bool* IMP_BOOL);
Grid FinalIndexRotationAgg(Grid* InitialG, int index, bool* IMP_BOOL);
int ExtractNumberOfPolymers(std::string filename);

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// =====       END OF MONTE CARLO MOVES      =====

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// extract topology from the topology file 
std::vector <std::string> ExtractContentFromFile(std::string filename); 


#endif 
