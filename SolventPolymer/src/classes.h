#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_ 





/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

class Particle{
public: 
    std::vector <int> coords;                                        // the coordinates of the particles
    std::string ptype; 
    int orientation; 
    bool operator<(const Particle& rhs)const{
        return coords < rhs.coords; 
    }


    // constructor 
    Particle(){};  // default constructor

    Particle (std::vector <int> crds, std::string type_, int orientation_): coords (crds), ptype (type_), orientation (orientation_){

    }

    // destructor 
    ~Particle(){

    }

    // Particle( const Particle &other); // copy constructor for some reason i dont know

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
This is the Grid Class. 
This is the master class. Everything cool happens to this guy over here. 
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/

class Grid{

public:
    std::vector <Polymer> PolymersInGrid;                 // all the polymers in the grid 
    std::vector <Particle> SolventInGrid;                // solvent molecules in the grid 
    const int x;                                          // length of x-edge of grid 
    const int y;                                          // length of y-edge of grid 
    const int z;                                          // length of z-edge of grid 
    const double kT;                                      // energy factor 
    const double Emm ;                                    // monomer-monomer interaction 
    const double Ess ;                                    // solvent-solvent interaction 
    const double Ems_n ;                                  // monomer-solvent when Not aligned 
    const double Ems_a ;                                  // monomer-solvent when Aligned

    double Energy; 
    std::map <std::vector <int>, Particle> OccupancyMap;     // a map that gives the particle given the location
    

    Grid(int xlen, int ylen, int zlen, double kT_, double Emm_, double Ess_, double Ems_n_, double Ems_a_): x (xlen), y (ylen), z (zlen), kT (kT_), Emm(Emm_), Ess (Ess_), Ems_n (Ems_n_), Ems_a (Ems_a_) {        // Constructor of class
        // this->instantiateOccupancyMap(); 
    };

    ~Grid(){                                    // Destructor of class 

    }; 

    // get the initial completed unoccupied map
    void instantiateOccupancyMap();      

    // update OccupancyMap 
    void updateOccupancyMap();      

    // plant the polymer from input file 
    void plantPolymersInGrid(std::string filename);  

    int ExtractNumberOfPolymers(std::string filename);

    std::vector <std::string> ExtractContentFromFile(std::string filename); 

    void ExtractPolymersFromFile(std::string filename);
    

    void CalculateEnergy(); 

    void dumpPositionsOfPolymers (int step); 


    // check validity of input coords
    bool checkValidityOfCoords(std::vector <int> v);
    bool checkForOverlaps(std::vector <Polymer> PolymerVector); 
    bool checkConnectivity(std::vector <Polymer> PolymerVector); 

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
    void KinkJump_MC(int polymer_index);

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // crank shaft 
    void CrankShaft(int polymer_index); 
    void CrankShaft_MC(int polymer_index); 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    // reptation 
    void ZeroToFinalReptation(int index);
    void ZeroToFinalReptation_MC(int index);

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    void FinalToZeroReptation(int index); 
    void FinalToZeroReptation_MC(int index); 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    void Reptation(int polymer_index); 
    void Reptation_MC(int polymer_index); 



    void MonteCarloExecuter(int move_index, int polymer_index);
    void TheElementaryGridEvolver();
}; 


Grid CreateGridObject(std::string positions, std::string topology);





#endif 