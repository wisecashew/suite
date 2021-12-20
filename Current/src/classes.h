#ifndef _MC_CLASSES_H_
#define _MC_CLASSES_H_ 




class Particle{
public: 
    std::vector <int> coords;                                        // the coordinates of the particles
    std::string ptype; 
    int orientation;

    bool operator<(const Particle& rhs)const{
        return coords < rhs.coords; 
    }

    bool operator==(const Particle& rhs){
        return std::tie(coords, ptype, orientation) == std::tie(rhs.coords, rhs.ptype, rhs.orientation);
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
    std::vector <Particle> SolventInGrid;                 // solvent molecules in the grid 
    int x;                                          // length of x-edge of grid 
    int y;                                          // length of y-edge of grid 
    int z;                                          // length of z-edge of grid 
    double kT;                                      // energy factor 
    double Emm ;                                    // monomer-monomer interaction 
    double Ess ;                                    // solvent-solvent interaction 
    double Ems_n ;                                  // monomer-solvent when Not aligned 
    double Ems_a ;                                  // monomer-solvent when Aligned

    double Energy; 
    std::map <std::vector <int>, Particle> OccupancyMap;     // a map that gives the particle given the location
    

    Grid(int xlen, int ylen, int zlen, double kT_, double Emm_, double Ess_, double Ems_n_, double Ems_a_): x (xlen), y (ylen), z (zlen), kT (kT_), Emm(Emm_), Ess (Ess_), Ems_n (Ems_n_), Ems_a (Ems_a_) {        // Constructor of class
        // this->instantiateOccupancyMap(); 
    };


    // Destructor of class 
    ~Grid(){                                    

    }; 

    // assignment operator that allows for a correct transfer of properties. Important to functioning of program. 
    Grid& operator=(Grid other){
        std::swap(PolymersInGrid, other.PolymersInGrid); 
        std::swap(SolventInGrid, other.SolventInGrid); 
        std::swap(Energy, other.Energy); 
        std::swap(OccupancyMap, other.OccupancyMap);
        return *this; 
    } 

    // get the initial completed unoccupied map
    void instantiateOccupancyMap();      

    // update OccupancyMap 
    void updateOccupancyMap();      

    // plant the polymer from input file 
    void plantPolymersInGrid(std::string filename);  

    // extract number of polymers from the topology file 
    int ExtractNumberOfPolymers(std::string filename);

    // extract topology from the topology file 
    std::vector <std::string> ExtractContentFromFile(std::string filename); 


    void ExtractPolymersFromFile(std::string filename);
    

    void CalculateEnergy(); 

    void dumpPositionsOfPolymers (int step, std::string filename="dumpfile.txt"); 


    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#
    double EnergyPredictor(Particle p1, Particle p2){
        if (p1.ptype == "monomer" && p2.ptype == "monomer"){
            return this->Emm;
        }
        else if (p1.ptype=="solvent" && p2.ptype == "solvent"){
            return this->Ess; 
        }
        else if (p1.ptype == "monomer" && p2.ptype == "solvent"){
            if (p1.orientation == p2.orientation){
                return this->Ems_a;
            }
            else {
                return this->Ems_n;
            }
        }
        else if (p1.ptype=="solvent" && p2.ptype == "monomer"){
            if (p1.orientation == p2.orientation){
                return this->Ems_a;
            }
            else {
                return this->Ems_n;
            }

        }
        else {
            std::cout << "There is a bad energetic interaction."; 
            exit(EXIT_FAILURE);
        }
    };
    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#

    // check validity of input coords
    bool checkValidityOfCoords(std::vector <int> v);
    bool checkForOverlaps(std::vector <Polymer> PolymerVector); 
    bool checkConnectivity(std::vector <Polymer> PolymerVector); 

    //~#~#~~#~#~#~#~#~~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~~~##~#~#~#~##~~#~#~#~#


    std::vector <Particle> ClusterParticleMaker();
    std::vector <Particle> ClusterMaker(std::vector <Particle> Particles, std::vector <Particle> final, int count=0); 
    

}; 


// 
Grid CreateGridObject(std::string positions, std::string topology);
Grid IsingFlip(Grid InitialG);
Grid ZeroIndexRotation(Grid InitialG, int index);
Grid FinalIndexRotation(Grid InitialG, int index);
Grid EndRotation(Grid InitialG, int index); 
Grid KinkJump(Grid InitialG, int index); 
Grid CrankShaft(Grid InitialG, int index); 
Grid ForwardReptation(Grid InitialG, int index); 
Grid BackwardReptation(Grid InitialG, int index);
Grid Reptation(Grid InitialG, int index);
Grid Translation(Grid InitialG, int index, std::vector <int> direction);  
Grid MoveChooser(Grid InitialG, bool v); 


#endif 