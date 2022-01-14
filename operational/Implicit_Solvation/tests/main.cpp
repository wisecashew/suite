#include <iostream> 
#include <vector> 
#include <string>
#include <algorithm> 
#include <map>
#include <random>
#include <chrono>
#include <regex> 
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <stdlib.h>


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
    Particle () {};

    // constructor 
    Particle (std::vector <int> crds, std::string type_, int orientation_): coords (crds), ptype (type_), orientation (orientation_){

    }

    // destructor 
    ~Particle(){

    }

    // print location of the particle 
    void printCoords(); 

};


int main(int argc, char* argv[]){


  std::map <std::vector <int>, Particle> OccupancyMap;

  Particle p1 ({0,0,0}, "monomer", 0); 
  Particle p2 ({0,0,0}, "monomer", 0); 
  std::cout << std::boolalpha;
  bool b = (p1==p2);
  std::cout << "Check equality of p1, p2: " << b << std::endl;
  //std::cout << p.ptype << std::endl;
  //OccupancyMap[p.coords] = p ;
  //std::cout << OccupancyMap[p.coords].orientation << std::endl;
}