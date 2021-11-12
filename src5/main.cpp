#include <iostream> 
#include <vector> 
#include <string> 
#include <map>
#include "classes.h"
#include "misc.h"

int main() {

    const int x{10}, y{10}, z{10}, kT{1}; 

    Grid G(x, y, z, kT); 

    int N = G.ExtractNumberOfPolymers("positions.txt");
    std::cout << "Number of polymers provided is " << N << std::endl;
    
    G.ExtractPolymersFromFile("positions.txt");
    G.PolymersInGrid.at(0).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    //G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    G.CalculateEnergy(); 
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;

    // =============================


    /*G.EndRotation_MC(0); 
    G.CalculateEnergy(); 

    G.PolymersInGrid.at(0).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    //G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;

    //==============================

    
    G.EndRotation_MC(0); 
    G.CalculateEnergy(); 

    G.PolymersInGrid.at(0).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;
    

    //==============================
    G.KinkJump_MC(0);
    G.CalculateEnergy(); 

    G.PolymersInGrid.at(0).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    //G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;
*/

    //==============================
    G.CrankShaft_MC(0);
    G.CalculateEnergy(); 

    
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    G.PolymersInGrid.at(0).printChainCoords();
    //G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;


    //==============================
    G.ZeroToFinalReptation_MC(0);
    G.CalculateEnergy(); 

    
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    G.PolymersInGrid.at(0).printChainCoords();
    //G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;

    //==============================
    G.FinalToZeroReptation_MC(0);
    G.CalculateEnergy(); 

    
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    G.PolymersInGrid.at(0).printChainCoords();
    //G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;

    for (int j=0; j< 3; j ++){
    std::cout << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~" << std::endl;
    //==============================
    G.Reptation_MC(0);
    G.CalculateEnergy(); 

    
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    G.PolymersInGrid.at(0).printChainCoords();
    //G.PolymersInGrid.at(1).printChainCoords();
    std::cout << "~=~=~+~+~+~+~+" << std::endl;
    std::cout << "Energy of the system is " << G.Energy << std::endl << std::endl;


    }

   



    return 0;

}