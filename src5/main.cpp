#include <iostream> 
#include <vector> 
#include <string> 
#include <map>
#include <chrono>
#include "classes.h"
#include "misc.h"

int main() {

    const int x{10}, y{10}, z{10}, kT{1}; 

    Grid G(x, y, z, kT); 

    int N = G.ExtractNumberOfPolymers("positions.txt");
    std::cout << "Number of polymers provided is " << N << std::endl;
    G.ExtractPolymersFromFile("positions.txt");
    std::cout << "Polymers extracted successfuly!" << std::endl;
    std::cout <<"Energy of the system is " << G.Energy << std::endl;
    //G.CalculateEnergy();



    int NumberOfMCMoves = 10; 

    auto start = std::chrono::high_resolution_clock::now(); 

    for (int i{0}; i < NumberOfMCMoves; i++){
        G.TheElementaryGridEvolver();
    }

    auto stop = std::chrono::high_resolution_clock::now(); 

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start); 

    std::cout << "Your function ran for " << duration.count() << "seconds." << std::endl;


    return 0;

}