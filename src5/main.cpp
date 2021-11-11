#include <iostream> 
#include <vector> 
#include <string> 
#include <map>
#include "classes.h"
#include "misc.h"

int main() {

    const int x{10}, y{10}, z{10}; 

    Grid G(x, y, z); 

    std::cout << G.OccupancyMap[{0,0,0}] << std::endl;

    int N = G.ExtractNumberOfPolymers("testfile.txt");
    std::cout << "NoP is " << N << std::endl;
    G.ExtractPolymersFromFile("testfile.txt");
    G.PolymersInGrid.at(0).printChainCoords();
    std::cout << "Occupancy of {1,2,3} is " << G.OccupancyMap[{1,2,3}] << std::endl;

    return 0;

}