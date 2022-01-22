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


int main(){ 
    std::vector <int> coords {1,2,3,4,5,6,7,8}; 
    int i = 5; 

    if (std::find( coords.begin(), coords.begin() + 4, i) != coords.begin()+5){
        printf("found %d\n", i);
    }
    else {
        printf("not found.\n");
    }

    return 0;

}