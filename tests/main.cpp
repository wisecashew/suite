#include <iostream> 
#include <vector> 
#include <string>
#include <algorithm> 
#include <map>
#include "particle.h"
#include <random>
#include <chrono>
int main() {


	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
	for (int i{0}; i<10; i++){

    std::cout << "rng is " << distribution(generator) << std::endl;

	}
/*
	std::vector <int> l {1, 2 , 3}, m { 2,3,3}, n {1,3,5};  

	int a {2}; 

	std::vector <std::vector <int>> l2 = {l, m, n}; 

	if ( std::find(l2.begin(), l2.end(), l) != l2.end() ){
		std::cout << "the item was found" << std::endl;
	} 
	else {
		std::cout << "the item was not found. " << std::endl;
	}

	std::vector <int> v {1,2}; 

	std::map <std::string, std::vector <int> > M; 
	M["Turf"] = v;

	for (int i: M["Turf"]){
		std::cout << i << std::endl;
	}

	std::map <Particle, std::vector <Particle>> pmaps; 

*/
	return 0; 
}
