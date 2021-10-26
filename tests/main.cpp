#include <iostream> 
#include <vector> 
#include <algorithm> 
#include "particle.h"

int main() {

	std::vector <int> l {1, 2 , 3}, m { 2,3,3}, n {1,3,5};  

	int a {2}; 

	std::vector <std::vector <int>> l2 = {l, m, n}; 

	if ( std::find(l2.begin(), l2.end(), l) != l2.end() ){
		std::cout << "the item was found" << std::endl;
	} 
	else {
		std::cout << "the item was not found. " << std::endl;
	}

	return 0; 
}