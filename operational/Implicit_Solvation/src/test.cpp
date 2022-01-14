#include <chrono> 
#include <random> 
#include <iostream>

int main() {

	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().count() ); 
	std::mt19937 generator (seed); 
	std::uniform_int_distribution<int> distribution (0, 10); 

	for (int i{0}; i<10; i++){
		std::cout << "rng is " << distribution(generator) << std::endl;
	}


	return 0; 
}
	
