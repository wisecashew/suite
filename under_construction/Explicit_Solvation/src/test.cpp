#include <iostream>
#include <array>
#include <algorithm> 
#include <vector> 

int main() {


	std::array <int,3> a1 = {1,0,0}, a2 = {2,0,0}, a3 = {3,0,0}, a4= {4,0,0}, a5={5,0,0}, a6={6,0,0};

	std::vector <std::array <int,3>> v1 = {a1,a2,a4,a5,a6};
	std::vector <std::array <int,3>> v2 = {a3,a4,a5};

	std::vector <std::array <int,3>> v3;

	std::set_intersection (v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(v3)); 

	for (auto a: v3){

		for (int i: a){
			std::cout << i << " | ";
		}
		std::cout << std::endl;

	}

	return 0;
}
