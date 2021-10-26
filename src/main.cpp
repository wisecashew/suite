#include <iostream> 
#include <vector> 
#include <string> 

#include "classes.h"
#include "misc.h"

int main(){

	int x_len {4}, y_len{4}, z_len{4}; 
	std::vector <int> seed{0,0,0};  // this is where the first particle will fall 

	// std::vector <Particle> gp; 
	int dop = 10;
	Grid Gp( x_len, y_len, z_len); 
	Gp.polymer_insertion(dop, seed);
	// Gp.print_loclist();

	std::vector <std::vector <int>> lat; 
	lat = create_lattice_pts(x_len, y_len, z_len);

	std::vector <std::vector <int>> solv_pts; 
	solv_pts = Gp.solvate(); 

	std::cout << "printing out solvated coordinates..." << std::endl;
	// for (auto v: solv_pts){
	//	print(v); 
	// }

	std::cout << "length of solvate coordinates is " << solv_pts.size() << std::endl; 
	std::cout << "length of polymer coordinates is " << std::endl;

	// std::cout << "Printing out lattice points: " << std::endl;

	// for (auto v: lat){
	// 	print(v);
	// }




	return 0;
}