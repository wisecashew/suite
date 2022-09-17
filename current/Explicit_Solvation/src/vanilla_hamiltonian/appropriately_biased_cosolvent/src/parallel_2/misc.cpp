#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include <random>
#include <chrono>
#include "misc.h"
#include "classes.h"
#include <sstream>    
#include <fstream> 
#include <regex>
#include <tuple> 
#include <array>
#include <iterator>
#include <utility>
#include <unordered_set>
#include <algorithm>
#include <omp.h>
#define NUM_THREADS 20

/* 
=============================================================================================
These are some objects I have defined which I tend to use often. 
Helpful definitions which are employed often in the context of the z=6 lattice I am using.
=============================================================================================
*/ 

const std::vector <int> ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1}; 	// unit directions 
const std::vector <std::vector <int>> drns = {ex, nex, ey, ney, ez, nez};  							// vector of unit directions 

std::array <int,3> ax     = {1,0,0}, ay     = {0,1,0} , az     = {0,0,1} , nx     = {-1,0,0} ,  ny     = {0,-1,0}, nz     =  {0,0,-1} ; 
std::array <int,3> axay   = {1,1,0}, axaz   = {1,0,1} , axny   = {1,-1,0}, axnz   = {1,0,-1} , nxay    = {-1,1,0}, nxaz   =  {-1,0,1} , nxny = {-1,-1,0} , nxnz = {-1,0,-1}; 
std::array <int,3> ayaz   = {0,1,1}, aynz   = {0,1,-1}, nyaz   = {0,-1,1}, nynz   = {0,-1,-1};  
std::array <int,3> axayaz = {1,1,1}, axaynz = {1,1,-1}, axnyaz = {1,-1,1}, axnynz = {1,-1,-1},  nxayaz = {-1,1,1}, nxaynz = {-1,1,-1}, nxnyaz = {-1,-1,1}, nxnynz = {-1,-1,-1}; 
std::array <std::array <int,3>, 26> adrns = { ax, ay, az, nx, ny, nz, axay, axaz, axny, axnz, nxay, nxaz, nxny, nxnz, ayaz, aynz, nyaz, nynz, axayaz, axnyaz, axaynz, axnynz, nxayaz, nxaynz, nxnyaz, nxnynz }; 

std::map <int, std::array<double,3>> Or2Dir = { {0, {1.0,0,0}}, {1, {0,1.0,0}}, {2, {0,0,1}}, {3, {-1,0,0}}, {4, {0,-1,0}}, {5, {0,0,-1}}, {6, {1.0/(std::sqrt(2)), 1.0/(std::sqrt(2)), 0}}, {7, {1.0/(std::sqrt(2)), 0, 1.0/(std::sqrt(2))}}, {8, {1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}}, {9, {1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}}, {10, {-1.0/(std::sqrt(2)),1.0/(std::sqrt(2)),0}}, {11, {-1.0/(std::sqrt(2)),0,1.0/(std::sqrt(2))}}, {12, {-1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}}, {13, {-1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}}, {14, {0,1.0/(std::sqrt(2)),1.0/(std::sqrt(2))}}, {15, {0,1.0/(std::sqrt(2)),-1.0/(std::sqrt(2))}}, {16, {0,-1.0/(std::sqrt(2)), 1.0/(std::sqrt(2))}}, {17, {0,-1.0/(std::sqrt(2)), -1.0/(std::sqrt(2))}}, {18, {1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {19, {1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {20, {1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {21, {1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {22, {-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {23, {-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {24, {-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {25, {-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}} };
std::map <std::array<double,3>, int> Dir2Or = { {{1.0,0,0}, 0}, {{0,1.0,0}, 1}, {{0,0,1}, 2}, {{-1,0,0}, 3}, {{0,-1,0}, 4}, {{0,0,-1}, 5}, {{1.0/(std::sqrt(2)), 1.0/(std::sqrt(2)), 0}, 6}, {{1.0/(std::sqrt(2)), 0, 1.0/(std::sqrt(2))}, 7}, {{1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}, 8}, {{1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}, 9}, {{-1.0/(std::sqrt(2)),1.0/(std::sqrt(2)),0}, 10}, {{-1.0/(std::sqrt(2)),0,1.0/(std::sqrt(2))}, 11}, {{-1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}, 12}, {{-1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}, 13}, {{0,1.0/(std::sqrt(2)),1.0/(std::sqrt(2))}, 14}, {{0,1.0/(std::sqrt(2)),-1.0/(std::sqrt(2))}, 15}, {{0,-1.0/(std::sqrt(2)), 1.0/(std::sqrt(2))}, 16}, {{0,-1.0/(std::sqrt(2)), -1.0/(std::sqrt(2))}, 17}, {{1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 18}, {{1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 19}, {{1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 20}, {{1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 21}, {{-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 22}, {{-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 23}, {{-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 24}, {{-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 25} };

//==============================================================================================
// impose periodic boundary conditions on vector 
//==============================================================================================

void impose_pbc(std::vector <int>* vect, int x_len, int y_len, int z_len){
	for (int i{0}; i<3; i++){
		if (i==0){
			(*vect)[i] = ((((*vect)[i])%x_len)+x_len)%x_len; 
		}
		else if (i==1){
			(*vect)[i] = ((((*vect)[i])%y_len)+y_len)%y_len; 
		}
		else {
			(*vect)[i] = ((((*vect)[i])%z_len)+z_len)%z_len; 
		}

	}
	
	return; 
}

void impose_pbc(std::array <int,3>* arr, int x_len, int y_len, int z_len) {
	for (int i{0}; i<3; ++i){
		if (i==0){
			(*arr)[i] = (((*arr)[i]%x_len)+x_len)%x_len; 
		}
		else if (i==1){
			(*arr)[i] = (((*arr)[i]%y_len)+y_len)%y_len; 	
		}
		else {
			(*arr)[i] = (((*arr)[i]%z_len)+z_len)%z_len; 
		}

	}

	return; 
}


void impose_pbc(std::array <double,3>* arr, int x_len, int y_len, int z_len) {
	for (int i{0}; i<3; ++i){
		if (i==0){
			(*arr)[i] = std::fmod( ( std::fmod( (*arr)[i], x_len) +x_len), x_len); 
		}
		else if (i==1){
			(*arr)[i] = std::fmod( ( std::fmod( (*arr)[i], y_len) +y_len), y_len); 
			// (*arr)[i] = (((*arr)[i]%y_len)+y_len)%y_len; 	
		}
		else {
			(*arr)[i] = std::fmod( ( std::fmod( (*arr)[i], z_len) +z_len), z_len); 
			// (*arr)[i] = (((*arr)[i]%z_len)+z_len)%z_len; 
		}

	}

	return; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
// impose a modified modulo behavior 
//==============================================================================================

int modified_modulo(int divident, int divisor){
	double midway = static_cast<double>(divisor/2); 
	
	
	if ( ((divident%divisor)+divisor)%divisor > midway){
		
		return ((divident%divisor)+divisor)%divisor-divisor; 
		
	}
	else {
		// std::cout << "result is " << result << std::endl;
		return (((divident%divisor)+divisor)%divisor);
	}

}

double modified_modulo ( double divident, int divisor){
	
	double midway = static_cast<double>(divisor/2); 
	
	if ( std::fmod(std::fmod(divident, divisor)+divisor,divisor) > midway){
		
		return std::fmod((std::fmod(divident,divisor)+divisor), divisor)-divisor; 
		
	}
	else {
		// std::cout << "result is " << result << std::endl;
		return std::fmod(std::fmod(divident,divisor)+divisor,divisor);
	}

}

// modifies array to be a legit direction, if it can be done 
// i.e. (9, 0, 0) becomes (-1, 0, 0) 

void modified_direction(std::array<int,3>* a, int x, int y, int z){
	
	for (int i{0}; i<3; ++i){
		switch (i){
			case (0):
				(*a)[i] = modified_modulo((*a)[i], x); 
				break; 
			case (1): 
				(*a)[i] = modified_modulo((*a)[i], y);
				break;
			case (2):
				(*a)[i] = modified_modulo((*a)[i], z); 
				break;
			}
	}

	return;
	
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
// random number generators 
//==============================================================================================

int rng_uniform(int start, int end){
	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (start, end);

    return distribution(generator); 

}

double rng_uniform(double start, double end){
	
	std::default_random_engine generator; 
	std::uniform_real_distribution <double> distribution (start, end); 
	generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	
	return distribution (generator); 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
// functions to access certain locations on the LATTICE object  
//==============================================================================================

int lattice_index ( std::array<int,3> location, int y, int z ){

	int idx = (location)[2]*(z*z) + (location)[1]*(y) + (location)[0];

	return idx; 

}

std::array <int,3> location ( int lattice_index, int x, int y, int z){

	int zcoord = lattice_index / (z*z);
	int ycoord = (lattice_index % (z*z))/ y; 
	int xcoord = ((lattice_index % (z*z) ) % y) % x; 

	std::array <int,3> loc = {xcoord, ycoord, zcoord}; 

	return loc; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
// functions to print data types. primarily vectors and arrays. 
//==============================================================================================

void print(std::vector <int> v, std::string c){
	for (int i: v){
		std::cout << i << " | ";
	}
	std::cout << c;
}

void print(std::vector <double> v, std::string c){
	for (double i: v){
		std::cout << i << " | ";
	}
	std::cout << c;
}

void print(std::array <int, 3> v, std::string c){
	for (int i: v){
		std::cout << i << " | ";
	}
	std::cout << c;
}

void print(std::array <int, 6> v, std::string c){
	for (int i: v){
		std::cout << i << " | ";
	}
	std::cout << c;
}

void print(std::array <double, 3> v, std::string c){
	for (double i: v){
		std::cout << i << " | ";
	}
	std::cout << c;
}

void print(std::array <double, 4> v, std::string c){
	for (double i: v){
		std::cout << i << " | ";
	}
	std::cout << c;
}

void print(std::array <double, 6> v, std::string c){
	for (double i: v){
		std::cout << i << " | ";
	}
	std::cout << c;
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
// functions to print data types. primarily vectors and arrays. 
//==============================================================================================

void print(std::vector <std::vector <int>>vv){
	for (std::vector <int> v: vv){
		print(v);
	}
}

void print(std::vector <std::vector <double>>vv){
	for (std::vector <double> v: vv){
		print(v);
	}
}

void print(std::vector <std::string>vv){
	for (std::string s: vv){
		std::cout << s << std::endl;
	}
}

void print(std::vector <Particle> pvec){
	for (Particle p: pvec){
		p.printCoords();
		// std::cout << p.ptype << std::endl;
	}
}

void print( std::array <std::array<int,3>,6> aa ){

	for ( std::array<int,3>& a: aa){
		print(a);
	}

}

void print ( std::vector <std::array<int,3>> aa ){

	for ( std::array <int,3>& a: aa){
		print (a);
	}

}

void print ( std::vector <Particle*>* LATTICE ){

	for ( Particle*& p: (*LATTICE) ){
		print( p->coords );
	}
}

void print ( std::vector <std::vector <std::array<int,3>>> V){

	for ( std::vector <std::array<int,3>>& v: V){
		for (std::array<int,3>& a: v){
			print (a); 
		}
	}
	return;
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
// functions to perform arithmetic with vectors and arrays 
//==============================================================================================

std::vector <int> add_vectors(std::vector <int>* v1, std::vector <int>* v2){
	size_t s = (*v1).size(); 
	std::vector <int> v3 (s,0); 
	for (int i{0}; i < static_cast<int>(s); i++){
		v3.at(i) = (*v2).at(i) + (*v1).at(i);
	}

	return v3; 
}

std::array <int,3> add_arrays(std::array <int,3>* a1, std::array <int,3>* a2){
	std::array<int, 3> a3; 

	for (int i{0}; i<3; ++i){
		a3[i] = (*a1)[i] + (*a2)[i]; 
	}

	return a3;

}

std::array <double,3> add_arrays(std::array <double,3>* a1, std::array <double,3>* a2){
	std::array<double, 3> a3; 

	for (int i{0}; i<3; ++i){
		a3[i] = (*a1)[i] + (*a2)[i]; 
	}

	return a3;

}

std::array <double,3> add_arrays(std::array <int,3>* a1, std::array <double,3>* a2){
	std::array<double, 3> a3; 

	for (int i{0}; i<3; ++i){
		a3[i] = (*a1)[i] + (*a2)[i]; 
	}

	return a3;

}

std::vector <int> subtract_vectors(std::vector <int>* v1, std::vector <int>* v2){
	size_t s = (*v1).size(); 
	std::vector <int> v3 (s,0); 
	for (int i{0}; i < static_cast<int>(s); i++){
		v3.at(i) = (*v1).at(i) -  (*v2).at(i);
	}

	return v3; 
}

std::array <int,3> subtract_arrays(std::array <int,3>* a1, std::array <int,3>* a2){
	std::array<int, 3> a3; 

	for (int i{0}; i<3; ++i){
		a3[i] = (*a1)[i] - (*a2)[i]; 
	}

	return a3;

}

std::array <double,3> subtract_arrays(std::array <double,3>* a1, std::array <double,3>* a2){
	std::array<double, 3> a3; 

	for (int i{0}; i<3; ++i){
		a3[i] = (*a1)[i] - (*a2)[i]; 
	}

	return a3;

}


std::array <double,3> scale_arrays ( double scalar, std::array <double,3>* array){

	std::array <double,3> arr_n;
	for (int i{0}; i<3; ++i){
		arr_n [i] = (*array)[i]*scalar; 
	}
	return arr_n;

}

std::array <double,3> scale_arrays ( double scalar, std::array <int,3>* array){

	std::array <double,3> arr_n;
	for (int i{0}; i<3; ++i){
		arr_n [i] = (*array)[i]*scalar; 
	}
	return arr_n;
	
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
// functions to perform geometric operations on arrays 
//==============================================================================================

double distance_between_points (std::array <double,3>* a1, std::array <double,3>* a2, int xlen, int ylen, int zlen){

	std::array <double,3> delta = subtract_arrays (a1, a2); 
	// std::cout << "delta array is "; print (delta);
	for ( int i{0}; i<3; ++i){
		if (i==0){
			delta[i] = modified_modulo ( delta[i], xlen );
		}
		else if ( i == 1 ){
			delta[i] = modified_modulo ( delta[i], ylen ); 
		}
		else if ( i == 2 ){
			delta[i] = modified_modulo ( delta[i], zlen ); 
		}
	}

	// impose_pbc ( &delta, xlen, ylen, zlen ); 

	return std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

double take_dot_product (int o1, int o2){

	// take the dot product 
	double dot_prod = 0;	
	std::array <double,3> d1 = Or2Dir [o1]; 
	std::array <double,3> d2 = Or2Dir [o2]; 

	for (int i{0}; i<3; ++i){
		dot_prod += d1[i]*d2[i];
	}

	return dot_prod; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
// function to check for avoidance in the walk 
//==============================================================================================

bool check_avoidance(const std::vector <int> to_check, const std::vector<std::vector <int>> loc_list){
	for (std::vector <int> v: loc_list){
		if (v==to_check){
			return false;
		}
	}
	return true;
}



//==============================================================================================
// executing self-avoiding random walk
//==============================================================================================

void sarw(std::vector<std::vector<int>>* loc_list, int dop){
	
	if (dop == 0){
		return; // once the degree of polymerization hits zero, you are done
	}; 

	// until then 
	// increment final vector in loc_list in a unit direction 
	std::vector <int> next(3,0); 
	for (auto v: drns){
		 
		next = add_vectors(&((*loc_list).at((*loc_list).size()-1)), &v);
		 
		if (check_avoidance(next, *loc_list)){
			dop--; // decrease dop now that a monomer unit has been added 
			(*loc_list).push_back(next); // add monomer to list 
			return sarw(loc_list, dop); 
		} 
	}

	std::cout << "no solution found for the self-avoiding random walk...";
	return; 
}


//==============================================================================================
// executing sarw with pbc
//==============================================================================================

void sarw_pbc(std::vector<std::vector<int>>* loc_list, int dop, int x_len, int y_len, int z_len){
	
	if (dop == 0){
		return; // once the degree of polymerization hits zero, you are done
	}; 

	// until then 
	// increment final vector in loc_list in a unit direction 
	std::vector <int> next(3,0); 
	for (auto v: drns){

		next = add_vectors(&((*loc_list).at((*loc_list).size()-1)), &v);
		
		impose_pbc(&next, x_len, y_len, z_len); 
		
		
		if (check_avoidance(next, *loc_list)){
			dop--; // decrease dop now that a monomer unit has been added 
			(*loc_list).push_back(next); // add monomer to list 
			return sarw_pbc(loc_list, dop, x_len, y_len, z_len); 
		} 
	}

	std::cout << "no solution found for the self-avoiding random walk...";
	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
// create_linked_list is a function that takes two lists and creates a well-ordered linked-list 
// between the two lists 
//==============================================================================================

void create_linked_list ( std::vector<std::array<int,3>> v1, std::vector<std::array<int,3>> v2, \
	std::vector <std::array<int,3>> link, std::vector <std::vector <std::array<int,3>>>* master_linked_list, \
	int beginning){


	if ((v1).size() == 0 ){
		if (link.size() != 0){
			(*master_linked_list).push_back ( link ); 
		}
		return; 
	}

	if (beginning){
		(link) = { (v1)[0], (v2)[0]}; 
		(v1).erase ((v1).begin()); 
		(v2).erase ((v2).begin()); 
		beginning = 0; 
	}

	std::vector<std::array<int,3>>::iterator it = std::find ( (v1).begin(), (v1).end(), (link)[ link.size()-1 ] );

	if ( it != (v1).end() ){

		int idx = it - (v1).begin(); 

		link.push_back ( (v1)[idx] );
		link.push_back ( (v2)[idx] );

		

		(v1).erase (v1.begin()+idx );
		(v2).erase (v2.begin()+idx );


		create_linked_list (v1, v2, link, master_linked_list, beginning); 

	}

	else {
		
		beginning = 1; 
		(*master_linked_list).push_back ( link );
		(link).clear(); 
		create_linked_list (v1, v2, link, master_linked_list, beginning); 

	}

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of create_linked_list
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
// functions to generate neighbor lists given the location of a particle 
//==============================================================================================

std::vector <std::vector <int>> obtain_ne_list(std::vector <int> loc, int x_len, int y_len, int z_len){
	std::vector <std::vector <int>> nl;
	nl.reserve(6); 
	for (std::vector <int> d: drns){
		std::vector <int> v = add_vectors(&loc, &d); 
		impose_pbc(&v, x_len, y_len, z_len); 
		nl.push_back(v);
	}

	std::sort(nl.begin(), nl.end() );
	nl.erase(std::unique (nl.begin(), nl.end() ), nl.end() );  

	return nl; 

}

std::array <std::array <int,3>, 26> obtain_ne_list(std::array <int,3> loc, int x_len, int y_len, int z_len){
	std::array <std::array <int,3>,26> nl;
	int i {0}; 
	for (std::array <int,3> d: adrns) {
		std::array <int,3>  a = add_arrays(&loc, &d); 
		impose_pbc(&a, x_len, y_len, z_len); 
		nl[i] = a;
		++i; 
	}

	// std::sort(nl.begin(), nl.end() );
	// nl.erase(std::unique (nl.begin(), nl.end() ), nl.end() );  

	return nl; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
// functions for acceptance of a position 
//==============================================================================================

bool acceptance(int dE, double kT){
	// prob = min(1, exp(-1/(kT)*dE))
	std::vector <double> P {1, std::exp(-dE/kT)}; 
	double prob = *std::min_element(P.begin(), P.end()) ; 


	// rng stuff 
	std::default_random_engine generator; 
	std::uniform_real_distribution <double> distribution (0.0, 1.0); 
	generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	double num = distribution (generator); 

	//std::cout << "acceptance has been called..." << std::endl; 
	// std::cout << "dE is " << dE << std::endl; 
	// std::cout << "probability of move is " << prob << std::endl; 
	// std::cout << "rng is " << num << std::endl;
	
	// std::cout << "random number is " << num << std::endl; 
	// std::cout << "probability is " << prob << std::endl;
	if (num <= prob){
		return true;
	}
	else{
		return false; 
	}

}

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: ExtractNumberOfPolymers 
//
// PARAMETERS: std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the positions file, and extracts the number of polymers from that file. 
// 
// DEPENDENCIES: subtract_vectors 
//
// THE CODE: 

int ExtractNumberOfPolymers(std::string filename){
    std::regex start("START"); 
    std::regex end ("END"); 
    std::ifstream myfile (filename); 

    if ( !myfile ){
        std::cerr << "File named " << filename << " could not be opened!" << std::endl; 
        exit(EXIT_FAILURE);
    }

    std::string myString; 
    // std::vector <std::string> StartContents; 
    // std::vector <std::string> EndContents; 

    int numberOfStarts {0}, numberOfEnds {0}; 

    // std::cout << "This is from ExtractNumberOfPolymers..." << std::endl;
    if (myfile.is_open()) {
        while ( std::getline(myfile, myString) ) {
            
            // std::cout << myString << std::endl; 
            // std::getline(myfile, myString); // pipe file's content into stream 

            if (std::regex_search(myString, start)){
                ++numberOfStarts; 
            }
            else if (std::regex_search(myString, end)){
                ++numberOfEnds; 
            }
            else if ( myString.empty()){
                std::cerr << "ERROR: Empty line found. Bad positions file. " << std::endl;
                exit (EXIT_FAILURE);
            }
        }
    }

    if (numberOfStarts==numberOfEnds){
        return numberOfStarts;
    }
    else {
        std::cerr << "Number of starts is not the same as number of ends. Bad input file." << std::endl;
        return 0;
    }

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of ExtractNumberOfPolymers. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: ExtractContentFromFile
//
// PARAMETERS: std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the positions file, and extracts all the data from the file 
// in the form of a vector of strings.  
// 
// DEPENDENCIES: none apart from the STL
//
// THE CODE: 

std::vector <std::string> ExtractContentFromFile(std::string filename){
    std::ifstream myfile (filename); 

    if ( !myfile ){
        std::cerr << "File named " << filename << " could not be opened!" << std::endl; 
        exit(EXIT_FAILURE);
    }

    std::string mystring; 
    std::vector <std::string> contents; 

    if (myfile.is_open() ){
        while ( std::getline(myfile, mystring) ) {
            // pipe file's content into stream 
            contents.push_back(mystring); 
        }
    }

    return contents; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of ExtractContentFromFile. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: NumberExtractor 
//
// PARAMETERS: std::string line_from_a_file 
//
// WHAT THE FUNCTION DOES: it looks at a string, and extracts all the numbers in the file. 
// used to get box dimensions and temperature and so on from the topology file. 
// 
// DEPENDENCIES: none apart from the STL. 
// 
// THE CODE: 

double NumberExtractor(std::string s){

	double info {0}; 
	std::stringstream ss (s); 
	std::string temp; 

	double found; 

	while (!(ss.eof()) ) {
		ss >> temp; 

		if (std::stringstream(temp) >> found){
			info = found; 
		}

		temp = "" ; //not sure what this is for; 
	}

	return info; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//                 End of NumberExtractor
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: ExtractTopologyFromFile 
//
// PARAMETERS: std::string name_of_file 
//
// WHAT THE FUNCTION DOES: it looks at a string, and extracts all the numbers in the file. 
// used to get box dimensions and temperature and so on from the topology file. 
// 
// DEPENDENCIES: none apart from the STL. 
// 
// THE CODE: 

/*
std::array <double,8> ExtractTopologyFromFile(std::string filename){
    
    std::array <double, 8> info_vec; 
    double info; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(filename); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), Emm_a ("Emm_a"), Emm_n ("Emm_n"), Ems_a ("Ems_a"), Ems_n ("Ems_n"), eof ("END OF FILE"); 
    //bool out_mat = true; 

    // print(contents);


    for (std::string s: contents){

    	if (std::regex_search(s, x)){
    		info = NumberExtractor(s); 
    		info_vec[0]=info; 
    		continue; 
    	}

    	else if (std::regex_search(s, y)){
    		double info = NumberExtractor(s); 
    		info_vec[1] = info; 
    		continue; 
    	}

    	else if (std::regex_search(s, z)){
    		double info = NumberExtractor(s); 
    		info_vec[2] = info; 
    		continue; 
    	}

		else if (std::regex_search(s, kT)){
    		double info = NumberExtractor(s); 
    		info_vec[3] = info ; 
    		continue; 
    	}

    	else if (std::regex_search (s, Emm_a)){
    		double info = NumberExtractor(s); 
    		info_vec[4] = info; 
    		continue; 
    	}

    	else if (std::regex_search (s, Emm_n)){
    		double info = NumberExtractor(s); 
    		info_vec[5] = info; 
    		continue; 
    	}

    	else if (std::regex_search (s, Ems_a)){
    		double info = NumberExtractor(s); 
    		info_vec[6] = info; 
    		continue; 
    	}

    	else if (std::regex_search (s, Ems_n)){

    		double info = NumberExtractor(s);
    		info_vec[7] = info; 
    		continue;
    	}

    	else if (std::regex_search(s, eof)){
    		// std::cout << "End of topology file." << std::endl;
    		break;
    	}

    	else {
    		std::cout << s << std::endl;
    		std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
    		exit(EXIT_FAILURE); 
    	}

    }


    return info_vec;

}
*/
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

std::array <double,13> ExtractTopologyFromFile(std::string filename){
    
    std::array <double, 13> info_vec; 
    double info; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(filename); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), Emm_a ("Emm_a"),\
    Emm_n ("Emm_n"), Ems1_a ("Ems1_a"), Ems1_n ("Ems1_n"), Ems2_a ("Ems2_a"), Ems2_n ("Ems2_n"), Es1s2_a("Es1s2_a"), Es1s2_n ("Es1s2_n"), \
    frac ("frac"), eof ("END OF FILE"); 
    //bool out_mat = true; 

    // print(contents);


    for (std::string s: contents){

    	if (std::regex_search(s, x)){
    		info = NumberExtractor(s); 
    		info_vec[0]=info; 
    		continue; 
    	}

    	else if (std::regex_search(s, y)){
    		double info = NumberExtractor(s); 
    		info_vec[1] = info; 
    		continue; 
    	}

    	else if (std::regex_search(s, z)){
    		double info = NumberExtractor(s); 
    		info_vec[2] = info; 
    		continue; 
    	}

		else if (std::regex_search(s, kT)){
    		double info = NumberExtractor(s); 
    		info_vec[3] = info ; 
    		continue; 
    	}

    	else if (std::regex_search (s, Emm_a)){
    		double info = NumberExtractor(s); 
    		info_vec[4] = info; 
    		continue; 
    	}

    	else if (std::regex_search (s, Emm_n)){
    		double info = NumberExtractor(s); 
    		info_vec[5] = info; 
    		continue; 
    	}

    	else if (std::regex_search (s, Ems1_a)){
    		double info = NumberExtractor(s); 
    		info_vec[6] = info; 
    		continue; 
    	}

    	else if (std::regex_search (s, Ems1_n)){

    		double info = NumberExtractor(s);
    		info_vec[7] = info; 
    		continue;
    	}
    	else if (std::regex_search (s, Ems2_a)){

    		double info = NumberExtractor(s);
    		info_vec[8] = info; 
    		continue;
    	}
    	else if (std::regex_search (s, Ems2_n)){

    		double info = NumberExtractor(s);
    		info_vec[9] = info; 
    		continue;
    	}
    	else if (std::regex_search (s, Es1s2_a)){

    		double info = NumberExtractor(s);
    		info_vec[10] = info; 
    		continue;
    	}
    	else if (std::regex_search (s, Es1s2_n)){

    		double info = NumberExtractor(s);
    		info_vec[11] = info; 
    		continue;
    	}
    	else if (std::regex_search (s, frac)) {
    		double info = NumberExtractor(s);
    		info_vec[12] = info;
    		continue; 
    	}

    	else if (std::regex_search(s, eof)){
    		// std::cout << "End of topology file." << std::endl;
    		break;
    	}

    	else {
    		std::cout << s << std::endl;
    		std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
    		exit(EXIT_FAILURE); 
    	}

    }


    return info_vec;

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// End of ExtractTopologyFromFile
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: makePolymer  
//
// PARAMETERS: locations of monomers, and type of monomers 
//
// WHAT THE FUNCTION DOES: it takes the positions and types provided by the user and makes the Polymer object. 
// 
// DEPENDENCIES: none apart from the STL. 
// 
// THE CODE: 


Polymer makePolymer(std::vector <std::array <int,3> > locations, std::string type_m){
	std::vector <int> pmer_spins; 
    short size_ = locations.size(); 
    for (short i=0; i<size_; i++){
        // unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
        std::random_device rd;
        std::mt19937 generator(rd () ); 
        // std::cout << "check..." << std::endl;
        std::uniform_int_distribution<int> distribution (0, 25); 
        pmer_spins.push_back(distribution(generator));
        
    }

    std::vector <Particle*> ptc_vec; 

    for (int i=0;i< size_ ; i++ ){
        Particle* p = new Particle (locations.at(i), type_m, pmer_spins.at(i)); 
        ptc_vec.push_back(p); 
    }

    Polymer pmer (size_, ptc_vec);

    return pmer; 
}

// WHAT THIS FUNCTION DOES:
// It takes a spins vector, which assigns a specific spin to each monomer unit. This is not the case for the above 
// function. 


Polymer makePolymer(std::vector <std::array <int,3> > locations, std::vector<int> pmer_spins, std::string type_m){

    std::vector <Particle*> ptc_vec; 
    short size_ = locations.size(); 

    for (int i=0; i < size_ ; i++ ){
        Particle* p = new Particle (locations.at(i), type_m, pmer_spins.at(i)); 
        ptc_vec.push_back(p); 
    }

    Polymer pmer (size_, ptc_vec);

    return pmer; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//				End of makePolymer 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//============================================================
//============================================================
// 
// NAME OF FUNCTION: ExtractPolymersFromTraj
//
// PARAMETERS: std::string trajectory, std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the positions file, and extracts all the data from the file 
// in the form of a vector of strings.  
// 
// DEPENDENCIES: ExtractContentFromFile, makePolymer, checkValidityOfCoords 
//
// THE CODE: 


std::vector <Particle*> ExtractLatticeFromRestart ( std::string rfile, int* step_num, int x, int y, int z ){

	std::vector <Particle*> LATTICE; 
	LATTICE.reserve (x*y*z);

	std::vector <std::string> contents = ExtractContentFromFile ( rfile );

	std::regex start ("FINAL STEP: "), end ("END"); 
	std::regex numbers ("[0-9]+"); 
	std::regex characters ("[a-z][0-9]+");

	int orientation   = -1; 
	std::string ptype = "x"; 
	int index         = -1; 
	std::smatch match; 


	for ( std::string& s: contents) {

		// send content into stringstream 

		if ( std::regex_search (s, start) ) {
			std::regex_search ( s, match, numbers ); 
			// std::cout << "match for step number is " << match[0] << std::endl;
			*step_num = std::stoi(match[0].str()); 
		}

		else if ( std::regex_search (s, end) ){
			break; 
		}

		else { 
			
			std::regex_search ( s, match, numbers );
			std::regex_token_iterator<std::string::iterator> rend; 
			std::regex_token_iterator<std::string::iterator> a ( s.begin(), s.end(), numbers );

			for ( int i=0; i<2; ++i ){
				if ( i == 0 ){
					// std::cout << "match for orientation is " << *a << std::endl;
					orientation = std::stoi ( *a );
					// std::cout << "match for orientation is " << *a << std::endl;
					*a++;  
				}
				else {
					*a++; 
					// std::cout << "match for index is " << *a << std::endl;
					index = std::stoi ( *a );
					// std::cout << "match for index is " << *a << std::endl;
				}
			}
 

			std::regex_search ( s, match, characters );
			ptype = match[0].str(); 

			// std::cout << "ptype is " << ptype << std::endl; 

			Particle* p_ptr = new Particle (location(index, x, y, z), ptype, orientation); 

			LATTICE.insert ( LATTICE.begin() + index, p_ptr); 

		}

	}

	std::cout << "Created lattice from file!" << std::endl;
	return LATTICE; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// End of ExtractPolymersFromRestart
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

std::vector <Polymer> ExtractPolymersFromTraj(std::string trajectory, std::string position, int final_step_num, int x, int y, int z){

    int NumberOfPolymers = ExtractNumberOfPolymers(position); 

    std::vector <Polymer> PolymerVector; 
    PolymerVector.reserve(NumberOfPolymers);

    std::vector <std::array <int,3>> locations; 
    std::vector <int> spins; 

    std::vector <std::string> contents = ExtractContentFromFile(trajectory); // this extracts every line of the file

    // std::vector <int> step_num_store; 
    std::vector <int> index_store; // need this guy to hold the index of the final set of coordinates. 

    // std::cout << "Is content being extracted?" << std::endl;

    
    bool step_bool {false}, start_bool {false}, end_bool {false}; 

    std::regex start ("START"), end ("END"), step ("step " + std::to_string(final_step_num) );
    std::regex step_generic ("step"); 
    std::regex reg_poly ("Dumping coordinates of Polymer"); 
    std::regex numbers ("[0-9]+"); 

    int startCount{0}, endCount{0}; 
    std::array <int,3> loc;
    // std::stringstream ss; 
    std::smatch match; 

    // std::cout << "final step number is " << final_step_num << std::endl;

    for (std::string& s: contents){
    	
    	if ( std::regex_search (s, step_generic) ){

    		std::regex_search ( s, match, numbers ); 
			std::regex_token_iterator<std::string::iterator> rend; 
			std::regex_token_iterator<std::string::iterator> a ( s.begin(), s.end(), numbers );

			// std::cout << "*a is " << std::stoi(*a) << std::endl;

			if ( std::stoi(*a) > final_step_num ){
				// std::cout << "*a is " << std::stoi(*a) << std::endl;
				// std::cout << "final_step_num = " << final_step_num << std::endl;				
				std::cerr << "\n\nYou coordinates file and restart file are not in sync. \nIt is probably because the restart file you chose is not for that simulation.\nPlease check them out. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}

    	}

    	if ( std::regex_search ( s, step) ) {
    		// std::cout << s << std::endl;
    		step_bool = true; 
    		continue; 
    	}

        // sending to the stringstream
        // ss (s); 

        if ( step_bool ){

	        if ( std::regex_search(s, start) ){
	            ++startCount;
	            start_bool = true; 
	            end_bool   = false; 
	            continue; 
	        }

	        else if (std::regex_search(s, end) ) {
	            ++endCount;
	            start_bool = false; 
	            end_bool   = false; 
	            step_bool  = false; 

	            Polymer pmer = makePolymer(locations, spins);
	            PolymerVector.push_back(pmer);
	            
	            locations.clear();
	            break;
	            // continue; 
	            
	        }

	        else if (start_bool == end_bool){
	            continue;
	        }

	        else{
	        	// std::cout << s << std::endl;
	            std::regex_search ( s, match, numbers ); 
				std::regex_token_iterator<std::string::iterator> rend; 
				std::regex_token_iterator<std::string::iterator> a ( s.begin(), s.end(), numbers );



				for ( int i=0; i<4; ++i ){
					if ( i ==0 ) {
						// std::cout << "x-coords is " << *a << std::endl; 
						loc[i] = std::stoi ( *a );
						*a++; 
					}
					else if ( i==1 ){
						// std::cout << "y-coord is " << *a << std::endl; 
						loc[i] = std::stoi ( *a );
						*a++; 
					}
					else if ( i==2 ){
						// std::cout << "z-coord is " << *a << std::endl; 
						loc[i] = std::stoi ( *a );
						*a++;
					}
					else if ( i==3 ){
						spins.push_back ( std::stoi(*a) ); 
					}
				}  
				
				// std::cout << "Location is "; print (loc);
	            
	            if (!checkValidityOfCoords(loc, x, y, z)){
	            	std::cerr << "Coordinates are out of bounds. Bad input file." << std::endl;
	            	exit(EXIT_FAILURE); 
	            }
	        
	            locations.push_back(loc); 
	            
	        }
    	}
    }

    if(!(checkForOverlaps(PolymerVector))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Overlap detected." << std::endl;
        exit(EXIT_FAILURE);  
    }

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // throw in a check for connectivity of polymer chains 
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    if (!(checkConnectivity(PolymerVector, x, y, z))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Monomer units are not adjacent to one another on the lattice." << std::endl;
        exit(EXIT_FAILURE); 
    }
    
    

    return PolymerVector; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			End of ExtractPolymersFromRestart
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			End of Extractor functions 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: SetUpLatticeFromScratch and SetUpLatticeFromRestart 
//
// PARAMETERS: Polymers, LATTICE, polymer coords file name 
//
// WHAT THE FUNCTION DOES: it takes the positions and types provided by the user and makes the Polymer object and LATTICE OBJECT. 
// 
// DEPENDENCIES: none apart from the STL. 
// 
// THE CODE: 

void SetUpLatticeFromScratch (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector<Particle*>* LATTICE, std::string positions, double frac, int x, int y, int z){

	(*Polymers) = ExtractPolymersFromFile (positions, x, y, z); 
	AddSolvent (LATTICE, x, y, z); 

	// populate the lattice 
	for (Polymer& pmer: (*Polymers)) {
		for (Particle*& p: pmer.chain){

			// now that I have my polymer coordinates, time to create the grand lattice 
			(*LATTICE).at(lattice_index (p->coords, y, z) ) = p; 
		
		}
	}

	AddCosolvent ( Cosolvent, LATTICE, frac, static_cast<int>((*Polymers)[0].chain.size()), x, y, z ); 

	return;
}


void SetUpLatticeFromRestart ( std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector<Particle*>* LATTICE, int* step_number, std::string lattice_file_read, std::string dfile, std::string positions, int x, int y, int z ){
	
	(*LATTICE)  = ExtractLatticeFromRestart ( lattice_file_read, step_number, x, y, z ); 
    (*Polymers) = ExtractPolymersFromTraj   ( dfile, positions, *step_number, x, y, z );

    for ( Polymer& pmer: (*Polymers)) {
        for ( Particle*& p: pmer.chain){

            if ( !( (*LATTICE) [ lattice_index(p->coords, y, z) ]->ptype == p->ptype && (*LATTICE) [ lattice_index(p->coords, y, z) ]->coords == p->coords && (*LATTICE) [ lattice_index(p->coords, y, z) ]->orientation == p->orientation) ) {
                std::cerr << "There is a problem with extraction for restart..." << std::endl;
                exit (EXIT_FAILURE); 
            } 
            (*LATTICE) [ lattice_index(p->coords, y, z) ] = p;
        }
    }

    for (int i{0}; i< x*y*z; ++i){
    	if ( (*LATTICE).at(i)->ptype == "s2" ){
    		(*Cosolvent).push_back((*LATTICE).at(i));
    	}
    }

    return; 

}


void BiasTheStart ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z ) {

    std::array <std::array<int,3>, 26> ne_list;
    std::vector <int> solvent_indices; 
    
    for ( Polymer& pmer: *Polymers ){
        for ( Particle*& p: pmer.chain ) {
            
            ne_list = obtain_ne_list (p->coords, x, y, z);
            for ( std::array<int,3>& ne: ne_list ){
                if ( (*LATTICE).at(lattice_index(ne, y, z))->ptype[0] == 's' ) {
                    solvent_indices.push_back (lattice_index(ne, y, z));
                }
            }
        }
    }  

    // get rid of duplicates 
    std::unordered_set <int> s (solvent_indices.begin(), solvent_indices.end() );
    solvent_indices.assign (s.begin(), s.end() ); 

    for ( Polymer& pmer: (*Polymers) ) {
        for ( Particle*& p: pmer.chain ) {
            p->orientation = 0;
        }
    }

    for (int i: solvent_indices){
        (*LATTICE)[i]->orientation = 0;
    }

    return;

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of SetUp
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: AddSolvent
//
// PARAMETERS: int x, int y, int z, std::vector <Particle*>* LATTICE 
// 
// WHAT THE FUNCTION DOES: it looks at the polymer vector and solvates the lattice with solvent molecules 
// 
// DEPENDENCIES: LATTICE
//
// THE CODE: 

void AddSolvent (std::vector <Particle*>* LATTICE, int x, int y, int z){

	int c_idx {-1};
	std::array <int,3> loc = {0,0,0}; 
	for (int k{0}; k<z; ++k){
		for (int j{0}; j<y; ++j){
			for (int i{0}; i<x; ++i){
				
				loc = {i, j, k};

				if (lattice_index(loc,y,z)-c_idx != 1){
					std::cerr << "Something is fucked in LATTICE creation." << std::endl;
					exit (EXIT_FAILURE);
				}
				c_idx = lattice_index (loc,y,z); 

				Particle* p_ptr = new Particle ( loc, "s1", rng_uniform (0, 25) ); 

				// std::cout << "loc is "; print(loc);
				// std::cout << "index in lattice is " << lattice_index (loc, y, z) << std::endl;

				(*LATTICE).insert( (*LATTICE).begin() + lattice_index(loc, y, z), p_ptr) ;
				

			}
		}
	}

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of AddSolvent
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of AddCosolvent
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void AddCosolvent (std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, double frac, int Nmonomer, int x, int y, int z){

	int nsol2         = std::floor ((x*y*z-Nmonomer)*frac); 
	std::cout << "Number of particles of cosolvent is " << nsol2 << "." << std::endl;

	std::vector <int> indices (x*y*z);
	std::iota (indices.begin(), indices.end(), 0); 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle ( indices.begin(), indices.end(), std::default_random_engine(seed) );
    int count = 0;
    int i = 0; 
	while ( count < nsol2 ) {

		if ( (*LATTICE).at( indices[i] )->ptype[0] == 'm' ){
            ; //  count += 1; 
		}
		else {
			(*LATTICE).at( indices[i] )->ptype = "s2"; 
			(*Cosolvent).push_back((*LATTICE).at( indices[i] )); 
            count += 1; 
		}
        i += 1; 
			
	}
	
	return; 

}



// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of AddCosolvent
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: MetropoloisAcceptance  
//
// PARAMETERS: Initial energy, final energy, and kT. The second function also takes a weight. 
//
// WHAT THE FUNCTION DOES: calculates metropolis acceptance based on energies and kT
// 
// DEPENDENCIES: random number generator. 
// 
// THE CODE: 

bool MetropolisAcceptance(double E1, double E2, double kT){

	double dE = E2-E1; 
	double prob = std::exp(-1/kT*dE); 
	double r = rng_uniform(0.0, 1.0); 
	// std::cout << "Probability is " << prob <<"." << std::endl;

	// std::cout << "E1 is " << E1 << std::endl;
	// std::cout << "E2 is " << E2 << std::endl;

	// std::cout << "Probability of acceptance is " << prob << "." << std::endl;
	// std::cout << "RNG is " << r << "." << std::endl;
	if (r < prob){
		return true; 
	}
	else {
		return false; 
	}

}

bool MetropolisAcceptance(double E1, double E2, double kT, double rweight){

	double dE = E2-E1; 
	double prob = std::exp(-1/kT*dE) * rweight; 
	double r = rng_uniform(0.0, 1.0); 
	// std::cout << "Probability is " << prob <<"." << std::endl;
    // std::cout << "rweight is " << rweight << "." << std::endl;
	// std::cout << "E1 is " << E1 << std::endl;
	// std::cout << "E2 is " << E2 << std::endl;
	// std::cout << "Probability of acceptance is " << prob << "." << std::endl;
	// std::cout << "RNG is " << r << "." << std::endl;
	if (r < prob){
		return true; 
	}
	else {
		return false; 
	}

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of MetropoliAcceptance 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: InputParser
//
// PARAMETERS: All the inputs provided by the user. 
//
// WHAT THE FUNCTION DOES: Inputs are collected and examined. If there is a bad input, this function will 
// exit and raise an error. Important function. 
// 
// DEPENDENCIES: none other than the STL. 
// 
// THE CODE: 

void InputParser(int dfreq, int max_iter, bool r, \
	std::string positions, std::string topology, std::string dfile, \
	std::string efile, std::string mfile, std::string stats_file, \
	std::string lattice_file_read){

	
	if (!r) {

	    if (dfreq == -1 || max_iter == -1) {
	        std::cerr << "ERROR: No value for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
	        exit (EXIT_FAILURE);
	    }
	    
	    if ( positions== "__blank__" || topology == "__blank__" || dfile== "__blank__" || efile == "__blank__" || \
	    	mfile == "__blank__" || stats_file == "__blank__" ){
	        std::cerr << "polymer coords file is " << positions <<",\ntopology is " << topology <<",\npolymer coordinate dump file is " << dfile << ",\nenergy dump file is " \
	        << efile << ",\norientation file is " << mfile << ",\nmove statistics file is " << stats_file << "." << std::endl;
	        std::cerr << "ERROR: No value for option p (polymer coordinate file) and/or\nfor option S (solvent coordinate file) and/or\n" <<
	        "for option t (energy and geometry file) and/or\nfor option o (name of output dump file) and/or\nfor option e (name of orientation file) and/or\n" <<
	        "for option s (name of move stats file) and/or\n for option u (name of energy dump file) was provided. Exiting..." << std::endl;
	        exit (EXIT_FAILURE);    
	    }
	    
	    // set up these files 
	    std::ofstream polymer_dump_file (dfile);
	    std::ofstream energy_dump_file (efile);
	    std::ofstream orientation_dump_file (mfile); 
	    std::ofstream statistics_dump_file (stats_file); 

	    if ( lattice_file_read != "__blank__" ){
	    	std::cerr << "Restart has not been requested. Do not provide a restart file to read. Exiting..." << std::endl;
	    	exit (EXIT_FAILURE);
	    } 
	    

	}     

	else {

		if ( dfreq == -1 || max_iter == -1 ){
			std::cerr << "ERROR: No value for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
	        exit (EXIT_FAILURE);	
		}

		if ( positions== "__blank__" || topology == "__blank__" || dfile== "__blank__" || efile == "__blank__" || \
			mfile == "__blank__" || stats_file == "__blank__" || lattice_file_read == "__blank__" ) {
			std::cerr << "polymer coords file is " << positions <<",\ntopology is " << topology <<",\npolymer coordinate dump file is " << dfile << ",\nenergy dump file is " \
	        << efile << ",\norientation file is " << mfile << ",\nmove statistics file is " << stats_file << ", " << \
	        "\nlattice file to read is " << lattice_file_read << "." << std::endl;
	        std::cerr << "ERROR: No value for option p (polymer coordinate file) and/or\n" << 
	        "for option t (energy and geometry file) and/or\nfor option o (name of output dump file) and/or\nfor option e (name of orientation file) and/or\n" <<
	        "for option s (name of move stats file) and/or\n for option u (name of energy dump file) was provided. Exiting..." << std::endl;
	        exit (EXIT_FAILURE);    

		}

	}


	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of InputParser  
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: checkValidityOfCoords
//
// PARAMETERS: std::vector <int> v, some attributes of the Grid Object
// WHAT THE FUNCTION DOES: It look at a singular location coordinate and checks if it is a "good" location, 
// and does not lie outside the box, or something unphysical like that. 

// DEPENDENCIES: No custom functions required apart from those defined previously in object Grid. 
//
// THE CODE: 

bool checkValidityOfCoords(std::array <int,3> v, int x, int y, int z){
    if (v.at(0)> x || v.at(0) < 0 || v.at(1)> y || v.at(1)<0 || v.at(2)> z || v.at(2)<0){
        return false;
    }
    else {
        return true;
    }
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of checkValidityOfCoords 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: checkForOverlaps
//
// PARAMETERS: std::vector <int> PolymerVector
// 
// WHAT THE FUNCTION DOES: It looks at a vector of polymers that are supposed to go into the Grid. If
// two (or more) monomer units occupy the same spot, it will freak out and exit out of compilation. 
// 
// DEPENDENCIES: No custom functions required apart from those defined previously in object Grid. 
//
// THE CODE: 


bool checkForOverlaps(std::vector <Polymer> Polymers){
    
    std::vector <std::array <int,3>> loc_list; 

    for (Polymer& pmer: Polymers){
        for (Particle*& p: pmer.chain){
            // check if element exists in vector 
                if (std::find(loc_list.begin(), loc_list.end(), p->coords) != loc_list.end() ){
                    std::cerr << "you have a repeated element." << std::endl;
                    // std::cout << "current element is: " << std::endl;
                    print(p->coords); 
                    print(loc_list);
                    return false; 
                    }
            
                else{
                    loc_list.push_back(p->coords);  
                }
            }
        }    
    
    std::cout << "Input file has no overlaps!" << std::endl;
    return true;

}

bool checkForOverlaps ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE){

	bool pmer_loop_flag = false; 
	bool particle_found_flag = false; 
	for ( int i{0}; i< static_cast<int>((*LATTICE).size()); ++i ){
		
		pmer_loop_flag = false; 
		particle_found_flag = false; 

		if ( (*LATTICE)[i]->ptype[0] == 'm' ){

			for ( Polymer& pmer: (*Polymers) ){
				for (Particle*& p: pmer.chain){
					if ( (*LATTICE)[i]->coords == p->coords){
						pmer_loop_flag = true; 
						break;
					}
				}
				if (pmer_loop_flag) {
					particle_found_flag = true; 
					break;
				}
			}

			if (particle_found_flag){
				continue;
			}
			else {
				std::cout << "Something is fucked. There is a random monomer floating around." << std::endl;
				exit (EXIT_FAILURE); 
			}
		}

	}

	return true;

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of checkForOverlaps. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: checkForSolventMonomerOverlap
//
// PARAMETERS: std::vector <int> PolymerVector
// 
// WHAT THE FUNCTION DOES: It looks at a vector of polymers that are supposed to go into the Grid. If
// two (or more) monomer units occupy the same spot, it will freak out and exit out of compilation. 
// 
// DEPENDENCIES: No custom functions required apart from those defined previously in object Grid. 
//
// THE CODE: 


bool checkForSolventMonomerOverlap(std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int y, int z){
    
    std::vector <std::array <int,3>> loc_list;

    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){
            // check if element exists in vector 
                if (std::find(loc_list.begin(), loc_list.end(), p->coords) != loc_list.end() ){
                    std::cerr << "you have a repeated element." << std::endl;
                    return false; 
                }
            
                else{
                    loc_list.push_back(p->coords);  
                }
        }
    }    
    
    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){

        	if ( (*LATTICE)[ lattice_index (p->coords, y, z)]->coords == p->coords ) {

        		if ( (*LATTICE)[ lattice_index(p->coords, y, z)]->ptype[0] == 's' ){
        			std::cerr << "Some kind of bad solvent-monomer overlap that has taken a place. A monomer is being represented by a solvent. Something's fucked." << std::endl;
        			std::cerr << "Location is: "; print (p->coords); 
        			std::cerr << "Type is: " << ((*LATTICE)[ lattice_index (p->coords, y, z)]->ptype) << std::endl; 
        			std::cerr << "Type is: " << p->ptype << std::endl; 
        			std::cerr << "Location on lattice is: "; print ((*LATTICE)[ lattice_index (p->coords, y, z)]->coords); 
        			exit (EXIT_FAILURE); 
        		}
    			else {
    				continue;
    			}

        	}

        	else {
        		std::cerr << "Something is wrong with the LATTICE map. Monomer unit lost. Something is fucked. " << std::endl;
        		std::cerr << "Output from *Polymers: "; print (p->coords);
        		std::cerr << "Output from *LATTICE: "; print ((*LATTICE)[ lattice_index (p->coords, y, z)]->coords);
        		exit(EXIT_FAILURE);
        	}

        }
    }

    std::cout << "Input file has no overlap between and internally amongst solvent and monomers!" << std::endl;
    return true;

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of checkForOverlaps. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: checkConnectivity
//
// PARAMETERS: std::vector <int> PolymerVector
// 
// WHAT THE FUNCTION DOES: It looks at a vector of polymers that are supposed to go into the Grid. If
// two adjacent monomer units in the polymer vector do not have adjacent coordinates, you have bad connectivity.  
// 
// DEPENDENCIES: subtract_vectors 
//
// THE CODE: 

bool checkConnectivity(std::vector <Polymer> Polymers, int x, int y, int z) {
    
    for (Polymer& pmer: Polymers){
        size_t length = pmer.chain.size(); 
        std::array <int,3> connection = {0,0,0}; 
        std::sort (adrns.begin(), adrns.end() ); 

        for (int i{1}; i<static_cast<int>(length); ++i){
            
            connection = subtract_arrays(&(pmer.chain[i]->coords), &(pmer.chain[i-1]->coords));
            impose_pbc(&connection, x, y, z);
            modified_direction ( &connection, x, y, z); 

            if ( binary_search ( adrns.begin(), adrns.end(), connection) ) {
            	continue;
            }
            else {
            	std::cerr << "Shit, you have bad connectivity inside one (or maybe more) polymers. Check input file." << std::endl;
                return false; 
            }

        }
    }

    std::cout << "Input polymers are well-connected!" << std::endl;
    return true;
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of checkConnectivity. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: checkPointersOnLattice
//
// PARAMETERS: LATTICE
// 
// WHAT THE FUNCTION DOES: It looks at each point on the lattice and makes sure that every lattice point has 
// a pointer with the same location as the lattice site
// 
// DEPENDENCIES: 
//
// THE CODE: 

bool checkPointersOnLattice (std::vector <Particle*>* LATTICE, int x, int y, int z){
	
	for ( int i{0}; i<x*y*z; ++i ) {

		if ( location(i, x, y, z) == (*LATTICE)[i]->coords ) {
			continue;
		}
		else {
			std::cout << "Problem with pointers." << std::endl;
			std::cout << "Bad location is "; print ( location (i, x, y, z) ); 
			std::cout << "LATTICE says "; print ((*LATTICE)[i]->coords);
			std::cerr << "Something is fucked. Pointer does not correspond to position. " << std::endl;
			exit (EXIT_FAILURE);
			return false;
		}

	}

	return true;

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of checkPointersOfLattice. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
/*
bool checkSolvationShells (std::vector <Polymer>* Polymers,  int y, int z){

	for ( Polymer& pmer: *Polymers){
		for ( Particle*& p: pmer.chain ){
			if ( std::find ( solvation_shells->begin(), solvation_shells->end(), lattice_index (p->coords, y, z ) ) != solvation_shells->end() ) {
				std::cerr << "Solvation shell has the index of monomer particle. Something's fucked." << std::endl; 
				exit (EXIT_FAILURE);
			}
		}
	}
	return true; 
}
*/

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of checksolvationshells. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: CheckStructures 
//
// PARAMETERS: LATTICE, Polymers
// 
// WHAT THE FUNCTION DOES: Grand function which runs all checks in one hit. 
//
// DEPENDENCIES: checkforoverlaps, checkconnectivity, checkpointers, checkoverlaps -- all the check functions 
//
// THE CODE: 


void CheckStructures (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector<Particle*>* LATTICE,  int x, int y, int z){

	std::cout << "Checking validity of coords..." << std::endl;
    std::cout << "checkForOverlaps says: " << checkForOverlaps(*Polymers) << "." << std::endl; 
    if (!checkForOverlaps(*Polymers)){
        std::cout << "Something is fucked up overlaps-wise in the polymer itself." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (! checkForOverlaps (Polymers, LATTICE)){
    	std::cout << "Random monomer floating in the lattice! Breaking out..." << std::endl; 
        exit (EXIT_FAILURE);
    }
    std::cout << "No random monomers floating around!!" << std::endl;

    if (!checkForSolventMonomerOverlap (Polymers, LATTICE, y, z) ){
    	std::cout << "Something is fucked up solvent-monomer overlaps-wise. " << std::endl; 
        exit(EXIT_FAILURE);
    }

    std::cout << "checkConnectivity says: " << checkConnectivity(*Polymers, x, y, z) << "." << std::endl;
    if (!checkConnectivity(*Polymers, x, y, z) ){
    	std::cout << "Something is fucked up connectivity-wise." << std::endl; 
        exit(EXIT_FAILURE);
    }

    if (checkPointersOnLattice (LATTICE, x, y, z) ){
        std::cout << "We good. LATTICE is in good shape." << std::endl; 
    }
    else {
        std::cerr <<"Something is fucked with pointers on LATTICE." << std::endl; 
        exit (EXIT_FAILURE);
    }


    for ( Particle*& p: *Cosolvent ){

    	if ( (*LATTICE).at( lattice_index(p->coords, y, z) )->ptype != "s2" && (*LATTICE).at( lattice_index(p->coords, y, z) )->ptype != p->ptype && (*LATTICE).at( lattice_index(p->coords, y, z) )->orientation != p->orientation && \
    	(*LATTICE).at( lattice_index(p->coords, y, z) )->coords != p->coords ) {
    		std::cerr << "Cosolvent and LATTICE do not align. Something's fucked." << std::endl;
    		exit (EXIT_FAILURE);
    	}

    }
    std::cout << "Cosolvent vector looks good!" << std::endl;

    return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of CheckStructures
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: ExtractPolymersFromFile
//
// PARAMETERS: std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the positions file, and extracts all the data from the file 
// in the form of a vector of strings.  
// 
// DEPENDENCIES: ExtractContentFromFile, makePolymer, checkValidityOfCoords 
//
// THE CODE: 


std::vector<Polymer> ExtractPolymersFromFile(std::string filename, int x, int y, int z){

    int NumberOfPolymers = ExtractNumberOfPolymers(filename); 

    std::vector <Polymer> PolymerVector; 
    PolymerVector.reserve(NumberOfPolymers);

    std::vector <std::array <int,3>> locations; 

    std::vector <std::string> contents = ExtractContentFromFile(filename); // this extracts every line of the file

    std::regex start ("START"), end ("END"); 


    int startCount{0}, endCount {0}; 
    
    
    for (std::string s: contents){
        
        
         
        std::stringstream ss(s); 
        if (std::regex_search(s, start) ){
            ++startCount;
            continue; 
        }

        else if (std::regex_search(s, end) ) {
            ++endCount;
            
            Polymer pmer = makePolymer(locations);
            PolymerVector.push_back(pmer);
            
            locations.clear();
            
        }

        else{
            std::array <int,3> loc;
            int j{0};  
            for (int i=0; ss>>i; ){
                
                loc[j] = i;
                ++j;

            }

            if (! checkValidityOfCoords(loc, x, y, z)){
            std::cerr << "Coordinates are out of bounds. Bad input file." << std::endl;
            exit(EXIT_FAILURE); 
            }
        
            locations.push_back(loc); 
            
        
        }
    }
    
    

    if(!(checkForOverlaps(PolymerVector))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Overlap detected." << std::endl;
        exit(EXIT_FAILURE);  
    }

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // throw in a check for connectivity of polymer chains 
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    if (!(checkConnectivity(PolymerVector, x, y, z))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Monomer units are not adjacent to one another on the lattice." << std::endl;
        exit(EXIT_FAILURE); 
    }
    
    

    return PolymerVector; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of ExtractPolymersFromFile. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
// 
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: MonomerReporter
//
// PARAMETERS: Polymers, to_check
//
// WHAT THE FUNCTION DOES: Calculates energy of the current Grid. Critical to correctly evolve system. 
// includes nearest neighbor interactions with and without directional effects.  
//
// DEPENDENCIES: obtain_ne_list 
//
// OPTIMIZATION OPPORTUNITY: I am double counting monomer-monomer interactions. This can possibly be avoided. 
//
// THE CODE: 

bool MonomerReporter (std::vector <Polymer>* Polymers, std::array <int,3>* to_check){


	for (Polymer& pmer: (*Polymers)) {
		for (Particle*& p: pmer.chain){ 

			if ((*to_check) == p->coords){
				return true;
			}
		}
	}

	return false; 

}

bool MonomerReporter (std::vector <Particle*>* LATTICE, std::array<int,3>* to_check, int y, int z){

	if ( (*LATTICE)[ lattice_index((*to_check), y, z) ]->ptype[0] == 'm' ){
		return true;
	}
	else {
		return false; 
	}

}

bool MonomerReporter (std::vector <Polymer>* Polymers, std::array <int,3>* check_1, std::array <int,3>* check_2){


	for (Polymer& pmer: (*Polymers)) {
      		for (Particle*& p: pmer.chain){

      			if (p->coords == (*check_1) || p->coords == (*check_2) ){	 
      				return true; 
      			}
      		}
      	}
    return false; 
}


bool MonomerReporter (std::vector <Particle*>* LATTICE, std::array <int,3>* check_1, std::array <int,3>* check_2, int y, int z){

	if ( (*LATTICE)[lattice_index((*check_1), y, z)]->ptype[0] == 'm' || (*LATTICE)[lattice_index((*check_2), y, z)]->ptype[0] == 'm'){
		return true;
	}

	return false; 

}

bool MonomerNeighborReporter ( std::vector <Polymer>* Polymers, std::array <int,3>* to_check, int x, int y, int z){

	std::array <int,3> diff = {0,0,0}; 
	for (Polymer& pmer: (*Polymers) ){
		for (Particle*& p: pmer.chain) {

			diff = subtract_arrays ( &(p->coords), to_check ); 
			modified_direction (&diff, x, y, z); 

			if ( std::find ( adrns.begin(), adrns.end(), diff ) != adrns.end() ) {
				return true;
			}

		} // particle loop 
	} // polymer loop 

	return false; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of MonomerReporter 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
// 
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: CalculateEnergy
//
// PARAMETERS: PolymerVector, SolventVector, x, y, z
// WHAT THE FUNCTION DOES: Calculates energy of the current Grid. Critical to correctly evolve system. 
// includes nearest neighbor interactions with and without directional effects.  
//
// DEPENDENCIES: obtain_ne_list 
//
// OPTIMIZATION OPPORTUNITY: I am double counting monomer-monomer interactions. This can possibly be avoided. 
//
// THE CODE: 

double CalculateEnergy (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array<double,8>* E, std::array<double,8>* contacts, int x, int y, int z) {
    
    double Energy {0.0};
    (*contacts) = {0,0,0,0,0,0,0,0}; 

    double             dot_product   = -2; 
    std::array <int,3> connvec       = {0,0,0}; 
    double             theta_1       = 0; 
    double             theta_2       = 0; 
    double             magnitude     = 0; 

    std::array <std::array <int,3>, 26> ne_list; 

    // run energy computations for every monomer bead 
    // m-m  = stacking interaction
    // m-s1 = stacking interaction 
    // m-s2 = isotropic interaction 

    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){
        	// std::cout << "ploc = "; print (p->coords); 
            ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list 
            
            // std::cout << "Particle loc is "; print (p->coords); 

            for ( std::array <int, 3>& loc: ne_list){
            	// std::cout << "l index = " << lattice_index (loc, y, z) << ", "; print(loc);
            	dot_product = take_dot_product (  p->orientation, (*LATTICE)[ lattice_index(loc, y, z) ]->orientation );
            	if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "m1"){
            		// m-m interactions

            		if (dot_product > 0.54){
            			Energy += 0.5* (*E)[0];
            			(*contacts)[0]   += 0.5;
            		}
            		else {
            			Energy += 0.5* (*E)[1];
            			(*contacts)[1]  += 0.5;
            		}
            	}
            	else if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s1" ){ 
            		// m-s1 interactions 
            		if (dot_product > 0.54){
            			Energy += (*E)[2];
            			(*contacts)[2] += 1;
            		}
            		else {
            			Energy += (*E)[3]; 
            			(*contacts)[3]  += 1;
            		}
            	}

            	else {
					// m-s2 interactions 
            		Energy += (*E)[4]; 
            		(*contacts)[4] += 1; 
            		(*contacts)[5] += 1; 

            	}
            }
        }
    }
    

    for ( Particle*& p: *Cosolvent ){

    	ne_list = obtain_ne_list ( p->coords, x, y, z ); 
    	for ( std::array <int,3>& loc: ne_list ){

    		if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "m1" || (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s2"){
    			continue; 
    		}
    		else {

    			connvec   = subtract_arrays ( &(*LATTICE)[lattice_index(loc, y, z)]->coords, &(p->coords) );
    			modified_direction ( &connvec, x, y, z); 
    			magnitude = std::sqrt (connvec[0]*connvec[0]+connvec[1]*connvec[1]+connvec[2]*connvec[2]); 
    			theta_1 = std::acos (take_dot_product ( Dir2Or[scale_arrays(1/magnitude, &connvec)], p->orientation ) ); 
    			theta_2 = std::acos (take_dot_product ( Dir2Or[scale_arrays(-1/magnitude, &connvec)], (*LATTICE)[lattice_index(loc, y, z)]->orientation ) ); 

    			if ( theta_1+theta_2 > M_PI/2 ){
					Energy += (*E)[7]; 
					(*contacts)[7] += 1; 
				}
				else {
					// locked or aligned
					Energy += (*E)[6]; 
					(*contacts)[6] += 1; 
				}

    		}

    	}

    }
    
    return Energy; 
}



double CalculateEnergy_parallel (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array<double,8>* E, std::array<double,8>* contacts, int x, int y, int z) {
    
    double Energy {0.0};
    (*contacts) = {0,0,0,0,0,0,0,0}; 
    int NCosolvent = static_cast<int>((*Cosolvent).size()); 
    /*
    double             dot_product   = -2; 
    // std::array <int,3> connvec       = {0,0,0}; 
    // double             theta_1       = 0; 
    // double             theta_2       = 0; 
    // double             magnitude     = 0; 

    std::array <std::array <int,3>, 26> ne_list; 

    // run energy computations for every monomer bead 
    // m-m  = stacking interaction
    // m-s1 = stacking interaction 
    // m-s2 = isotropic interaction 

    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){
        	// std::cout << "ploc = "; print (p->coords); 
            ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list 
            
            // std::cout << "Particle loc is "; print (p->coords); 

            for ( std::array <int, 3>& loc: ne_list){
            	// std::cout << "l index = " << lattice_index (loc, y, z) << ", "; print(loc);
            	dot_product = take_dot_product (  p->orientation, (*LATTICE)[ lattice_index(loc, y, z) ]->orientation );
            	if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "m1"){
            		// m-m interactions

            		if (dot_product > 0.54){
            			Energy += 0.5* (*E)[0];
            			(*contacts)[0]   += 0.5;
            		}
            		else {
            			Energy += 0.5* (*E)[1];
            			(*contacts)[1]  += 0.5;
            		}
            	}
            	else if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s1" ){ 
            		// m-s1 interactions 
            		if (dot_product > 0.54){
            			Energy += (*E)[2];
            			(*contacts)[2] += 1;
            		}
            		else {
            			Energy += (*E)[3]; 
            			(*contacts)[3]  += 1;
            		}
            	}

            	else {
					// m-s2 interactions 
            		Energy += (*E)[4]; 
            		(*contacts)[4] += 1; 
            		(*contacts)[5] += 1; 

            	}
            }
        }
    }
    */     
    // start the set up for parallelized energy computation ... 
    
    omp_set_num_threads (NUM_THREADS); 
    // these need to be updated for the whole thing to work! 
    double energy_chunks[NUM_THREADS]; 
    double contacts_aligned [NUM_THREADS]; 
    double contacts_naligned [NUM_THREADS]; 
    std::array <double,8> c_contacts = *contacts; 

    #pragma omp parallel 
    { // start of pragma
    int id, dummy_idx, nthrds;
    id = omp_get_thread_num();  
    nthrds = omp_get_num_threads();  
    
    // reinitializing a bunch of new variables for locking interaction 
    std::array <int, 3> connvec_p; 
    double magnitude_p {0}; 
    double theta_1_p   {0}; 
    double theta_2_p   {0}; 
    std::array <double,8> E_ = *E;  
    double dot_product_p   = -2; 
    std::array <std::array <int,3>, 26> ne_list; 
    
    if (id == 0) {
    // do the polymer loop 
        energy_chunks [0] = 0.0; 
        for ( Polymer& pmer: (*Polymers) ) {
            for ( Particle*& p: pmer.chain ) {
                ne_list = obtain_ne_list(p->coords, x, y, z); 
                for ( std::array <int,3>& loc: ne_list ) {
            	    dot_product_p = take_dot_product (  p->orientation, (*LATTICE)[ lattice_index(loc, y, z) ]->orientation );
            	    if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "m1"){
            		    // m-m interactions

            		    if (dot_product_p > 0.54){
            			    energy_chunks[0] += 0.5* (E_)[0];
            			    (c_contacts)[0]   += 0.5;
            		    }
            		    else {
            			    energy_chunks[0] += 0.5* (E_)[1];
            			    (c_contacts)[1]  += 0.5;
            		    }
            	    }
            	    else if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s1" ){ 
            		    // m-s1 interactions 
            		    if (dot_product_p > 0.54){
            			    energy_chunks[0] += E_[2];
            			    (c_contacts)[2] += 1;
            		    }
            		    else {
            			    energy_chunks[0] += E_[3]; 
            			    (c_contacts)[3]  += 1;
            		    }
            	    }

            	    else {
					    // m-s2 interactions 
            		    energy_chunks[0] += (E_)[4]; 
            		    (c_contacts)[4] += 1; 
            		    (c_contacts)[5] += 1; 
            	    }
                }
            }
        }
    }

    else {
        // end of reinitialization
            // std::cout << id << std::endl;
        for (dummy_idx=id-1, energy_chunks[id]=0.0, contacts_naligned[id-1]=0, contacts_aligned[id-1]=0; dummy_idx<NCosolvent; dummy_idx = dummy_idx+nthrds-1) {
            ne_list = obtain_ne_list ((*Cosolvent)[dummy_idx]->coords, x, y, z);
            for ( std::array <int,3>& loc: ne_list ) {
            
                if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "m1" || (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s2" ) {
                    continue; 
                }
                else {
    			    connvec_p   = subtract_arrays ( &(*LATTICE)[lattice_index(loc, y, z)]->coords, &( (*Cosolvent)[dummy_idx]->coords) );
    			    modified_direction ( &connvec_p, x, y, z); 
    			    magnitude_p = std::sqrt (connvec_p[0]*connvec_p[0]+connvec_p[1]*connvec_p[1]+connvec_p[2]*connvec_p[2]); 
    			    theta_1_p = std::acos (take_dot_product ( Dir2Or[scale_arrays(1/magnitude_p, &connvec_p)] , (*Cosolvent)[dummy_idx]->orientation ) ); 
    			    theta_2_p = std::acos (take_dot_product ( Dir2Or[scale_arrays(-1/magnitude_p, &connvec_p)], (*LATTICE)[lattice_index(loc, y, z)]->orientation ) ); 
                
                    if ( theta_1_p + theta_2_p > M_PI/2 ) {
                        energy_chunks [id] += E_[7]; // naligned_energy;
                        contacts_naligned [id-1] += 1;
                    }
                    else {
                        energy_chunks [id] += E_[6]; // aligned_energy;
                        contacts_aligned [id-1] += 1; 
                    }
                }
            }
        }
    }
    } // end of pragma
    for (int i{0}; i<6; ++i) {
        (*contacts)[i] = c_contacts[i]; 
    }
    
    for ( int i{0}; i<NUM_THREADS; ++i) {
        Energy += energy_chunks[i];
    }
    for (int i{0}; i<NUM_THREADS-1; ++i){
        (*contacts)[6] += contacts_aligned[i];
        (*contacts)[7] += contacts_naligned[i];
    }
    
    return Energy; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of calculateEnergy. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: dumpPositionsOfPolymers 
//
// PARAMETERS: (int step, std::string filename), and some attributes present in the Grid Object 
// 'step' is the current time step we are at. This is an integer which is likely defined in the driver code.  
// 'filename' is the file to which I am going to print out coordinate information about the polymer. 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers in the current Grid object \, and prints them 
// out to a text file.  
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 


void dumpPositionsOfPolymers (std::vector <Polymer>* Polymers, int step, std::string filename){
    // std::ostringstream os; 
    std::ofstream dump_file(filename, std::ios::app); 
    dump_file <<"Dumping coordinates at step " << step << ".\n";
    
 
    // dump_file <<"Dumping coordinates at step " << step << ".\n";
    int count = 0; 
    for (Polymer& pmer: *Polymers){
        
        dump_file <<"Dumping coordinates of Polymer # " << count << ".\n";
        dump_file<<"START" << "\n";
        for (Particle*& p: pmer.chain){
            for (int i: p->coords){
                dump_file << i << " | "; 
            }
            dump_file << p->orientation << " | ";
            dump_file << "\n"; 
        }
        ++count ; 
        dump_file <<"END"<<"\n";
    }
    
    dump_file << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n";

    return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of dumpPositionsOfPolymer. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: dumpPositionsOfSolvent 
//
// PARAMETERS: (int step, std::string filename), and some attributes present in the Grid Object 
// 'step' is the current time step we are at. This is an integer which is likely defined in the driver code.  
// 'filename' is the file to which I am going to print out coordinate information about the polymer. 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers in the current Grid object \, and prints them 
// out to a text file.  
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 


void dumpPositionOfSolvent(std::vector <Particle*>* LATTICE, int step, std::string filename){

	std::ofstream dump_file(filename, std::ios::app); 
    dump_file <<"Dumping coordinates at step " << step << ".\n";
    
    // dump_file <<"Dumping coordinates at step " << step << ".\n";
    int count = 0; 
        
    dump_file <<"Dumping coordinates of solvent # " << count << ".\n";
    dump_file<<"START" << "\n";
    
    for (Particle*& p: *LATTICE){
    	// std::cout << "ptype is " << (*p).ptype << std::endl;
    	if (p->ptype[0] == 's'){
    		dump_file<<"Orientation: " << p->orientation <<", ";
        	for (int i: p->coords){
            	dump_file << i << " | "; 
        	}
        	dump_file << "\n"; 
    	}
    }
    ++count ; 
    dump_file <<"END"<<"\n";
    
    
    dump_file << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n";
    dump_file.close();

    return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of dumpPositionOfSolvent. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//============================================================
//============================================================
//
// NAME OF FUNCTION: dumpEnergy
//
// PARAMETERS: (int step, std::string filename), and some attributes present in the Grid Object 
// 'step' is the current time step we are at. This is an integer which is likely defined in the driver code.  
// 'filename' is the file to which I am going to print out coordinate information about the polymer. 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers in the current Grid object \, and prints them 
// out to a text file.  
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 


void dumpEnergy (double sysEnergy, int step, std::array<double,8>* contacts, std::string filename) {
    std::ofstream dump_file(filename, std::ios::app); 
    // std::ostringstream os; 
    
    // (mm_a + mm_n)     0 = mm_a  , 1 = mm_n
    // (ms1_a+ms1_n)     2 = ms1_a , 3 = ms1_n 
    // (ms2_a+ms2_n)     4 = ms2_a , 5 = ms2_n 
    // (s1s2_a+s1s2_n)   6 = s1s2_a, 7 = s1s2_n

    dump_file << sysEnergy << " | " \
			<< (*contacts)[0]+(*contacts)[1] << " | " << (*contacts)[0] << " | " << (*contacts)[1] << " | " \
			<< (*contacts)[2]+(*contacts)[3] << " | " << (*contacts)[2] << " | " << (*contacts)[3] << " | " \
			<< (*contacts)[4]+(*contacts)[5] << " | " << (*contacts)[4] << " | " << (*contacts)[5] << " | " \
			<< (*contacts)[6]+(*contacts)[7] << " | " << (*contacts)[6] << " | " << (*contacts)[7] << " | " << step << "\n";

    return; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of dumpEnergy
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//============================================================
//============================================================
//
// NAME OF FUNCTION: dumpMoveStatistics 
//
// PARAMETERS: (int step, std::string filename), and some attributes present in the Grid Object 
// 'step' is the current time step we are at. This is an integer which is likely defined in the driver code.  
// 'filename' is the file to which I am going to print out coordinate information about the polymer. 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers in the current Grid object \, and prints them 
// out to a text file.  
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 

void dumpMoveStatistics (std::array <int,9>* attempts, std::array <int,9>* acceptances, int step, std::string stats_file) {
    
    std::ofstream dump_file (stats_file, std::ios::out); 
    dump_file << "At step " << step << "...\n";

    dump_file << "End rotations without bias         - attempts: " << (*attempts)[0] <<", acceptances: " << (*acceptances)[0] << ", acceptance fraction: " << static_cast<double>((*acceptances)[0])/static_cast<double>((*attempts)[0]) << std::endl; 
    dump_file << "Reptation without bias             - attempts: " << (*attempts)[1] <<", acceptances: " << (*acceptances)[1] << ", acceptance fraction: " << static_cast<double>((*acceptances)[1])/static_cast<double>((*attempts)[1]) << std::endl; 
    dump_file << "Chain regrowth with overlap bias   - attempts: " << (*attempts)[2] <<", acceptances: " << (*acceptances)[2] << ", acceptance fraction: " << static_cast<double>((*acceptances)[2])/static_cast<double>((*attempts)[2]) << std::endl; 
    dump_file << "Chain regrowth with ori flip       - attempts: " << (*attempts)[3] <<", acceptances: " << (*acceptances)[2] << ", acceptance fraction: " << static_cast<double>((*acceptances)[3])/static_cast<double>((*attempts)[3]) << std::endl; 
    dump_file << "Solvent flips without bias         - attempts: " << (*attempts)[4] <<", acceptances: " << (*acceptances)[3] << ", acceptance fraction: " << static_cast<double>((*acceptances)[4])/static_cast<double>((*attempts)[4]) << std::endl;
    dump_file << "Solvation shell flip with bias     - attempts: " << (*attempts)[5] <<", acceptances: " << (*acceptances)[5] << ", acceptance fraction: " << static_cast<double>((*acceptances)[5])/static_cast<double>((*attempts)[5]) << std::endl;
    dump_file << "Polymer flips                      - attempts: " << (*attempts)[6] <<", acceptances: " << (*acceptances)[6] << ", acceptance fraction: " << static_cast<double>((*acceptances)[6])/static_cast<double>((*attempts)[6]) << std::endl;
    dump_file << "Solvent exchange with bias         - attempts: " << (*attempts)[7] <<", acceptances: " << (*acceptances)[7] << ", acceptance fraction: " << static_cast<double>((*acceptances)[7])/static_cast<double>((*attempts)[7]) << std::endl;

    return;
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of dumpMoveStatistics
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//============================================================
//============================================================
//
// NAME OF FUNCTION: dumpOrientations 
//
// PARAMETERS: (std::vector <Polymer>* PV, std::vector <Particle>* SV, int step, std::string filename), and some attributes present in the Grid Object 
// 'step' is the current time step we are at. This is an integer which is likely defined in the driver code.  
// 'filename' is the file to which I am going to print out coordinate information about the polymer. 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers and solvents and throws out neighboring orientations  
//   
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 
// 
// dumping orientations of solvent around monomer segments

void dumpOrientation( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int step, std::string filename, int x, int y, int z ) {
    // std::cout<< "just inside dumpO..."<<std::endl; 
    std::ofstream dump_file (filename, std::ios::app); 
    // std::cout << "step = " << step << ".\n";
    dump_file << "START for Step " << step << ".\n";
    std::vector <int> solvent_indices; 
    // std::pair <char, int> properties (' ' , -1); 
    for ( Polymer& pmer: (*Polymers) ) {
        for ( Particle*& p: pmer.chain ) {
            
            dump_file << p->orientation << " | ";
            std::array <std::array<int,3>, 26> ne_list = obtain_ne_list (p->coords, x, y, z) ;
            
            for ( std::array <int,3>& ne: ne_list) {
                
                // std::cout << "Reported~\n"; 
                
                if ( (*LATTICE)[ lattice_index(ne, y, z) ]->ptype[0] == 's' && std::find( solvent_indices.begin(), solvent_indices.end(), lattice_index(ne, y, z)) == solvent_indices.end()  ){
                    dump_file << ((*LATTICE)[ lattice_index(ne, y, z) ])->orientation << " | ";  
                } 
            }
            dump_file << "\n"; 
        } 
    }
    
    dump_file << "END. \n";
    return;

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of dumpOrientation 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//============================================================
//============================================================
//
// NAME OF FUNCTION: dumpLATTICE
//
// PARAMETERS: (std::vector <Polymer>* LATTICE, int step, std::string filename). 
//
// WHAT THE FUNCTION DOES: It takes the coordinates of the polymers and solvents and dumps all of it in a file.   
//
// DEPENDENCIES: No custom function required. All can be done from C++ STL. 
//
// THE CODE: 
// 

void dumpLATTICE ( std::vector <Particle*> *LATTICE, int step, int y, int z, std::string filename ){

	std::ofstream dump_file ( filename, std::ios::out ); 
	dump_file << "FINAL STEP: " << step << ".\n"; 
	for ( Particle*& p: (*LATTICE) ){
		dump_file << p->orientation << ", " << p->ptype << ", " << lattice_index(p->coords, y, z) << "\n"; 
	}

	dump_file << "END. \n";
	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of dumpLATTICE
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: TailRotation_UNBIASED and TailRotation_BIASED
//
// PARAMETERS: 
// 
// WHAT THE FUNCTION DOES: it will perform a rotation of the monomer at the first index (the tail)
// of the polymer. As it stands, this function only rotates one molecule at the tail. It does so in an unbiased fashion. 
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully.    
//
// THE CODE: 

void TailRotation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {


	std::array <double,8> c_contacts = *contacts;

    // get the neighborlist of particle at index 1 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[0]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[1]->coords;

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

	int choice = rng_uniform(0, 25); 

	// make the change on the lattice 
	// make the change only to the SOLVENT site... 

	if ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] == 's' ){

		(*Polymers)[index].chain[0]->coords = ne_list[choice]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 
		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[0]; 
	}
	else {
		*IMP_BOOL = false; 
		return; 
	}

	double energy_n = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z); 

	if ( MetropolisAcceptance (*sysEnergy, energy_n, temperature)){
		*sysEnergy        = energy_n; 
		*contacts         = c_contacts; 
	}

	else {

		// revert the polymer 
		(*Polymers)[index].chain[0]->coords = loc_0; 
		(*LATTICE) [ lattice_index (loc_0, y, z) ]->coords  = ne_list[choice];

		// do the switch 
		// take the pointer of the solvent and put it where the monomer was on the lattice 

		(*LATTICE) [ lattice_index(ne_list[choice], y, z) ] = (*LATTICE)[ lattice_index (loc_0, y, z) ];
		(*LATTICE) [ lattice_index(loc_0, y, z) ]           = (*Polymers)[index].chain[0];  
		*IMP_BOOL = false; 
	}
	
	return; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRotation_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void TailRotation_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

    // get the neighborlist of particle at index 1 

	std::array <double,8> c_contacts = *contacts; 

    std::array <int,3> loc_0 = (*Polymers)[index].chain[0]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[1]->coords;
    
    double w2 {0}; // acceptance criteria weights 
    
    // first check if tail rotation can be performed at all... 
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

	// find all spots that are available 
	std::vector <std::array <int,3>> idx_v; 
	idx_v.reserve(26); 

	for (std::array <int,3>& to_rot: ne_list){
		if ( to_rot == loc_0){
			continue; 
		}

		else if ( to_rot == (*Polymers)[index].chain[2]->coords ) {
			continue; 
		}

		else if ( (*LATTICE)[lattice_index(to_rot, y, z) ]->ptype[0] != 'm'  ){
			idx_v.push_back( to_rot );  
		}
	}

    if ( static_cast<int>(idx_v.size()) == 0 ){
    	// if nothing can be done, just get out. 
    	*IMP_BOOL = false;
    	return; 
    }

    // otherwise, let's fuckin go
	
	// prepare the initial boltzmann sampling... 
	// get the solvent neighbors in the current location 

    // save current state... 
    // initial location is loc_0 
    // instantiate a variable that contains all energies and contacts 

    std::vector <double> energies; 
    energies.reserve(26); 

    std::vector <std::array<double,8>> contacts_store; 
    contacts_store.reserve(26); 

    std::vector <std::vector <int>> solvation_shell_store; 
    solvation_shell_store.reserve(26); 

    for (std::array<int,3>& rot: idx_v) {

    	(*Polymers)[index].chain[0]->coords = rot;
    	(*LATTICE)[ lattice_index (rot, y, z) ]->coords = loc_0; 

    	// do the switch 
    	(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (rot, y, z) ]; 
		(*LATTICE)[ lattice_index (rot  , y, z) ]    = (*Polymers)[index].chain[0]; 

		// the LATTICE is in a new configuration! 
		// calculate energy of new configutation... 

		energies.push_back( CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z) );
		contacts_store.push_back(c_contacts);

		// now that I have the energy of the new configuration, time to revert back to the old configuration... 
		(*Polymers)[index].chain[0]->coords = loc_0; 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]->coords = rot;

		// undo the switch 
		(*LATTICE)[ lattice_index (rot   , y, z) ]    = (*LATTICE)[ lattice_index (loc_0, y, z) ]; 
		(*LATTICE)[ lattice_index (loc_0 , y, z) ]    = (*Polymers)[index].chain[0]; 

    }

    idx_v.push_back (loc_0);
    energies.push_back (*sysEnergy); 
    contacts_store.push_back (*contacts);

    // now that i have the energies of each config, i can calculate boltzmann weights 
    for ( double e: energies ){
    	w2 += std::exp (-1/temperature*e); 
    }

    // choose a location now that boltzmann weights have been calculated. 
    double r    = rng_uniform (0.0, 1.0); 
    double rsum = 0; 
    int count   = 0; 

    for ( double e: energies ){
    	rsum += std::exp (-1/temperature*e)/w2; 
    	if ( r < rsum ){
    		break;
    	}
    	count += 1; 
    }

    // get system to new state 
    (*Polymers)[ index ].chain[0]->coords = idx_v[count]; 
	(*LATTICE) [ lattice_index (idx_v[count], y, z) ]->coords = loc_0;

	// make the switch 
	(*LATTICE)[ lattice_index (loc_0       , y, z) ]    = (*LATTICE) [ lattice_index (idx_v[count], y, z) ]; 
	(*LATTICE)[ lattice_index (idx_v[count], y, z) ]    = (*Polymers)[ index ].chain[0];
	
	*sysEnergy = energies[count];
    *contacts  = contacts_store[count];
	
	return; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRotation_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRotation
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: HeadRotation_UNBIASED and HeadRotation_BIASED
//
// PARAMETERS: index of a polymer to perform ZeroIndexRotation, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a rotation of the monomer at the first index (the tail)
// of the polymer. As it stands, this function only rotates one molecule at the tail.
//
// PLANNED EXTENSION: Multiple molecules at the tail need to be rotated.    
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully.    
//
// THE CODE: 

void HeadRotation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

    // get the neighborlist of particle at index 1 
    int dop = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2]->coords; 

    std::array <double,8> c_contacts         = *contacts;

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

	int choice = rng_uniform(0, 25); 

	// made the change to the polymer
	// make the change on the lattice 
	// make the change only to the SOLVENT site... 

	if ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] == 's' ){
		(*Polymers)[index].chain[dop-1]->coords = ne_list[choice]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 

		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[dop-1]; 
	}
	else {
		*IMP_BOOL = false; 
		return; 
	}

	double energy_n = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z); 

	if ( MetropolisAcceptance (*sysEnergy, energy_n, temperature)){
		*sysEnergy = energy_n; 
		*contacts  = c_contacts; 
	}

	else {
		
		// revert the polymer 
		(*Polymers)[index].chain[dop-1]->coords = loc_0; 
		(*LATTICE) [ lattice_index (loc_0, y, z) ]->coords  = ne_list[choice]; 
		
		// do the switch 
		// take the pointer of the solvent and put it where the monomer was on the lattice 
		
		(*LATTICE) [ lattice_index(ne_list[choice], y, z) ] = (*LATTICE)[ lattice_index (loc_0, y, z) ];
		(*LATTICE) [ lattice_index(loc_0, y, z) ]           = (*Polymers)[index].chain[dop-1];  
		*IMP_BOOL = false; 
	}

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRotation_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void HeadRotation_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts,\
	bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int index, int x, int y, int z){

    // get the neighborlist of particle at index 1 
	int dop = (*Polymers)[index].deg_poly; 
	std::array <double,8> c_contacts     = *contacts; 

    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2]->coords;
    
    double w2 {0}; // acceptance criteria weights 
    
    // first check if tail rotation can be performed at all... 
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

	// find all spots that are available 
	std::vector <std::array <int,3>> idx_v; 
	idx_v.reserve(26);

	for (std::array <int,3>& to_rot: ne_list){
		if ( to_rot == loc_0){
			continue; 
		}

		else if ( to_rot == (*Polymers)[index].chain[dop-3]->coords ) {
			continue; 
		}

		else if ( (*LATTICE)[lattice_index(to_rot, y, z) ]->ptype[0] != 'm'  ){
			idx_v.push_back( to_rot );  
		}
	}

    if ( static_cast<int>(idx_v.size()) == 0 ){
    	// if nothing can be done, just get out. 
    	*IMP_BOOL = false;
    	return; 
    }

    // otherwise, let's fuckin go
	
	// prepare the initial boltzmann sampling... 
	// get the solvent neighbors in the current location 

    // save current state... 
    // initial location is loc_0 
    // instantiate a variable that contains all energies and contacts 

    std::vector <double> energies; 
    energies.reserve(26); 
    
    std::vector <std::array<double,8>> contacts_store; 
    contacts_store.reserve(26); 
   
    for (std::array<int,3>& rot: idx_v) {
    	(*Polymers)[index].chain[dop-1]->coords = rot;
    	(*LATTICE)[ lattice_index (rot, y, z) ]->coords = loc_0; 
    	// (*(*LATTICE)[lattcei_index]).coords 

    	// do the switch 
    	(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (rot, y, z) ]; 
		(*LATTICE)[ lattice_index (rot  , y, z) ]    = (*Polymers)[index].chain[dop-1]; 

		// the LATTICE is in a new configuration! 
		// calculate energy of new configutation... 

		energies.push_back( CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z) );
		contacts_store.push_back(c_contacts);

		// now that I have the energy of the new configuration, time to revert back to the old configuration... 
		(*Polymers)[index].chain[dop-1]->coords = loc_0; 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]->coords = rot;

		// undo the switch 
		(*LATTICE)[ lattice_index (rot   , y, z) ]    = (*LATTICE)[ lattice_index (loc_0, y, z) ]; 
		(*LATTICE)[ lattice_index (loc_0 , y, z) ]    = (*Polymers)[index].chain[dop-1]; 

    }

    idx_v.push_back (loc_0);
    energies.push_back (*sysEnergy); 
    contacts_store.push_back (*contacts);

    // now that i have the energies of each config, i can calculate boltzmann weights 
    for ( double e: energies ){
    	w2 += std::exp (-1/temperature*e); 
    }

    // choose a location now that boltzmann weights have been calculated. 
    double r    = rng_uniform (0.0, 1.0); 
    double rsum = 0; 
    int count   = 0;

    for ( double e: energies ){
    	rsum += std::exp (-1/temperature*e)/w2; 
    	if ( r < rsum ){
    		break;
    	}
    	count += 1; 
    }
    std::cout << "New location = "; print (idx_v[count]);
    std::cout << "Count = " << count << std::endl;
    std::cout << "energy = " << energies[count] << std::endl;
   
    // get system to new state 
    (*Polymers)[ index ].chain[dop-1]->coords = idx_v[count]; 
	(*LATTICE) [ lattice_index (idx_v[count], y, z) ]->coords = loc_0;

	// make the switch 
	(*LATTICE)[ lattice_index (loc_0       , y, z) ]    = (*LATTICE) [ lattice_index (idx_v[count], y, z) ]; 
	(*LATTICE)[ lattice_index (idx_v[count], y, z) ]    = (*Polymers)[ index ].chain[dop-1];
	
   	*sysEnergy = energies[count];
    *contacts  = contacts_store[count];

	// CheckStructures(Polymers, LATTICE, x, y, z);
	std::cerr << "Everything seems to be okay!" << std::endl;
	// exit(EXIT_SUCCESS);

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRotation_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRotation
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: EndRotation
//
// PARAMETERS: index of polymer to perform EndRotation 
// 
// WHAT THE FUNCTION DOES: it will perform a rotation of the monomer at the terminal index (the tail or head)
// of the polymer. As it stands, this function only rotates one molecule at the tail or head. The choice of head 
// or tail rotation comes from a random distribution. 
//
// PLANNED EXTENSION: Multiple molecules at the termini need to be rotated.    
//
// DEPENDENCIES: ZeroIndexRotation, FinalIndexRotation
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void EndRotation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int index, int x, int y, int z){


	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        TailRotation_UNBIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        HeadRotation_UNBIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of EndRotation_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void EndRotation_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int index, int x, int y, int z){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = 1; // distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        TailRotation_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        HeadRotation_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of EndRotation_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//				End of EndRotation 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTIONS: 
// 1. forward_reptation_with_tail_biting
// 2. forward_reptation_without_tail_biting 
// 3. backward_repation_wtih_head_butting 
// 4. backward_reptation_without_head_butting 
//
// PARAMATERS: std::vector <Polymer>* Polymer, std::vector <Particle*>* LATTICE, std::array <int,3>* to_slither, 
// int deg_poly, int index, int y, int z, nullptr
//
// WHAT THE FUNCTIONS DO: this function performs reptation appropriately, making sure that the state is reached without 
// corrupting the system. 
//
// THE CODE: 

void forward_reptation_with_tail_biting (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <int,3>* to_slither, int deg_poly, int index, int y, int z){

	for (int i{0}; i<deg_poly; ++i){
		// push everything forward 
		// starting from the head 
		if ( i == deg_poly-1 ){
			(*Polymers)[index].chain[deg_poly-1]->coords   = *to_slither; 
			(*LATTICE)[ lattice_index (*to_slither, y, z) ] = (*Polymers)[index].chain[deg_poly-1]; 
		}
		else {
			(*Polymers)[index].chain[i]->coords = (*Polymers)[index].chain[i+1]->coords; 
			(*LATTICE)[ lattice_index ((*Polymers)[index].chain[i]->coords, y, z) ] = (*Polymers)[index].chain[i]; 
		}
	}

	return;
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of forward_reptation_with_tail_biting
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void forward_reptation_without_tail_biting (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* loc_0, int deg_poly, int index, int y, int z){

	for (int i{0}; i<deg_poly; ++i){
		if (i == deg_poly-1){
			(*Polymers)[index].chain[deg_poly-1]->coords = *to_slither; 
			(*LATTICE)[ lattice_index (*to_slither, y, z) ] = (*Polymers)[index].chain[deg_poly-1]; 
		}
		else {
			(*Polymers)[index].chain[i]->coords = (*Polymers)[index].chain[i+1]->coords; 
			(*LATTICE)[ lattice_index ((*Polymers)[index].chain[i]->coords, y, z) ] = (*Polymers)[index].chain[i]; 
		} 
	}

	// put the solvent molecule back 
	(*LATTICE)[ lattice_index (*loc_0, y, z) ] = tmp; 
	(*LATTICE)[ lattice_index (*loc_0, y, z) ]->coords = *loc_0;

	return; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of forward_reptation_without_tail_biting
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void backward_reptation_with_head_butting (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <int,3>* to_slither, int deg_poly, int index, int y, int z){

	for (int i{0}; i<deg_poly; ++i){
		// push everything forward 
		// starting from the head 
		if ( deg_poly-1-i == 0 ){
			(*Polymers)[index].chain[deg_poly-1-i]->coords   = *to_slither; // to_slither = locf 
			(*LATTICE)[ lattice_index (*to_slither, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i]; 
		}
		else {
			(*Polymers)[index].chain[deg_poly-1-i]->coords = (*Polymers)[index].chain[deg_poly-1-i-1]->coords; 
			(*LATTICE)[ lattice_index ((*Polymers)[index].chain[deg_poly-1-i]->coords, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i]; 
		}
	}

	return;
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of backward_reptation_with_head_butting
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void backward_reptation_without_head_butting (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, Particle* tmp, std::array <int,3>* to_slither, std::array <int,3>* locf, int deg_poly, int index, int y, int z){

	for (int i{0}; i<deg_poly; ++i){
		// push everything forward 
		// starting from the head 
		if ( deg_poly-1-i == 0 ){
			(*Polymers)[index].chain[deg_poly-1-i]->coords   = *to_slither; // to_slither = locf 
			(*LATTICE)[ lattice_index (*to_slither, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i]; 
		}
		else {
			(*Polymers)[index].chain[deg_poly-1-i]->coords = (*Polymers)[index].chain[deg_poly-1-i-1]->coords; 
			(*LATTICE)[ lattice_index ((*Polymers)[index].chain[deg_poly-1-i]->coords, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i]; 
		}
	}

	(*LATTICE)[lattice_index(*locf, y, z)] = tmp; 
	(*LATTICE)[lattice_index(*locf, y, z)]->coords = *locf; 

	return;
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of backward_reptation_without_head_butting
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//				End of reptation with or without collisions with another monomer 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: ForwardReptation_UNBIASED 
//
// PARAMETERS: polymer and lattice objects, along with energy, box dimensions, and topology of system
// 
// WHAT THE FUNCTION DOES: reptates a polymer forward
// 
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void ForwardReptation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	int deg_poly = (*Polymers)[index].deg_poly; 
	std::array <int,3> loc0 = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3> locf = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,8> c_contacts = *contacts; 

	// first check if tail rotation can be performed at all 
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( locf, x, y, z ); 
    int choice = rng_uniform (0, 25);	

	std::array <int,3> to_slither = ne_list [choice]; 

	Particle* tmp {nullptr}; 
	
	if ( to_slither == loc0 ){
		forward_reptation_with_tail_biting (Polymers, LATTICE, &to_slither, deg_poly, index, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		tmp = (*LATTICE)[ lattice_index (to_slither, y, z) ]; 
		forward_reptation_without_tail_biting (Polymers, LATTICE, tmp, &to_slither, &loc0, deg_poly, index, y, z);
	}

	else {
		*IMP_BOOL = false; 
		return; 
	}

	// calculate energy of current state 
	double energy_n = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z);

	if ( MetropolisAcceptance (*sysEnergy, energy_n, temperature) ){

		*sysEnergy = energy_n;
		*contacts  = c_contacts; 

	}
	else {
		// revert back to old state 
		if ( to_slither == loc0 ) {
			backward_reptation_with_head_butting (Polymers, LATTICE, &to_slither, deg_poly, index, y, z); 
		}
		else {
			tmp = (*LATTICE)[ lattice_index (loc0, y, z) ]; 
			backward_reptation_without_head_butting (Polymers, LATTICE, tmp, &loc0, &to_slither, deg_poly, index, y, z); 

		}
		*IMP_BOOL = false; 
	}

	return; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of ForwardReptation. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: BackwardReptation_UNBIASED 
//
// PARAMETERS: polymer and lattice objects, along with energy, box dimensions, and topology of system
// 
// WHAT THE FUNCTION DOES: reptates a polymer backward
// 
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void BackwardReptation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	int deg_poly = (*Polymers)[index].deg_poly; 
	std::array <int,3> loc0 = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3> locf = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,8> c_contacts = *contacts; 

	// first check if tail rotation can be performed at all 
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( loc0, x, y, z ); 
    int choice = rng_uniform (0, 25);	

	std::array <int,3> to_slither = ne_list[choice]; 

	Particle* tmp {nullptr}; 
	
	if ( to_slither == locf ){
		backward_reptation_with_head_butting (Polymers, LATTICE, &to_slither, deg_poly, index, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		tmp = (*LATTICE)[ lattice_index (to_slither, y, z) ]; 
		backward_reptation_without_head_butting (Polymers, LATTICE, tmp, &to_slither, &locf, deg_poly, index, y, z);
	}

	else {
		*IMP_BOOL = false; 
		return; 
	}

	// calculate energy of current state 
	double energy_n = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z);

	if ( MetropolisAcceptance (*sysEnergy, energy_n, temperature) ){

		*sysEnergy       = energy_n;
		*contacts        = c_contacts; 

	}
	else {
		// revert back to old state 
		if ( to_slither == locf ) {
			forward_reptation_with_tail_biting (Polymers, LATTICE, &to_slither, deg_poly, index, y, z); 
		}
		else {
			tmp = (*LATTICE)[ lattice_index (locf, y, z) ]; 
			forward_reptation_without_tail_biting (Polymers, LATTICE, tmp, &locf, &to_slither, deg_poly, index, y, z); 
		}
		*IMP_BOOL = false; 
	}

	return; 


}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of BackwardReptation. 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
// 
// NAME OF FUNCTION: Reptation_UNBIASED 
//
// PARAMETERS: index of a polymer to reptate forward or backward
// 
// WHAT THE FUNCTION DOES: reptates a polymer forwards or backwards 
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 


void Reptation_UNBIASED (std::vector<Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = 1; // distribution(generator); 

    if (num==0){
        // std::cout << "Backward reptation only!" << std::endl;
        BackwardReptation_UNBIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
        return; 
    }
    else {
        // std::cout << "Forward reptation!" << std::endl;
        ForwardReptation_UNBIASED  (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of Reptation_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
//
// 
// NAME OF FUNCTION: SolventFlip_UNBIASED, SolvationShellFlip_BIASED, PolymerFlip_BIASED, PolymerFlip_UNBIASED
//
// PARAMETERS: LATTICE, Polymers, and everything that goes into calculating energies
// 
// WHAT THE FUNCTION DOES: randomly picks a region of space, and perturbs the orientation of the solvent molecule 
// in that region. this function is a bit aggressive, imo. 
//
// DEPENDENCIES: metropolis, calculateenergy 
//
// THE CODE: 

void SolventFlip_UNBIASED ( std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array<double,8>* E, std::array<double,8>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z ){

	// number of sites to flip 
	int deg_poly = static_cast<int> ( (*Polymers)[index].chain.size() );
	int nflips = rng_uniform (1, static_cast<int>(std::ceil(deg_poly/2.0*deg_poly/2.0*deg_poly/2.0) ) ); 
	std::array <double,8> c_contacts         = *contacts; 


	std::vector <int> indices;
	std::vector <int> oris;  
	indices.reserve(nflips); 
	oris.reserve(nflips);

	int r_idx {-1}; 
	int s {0}; 

	for (int i{0}; i<nflips; ++i){
		// generate an index 
		r_idx = rng_uniform (0, x*y*z-1); 
		if ( (*LATTICE)[r_idx]->ptype[0] == 's' ){

			std::vector<int>::iterator it = std::find (indices.begin(), indices.end(), r_idx);
			// if it is a solvent, perturb orientation 
			if ( it == indices.end() ){
				indices.push_back (r_idx); 
				oris.push_back ((*LATTICE)[r_idx]->orientation); 
				// generate an orientation 
				(*LATTICE)[r_idx]->orientation = rng_uniform(0, 25); 
				
			}
			else {
				(*LATTICE)[r_idx]->orientation = rng_uniform(0, 25); 
			}
		}
		else {
			*IMP_BOOL = false; 
			break;
		}

	}
	// std::cout << "Does it reach here 2" << std::endl;
	if (*IMP_BOOL){
		// calculate the energy 
		double energy = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z); 

		// perform the metropolis acceptance 
		if ( MetropolisAcceptance (*sysEnergy, energy, temperature) ){
			*sysEnergy = energy;
			*contacts  = c_contacts; 
		}
		else {
			*IMP_BOOL = false; 
			// reverse all the perturbations performed 
			s = static_cast<int>(indices.size());
			for ( int i{0}; i < s; ++i ){
				(*LATTICE)[ indices[i] ]->orientation = oris[i];
			}
		}
	}

	else {
		// reverse all the perturbations performed 
		s = static_cast<int>(indices.size());
		for ( int i{0}; i < s; ++i ){
			(*LATTICE)[ indices[i] ]->orientation = oris[i]; 
		}
		// and that should do it. 
	}

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of SolventFlip_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of FirstSolvationShellFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void SolvationShellFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z ){

	std::set <int> solvation_shell_set; 
	std::array <std::array<int,3>, 26> ne_list; 

	// get the first solvation shell 
	for ( Polymer& pmer: *Polymers){
		for ( Particle*& p: pmer.chain ){

			ne_list = obtain_ne_list ( p->coords, x, y, z); 
			for ( std::array <int,3>& loc: ne_list ){
				if ( (*LATTICE)[ lattice_index (loc, y, z) ]->ptype[0] == 's' ){
					solvation_shell_set.insert (lattice_index (loc, y, z)); 
				}
			}
		}
	}


	std::vector <int> solvation_shell_indices (solvation_shell_set.begin(), solvation_shell_set.end()); 

	// std::cout << "Is this being hit1?" << std::endl;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// std::iota (solvation_shell_indices.begin(), solvation_shell_indices.end(), 0); 
	std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));
	// set up of solvation shell has been performed 

	// start setting up the flipping process 
    int nflip = rng_uniform(1, static_cast<int>(solvation_shell_indices.size() ) ); 			// number of sites to be flipped 

    std::vector <int> old_ori; 					// vector to hold old orientations 
    std::vector <int> new_ori; 					// vector to hold new orientations 
    old_ori.reserve(nflip); 					// reserving... 
    new_ori.reserve(nflip); 					// reserving... 

    // energy store for boltzmann sampling 
    // instantiating a bunch of variables for the process 
    int                                 ntest            = 5; 
    std::array <double,5>               energies         = {0,0,0,0,0}; 
    std::array <double,5>               boltzmann        = {0,0,0,0,0};
    std::array <int,5>                  orientations     = {0,0,0,0,0}; 
    double                              rboltzmann       = 0;  
    double                              frontflow_energy = 0; 
    double                              prob_o_to_n      = 1; 
	std::array<std::array<double,8>,5>  contacts_store   = {*contacts, *contacts, *contacts, *contacts, *contacts}; 
    std::array <double,8>               c_contacts1      = *contacts; 
    double                              Emin             = 0; 

    double rng     = 0; // rng_uniform (0.0, 1.0); 
	double rsum    = 0; 
	int    e_idx   = 0; 

	// std::cout << "solvent_indices = "; print((*solvation_shells));

    // loop over all solvent_indices
    for ( int i{0}; i < nflip; ++i ) {

    	rboltzmann = 0; 
    	old_ori.push_back( (*LATTICE)[ solvation_shell_indices[i] ]->orientation );

    	for ( int j{0}; j < ntest; ++j ){
    		
    		(*LATTICE)[ solvation_shell_indices[i] ]->orientation = rng_uniform (0, 25); 
    		orientations [j]    = (*LATTICE) [ solvation_shell_indices [i] ]->orientation; 
    		energies [j]        = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts1, x, y, z); 
    		contacts_store [j]  = c_contacts1; 

    	}

    	// std::cout << "Energies are "; print(energies); 
		Emin = *std::min_element ( energies.begin(), energies.end() ); 

		for (int k{0}; k < ntest; ++k){
			boltzmann[k] = std::exp (-1/temperature*( energies[k] - Emin ) ); 
			rboltzmann  += boltzmann[k]; 
		}


		// std::cout << "normalization = " << rboltzmann << std::endl;
		// std::cout << "rng = " << rng_acc << std::endl;

		rng     = rng_uniform (0.0, 1.0); 
		rsum    = 0; 
		e_idx   = 0; 

		for (int j{0}; j < ntest; ++j){
			rsum += boltzmann[j]/rboltzmann; 
			if ( rng < rsum ) {
				e_idx = j; 
				break; 
			}	
		}

		// make the jump to the new state 
		new_ori.push_back (orientations[e_idx]); 
		(*LATTICE)[ solvation_shell_indices[i] ]->orientation = orientations[e_idx]; 
		prob_o_to_n *= boltzmann[e_idx]/rboltzmann; 
		 
    }
    
    frontflow_energy  = energies       [e_idx];
	c_contacts1       = contacts_store [e_idx];


    // figure out the backflow energy 
    double prob_n_to_o     = 1; 
    double backflow_energy = 1; 
    std::array <double,8> c_contacts2 = {0,0,0,0,0,0,0,0};

    for ( int i{0}; i < nflip; ++i ){

    	rboltzmann = 0; 
    	(*LATTICE) [ solvation_shell_indices[i] ]->orientation = old_ori[i]; 
    	energies[0] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z);
    	contacts_store[0] = c_contacts2; 

    	for ( int j{1}; j < ntest; ++j ){
    		(*LATTICE)[ solvation_shell_indices[i] ]->orientation = rng_uniform(0, 25); 
    		energies[j] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z);
    		contacts_store[j] = c_contacts2; 
    	}

		// std::cout << "Energies are "; print(energies); 

		Emin = *std::min_element ( energies.begin(), energies.end() ); 
		// std::cout << "Emin = " << Emin << std::endl; 

		for (int k{0}; k < ntest; ++k){
			boltzmann[k] = std::exp(-1/temperature*( energies[k] - Emin ) ); 
			rboltzmann  += boltzmann[k]; 
		}
		prob_n_to_o     *= boltzmann[0]/rboltzmann; 
		
		// make the jump to the old state 
		(*LATTICE) [ solvation_shell_indices [i] ]->orientation = old_ori[i]; 

    }

    backflow_energy  = energies[0];
    c_contacts2      = contacts_store[0]; 

    if ( backflow_energy != *sysEnergy || c_contacts2 != *contacts ){
    	std::cout << "Something is fucked. Energies do not match." << std::endl;
    	std::cout << "backflow_energy = " << backflow_energy << ", sysEnergy = " << *sysEnergy << std::endl;
    	std::cout << "c_contacts = "; print (c_contacts2, ", "); std::cout << "contacts = "; print (*contacts); 
    	exit(EXIT_FAILURE); 
    }

    // check the acceptance criterion 

	double rng_acc = rng_uniform (0.0, 1.0); 
	if ( rng_acc < std::exp (-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n  ) {
		// if accepted, return to the new orientations 
		
		for ( int j{0}; j < nflip; ++j ) {
			(*LATTICE)[ solvation_shell_indices[j] ]->orientation = new_ori[j]; 
		}

		double en = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z); 
		
		if ( en != frontflow_energy || c_contacts2 != c_contacts1 ){
    		std::cout << "Something is fucked. Energies do not match." << std::endl;
    		std::cout << "en = " << en << ", frontflow energy = " << frontflow_energy << std::endl;
    		std::cout << "c_contacts2 = "; print (c_contacts2, ", "); std::cout << "c_contacts1 = "; print (c_contacts1); 
    		exit(EXIT_FAILURE); 
    	}

		*sysEnergy = frontflow_energy; 
		*contacts  = c_contacts1; 

	}
	else {
		*IMP_BOOL = false; 
	}
	
	// CheckStructures (Polymers, Cosolvent, LATTICE, x, y, z);
	return;
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of FirstSolvationShellFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of SecondSolvationShellFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of PolymerFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void PolymerFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array<double,8>* E, std::array<double,8>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z ){

	int deg_poly = (*Polymers)[0].deg_poly; 
	std::vector <int> polymer_indices (deg_poly);
	std::iota (polymer_indices.begin(), polymer_indices.end(), 0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (polymer_indices.begin(), polymer_indices.end(), std::default_random_engine(seed)); 

	// std::cout << "polymer indices = "; print (polymer_indices); 

	int ntest    = 5; 
	int nflip    = (deg_poly==1) ? 1 : rng_uniform (1, deg_poly-1); // if deg poly is 1, nflip = 1, otherwise do rng(1, degpoly-1)
	// std::cout << "nflip = " << nflip << std::endl;

	std::vector <int> old_ori; 
	std::vector <int> new_ori;
	old_ori.reserve(nflip);
	new_ori.reserve(nflip);

	std::array<double,5>               energies           = {0,0,0,0,0}; 
	std::array<double,5>               boltzmann          = {0,0,0,0,0};
	std::array<int,5>                  orientations       = {0,0,0,0,0};
	std::array<std::array<double,8>,5> contacts_store     = {*contacts, *contacts, *contacts, *contacts, *contacts}; 
	double                             rboltzmann         = 0;
	double                             frontflow_energy   = 0; 
	double                             prob_o_to_n        = 1; 
	std::array  <double,8>             c_contacts1        = *contacts; 
	double                             Emin               = 0; 

	double rng   = 0; 
	double rsum  = 0;
	int    e_idx = 0; 

	// loop over all solvent_indices 
	for ( int i{0}; i < nflip; ++i ){

		rboltzmann = 0; 
		old_ori.push_back ( (*Polymers)[index].chain[ polymer_indices[i] ]->orientation ); 

		for ( int j{0}; j < ntest; ++j ){

			(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = rng_uniform (0, 25); 
			orientations[j]   = (*Polymers)[index].chain[ polymer_indices[i] ]->orientation;
			energies[j]       = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts1, x, y, z); 
			contacts_store[j] = c_contacts1;   

		}

		Emin = *std::min_element ( energies.begin(), energies.end() ); 

		for ( int k{0}; k < ntest; ++k ){
			boltzmann[k] = std::exp (-1/temperature * (energies[k] - Emin) );
			rboltzmann  += boltzmann [k]; 
		}

		// std::cout << "normalization = " << rboltzmann << std::endl;
		// std::cout << "rng = " << rng_acc << std::endl;

		rng   = rng_uniform (0.0, 1.0); 
		rsum  = 0; 
		e_idx = 0;

		for ( int j{0}; j<5; ++j){
			rsum += boltzmann[j]/rboltzmann;
			if ( rng < rsum ){
				e_idx = j; 
				break; 
			}
		}

		// make the jump to the new state 
		new_ori.push_back (orientations[e_idx]); 
		(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = orientations[e_idx]; 
		prob_o_to_n *= boltzmann[e_idx]/rboltzmann; 

	}

	frontflow_energy = energies[e_idx]; 
	c_contacts1      = contacts_store[e_idx]; 

	// std::cout << "Starting backflow... " << std::endl;

	// figure out the backflow energy 
	double prob_n_to_o     = 1; 
	double backflow_energy = 1;
	std::array <double,8> c_contacts2 = {0,0,0,0,0,0,0,0}; 

	for ( int i{0}; i < nflip; ++i ){

		// std::cout << "i = " << i << std::endl;
		rboltzmann = 0; 
		(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = old_ori[i]; 
		energies[0] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z); 
		contacts_store[0] = c_contacts2; 

		for (int j{1}; j<ntest; ++j){
			(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = rng_uniform (0, 25); 
			energies[j] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z); 
			contacts_store[j] = c_contacts2; 
		}

		Emin = *std::min_element ( energies.begin(), energies.end() ); 

		for (int k{0}; k < ntest; ++k){
			boltzmann [k] = std::exp (-1/temperature * (energies[k] - Emin ) );
			rboltzmann   += boltzmann[k]; 
		}
		prob_n_to_o      *= boltzmann[0]/rboltzmann; 

		// make the jump to the old state 
		(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = old_ori[i]; 

	}

	// std::cout << "Done with backflow! " << std::endl;

	backflow_energy = energies[0]; 
	c_contacts2     = contacts_store[0]; 

	if ( backflow_energy != *sysEnergy || c_contacts2 != *contacts ){
    	std::cout << "Something is fucked. Energies do not match." << std::endl;
    	std::cout << "backflow_energy = " << backflow_energy << ", sysEnergy = " << *sysEnergy << std::endl;
    	std::cout << "c_contacts = "; print (c_contacts2, ", "); std::cout << "contacts = "; print (*contacts); 
    	exit(EXIT_FAILURE); 
    }

    // if ( c_solvation_shells != *solvation_shells ){
    	// std::cerr << "Something is fucked. Solvation shell should not be affected by polymer flips." << std::endl;
    	// exit (EXIT_FAILURE); 
    // }

    // check the acceptance criterion 

    double rng_acc = rng_uniform (0.0, 1.0); 

    if ( rng_acc < std::exp (-1/temperature * (frontflow_energy - *sysEnergy ) ) ){
    	// if accepted, return to the new orientations 
    	for (int j{0}; j < nflip; ++j){
    		(*Polymers)[index].chain[polymer_indices[j]]->orientation = new_ori[j];
    	}

		double en = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z); 
		
		if ( en != frontflow_energy || c_contacts2 != c_contacts1 ){
    		std::cout << "Something is fucked. Energies do not match." << std::endl;
    		std::cout << "en = " << en << ", frontflow energy = " << frontflow_energy << std::endl;
    		std::cout << "c_contacts2 = "; print (c_contacts2, ", "); std::cout << "c_contacts1 = "; print (c_contacts1); 
    		exit(EXIT_FAILURE); 
    	}

		*sysEnergy = frontflow_energy; 
		*contacts  = c_contacts1; 

    }
    else {
    	*IMP_BOOL = false; 
    }

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of PolymerFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void PolymerFlip_UNBIASED ( std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array<double,8>* E, std::array<double,8>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, int p_index, int x, int y, int z ){

	// number of monomers to flip 
	int deg_poly = static_cast<int> ( (*Polymers)[p_index].chain.size() );
	int nflips = rng_uniform (1, static_cast<int>(deg_poly/2)); 
	std::array <double,8> c_contacts = *contacts; 

	std::vector <int> indices;
	std::vector <int> oris;  
	indices.reserve(nflips); 
	oris.reserve(nflips);	

	int r_idx {-1}; 
	std::vector<int>::iterator it;

	for (int i{0}; i<nflips; ++i){

		r_idx = rng_uniform (0, deg_poly-1); 
		it = std::find (indices.begin(), indices.end(), r_idx);

		if ( it == indices.end() ){
			indices.push_back (r_idx);
			oris.push_back ((*Polymers)[p_index].chain[r_idx]->orientation); 
			(*Polymers)[p_index].chain[r_idx]->orientation = rng_uniform(0, 25); 
		}
		else {
			(*Polymers)[p_index].chain[r_idx]->orientation = rng_uniform(0,25); 
		}

	}

	double energy = CalculateEnergy_parallel( Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z); 
	
	if ( MetropolisAcceptance (*sysEnergy, energy, temperature) ){
		
		*sysEnergy = energy;
		*contacts  = c_contacts; 

	}
	else {
		*IMP_BOOL = false; 
		// reverse all the perturbations performed 
		int s = static_cast<int>(indices.size());
		for ( int i{0}; i < s; ++i ){
			(*Polymers)[p_index].chain[ indices[i] ]->orientation = oris[i]; 
		}
	}

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of PolymerFlip_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			End of unbiased orientation flipping moves 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: ChainRegrowth_BIASED, HeadRegrowth_BIASED, TailRegrowth_BIASED
//
// PARAMETERS: Polymers, LATTICE, and all the stuff that goes into calculating energy 
// 
// WHAT THE FUNCTION DOES: Performs a chain regrowth on a polymer. 
//
// DEPENDENCIES: metropolis, calculateenergy 
//
// THE CODE: 

void ChainRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int p_index, int x, int y, int z){


	int deg_poly = (*Polymers)[p_index].deg_poly; 
	int m_index  = rng_uniform (1, deg_poly-2); 


	std::array <double,8> c1_contacts     = *contacts; 
	std::array <double,8> c2_contacts     = *contacts; 

	double prob_o_to_n {1}; 
	double prob_n_to_o {1}; 
	double frontflow_energy {*sysEnergy}; 
	double backflow_energy  {0}; 
	int    recursion_depth  {0}; 

	// old_cut contains the old positions of thr monomer 
	std::vector <std::array<int,3>> old_cut; 
	std::vector <std::array<int,3>> new_cut;
	old_cut.reserve (deg_poly);
	new_cut.reserve (deg_poly); 
	 

	// regrowth will be in the direction where there are fewer monomer
	// if m_index/deg_poly is lesser than 0.5, growth is false, otherwise growth is true 
	int growth {-1}; 

	if ( deg_poly % 2 == 0 ){
		growth = (0.5 >= (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
	}
	else {
		if ( 0.5 == (m_index+1)/static_cast<double>(deg_poly+1) ){
			growth = rng_uniform (0, 1);
		}
		else {
			growth = (0.5 > (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
		}
	}
	
	if ( growth ){
		
		// std::cout << "Head regrowth... " << std::endl;
		// std::cout << "OLD ENERGY = " << *sysEnergy << std::endl << std::endl;
		// get old_cut 
		
		for (int i {m_index+1}; i<deg_poly; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
		}

		// regrow the polymer frontwards
		HeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		

		for (int i {m_index+1}; i<deg_poly; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
		}

		if ( old_cut == new_cut ){
			// std::cout << "------------------------------------------------------------------------"  << std::endl;
			// std::cout << "POLYMER CONFIGURATION WAS NOT CHANGED. RETURN BACK TO MAIN CONFIGURATION." << std::endl;
			// std::cout << "------------------------------------------------------------------------"  << std::endl << std::endl;
			return; 
		}

		if ( !(*IMP_BOOL) ){
			// std::cout << "ALL BLOCKS! REVERTING!" << std::endl;
			// revert back to the original state. 

			acceptance_after_head_regrowth ( LATTICE, &new_cut, &old_cut, y, z );

			return; 
		}

		backflow_energy = frontflow_energy; 
		BackFlowFromHeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

		if ( *sysEnergy != backflow_energy || c2_contacts != *contacts ){
			std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
			std::cout << "*sysEnergy = " << *sysEnergy << ", backflow_energy = " << backflow_energy << "." << std::endl;
			std::cout << "c2_contacts = "; print (c2_contacts, ", "); std::cout << "*contacts = "; print(*contacts);
			std::cout << "Shit's fucked." << std::endl;
			exit(EXIT_FAILURE);
		}

		// check acceptance criterion
		double rng_acc = rng_uniform (0.0, 1.0); 
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy- *sysEnergy)) * prob_n_to_o/prob_o_to_n ){

			acceptance_after_head_regrowth (LATTICE, &old_cut, &new_cut, y, z); 
			*sysEnergy = frontflow_energy; 
			*contacts  = c1_contacts;

		}

		else {
			*IMP_BOOL = false; 
		}

	}
	else {

		// std::cout << "Tail regrowth..." << std::endl;
		// get old cut 
		for ( int i{0}; i<m_index; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords);
		}

		// regrow the polymer backwards
		TailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 

		for (int i {0}; i<m_index; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
		}

		if ( old_cut == new_cut ){
			// std::cout << "------------------------------------------------------------------------"  << std::endl;
			// std::cout << "POLYMER CONFIGURATION WAS NOT CHANGED. RETURN BACK TO MAIN."               << std::endl;
			// std::cout << "------------------------------------------------------------------------"  << std::endl << std::endl;
			return; 
		}

		
		if ( !(*IMP_BOOL) ){

			acceptance_after_tail_regrowth ( LATTICE, &new_cut, &old_cut, y, z); 
			
			return; 
		}

		backflow_energy = frontflow_energy; 

		// std::cout << "BEGIN BACK FLOW! " << std::endl;

		BackFlowFromTailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

		if ( *sysEnergy != backflow_energy || c2_contacts != *contacts ){
			std::cout << "Energies are bad, or contacts are not right." << std::endl;
			std::cout << "*sysEnergy = " << *sysEnergy << ", backflow_energy = " << backflow_energy << "." << std::endl;
			std::cout << "c2_contacts = "; print (c2_contacts, ", "); std::cout << "*contacts = "; print(*contacts);
			std::cout << "Shit's fucked." << std::endl;
			exit(EXIT_FAILURE);
		}

		// check acceptance criterion 
		double rng_acc = rng_uniform (0.0, 1.0); 
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy- *sysEnergy)) * prob_n_to_o/prob_o_to_n ){
			// accept new cut 
			// perform swaps 
			acceptance_after_tail_regrowth (LATTICE, &old_cut, &new_cut, y, z); 

			*sysEnergy = frontflow_energy; 
			*contacts  = c1_contacts;
		}
		else {
			*IMP_BOOL = false; 
		}

	}

	return; 
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of ChainRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of HeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void HeadRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	std::array <double,8> current_contacts = *contacts; 

	std::array <double,5> energies; 
	std::array <std::array<double,8>,5> contacts_store; 
	std::array <double,5> boltzmann; 

	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 
	
	// start attempting jumps 

	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 


	while ( idx_counter < 5 ){

		// std::cout << "ptype = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;

		if ( ne_list[idx_counter] == loc_m ){ 
			// current position
			energies[idx_counter]       = *frontflow_energy; 
			contacts_store[idx_counter] = current_contacts;
			block_counter += 1; 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			// check which region of the polymer the swap is taking place 

			for ( int u{0}; u<deg_poly; ++u ) {
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ) {
					self_swap_idx = u; 
					break; 
				}
			}

			// if it is with other head units, do the swap 
			// if not, discourage it 
			// std::cout << "self_swap_idx = " << self_swap_idx << std::endl; 

			if ( self_swap_idx < m_index ){

				// maintain current state, and sample another state 
				// maintain_idx.push_back(idx_counter); 

				energies[idx_counter] = 1e+08; // very unfavorable state 
				contacts_store[idx_counter] = {-1,-1,-1,-1}; 
				block_counter += 1; 
			}
			else {
				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the energy
				energies [idx_counter] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E,contacts, x, y, z); 
				contacts_store[idx_counter] = *contacts; 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];
			}
			
		}
		else {
			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

			// get the energy
			energies [idx_counter] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
			contacts_store[idx_counter] = *contacts; 

			// revert back to original structure 
			(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
			(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
			(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];
		}

		idx_counter += 1; 
		
	}

	if ( block_counter == 5 ){
		*IMP_BOOL = false; 
		return; 
	}

	// now that i have all the energies, and boltzmann weights, i can choose a configuration 
	// std::cout << "Energies are "; print(energies); 

	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<5; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	double rng_acc = rng_uniform (0.0, 1.0); 
	double rsum    = 0; 
	int    e_idx   = 0; 
	
	// std::cout << "normalization = " << rboltzmann << std::endl;
	// std::cout << "rng = " << rng_acc << std::endl;

	for (int j{0}; j<5; ++j){
		rsum += boltzmann[j]/rboltzmann; 
		if ( rng_acc < rsum ){
			e_idx = j; 
			break; 
		}	
	}
	
	// now that I have chosen a configuration, go with it
	*prob_o_to_n     *= boltzmann[e_idx]/rboltzmann; 
	*frontflow_energy = energies[e_idx];
	*contacts         = contacts_store [e_idx]; 
	// std::cout << "e_idx = " << e_idx << std::endl;
	// std::cout << "maintain_idx = "; print(maintain_idx);
	// std::cout << "frontflow_energy is " << *frontflow_energy << std::endl; 
	// std::cout << "position chosen = "; print (ne_list[ e_idx ]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(o->n) = " << *prob_o_to_n << std::endl;
	// std::cout << "------------------------------" << std::endl;
	// std::cout << std::endl;


	// do the swap again
	(*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords           = ne_list[e_idx];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];

	// else { do nothing and maintain structure, because suggested index is in the maintain index vector }
	HeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of BackFlowFromHeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	
	std::array <double,5> energies; 
	std::array <double,5> boltzmann; 
	std::array <double,8> current_contacts = *contacts; 
	std::array <std::array<double,8>,5> contacts_store; 


	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; // this is the key item of interest

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[recursion_depth]) - ne_list.begin() ; 

	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	// std::cout << "shuffled ne_list is "; 
	// for (std::array <int,3>& n: ne_list){
		// print(n, ", ");
	// }
	// std::cout << std::endl;

	// i now have a vector which has the back peddling step at position index 0 

	// start attempting jumps 
	int idx_counter         = 0; 
	int self_swap_idx 		= -1; 

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);
		if ( ne_list[idx_counter] == loc_m ){
			energies[idx_counter]       = *backflow_energy;
			contacts_store[idx_counter] = current_contacts;
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list [idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}
			}

			if ( self_swap_idx < m_index ){
				// maintain current state 
				// std::cout << "self_swap_idx = " << self_swap_idx << std::endl;
				// std::cout << "maintain current state..." << std::endl;
				// std::cout << "Run into an undoable swap during backflow. " << std::endl; 
				// std::cout << "Selected position is "; print (ne_list[idx_counter]);
				energies [idx_counter] = 1e+08; 
				contacts_store[idx_counter] = {-1,-1,-1,-1}; 
			}
			else {
				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap 
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
				
				// get the energy
				energies [idx_counter] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = *contacts; 
				// std::cout << "current_contacts = "; print(current_contacts); 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];
			}
		}
		else {
			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index+1];

			// get the energy 
			energies [idx_counter] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
			contacts_store [idx_counter] = *contacts; 

			// revert back to original structure 
			(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
			(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
			(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];
		}

		idx_counter += 1;

	}

	// std::cout << "Energies are "; print(energies); 
	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<5; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	// std::cout << "normalization = " << rboltzmann << std::endl;

	*prob_n_to_o      = (*prob_n_to_o)*boltzmann[0]/rboltzmann;
	*backflow_energy  = energies[0]; 
	*contacts         = contacts_store [0];

	// std::cout << "position chosen = "; print (ne_list[0]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(n->o) = " << *prob_n_to_o << std::endl;
	// std::cout << "------------------------------" << std::endl;


	// do the swap again
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords           = ne_list[0];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];

	// else { do nothing, MAINTAIN, because suggested index is in the maintain index vector }
	BackFlowFromHeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, old_cut, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of BackflowHeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of TailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void TailRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	std::array <double,8> current_contacts = *contacts; 
	std::array <double,5> energies; 
	std::array <std::array<double,8>,5> contacts_store; 
	std::array <double,5> boltzmann; 

	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// std::cout << "shuffled ne_list is "; 
	// for (std::array <int,3>& n: ne_list){
		// print(n, " : ");
	// }
	// std::cout << std::endl;
	

	// start attempting jumps 
	int block_counter       = 0; 
	int idx_counter         = 0; 
	int self_swap_idx       = -1; 

	while ( idx_counter < 5 ){

		// std::cout << "ptype = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;

		if ( ne_list[idx_counter] == loc_m ){
			energies[idx_counter]                = *frontflow_energy;
			contacts_store[idx_counter]          = current_contacts;
			block_counter += 1; 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			// check which region of the polymer the swap is taking place 

			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ){
					self_swap_idx = u; 
					break; 
				}				
			}

			// if it is with other tail units, do the swap 
			// if not, discourage it 

			// std::cout << "self_swap_idx = " << self_swap_idx << std::endl; 
			if ( self_swap_idx > m_index ){
				// std::cout << "maintain current state..." << std::endl;
				// std::cout << "You have run into an undoable swap. " << std::endl; 
				// std::cout << "Selected position is "; print (ne_list[idx_counter]);	
				energies[idx_counter] = 1e+08; 
				contacts_store[idx_counter] = {-1,-1,-1,-1}; 
				block_counter += 1; 
			}
			else {
				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the energy
				energies [idx_counter] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = *contacts; 

				// revert back to original structure 
				(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords   = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
			}

		}
		else {
			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the energy
			energies [idx_counter] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
			contacts_store[idx_counter] = *contacts; 

			// revert back to original structure 
			(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
			(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
			(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
		}

		idx_counter += 1;

	}

	if ( block_counter == 5 ){
		*IMP_BOOL = false;
		return; 
	}

	// now that i have all the energies, and boltzmann weights, i can choose a configuration 
	// std::cout << "Energies are "; print(energies); 

	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<5; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	double rng_acc = rng_uniform (0.0, 1.0); 
	double rsum    = 0; 
	int e_idx      = 0; 
	// std::cout << "normalization = " << rboltzmann << std::endl;
	// std::cout << "rng = " << rng_acc << std::endl;

	for (int j{0}; j<5; ++j){
		rsum += boltzmann[j]/rboltzmann; 
		if ( rng_acc < rsum ) {
			e_idx = j; 
			break; 
		}	
	}

	// now that I have chosen a configuration, go with it
	*prob_o_to_n      = (*prob_o_to_n) * boltzmann[e_idx]/rboltzmann; 
	*frontflow_energy = energies[e_idx];
	*contacts         = contacts_store [e_idx];

	// std::cout << "e_idx = " << e_idx << std::endl;
	// std::cout << "frontflow_energy is " << *frontflow_energy << std::endl; 
	// std::cout << "position chosen = "; print (ne_list[e_idx]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(o->n) = " << *prob_o_to_n << std::endl;
	// std::cout << "------------------------------" << std::endl;
	// std::cout << std::endl; 

	// if { e_idx is not in the maintain index vector, perform the swap}
		
	// do the swap again 	
	(*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords           = ne_list[e_idx];
	
	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
	

	// else { do nothing and maintain structure, because the suggested index is in the the maintain index vector }
	TailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of BackFlowFromTailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromTailRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::array<double,8>* E, std::array <double,8>* contacts, \
	bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, \
	int p_index, int m_index, int recursion_depth, int x, int y, int z){

	if (m_index == 0){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	
	// generate an array for energies 
	
	std::array <double,5> energies; 
	std::array <double,5> boltzmann; 

	std::array <double,8> current_contacts = *contacts;
	std::array <std::array<double,8>,5> contacts_store; 

	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 
	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[m_index-1]) - ne_list.begin() ; 

	// make sure old_cut is present at position 0 
	// crit_idx = 1; 
	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	// std::cout << "shuffled ne_list is "; 
	// for (std::array <int,3>& n: ne_list){
		// print(n, ", ");
	// }
	// std::cout << std::endl;

	// i now have a vector which has the back peddling position at position index 0 

	// start attempting jumps 
	int idx_counter         = 0 ; 
	int self_swap_idx		= -1;

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);
		if ( ne_list [idx_counter] == loc_m ){
			energies[idx_counter]               = *backflow_energy;
			contacts_store[idx_counter]         = current_contacts;
			// std::cout << "current_contacts = "; print(current_contacts); 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){
			
			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}				
			}

			if ( self_swap_idx > m_index ){
				// maintain current state
				// std::cout << "self_swap_idx = " << self_swap_idx << std::endl;
				// std::cout << "maintain current state..." << std::endl;
				// std::cout << "Run into an undoable swap during backflow. " << std::endl; 
				// std::cout << "Selected position is "; print (ne_list[idx_counter]);
				energies[idx_counter] = 1e+08; 
				contacts_store[idx_counter] = {-1,-1,-1,-1}; 
			}

			else {
				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the energy
				energies [idx_counter]               = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter]          = *contacts; 
				// std::cout << "current_contacts = "; print(current_contacts); 

				// revert back to original structure 
				(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
			}

		}
		else {
			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords                 = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ]                = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the energy
			energies [idx_counter]               = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
			contacts_store[idx_counter]          = *contacts; 
			// std::cout << "current_contacts = "; print(current_contacts); 

			// revert back to original structure 
			(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
			(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
			(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index-1];
		}

		idx_counter += 1; 
		
	}

	// std::cout << "Energies are "; print(energies); 
	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<5; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	// std::cout << "normalization = " << rboltzmann << std::endl;
	
	*prob_n_to_o     *= boltzmann[0]/rboltzmann;
	*backflow_energy  = energies[0]; 
	*contacts         = contacts_store[0];

	// std::cout << "position chosen = "; print (ne_list[0]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(n->o) = " << *prob_n_to_o << std::endl;
	// std::cout << "------------------------------" << std::endl;


	// do the swap again to the positions dictated by old_cut  
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords           = ne_list[0];
	
	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

	BackFlowFromTailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, old_cut, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

	return; 

}



// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of all biased chain regrowth moves 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

//==============================================================================================
//==============================================================================================
//
// 
// NAME OF FUNCTION: ChainRegrowth_UNBIASED, HeadRegrowth_UNBIASED, TailRegrowth_UNBIASED
//
// PARAMETERS: Polymers, LATTICE, and all the stuff that goes into calculating energy 
// 
// WHAT THE FUNCTION DOES: Performs a chain regrowth on a polymer. 
//
// DEPENDENCIES: metropolis, calculateenergy 
//
// THE CODE: 

void ChainRegrowth_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int p_index, int x, int y, int z){

	// choose an index to regrow from 
	int deg_poly = (*Polymers)[p_index].deg_poly; 
	int m_index  = rng_uniform (1, deg_poly-2); 
	// std::cout << "m_index = " << m_index << std::endl;

	std::array <double,8> c1_contacts = *contacts; 

	double frontflow_energy {*sysEnergy}; 
	
	// old_cut contains the old positions of thr monomer 
	std::vector <std::array<int,3>> old_cut; 
	std::vector <std::array<int,3>> new_cut; 
	old_cut.reserve(deg_poly);
	new_cut.reserve(deg_poly);
	 

	// regrowth will be in the direction where there are fewer monomer
	// if m_index/deg_poly is lesser than 0.5, growth is false, otherwise growth is true 
	int growth {-1}; 

	if ( deg_poly % 2 == 0 ){
		growth = (0.5 >= (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
	}
	else {
		if ( 0.5 == (m_index+1)/static_cast<double>(deg_poly+1) ){
			growth = rng_uniform (0, 1);
		}
		else {
			growth = (0.5 > (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
		}
	}
	
	if ( growth ){

		for (int i {m_index+1}; i<deg_poly; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
		}

		// regrow the polymer frontwards
		HeadRegrowth_UNBIASED (Polymers, Cosolvent, LATTICE, IMP_BOOL, deg_poly, p_index, m_index, x, y, z); 

		for (int i {m_index+1}; i<deg_poly; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
		}

		// if the move could not be completed, i need to back track 
		if ( !(*IMP_BOOL) ){

			acceptance_after_head_regrowth (LATTICE, &new_cut, &old_cut, y, z); 
			return; 
		}

		frontflow_energy = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c1_contacts, x, y, z); 

		// check acceptance criterion
		double rng_acc = rng_uniform (0.0, 1.0); 
		
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy - *sysEnergy)) ){
			
			*sysEnergy = frontflow_energy; 
			*contacts  = c1_contacts; 
		}
		else {
			// since move was not accepted, need to get the original configuration 
			*IMP_BOOL = false; 
			acceptance_after_head_regrowth (LATTICE, &new_cut, &old_cut, y, z); 
		}
	}
	else {

		// std::cout << "Tail regrowth... " << std::endl;
		for ( int i{0}; i<m_index; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords);
		}

		// regrow the polymer backwards
		TailRegrowth_UNBIASED (Polymers, Cosolvent, LATTICE, IMP_BOOL, deg_poly, p_index, m_index, x, y, z); 

		for ( int i{0}; i<m_index; ++i) {
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords); 
		}
		
		// if the move could not be completed, i need to back track 
		if ( !(*IMP_BOOL) ){

			acceptance_after_tail_regrowth (LATTICE, &new_cut, &old_cut, y, z);
			return; 
		}

		frontflow_energy = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c1_contacts, x, y, z); 
		
		double rng_acc = rng_uniform (0.0, 1.0); 

		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy - *sysEnergy)) ){
			
			*sysEnergy = frontflow_energy; 
			*contacts  = c1_contacts; 
		}
		else {
			// since move was not accepted, need to get the original configuration 
			*IMP_BOOL = false; 
			acceptance_after_tail_regrowth (LATTICE, &new_cut, &old_cut, y, z); 
		}

	}

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of ChainRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of HeadRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void HeadRegrowth_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	bool* IMP_BOOL, int deg_poly, int p_index, \
	int m_index, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	int rint = rng_uniform (0, 25); 

	// std::cout << "Suggested location = "; print(ne_list[rint]);  

	if ( (*LATTICE)[ lattice_index (ne_list[rint], y, z) ]->ptype[0]=='s' ){

		(*LATTICE)[lattice_index(ne_list[rint], y, z)]->coords = loc_m;
		(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[rint];
		
		// perform the swap (since coords were changed, this swap works)
		(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[rint], y, z)];
		(*LATTICE)[ lattice_index (ne_list[rint], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
	}

	else if ( (*LATTICE)[ lattice_index (ne_list[rint], y, z) ]->ptype[0]=='m' ){

		int self_swap_idx = -1; 
		// check which region of the polymer the swap is taking place 
		for ( int u{0}; u<deg_poly; ++u ){
			if ( (*Polymers)[p_index].chain[u]->coords == ne_list[rint] ){
				self_swap_idx = u; 
				break; 
			}
		}

		// if it is with other head units do the swap
		// if not, stop the growth 

		if ( self_swap_idx < m_index+1 ){
			*IMP_BOOL = false;
			return;
		}
		else {
			// perform the swap 
			(*LATTICE)[lattice_index(ne_list[rint], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[rint];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[rint], y, z)];
			(*LATTICE)[ lattice_index (ne_list[rint], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
		}

	}

	else {
		*IMP_BOOL = false; 
		return; 
	}
	
	HeadRegrowth_UNBIASED (Polymers, Cosolvent, LATTICE, IMP_BOOL, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of TailRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void TailRegrowth_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	bool* IMP_BOOL, int deg_poly, int p_index, \
	int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true;
		return; 
	}

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	int rint = rng_uniform (0, 25); 

	// std::cout << "Suggested location = "; print(ne_list[rint]); 

	if ( (*LATTICE)[ lattice_index (ne_list[rint], y, z) ]->ptype[0]=='s' ){

		(*LATTICE)[lattice_index(ne_list[rint], y, z)]->coords = loc_m;
		(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[rint];
		
		// perform the swap (since coords were changed, this swap works)
		(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[rint], y, z)];
		(*LATTICE)[ lattice_index (ne_list[rint], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
	}

	else if ( (*LATTICE)[ lattice_index (ne_list[rint], y, z) ]->ptype[0]=='m' ){

		int self_swap_idx = -1; 
		// check which region of the polymer with the swap is taking place 
		for ( int u{0}; u<deg_poly; ++u ){
			if ( (*Polymers)[p_index].chain[u]->coords == ne_list[rint] ){
				self_swap_idx = u; 
				break; 
			}
		}

		// if it is with other tail units do the swap 
		// if not, stop the growth 

		if ( self_swap_idx > m_index-1){
			*IMP_BOOL = false; 
			return;
		}
		else {
			// perform the swap 
			(*LATTICE)[lattice_index(ne_list[rint], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[rint];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[rint], y, z)];
			(*LATTICE)[ lattice_index (ne_list[rint], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

		}

	}

	else {
		*IMP_BOOL = false; 
		return; 
	}
	
	TailRegrowth_UNBIASED (Polymers, Cosolvent, LATTICE, IMP_BOOL, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
//
// 
// NAME OF FUNCTION: acceptance_after_head_regrowth, acceptance_after_tail_regrowth
//
// PARAMETERS: Polymers, LATTICE, and all the stuff that goes into calculating energy 
// 
// WHAT THE FUNCTION DOES: Performs a chain regrowth on a polymer. 
//
// DEPENDENCIES: metropolis, calculateenergy 
//
// THE CODE: 

void acceptance_after_head_regrowth (std::vector <Particle*>* LATTICE, \
	std::vector <std::array<int,3>>* old_cut, std::vector <std::array<int,3>>* new_cut, int y, int z){


	std::vector< std::vector<std::array<int,3>>> master_linked_list;
	std::vector <std::array<int,3>> link; 

	create_linked_list (*old_cut, *new_cut, link, &master_linked_list, 1);

	int L; 
	Particle* tmp_par_ptr  {nullptr}; 
	Particle* tmp_par_ptr_ {nullptr}; 

	for ( std::vector<std::array<int,3>>& linked_list: master_linked_list ){
		
		L = static_cast<int> (linked_list.size()) ; 
		
		// check if the final object location in link goes to a solvent
		// then do a cyclical transition, starting from the end index and make it to the first index 
		if ( (*LATTICE)[ lattice_index(linked_list[L-1], y, z) ]->ptype[0] == 's' ){

			tmp_par_ptr = (*LATTICE) [ lattice_index (linked_list.back(), y, z) ]; 
			// go backwards 
			for (int i{0}; i < L ; i=i+2 ){

				(*LATTICE)[lattice_index (linked_list[L-2-i], y, z)]->coords = linked_list[L-1-i]; 
				(*LATTICE)[lattice_index (linked_list[L-1-i], y, z)] = (*LATTICE)[lattice_index( linked_list[L-2-i], y, z )]; 
				
			}
			tmp_par_ptr->coords = linked_list[0]; 
			(*LATTICE)[lattice_index (linked_list[0], y, z)] = tmp_par_ptr; 

		}

		else { // when it is a circulation of monomers 

			for ( int i{0}; i < L; i=i+2){
			
				if ( i == 0 ) {
					// store info about linked_list[1]... and consequently also linked_list[2]
					tmp_par_ptr = (*LATTICE)[lattice_index ( linked_list[i+1], y, z) ];

					(*LATTICE)[ lattice_index ( linked_list[i],   y, z) ]->coords = linked_list[i+1]; 
					(*LATTICE)[ lattice_index ( linked_list[i+1], y, z) ] = (*LATTICE)[lattice_index ( linked_list[i], y, z) ]; 
			
				}
				else {

					tmp_par_ptr->coords = linked_list[i+1]; 
					tmp_par_ptr_ = (*LATTICE)[ lattice_index (linked_list[i+1], y, z) ]; 
					(*LATTICE)[ lattice_index (linked_list[i+1], y, z) ] = tmp_par_ptr;
					tmp_par_ptr = tmp_par_ptr_; 
					
				}
						
			}

		}

	}

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of acceptance_after_head_regrowth
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~



// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of acceptance_after_tail_regrowth
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void acceptance_after_tail_regrowth (std::vector <Particle*>* LATTICE, \
	std::vector <std::array<int,3>>* old_cut, std::vector <std::array <int,3>>* new_cut, int y, int z){

	std::vector< std::vector<std::array<int,3>>> master_linked_list;
	std::vector <std::array<int,3>> link; 


	create_linked_list (*old_cut, *new_cut, link, &master_linked_list, 1);

	int L; // to store length of linked_list 

	Particle* tmp_par_ptr  {nullptr}; // to help with swaps of pointers
	Particle* tmp_par_ptr_ {nullptr}; // this guy as well 

	for ( std::vector <std::array<int,3>>& linked_list: master_linked_list ) {

		L = static_cast<int> (linked_list.size()) ;

		if ( (*LATTICE)[ lattice_index(linked_list[L-1], y, z) ]->ptype[0] == 's' ){

			tmp_par_ptr = (*LATTICE) [ lattice_index (linked_list.back(), y, z) ]; 
			// go backwards 
			for (int i{0}; i < L ; i=i+2 ){

				(*LATTICE)[lattice_index (linked_list[L-2-i], y, z)]->coords = linked_list[L-1-i]; 
				(*LATTICE)[lattice_index (linked_list[L-1-i], y, z)] = (*LATTICE)[lattice_index( linked_list[L-2-i], y, z )]; 
				
			}
			tmp_par_ptr->coords = linked_list[0]; 
			(*LATTICE)[lattice_index (linked_list[0], y, z)] = tmp_par_ptr; 

		}

		else { // when it is a circulation of monomers 

			for ( int i{0}; i < L; i=i+2){
				
				if ( i == 0 ) {
					// store info about linked_list[1]... and consequently also linked_list[2]
					tmp_par_ptr = (*LATTICE)[lattice_index ( linked_list[i+1], y, z) ];

					(*LATTICE)[ lattice_index ( linked_list[i],   y, z) ]->coords = linked_list[i+1]; 
					(*LATTICE)[ lattice_index ( linked_list[i+1], y, z) ] = (*LATTICE)[lattice_index ( linked_list[i], y, z) ]; 
				
				}
				else {

					tmp_par_ptr->coords = linked_list[i+1]; 
					tmp_par_ptr_ = (*LATTICE)[ lattice_index (linked_list[i+1], y, z) ]; 
					(*LATTICE)[ lattice_index (linked_list[i+1], y, z) ] = tmp_par_ptr;
					tmp_par_ptr = tmp_par_ptr_; 
					
				}
						
			}

		}	

	}

	return;
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of acceptance_after_tail_regrowth
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of TailRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void TailRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	// std::array <double,4> current_contacts = *contacts; 
	std::array <double,25> energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25> boltzmann; 
	std::array <int,25>    orientations;

	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 
	int original_ori         = (*Polymers)[p_index].chain[m_index-1]->orientation; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	/*
	std::cout << "shuffled ne_list for head regrowth = "; 
	for (std::array<int,3>& v: ne_list){
		print (v, ", ");
	}
	*/ 

	// start attempting jumps 
	int block_counter       =  0; 
	int idx_counter         =  0; 
	int self_swap_idx       = -1; 

	while ( idx_counter < 5 ){

		// std::cout << "ptype = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;

		if ( ne_list[idx_counter] == loc_m ){

			for ( int i{0}; i<5; ++i){

				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
				orientations[5*idx_counter+i] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation; 
		
				energies    [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				
				contacts_store[5*idx_counter+i] = *contacts;

				// go back to the original state
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 

			}
			 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			// check which region of the polymer the swap is taking place 

			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ){
					self_swap_idx = u; 
					break; 
				}				
			}

			// if it is with other tail units, do the swap 
			// if not, discourage it 

			// std::cout << "self_swap_idx = " << self_swap_idx << std::endl; 
			for ( int i{0}; i < 5; ++i ){
				
				if ( self_swap_idx > m_index ){
					// std::cout << "maintain current state..." << std::endl;
					// std::cout << "You have run into an undoable swap. " << std::endl; 
					// std::cout << "Selected position is "; print (ne_list[idx_counter]);	
					orientations[5*idx_counter+i] = original_ori; 
					energies[5*idx_counter+i] = 1e+08; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1}; 
					block_counter += 1; 
				}
				else {

					// (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
					(*Polymers)[p_index].chain[m_index-1]->orientation = rng_uniform(0,25);
					orientations[5*idx_counter + i] = (*Polymers)[p_index].chain[m_index-1]->orientation; 

					// prep the swap 
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

					// get the energy
					energies [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
					contacts_store[5*idx_counter+i] = *contacts; 

					// revert back to original structure 
					(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index-1]->coords   = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];

					// reset orientation 
					(*Polymers)[p_index].chain[m_index-1]->orientation = original_ori; 
				}
			}
		}
		else {

			for ( int i{0}; i < 5; ++i ){

				(*Polymers)[p_index].chain[m_index-1]->orientation = rng_uniform(0,25);
				orientations[5*idx_counter + i] = (*Polymers)[p_index].chain[m_index-1]->orientation; 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the energy
				energies [5*idx_counter + i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[5*idx_counter + i] = *contacts; 

				// revert back to original structure 
				(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];

				// reset orientation 
				(*Polymers)[p_index].chain[m_index-1]->orientation = original_ori; 

			}
		}

		idx_counter += 1;

	}

	if ( block_counter == 25 ){
		*IMP_BOOL = false;
		return; 
	}

	// std::cout << std::endl;

	// std::cout << "contacts store = "; 
	/*
	for (std::array <double,4>& n: contacts_store){
		print(n, ", ");
	}
	std::cout << std::endl;

	// now that i have all the energies, and boltzmann weights, i can choose a configuration 
	std::cout << "Energies are "; print(energies); 
	*/ 

	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<25; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	double rng_acc = rng_uniform (0.0, 1.0); 
	double rsum    = 0; 
	int e_idx      = 0; 
	// std::cout << "normalization = " << rboltzmann << std::endl;
	// std::cout << "rng = " << rng_acc << std::endl;

	for (int j{0}; j<25; ++j){
		rsum += boltzmann[j]/rboltzmann; 
		if ( rng_acc < rsum ) {
			e_idx = j; 
			break; 
		}	
	}

	// now that I have chosen a configuration, go with it
	*prob_o_to_n      = (*prob_o_to_n) * boltzmann[e_idx]/rboltzmann; 
	*frontflow_energy = energies[e_idx];
	*contacts         = contacts_store [e_idx];
	// std::cout << "e_idx = " << e_idx << std::endl;
	// std::cout << "frontflow_energy is " << *frontflow_energy << std::endl; 
	// std::cout << "position chosen = "; print (ne_list[e_idx/5]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(o->n) = " << *prob_o_to_n << std::endl;
	// std::cout << "------------------------------" << std::endl;
	// std::cout << std::endl; 

	// std::cout << "Orientations suggested = "; print (orientations); 
	// std::cout << "Assigned orientation = " << orientations[e_idx] << std::endl;

	// if { e_idx is not in the maintain index vector, perform the swap}
		
	// do the swap again 	

	(*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords             = ne_list[e_idx/5];
	(*Polymers)[p_index].chain[m_index-1]->orientation        = orientations[e_idx];
	
	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx/5], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
	
	// else { do nothing and maintain structure, because the suggested index is in the the maintain index vector }
	// std::cout << "m_index = " << m_index << std::endl;
	// (*Polymers)[0].printChainCoords(); 

	TailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void BackFlowFromTailRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::vector <int>* old_ori, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	
	std::array <double,25> energies; 
	std::array <double,25> boltzmann; 
	// std::array <double,4>  current_contacts = *contacts; 
	std::array <std::array<double,8>,25> contacts_store; 

	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; // this is the key item of interest
	int original_ori         = (*Polymers)[p_index].chain[m_index-1]->orientation; 

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[m_index-1]) - ne_list.begin() ; 

	std::array <int,25> test_ori; 

	for (int i{0}; i<25; ++i){
		if ( i == 0 ){
			test_ori[i] = (*old_ori)[m_index-1];
		}
		else {
			test_ori[i] = rng_uniform(0,25); 
		}
	}

	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	/*
	std::cout << "shuffled ne_list is "; 
	for (std::array <int,3>& n: ne_list){
		print(n, ", ");
	}
	std::cout << std::endl;
	*/ 

	// i now have a vector which has the back peddling step at position index 0 

	// start attempting jumps 
	int idx_counter         =  0; 
	int self_swap_idx 		= -1; 

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);
		if ( ne_list[idx_counter] == loc_m ){

			for (int i{0}; i < 5; ++i ){
				// std::cout << "monomer self-swap: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index-1]->orientation = test_ori[5*idx_counter+i];
				energies[5*idx_counter+i]       = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z);
				contacts_store[5*idx_counter+i] = *contacts;
			}

			(*Polymers)[p_index].chain[m_index-1]->orientation = original_ori; 
			
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			for ( int u{0}; u<deg_poly; ++u) {
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list [idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}
			}

			if ( self_swap_idx > m_index-1 ){
					
				// std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
				// std::cout << "self_swap_idx = " << self_swap_idx << std::endl;				
				// std::cout << "Selected position is "; print (ne_list[idx_counter]);
				for (int i{0}; i<5; ++i){
					(*Polymers)[p_index].chain[m_index-1]->orientation = test_ori[5*idx_counter+i];
					energies [5*idx_counter+i]      = 1e+08; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1}; 
				}
			}
			else {
				for (int i{0}; i<5; ++i){

					// std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
					// prep the swap 
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
					
					// perform the swap 
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
					(*Polymers)[p_index].chain[m_index-1]->orientation     = test_ori[5*idx_counter+i];

					// get the energy
					energies [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
					contacts_store[5*idx_counter+i] = *contacts; 
					// std::cout << "current_contacts = "; print(current_contacts); 

					// revert back to original structure 
					(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
					(*Polymers)[p_index].chain[m_index-1]->orientation       = original_ori; 
				}
			}
		}
		else {
			for (int i{0}; i<5; ++i){
				// prep the swap 
				// std::cout << "solventswap: idx_counter = " << idx_counter << std::endl;
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords  = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index-1];
				(*Polymers)[p_index].chain[m_index-1]->orientation     = test_ori[5*idx_counter+i];

				// get the energy 
				energies [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store [5*idx_counter+i] = *contacts;

				// std::cout << "current_contacts = "; print(current_contacts); 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
				(*Polymers)[p_index].chain[m_index-1]->orientation       = original_ori; 
			}
		}

		idx_counter += 1;

	}

	// std::cout << std::endl;
	/*
	std::cout << "contact_store = ";
	for ( std::array <double,4>& v: contacts_store) {
		print(v, ", ");
	}

	std::cout << "Energies are "; print(energies); 
	*/ 
	
	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<25; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	// std::cout << "normalization = " << rboltzmann << std::endl;

	*prob_n_to_o      = (*prob_n_to_o)*boltzmann[0]/rboltzmann;
	*backflow_energy  = energies[0]; 
	*contacts         = contacts_store [0];

	// std::cout << "position chosen = "; print (ne_list[0]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(n->o) = " << *prob_n_to_o << std::endl;
	// std::cout << "------------------------------" << std::endl;


	// do the swap again
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords           = ne_list[0];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ]  = (*Polymers)[p_index].chain[m_index-1];
	(*Polymers)[p_index].chain[m_index-1]->orientation = test_ori[0]; 

	// else { do nothing, MAINTAIN, because suggested index is in the maintain index vector }
	BackFlowFromTailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, old_cut, old_ori, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			Start of BackFlowFromTailRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			Start of HeadRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void HeadRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	// std::array <double,4> current_contacts = *contacts; 
	std::array <double,25> energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25> boltzmann; 
	std::array <int,25> orientations; 

	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; 
	int original_ori         = (*Polymers)[p_index].chain[m_index+1]->orientation; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	/*
	std::cout << "shuffled ne_list for head regrowth = "; 
	for (std::array<int,3>& v: ne_list){
		print (v, ", ");
	}
	std::cout << std::endl;
	*/ 
	// start attempting jumps 

	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 


	while ( idx_counter < 5 ){ // 5 locations to test 

		// std::cout << "ptype = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;

		if ( ne_list[idx_counter] == loc_m ){ 

			for (int i{0}; i < 5; ++i ){

				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
				orientations[ 5*idx_counter + i ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation;

				energies[ 5*idx_counter + i ] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 

				contacts_store[ 5*idx_counter + i ] = *contacts;

				//go back to the original state 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 

			}
			// block_counter += 1; 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			// check which region of the polymer the swap is taking place 

			for ( int u{0}; u<deg_poly; ++u ){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ){
					self_swap_idx = u; 
					break; 
				}
			}

			for (int i{0}; i < 5; ++i ){

				// if it is with other head units, do the swap 
				// if not, discourage it 

				// std::cout << "self_swap_idx = " << self_swap_idx << std::endl; 
				if ( self_swap_idx < m_index ){
					// maintain current state, and sample another state 
					// maintain_idx.push_back(idx_counter); 
					// std::cout << "In the bad zone..." << std::endl;
					energies[5*idx_counter+i] = 1e+08; // very unfavorable state 
					orientations[5*idx_counter+i] = original_ori; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1}; 
					block_counter += 1; 
				}
				else {

					// perturb orientation 
					// (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
					(*Polymers)[p_index].chain[m_index+1]->orientation = rng_uniform(0,25); 
					orientations[5*idx_counter+i] = (*Polymers)[p_index].chain[m_index+1]->orientation; 

					// prep the swap 
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
					

					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

					// get the energy
					energies [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
					contacts_store[5*idx_counter+i] = *contacts; 

					// revert back to original structure 
					(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]         			 = (*Polymers)[p_index].chain[m_index+1];

					// reset orientation
					(*Polymers)[p_index].chain[m_index+1]->orientation = original_ori; 

				}
			}
			
		}
		else {

			for (int i{0}; i < 5; ++i ){ // 5 orientations to test 

				(*Polymers)[p_index].chain[m_index+1]->orientation = rng_uniform(0,25); 
				orientations[5*idx_counter+i] = (*Polymers)[p_index].chain[m_index+1]->orientation; 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the energy
				energies [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[5*idx_counter+i] = *contacts; 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];

				// reset orientation
				(*Polymers)[p_index].chain[m_index+1]->orientation = original_ori; 

			}
		
		}

		idx_counter += 1; 
		
	}

	/*
	std::cout << "contacts store = "; 
	for (std::array <double,4>& n: contacts_store){
		print(n, ", ");
	}
	std::cout << std::endl;

	std::cout << "Energies = "; print(energies); 
	*/ 

	if ( block_counter == 25 ){
		*IMP_BOOL = false; 
		return; 
	}

	// now that i have all the energies, and boltzmann weights, i can choose a configuration 
	// std::cout << "Energies are "; print(energies); 

	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<25; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	double rng_acc = rng_uniform (0.0, 1.0); 
	double rsum    = 0; 
	int    e_idx   = 0; 
	
	// std::cout << "normalization = " << rboltzmann << std::endl;
	// std::cout << "rng = " << rng_acc << std::endl;

	for (int j{0}; j<25; ++j){
		rsum += boltzmann[j]/rboltzmann; 
		if ( rng_acc < rsum ){
			e_idx = j; 
			break; 
		}	
	}
	
	// now that I have chosen a configuration, go with it
	*prob_o_to_n     *= boltzmann[e_idx]/rboltzmann; 
	*frontflow_energy = energies[e_idx];
	*contacts         = contacts_store [e_idx]; 

	// std::cout << "e_idx = " << e_idx << std::endl;
	// std::cout << "maintain_idx = "; print(maintain_idx);
	// std::cout << "frontflow_energy is " << *frontflow_energy << std::endl; 
	// std::cout << "position chosen = "; print (ne_list[ e_idx/5 ]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(o->n) = " << *prob_o_to_n << std::endl;
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "e_idx/5 = " << e_idx/5 << std::endl;
	// std::cout << std::endl;

	// std::cout << "Orientations suggested = "; print (orientations); 
	// std::cout << "Assigned orientation = " << orientations[e_idx] << std::endl;
	
	// do the swap again
	
	(*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords             = ne_list[e_idx/5];
	(*Polymers)[p_index].chain[m_index+1]->orientation        = orientations[e_idx];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx/5], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];

	// else { do nothing and maintain structure, because suggested index is in the maintain index vector }
	HeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			end of HeadRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			Start of BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::vector <int>* old_ori, std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	
	std::array <double,25> energies; 
	std::array <double,25> boltzmann; 
	// std::array <double,4> current_contacts = *contacts; 
	std::array <std::array<double,8>,25> contacts_store; 

	double rboltzmann = 0; // running sum for boltzmann weights 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; // this is the key item of interest
	int original_ori         = (*Polymers)[p_index].chain[m_index+1]->orientation; 

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[recursion_depth]) - ne_list.begin() ; 
	// std::cout << "old_ori = "; print(*old_ori);

	std::array <int,25> test_ori; 
	for (int i{0}; i<25; ++i){
		if ( i == 0 ){
			test_ori[i] = (*old_ori)[recursion_depth];
		}
		else {
			test_ori[i] = rng_uniform(0,25); 
		}
	}

	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	/*
	std::cout << "shuffled ne_list is "; 
	for (std::array <int,3>& n: ne_list){
		print(n, ", ");
	}
	std::cout << std::endl;
	*/
	// i now have a vector which has the back peddling step at position index 0 
	// std::cout << "test_ori = "; print (test_ori);

	// start attempting jumps 
	int idx_counter         = 0; 
	int self_swap_idx 		= -1; 

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);
		if ( ne_list[idx_counter] == loc_m ){

			for (int i{0}; i < 5; ++i ){
				// std::cout << "self-place: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index+1]->orientation = test_ori[5*idx_counter+i]; 
				energies[5*idx_counter+i]       = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z);
				contacts_store[5*idx_counter+i] = *contacts;
			}

			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 
			

			// std::cout << "current_contacts = "; print(current_contacts); 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list [idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}
			}

			if ( self_swap_idx < m_index ){
				
				// std::cout << "Selected position is "; print (ne_list[idx_counter]);
				for (int i{0}; i<5; ++i){
					// std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
					// std::cout << "self_swap_idx = " << self_swap_idx << std::endl;
					(*Polymers)[p_index].chain[m_index+1]->orientation = test_ori[5*idx_counter+i];
					energies [5*idx_counter+i]      = 1e+08; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1}; 
				}
			}
			else {
				for (int i{0}; i<5; ++i){
					// std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
					
					// prep the swap 
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
					
					// perform the swap 
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
					(*Polymers)[p_index].chain[m_index+1]->orientation     = test_ori[5*idx_counter+i];

					// get the energy
					energies [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
					contacts_store[5*idx_counter+i] = *contacts;

					// std::cout << "current_contacts = "; print(current_contacts); 

					// revert back to original structure 
					(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];
					(*Polymers)[p_index].chain[m_index+1]->orientation       = original_ori; 
				}
			}
		}
		else {
			for (int i{0}; i<5; ++i){
				// prep the swap 
				// std::cout << "solventswap: idx_counter = " << idx_counter << std::endl;
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index+1];
				(*Polymers)[p_index].chain[m_index+1]->orientation     = test_ori[5*idx_counter+i];

				// get the energy 
				energies [5*idx_counter+i] = CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store [5*idx_counter+i] = *contacts;
				// std::cout << "current_contacts = "; print(current_contacts); 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];
				(*Polymers)[p_index].chain[m_index+1]->orientation       = original_ori; 
			}
		}

		idx_counter += 1;

	}

	// std::cout << std::endl;

	// std::cout << "contact_store = ";
	// for ( std::array <double,4>& v: contacts_store) {
		// print(v, ", ");
	// }
	// std::cout << "Energies are "; print(energies); 
	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<25; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	// std::cout << "normalization = " << rboltzmann << std::endl;

	*prob_n_to_o      = (*prob_n_to_o)*boltzmann[0]/rboltzmann;
	*backflow_energy  = energies[0]; 
	*contacts         = contacts_store [0];
	// std::cout << "position chosen = "; print (ne_list[0]);
	// std::cout << "------------------------------" << std::endl;
	// std::cout << "p(n->o) = " << *prob_n_to_o << std::endl;
	// std::cout << "------------------------------" << std::endl;


	// do the swap again
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords           = ne_list[0];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];
	(*Polymers)[p_index].chain[m_index+1]->orientation = test_ori[0]; 

	// else { do nothing, MAINTAIN, because suggested index is in the maintain index vector }
	BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, old_cut, old_ori, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			End of BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// Start of ChainRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void ChainRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int p_index, int x, int y, int z){


	int deg_poly = (*Polymers)[p_index].deg_poly; 
	int m_index  = rng_uniform (1, deg_poly-2); 

	std::array <double,8> c1_contacts = *contacts; 
	std::array <double,8> c2_contacts = *contacts; 

	double prob_o_to_n {1}; 
	double prob_n_to_o {1}; 
	double frontflow_energy {*sysEnergy}; 
	double backflow_energy  {0}; 
	int    recursion_depth  {0}; 
	// old_cut contains the old positions of thr monomer 
	std::vector <std::array<int,3>> old_cut; 
	std::vector <std::array<int,3>> new_cut;
	std::vector <int> old_ori;
	std::vector <int> new_ori;
	old_cut.reserve(deg_poly);
	old_ori.reserve(deg_poly);
	new_cut.reserve(deg_poly);
	new_ori.reserve(deg_poly);  

	// regrowth will be in the direction where there are fewer monomer
	// if m_index/deg_poly is lesser than 0.5, growth is false, otherwise growth is true 
	int growth {-1}; 

	if ( deg_poly % 2 == 0 ){
		growth = (0.5 >= (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
	}
	else {
		if ( 0.5 == (m_index+1)/static_cast<double>(deg_poly+1) ){
			growth = rng_uniform (0, 1);
		}
		else {
			growth = (0.5 > (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
		}
	}
	
	if ( growth ){
		
		// std::cout << "Head regrowth... " << std::endl;
		// std::cout << "OLD ENERGY = " << *sysEnergy << std::endl << std::endl;
		// std::cout << "initial orientation of final bead = " << (*Polymers)[p_index].chain[3]->orientation << std::endl;
		// get old_cut 
		for (int i {m_index+1}; i<deg_poly; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
			old_ori.push_back ((*Polymers)[p_index].chain[i]->orientation) ; 
		}

		// std::cout << "old orientation = "; print (old_ori);
		// std::cout << "printing out neighbors with orientations..." << std::endl;

		/*
		t_ne_list = obtain_ne_list ( (*Polymers)[p_index].chain[deg_poly-1]->coords, x, y, z ); 

		std::cout << "orientation dump - pre regrowth... " << std::endl;
		for ( std::array <int,3>& t: t_ne_list){
			print (t, ", "); std::cout << "orientation = " << (*LATTICE)[ lattice_index(t, y, z) ]->orientation << std::endl;
		}
		*/ 

		// regrow the polymer frontwards
		HeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		
		// CheckStructures (Polymers, LATTICE, &c1_solvation_shells, x, y, z); 

		/*
		std::cout << "orientation dump - post regrowth... " << std::endl;
		t_ne_list = obtain_ne_list ( (*Polymers)[p_index].chain[deg_poly-1]->coords, x, y, z ); 

		for ( std::array <int,3>& t: t_ne_list){
			print (t, ", "); std::cout << "orientation = " << (*LATTICE)[ lattice_index(t, y, z) ]->orientation << std::endl;
		}


		std::cout << "final orientation of final bead = " << (*Polymers)[p_index].chain[3]->orientation << std::endl;
		*/ 

		for (int i {m_index+1}; i<deg_poly; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
			new_ori.push_back ((*Polymers)[p_index].chain[i]->orientation) ;
		}

		// std::cout << "new orientation = "; print (new_ori);

		if ( old_cut == new_cut && old_ori == new_ori){
			// std::cout << "-------------------------------------------------------------------------"  << std::endl;
			// std::cout << "POLYMER CONFIGURATION WAS NOT CHANGED. RETURN BACK TO MAIN CONFIGURATION." << std::endl;
			// std::cout << "-------------------------------------------------------------------------"  << std::endl << std::endl;
			return; 
		}

		if ( !(*IMP_BOOL) ){
			// std::cout << "ALL BLOCKS! REVERTING!" << std::endl;
			// revert back to the original state. 
			// this is actually a rejection move 

			acceptance_after_head_regrowth ( LATTICE, &new_cut, &old_cut, y, z );

			// change orientations of polymer bead to old 
			for (int i{m_index+1}; i<deg_poly; ++i){
				(*Polymers)[p_index].chain[i]->orientation = old_ori[i-m_index-1]; 
			}

			return; 
		}

		// std::cout << "frontflow energy = " << frontflow_energy << std::endl;
		// std::cout << "front flow contacts are "; print (c1_contacts);
		// std::cout << "Coordinates of perturbed polymer: " << std::endl;
		// (*Polymers)[p_index].printChainCoords(); 

		backflow_energy = frontflow_energy; 
		// std::cout << "Begin backflow... " << std::endl;
		
		BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, &old_ori, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

		// CheckStructures (Polymers, LATTICE, &c2_solvation_shells, x, y, z); 

		/*
		t_ne_list = obtain_ne_list ( (*Polymers)[p_index].chain[deg_poly-1]->coords, x, y, z ); 
		std::cout << "orientation dump - post backflow... " << std::endl;
		for ( std::array <int,3>& t: t_ne_list){
			print (t, ", "); std::cout << "orientation = " << (*LATTICE)[ lattice_index(t, y, z) ]->orientation << std::endl;
		}

		std::cout << "post backflow orientation of final bead = " << (*Polymers)[p_index].chain[deg_poly-1]->orientation << std::endl;
		

		std::cout << "Coordinates of polymer after backflow: " << std::endl;
		for (int j{0}; j< int((*Polymers)[0].chain.size()); ++j){
            print((*Polymers)[0].chain[j]->coords, ", "); std::cout << "o = " << (*Polymers)[0].chain[j]->orientation << std::endl;
        }
		*/ 

		if ( *sysEnergy != backflow_energy || c2_contacts != *contacts){
			std::cout << "Energies are bad, or contacts are not right." << std::endl;
			std::cout << "*sysEnergy = " << *sysEnergy << ", backflow_energy = " << backflow_energy << "." << std::endl;
			std::cout << "c2_contacts = "; print (c2_contacts, ", "); std::cout << "*contacts = "; print(*contacts);
			std::cout << "Shit's fucked." << std::endl;
			exit(EXIT_FAILURE);
		}

		// check acceptance criterion
		double rng_acc = rng_uniform (0.0, 1.0); 
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n ){

			acceptance_after_head_regrowth (LATTICE, &old_cut, &new_cut, y, z); 
			for (int i{m_index+1}; i<deg_poly; ++i){
				(*Polymers)[p_index].chain[i]->orientation = new_ori[i-m_index-1]; 
				// std::cout << "new_ori[" << i-m_index-1 << "] = " << new_ori[i-m_index-1] << std::endl;
			}
			*sysEnergy = frontflow_energy; 
			*contacts  = c1_contacts;

		}

		else {
			*IMP_BOOL = false; 
		}

	}
	else {
		

		// std::cout << "Initiate tail regrowth..." << std::endl;
		// get old cut 
		for ( int i{0}; i<m_index; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords);
			old_ori.push_back ((*Polymers)[p_index].chain[i]->orientation);
		}
		// std::cout << "old orientation = "; print (old_ori);
		// regrow the polymer backwards
		
		/*
		t_ne_list = obtain_ne_list ( (*Polymers)[p_index].chain[m_index-1]->coords, x, y, z ); 

		std::cout << "orientation dump - pre regrowth... " << std::endl;
		for ( std::array <int,3>& t: t_ne_list){
			print (t, ", "); std::cout << "orientation = " << (*LATTICE)[ lattice_index(t, y, z) ]->orientation << std::endl;
		}
		*/ 

		TailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		/*
		std::cout << "Performed regrowth... " << std::endl;
		for (int j{0}; j< int((*Polymers)[0].chain.size()); ++j){
            print((*Polymers)[0].chain[j]->coords, ", "); std::cout << "o = " << (*Polymers)[0].chain[j]->orientation << std::endl;
        }
		*/
		// CheckStructures (Polymers, LATTICE, &c1_solvation_shells, x, y, z); 
		/*
		t_ne_list = obtain_ne_list ( (*Polymers)[p_index].chain[m_index-1]->coords, x, y, z ); 

		std::cout << "orientation dump - post regrowth... " << std::endl;
		for ( std::array <int,3>& t: t_ne_list){
			print (t, ", "); std::cout << "orientation = " << (*LATTICE)[ lattice_index(t, y, z) ]->orientation << std::endl;
		}

		std::cout << "Performed lattice checks..." << std::endl; 
		*/ 

		for (int i {0}; i<m_index; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
			new_ori.push_back ((*Polymers)[p_index].chain[i]->orientation) ;
		}

		// std::cout << "new orientation = "; print (new_ori);

		// std::cout << "created new_cut and new_ori!" << std::endl; 

		if ( old_cut == new_cut && old_ori == new_ori ){
			// std::cout << "------------------------------------------------------------------------"  << std::endl;
			// std::cout << "POLYMER CONFIGURATION WAS NOT CHANGED. RETURN BACK TO MAIN."               << std::endl;
			// std::cout << "------------------------------------------------------------------------"  << std::endl << std::endl;
			return; 
		}

		
		if ( !(*IMP_BOOL) ){
			// std::cout << "ALL BLOCKS! REVERTING!" << std::endl;
			acceptance_after_tail_regrowth ( LATTICE, &new_cut, &old_cut, y, z); 
			
			for (int i{0}; i<m_index; ++i){
				(*Polymers)[p_index].chain[i]->orientation = old_ori[i]; 
			}
			
			return; 
		}

		backflow_energy = frontflow_energy; 

		// std::cout << "BEGIN BACK FLOW for Tail! " << std::endl;

		BackFlowFromTailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, &old_ori, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

		// t_ne_list = obtain_ne_list ( (*Polymers)[p_index].chain[m_index-1]->coords, x, y, z ); 

		/*
		std::cout << "orientation dump - post backflow... " << std::endl;
		for ( std::array <int,3>& t: t_ne_list){
			print (t, ", "); std::cout << "orientation = " << (*LATTICE)[ lattice_index(t, y, z) ]->orientation << std::endl;
		}

		std::cout << "Coordinates of polymer post backflow: " << std::endl;
		for (int j{0}; j< int((*Polymers)[0].chain.size()); ++j){
            print((*Polymers)[0].chain[j]->coords, ", "); std::cout << "o = " << (*Polymers)[0].chain[j]->orientation << std::endl;
        }

		std::cout << "Performed BACK FLOW for Tail! " << std::endl;
		*/ 

		if ( *sysEnergy != backflow_energy || c2_contacts != *contacts ){
			std::cout << "Energies are bad, or contacts are not right." << std::endl;
			std::cout << "*sysEnergy = " << *sysEnergy << ", backflow_energy = " << backflow_energy << "." << std::endl;
			std::cout << "c2_contacts = "; print (c2_contacts, ", "); std::cout << "*contacts = "; print(*contacts);
			std::cout << "Shit's fucked." << std::endl;
			exit(EXIT_FAILURE);
		}

		// CheckStructures (x, y, z, Polymers, LATTICE); 

		// check acceptance criterion 
		double rng_acc = rng_uniform (0.0, 1.0); 
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy- *sysEnergy)) * prob_n_to_o/prob_o_to_n ){
			// accept new cut 
			// perform swaps 
			acceptance_after_tail_regrowth (LATTICE, &old_cut, &new_cut, y, z); 

			for (int i{0}; i<m_index; ++i){
				(*Polymers)[p_index].chain[i]->orientation = new_ori[i]; 
			}

			*sysEnergy        = frontflow_energy; 
			*contacts         = c1_contacts;
		}
		else {
			*IMP_BOOL = false; 
		}

	}

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			   End of ChainRegrowthPlusOrientationalFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ##############################################################################################

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of end rotation with solvent flips with a bias 
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

/*
void TailRotationWithSolventFlips_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int p_index, int x, int y, int z){

	// get the neighborlist of particles at index 1 
	std::array <int,3> loc_0 = (*Polymers)[p_index].chain[0]->coords;
	std::array <int,3> loc_1 = (*Polymers)[p_index].chain[1]->coords; 

	std::array <std::array <int,3>, 26> ne_list = obtain_ne_list (loc_1, x, y, z); 

	// std::array <double,4> current_contacts = *contacts; 
	std::array <double,5>               energies; 
	std::array <std::array<double,4>,5> contacts_store; 
	std::array <double,5>               boltzmann;
	std::array <int,5>				    old_orientations; 
	std::array <int,5>				    new_orientations; 
	double rboltzmann   = 0; 
	double prob_o_to_n  = 1; 
	double prob_n_to_o  = 1; 
	
	double Emin  = 0; 
	double rng   = 0; 
	double rsum  = 0; 
	int    e_idx = 0; 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// start checking your choices
	
	int block_counter = 0 ;
	int self_swap_idx = -1; 

	for (int i{0}; i < 5; ++i){

		if ( (*LATTICE)[ lattice_index (ne_list[i], y, z) ]->ptype[0] == 's' ){

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[i], y, z)]->coords = loc_0;
			(*Polymers)[p_index].chain[0]->coords               = ne_list[i];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_0, y, z) ] = (*LATTICE)[lattice_index(ne_list[i], y, z)];
			(*LATTICE)[ lattice_index (ne_list[i], y, z) ] = (*Polymers)[p_index].chain[0];


			// start evaluating boltzmann factors for each position, much like for regular tail rotation 
			energies[i] = CalculateEnergy (Polymers, LATTICE, E, &(contacts_store[i]), x, y, z);  

			// revert the polymer to original structure 
			(*LATTICE)[lattice_index(loc_0, y, z)]->coords = ne_list[i];
			(*Polymers)[p_index].chain[0]->coords  = loc_0;
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[lattice_index(ne_list[i], y, z)]    = (*LATTICE)[lattice_index (loc_0, y, z)];
			(*LATTICE)[lattice_index(loc_0, y, z)]         = (*Polymers)[p_index].chain[0];	


		}

		else {
			// do nothing 
			energies[i] = 1e+08;        // very unfavorable state 
			contacts_store[i] = {-1,-1,-1,-1}; 
			block_counter += 1; 

		}

	}

	Emin = *std::min_element ( energies.begin(), energies.end() ); 

	for ( int k{0}; k<5; ++k ){
		boltzmann [k] = std::exp (-1/temperature* (energies[k] - Emin) ); 
		rboltzmann   += boltzmann [k]; 
	}

	rng   = rng_uniform (0.0, 1.0); 
	// pick a location 
	for ( int j{0}; j<5; ++j){
		rsum += boltzmann[j]/rboltzmann;
		if ( rng < rsum ){
			e_idx = j; 
			break; 
		}
	}

	// prep the swap 
	(*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_0;
	(*Polymers)[p_index].chain[0]->coords               = ne_list[e_idx];
	
	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_0, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ] = (*Polymers)[p_index].chain[0];

	prob_o_to_n *= boltzmann[e_idx]/rboltzmann; 

	// now perform the solvent flips 
	new_orientations = {0,0,0,0,0}; 
	std::array <int,5> solvent_indices = {0,0,0,0,0}; 
	std::array <int,4> new_contacts = {0,0,0,0};
	double new_energy = 0; 

	SolventFlip_BIASED (Polymers, LATTICE, E, &new_contacts, &solvent_indices, &new_orientations, IMP_BOOL, sysEnergy, &new_energy, &prob_o_to_n, &prob_n_to_o, temperature, 1, p_index, x, y, z); 

	// now do the backtracking wrt monomer swinging 
	// *contacts = new_contacts; // 


	// I AM NOT SURE WHAT IS GOING ON HERE 
	// int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*Polymers)[p_index].chain[0]->coords) - ne_list.begin() ; 
	// std::array <int,3> tmp = ne_list[0]; 
	// ne_list[0]             = ne_list[ne_idx];
	// ne_list[ne_idx]        = tmp; 

	return; 
}
*/

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of acceptance_after_tail_regrowth
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
/*
void SolventFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::array <double,4>* E, std::array <double,4>* new_contacts, std::array<int,5>* solvent_indices, std::array <int,5>* new_orientations, \
	bool* IMP_BOOL, double* sysEnergy, double* new_energy, \
	double* prob_o_to_n, double* prob_n_to_o, double temperature, \
	int m_index, int p_index,  int x, int y, int z){

	std::array <std::array <int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 
	
	// to get rid of! 
	// temporary variables, to delete later 
	std::array <double,4> c_contacts = *new_contacts; 

	// #################################################

	std::array <std::array <double,4>, 5> contacts_store_1;
	std::array <int,5> old_orientations;				// need this for the unflipping process
	std::array <int,5> suggested_orientations; 			// need this for holding orientations before assigning them 
	std::array <double,5> energies; 
	std::array <double,5> boltzmann;
	// std::array <int,5> new_orientations; 				// holds the best possible orientations 
	double Emin       = 0; 
	double rboltzmann = 0; 
	double rng        = 0; 
	int    rsum       = 0;
	int    e_idx      = -1;

	for (int i{0}; i<5; ++i){

		(*solvent_indices)[i] = lattice_index (ne_list[i], y, z); // assignment, in case this move gets accepted, i need this information 
		old_orientations[i]  = (*LATTICE)[ lattice_index (ne_list[i], y, z) ]->orientation; 

		for (int j{0}; j<5; ++j){
			(*LATTICE)[ lattice_index (ne_list[i], y, z) ]->orientation = rng_uniform(0, 25); 
			energies[j] = CalculateEnergy (Polymers, LATTICE, E, &(contacts_store_1[i]), x, y, z); 
			suggested_orientations [j] = (*LATTICE)[ lattice_index (ne_list[i], y, z) ]->orientation; 
		}

		// revert back to old orientation for index i 
		// (*LATTICE)[ lattice_index (ne_list[i], y, z) ]->orientation = old_orientations[i]; 

		// calculate boltzmann weights 
		Emin = *std::min_element (energies.begin(), energies.end()); 
		rboltzmann = 0; 
		rsum = 0; 

		for ( int k{0}; k<5; ++k ){
			boltzmann [k] = std::exp (-1/temperature* (energies[k] - Emin) ); 
			rboltzmann   += boltzmann [k]; 
		}

		rng   = rng_uniform (0.0, 1.0); 
	
		// pick a location 
		for ( int k{0}; k<5; ++k){
			rsum += boltzmann[k]/rboltzmann;
			if ( rng < rsum ){
				e_idx = k; 
				break; 
			}
		}
		*prob_o_to_n *= boltzmann[e_idx]/rboltzmann; 

		// make the jump 
		(*LATTICE)[ lattice_index (ne_list[i], y, z) ]->orientation = suggested_orientations[e_idx]; 
		(*new_orientations) [i] = suggested_orientations[e_idx]; 
		
	}

	(*new_contacts) = contacts_store_1[e_idx]; 
	(*new_energy)   = energies[e_idx]; 

	// i now have prob_o_to_n 
	// time to get prob_n_to_o 

	for (int i{0}; i < 5; ++i){

		(*LATTICE)[ lattice_index(ne_list[i], y, z) ]->orientation = old_orientations[i]; 
		energies[0] = CalculateEnergy (Polymers, LATTICE, E, &(contacts_store_1[0]), x, y, z);

		for (int j{1}; j < 5; ++j){
			(*LATTICE)[ lattice_index (ne_list[i], y, z) ]->orientation = rng_uniform(0, 25); 
			energies[j] = CalculateEnergy (Polymers, LATTICE, E, &(contacts_store_1[j]), x, y, z); 
		}

		// calculate boltzmann weights 
		Emin = *std::min_element (energies.begin(), energies.end()); 
		rboltzmann = 0; 

		for ( int k{0}; k<5; ++k ){
			boltzmann [k] = std::exp (-1/temperature * (energies[k] - Emin ) );
			rboltzmann   += boltzmann [k]; 
		}

		*prob_n_to_o *= boltzmann[0]/rboltzmann;

	}

	// TO GET RID OF
	// can get rid of sysenergy and c_contacts

	if ( *sysEnergy != energies[0] || c_contacts != contacts_store_1[0] ){
		std::cout << "Energies are bad, or contacts are not right. " << std::endl;
		std::cout << "*sysEnergy = " << *sysEnergy << ", energies[0] = " << energies[0] << "." << std::endl;
		std::cout << "*contacts = "; print(*new_contacts, ", "); std::cout << "contacts_store_1[0] = "; print (contacts_store_1[0]); 
		exit(EXIT_FAILURE);
	}

	return; 

}
*/

/*
void HeadRotationWithSolventFlips_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::array <double,4>* E, std::array <double,4>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int p_index, int x, int y, int z);

void EndRotationWithSolventFlips_BIASED(std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::array <double,4>* E, std::array <double,4>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int p_index, int x, int y, int z){

	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = 1; // distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        TailRotationWithSolventFlips_BIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        HeadRotationWithSolventFlips_BIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }


	return;
} 

void RegrowthWithSolventFlips_BIASED(); 
*/


void SolventExchange_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::array <double,8>* E, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int x, int y, int z) {

	// int deg_poly = static_cast<int>( (*Polymers)[0].chain.size() ); 

	
	// solvation_shell_indices.reserve (26*26*deg_poly); 
	std::set <int> solvation_shell_set; 
	std::array <std::array<int,3>, 26> ne_list, ne_list_; 

	// get the first solvation shell 
	for ( Polymer& pmer: *Polymers){
		for ( Particle*& p: pmer.chain ){

			ne_list = obtain_ne_list ( p->coords, x, y, z); 
			for ( std::array <int,3>& loc: ne_list ){
				if ( (*LATTICE)[ lattice_index (loc, y, z) ]->ptype[0] == 's' ){
					solvation_shell_set.insert (lattice_index (loc, y, z)); 
					ne_list_ = obtain_ne_list (loc, x, y, z);
					for ( std::array <int,3>& loc_: ne_list_){
						if ( (*LATTICE)[ lattice_index (loc_, y, z) ]->ptype[0] == 's' ){
							solvation_shell_set.insert ( lattice_index (loc_, y, z) ); 
						}
					}
				}
			}
		}
	}

	std::vector <int> solvation_shell_indices (solvation_shell_set.begin(), solvation_shell_set.end()); 

	// std::cout << "Is this being hit1?" << std::endl;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// std::iota (solvation_shell_indices.begin(), solvation_shell_indices.end(), 0); 
	std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));

	int nexchange =  1; 
	int exc_idx   = -1; 

    int                                 ntest            = 5; 
    std::array <double,5>               energies         = {0,0,0,0,0}; 
    std::array <double,5>               boltzmann        = {0,0,0,0,0};
    std::array <int,5>                  orientations     = {0,0,0,0,0}; 
    double                              prob_o_to_n      = 1; 
    double                              prob_n_to_o      = 1; 
    double                              frontflow_energy = 0; 
	std::array<std::array<double,8>,5>  contacts_store   = {*contacts, *contacts, *contacts, *contacts, *contacts}; 
    std::array <double,8>               c_contacts1      = *contacts; 
    std::array <double,8>				c_contacts2      = *contacts; 
    double                              Emin             = 0; 
    double                              rng_acc          = 0;
    double                              rng              = 0; 
    double                              rsum             = 0; 
    double                              rboltzmann       = 0; 
    int                                 e_idx            = 0; 
    

	Particle* tmp_par_ptr {nullptr};
	std::vector <int> old_ori; 
	std::vector <int> new_ori; 
	old_ori.reserve(nexchange);  
	new_ori.reserve(nexchange);

	// std::vector <int> exc_ids;
	// exc_ids.reserve(nexchange);  


	exc_idx = rng_uniform (0, x*y*z-1); 
	int my_idx = 0; 

	// switch positions, of exc idx and solvation_shell[i], then perturb orientation of exc_idx on switch 

	if ( (*LATTICE)[exc_idx]->ptype[0] == 's' ){
		// made a copy of the particle in the solvation shell 
		// swap particles 

		tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[ my_idx ] ]; 
		old_ori.push_back ((*LATTICE)[ exc_idx ]->orientation); 
		
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ] = (*LATTICE)[ exc_idx ]; 
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ]->coords = location ( (solvation_shell_indices)[ my_idx ], x, y, z); 

		(*LATTICE)[ exc_idx ] = tmp_par_ptr; 
		(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); 

		// flip the particle newly added to the solvation shell 
		for (int j{0}; j < ntest; ++j) { 

			(*LATTICE)[ solvation_shell_indices[ my_idx ] ]->orientation = rng_uniform (0, 25); 
			orientations   [ j ] = (*LATTICE)[ (solvation_shell_indices)[ my_idx ] ]->orientation; 
			energies       [ j ] = CalculateEnergy_parallel ( Polymers, Cosolvent, LATTICE, E, &c_contacts1, x, y, z ); 
			contacts_store [ j ] = c_contacts1; 
		}
		

	}

	else if ( (*LATTICE)[exc_idx]->ptype == "m1" ){
		(*IMP_BOOL) = false; 
		return; 
	}

	// std::cout << "Is this being hit3?" << std::endl;

	Emin = *std::min_element ( energies.begin(), energies.end() ); 
	rboltzmann = 0; 

	for (int k{0}; k<5; ++k){
		boltzmann[k] = std::exp ( -1/temperature*( energies[k] - Emin )  );
		rboltzmann  += boltzmann[k];
	}

	rng       = rng_uniform (0.0, 1.0); 
	rsum      = 0; 
	e_idx     = 0; 

	for (int j{0}; j<5; ++j){
		rsum += boltzmann[j]/rboltzmann; 
		if ( rng < rsum ) {
			e_idx = j;
			break; 
		}	
	}

	prob_o_to_n = boltzmann[e_idx]/rboltzmann; 

	// once I know which orientation is good, assign it 
	(*LATTICE)[ solvation_shell_indices[ my_idx ] ]->orientation = orientations[e_idx]; 
	new_ori.push_back (orientations[e_idx]); // the chosen orientation


	frontflow_energy     = energies[e_idx]; // the chosen energy
	c_contacts1          = contacts_store [e_idx]; // the chosen contacts

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// reversing the flow 

	// made a copy of the particle in the solvation shell 
	// swap particles 

	tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[my_idx] ];
	
	(*LATTICE)[ solvation_shell_indices [my_idx] ]        = (*LATTICE)[ exc_idx ]; 
	(*LATTICE)[ solvation_shell_indices[my_idx] ]->coords = location ( solvation_shell_indices[my_idx], x, y, z); 

	(*LATTICE)[ exc_idx ]         = tmp_par_ptr; 
	(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); // this particle is the particle that was in the solvation shell 

	// flip the particle sent away 
	(*LATTICE) [ exc_idx ]->orientation = old_ori[0]; 
	energies [ 0 ]       = CalculateEnergy_parallel ( Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z ); 

	if ( energies[0] != *sysEnergy || c_contacts2 != *contacts ){
    	std::cout << "Something is fucked. Energies do not match." << std::endl;
    	std::cout << "backflow_energy = " << energies[0] << ", sysEnergy = " << *sysEnergy << std::endl;
    	std::cout << "c_contacts2 = "; print (c_contacts2, ", "); std::cout << "*contacts = "; print (*contacts); 
    	exit(EXIT_FAILURE); 
    }

	for (int j{1}; j < ntest; ++j){

		(*LATTICE)[ exc_idx ]->orientation = rng_uniform (0, 25); 
		energies  [ j ]      = CalculateEnergy_parallel ( Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z ); 

	}

	for (int k{0}; k < ntest; ++k){
		boltzmann[k] = std::exp ( -1/temperature*( energies[k] - Emin )  );
		rboltzmann  += boltzmann[k];
	}

	prob_n_to_o = boltzmann[0]/rboltzmann; 
	rng_acc = rng_uniform (0.0, 1.0); 

	if ( rng_acc < std::exp ( -1/temperature * (frontflow_energy - *sysEnergy) ) * prob_n_to_o/prob_o_to_n ) {

		// make the swap again 
		tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[my_idx] ];
		
		(*LATTICE)[ solvation_shell_indices[my_idx] ] = (*LATTICE)[ exc_idx ]; 
		(*LATTICE)[ solvation_shell_indices[my_idx] ]->coords = location ( solvation_shell_indices[my_idx], x, y, z); 
		(*LATTICE)[ solvation_shell_indices[my_idx] ]->orientation = new_ori[0]; 

		(*LATTICE)[ exc_idx ] = tmp_par_ptr; 
		(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); 
		*sysEnergy = frontflow_energy;
		*contacts  = c_contacts1;
	}
	else {
		*IMP_BOOL = false;
		// get the particle back to its original orientation 
		(*LATTICE) [ exc_idx ]->orientation = old_ori[0]; 
	}	


	return; 

}


//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: PerturbSystem_UNBIASED, PerturbSystem_BIASED
//
// PARAMETERS: Polymers, LATTICE, Energies, contacts, attempts, IMP_BOOL, v, currentenergy, temperature, 
// move number, box dimensions  
//
// WHAT THE FUNCTION DOES: This is the main simulation function. It takes all of these inputs and perturbs the system. 
// This is the heart of the monte carlo integrator.  
// 
// DEPENDENCIES: every move used above. 
// 
// THE CODE: 

void PerturbSystem_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, std::array <int,9>* attempts, \
	bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, \
	int* move_number, int x, int y, int z){

	int index = 0; // rng_uniform (0, static_cast<int>((*Polymers).size()-1) ); 
	int r     = 2; // rng_uniform (0, 4); 

	switch (r) {

		case(0):
			if (v) {
				std::cout << "Performing end rotations..." << std::endl; 
			}
			EndRotation_UNBIASED  (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = 0; 
			(*attempts)[0] += 1; 
			break; 

		case (1):
			if (v) {
				std::cout << "Performing reptation..." << std::endl; 
			}
			Reptation_UNBIASED    (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = 1;
			(*attempts)[1] += 1; 
			break;

		case (2):
			if (v) {
				std::cout << "Performing chain regrowth..." << std::endl;
			}
			ChainRegrowth_UNBIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = 2; 
			(*attempts)[2] += 1;
			break; 

		case(3):
			if (v) {
				std::cout << "Performing solvent flips..." << std::endl;
			}
			SolventFlip_UNBIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = 3; 
			(*attempts)[3] += 1; 
			break;

		case (4):
			if (v) {
				std::cout << "Performing monomer flips..." << std::endl;
			}
			PolymerFlip_UNBIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = 4;
			(*attempts)[4] += 1; 
			break; 

	}

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of PerturbSystem_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void PerturbSystem_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::array <double,8>* E, std::array <double,8>* contacts, std::array <int,9>* attempts, \
	bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, \
	int* move_number, int x, int y, int z) {

	int index = 0; // rng_uniform (0, static_cast<int>((*Polymers).size()-1) ); 
	int r     = rng_uniform (0, 7); 
	// std::array <double,4> c_contacts = *contacts; 

	switch (r) {

		case(0):
			if (v) {
				std::cout << "Performing end rotations..." << std::endl; 
			}
			EndRotation_UNBIASED  (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break; 

		case (1):
			if (v) {
				std::cout << "Performing reptation..." << std::endl; 
			}
			Reptation_UNBIASED    (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (2):
			if (v) {
				std::cout << "Performing biased chain regrowth..." << std::endl;
			}
			ChainRegrowth_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

		case (3):
			if (v) {
				std::cout << "Performing biased chain regrowth with orientation flip..." << std::endl; 
			}
			ChainRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

		case (4):
			if (v) {
				std::cout << "Performing solvent flips..." << std::endl;
			}
			SolventFlip_UNBIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break;

		case (5):
			if (v) {
				std::cout << "Performing a biased solvation shell flip..." << std::endl;
			}
			SolvationShellFlip_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1;
			break;

		case (6):
			if (v) {
				std::cout << "Performing a biased polymer flip..." << std::endl;
			}
			PolymerFlip_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (7): 
			if (v) {
				std::cout << "Performing an identity swap... " << std::endl;
			}
			SolventExchange_BIASED (Polymers, Cosolvent, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z);
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;			

	}

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of PerturbSystem_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of PerturbSystem
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
