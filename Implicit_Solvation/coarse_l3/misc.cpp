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
std::array <std::array <int,3>, 26> adrns  = { ax, ay, az, nx, ny, nz, axay, axaz, axny, axnz, nxay, nxaz, nxny, nxnz, ayaz, aynz, nyaz, nynz, axayaz, axnyaz, axaynz, axnynz, nxayaz, nxaynz, nxnyaz, nxnynz }; 
std::array <std::array <int,3>, 98> nadrns = { { { 2, 0, 0 }, { 2, 1, 0 }, { 2, -1, 0 }, { 2, 0, 1 }, { 2, 0, -1 }, { 2, 1, 1 }, { 2, -1, 1 }, { 2, 1, -1 }, { 2, -1, -1 }, { -2, 0, 0 }, { -2, 1, 0 }, \
{ -2, -1, 0 }, { -2, 0, 1 }, { -2, 0, -1 }, { -2, 1, 1 }, { -2, -1, 1 }, { -2, 1, -1 }, { -2, -1, -1 }, { 0, 2, 0 }, { 1, 2, 0 }, { -1, 2, 0 }, { 0, 2, 1 }, { 0, 2, -1 }, { 1, 2, 1 }, { -1, 2, 1 }, \
{ 1, 2, -1 }, { -1, 2, -1 }, { 0, -2, 0 }, { 1, -2, 0 }, { -1, -2, 0 }, { 0, -2, 1 }, { 0, -2, -1 }, { 1, -2, 1 }, { -1, -2, 1 }, { 1, -2, -1 }, { -1, -2, -1 }, { 0, 0, 2 }, { 0, 1, 2 }, { 0, -1, 2 }, \
{ 1, 0, 2 }, { -1, 0, 2 }, { 1, 1, 2 }, { -1, 1, 2 }, { 1, -1, 2 }, { -1, -1, 2 }, { 0, 0, -2 }, { 0, 1, -2 }, { 0, -1, -2 }, { 1, 0, -2 }, { -1, 0, -2 }, { 1, 1, -2 }, { -1, 1, -2 }, { 1, -1, -2 }, \
{ -1, -1, -2 }, { 2, 2, 0 }, { 2, 2, 1 }, { 2, 2, -1 }, { -2, 2, 0 }, { -2, 2, 1 }, { -2, 2, -1 }, { 2, -2, 0 }, { 2, -2, 1 }, { 2, -2, -1 }, { -2, -2, 0 }, { -2, -2, 1 }, { -2, -2, -1 }, { 0, 2, 2 },\
{ 1, 2, 2 }, { -1, 2, 2 }, { 0, -2, 2 }, { 1, -2, 2 }, { -1, -2, 2 }, { 0, 2, -2 }, { 1, 2, -2 }, { -1, 2, -2 }, { 0, -2, -2 }, { 1, -2, -2 }, { -1, -2, -2 }, { 2, 0, 2 }, { 2, 1, 2 }, { 2, -1, 2 }, \
{ -2, 0, 2 }, { -2, 1, 2 }, { -2, -1, 2 }, { 2, 0, -2 }, { 2, 1, -2 }, { 2, -1, -2 }, { -2, 0, -2 }, { -2, 1, -2 }, { -2, -1, -2 }, { 2, 2, 2 }, { -2, 2, 2 }, { 2, -2, 2 }, { 2, 2, -2 }, { -2, -2, 2 }, { -2, 2, -2 }, \
{ 2, -2, -2 }, { -2, -2, -2 } } }; 



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

//==============================================================================================
// functions to perform arithmetic with vectors and arrays 
//==============================================================================================

template <typename T> 
T add (T a, T b){
	
	T c;
	size_t s = a.size(); 
	for (size_t i{0}; i < s; ++i) {
		c[i] = a[i] + b[i];
	}

	return c;

}

template <typename T> 
T add (T* a, T* b){
	
	T c;
	size_t s = (*a).size(); 
	for (size_t i{0}; i < s; ++i) {
		c[i] = (*a)[i] + (*b)[i];
	}

	return c;

}


template <typename T> 
T subtract (T a, T b){
	
	T c;
	size_t s = (a).size(); 
	for (size_t i{0}; i < s; ++i) {
		c[i] = (a)[i] - (b)[i];
	}

	return c;

}

template <typename T> 
T subtract (T* a, T* b){
	
	T c;
	size_t s = (*a).size(); 
	for (size_t i{0}; i < s; ++i) {
		c[i] = (*a)[i] - (*b)[i];
	}

	return c;

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
		arr_n [i] = static_cast<double>((*array)[i])*scalar; 
	}
	return arr_n;
	
}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
// functions to perform geometric operations on arrays 
//==============================================================================================

double distance_between_points (std::array <double,3>* a1, std::array <double,3>* a2, int xlen, int ylen, int zlen){

	std::array <double,3> delta = subtract (a1, a2); 
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

double take_dot_product (std::array <double,3> o1, std::array <double,3> o2){

	// take the dot product 
	double dot_prod = 0;	

	for (int i{0}; i<3; ++i){
		dot_prod += o1[i]*o2[i];
	}

	return dot_prod; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
// function to check for avoidance in the walk 
//==============================================================================================

bool check_avoidance (const std::vector <int> to_check, const std::vector<std::vector <int>> loc_list){
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

void sarw (std::vector<std::vector<int>>* loc_list, int dop){
	
	if (dop == 0){
		return; // once the degree of polymerization hits zero, you are done
	}; 

	// until then 
	// increment final vector in loc_list in a unit direction 
	std::vector <int> next(3,0); 
	for (auto v: drns){
		 
		next = add (&((*loc_list).at((*loc_list).size()-1)), &v);
		 
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

void sarw_pbc (std::vector<std::vector<int>>* loc_list, int dop, int x_len, int y_len, int z_len){
	
	if (dop == 0){
		return; // once the degree of polymerization hits zero, you are done
	}; 

	// until then 
	// increment final vector in loc_list in a unit direction 
	std::vector <int> next(3,0); 
	for (auto v: drns){

		next = add (&((*loc_list).at((*loc_list).size()-1)), &v);
		
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

std::vector <std::vector <int>> obtain_ne_list(std::vector <int> loc, int x_len, int y_len, int z_len) {

	std::vector <std::vector <int>> nl;
	nl.reserve(6); 
	std::vector <int> v;
	for (std::vector <int> d: drns){
		v = add (&loc, &d); 
		impose_pbc(&v, x_len, y_len, z_len); 
		nl.push_back(v);
	}

	std::sort(nl.begin(), nl.end() );
	nl.erase(std::unique (nl.begin(), nl.end() ), nl.end() );  

	return nl; 

}

std::array <std::array <int,3>, 26> obtain_ne_list(std::array <int,3> loc, int x_len, int y_len, int z_len) {

	std::array <std::array <int,3>,26> nl;
	int i {0}; 
	std::array <int,3>  a = {0,0,0};
	for (std::array <int,3> d: adrns) {
		a = add (&loc, &d); 
		impose_pbc (&a, x_len, y_len, z_len); 
		nl[i] = a;
		++i; 
	}

	return nl; 

}


std::array <std::array <int,3>,98> obtain_next_ne_list (std::array <int,3> loc, int x_len, int y_len, int z_len) {

	std::array <int,3> a = {0,0,0};
	std::array <std::array <int,3>,98> next_nl; 
	int i{0};

	for ( std::array <int,3>& d1: nadrns ) {
		a = add (&loc, &d1);
		impose_pbc (&a, x_len, y_len, z_len);
		next_nl[i] = a;
		++i; 
	}

	if (i != 98) {
		std::cout << "Problem with next ne list." << std::endl;
		std::cout << "loc = "; print (loc);
		std::cout << "next neighbor list: " << std::endl;
		print_nested (next_nl);
		exit (EXIT_FAILURE);
	}

	
	return next_nl;

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


std::array <double,6> ExtractTopologyFromFile(std::string filename){
    
    std::array <double, 6> info_vec; 
    double info; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(filename); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), Emm ("Emm_1"),\
    Ems ("Emm_2"), eof ("END OF FILE"); 
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

    	else if (std::regex_search (s, Emm)){
    		double info = NumberExtractor(s); 
    		info_vec[4] = info; 
    		continue; 
    	}

    	else if (std::regex_search (s, Ems)){
    		double info = NumberExtractor(s); 
    		info_vec[5] = info; 
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

void SetUpLatticeFromScratch (std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, std::string positions, int x, int y, int z){

	(*Polymers) = ExtractPolymersFromFile (positions, x, y, z); 

	AddVoid (LATTICE, x, y, z);

	// populate the lattice 
	for (Polymer& pmer: (*Polymers)) {
		for (Particle*& p: pmer.chain){
			// now that I have my polymer coordinates, time to create the grand lattice 
			(*LATTICE).at(lattice_index (p->coords, y, z) ) = p; 
		}
	}
	
	return;
}


void SetUpLatticeFromRestart ( std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, int* step_number, std::string lattice_file_read, std::string dfile, std::string positions, int x, int y, int z ){
	
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

    return; 

}


void BiasTheStart ( std::vector <Polymer>* Polymers) {


    for ( Polymer& pmer: (*Polymers) ) {
        for ( Particle*& p: pmer.chain ) {
            p->orientation = 0;
        }
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

void AddVoid (std::vector <Particle*>* LATTICE, int x, int y, int z) {

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

				Particle* p_ptr = new Particle ( loc, "void", 0 ); 

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


void AddCosolvent (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, bool solvation, double frac, int Nmonomer, int x, int y, int z) {

	int nsol2         = std::floor ((x*y*z-Nmonomer)*frac); 
	std::cout << "Number of particles of cosolvent is " << nsol2 << "." << std::endl;
	std::vector <int> indices (x*y*z);
	std::iota (indices.begin(), indices.end(), 0); 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// get indices to loop through 
	std::shuffle ( indices.begin(), indices.end(), std::default_random_engine(seed) );
	int count = 0;
	int i = 0; 
	// find solvation shell indices 
	if ( solvation ) {
		std::set <int> solvation_shell_set;
		std::array <std::array<int,3>, 26> ne_list; 
		// get the first solvation shell 
		// auto start = std::chrono::high_resolution_clock::now(); 
		for ( Polymer& pmer: *Polymers){
			for ( Particle*& p: pmer.chain ){
				ne_list = obtain_ne_list ( p->coords, x, y, z); 
				for ( std::array <int,3>& loc: ne_list ){
					if ( (*LATTICE).at(lattice_index (loc, y, z))->ptype[0] == 's' ){
						// std::cout << "loc = "; print (loc); 
						solvation_shell_set.insert (lattice_index (loc, y, z)); 
					}
				}
			}
		}
		std::vector <int> solvation_shell_indices (solvation_shell_set.begin(), solvation_shell_set.end());
		// std::cout << "#ss = " << solvation_shell_indices.size() << std::endl;
		while ( (count < nsol2) && (count < static_cast <int> (solvation_shell_indices.size()) ) ) {
			// std::cout << "solvation_shell_indices[" << i << "] = " << solvation_shell_indices[i] << std::endl;
			// std::cout << "loc = "; print ( location(solvation_shell_indices[i], x, y, z) ); 
			// std::cout << "ptype = " << (*LATTICE).at( solvation_shell_indices.at(i) )->ptype << std::endl; 
			(*LATTICE).at( solvation_shell_indices.at(i) )->ptype = "s2"; 
			// std::cout << "ptype = " << (*LATTICE).at( solvation_shell_indices.at(i) )->ptype << std::endl; 
			(*Cosolvent).push_back ( (*LATTICE).at( solvation_shell_indices.at(i) ) );
			count += 1;
			i     += 1;
		}
	}
	// std::cout << "count = " << count << std::endl;
	i = 0; 
	while ( count < nsol2 ) {

		if ( (*LATTICE).at( indices[i] )->ptype[0] == 'm' ){
			; 
		}
		else if ( (*LATTICE).at(indices[i])->ptype == "s2" ) {
			;
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
				std::cout << "Location = "; print ((*LATTICE)[i]->coords);
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

        		if ( (*LATTICE)[ lattice_index(p->coords, y, z)]->ptype[0] == 'v' ){
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
            
            connection = subtract(&(pmer.chain[i]->coords), &(pmer.chain[i-1]->coords));
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


void CheckStructures (std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE,  int x, int y, int z){

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

			diff = subtract ( &(p->coords), to_check ); 
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

double CalculateEnergy (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array<double,nedim>* E, std::array<double,ncdim>* contacts, int x, int y, int z) {
    
    double Energy {0.0};
    (*contacts) = {0,0}; 
    std::array <std::array <int,3>, 26> ne_list; 
    std::array <std::array <int,3>, 98> next_ne_list;
    // run energy computations for every monomer bead 
    // m-m    = isotropic interaction
    // m-void = isotropic interaction 
    
    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){
        	// first neighbor list 
            ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list 
            for ( std::array <int, 3>& loc: ne_list){ 
            	switch ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype[0] ) {
            		case 'm':
            			Energy += 0.5 * (*E)[0];
            			(*contacts)[0] += 0.5; 
            			break;

            		case 'v':
            			(*contacts)[2] += 1;
            			break;
            	}
            }

            // second neighbor list 
            //std::cout << "Inside neighbor energy..." << std::endl;
            next_ne_list = obtain_next_ne_list (p->coords, x, y, z) ; 
            for ( std::array <int,3>& loc: next_ne_list) {
            	switch ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype[0] ) {
            		case 'm':
            			Energy += 0.5 * (*E)[1];
            			(*contacts)[1] += 0.5; 
            			break;

            		case 'v':
            			break;
            	}            	
            }
        }
    }
    
    return Energy; 
}

//==============================================================================================
//==============================================================================================


double NeighborEnergy ( std::vector <Particle*>* LATTICE, std::array <double,nedim>* E, std::array <double,ncdim>* contacts, int ss_index, int x, int y, int z){

	// reinitialize contacts 
	(*contacts) = {0,0,0}; 

	// information of particle whose neighbor list is being calculated 
	std::array <std::array<int,3>, 98> next_ne_list = obtain_next_ne_list ( (*LATTICE)[ss_index]->coords, x, y, z );
	std::array <std::array<int,3>, 26> ne_list      = obtain_ne_list ( (*LATTICE)[ss_index]->coords, x, y, z ); 
	std::string particle_type                       = (*LATTICE).at(ss_index)->ptype; 
	
	double Ei       = 0; 

	// std::cout << "Coords of particle of interest = "; print ( ploc );
	// start of neighbor loop 

	switch (particle_type[0]) {

		case 'm':
			for ( std::array <int,3>& loc: ne_list ) {

				// if neighbor is m1
				switch ((*LATTICE)[lattice_index(loc, y, z)]->ptype[0]) {

					case 'm':
						Ei += (*E)[0];
						(*contacts)[0] += 1;
						break;

					case 'v':
						(*contacts)[2] += 1; 
						break; 
				} 
				
			} // end of neighbor loop 

			for ( std::array <int,3>& loc: next_ne_list ) {
				// if neighbor is m1
				switch ((*LATTICE)[lattice_index(loc, y, z)]->ptype[0]) {

					case 'm':
						Ei += (*E)[1];
						(*contacts)[1] += 1;
						break;

					case 'v':
						break;

				} 
				
			} // end of next neighbor loop			



			break; 

		case 'v':
			for ( std::array <int,3>& loc: ne_list) {
				// if neighbor is m1
				switch ((*LATTICE)[lattice_index(loc, y, z)]->ptype[0]) {
					case 'm':
						(*contacts)[2] += 1;
						break;

					case 'v': 
						break; 
				} 
			} // end of neighbor loop
			// dont worry about next neighbor loop because m-v at next neighbor is not a real thing

			break;

	}

	return Ei; 

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
		if (p->ptype[0] == 'v'){
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


void dumpEnergy (double sysEnergy, int step, std::array<double,ncdim>* contacts, std::string filename) {
    std::ofstream dump_file(filename, std::ios::app); 
    // std::ostringstream os; 

    dump_file << sysEnergy << " | " << (*contacts)[0] << " | " << (*contacts)[1] << " | " << (*contacts)[2] << " | " << step << "\n";

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

void dumpMoveStatistics (std::array <int,3>* attempts, std::array <int,3>* acceptances, int step, std::string stats_file) {
    
    std::ofstream dump_file (stats_file, std::ios::out); 
    dump_file << "At step " << step << "...\n";
    dump_file << "End rotations without bias         - attempts: " << (*attempts)[0] <<", acceptances: " << (*acceptances)[0] << ", acceptance fraction: " << static_cast<double>((*acceptances)[0])/static_cast<double>((*attempts)[0]) << std::endl; 
    dump_file << "Reptation without bias             - attempts: " << (*attempts)[1] <<", acceptances: " << (*acceptances)[1] << ", acceptance fraction: " << static_cast<double>((*acceptances)[1])/static_cast<double>((*attempts)[1]) << std::endl; 
    dump_file << "Chain regrowth with overlap bias   - attempts: " << (*attempts)[2] <<", acceptances: " << (*acceptances)[2] << ", acceptance fraction: " << static_cast<double>((*acceptances)[2])/static_cast<double>((*attempts)[2]) << std::endl; 

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
                
                if ( (*LATTICE)[ lattice_index(ne, y, z) ]->ptype[0] == 'v' && std::find( solvent_indices.begin(), solvent_indices.end(), lattice_index(ne, y, z)) == solvent_indices.end()  ){
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

	std::ofstream dump_file ( filename, std::ios::app ); 
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

void TailRotation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array<double,nedim>* E, std::array<double,ncdim>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

	
    // get the neighborlist of particle at index 1 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[0]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[1]->coords;

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

    std::array <double,ncdim> cm, cm_n;
    std::array <double,ncdim> cs, cs_n;
    std::array <double,ncdim> final_contacts; 
    std::array <double,ncdim> c_contacts = *contacts;

	int choice = rng_uniform(0, 25); 

	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	double Es = 0, Es_n = 0;
	double Em = 0, Em_n = 0; 
	double Ef = 0; 

	switch ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] ){

		case 'v':
			// find the energetic interaction for the solvent molecule 
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index (loc_0, y, z), x, y, z);  
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index (ne_list[choice], y, z), x, y, z);  

			(*Polymers)[index].chain[0]->coords = ne_list[choice]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 
			
			// do the switch 
			// take the pointer of the solvent, and put it where the monomer was on the lattice 
			(*LATTICE)[ lattice_index (loc_0, y, z) ]           = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[0]; 

			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);	
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index (loc_0, y, z), x, y, z);	
			break;	

		case 'm':
			*IMP_BOOL = false;
			return;
			break; 
	
	}

	Ef = *sysEnergy - (Em+Es) + (Em_n+Es_n); 
	final_contacts = add(subtract(*contacts, add(cs, cm)), add(cs_n, cm_n)); 

	// DELETE THIS LATER 
	double energy_n = CalculateEnergy (Polymers, LATTICE, E, &c_contacts, x, y, z); 
	if (Ef != energy_n || final_contacts != c_contacts) {
		std::cout << "Either energy or contacts is messed up in end rotation... " << std::endl;
		std::cout << "Ef = " << Ef << ", energy_n = " << energy_n << ". " << std::endl; 
		std::cout << "final_contacts = "; print (final_contacts, ", "); std::cout << "c_contacts = "; print(c_contacts);
		exit(EXIT_FAILURE);
	} 

	// DELETE ABOVE LATER 
	if ( MetropolisAcceptance (*sysEnergy, Ef, temperature)){
		*sysEnergy        = Ef; 
		*contacts         = final_contacts; 
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

//==============================================================================================
//==============================================================================================

void TailRotation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array<double,nedim>* E, std::array<double,ncdim>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

	
    // get the neighborlist of particle at index 1 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[0]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[1]->coords;

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

    std::array <double,ncdim> cs, cs_n; 
    std::array <double,ncdim> cm, cm_n; 
    std::array <double,ncdim> final_contacts; 

	int choice = rng_uniform(0, 25); 

	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	double Es   = 0, Es_n = 0;
	double Em   = 0, Em_n = 0; 
	double Ef   = 0; 

	switch ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] ){

		case 'v':
			// find the energetic interaction for the solvent molecule 
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index (loc_0, y, z), x, y, z);  
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index (ne_list[choice], y, z), x, y, z);

			(*Polymers)[index].chain[0]->coords = ne_list[choice]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 
			
			// do the switch 
			// take the pointer of the solvent, and put it where the monomer was on the lattice 
			(*LATTICE)[ lattice_index (loc_0, y, z) ]           = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[0]; 

			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index (loc_0, y, z), x, y, z);
			break;	

		case 'm':
			*IMP_BOOL = false;
			return;
			break; 
	
	}

	Ef = *sysEnergy - (Em+Es) + (Em_n+Es_n); 
	final_contacts = add( subtract(*contacts, add(cm, cs)), add(cm_n, cs_n) );


	if ( MetropolisAcceptance (*sysEnergy, Ef, temperature)){
		*sysEnergy        = Ef; 
		*contacts         = final_contacts; 
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

void HeadRotation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array<double,nedim>* E, std::array<double,ncdim>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

    // get the neighborlist of particle at index 1 
    int dop = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2]->coords; 

    std::array <double,ncdim> c_contacts         = *contacts;
    std::array <double,ncdim> cs, cs_n;
    std::array <double,ncdim> cm, cm_n;
    std::array <double,ncdim> final_contacts; 

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

	int choice = rng_uniform(0, 25); 

	// made the change to the polymer
	// make the change on the lattice 
	double Es = 0, Es_n = 0;
	double Em = 0, Em_n = 0; 
	double Ef = 0; 

	switch ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] ){
		
		case 'v':
			// find the energetic interactions for the solvent molecule 
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index (loc_0, y, z), x, y, z);  
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index (ne_list[choice], y, z), x, y, z);  

			(*Polymers)[index].chain[dop-1]->coords = ne_list[choice]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 

			// do the switch 
			// take the pointer of the solvent, and put it where the monomer was on the lattice 
			(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[dop-1]; 

			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index (loc_0, y, z), x, y, z);
			break;

		case 'm':
			*IMP_BOOL = false; 
			return; 
			break;
	}

	Ef = *sysEnergy - (Em+Es) + (Em_n+Es_n);
	final_contacts = add( subtract(*contacts, add(cm, cs)), add(cm_n, cs_n) ); 

	// DELETE THIS LATER 
	double energy_n = CalculateEnergy (Polymers, LATTICE, E, &c_contacts, x, y, z); 
	if (Ef != energy_n || final_contacts != c_contacts) {
		std::cout << "Either energy or contacts is messed up in end rotation... " << std::endl; 
		std::cout << "Ef = " << Ef << ", energy_n = " << energy_n << ". " << std::endl; 
		std::cout << "final_contacts = "; print (final_contacts, ", "); std::cout << "c_contacts = "; print(c_contacts); 
		exit (EXIT_FAILURE); 
	}
	// DELETE ABOVE LATER 

	if ( MetropolisAcceptance (*sysEnergy, Ef, temperature) ){
		// std::cout << "accepted." << std::endl;
		*sysEnergy = Ef; 
		*contacts  = final_contacts; 
	}

	else {
		// std::cout << "rejected." << std::endl;
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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void HeadRotation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array<double,nedim>* E, std::array<double,ncdim>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

    // get the neighborlist of particle at index 1 
    int dop = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2]->coords; 

    // these are contacts arrays which are necessary to check validity of code 
    std::array <double,ncdim> cs, cs_n;
    std::array <double,ncdim> cm, cm_n;
    std::array <double,ncdim> final_contacts; 

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 
	int choice = rng_uniform(0, 25); 

	// made the change to the polymer
	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	double Es = 0, Es_n = 0; 
	double Em = 0, Em_n = 0; 
	double Ef = 0; 

	switch ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] ) {
		
		case 'v':
		// find the energetic interactions for the solvent molecule 
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index (loc_0, y, z), x, y, z);  
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index (ne_list[choice], y, z), x, y, z);  

			(*Polymers)[index].chain[dop-1]->coords = ne_list[choice]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 

			// do the switch 
			// take the pointer of the solvent, and put it where the monomer was on the lattice 
			(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
			(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[dop-1]; 

			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index (loc_0, y, z), x, y, z);
			break;

		case 'm':
			*IMP_BOOL = false;
			return;
			break;
	}
	

	Ef = *sysEnergy - (Em+Es) + (Em_n+Es_n);
	final_contacts = add( subtract(*contacts, add(cm, cs)), add(cm_n, cs_n) );

	if ( MetropolisAcceptance (*sysEnergy, Ef, temperature) ){
		// std::cout << "accepted." << std::endl;
		*sysEnergy = Ef; 
		*contacts  = final_contacts; 
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

void EndRotation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int index, int x, int y, int z){

	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        std::cout << "Zero index rotation!" << std::endl;
        TailRotation_UNBIASED_debug (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    else {
        std::cout << "Final index rotation!" << std::endl;
        HeadRotation_UNBIASED_debug (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void EndRotation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int index, int x, int y, int z){

	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    switch (num) {
    case (0):
        // std::cout << "Zero index rotation!" << std::endl;
        TailRotation_UNBIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        break; 
    case (1):
        // std::cout << "Final index rotation!" << std::endl;
        HeadRotation_UNBIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        break; 
    }
    return;
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of EndRotation_UNBIASED
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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void forward_reptation_with_tail_biting_new (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,nedim>* E, std::array <double,ncdim>* contacts, double* frontflow_energy, int deg_poly, int index, int x, int y, int z){

	// start performing exchanges 
	// doubles for energy transfer 
	double Em1         = 0;
	double Em2         = 0; 
	double Em1_n       = 0; 
	double Em2_n       = 0; 
	double Esys        = *frontflow_energy; 

	// ARRAY INITIALIZATION BLOCK 
	std::array <double,ncdim> cm1               = {0,0,0};
	std::array <double,ncdim> cm2               = {0,0,0};
	std::array <double,ncdim> cm1_n             = {0,0,0};
	std::array <double,ncdim> cm2_n             = {0,0,0};
	std::array <double,ncdim> current_contacts  = *contacts;
	std::array <int,3>    loc_m1                = {0,0,0}; 
	std::array <int,3>    loc_m2                = {0,0,0}; 
	// ARRAY INITIALIZATION BLOCK 

	for (int i{0}; i<deg_poly-1; ++i){
		
		loc_m1 = (*Polymers)[index].chain[i]->coords;
		loc_m2 = (*Polymers)[index].chain[i+1]->coords;
		
		// find energies 
		Em1 = NeighborEnergy (LATTICE, E, &cm1, lattice_index (loc_m1, y, z), x, y, z);
		Em2 = NeighborEnergy (LATTICE, E, &cm2, lattice_index (loc_m2, y, z), x, y, z);

		(*LATTICE)[ lattice_index(loc_m1, y, z) ] = (*LATTICE)[ lattice_index (loc_m2, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ] = (*Polymers)[index].chain[i]; 
		(*LATTICE)[ lattice_index(loc_m1, y, z) ]->coords = loc_m1; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ]->coords = loc_m2; 

		Em1_n = NeighborEnergy (LATTICE, E, &cm1_n, lattice_index (loc_m1, y, z), x, y, z);
		Em2_n = NeighborEnergy (LATTICE, E, &cm2_n, lattice_index (loc_m2, y, z), x, y, z);

		current_contacts = add ( subtract (current_contacts, add (cm1, cm2)), add (cm1_n, cm2_n) ); 

		Esys = Esys - (Em1+Em2) + (Em1_n+Em2_n); 

	}

	*frontflow_energy = Esys; 
	*contacts         = current_contacts; 

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void forward_reptation_without_tail_biting_new (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,nedim>* E, std::array <double,ncdim>* contacts, \
	std::array <int,3>* to_slither, double* frontflow_energy, int deg_poly, int index, int x, int y, int z){

	// start performing exchanges 
	// doubles for energy transfer 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 

	std::array <double,ncdim> cm               = {0,0,0};
	std::array <double,ncdim> cm_n             = {0,0,0};
	std::array <double,ncdim> cs               = {0,0,0};
	std::array <double,ncdim> cs_n             = {0,0,0};
	std::array <double,ncdim> current_contacts = *contacts;
	std::array <int,3>    loc_s            = *to_slither; 
	std::array <int,3>    loc_m            = {0,0,0}; 

	for (int i{0}; i<deg_poly; ++i){
		
		loc_s = *to_slither; 
		loc_m = (*Polymers)[index].chain[deg_poly-1-i]->coords;
		
		// find energies 
		Em = NeighborEnergy (LATTICE, E, &cm, lattice_index (loc_m, y, z), x, y, z);
		Es = NeighborEnergy (LATTICE, E, &cs, lattice_index (loc_s, y, z), x, y, z);

		(*LATTICE)[ lattice_index(loc_m, y, z) ] = (*LATTICE)[ lattice_index (loc_s, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i]; 
		(*LATTICE)[ lattice_index(loc_m, y, z) ]->coords = loc_m; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ]->coords = loc_s; 

		Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index (loc_m, y, z), x, y, z);
		Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index (loc_s, y, z), x, y, z);

		current_contacts = add ( subtract (current_contacts, add (cm, cs)), add (cm_n, cs_n) ); 

		Esys = Esys - (Em+Es) + (Em_n+Es_n); 

		// the particle is now at the old position of the monomer bead 
		*to_slither = loc_m;

	}

	*frontflow_energy = Esys; 
	*contacts         = current_contacts; 

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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void backward_reptation_with_head_butting_new (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, std::array <double,nedim>* E, std::array <double,ncdim>* contacts, \
	double* frontflow_energy, int deg_poly, int index, int x, int y, int z){

	// start performing exchanges 
	// doubles for energy transfer 
	double Em1         = 0;
	double Em2         = 0; 
	double Em1_n       = 0; 
	double Em2_n       = 0; 
	double Esys        = *frontflow_energy; 

	std::array <double,ncdim> cm1                = {0,0,0};
	std::array <double,ncdim> cm2                = {0,0,0};
	std::array <double,ncdim> cm1_n              = {0,0,0};
	std::array <double,ncdim> cm2_n              = {0,0,0};
	std::array <double,ncdim> current_contacts   = *contacts;
	std::array <int,3>    loc_m1             = {0,0,0}; 
	std::array <int,3>    loc_m2             = {0,0,0}; 

	for (int i{0}; i<deg_poly-1; ++i){
		
		loc_m1 = (*Polymers)[index].chain[deg_poly-(i+1)]->coords;
		loc_m2 = (*Polymers)[index].chain[deg_poly-(i+2)]->coords;
		
		// find energies 
		Em1 = NeighborEnergy (LATTICE, E, &cm1, lattice_index (loc_m1, y, z), x, y, z);
		Em2 = NeighborEnergy (LATTICE, E, &cm2, lattice_index (loc_m2, y, z), x, y, z);

		(*LATTICE)[ lattice_index(loc_m1, y, z) ] = (*LATTICE)[ lattice_index (loc_m2, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ] = (*Polymers)[index].chain[deg_poly-(i+1)]; 
		(*LATTICE)[ lattice_index(loc_m1, y, z) ]->coords = loc_m1; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ]->coords = loc_m2; 

		Em1_n = NeighborEnergy (LATTICE, E, &cm1_n, lattice_index (loc_m1, y, z), x, y, z);
		Em2_n = NeighborEnergy (LATTICE, E, &cm2_n, lattice_index (loc_m2, y, z), x, y, z);

		current_contacts = add ( subtract (current_contacts, add (cm1, cm2)), add (cm1_n, cm2_n) ); 

		Esys = Esys - (Em1+Em2) + (Em1_n+Em2_n); 

	}
	*frontflow_energy = Esys; 
	*contacts         = current_contacts;

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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void backward_reptation_without_head_butting_new (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, std::array <int,3>* to_slither, double* frontflow_energy, int deg_poly, int index, int x, int y, int z) {

	// start performing exchanges 
	// doubles for energy transfer 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 

	// ARRAY INITIALIZATION BLOCK

	std::array <double,ncdim> cm               = {0,0};
	std::array <double,ncdim> cs               = {0,0};
	std::array <double,ncdim> cm_n             = {0,0};
	std::array <double,ncdim> cs_n             = {0,0};
	std::array <double,ncdim> current_contacts = *contacts;
	std::array <int,3>    loc_s            = *to_slither; 
	std::array <int,3>    loc_m            = {0,0,0}; 

	// ARRAY INITIALIZATION BLOCK 

	for (int i{0}; i<deg_poly; ++i){
		
		loc_s = *to_slither; 
		loc_m = (*Polymers)[index].chain[i]->coords;
		
		// find energies 
		Es = NeighborEnergy (LATTICE, E, &cs, lattice_index (loc_s, y, z), x, y, z);
		Em = NeighborEnergy (LATTICE, E, &cm, lattice_index (loc_m, y, z), x, y, z);

		(*LATTICE)[ lattice_index(loc_m, y, z) ] = (*LATTICE)[ lattice_index (loc_s, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ] = (*Polymers)[index].chain[i]; 
		(*LATTICE)[ lattice_index(loc_m, y, z) ]->coords = loc_m; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ]->coords = loc_s; 

		Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index (loc_s, y, z), x, y, z);
		Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index (loc_m, y, z), x, y, z);

		current_contacts = add ( subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 

		Esys = Esys - (Es+Em) + (Es_n+Em_n); 

		// the particle is now at the old position of the monomer bead 
		*to_slither = loc_m;

	}

	*frontflow_energy = Esys; 
	*contacts         = current_contacts; 

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

void ForwardReptation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	std::cout << "In forward reptation unbiased." <<std::endl;

	int                   deg_poly             = (*Polymers)[index].deg_poly; 
	std::array <int,3>    loc0                 = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>    locf                 = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,ncdim> current_contacts = *contacts; 
	std::array <double,ncdim> copy_contacts    = *contacts; 
	double                Esys                 = *sysEnergy; 


	// first check if tail rotation can be performed at all 	
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( locf, x, y, z ); 
    int choice = rng_uniform (0, 25);	

	std::array <int,3> to_slither      = ne_list [choice]; 
	std::array <int,3> to_slither_copy = to_slither; 
	Particle* tmp {nullptr}; 
	

	std::cout << "Suggested location is "; print (ne_list[choice]);
	std::cout << "to_slither = "; print (to_slither);
	if ( to_slither == loc0 ){
		std::cout << "With tail biting..." << std::endl;
		// if you are performing tail biting, you need neighborhood information about all monomer beads 
		forward_reptation_with_tail_biting_new    (Polymers, LATTICE, E, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 'v' ){
		// IMPORTANT STEP 
		std::cout << "Without tail biting... " << std::endl;
		forward_reptation_without_tail_biting_new (Polymers, LATTICE, E, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
	}

	else {
		std::cout << "Rejected due to collision." << std::endl;
		*IMP_BOOL = false; 
		return; 
	}
	(*Polymers)[0].printChainCoords(); 
	std::cout << "First check structures..." << std::endl;
	CheckStructures (Polymers, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl; 
	// calculate energy of current state 
	// delete later 

	double energy_n = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z);
	if ( Esys != energy_n || current_contacts != copy_contacts){
		std::cout << "Either energy or contacts is messed up in forward reptation... " << std::endl; 
		std::cout << "Esys = " << Esys << ", energy_n = " << energy_n << ". " << std::endl; 
		std::cout << "current_contacts = "; print (current_contacts, ", "); std::cout << "copy_contacts = "; print(copy_contacts); 
		exit (EXIT_FAILURE); 
	}
	// delete above 

	if ( MetropolisAcceptance (*sysEnergy, energy_n, temperature) ){
		std::cout << "Accepted!" << std::endl;
		*sysEnergy = Esys;
		*contacts  = current_contacts; 

	}
	else {
		// revert back to old state 
		std::cout << "Energy = " << Esys << std::endl;
		std::cout << "Rejected..." << std::endl;
		if ( to_slither_copy == loc0 ) {
			std::cout << "with head_butting... " << std::endl;
			backward_reptation_with_head_butting (Polymers, LATTICE, &to_slither_copy, deg_poly, index, y, z); 
		}
		else {
			// IMPORTANT STEP 
			std::cout << "without head_butting... " << std::endl;
			tmp = (*LATTICE)[ lattice_index (loc0, y, z) ]; 
			backward_reptation_without_head_butting (Polymers, LATTICE, tmp, &loc0, &to_slither_copy, deg_poly, index, y, z); 

		}
		*IMP_BOOL = false; 
	}

	std::cout << "Second check structures..." << std::endl;
	CheckStructures (Polymers, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl; 

	return; 
}

//==============================================================================================
//==============================================================================================

void ForwardReptation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	// std::cout << "In forward reptation unbiased." <<std::endl;

	int                       deg_poly         = (*Polymers)[index].deg_poly; 
	std::array <int,3>        loc0             = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>        locf             = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,ncdim> current_contacts = *contacts; 
	double                    Esys             = *sysEnergy; 


	// first check if tail rotation can be performed at all 	
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( locf, x, y, z ); 
    int choice = rng_uniform (0, 25);	

	std::array <int,3> to_slither      = ne_list [choice]; 
	std::array <int,3> to_slither_copy = to_slither; 
	Particle* tmp {nullptr}; 
	

	// std::cout << "Suggested location is "; print (ne_list[choice]);
	// std::cout << "to_slither = "; print (to_slither);
	if ( to_slither == loc0 ){
		// std::cout << "With tail biting..." << std::endl;
		// if you are performing tail biting, you need neighborhood information about all monomer beads 
		forward_reptation_with_tail_biting_new    (Polymers, LATTICE, E, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 'v' ){
		// IMPORTANT STEP 
		// std::cout << "Without tail biting... " << std::endl;
		forward_reptation_without_tail_biting_new (Polymers, LATTICE, E, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
	}

	else {
		*IMP_BOOL = false; 
		return; 
	}

	// calculate energy of current state 
	// delete later 

	
	if ( MetropolisAcceptance (*sysEnergy, Esys, temperature) ){
		// std::cout << "Accepted!" << std::endl;
		*sysEnergy = Esys;
		*contacts  = current_contacts; 

	}
	else {
		// revert back to old state 
		// std::cout << "Rejected..." << std::endl;
		if ( to_slither_copy == loc0 ) {
			// std::cout << "with head_butting... " << std::endl;
			backward_reptation_with_head_butting (Polymers, LATTICE, &to_slither_copy, deg_poly, index, y, z); 
		}
		else {
			// IMPORTANT STEP 
			// std::cout << "without head_butting... " << std::endl;
			tmp = (*LATTICE)[ lattice_index (loc0, y, z) ]; 
			backward_reptation_without_head_butting (Polymers, LATTICE, tmp, &loc0, &to_slither_copy, deg_poly, index, y, z); 

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

void BackwardReptation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	std::cout << "In backward reptation unbiased." <<std::endl;

	int                   deg_poly         = (*Polymers)[index].deg_poly; 
	std::array <int,3>    loc0             = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>    locf             = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,ncdim> current_contacts = *contacts; 
	std::array <double,ncdim> copy_contacts    = *contacts;
	double                Esys             = *sysEnergy; 

	// first check if tail rotation can be performed at all 
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( loc0, x, y, z ); 
    int choice = rng_uniform (0, 25);	

	std::array <int,3> to_slither = ne_list[choice]; 
	std::array <int,3> to_slither_copy = to_slither; 
	Particle* tmp {nullptr}; 

	std::cout << "Suggested location is "; print (ne_list[choice]);
	std::cout << "to_slither = "; print (to_slither);

	if ( to_slither == locf ){
		std::cout << "With head butting..." << std::endl;
		backward_reptation_with_head_butting_new (Polymers, LATTICE, E, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 'v' ){
		std::cout << "Without head butting... " << std::endl;
		backward_reptation_without_head_butting_new (Polymers, LATTICE, E, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
	}

	else {
		std::cout << "Rejected from collision." << std::endl;
		*IMP_BOOL = false; 
		return; 
	}

	(*Polymers)[0].printChainCoords();

	std::cout << "First check structures..." << std::endl;
	CheckStructures (Polymers, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl;

	// delete later 
	
	double energy_n = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z);
	if ( Esys != energy_n || current_contacts != copy_contacts){
		std::cout << "Either energy or contacts is messed up in forward reptation... " << std::endl; 
		std::cout << "Esys = " << Esys << ", energy_n = " << energy_n << ". " << std::endl; 
		std::cout << "current_contacts = "; print (current_contacts, ", "); std::cout << "copy_contacts = "; print(copy_contacts); 
		exit (EXIT_FAILURE); 
	}
	
	// delete above 
	if ( MetropolisAcceptance (*sysEnergy, Esys, temperature) ){
		std::cout << "Accepted!" << std::endl;
		*sysEnergy       = Esys;
		*contacts        = current_contacts; 

	}
	else {
		// revert back to old state 
		std::cout << "Energy = " << Esys << std::endl;
		std::cout << "Rejected..." << std::endl;
		if ( to_slither_copy == locf ) {

			std::cout << "with head_butting... " << std::endl;
			forward_reptation_with_tail_biting (Polymers, LATTICE, &to_slither_copy, deg_poly, index, y, z); 
		}
		else {
			std::cout << "without head_butting... " << std::endl;
			tmp = (*LATTICE)[ lattice_index (locf, y, z) ]; 
			forward_reptation_without_tail_biting (Polymers, LATTICE, tmp, &locf, &to_slither_copy, deg_poly, index, y, z); 
		}
		*IMP_BOOL = false; 
	}

	std::cout << "Second check structures..." << std::endl;
	CheckStructures (Polymers, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl;

	return; 
}

//==============================================================================================
//==============================================================================================


void BackwardReptation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	// std::cout << "In backward reptation unbiased." <<std::endl;

	int                       deg_poly         = (*Polymers)[index].deg_poly; 
	std::array <int,3>        loc0             = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>        locf             = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,ncdim> current_contacts = *contacts; 
	double                    Esys             = *sysEnergy; 

	// first check if tail rotation can be performed at all 
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( loc0, x, y, z ); 
    int choice = rng_uniform (0, 25);	

	std::array <int,3> to_slither = ne_list[choice]; 
	std::array <int,3> to_slither_copy = to_slither; 
	Particle* tmp {nullptr}; 

	// std::cout << "Suggested location is "; print (ne_list[choice]);
	// std::cout << "to_slither = "; print (to_slither);

	if ( to_slither == locf ){
		// std::cout << "With head butting..." << std::endl;
		backward_reptation_with_head_butting_new (Polymers, LATTICE, E, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 'v' ){
		// std::cout << "Without head butting... " << std::endl;
		backward_reptation_without_head_butting_new (Polymers, LATTICE, E, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
	}

	else {
		*IMP_BOOL = false; 
		return; 
	}

	if ( MetropolisAcceptance (*sysEnergy, Esys, temperature) ){
		// std::cout << "Accepted!" << std::endl;
		*sysEnergy       = Esys;
		*contacts        = current_contacts; 

	}
	else {
		// revert back to old state 
		// std::cout << "Rejected..." << std::endl;
		if ( to_slither_copy == locf ) {
			// std::cout << "with head_butting... " << std::endl;
			forward_reptation_with_tail_biting (Polymers, LATTICE, &to_slither_copy, deg_poly, index, y, z); 
		}
		else {
			// std::cout << "without head_butting... " << std::endl;
			tmp = (*LATTICE)[ lattice_index (locf, y, z) ]; 
			forward_reptation_without_tail_biting (Polymers, LATTICE, tmp, &locf, &to_slither_copy, deg_poly, index, y, z); 
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


void Reptation_UNBIASED_debug (std::vector<Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 

    if (num==0){
        std::cout << "Backward reptation only!" << std::endl;
        BackwardReptation_UNBIASED_debug (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
        return; 
    }
    else {
        std::cout << "Forward reptation!" << std::endl;
        ForwardReptation_UNBIASED_debug  (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    
}

//==============================================================================================
//==============================================================================================

void Reptation_UNBIASED (std::vector<Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 

    switch (num){
    case (0):
        // std::cout << "Backward reptation only!" << std::endl;
        BackwardReptation_UNBIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
        break;
         
    case (1):
        // std::cout << "Forward reptation!" << std::endl;
        ForwardReptation_UNBIASED  (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        break;
    }
    
    return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of Reptation_UNBIASED
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

void ChainRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int p_index, int x, int y, int z){


	int deg_poly = (*Polymers)[p_index].deg_poly; 
	int m_index  = rng_uniform (1, deg_poly-2); 


	std::array <double,ncdim> c1_contacts     = *contacts; 
	std::array <double,ncdim> c2_contacts     = *contacts; 

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
			old_cut.push_back ((*Polymers).at(p_index).chain.at(i)->coords) ;
		}

		// regrow the polymer frontwards
		// std::cout << "Regrowth time!" << std::endl;
		HeadRegrowth_BIASED (Polymers, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		// std::cout << "Worked!" << std::endl;
		
		for (int i {m_index+1}; i<deg_poly; ++i){
			new_cut.push_back ((*Polymers).at(p_index).chain.at(i)->coords) ;
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

		c2_contacts = c1_contacts; 
		// std::cout << "c2_contacts = "; print (c2_contacts);
		backflow_energy = frontflow_energy; 

		// std::cout << "-------------------" << std::endl;
		// std::cout << "Backflow time!" << std::endl;
		BackFlowFromHeadRegrowth_BIASED (Polymers, LATTICE, &old_cut, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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
		TailRegrowth_BIASED (Polymers, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 

		for (int i {0}; i<m_index; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
		}

		if ( old_cut == new_cut ){

			return; 
		}

		
		if ( !(*IMP_BOOL) ) {
			acceptance_after_tail_regrowth ( LATTICE, &new_cut, &old_cut, y, z); 
			
			return; 
		}

		c2_contacts = c1_contacts; 
		backflow_energy = frontflow_energy; 

		// std::cout << "BEGIN BACK FLOW! " << std::endl;

		BackFlowFromTailRegrowth_BIASED (Polymers, LATTICE, &old_cut, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

		// check acceptance criterion 
		double rng_acc = rng_uniform (0.0, 1.0); 
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n ) {

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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void ChainRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int p_index, int x, int y, int z){


	int deg_poly = (*Polymers)[p_index].deg_poly; 
	std::cout << "deg_poly = " << deg_poly << std::endl;
	int m_index  = rng_uniform (deg_poly/2, deg_poly-2); 


	std::array <double,ncdim> c1_contacts     = *contacts; 
	std::array <double,ncdim> c2_contacts     = *contacts; 

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
		
		std::cout << "Head regrowth... " << std::endl;
		// std::cout << "OLD ENERGY = " << *sysEnergy << std::endl << std::endl;
		// get old_cut 
		
		for (int i {m_index+1}; i<deg_poly; ++i){
			old_cut.push_back ((*Polymers).at(p_index).chain.at(i)->coords) ;
		}

		// regrow the polymer frontwards
		// std::cout << "Regrowth time!" << std::endl;
		
		HeadRegrowth_BIASED_debug (Polymers, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		
		// std::cout << "Worked!" << std::endl;
		
		for (int i {m_index+1}; i<deg_poly; ++i){
			new_cut.push_back ((*Polymers).at(p_index).chain.at(i)->coords) ;
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

		c2_contacts = c1_contacts; 
		backflow_energy = frontflow_energy; 

		std::cout << "-------------------" << std::endl;
		std::cout << "Backflow time!" << std::endl;
		BackFlowFromHeadRegrowth_BIASED_debug (Polymers, LATTICE, &old_cut, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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

		std::cout << "Tail regrowth..." << std::endl;
		// get old cut 
		for ( int i{0}; i<m_index; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords);
		}

		// regrow the polymer backwards
		TailRegrowth_BIASED_debug (Polymers, LATTICE, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 

		for (int i {0}; i<m_index; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
		}

		if ( old_cut == new_cut ){
			// std::cout << "------------------------------------------------------------------------"  << std::endl;
			// std::cout << "POLYMER CONFIGURATION WAS NOT CHANGED. RETURN BACK TO MAIN."               << std::endl;
			// std::cout << "------------------------------------------------------------------------"  << std::endl << std::endl;
			return; 
		}

		
		if ( !(*IMP_BOOL) ) {
			acceptance_after_tail_regrowth ( LATTICE, &new_cut, &old_cut, y, z); 
			
			return; 
		}

		c2_contacts = c1_contacts; 
		backflow_energy = frontflow_energy; 

		// std::cout << "BEGIN BACK FLOW! " << std::endl;

		BackFlowFromTailRegrowth_BIASED_debug (Polymers, LATTICE, &old_cut, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

		if ( *sysEnergy != backflow_energy || c2_contacts != *contacts ){
			std::cout << "Energies are bad, or contacts are not right." << std::endl;
			std::cout << "*sysEnergy = " << *sysEnergy << ", backflow_energy = " << backflow_energy << "." << std::endl;
			std::cout << "c2_contacts = "; print (c2_contacts, ", "); std::cout << "*contacts = "; print(*contacts);
			std::cout << "Shit's fucked." << std::endl;
			exit(EXIT_FAILURE);
		}

		// check acceptance criterion 
		double rng_acc = rng_uniform (0.0, 1.0); 
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n ) {

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

void HeadRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	std::array <double,ncdim>               current_contacts = *contacts; 

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>               energies; 
	std::array <std::array<double,ncdim>,5> contacts_store; 
	std::array <double,5>               boltzmann; 
	std::array <int,3>                  loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 
	
	// doubles for energy transfer 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 

	// ARRAY INITIALIZATION BLOCK
	std::array <double,ncdim> cm   = {0,0};
	std::array <double,ncdim> cs   = {0,0};
	std::array <double,ncdim> cm_n = {0,0};
	std::array <double,ncdim> cs_n = {0,0};
	// ARRAY INITIALIZATION BLOCK

	// start attempting jumps 

	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 


	while ( idx_counter < 5 ){

		if ( ne_list[idx_counter] == loc_m ){ 
			// current position
			energies[idx_counter]       = Esys; 
			contacts_store[idx_counter] = current_contacts;
			block_counter += 1; 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

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
				contacts_store[idx_counter] = {-1,-1}; 
				block_counter += 1; 
			}
			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the new energies
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
				contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];
			}
			
		}
		else {
			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

			// get the new energies
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			// std::cout << "Es_n = " << Es_n << std::endl; 
			// std::cout << "Em_n = " << Em_n << std::endl;

			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	

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

	// do the swap again
	(*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords           = ne_list[e_idx];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];

	// else { do nothing and maintain structure, because suggested index is in the maintain index vector }
	HeadRegrowth_BIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void HeadRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z) {


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	std::array <double,ncdim>               current_contacts = *contacts; 
	std::array <double,ncdim>				copy_contacts    = *contacts; // this is only for the checks with brute force calculation

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>                   energies; 
	std::array <std::array<double,ncdim>,5> contacts_store; 
	std::array <double,5>                   boltzmann; 
	std::array <int,3>                      loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 
	
	// doubles for energy transfer 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	double Etest      = 0;

	// ARRAY INITIALIZATION BLOCK 
	std::array <double,ncdim> cm   = {0,0,0};
	std::array <double,ncdim> cs   = {0,0,0};
	std::array <double,ncdim> cm_n = {0,0,0};
	std::array <double,ncdim> cs_n = {0,0,0};
	// ARRAY INITIALIZATION BLOCK 

	// start attempting jumps 

	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 


	while ( idx_counter < 5 ){
		/*
		std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << std::endl;
		std::cout << "idx_counter = " << idx_counter << std::endl;
		std::cout << "ptype  = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;
		std::cout << "coords = "; print ((*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->coords);
		std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% " << std::endl;
		*/ 
		if ( ne_list[idx_counter] == loc_m ){ 
			// current position
			energies[idx_counter]       = Esys; 
			contacts_store[idx_counter] = current_contacts;
			block_counter += 1; 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

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
				// ARRAY INITIALIZATION BLOCK
				contacts_store[idx_counter] = {-1,-1,-1}; 
				// ARRAY INITIALIZATION BLOCK 
				block_counter += 1; 
			}
			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the new energies
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
				contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 
				
				// get the energy
				// delete later 
				
				Etest  = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
				if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				 
				// delete above later 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];
			}
			
		}
		else {
			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

			// get the new energies
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			// std::cout << "Es_n = " << Es_n << std::endl; 
			// std::cout << "Em_n = " << Em_n << std::endl;

			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	

			// get the energy
			// delete later 
			
			Etest  = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
			if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
				std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
				std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
				std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
				std::cout << "lattice index = " << lattice_index( ne_list[idx_counter], y, z ) << ", loc = "; print (ne_list[idx_counter]); 
				std::cout << "lattice index = " << lattice_index( loc_m, y, z ) << ", loc = "; print (loc_m);
				std::cout << "Shit's fucked." << std::endl;
				dumpLATTICE(LATTICE, 1, y, z, "debug.lattice");
				dumpPositionsOfPolymers (Polymers, 1, "debug.coords");
				exit(EXIT_FAILURE);
			}
			
			// delete above 

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


	// do the swap again
	(*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords           = ne_list[e_idx];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];

	// else { do nothing and maintain structure, because suggested index is in the maintain index vector }
	HeadRegrowth_BIASED_debug (Polymers, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of BackFlowFromHeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	double                                   rboltzmann       = 0; // running sum for boltzmann weights 
	std::array <double,5>                    energies         = {0,0,0,0,0}; 
	std::array <double,5>                    boltzmann        = {0,0,0,0,0};  
	std::array <double,ncdim>                current_contacts = *contacts; 
	std::array <std::array<double,ncdim>,5>  contacts_store; 

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


	// these are doubles for energies 
	double Em      = 0;
	double Es      = 0; 
	double Em_n    = 0; 
	double Es_n    = 0; 
	double Esys    = *backflow_energy; 

	// ARRAY INITIALIZATION BLOCK 
	std::array <double,ncdim> cm   = {0,0,0};
	std::array <double,ncdim> cs   = {0,0,0};
	std::array <double,ncdim> cm_n = {0,0,0};
	std::array <double,ncdim> cs_n = {0,0,0};
	// ARRAY INITIALIZATION BLOCK 

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
					break; 
				}
			}

			if ( self_swap_idx < m_index ){
				energies [idx_counter] = 1e+08; 
				contacts_store[idx_counter] = {-1,-1}; 
			}
			else {

				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap 
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
				
				// get the new energies
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em) + (Es_n+Em_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];
			}
		}

		else {

			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index+1];

			// get the energy 
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	

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
	BackFlowFromHeadRegrowth_BIASED (Polymers, LATTICE, old_cut, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::vector <std::array<int,3>>* old_cut, std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z) {


	if (m_index == deg_poly-1){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	
	std::array <double,5>                   energies; 
	std::array <double,5>                   boltzmann; 
	std::array <double,ncdim>               copy_contacts    = *contacts; 
	std::array <double,ncdim>               current_contacts = *contacts; 
	std::array <std::array<double,ncdim>,5> contacts_store; 
	double                                  rboltzmann       = 0; // running sum for boltzmann weights 

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


	// these are doubles for energies 
	double Em      = 0;
	double Es      = 0; 
	double Em_n    = 0; 
	double Es_n    = 0; 
	double Esys    = *backflow_energy; 
	double Etest   = 0;


	// ARRAY INITIALIZATION BLOCK 

	std::array <double,ncdim> cm   = {0,0,0};
	std::array <double,ncdim> cs   = {0,0,0};
	std::array <double,ncdim> cm_n = {0,0,0};
	std::array <double,ncdim> cs_n = {0,0,0};

	// ARRAY INITIALIZATION BLOCK

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
				contacts_store[idx_counter] = {-1,-1}; 
			}
			else {

				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap 
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
				
				// get the new energies
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em) + (Es_n+Em_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 
				
				// delete later 
				
				Etest  = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
				if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				
				// delete above 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]         = (*Polymers)[p_index].chain[m_index+1];
			}
		}

		else {

			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index+1];

			// get the energy 
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	
			
			// delete later
			
			Etest  = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
			if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
				std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
				std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
				std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
				std::cout << "lattice index = " << lattice_index( ne_list[idx_counter], y, z ) << ", loc = "; print (ne_list[idx_counter]); 
				std::cout << "lattice index = " << lattice_index( loc_m, y, z ) << ", loc = "; print (loc_m);
				std::cout << "Shit's fucked." << std::endl;
				// dumpLATTICE(LATTICE, 1, y, z, "debug.lattice");
				// dumpPositionsOfPolymers (Polymers, 1, "debug.coords");
				exit(EXIT_FAILURE);
			}
			
			// delete above later

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
	BackFlowFromHeadRegrowth_BIASED_debug (Polymers, LATTICE, old_cut, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of BackflowHeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of TailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void TailRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	std::array <double,ncdim>               current_contacts = *contacts; 

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>                   energies; 
	std::array <std::array<double,ncdim>,5> contacts_store; 
	std::array <double,5>                   boltzmann; 
	std::array <int,3>                      loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// doubles for energy transfer 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	
	// contacts for energy transfer 
	// ARRAY INITIALIZATION BLOCK 
	std::array <double,ncdim> cm   = {0,0,0};
	std::array <double,ncdim> cs   = {0,0,0};
	std::array <double,ncdim> cm_n = {0,0,0};
	std::array <double,ncdim> cs_n = {0,0,0};
	// ARRAY INITIALIZATION BLOCK 

	// start attempting jumps 
	int block_counter          = 0; 
	int idx_counter            = 0; 
	int self_swap_idx          = -1; 


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

			if ( self_swap_idx > m_index ){
				energies[idx_counter] = 1e+08; 
				// ARRAY INITIALIZATION BLOCK
				contacts_store[idx_counter] = {-1,-1,-1}; 
				// ARRAY INITIALIZATION BLOCK 
				block_counter += 1; 
			}
			else {
				// get the initial neighboring energies 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z); 
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the energy
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
				contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 

				// revert back to original structure 
				(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords   = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
			}

		}
		else {
			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the new energies
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			// run the computation 
			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	

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

	// if { e_idx is not in the maintain index vector, perform the swap}
		
	// do the swap again 	
	(*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords           = ne_list[e_idx];
	
	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
	

	// else { do nothing and maintain structure, because the suggested index is in the the maintain index vector }
	TailRegrowth_BIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void TailRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	std::array <double,ncdim>               current_contacts = *contacts; 
	std::array <double,ncdim>				copy_contacts    = *contacts; // this is only for the checks with brute force calculation

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>                   energies; 
	std::array <std::array<double,ncdim>,5> contacts_store; 
	std::array <double,5>                   boltzmann; 
	std::array <int,3>                      loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// doubles for energy transfer 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	double Etest      = 0;
	
	// contacts for energy transfer 
	// ARRAY INITIALIZATION BLOCK
	std::array <double,ncdim> cm   = {0,0,0};
	std::array <double,ncdim> cs   = {0,0,0};
	std::array <double,ncdim> cm_n = {0,0,0};
	std::array <double,ncdim> cs_n = {0,0,0};
	// ARRAY INITIALIZATION BLOCK 

	// start attempting jumps 
	int block_counter          = 0; 
	int idx_counter            = 0; 
	int self_swap_idx          = -1; 


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
				contacts_store[idx_counter] = {-1,-1,-1}; 
				block_counter += 1; 
			}
			else {
				// get the initial neighboring energies 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z); 
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the energy
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				energies [idx_counter] = Esys - (Es+Em) + (Es_n + Em_n); 
				contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 
				
				// get the energy 
				// delete below
				
				Etest = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
				if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				 
				// delete above 

				// revert back to original structure 
				(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords   = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
			}

		}
		else {
			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the new energies
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			// run the computation 
			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	

			// delete this later 
			
			Etest  = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
			if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
				std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
				std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
				std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
				std::cout << "lattice index = " << lattice_index( ne_list[idx_counter], y, z ) << ", loc = "; print (ne_list[idx_counter]); 
				std::cout << "lattice index = " << lattice_index( loc_m, y, z ) << ", loc = "; print (loc_m);
				std::cout << "Shit's fucked." << std::endl;
				dumpLATTICE(LATTICE, 1, y, z, "debug.lattice");
				dumpPositionsOfPolymers (Polymers, 1, "debug.coords");
				exit(EXIT_FAILURE);
			}
			
			// delete above later 

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
	TailRegrowth_BIASED_debug (Polymers, LATTICE, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of BackFlowFromTailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromTailRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::array<double,nedim>* E, std::array <double,ncdim>* contacts, \
	bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, \
	int p_index, int m_index, int recursion_depth, int x, int y, int z){

	if (m_index == 0){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	
	// generate an array for energies 
	
	std::array <double,5>               energies; 
	std::array <double,5>               boltzmann; 
	std::array <double,ncdim>               current_contacts = *contacts; 
	std::array <std::array<double,ncdim>,5> contacts_store; 
	double                              rboltzmann       = 0; // running sum for boltzmann weights  

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; // key item of interest 

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 
	
	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[m_index-1]) - ne_list.begin() ; 

	// make sure old_cut is present at position 0 
	// crit_idx = 1; 
	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	// these are doubles for energies 
	double Em   = 0;
	double Es   = 0; 
	double Em_n = 0; 
	double Es_n = 0; 
	double Esys = *backflow_energy; 

	// ARRAY INITIALIZATION BLOCK 
	std::array <double,ncdim> cm   = {0,0,0};
	std::array <double,ncdim> cs   = {0,0,0};
	std::array <double,ncdim> cm_n = {0,0,0};
	std::array <double,ncdim> cs_n = {0,0,0};
	// ARRAY INITIALIZATION BLOCK 

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
				energies[idx_counter] = 1e+08; 
				contacts_store[idx_counter] = {-1,-1,-1}; 
			}

			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 				

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the new energies	
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em) + (Es_n+Em_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 

				// revert back to original structure 
				(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
			}

		}
		else {
			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords                 = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ]                = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the energy 
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	

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

	BackFlowFromTailRegrowth_BIASED (Polymers, LATTICE, old_cut, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromTailRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::array<double,nedim>* E, std::array <double,ncdim>* contacts, \
	bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, \
	int p_index, int m_index, int recursion_depth, int x, int y, int z){

	if (m_index == 0) {
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	
	// generate an array for energies 
	
	std::array <double,5>                   energies; 
	std::array <double,5>                   boltzmann; 
	std::array <double,ncdim>               copy_contacts    = *contacts; 
	std::array <double,ncdim>               current_contacts = *contacts; 
	std::array <std::array<double,ncdim>,5>     contacts_store; 
	double                                  rboltzmann       = 0; // running sum for boltzmann weights  

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; // key item of interest 

	// generate possible locations to jump to 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 
	
	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle   (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	int ne_idx = std::find (ne_list.begin(), ne_list.end(), (*old_cut)[m_index-1]) - ne_list.begin() ; 

	// make sure old_cut is present at position 0 
	// crit_idx = 1; 
	std::array <int,3> tmp = ne_list[0]; 
	ne_list[0]             = ne_list[ne_idx];
	ne_list[ne_idx]        = tmp; 

	// these are doubles for energies 
	double Em   = 0;
	double Es   = 0; 
	double Em_n = 0; 
	double Es_n = 0; 
	double Esys = *backflow_energy; 
	double Etest = 0;

	// ARRAY INITIALIZATION BLOCK
	std::array <double,ncdim> cm   = {0,0,0};
	std::array <double,ncdim> cs   = {0,0,0};
	std::array <double,ncdim> cm_n = {0,0,0};
	std::array <double,ncdim> cs_n = {0,0,0};
	// ARRAY INITIALIZATION BLOCK 


	// start attempting jumps 
	int idx_counter         = 0 ; 
	int self_swap_idx		= -1;

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);
		if ( ne_list [idx_counter] == loc_m ){
			energies[idx_counter]               = *backflow_energy;
			contacts_store[idx_counter]         = current_contacts;
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
				energies[idx_counter] = 1e+08; 
				// ARRAY INITIALIZATION BLOCK 
				contacts_store[idx_counter] = {-1,-1,-1}; 
				// ARRAY INITIALIZATION BLOCK 
			}

			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 				

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the new energies	
				Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z); 

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em) + (Es_n+Em_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 

				// delete later 
				
				Etest  = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
				if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				
				// delete above 

				// revert back to original structure 
				(*LATTICE) [lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];
			}

		}
		else {
			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergy (LATTICE, E, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergy (LATTICE, E, &cm, lattice_index(loc_m, y, z), x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords                 = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ]                = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the energy 
			Es_n = NeighborEnergy (LATTICE, E, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergy (LATTICE, E, &cm_n, lattice_index(loc_m, y, z), x, y, z);

			energies       [ idx_counter ] = Esys - (Es+Em) + (Es_n+Em_n); 
			contacts_store [ idx_counter ] = add (subtract (current_contacts, add (cs, cm)), add (cs_n, cm_n) ); 	

			// delete later
			
			Etest  = CalculateEnergy (Polymers, LATTICE, E, &copy_contacts, x, y, z); 
			if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ){
				std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
				std::cout << "Etest = " << Etest << ", energies [" << idx_counter <<"] = " << energies[idx_counter] << "." << std::endl;
				std::cout << "contacts_store = "; print (contacts_store[idx_counter], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
				std::cout << "lattice index = " << lattice_index( ne_list[idx_counter], y, z ) << ", loc = "; print (ne_list[idx_counter]); 
				std::cout << "lattice index = " << lattice_index( loc_m, y, z ) << ", loc = "; print (loc_m);
				std::cout << "Shit's fucked." << std::endl;
				// dumpLATTICE(LATTICE, 1, y, z, "debug.lattice");
				// dumpPositionsOfPolymers (Polymers, 1, "debug.coords");
				exit(EXIT_FAILURE);
			}
			
			// delete above later

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

	BackFlowFromTailRegrowth_BIASED_debug (Polymers, LATTICE, old_cut, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

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

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of ChainRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of HeadRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of TailRegrowth_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


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
		if ( (*LATTICE)[ lattice_index(linked_list[L-1], y, z) ]->ptype[0] == 'v' ){

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

		if ( (*LATTICE)[ lattice_index(linked_list[L-1], y, z) ]->ptype[0] == 'v' ){

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



////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of acceptance_after_tail_regrowth
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


//==============================================================================================
//==============================================================================================
//
// NAME OF FUNCTION: PerturbSystem_BIASED
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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of PerturbSystem_UNBIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void PerturbSystem_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, std::array <int,3>* attempts, \
	bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, \
	int* move_number, int x, int y, int z) {

	int index = 0; // rng_uniform (0, static_cast<int>((*Polymers).size()-1) ); 
	int r     = rng_uniform (0, 2); 
	// std::array <double,4> c_contacts = *contacts; 

	switch (r) {

		case (0):
			if (v) {
				std::cout << "Performing end rotations..." << std::endl; 
			}
			EndRotation_UNBIASED  (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break; 

		case (1):
			if (v) {
				std::cout << "Performing reptation..." << std::endl; 
			}
			Reptation_UNBIASED    (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (2):
			if (v) {
				std::cout << "Performing biased chain regrowth..." << std::endl;
			}
			ChainRegrowth_BIASED (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

	}

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of PerturbSystem_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void PerturbSystem_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::array <double,nedim>* E, std::array <double,ncdim>* contacts, std::array <int,3>* attempts, \
	bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, \
	int* move_number, int x, int y, int z) {

	int index = 0; // rng_uniform (0, static_cast<int>((*Polymers).size()-1) ); 
	int r     = rng_uniform(0, 2); 
	switch (r) {

		case (0):
			if (v) {
				std::cout << "Performing end rotations..." << std::endl; 
			}
			EndRotation_UNBIASED_debug  (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break; 

		case (1):
			if (v) {
				std::cout << "Performing reptation..." << std::endl; 
			}
			Reptation_UNBIASED_debug    (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (2):
			if (v) {
				std::cout << "Performing biased chain regrowth..." << std::endl;
			}
			ChainRegrowth_BIASED_debug (Polymers, LATTICE, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

	}

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of PerturbSystem
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
