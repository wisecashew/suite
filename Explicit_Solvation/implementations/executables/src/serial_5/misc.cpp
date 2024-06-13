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
std::array <std::array <int,3>, 26> adrns = { ax, ay, az, nx, ny, nz, axay, axaz, axny, axnz, nxay, nxaz, nxny, nxnz, ayaz, aynz, nyaz, nynz, axayaz, axnyaz, axaynz, axnynz, nxayaz, nxaynz, nxnyaz, nxnynz }; 

std::map <int, std::array<double,3>> Or2Dir = { {0, {1.0,0,0}}, {1, {0,1.0,0}}, {2, {0,0,1}}, {3, {-1,0,0}}, {4, {0,-1,0}}, {5, {0,0,-1}}, {6, {1.0/(std::sqrt(2)), 1.0/(std::sqrt(2)), 0}}, {7, {1.0/(std::sqrt(2)), 0, 1.0/(std::sqrt(2))}}, {8, {1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}}, {9, {1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}}, {10, {-1.0/(std::sqrt(2)),1.0/(std::sqrt(2)),0}}, {11, {-1.0/(std::sqrt(2)),0,1.0/(std::sqrt(2))}}, {12, {-1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}}, {13, {-1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}}, {14, {0,1.0/(std::sqrt(2)),1.0/(std::sqrt(2))}}, {15, {0,1.0/(std::sqrt(2)),-1.0/(std::sqrt(2))}}, {16, {0,-1.0/(std::sqrt(2)), 1.0/(std::sqrt(2))}}, {17, {0,-1.0/(std::sqrt(2)), -1.0/(std::sqrt(2))}}, {18, {1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {19, {1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {20, {1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {21, {1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {22, {-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {23, {-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {24, {-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {25, {-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}} };
std::map <std::array<double,3>, int> Dir2Or = { {{1.0,0,0}, 0}, {{0,1.0,0}, 1}, {{0,0,1}, 2}, {{-1,0,0}, 3}, {{0,-1,0}, 4}, {{0,0,-1}, 5}, {{1.0/(std::sqrt(2)), 1.0/(std::sqrt(2)), 0}, 6}, {{1.0/(std::sqrt(2)), 0, 1.0/(std::sqrt(2))}, 7}, {{1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}, 8}, {{1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}, 9}, {{-1.0/(std::sqrt(2)),1.0/(std::sqrt(2)),0}, 10}, {{-1.0/(std::sqrt(2)),0,1.0/(std::sqrt(2))}, 11}, {{-1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}, 12}, {{-1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}, 13}, {{0,1.0/(std::sqrt(2)),1.0/(std::sqrt(2))}, 14}, {{0,1.0/(std::sqrt(2)),-1.0/(std::sqrt(2))}, 15}, {{0,-1.0/(std::sqrt(2)), 1.0/(std::sqrt(2))}, 16}, {{0,-1.0/(std::sqrt(2)), -1.0/(std::sqrt(2))}, 17}, {{1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 18}, {{1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 19}, {{1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 20}, {{1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 21}, {{-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 22}, {{-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 23}, {{-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}, 24}, {{-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}, 25} };

const char* ws = " \t\n\r\f\v";

//==============================================================================================
// set of functions to trim strings 
//==============================================================================================

// trim from end of string (right)
inline std::string& rtrim(std::string& s, const char* t = ws)
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from beginning of string (left)
inline std::string& ltrim(std::string& s, const char* t = ws)
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from both ends of string (right then left)
inline std::string& trim(std::string& s, const char* t = ws)
{
    return ltrim(rtrim(s, t), t);
}


std::vector<std::string> split (const std::string &s, char delim) {
	std::vector<std::string> result;
	std::stringstream ss (s);
	std::string item;

	while (getline (ss, item, delim)) {
		result.push_back (item);
	}

	return result;
}



//==============================================================================================
// impose periodic boundary conditions on vector 
//==============================================================================================

void impose_pbc(std::vector <int>* vect, int x_len, int y_len, int z_len){
	for (int i{0}; i<3; i++) {
		switch (i) {
		case (0):
			(*vect)[i] = ((((*vect)[i])%x_len)+x_len)%x_len; 
			break;

		case (1):
			(*vect)[i] = ((((*vect)[i])%y_len)+y_len)%y_len; 
			break;

		case (2):
			(*vect)[i] = ((((*vect)[i])%z_len)+z_len)%z_len; 
			break;
		}
	}
	return; 
}

void impose_pbc(std::array <int,3>* arr, int x_len, int y_len, int z_len) {
	for (int i{0}; i<3; ++i){
		switch (i) {
			case (0):
				(*arr)[i] = (((*arr)[i]%x_len)+x_len)%x_len; 
				break;

			case (1):
				(*arr)[i] = (((*arr)[i]%y_len)+y_len)%y_len;
				break; 

			case (2):
				(*arr)[i] = (((*arr)[i]%z_len)+z_len)%z_len; 
				break; 
		}
	}

	return; 
}


void impose_pbc(std::array <double,3>* arr, int x_len, int y_len, int z_len) {
	for (int i{0}; i<3; ++i){
		switch (i) {
			case (0):
				(*arr)[i] = std::fmod( ( std::fmod( (*arr)[i], x_len) +x_len), x_len); 
				break;
			case (1):
				(*arr)[i] = std::fmod( ( std::fmod( (*arr)[i], y_len) +y_len), y_len); 
				break; 
			
			case (2):
				(*arr)[i] = std::fmod( ( std::fmod( (*arr)[i], z_len) +z_len), z_len);
				break; 
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

std::array <double,8> add_arrays(std::array <double,8>* a1, std::array <double,8>* a2){
	std::array<double, 8> a3; 

	for (int i{0}; i<8; ++i){
		a3[i] = (*a1)[i] + (*a2)[i]; 
	}

	return a3;

}

std::array <double,8> add_arrays (std::array<double,8> a1, std::array <double,8> a2) {
	
	std::array<double, 8> a3; 
	for (int i{0}; i<8; ++i){
		a3[i] = (a1)[i] + (a2)[i]; 
	}

	return a3;

}

std::array <double,8> subtract_arrays (std::array<double,8> a1, std::array <double,8> a2) {
	
	std::array<double, 8> a3; 
	for (int i{0}; i<8; ++i){
		a3[i] = (a1)[i] - (a2)[i]; 
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

std::array <int,8> subtract_arrays(std::array <int,8>* a1, std::array <int,8>* a2){
	std::array<int, 8> a3; 
	for (int i{0}; i<8; ++i){
		a3[i] = (*a1)[i] - (*a2)[i]; 
	}
	return a3;
}

std::array <double,8> subtract_arrays(std::array <double,8>* a1, std::array <double,8>* a2){
	std::array<double, 8> a3; 

	for (int i{0}; i<8; ++i){
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

	std::array <double,3> delta = subtract_arrays (a1, a2); 
	// std::cout << "delta array is "; print (delta);
	for ( int i{0}; i<3; ++i){
		switch (i) {
		case (0):
			delta[i] = modified_modulo ( delta[i], xlen );
			break; 
		case (1):
			delta[i] = modified_modulo ( delta[i], ylen ); 
			break;
		case (2):
			delta[i] = modified_modulo ( delta[i], zlen ); 
			break; 
		}
	}

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


std::array <double,13> ExtractTopologyFromFile( std::string filename, std::map <std::pair <std::string, std::string>, std::tuple<std::string, double, double, int, int> >* InteractionMap ){
    
    std::array <double, 13> info_vec; 
    double info; 
    std::array <int,17> input_hit_check = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::pair  <std::string, std::string> p; 
    std::tuple <std::string, double, double, int, int> t; 
    std::string interaction_type; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(filename); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), Emm_a ("Emm_a"),\
    Emm_n ("Emm_n"), Ems1_a ("Ems1_a"), Ems1_n ("Ems1_n"), Ems2_a ("Ems2_a"), \
    Ems2_n ("Ems2_n"), Es1s2_a("Es1s2_a"), Es1s2_n ("Es1s2_n"), \
    m1_m1 ("m1-m1"), m1_s1 ("m1-s1"), m1_s2 ("m1-s2"), s1_s2 ("s1-s2"),\
    frac ("frac"), eof ("END OF FILE"); 

    for (std::string s: contents){

    	if (std::regex_search(s, x)){
    		info = NumberExtractor(s); 
    		info_vec[0]=info; 
    		input_hit_check[0] += 1;
    		if (input_hit_check [0] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue; 
    	}

    	else if (std::regex_search(s, y)){
    		info = NumberExtractor(s); 
    		info_vec[1] = info; 
    		input_hit_check[1] += 1;
    		if (input_hit_check [1] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue; 
    	}

    	else if (std::regex_search(s, z)){
    		info = NumberExtractor(s); 
    		info_vec[2] = info; 
    		input_hit_check[2] += 1;
    		if (input_hit_check [2] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue; 
    	}

		else if (std::regex_search(s, kT)){
    		info = NumberExtractor(s); 
    		info_vec[3] = info ; 
    		input_hit_check[3] += 1;
    		if (input_hit_check [3] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue; 
    	}

    	else if (std::regex_search(s, m1_m1)){
    		interaction_type = trim (split(s, ':')[1]);
    		if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
    			p = std::make_pair ("m1", "m1");
    			t = std::make_tuple(interaction_type, 0, 0, 0, 1);
    			(*InteractionMap)[p] = t;
    		}
    		else {
    			std::cout << "\"" << s << "\"" << std::endl;
    			std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
    			exit(EXIT_FAILURE);     			
    		}
    		input_hit_check[4] += 1;
    		if (input_hit_check [4] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}

    	else if (std::regex_search(s, m1_s1)){

    		interaction_type = trim (split(s, ':')[1]);
    		if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
    			p = std::make_pair ("m1", "s1");
    			(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 2, 3);
    			p = std::make_pair ("s1", "m1") ;
    			(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 2, 3);
    		}
    		else {
    			std::cout << "\"" << s << "\"" << std::endl;
    			std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
    			exit(EXIT_FAILURE); 
    		}
    		input_hit_check[5] += 1;
    		if (input_hit_check [5] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}

    	else if (std::regex_search(s, m1_s2)){
    		interaction_type = trim (split(s, ':')[1]);
    		if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
	    		p = std::make_pair ("m1", "s2");
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 4, 5);
	    		p = std::make_pair ("s2", "m1");
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 4, 5);	
	    	}
	    	else {
    			std::cout << "\"" << s << "\"" << std::endl;
    			std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
    			exit(EXIT_FAILURE); 	    		
	    	}
	    	input_hit_check[6] += 1;
	    	if (input_hit_check [6] > 1){
	    		std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
	    	continue;
    	}

    	else if (std::regex_search(s, s1_s2)){
    		interaction_type = trim (split(s, ':')[1]);
    		if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
	    		p = std::make_pair ("s1", "s2");
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 6, 7);
	    		p = std::make_pair ("s2", "s1");
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 6, 7);	
    		}
    		else {
    			std::cout << "\"" << s << "\"" << std::endl;
    			std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
    			exit(EXIT_FAILURE); 	    
    		}
    		input_hit_check[7] += 1;
    		if (input_hit_check [7] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}

    	else if (std::regex_search (s, Emm_a)){
    		info = NumberExtractor(s); 
    		info_vec[4] = info; 
    		input_hit_check[16] += 1;
    		if (input_hit_check [16] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue; 
    	}

    	else if (std::regex_search (s, Emm_n)){
    		info = NumberExtractor(s); 
    		info_vec[5] = info; 
    		input_hit_check[8] += 1;
    		if (input_hit_check [8] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue; 
    	}

    	else if (std::regex_search (s, Ems1_a)){
    		info = NumberExtractor(s); 
    		info_vec[6] = info; 
    		input_hit_check[9] += 1;
    		if (input_hit_check [9] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue; 
    	}

    	else if (std::regex_search (s, Ems1_n)){

    		info = NumberExtractor(s);
    		info_vec[7] = info; 
    		input_hit_check[10] += 1;
    		if (input_hit_check [10] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}
    	else if (std::regex_search (s, Ems2_a)){

    		info = NumberExtractor(s);
    		info_vec[8] = info; 
    		input_hit_check[11] += 1;
    		if (input_hit_check [11] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}
    	else if (std::regex_search (s, Ems2_n)){

    		info = NumberExtractor(s);
    		info_vec[9] = info; 
    		input_hit_check[12] += 1;
    		if (input_hit_check [12] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}
    	else if (std::regex_search (s, Es1s2_a)){

    		info = NumberExtractor(s);
    		info_vec[10] = info; 
    		input_hit_check[13] += 1;
    		if (input_hit_check [13] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}
    	else if (std::regex_search (s, Es1s2_n)){

    		info = NumberExtractor(s);
    		info_vec[11] = info; 
    		input_hit_check[14] += 1;
    		if (input_hit_check [14] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		continue;
    	}
    	else if (std::regex_search (s, frac)) {
    		info = NumberExtractor(s);
    		info_vec[12] = info;
    		input_hit_check[15] += 1;
    		if (input_hit_check [15] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
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

    p = std::make_pair ("m1", "m1");
    std::get<1>((*InteractionMap)[p]) = info_vec[4];
    std::get<2>((*InteractionMap)[p]) = info_vec[5];

    p = std::make_pair ("m1", "s1");
	std::get<1>((*InteractionMap)[p]) = info_vec[6];
    std::get<2>((*InteractionMap)[p]) = info_vec[7];

    p = std::make_pair ("s1", "m1");
    std::get<1>((*InteractionMap)[p]) = info_vec[6];
    std::get<2>((*InteractionMap)[p]) = info_vec[7];

    p = std::make_pair ("m1", "s2");
    std::get<1>((*InteractionMap)[p]) = info_vec[8];
    std::get<2>((*InteractionMap)[p]) = info_vec[9];

	p = std::make_pair ("s2", "m1");
    std::get<1>((*InteractionMap)[p]) = info_vec[8];
    std::get<2>((*InteractionMap)[p]) = info_vec[9];

	p = std::make_pair ("s2", "s1");
    std::get<1>((*InteractionMap)[p]) = info_vec[10];
    std::get<2>((*InteractionMap)[p]) = info_vec[11];

	p = std::make_pair ("s1", "s2");
    std::get<1>((*InteractionMap)[p]) = info_vec[10];
    std::get<2>((*InteractionMap)[p]) = info_vec[11];


    return info_vec;

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// End of ExtractTopologyFromFile
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


double NeighborEnergetics ( std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double, 8>* contacts, int ss_index, int x, int y, int z ){

	double   Ei = 0; 
	(*contacts) = {0, 0, 0, 0, 0, 0, 0, 0};
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ( (*LATTICE)[ss_index]->coords, x, y, z); 

	for ( std::array <int,3>& loc: ne_list) {
		ParticlePairEnergyContribution ((*LATTICE)[ss_index], (*LATTICE)[lattice_index (loc, y, z)], InteractionMap, &Ei, contacts, x, y, z);
	}

	return Ei; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void ParticlePairEnergyContribution (Particle* p1, Particle* p2, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	double* pair_energy, std::array <double,8>* contacts, int x, int y, int z) {

	double dot_product = 0;
	double theta_1     = 0;
	double theta_2     = 0;
	double magnitude   = 0;

	std::array <int,3> connvec = {0, 0, 0};

	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 

	switch ( std::get<0>((*InteractionMap)[particle_pair])[0] ) {

		case 'i':
			(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 1; 
			*pair_energy += std::get<1>((*InteractionMap)[particle_pair]);
			break;

		case 'p':
				dot_product = take_dot_product ( p1->orientation, p2->orientation );
				if (dot_product > 0.54){
					(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 1; 
					*pair_energy += std::get<1>((*InteractionMap)[particle_pair]); 
				}
				else {
					(*contacts)[std::get<4>((*InteractionMap)[particle_pair])] += 1; 
					*pair_energy += std::get<2>((*InteractionMap)[particle_pair]); 
				}	
				break;

		case 'a':
			connvec = subtract_arrays ( &(p2->coords), &(p1->coords) );
			modified_direction ( &connvec, x, y, z); 
			magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
			theta_1   = std::acos (take_dot_product (scale_arrays (1/magnitude , &connvec), Or2Dir[p1->orientation] ) );
			theta_2   = std::acos (take_dot_product (scale_arrays (-1/magnitude, &connvec), Or2Dir[p2->orientation] ) );

			if ( theta_1 + theta_2 > M_PI/2 ){

				(*contacts)[std::get<4>((*InteractionMap)[particle_pair])] += 1; 
				*pair_energy += std::get<2>((*InteractionMap)[particle_pair]); 
			}
			else {
				(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 1; 
				*pair_energy += std::get<1>((*InteractionMap)[particle_pair]); 					
			}
			break;

	}

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

double IsolatedPairParticleInteraction (Particle* p1, Particle* p2, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	int* c_idx, int x, int y, int z) {

	std::array <int,3> connvec = subtract_arrays ( &(p2->coords), &(p1->coords) ); 
	modified_direction (&connvec, x, y, z); 

	if ( std:: find (adrns.begin(), adrns.end(), connvec) == adrns.end() ){
		*c_idx = 0; 
		return 0;
	}

	double pE          = 0;
	double dot_prod    = 0;
	double magnitude   = 0;
	double theta_1     = 0;
	double theta_2     = 0;
	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 

	switch (std::get<0>( (*InteractionMap)[particle_pair])[0] ) {
		case 'i':
			pE     = std::get<1>((*InteractionMap)[particle_pair]); 
			*c_idx = std::get<3>((*InteractionMap)[particle_pair]);
		break;

		case 'p':
			dot_prod  = take_dot_product (p1->orientation, p2->orientation); 
			if (dot_prod > 0.54) {
				pE = std::get<1>((*InteractionMap)[particle_pair]);
				*c_idx = std::get<3>((*InteractionMap)[particle_pair]);
			}
			else {
				pE = std::get<2>((*InteractionMap)[particle_pair]);
				*c_idx = std::get<4>((*InteractionMap)[particle_pair]);	
			}
		break;

		case 'a':
			magnitude = std::sqrt ( connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2] );
			theta_1   = std::acos (take_dot_product (scale_arrays (1/magnitude, &connvec),  Or2Dir[p1->orientation]) ); 
			theta_2   = std::acos (take_dot_product (scale_arrays (-1/magnitude, &connvec), Or2Dir[p2->orientation]) );

			if (theta_1 + theta_2 > M_PI/2) {
				pE = std::get<2>((*InteractionMap)[particle_pair]);
				*c_idx = std::get<4>((*InteractionMap)[particle_pair]);  
			}
			else {
				pE = std::get<1>((*InteractionMap)[particle_pair]);
				*c_idx = std::get<3>((*InteractionMap)[particle_pair]); 
			}		 
		break; 
	}

	return pE; 
	
}


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

	int final_step = 0;
	bool start_recording = false; 

	// find the last step dumped in lattice file
	for ( std::string& s: contents) {

		if ( std::regex_search (s, start) ) {
			std::regex_search (s, match, numbers);
			final_step = std::stoi (match[0].str());
		}

	}

	*step_num = final_step;

	for ( std::string& s: contents) {

		// send content into stringstream 

		if ( std::regex_search (s, start) ) {
			std::regex_search ( s, match, numbers ); 
			if ( final_step == std::stoi(match[0].str()) ) {
				start_recording = true;
			}
		}

		else if ( std::regex_search (s, end) && start_recording ){
			break; 
		}

		else if (start_recording) {
			std::regex_search ( s, match, numbers );
			std::regex_token_iterator<std::string::iterator> rend; 
			std::regex_token_iterator<std::string::iterator> a ( s.begin(), s.end(), numbers );

			for ( int i=0; i<2; ++i ){
				if ( i == 0 ){
					orientation = std::stoi ( *a );
					*a++;  
				}
				else {
					*a++; 
					index = std::stoi ( *a );
				}
			}
 

			std::regex_search ( s, match, characters );
			ptype = match[0].str(); 

			// std::cout << "ptype is " << ptype << std::endl; 

			Particle* p_ptr = new Particle (location(index, x, y, z), ptype, orientation); 

			LATTICE.insert ( LATTICE.begin() + index, p_ptr); 

		}

	}

	if (!start_recording) {
		std::cout <<"Something is wrong with the restart!" << std::endl;
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

void SetUpLatticeFromScratch (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector<Particle*>* LATTICE, std::string positions, bool solvation, double frac, int x, int y, int z){

	(*Polymers) = ExtractPolymersFromFile (positions, x, y, z); 
	
	AddSolvent (LATTICE, x, y, z); 

	// populate the lattice 
	for (Polymer& pmer: (*Polymers)) {
		for (Particle*& p: pmer.chain){
			// now that I have my polymer coordinates, time to create the grand lattice 
			(*LATTICE).at(lattice_index (p->coords, y, z) ) = p; 
		}
	}
	
	AddCosolvent (Polymers, Cosolvent, LATTICE, solvation, frac, static_cast<int>((*Polymers)[0].chain.size()), x, y, z ); 
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

    for (int i{0}; i < x*y*z; ++i){
    	if ( (*LATTICE).at(i)->ptype == "s2" ){
    		(*Cosolvent).push_back((*LATTICE).at(i));
    	}
    }

    return; 

}


void AlignTheSolvationShell ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z ) {

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

void AlignTheLattice ( std::vector <Particle*>* LATTICE ) {

	for ( Particle*& p: (*LATTICE) ) {
		p->orientation = 0;
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


void AddCosolvent (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, bool solvation, double frac, int Nmonomer, int x, int y, int z) {
	std::cout << "\n--------------------------------------------------------------------\n" << std::endl;
	std::cout << "Begin adding cosolvent... \n"; 
	int nsol2         = std::floor ((x*y*z-Nmonomer)*frac); 
	std::cout << "Number of particles of cosolvent being added is " << nsol2 << "." << std::endl;
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

void InputParser(int dfreq, int lfreq, int max_iter, bool r, \
	std::string positions, std::string topology, std::string dfile, \
	std::string efile, std::string mfile, std::string stats_file, \
	std::string lattice_file_read){

	
	if (!r) {

		if (dfreq == -1 || lfreq == -1 || max_iter == -1) {
			std::cerr << "ERROR: No value for option f (frequency of dumping) (" << dfreq << ") and/or for option l (ell) (" << lfreq << ") and/or for option M (maximum number of moves to be performed) (" << max_iter << ") was provided. Exiting..." << std::endl;
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

		if ( dfreq == -1 || lfreq == -1 || max_iter == -1 ){
			std::cerr << "ERROR: No value for option f (frequency of dumping) and/or for option l (ell) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
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
	std::cout << "Input file has no overlaps! ";
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

    std::cout << "Input polymers are well-connected! \n";
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
    std::cout << std::boolalpha << "checkForOverlaps says: " << checkForOverlaps(*Polymers) << "." << std::endl; 
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

    std::cout << std::boolalpha << "checkConnectivity says: " << checkConnectivity(*Polymers, x, y, z) << "." << std::endl;
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
    // auto start = std::chrono::high_resolution_clock::now();

    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){
            ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list 

            for ( std::array <int, 3>& loc: ne_list){
            	dot_product = take_dot_product (  p->orientation, (*LATTICE)[ lattice_index(loc, y, z) ]->orientation );
            	if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "m1" ) {
            		// m-m interactions
            		if (dot_product > 0.54) {
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
			if ( dot_product > 0.54 ) {
				Energy += (*E)[2];
				(*contacts)[2] += 1;
			}
			else {
				Energy += (*E)[3];
				(*contacts)[3] += 1;
			}
		}
		else {
			Energy += (*E)[4];
			(*contacts)[4] += 1;
			// (*contacts)[5] += 1;
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
				// std::cout << "location of s1 neighbor is "; print ( loc, ", " );
				connvec   = subtract_arrays ( &(*LATTICE)[lattice_index(loc, y, z)]->coords, &(p->coords) );
				modified_direction ( &connvec, x, y, z); 
				magnitude = std::sqrt (connvec[0]*connvec[0]+connvec[1]*connvec[1]+connvec[2]*connvec[2]); 
				theta_1 = std::acos (take_dot_product ( scale_arrays( 1/magnitude, &connvec), Or2Dir[p->orientation] ) ); 
				theta_2 = std::acos (take_dot_product ( scale_arrays(-1/magnitude, &connvec), Or2Dir[(*LATTICE)[lattice_index(loc, y, z)]->orientation] ) ); 

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

    // stop = std::chrono::high_resolution_clock::now(); 
    // duration = std::chrono::duration_cast <std::chrono::microseconds> (stop-start); 
    // std::cout << "solvent energy computation took " << duration.count() << " microseconds. " << std::endl;
    
    return Energy; 
}

//==============================================================================================
//==============================================================================================

double CalculateEnergyRevamped (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple<std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* contacts, int x, int y, int z) {
    
    double Energy {0.0};
    (*contacts) = {0,0,0,0,0,0,0,0}; 

    std::array <std::array <int,3>, 26> ne_list; 

    // run energy computations for every monomer bead 
    // auto start = std::chrono::high_resolution_clock::now();

    for (Polymer& pmer: (*Polymers)) {
		for (Particle*& p: pmer.chain){
            ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list
            for ( std::array <int,3>& loc: ne_list) {
            	PairInteractionForRevampedEnergyCalc (p, (*LATTICE)[ lattice_index(loc, y, z) ], InteractionMap, &Energy, contacts, x, y, z);
        	}
        }
	}
    
	for ( Particle*& p: *Cosolvent ){
		// std::cout << "Entered cosolvent loop. " << std::endl;
		ne_list = obtain_ne_list ( p->coords, x, y, z );
		for ( std::array <int,3>& loc: ne_list ){
			if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "m1" || (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s2"){
				continue; 
			}
			ParticlePairEnergyContribution (p, (*LATTICE)[ lattice_index(loc, y, z) ], InteractionMap, &Energy, contacts, x, y, z);
		}
	}

	return Energy; 
}


void PairInteractionForRevampedEnergyCalc (Particle* p1, Particle* p2, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	double* pair_energy, std::array <double,8>* contacts, int x, int y, int z) {

	double dot_product = 0;
	double theta_1     = 0;
	double theta_2     = 0;
	double magnitude   = 0;

	std::array <int,3> connvec = {0, 0, 0};

	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 

	std::string interaction_type = std::get<0>((*InteractionMap)[particle_pair]);

	if (interaction_type == "isotropic") {

		if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
			// std::cout << "isotropic m1-m1..." << std::endl;
			// std::cout << "particle_pair = " << particle_pair.first << ", " << particle_pair.second << std::endl;
			(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 0.5; 
			*pair_energy += std::get<1>((*InteractionMap)[particle_pair])*0.5; 
		}
		else {
			// std::cout << "particle_pair = " << particle_pair.first << ", " << particle_pair.second << std::endl;
			(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 1; 	
			*pair_energy += std::get<1>((*InteractionMap)[particle_pair]); 
		}

	}
	else if (interaction_type == "parallel") {

		dot_product = take_dot_product ( p1->orientation, p2->orientation );
		if (dot_product > 0.54){
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				// std::cout << "aligned parallel m1-m1..." << std::endl;
				(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 0.5; 
				*pair_energy += std::get<1>((*InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 1; 	
				*pair_energy += std::get<1>((*InteractionMap)[particle_pair]); 
			}	
		}
		else {
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				// std::cout << "misaligned parallel m1-m1..." << std::endl;
				(*contacts)[std::get<4>((*InteractionMap)[particle_pair])] += 0.5; 
				*pair_energy += std::get<2>((*InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(*contacts)[std::get<4>((*InteractionMap)[particle_pair])] += 1; 
				*pair_energy += std::get<2>((*InteractionMap)[particle_pair]); 
			}
		}

	}

	else if (interaction_type == "antiparallel") {
		connvec = subtract_arrays ( &(p2->coords), &(p1->coords) );
		modified_direction (&connvec, x, y, z); 
		magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
		theta_1   = std::acos (take_dot_product (scale_arrays (1/magnitude , &connvec), Or2Dir[p1->orientation] ) );
		theta_2   = std::acos (take_dot_product (scale_arrays (-1/magnitude, &connvec), Or2Dir[p2->orientation] ) );

		if ( theta_1 + theta_2 > M_PI/2 ){
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ) {
				// std::cout << "misaligned antiparallel m1-m1..." << std::endl;
				(*contacts)[std::get<4>((*InteractionMap)[particle_pair])] += 0.5; 
				*pair_energy += std::get<2>((*InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(*contacts)[std::get<4>((*InteractionMap)[particle_pair])] += 1; 	
				*pair_energy += std::get<2>((*InteractionMap)[particle_pair]); 				
			}
		}
		else {
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				// std::cout << "aligned antiparallel m1-m1..." << std::endl;
				(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 0.5; 
				*pair_energy += std::get<1>((*InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(*contacts)[std::get<3>((*InteractionMap)[particle_pair])] += 1; 
				*pair_energy += std::get<1>((*InteractionMap)[particle_pair]); 
			}			
		}
	}

	return; 

}


//==============================================================================================
//==============================================================================================



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
    dump_file << "Chain regrowth with ori flip       - attempts: " << (*attempts)[3] <<", acceptances: " << (*acceptances)[3] << ", acceptance fraction: " << static_cast<double>((*acceptances)[3])/static_cast<double>((*attempts)[3]) << std::endl; 
    dump_file << "Solvent flips without bias         - attempts: " << (*attempts)[4] <<", acceptances: " << (*acceptances)[4] << ", acceptance fraction: " << static_cast<double>((*acceptances)[4])/static_cast<double>((*attempts)[4]) << std::endl;
    dump_file << "Solvation shell flip with bias     - attempts: " << (*attempts)[5] <<", acceptances: " << (*acceptances)[5] << ", acceptance fraction: " << static_cast<double>((*acceptances)[5])/static_cast<double>((*attempts)[5]) << std::endl;
    dump_file << "Polymer flips                      - attempts: " << (*attempts)[6] <<", acceptances: " << (*acceptances)[6] << ", acceptance fraction: " << static_cast<double>((*acceptances)[6])/static_cast<double>((*attempts)[6]) << std::endl;
    dump_file << "Solvent exchange with bias         - attempts: " << (*attempts)[7] <<", acceptances: " << (*acceptances)[7] << ", acceptance fraction: " << static_cast<double>((*acceptances)[7])/static_cast<double>((*attempts)[7]) << std::endl;
    dump_file << "Solvent exchange without bias      - attempts: " << (*attempts)[8] <<", acceptances: " << (*acceptances)[8] << ", acceptance fraction: " << static_cast<double>((*acceptances)[8])/static_cast<double>((*attempts)[8]) << std::endl;

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

void dumpSolvation( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int step, std::string filename, int x, int y, int z ) {
	std::ofstream dump_file (filename, std::ios::app); 
	std::set <int> solvent_indices  ;
	std::set <int> cosolvent_indices;
	
	for ( Polymer& pmer: (*Polymers) ) {
		for ( Particle*& p: pmer.chain ) {
			std::array <std::array<int,3>, 26> ne_list = obtain_ne_list (p->coords, x, y, z);
			for ( std::array <int,3>& ne: ne_list) {
				if ((*LATTICE)[lattice_index(ne, y, z)]->ptype[1] == '1' && (*LATTICE)[lattice_index(ne, y, z)]->ptype[0]=='s'){
					solvent_indices.insert(lattice_index(ne, y, z));
				}
				else if ((*LATTICE)[lattice_index(ne, y, z)]->ptype[1] == '2' && (*LATTICE)[lattice_index(ne, y, z)]->ptype[0]=='s'){
					cosolvent_indices.insert(lattice_index(ne, y, z));
				}
			}
		}
	}

	int total_solvent_particles   {static_cast<int>(solvent_indices.size())};
	int total_cosolvent_particles {static_cast<int>(cosolvent_indices.size())};
	dump_file << total_solvent_particles + total_cosolvent_particles << " | " << total_solvent_particles << " | " << total_cosolvent_particles << " | " << step << "\n";

	return;

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of dumpSolvation 
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

void TailRotation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

	
	// get the neighborlist of particle at index 1 
	std::cout << "*E[0] = " << (*E)[0] << std::endl;
	std::array <int,3> loc_0 = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3> loc_1 = (*Polymers)[index].chain[1]->coords;

	std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

	std::array <double,8> cs             = {0, 0, 0, 0, 0, 0, 0, 0};  
	std::array <double,8> cm             = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> cs_n           = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> cm_n           = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> final_contacts = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> c_contacts     = {0, 0, 0, 0, 0, 0, 0, 0};

	int choice = rng_uniform(0, 25); 

	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	double Es = 0, Es_n = 0; 
	double Em = 0, Em_n = 0; 
	double Ef = 0; 
	double Epair   = 0; 
	double Epair_n = 0; 

	int c_idx   = 0;
	int c_idx_n = 0;


	if ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] == 's' ){

		// find the energetic interaction for the solvent molecule 
		Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index (ne_list[choice], y, z), x, y, z);
		Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index (loc_0, y, z), x, y, z);  
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx, x, y, z);

		(*Polymers)[index].chain[0]->coords = ne_list[choice]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 
		
		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]           = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[0]; 

		
		Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);
		Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index (loc_0, y, z), x, y, z);
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx_n, x, y, z);  

	}
	else {
		*IMP_BOOL = false; 
		return; 
	}
	

	Ef = *sysEnergy - (Es + Em - Epair) + (Em_n + Es_n - Epair_n); 
	final_contacts = add_arrays ( subtract_arrays(*contacts, add_arrays (cs, cm)), add_arrays(cs_n, cm_n) );
	final_contacts[c_idx]   += 1;
	final_contacts[c_idx_n] -= 1;

	// DELETE THIS LATER 
	double energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &c_contacts, x, y, z); 
	if (Ef != energy_n || final_contacts != c_contacts) {
		std::cout << "The energy or contacts have become messed up in end rotation... " << std::endl;
		std::cout << "Ef (calc) = " << Ef << ", energy_n (true) = " << energy_n << ". " << std::endl; 
		std::cout << "final_contacts (calc) = "; print (final_contacts, ", "); std::cout << "c_contacts (true) = "; print(c_contacts);
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
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

	// get the neighborlist of particle at index 1 
	std::array <int,3> loc_0 = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3> loc_1 = (*Polymers)[index].chain[1]->coords;

	std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

	std::array <double,8> cs             = {0, 0, 0, 0, 0, 0, 0, 0};  
	std::array <double,8> cm             = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> cs_n           = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> cm_n           = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> final_contacts = {0, 0, 0, 0, 0, 0, 0, 0}; 

	int choice = rng_uniform(0, 25); 

	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	double Es      = 0; 
	double Es_n    = 0; 
	double Em      = 0; 
	double Em_n    = 0; 
	double Ef      = 0; 
	double Epair   = 0; 
	double Epair_n = 0; 
	int    c_idx   = 0;
	int    c_idx_n = 0;


	if ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] == 's' ){

		// find the energetic interaction for the solvent molecule 
		Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[choice], y, z), x, y, z);
		Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_0, y, z), x, y, z);
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx, x, y, z);

		// setting up the switch 
		(*Polymers)[index].chain[0]->coords = ne_list[choice]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 
		
		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]           = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[0]; 

		// get new neighbor energetics 
		Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);
		Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index (loc_0, y, z), x, y, z);  
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx_n, x, y, z);

	}
	else {
		*IMP_BOOL = false; 
		return; 
	}

	Ef = *sysEnergy - (Es + Em - Epair) + (Em_n + Es_n - Epair_n); 
	final_contacts = add_arrays( subtract_arrays(*contacts, add_arrays (cs, cm)), add_arrays(cs_n, cm_n) );
	final_contacts[c_idx]   += 1; 
	final_contacts[c_idx_n] -= 1; 

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

		energies.push_back( CalculateEnergy (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z) );
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

void HeadRotation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

    // get the neighborlist of particle at index 1 
    std::cout << "*E[0] = " << (*E)[0] << std::endl;
    int dop = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2]->coords; 

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 

    std::array <double,8> cs             = {0, 0, 0, 0, 0, 0, 0, 0};  
    std::array <double,8> cm             = {0, 0, 0, 0, 0, 0, 0, 0}; 
    std::array <double,8> cs_n           = {0, 0, 0, 0, 0, 0, 0, 0}; 
    std::array <double,8> cm_n           = {0, 0, 0, 0, 0, 0, 0, 0}; 
    std::array <double,8> final_contacts = {0, 0, 0, 0, 0, 0, 0, 0}; 
	std::array <double,8> c_contacts     = {0, 0, 0, 0, 0, 0, 0, 0};

	int choice = rng_uniform(0, 25); 

	// made the change to the polymer
	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	double Es = 0, Es_n = 0; 
	double Em = 0, Em_n = 0; 
	double Ef = 0; 
	double Epair   = 0;
	double Epair_n = 0;

	int c_idx   = 0; 
	int c_idx_n = 0;

	if ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] == 's' ){
		
		// find the energetic interactions for the solvent molecule 
		Es      = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index (ne_list[choice], y, z), x, y, z);
		Em      = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index (loc_0, y, z), x, y, z);  
		Epair   = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx, x, y, z);

		(*Polymers)[index].chain[dop-1]->coords = ne_list[choice]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 

		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[dop-1]; 

		// 
		Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);
		Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index (loc_0, y, z), x, y, z);  
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx_n, x, y, z);
	
	}
	else {
		*IMP_BOOL = false; 
		return; 
	}

	Ef = *sysEnergy - (Es + Em - Epair) + (Em_n + Es_n - Epair_n);
	final_contacts = add_arrays( subtract_arrays(*contacts, add_arrays (cs, cm)), add_arrays(cs_n, cm_n) );
	final_contacts[c_idx]   += 1;
	final_contacts[c_idx_n] -= 1;
	
	double energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &c_contacts, x, y, z); 
	if (Ef != energy_n || final_contacts != c_contacts) {
		std::cout << "Either energy or contacts is messed up in end rotation... " << std::endl; 
		std::cout << "Ef (calculated) = " << Ef << ", energy_n (real) = " << energy_n << ". " << std::endl; 
		std::cout << "final_contacts (calculated)= "; print (final_contacts, ", "); std::cout << "c_contacts (real)= "; print(c_contacts); 
		exit (EXIT_FAILURE); 
	}
	

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
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap,\
	std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int index, int x, int y, int z) {

    // get the neighborlist of particle at index 1 
    int dop = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2]->coords; 

    // these are contacts arrays which are necessary to check validity of code 
    std::array <double,8> cs             = {0,0,0,0,0,0,0,0};
    std::array <double,8> cm             = {0,0,0,0,0,0,0,0};
    std::array <double,8> cs_n           = {0,0,0,0,0,0,0,0};
    std::array <double,8> cm_n           = {0,0,0,0,0,0,0,0};
    std::array <double,8> final_contacts = {0,0,0,0,0,0,0,0}; 

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z); 
	int choice = rng_uniform(0, 25); 

	// made the change to the polymer
	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	int    c_idx   = 0;
	int    c_idx_n = 0;
	double Es      = 0; 
	double Es_n    = 0; 
	double Em      = 0; 
	double Em_n    = 0; 
	double Ef      = 0; 
	double Epair   = 0;
	double Epair_n = 0;

	if ( (*LATTICE)[lattice_index (ne_list[choice], y, z)]->ptype[0] == 's' ){
		
		// find the energetic interactions for the solvent molecule 
		Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index (ne_list[choice], y, z), x, y, z);
		Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index (loc_0, y, z), x, y, z);  
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx, x, y, z);

		(*Polymers)[index].chain[dop-1]->coords = ne_list[choice]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]->coords = loc_0; 

		// do the switch 
		// take the pointer of the solvent, and put it where the monomer was on the lattice 
		(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (ne_list[choice], y, z)]; 
		(*LATTICE)[ lattice_index (ne_list[choice], y, z)]  = (*Polymers)[index].chain[dop-1]; 

		Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index (ne_list[choice], y, z), x, y, z);
		Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index (loc_0, y, z), x, y, z);  
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[choice], y, z)], (*LATTICE)[lattice_index (loc_0, y, z)], InteractionMap, &c_idx_n, x, y, z);
	
	}
	else {
		*IMP_BOOL = false; 
		return; 
	}

	Ef = *sysEnergy - (Es + Em - Epair) + (Em_n + Es_n - Epair_n);
	final_contacts = add_arrays( subtract_arrays(*contacts, add_arrays (cs, cm)), add_arrays(cs_n, cm_n) );
	final_contacts[c_idx]   += 1;
	final_contacts[c_idx_n] -= 1;

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

		energies.push_back( CalculateEnergy (Polymers, Cosolvent, LATTICE, E, &c_contacts, x, y, z) );
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

void EndRotation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int index, int x, int y, int z){

	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        std::cout << "Zero index rotation!" << std::endl;
        TailRotation_UNBIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    else {
        std::cout << "Final index rotation!" << std::endl;
        HeadRotation_UNBIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void EndRotation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap,\
	std::array <double,8>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, \
	int index, int x, int y, int z){

	unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
	std::mt19937 generator(seed); 
	std::uniform_int_distribution<int> distribution (0,1); 
	int num = distribution(generator); 
	switch (num) {
		case (0):
			// std::cout << "Zero index rotation!" << std::endl;
			TailRotation_UNBIASED (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			break; 
		case (1):
			// std::cout << "Final index rotation!" << std::endl;
			HeadRotation_UNBIASED (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			break;
	}

	return;

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
    int num = distribution(generator); 
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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void forward_reptation_with_tail_biting_new (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, double* frontflow_energy, int deg_poly, int index, int x, int y, int z){

	// start performing exchanges 
	// doubles for energy transfer 
	int    c_idx       = 0;
	int    c_idx_n     = 0; 
	double Em1         = 0;
	double Em2         = 0; 
	double Em1_n       = 0; 
	double Em2_n       = 0; 
	double Esys        = *frontflow_energy; 
	double Epair       = 0;
	double Epair_n     = 0; 

	std::array <double,8> cm1               = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2               = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm1_n             = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2_n             = {0,0,0,0,0,0,0,0};
	std::array <double,8> current_contacts  = *contacts;
	std::array <int,3>    loc_m1             = {0,0,0}; 
	std::array <int,3>    loc_m2             = {0,0,0}; 

	for (int i{0}; i<deg_poly-1; ++i){
		
		loc_m1 = (*Polymers)[index].chain[i]->coords;
		loc_m2 = (*Polymers)[index].chain[i+1]->coords;
		
		// find energies 
		Em1   = NeighborEnergetics (LATTICE, InteractionMap, &cm1, lattice_index (loc_m1, y, z), x, y, z);
		Em2   = NeighborEnergetics (LATTICE, InteractionMap, &cm2, lattice_index (loc_m2, y, z), x, y, z);
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m1, y, z)], (*LATTICE)[lattice_index (loc_m2, y, z)], InteractionMap, &c_idx, x, y, z);

		(*LATTICE)[ lattice_index(loc_m1, y, z) ] = (*LATTICE)[ lattice_index (loc_m2, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ] = (*Polymers)[index].chain[i]; 
		(*LATTICE)[ lattice_index(loc_m1, y, z) ]->coords = loc_m1; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ]->coords = loc_m2; 

		Em1_n   = NeighborEnergetics (LATTICE, InteractionMap, &cm1_n, lattice_index (loc_m1, y, z), x, y, z);
		Em2_n   = NeighborEnergetics (LATTICE, InteractionMap, &cm2_n, lattice_index (loc_m2, y, z), x, y, z);
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m1, y, z)], (*LATTICE)[lattice_index (loc_m2, y, z)], InteractionMap, &c_idx_n, x, y, z);

		current_contacts = add_arrays ( subtract_arrays (current_contacts, add_arrays (cm1, cm2)), add_arrays (cm1_n, cm2_n) ); 
		current_contacts[c_idx]   += 1;
		current_contacts[c_idx_n] -= 1;
		Esys = Esys - (Em1+Em2-Epair) + (Em1_n+Em2_n-Epair_n); 

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

void forward_reptation_without_tail_biting_new (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, std::array <int,3>* to_slither, double* frontflow_energy, int deg_poly, int index, int x, int y, int z){

	// start performing exchanges 
	// doubles for energy transfer 
	int    c_idx      = 0;
	int    c_idx_n    = 0;
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	double Epair      = 0;
	double Epair_n    = 0;

	std::array <double,8> cm               = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs               = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n             = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n             = {0,0,0,0,0,0,0,0};
	std::array <double,8> current_contacts = *contacts;
	std::array <int,3>    loc_s            = *to_slither; 
	std::array <int,3>    loc_m            = {0,0,0}; 

	for (int i{0}; i<deg_poly; ++i){
		
		loc_s = *to_slither; 
		loc_m = (*Polymers)[index].chain[deg_poly-1-i]->coords;
		
		// find energies 
		Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index (loc_s, y, z), x, y, z);
		Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index (loc_m, y, z), x, y, z);
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_s, y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx, x, y, z);

		(*LATTICE)[ lattice_index(loc_m, y, z) ] = (*LATTICE)[ lattice_index (loc_s, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i]; 
		(*LATTICE)[ lattice_index(loc_m, y, z) ]->coords = loc_m; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ]->coords = loc_s; 

		Es_n    = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index (loc_s, y, z), x, y, z);
		Em_n    = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index (loc_m, y, z), x, y, z);
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_s, y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z);

		current_contacts = add_arrays ( subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
		current_contacts[c_idx]   += 1;
		current_contacts[c_idx_n] -= 1;

		Esys = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 

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

void backward_reptation_with_head_butting_new (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, double* frontflow_energy, int deg_poly, int index, int x, int y, int z){

	// start performing exchanges 
	// doubles for energy transfer 
	int    c_idx       = 0;
	int    c_idx_n     = 0;
	double Em1         = 0;
	double Em2         = 0; 
	double Em1_n       = 0; 
	double Em2_n       = 0; 
	double Esys        = *frontflow_energy; 
	double Epair       = 0;
	double Epair_n     = 0;

	std::array <double,8> cm1                = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2                = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm1_n              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm2_n              = {0,0,0,0,0,0,0,0};
	std::array <double,8> current_contacts   = *contacts;
	std::array <int,3>    loc_m1             = {0,0,0}; 
	std::array <int,3>    loc_m2             = {0,0,0}; 

	for (int i{0}; i<deg_poly-1; ++i){
		
		loc_m1 = (*Polymers)[index].chain[deg_poly-(i+1)]->coords;
		loc_m2 = (*Polymers)[index].chain[deg_poly-(i+2)]->coords;
		
		// find energies 
		Em1 = NeighborEnergetics (LATTICE, InteractionMap, &cm1, lattice_index (loc_m1, y, z), x, y, z);
		Em2 = NeighborEnergetics (LATTICE, InteractionMap, &cm2, lattice_index (loc_m2, y, z), x, y, z);
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m1, y, z)], (*LATTICE)[lattice_index (loc_m2, y, z)], InteractionMap, &c_idx, x, y, z);


		(*LATTICE)[ lattice_index(loc_m1, y, z) ] = (*LATTICE)[ lattice_index (loc_m2, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ] = (*Polymers)[index].chain[deg_poly-(i+1)]; 
		(*LATTICE)[ lattice_index(loc_m1, y, z) ]->coords = loc_m1; 
		(*LATTICE)[ lattice_index(loc_m2, y, z) ]->coords = loc_m2; 

		Em1_n = NeighborEnergetics (LATTICE, InteractionMap, &cm1_n, lattice_index (loc_m1, y, z), x, y, z);
		Em2_n = NeighborEnergetics (LATTICE, InteractionMap, &cm2_n, lattice_index (loc_m2, y, z), x, y, z);
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m1, y, z)], (*LATTICE)[lattice_index (loc_m2, y, z)], InteractionMap, &c_idx_n, x, y, z);

		current_contacts = add_arrays ( subtract_arrays (current_contacts, add_arrays (cm1, cm2)), add_arrays (cm1_n, cm2_n) ); 
		current_contacts[c_idx]   += 1;
		current_contacts[c_idx_n] -= 1;

		Esys = Esys - (Em1+Em2-Epair) + (Em1_n+Em2_n-Epair_n); 

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
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, std::array <int,3>* to_slither, double* frontflow_energy, int deg_poly, int index, int x, int y, int z) {

	// start performing exchanges 
	// doubles for energy transfer 
	int    c_idx      = 0;
	int    c_idx_n    = 0;
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	double Epair      = 0;
	double Epair_n    = 0;

	std::array <double,8> cm               = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs               = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n             = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n             = {0,0,0,0,0,0,0,0};
	std::array <double,8> current_contacts = *contacts;
	std::array <int,3>    loc_s            = *to_slither; 
	std::array <int,3>    loc_m            = {0,0,0}; 

	for (int i{0}; i<deg_poly; ++i){
		
		loc_s = *to_slither; 
		loc_m = (*Polymers)[index].chain[i]->coords;
		
		// find energies 
		Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index (loc_s, y, z), x, y, z);
		Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index (loc_m, y, z), x, y, z);
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (loc_s, y, z)], InteractionMap, &c_idx, x, y, z);

		(*LATTICE)[ lattice_index(loc_m, y, z) ] = (*LATTICE)[ lattice_index (loc_s, y, z) ]; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ] = (*Polymers)[index].chain[i]; 
		(*LATTICE)[ lattice_index(loc_m, y, z) ]->coords = loc_m; 
		(*LATTICE)[ lattice_index(loc_s, y, z) ]->coords = loc_s; 

		Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index (loc_s, y, z), x, y, z);
		Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index (loc_m, y, z), x, y, z);
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (loc_s, y, z)], InteractionMap, &c_idx_n, x, y, z);

		current_contacts = add_arrays ( subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
		current_contacts[c_idx]   += 1;
		current_contacts[c_idx_n] -= 1;

		Esys = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 

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

void ForwardReptation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	std::cout << "In forward reptation unbiased." <<std::endl;
	std::cout << "*E[0] = " << (*E)[0] << std::endl;
	int                   deg_poly         = (*Polymers)[index].deg_poly; 
	std::array <int,3>    loc0             = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>    locf             = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,8> current_contacts = *contacts; 
	std::array <double,8> copy_contacts    = *contacts; 
	double                Esys             = *sysEnergy; 


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
		forward_reptation_with_tail_biting_new    (Polymers, LATTICE, InteractionMap, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		// IMPORTANT STEP 
		std::cout << "Without tail biting... " << std::endl;
		forward_reptation_without_tail_biting_new (Polymers, LATTICE, InteractionMap, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
	}

	else {
		*IMP_BOOL = false; 
		return; 
	}
	(*Polymers)[0].printChainCoords(); 
	std::cout << "First check structures..." << std::endl;
	CheckStructures (Polymers, Cosolvent, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl; 
	// calculate energy of current state 
	// delete later 

	double energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z);
	if ( Esys != energy_n || current_contacts != copy_contacts){
		std::cout << "Either energy or contacts is messed up in forward reptation... " << std::endl; 
		std::cout << "Esys (calc) = " << Esys << ", energy_n (real) = " << energy_n << ". " << std::endl; 
		std::cout << "current_contacts (calc) = "; print (current_contacts, ", "); std::cout << "copy_contacts (real) = "; print(copy_contacts); 
		exit (EXIT_FAILURE); 
	}
	// delete above 

	if ( MetropolisAcceptance (*sysEnergy, energy_n, temperature) ) {
		std::cout << "Accepted!" << std::endl;
		*sysEnergy = Esys;
		*contacts  = current_contacts; 
	}

	else {
		// revert back to old state 
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
	CheckStructures (Polymers, Cosolvent, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl; 

	return; 
}

//==============================================================================================
//==============================================================================================

void ForwardReptation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	// std::cout << "In forward reptation unbiased." <<std::endl;

	int                   deg_poly         = (*Polymers)[index].deg_poly; 
	std::array <int,3>    loc0             = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>    locf             = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,8> current_contacts = *contacts; 
	double                Esys             = *sysEnergy; 


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
		forward_reptation_with_tail_biting_new    (Polymers, LATTICE, InteractionMap, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		// IMPORTANT STEP 
		// std::cout << "Without tail biting... " << std::endl;
		forward_reptation_without_tail_biting_new (Polymers, LATTICE, InteractionMap, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
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

void BackwardReptation_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z) {

	std::cout << "In backward reptation unbiased." <<std::endl;
	std::cout << "*E[0] = " << (*E)[0] << std::endl;

	int                   deg_poly         = (*Polymers)[index].deg_poly; 
	std::array <int,3>    loc0             = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>    locf             = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,8> current_contacts = *contacts; 
	std::array <double,8> copy_contacts    = *contacts;
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
		backward_reptation_with_head_butting_new (Polymers, LATTICE, InteractionMap, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		std::cout << "Without head butting... " << std::endl;
		backward_reptation_without_head_butting_new (Polymers, LATTICE, InteractionMap, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
	}

	else {
		*IMP_BOOL = false; 
		return; 
	}

	(*Polymers)[0].printChainCoords();

	std::cout << "First check structures..." << std::endl;
	CheckStructures (Polymers, Cosolvent, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl;

	// delete later 
	
	double energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z);
	if ( Esys != energy_n || current_contacts != copy_contacts){
		std::cout << "Either energy or contacts is messed up in forward reptation... " << std::endl; 
		std::cout << "Esys (real) = " << Esys << ", energy_n (calc) = " << energy_n << ". " << std::endl; 
		std::cout << "current_contacts (real) = "; print (current_contacts, ", "); std::cout << "copy_contacts (calc) = "; print(copy_contacts); 
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
	CheckStructures (Polymers, Cosolvent, LATTICE, x, y, z);
	std::cout << "PASSED!" << std::endl;

	return; 
}

//==============================================================================================
//==============================================================================================


void BackwardReptation_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

	// std::cout << "In backward reptation unbiased." <<std::endl;

	int                   deg_poly         = (*Polymers)[index].deg_poly; 
	std::array <int,3>    loc0             = (*Polymers)[index].chain[0]->coords; 
	std::array <int,3>    locf             = (*Polymers)[index].chain[deg_poly-1]->coords; 
	std::array <double,8> current_contacts = *contacts; 
	double                Esys             = *sysEnergy; 

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
		backward_reptation_with_head_butting_new    (Polymers, LATTICE, InteractionMap, &current_contacts, &Esys, deg_poly, index, x, y, z); 
	}
	else if ( (*LATTICE)[ lattice_index (to_slither, y, z) ]->ptype[0] == 's' ){
		// std::cout << "Without head butting... " << std::endl;
		backward_reptation_without_head_butting_new (Polymers, LATTICE, InteractionMap, &current_contacts, &to_slither, &Esys, deg_poly, index, x, y, z);
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


void Reptation_UNBIASED_debug (std::vector<Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int index, int x, int y, int z){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 

    if (num==0){
        std::cout << "Backward reptation only!" << std::endl;
        BackwardReptation_UNBIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
        return; 
    }
    else {
        std::cout << "Forward reptation only!" << std::endl;
        ForwardReptation_UNBIASED_debug  (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
        return; 
    }
    
}

//==============================================================================================
//==============================================================================================

void Reptation_UNBIASED (std::vector<Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 

    if (num==0){
        // std::cout << "Backward reptation only!" << std::endl;
        BackwardReptation_UNBIASED (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
        return; 
    }
    else {
        // std::cout << "Forward reptation!" << std::endl;
        ForwardReptation_UNBIASED  (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
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

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void SolventFlip_UNBIASED_debug ( std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* E, std::array<double,8>* contacts, \
	bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z ){

	// number of sites to flip 
	std::cout << "*E[0] = " << (*E)[0] << std::endl;
	int deg_poly = static_cast<int>((*Polymers)[index].chain.size());
	int nflips = rng_uniform (1, static_cast<int>(std::ceil(deg_poly/4.0*deg_poly/4.0*deg_poly/4.0) ) ); 
	
	std::array <double,8> copy_contacts            = *contacts; 
	std::array <double,8> current_contacts         = *contacts; 


	std::vector <int> indices;
	std::vector <int> oris;  
	indices.reserve(nflips); 
	oris.reserve(nflips);

	int r_idx {-1}; 
	int s {0}; 

	std::array <double,8> cs; 
	std::array <double,8> cs_n; 

	double Es   = 0;
	double Es_n = 0;
	double Ef   = *sysEnergy; 
	for (int i{0}; i<nflips; ++i){
		// generate an index 
		r_idx = rng_uniform (0, x*y*z-1); 
		if ( (*LATTICE)[r_idx]->ptype[0] == 's' ){
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, r_idx, x, y, z); 
			std::vector<int>::iterator it = std::find (indices.begin(), indices.end(), r_idx);
			// if it is a solvent, perturb orientation 
			if ( it == indices.end() ){ 
				indices.push_back (r_idx); 
				oris.push_back ((*LATTICE)[r_idx]->orientation); 
				// generate an orientation  
				(*LATTICE)[r_idx]->orientation = rng_uniform(0, 25); 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, r_idx, x, y, z); 
				current_contacts = add_arrays( subtract_arrays(current_contacts, cs), cs_n ); 
				Ef += Es_n - Es; 
			} 
			else { 
				(*LATTICE)[r_idx]->orientation = rng_uniform(0, 25); 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, r_idx, x, y, z); 
				current_contacts = add_arrays( subtract_arrays(current_contacts, cs), cs_n ); 
				Ef += Es_n - Es; 
			}
		}
		else {
			*IMP_BOOL = false; 
			break;
		}
	}
	// auto stop  = std::chrono::high_resolution_clock::now(); 
	// auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 
	// std::cout << "T(swap) = " << duration.count() << std::endl;
	// std::cout << "Does it reach here 2" << std::endl;
	if (*IMP_BOOL){
		// calculate the energy 
		// delete later 
		// start  = std::chrono::high_resolution_clock::now(); 
		double energy = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
		// stop  = std::chrono::high_resolution_clock::now(); 
		// duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 
		// std::cout << "T(reg) = " << duration.count() << std::endl;

		// std::cout << "Ef = " << Ef << ", energy = " << energy << ". " << std::endl; 
		if (Ef != energy || current_contacts != copy_contacts) {
			std::cout << "Either energy or contacts is messed up in end rotation... " << std::endl;
			std::cout << "Ef = " << Ef << ", energy = " << energy << ". " << std::endl; 
			std::cout << "final_contacts = "; print (current_contacts, ", "); std::cout << "copy_contacts = "; print(copy_contacts);
			exit(EXIT_FAILURE);
		} 
	
		// delete later 

		// perform the metropolis acceptance 
		if ( MetropolisAcceptance (*sysEnergy, Ef, temperature) ){
			*sysEnergy = Ef;
			*contacts  = current_contacts; 
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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void SolventFlip_UNBIASED ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z ){

	// number of sites to flip 
	int deg_poly = static_cast<int> ( (*Polymers)[index].chain.size() );
	int nflips = rng_uniform (1, static_cast<int>(std::ceil(deg_poly/2.0*deg_poly/2.0*deg_poly/2.0) ) ); 
	
	std::array <double,8> current_contacts         = *contacts; 

	std::vector <int> indices;
	std::vector <int> oris;  
	indices.reserve(nflips); 
	oris.reserve(nflips);

	int r_idx {-1}; 

	std::array <double,8> cs; 
	std::array <double,8> cs_n; 

	double Es   = 0;
	double Es_n = 0;
	double Ef   = *sysEnergy; 

	for (int i{0}; i<nflips; ++i){
		// generate an index 
		r_idx = rng_uniform (0, x*y*z-1); 
		if ( (*LATTICE)[r_idx]->ptype[0] == 's' ){
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, r_idx, x, y, z); 
			std::vector<int>::iterator it = std::find (indices.begin(), indices.end(), r_idx);
			// if it is a solvent, perturb orientation 
			if ( it == indices.end() ){ 
				indices.push_back (r_idx); 
				oris.push_back ((*LATTICE)[r_idx]->orientation); 
				// generate an orientation  
				(*LATTICE)[r_idx]->orientation = rng_uniform(0, 25); 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, r_idx, x, y, z); 
				current_contacts = add_arrays( subtract_arrays(current_contacts, cs), cs_n ); 
				Ef += Es_n - Es; 
			} 
			else { 
				(*LATTICE)[r_idx]->orientation = rng_uniform(0, 25); 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, r_idx, x, y, z); 
				current_contacts = add_arrays( subtract_arrays(current_contacts, cs), cs_n ); 
				Ef += Es_n - Es; 

			}
		}
		else {
			*IMP_BOOL = false; 
			break;
		}

	}
	// std::cout << "Does it reach here 2" << std::endl;
	if (*IMP_BOOL){

		// perform the metropolis acceptance 
		if ( MetropolisAcceptance (*sysEnergy, Ef, temperature) ) {
			*sysEnergy = Ef;
			*contacts  = current_contacts; 
		}
		else {
			*IMP_BOOL = false; 
			// reverse all the perturbations performed 
			for ( int i{0}; i < static_cast<int>(indices.size()); ++i ) {
				(*LATTICE)[ indices[i] ]->orientation = oris[i];
			}
		}
	}

	else {
		// reverse all the perturbations performed 
		for ( int i{0}; i < static_cast<int>(indices.size()); ++i ) {
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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of FirstSolvationShellFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void SolvationShellFlip_BIASED_remake2 (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z ){

	std::set <int> solvation_shell_set; 
	std::array <std::array<int,3>, 26> ne_list; 
	// int dop = static_cast<int>((*Polymers)[0].chain.size() ); 

	// get the first solvation shell 
	// auto start = std::chrono::high_resolution_clock::now(); 
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
	int nflip = rng_uniform(1, static_cast<int>(solvation_shell_indices.size()/2 ) ); 			// number of sites to be flipped 

	// std::cout << "Is this being hit1?" << std::endl;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// std::iota (solvation_shell_indices.begin(), solvation_shell_indices.end(), 0); 
	std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));
    
	std::vector <int> old_ori; 					// vector to hold old orientations 
	std::vector <int> new_ori; 					// vector to hold new orientations 
	old_ori.reserve(nflip); 					// reserving... 
	new_ori.reserve(nflip); 					// reserving... 

	// energy store for boltzmann sampling 
	// instantiating a bunch of variables for the process 
	int                                 ntest              = 5; 
	std::array <double,5>               energies           = {0,0,0,0,0}; 
	std::array <double,5>               boltzmann          = {0,0,0,0,0};
	std::array <int,5>                  orientations       = {0,0,0,0,0}; 
	double                              rboltzmann         = 0;  
	double                              frontflow_energy   = 0; 
	double                              prob_n_to_o        = 1; 
	double                              prob_o_to_n        = 1; 
	std::array<std::array<double,8>,5>  contacts_store     = {*contacts, *contacts, *contacts, *contacts, *contacts}; 
	std::array <double,8>               contacts_sys       = *contacts; 
	std::array <double,8>               contacts_i         = {0,0,0,0,0,0,0,0}; 
	std::array <double,8>               contacts_pert      = {0,0,0,0,0,0,0,0}; 
	std::array <double,8>               frontflow_contacts = {0,0,0,0,0,0,0,0};
	double                              Emin               = 0; 

	double rng     = 0; // rng_uniform (0.0, 1.0); 
	double rsum    = 0; 
	int    e_idx   = 0; 


	double Esys    = *sysEnergy; 
	double Ei      = 0; 
	double Epert   = 0; 

	// std::cout << "solvent_indices = "; print((*solvation_shells));
	// loop over all solvent_indices
	for ( int i{0}; i < nflip; ++i ) {
		// get the flip index 
		rboltzmann = 0; 
		old_ori.push_back( (*LATTICE)[ solvation_shell_indices[i] ]->orientation );
		// find the neighboring interaction energies 
		Ei = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, solvation_shell_indices[i], x, y, z ); 

		for ( int j{0}; j < ntest; ++j ){

			(*LATTICE)[ solvation_shell_indices[i] ]->orientation = rng_uniform (0, 25); 
			orientations [j]    = (*LATTICE) [ solvation_shell_indices [i] ]->orientation; 
			Epert               = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_pert, solvation_shell_indices[i], x, y, z); // Epert = energy after perturbation
			energies [j]        = Esys - Ei + Epert; 
			contacts_store [j]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
			contacts_store [j]  = add_arrays ( &contacts_store[j], &contacts_pert ); 

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
		Esys         = energies[e_idx];
		contacts_sys = contacts_store[e_idx];
	}

	frontflow_energy   = energies[e_idx]; 
	frontflow_contacts = contacts_store[e_idx];  

	// figure out the backflow energy 


	for ( int i{0}; i < nflip; ++i ){

		Ei = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, solvation_shell_indices[i], x, y, z ); 
		rboltzmann = 0; 

		// first iteration 
		(*LATTICE) [ solvation_shell_indices[i] ]->orientation = old_ori[i]; 
		Epert               = NeighborEnergetics (LATTICE, InteractionMap, &contacts_pert, solvation_shell_indices[i], x, y, z); 
		energies[0]         = Esys - Ei + Epert; // CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z);
		contacts_store [0]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
		contacts_store [0]  = add_arrays ( &contacts_store[0], &contacts_pert ); 

		for ( int j{1}; j < ntest; ++j ){

			(*LATTICE)[ solvation_shell_indices[i] ]->orientation = rng_uniform(0, 25); 
			Epert       = NeighborEnergetics (LATTICE, InteractionMap, &contacts_pert, solvation_shell_indices[i], x, y, z); 
			energies[j] = Esys - Ei + Epert; 
			contacts_store [j]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
			contacts_store [j]  = add_arrays      ( &contacts_store[j], &contacts_pert );

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
		Esys = energies[0]; 
		contacts_sys = contacts_store[0];
	}

	// check the acceptance criterion 
	double rng_acc = rng_uniform (0.0, 1.0); 
	if ( rng_acc < std::exp (-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n  ) {
		// if accepted, return to the new orientations 
		for ( int j{0}; j < nflip; ++j ) {
			(*LATTICE)[ solvation_shell_indices[j] ]->orientation = new_ori[j]; 
		}

		*sysEnergy = frontflow_energy; 
		*contacts  = frontflow_contacts; 

	}
	else {
		*IMP_BOOL = false; 
	}
	
	return;

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void SolvationShellFlip_BIASED_remake2_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int x, int y, int z ){

    std::set <int> solvation_shell_set; 
    std::array <std::array<int,3>, 26> ne_list; 
    std::cout << "*E[0] = " << (*E)[0] << std::endl;
    // int dop = static_cast<int>((*Polymers)[0].chain.size() ); 

    // get the first solvation shell 
    // auto start = std::chrono::high_resolution_clock::now(); 
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
    int nflip = rng_uniform(1, static_cast<int>(solvation_shell_indices.size()/2 ) ); 			// number of sites to be flipped 

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));
    
    std::vector <int> old_ori; 					// vector to hold old orientations 
    std::vector <int> new_ori; 					// vector to hold new orientations 
    old_ori.reserve(nflip); 					// reserving... 
    new_ori.reserve(nflip); 					// reserving... 

    // energy store for boltzmann sampling 
    // instantiating a bunch of variables for the process 
    int                                 ntest              = 5; 
    std::array <double,5>               energies           = {0,0,0,0,0}; 
    std::array <double,5>               boltzmann          = {0,0,0,0,0};
    std::array <int,5>                  orientations       = {0,0,0,0,0}; 
    double                              rboltzmann         = 0;  
    double                              frontflow_energy   = 0; 
    double                              backflow_energy    = 1; 
    double                              prob_n_to_o        = 1; 
    double                              prob_o_to_n        = 1; 
    std::array<std::array<double,8>,5>  contacts_store     = {*contacts, *contacts, *contacts, *contacts, *contacts}; 
    std::array <double,8>               contacts_sys       = *contacts; 
    std::array <double,8>               contacts_i         = {0,0,0,0,0,0,0,0}; 
    std::array <double,8>               contacts_pert      = {0,0,0,0,0,0,0,0}; 
    std::array <double,8>               frontflow_contacts = {0,0,0,0,0,0,0,0};
    std::array <double,8>               backflow_contacts  = {0,0,0,0,0,0,0,0};
    std::array <double,8>			    c_contacts2        = {0,0,0,0,0,0,0,0}; 
    std::array <double,8>               contacts_test      = {0,0,0,0,0,0,0,0}; 
    double                              Emin               = 0; 

    double rng     = 0; // rng_uniform (0.0, 1.0); 
    double rsum    = 0; 
    int    e_idx   = 0; 

    
    
    double Esys    = *sysEnergy; 
    double Ei      = 0; 
    double Epert   = 0; 
    double E_test = 0; 

    // std::cout << "solvent_indices = "; print((*solvation_shells));
    // loop over all solvent_indices
    for ( int i{0}; i < nflip; ++i ) {
        // get the flip index 
    	rboltzmann = 0; 
    	old_ori.push_back( (*LATTICE)[ solvation_shell_indices[i] ]->orientation );
    	// find the neighboring interaction energies 
    	Ei = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, solvation_shell_indices[i], x, y, z ); 

    	for ( int j{0}; j < ntest; ++j ){
    	    
    	    (*LATTICE)[ solvation_shell_indices[i] ]->orientation = rng_uniform (0, 25); 
    	    orientations [j]    = (*LATTICE) [ solvation_shell_indices [i] ]->orientation; 
    	    Epert               = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_pert, solvation_shell_indices[i], x, y, z); // Epert = energy after perturbation
    	    energies [j]        = Esys - Ei + Epert; 
    	    contacts_store [j]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
    	    contacts_store [j]  = add_arrays ( &contacts_store[j], &contacts_pert ); 
    	    
    	    // DELETE THIS LATER 
    	    
    	    E_test         = CalculateEnergyRevamped ( Polymers, Cosolvent, LATTICE, InteractionMap, &contacts_test, x, y, z); 
    	    if ( E_test != energies[j] || contacts_store[j] != contacts_test  ){
    	    	std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
    			std::cout << "E_test = " << E_test << ", energies[j] = " << energies[j] << std::endl;
    			std::cout << "contacts_store[j] = "; print (contacts_store[j], ", "); std::cout << "contacts_test = "; print (contacts_test); 
    			exit(EXIT_FAILURE); 
    	    }
    	    
    	    // DELETE ABOVE 

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
		Esys         = energies[e_idx];
		contacts_sys = contacts_store[e_idx];
    }
    
    frontflow_energy   = energies[e_idx]; 
    frontflow_contacts = contacts_store[e_idx];  

    // figure out the backflow energy 
    

    for ( int i{0}; i < nflip; ++i ){

    	Ei = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, solvation_shell_indices[i], x, y, z ); 
    	rboltzmann = 0; 

    	// first iteration 
    	(*LATTICE) [ solvation_shell_indices[i] ]->orientation = old_ori[i]; 
    	Epert               = NeighborEnergetics (LATTICE, InteractionMap, &contacts_pert, solvation_shell_indices[i], x, y, z); 
    	energies[0]         = Esys - Ei + Epert; 
    	contacts_store [0]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
    	contacts_store [0]  = add_arrays ( &contacts_store[0], &contacts_pert ); 

    	// DELETE THIS LATER 
    	
	    E_test         = CalculateEnergyRevamped ( Polymers, Cosolvent, LATTICE, InteractionMap, &contacts_test, x, y, z); 
	    if ( E_test != energies[0] || contacts_store[0] != contacts_test  ){
	    	std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
			std::cout << "E_test = " << E_test << ", energies[j] = " << energies[0] << std::endl;
			std::cout << "contacts_store[j] = "; print (contacts_store[0], ", "); std::cout << "contacts_test = "; print (contacts_test); 
			exit(EXIT_FAILURE); 
	    }
	     
    	// DELETE ABOVE 


    	for ( int j{1}; j < ntest; ++j ){

    		(*LATTICE)[ solvation_shell_indices[i] ]->orientation = rng_uniform(0, 25); 
    		Epert       = NeighborEnergetics (LATTICE, InteractionMap, &contacts_pert, solvation_shell_indices[i], x, y, z); 
    		energies[j] = Esys - Ei + Epert; // CalculateEnergy_parallel (Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z);
    	    contacts_store [j]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
    	    contacts_store [j]  = add_arrays      ( &contacts_store[j], &contacts_pert );
    	    
    	    // DELETE THIS LATER 
    	    E_test         = CalculateEnergyRevamped ( Polymers, Cosolvent, LATTICE, InteractionMap, &contacts_test, x, y, z); 
    	    if ( E_test != energies[j] || contacts_store[j] != contacts_test  ){
    	    	std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
    			std::cout << "E_test = " << E_test << ", energies[j] = " << energies[j] << std::endl;
    			std::cout << "contacts_store[j] = "; print (contacts_store[j], ", "); std::cout << "contacts_test = "; print (contacts_test); 
    			exit(EXIT_FAILURE); 
    	    }
    	    // DELETE ABOVE 

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
		Esys = energies[0]; 
		contacts_sys = contacts_store[0];
    }

    // THIS CAN BE DELETED 

    backflow_energy   = energies[0];
    backflow_contacts = contacts_store[0]; 
    if ( backflow_energy != *sysEnergy || backflow_contacts != *contacts ){
    	std::cout << "Something is fucked. Energies do not match." << std::endl;
    	std::cout << "backflow_energy = " << backflow_energy << ", sysEnergy = " << *sysEnergy << std::endl;
    	std::cout << "c_contacts = "; print (backflow_contacts, ", "); std::cout << "contacts = "; print (*contacts); 
    	exit(EXIT_FAILURE); 
    }
    
    // ABOVE CAN BE DELETED 

    // check the acceptance criterion 

	double rng_acc = rng_uniform (0.0, 1.0); 
	if ( rng_acc < std::exp (-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n  ) {
		// if accepted, return to the new orientations 
		for ( int j{0}; j < nflip; ++j ) {
			(*LATTICE)[ solvation_shell_indices[j] ]->orientation = new_ori[j]; 
		}

		// THIS CAN BE DELETED 
		
		double en = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &c_contacts2, x, y, z); 
		
		if ( en != frontflow_energy || c_contacts2 != frontflow_contacts ){
    		std::cout << "Something is fucked. Energies do not match." << std::endl;
    		std::cout << "en = " << en << ", frontflow energy = " << frontflow_energy << std::endl;
    		std::cout << "c_contacts2 = "; print (c_contacts2, ", "); std::cout << "frontflow_contacts = "; print (frontflow_contacts); 
    		exit(EXIT_FAILURE); 
    	}
		
		// ABOVE CAN BE DELETED 
		*sysEnergy = frontflow_energy; 
		*contacts  = frontflow_contacts; 

	}
	else {
		*IMP_BOOL = false; 
	}
	
	// CheckStructures (Polymers, Cosolvent, LATTICE, x, y, z);
	return;
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of SecondSolvationShellFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of PolymerFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void PolymerFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z ){

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
	std::array<double,8>               contacts_sys		  = *contacts; 
	std::array<double,8>			   contacts_i         = {0,0,0,0,0,0,0,0};
	std::array<double,8>			   contacts_pert      = {0,0,0,0,0,0,0,0};
	std::array<double,8>			   frontflow_contacts = {0,0,0,0,0,0,0,0};
	double                             rboltzmann         = 0;
	double                             frontflow_energy   = 0; 
	double                             prob_o_to_n        = 1; 
	double 							   prob_n_to_o        = 1; 
	double                             Emin               = 0; 


	double rng   = 0; 
	double rsum  = 0;
	int    e_idx = 0; 

	double Esys            = *sysEnergy; 
	double Ei              = 0; 
	double Epert           = 0; 
	int    m_lattice_idx   = -1; 

	// loop over all solvent_indices 
	for ( int i{0}; i < nflip; ++i ){

		rboltzmann        = 0; 
		old_ori.push_back ( (*Polymers)[index].chain[ polymer_indices[i] ]->orientation ); 
		m_lattice_idx     = lattice_index((*Polymers)[index].chain[polymer_indices[i]]->coords, y, z);
		Ei                = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, m_lattice_idx, x, y, z); 

		for ( int j{0}; j < ntest; ++j ){

			(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = rng_uniform (0, 25); 
			orientations[j]   = (*Polymers)[index].chain[ polymer_indices[i] ]->orientation;
			Epert             = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_pert, m_lattice_idx, x, y, z );
			energies[j]       = Esys - Ei + Epert; 
			contacts_store[j] = subtract_arrays (&contacts_sys, &contacts_i); 
			contacts_store[j] = add_arrays      (&contacts_store[j], &contacts_pert); 

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
		Esys         = energies[e_idx]; 
		contacts_sys = contacts_store [e_idx]; 
	}

	frontflow_energy        = energies[e_idx]; 
	frontflow_contacts      = contacts_store[e_idx]; 

	// std::cout << "Starting backflow... " << std::endl;

	// figure out the backflow energy 
	// double backflow_energy = 1;

	for ( int i{0}; i < nflip; ++i ){

		rboltzmann = 0; 
		m_lattice_idx = lattice_index((*Polymers)[index].chain[polymer_indices[i]]->coords, y, z);
		Ei = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, m_lattice_idx, x, y, z); 

		(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = old_ori[i]; 
		Epert       = NeighborEnergetics (LATTICE, InteractionMap, &contacts_pert, m_lattice_idx, x, y, z); 
		energies[0] = Esys - Ei + Epert; 
		contacts_store [0]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
		contacts_store [0]  = add_arrays      ( &contacts_store[0], &contacts_pert );  

		for (int j{1}; j<ntest; ++j){
			(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = rng_uniform (0, 25); 
			Epert       = NeighborEnergetics  (LATTICE, InteractionMap, &contacts_pert, m_lattice_idx, x, y, z); 
			energies[j] = Esys - Ei + Epert; 
			contacts_store [j]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
			contacts_store [j]  = add_arrays      ( &contacts_store[j], &contacts_pert ); 

		}

		Emin = *std::min_element ( energies.begin(), energies.end() ); 

		for (int k{0}; k < ntest; ++k){
			boltzmann [k] = std::exp (-1/temperature * (energies[k] - Emin ) );
			rboltzmann   += boltzmann[k]; 
		}
		prob_n_to_o      *= boltzmann[0]/rboltzmann; 

		// make the jump to the old state 
		(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = old_ori[i]; 
		Esys = energies[0]; 
		contacts_sys = contacts_store[0]; 

	}

	// check the acceptance criterion 

	double rng_acc = rng_uniform (0.0, 1.0); 

	if ( rng_acc < std::exp (-1/temperature * (frontflow_energy - *sysEnergy ) ) ){
		// if accepted, return to the new orientations 
		for (int j{0}; j < nflip; ++j){
			(*Polymers)[index].chain[polymer_indices[j]]->orientation = new_ori[j];
		}

		*sysEnergy = frontflow_energy; 
		*contacts  = frontflow_contacts;  

	}
	else {
		*IMP_BOOL = false; 
	}

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void PolymerFlip_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* E, std::array<double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, double temperature, int index, int x, int y, int z ) {

	int deg_poly = (*Polymers)[0].deg_poly; 
	std::vector <int> polymer_indices (deg_poly);
	std::iota (polymer_indices.begin(), polymer_indices.end(), 0);

	std::cout << "*E[0] = " << (*E)[0] << std::endl;

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
	std::array<double,8>               contacts_sys		  = *contacts; 
	std::array<double,8>			   contacts_i         = {0,0,0,0,0,0,0,0};
	std::array<double,8>			   contacts_pert      = {0,0,0,0,0,0,0,0};
	std::array<double,8>			   frontflow_contacts = {0,0,0,0,0,0,0,0};
	std::array<double,8>			   backflow_contacts  = {0,0,0,0,0,0,0,0};
	std::array<double,8>			   contacts_test      = {0,0,0,0,0,0,0,0};
	double                             rboltzmann         = 0;
	double                             frontflow_energy   = 0; 
	double                             backflow_energy    = 0; 
	double                             prob_o_to_n        = 1; 
	double 							   prob_n_to_o        = 1; 
	double                             Emin               = 0; 


	double rng   = 0; 
	double rsum  = 0;
	int    e_idx = 0; 

	double Esys    = *sysEnergy; 
	double Ei      = 0; 
	double Epert   = 0; 
	double E_test  = 0; 
	int    m_lattice_idx   = -1; 
	// loop over all solvent_indices 
	for ( int i{0}; i < nflip; ++i ){

		rboltzmann = 0; 
		old_ori.push_back ( (*Polymers)[index].chain[ polymer_indices[i] ]->orientation ); 
		m_lattice_idx = lattice_index((*Polymers)[index].chain[polymer_indices[i]]->coords, y, z);
		Ei    = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, m_lattice_idx, x, y, z); 

		for ( int j{0}; j < ntest; ++j ){

			(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = rng_uniform (0, 25); 
			orientations[j]   = (*Polymers)[index].chain[ polymer_indices[i] ]->orientation;
			Epert             = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_pert, m_lattice_idx, x, y, z );
			energies[j]       = Esys - Ei + Epert; 
			contacts_store[j] = subtract_arrays (&contacts_sys, &contacts_i); 
			contacts_store[j] = add_arrays      (&contacts_store[j], &contacts_pert); 

			// DELETE THIS LATER 
		
			E_test         = CalculateEnergyRevamped ( Polymers, Cosolvent, LATTICE, InteractionMap, &contacts_test, x, y, z); 
			if ( E_test != energies[j] || contacts_store[j] != contacts_test  ){
				std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
				std::cout << "E_test = " << E_test << ", energies[j] = " << energies[j] << std::endl;
				std::cout << "contacts_store[j] = "; print (contacts_store[j], ", "); std::cout << "contacts_test = "; print (contacts_test); 
				exit(EXIT_FAILURE); 
			}
		 
			// DELETE ABOVE 

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
		Esys         = energies[e_idx]; 
		contacts_sys = contacts_store [e_idx]; 
	}

	frontflow_energy        = energies[e_idx]; 
	frontflow_contacts      = contacts_store[e_idx]; 

	// std::cout << "Starting backflow... " << std::endl;

	// figure out the backflow energy 
	// double backflow_energy = 1;

	for ( int i{0}; i < nflip; ++i ){

		rboltzmann = 0; 
		m_lattice_idx = lattice_index((*Polymers)[index].chain[polymer_indices[i]]->coords, y, z);
		Ei = NeighborEnergetics ( LATTICE, InteractionMap, &contacts_i, m_lattice_idx, x, y, z); 

		(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = old_ori[i]; 
		Epert               = NeighborEnergetics  (LATTICE, InteractionMap, &contacts_pert, m_lattice_idx, x, y, z); 
		energies       [0]  = Esys - Ei + Epert; 
		contacts_store [0]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
		contacts_store [0]  = add_arrays ( &contacts_store[0], &contacts_pert );  

		// DELETE THIS LATER 
		
		E_test         = CalculateEnergyRevamped ( Polymers, Cosolvent, LATTICE, InteractionMap, &contacts_test, x, y, z); 
		if ( E_test != energies[0] || contacts_store[0] != contacts_test  ){
			std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
			std::cout << "E_test = " << E_test << ", energies[j] = " << energies[0] << std::endl;
			std::cout << "contacts_store[j] = "; print (contacts_store[0], ", "); std::cout << "contacts_test = "; print (contacts_test); 
			exit(EXIT_FAILURE); 
		}
	
	// DELETE ABOVE

		for (int j{1}; j<ntest; ++j){
			(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = rng_uniform (0, 25); 
			Epert       = NeighborEnergetics  (LATTICE, InteractionMap, &contacts_pert, m_lattice_idx, x, y, z); 
			energies[j] = Esys - Ei + Epert; 
			contacts_store [j]  = subtract_arrays ( &contacts_sys, &contacts_i ); 
			contacts_store [j]  = add_arrays      ( &contacts_store[j], &contacts_pert ); 
			// DELETE THIS LATER 
			
			E_test         = CalculateEnergyRevamped ( Polymers, Cosolvent, LATTICE, InteractionMap, &contacts_test, x, y, z); 
			if ( E_test != energies[j] || contacts_store[j] != contacts_test  ){
				std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
				std::cout << "E_test = " << E_test << ", energies[j] = " << energies[j] << std::endl;
				std::cout << "contacts_store[j] = "; print (contacts_store[j], ", "); std::cout << "contacts_test = "; print (contacts_test); 
				exit(EXIT_FAILURE); 
			}
			// DELETE ABOVE 
		}

		Emin = *std::min_element ( energies.begin(), energies.end() ); 

		for (int k{0}; k < ntest; ++k){
			boltzmann [k] = std::exp (-1/temperature * (energies[k] - Emin ) );
			rboltzmann   += boltzmann[k]; 
		}
		prob_n_to_o      *= boltzmann[0]/rboltzmann; 

		// make the jump to the old state 
		(*Polymers)[index].chain[ polymer_indices[i] ]->orientation = old_ori[i]; 
		Esys = energies[0]; 
		contacts_sys = contacts_store[0]; 

	}

	// std::cout << "Done with backflow! " << std::endl;
	// DELETE LATER
	backflow_energy   = energies[0]; 
	backflow_contacts = contacts_store[0]; 

	if ( backflow_energy != *sysEnergy || backflow_contacts != *contacts ){
    	std::cout << "Something is fucked. Energies do not match." << std::endl;
    	std::cout << "backflow_energy = " << backflow_energy << ", sysEnergy = " << *sysEnergy << std::endl;
    	std::cout << "backflow_contacts = "; print (backflow_contacts, ", "); std::cout << "contacts = "; print (*contacts); 
    	exit(EXIT_FAILURE); 
    }
    // DELETE ABOVE LATER

    // check the acceptance criterion 

    double rng_acc = rng_uniform (0.0, 1.0); 

    if ( rng_acc < std::exp (-1/temperature * (frontflow_energy - *sysEnergy ) ) ){
    	// if accepted, return to the new orientations 
    	for (int j{0}; j < nflip; ++j){
    		(*Polymers)[index].chain[polymer_indices[j]]->orientation = new_ori[j];
    	}

    	// DELETE LATER
    	
		double en = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &contacts_test, x, y, z); 
		
		if ( en != frontflow_energy || contacts_test != frontflow_contacts ){
    		std::cout << "Something is fucked. Energies do not match." << std::endl;
    		std::cout << "en = " << en << ", frontflow energy = " << frontflow_energy << std::endl;
    		std::cout << "contacts_test = "; print (contacts_test, ", "); std::cout << "frontflow_contacts = "; print (frontflow_contacts); 
    		exit(EXIT_FAILURE); 
    	}
    	 
    	// DELETE LATER

		*sysEnergy = frontflow_energy; 
		*contacts  = frontflow_contacts;  

    }
    else {
    	*IMP_BOOL = false; 
    }

	return; 
}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of PolymerFlip_BIASED
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

void ChainRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, \
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
			old_cut.push_back ((*Polymers).at(p_index).chain.at(i)->coords) ;
		}

		// regrow the polymer frontwards
		// std::cout << "Regrowth time!" << std::endl;
		HeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
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
		BackFlowFromHeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, InteractionMap, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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
		TailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 

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

		BackFlowFromTailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, InteractionMap, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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

void ChainRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int p_index, int x, int y, int z) {


	int deg_poly = (*Polymers)[p_index].deg_poly; 
	std::cout << "deg_poly = " << deg_poly << std::endl;
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
		
		std::cout << "Head regrowth... " << std::endl;
		// std::cout << "OLD ENERGY = " << *sysEnergy << std::endl << std::endl;
		// get old_cut 
		
		for (int i {m_index+1}; i<deg_poly; ++i){
			old_cut.push_back ((*Polymers).at(p_index).chain.at(i)->coords) ;
		}

		// regrow the polymer frontwards
		// std::cout << "Regrowth time!" << std::endl;
		
		HeadRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		
		// std::cout << "Worked!" << std::endl;
		
		for (int i {m_index+1}; i<deg_poly; ++i){
			new_cut.push_back ((*Polymers).at(p_index).chain.at(i)->coords) ;
		}

		if ( old_cut == new_cut ){
			return; 
		}

		if ( !(*IMP_BOOL) ){
			acceptance_after_head_regrowth ( LATTICE, &new_cut, &old_cut, y, z );
			return; 
		}

		c2_contacts = c1_contacts; 
		// std::cout << "c2_contacts = "; print (c2_contacts);
		backflow_energy = frontflow_energy; 

		std::cout << "-------------------" << std::endl;
		std::cout << "Backflow time!" << std::endl;
		BackFlowFromHeadRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, &old_cut, InteractionMap, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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
		TailRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 

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

		BackFlowFromTailRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, &old_cut, InteractionMap, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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

void HeadRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	std::array <double,8>               current_contacts = *contacts; 

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>               energies; 
	std::array <std::array<double,8>,5> contacts_store; 
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
	double Epair      = 0;
	double Epair_n    = 0;

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 

	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 
	int c_idx               = 0;
	int c_idx_n             = 0;

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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				block_counter += 1; 
			}
			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx, x, y, z);

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the new energies
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z);

				energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
				contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [ idx_counter ][ c_idx ] += 1;
				contacts_store [ idx_counter ][ c_idx_n ] -= 1; 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];
			}
			
		}
		else {
			// std::cout << "solvent swap time." << std::endl;
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx, x, y, z);

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

			// get the new energies
			Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z);

			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 	
			contacts_store [ idx_counter ][ c_idx ]   += 1;
			contacts_store [ idx_counter ][ c_idx_n ] -= 1;

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
	HeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void HeadRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z) {


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	std::array <double,8>               current_contacts = *contacts; 
	std::array <double,8>				copy_contacts    = *contacts; // this is only for the checks with brute force calculation

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>               energies; 
	std::array <std::array<double,8>,5> contacts_store; 
	std::array <double,5>               boltzmann; 
	std::array <int,3>                  loc_m = (*Polymers)[p_index].chain[m_index+1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 
	
	// doubles for energy transfer 
	int    c_idx      = 0;
	int    c_idx_n    = 0;

	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	double Etest      = 0;
	double Epair      = 0;
	double Epair_n    = 0; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				block_counter += 1; 
			}
			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx, x, y, z);

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the new energies
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z);

				energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
				contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [ idx_counter ][ c_idx   ] += 1;
				contacts_store [ idx_counter ][ c_idx_n ] -= 1;
				
				// get the energy
				// delete later 
				
				Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
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
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx, x, y, z);

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

			// get the new energies
			Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index (loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z);

			// std::cout << "Es_n = " << Es_n << std::endl; 
			// std::cout << "Em_n = " << Em_n << std::endl;

			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 	
			contacts_store [ idx_counter ][ c_idx   ] += 1;
			contacts_store [ idx_counter ][ c_idx_n ] -= 1;			

			// get the energy
			// delete later 
			
			Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
			if ( Etest != energies [idx_counter] || contacts_store[idx_counter] != copy_contacts ) {
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
	HeadRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of HeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of BackFlowFromHeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	double                               rboltzmann       = 0; // running sum for boltzmann weights 
	std::array <double,5>                energies; 
	std::array <double,5>                boltzmann; 
	std::array <double,8>                current_contacts = *contacts; 
	std::array <std::array<double,8>,5>  contacts_store; 

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
	int    c_idx   = 0;
	int    c_idx_n = 0;
	double Em      = 0;
	double Es      = 0; 
	double Em_n    = 0; 
	double Es_n    = 0; 
	double Esys    = *backflow_energy; 
	double Epair   = 0;
	double Epair_n = 0;

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};


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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
			}
			else {

				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap 
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
				
				// get the new energies
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z); 

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store[idx_counter][c_idx  ] += 1;
				contacts_store[idx_counter][c_idx_n] -= 1;

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
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index+1];

			// get the energy 
			Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z); 

			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 	
			contacts_store [idx_counter][c_idx  ] += 1;
			contacts_store [idx_counter][c_idx_n] -= 1;


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

	
	*prob_n_to_o      = (*prob_n_to_o)*boltzmann[0]/rboltzmann;
	*backflow_energy  = energies[0]; 
	*contacts         = contacts_store [0];

	// do the swap again
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords           = ne_list[0];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];

	// else { do nothing, MAINTAIN, because suggested index is in the maintain index vector }
	BackFlowFromHeadRegrowth_BIASED (Polymers, Cosolvent, LATTICE, old_cut, InteractionMap, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector <std::array<int,3>>* old_cut, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z) {


	if (m_index == deg_poly-1){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	
	std::array <double,5>               energies; 
	std::array <double,5>               boltzmann; 
	std::array <double,8>               copy_contacts    = *contacts; 
	std::array <double,8>               current_contacts = *contacts; 
	std::array <std::array<double,8>,5> contacts_store; 
	double                              rboltzmann       = 0; // running sum for boltzmann weights 

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
	int    c_idx   = 0;
	int    c_idx_n = 0;
	double Em      = 0;
	double Es      = 0; 
	double Em_n    = 0; 
	double Es_n    = 0; 
	double Esys    = *backflow_energy; 
	double Etest   = 0;
	double Epair   = 0;
	double Epair_n = 0;

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};


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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
			}
			else {

				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z); 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap 
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
				
				// get the new energies
				Es_n    = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n    = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z); 

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store[idx_counter][c_idx  ] += 1;
				contacts_store[idx_counter][c_idx_n] -= 1;
				
				// delete later 
				
				Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
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
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z); 

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index+1];

			// get the energy 
			Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z); 

			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 	
			contacts_store[idx_counter][c_idx  ] += 1;
			contacts_store[idx_counter][c_idx_n] -= 1;

			// delete later
			Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
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
	BackFlowFromHeadRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, old_cut, InteractionMap, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of BackflowHeadRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of TailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void TailRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	std::array <double,8>               current_contacts = *contacts; 

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>               energies; 
	std::array <std::array<double,8>,5> contacts_store; 
	std::array <double,5>               boltzmann; 
	std::array <int,3>                  loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// doubles for energy transfer 
	int    c_idx      = 0;
	int    c_idx_n    = 0; 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	double Epair      = 0;
	double Epair_n    = 0;
	
	// contacts for energy transfer 
	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				block_counter += 1; 
			}
			else {
				// get the initial neighboring energies 
				Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z); 
				Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the energy
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

				// do the calc
				energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n + Em_n-Epair_n); 
				contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [ idx_counter ][c_idx  ] += 1;
				contacts_store [ idx_counter ][c_idx_n] -= 1;

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
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the new energies
			Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

			// run the computation 
			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
			contacts_store [ idx_counter ][c_idx  ] += 1;
			contacts_store [ idx_counter ][c_idx_n] -= 1;



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
	TailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void TailRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	std::array <double,8>               current_contacts = *contacts; 
	std::array <double,8>				copy_contacts    = *contacts; // this is only for the checks with brute force calculation

	// current contacts will be perturbed, and eventually merged with *contacts. 

	std::array <double,5>               energies; 
	std::array <std::array<double,8>,5> contacts_store; 
	std::array <double,5>               boltzmann; 
	std::array <int,3>                  loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 

	// generate possible locations to jump to. 
	std::array <std::array<int,3>, 26> ne_list = obtain_ne_list ((*Polymers)[p_index].chain[m_index]->coords, x, y, z); 

	// randomly select five of them 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

	// doubles for energy transfer 
	int    c_idx      = 0;
	int    c_idx_n    = 0;
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0; 
	double Esys       = *frontflow_energy; 
	double Etest      = 0;
	double Epair      = 0;
	double Epair_n    = 0; 
	
	// contacts for energy transfer 
	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				block_counter += 1; 
			}
			else {
				// get the initial neighboring energies 
				Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z); 
				Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the energy
				Es_n    = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n    = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

				energies [idx_counter] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
				contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [ idx_counter ][c_idx]   += 1;
				contacts_store [ idx_counter ][c_idx_n] -= 1; 
				
				// get the energy 
				// delete below
				
				Etest = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
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

			Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the new energies
			Es_n    = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n    = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

			// run the computation 
			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 	
			contacts_store [ idx_counter ][c_idx]   += 1;
			contacts_store [ idx_counter ][c_idx_n] -= 1;

			// delete this later 
			
			Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
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
	
	// if { e_idx is not in the maintain index vector, perform the swap}
		
	// do the swap again 	
	(*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords           = ne_list[e_idx];
	
	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
	

	// else { do nothing and maintain structure, because the suggested index is in the the maintain index vector }
	TailRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             Start of BackFlowFromTailRegrowth_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromTailRegrowth_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::vector<std::array<int,3>>* old_cut, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, \
	int p_index, int m_index, int recursion_depth, int x, int y, int z){

	if (m_index == 0){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	
	// generate an array for energies 
	
	std::array <double,5>               energies; 
	std::array <double,5>               boltzmann; 
	std::array <double,8>               current_contacts = *contacts; 
	std::array <std::array<double,8>,5> contacts_store; 
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
	int    c_idx   = 0;
	int    c_idx_n = 0;
	double Em      = 0;
	double Es      = 0; 
	double Em_n    = 0; 
	double Es_n    = 0; 
	double Esys    = *backflow_energy; 
	double Epair   = 0;
	double Epair_n = 0;

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};


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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
			}

			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 				
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the new energies	
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store[idx_counter][c_idx  ] += 1;
				contacts_store[idx_counter][c_idx_n] -= 1;				

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
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords                 = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ]                = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the energy 
			Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 	
			contacts_store[idx_counter][c_idx  ] += 1;
			contacts_store[idx_counter][c_idx_n] -= 1;				

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

	BackFlowFromTailRegrowth_BIASED (Polymers, Cosolvent, LATTICE, old_cut, InteractionMap, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromTailRegrowth_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, std::vector<std::array<int,3>>* old_cut, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array<double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, double* backflow_energy, double temperature, int deg_poly, \
	int p_index, int m_index, int recursion_depth, int x, int y, int z){

	if (m_index == 0) {
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	
	// generate an array for energies 
	
	std::array <double,5>               energies; 
	std::array <double,5>               boltzmann; 
	std::array <double,8>               copy_contacts    = *contacts; 
	std::array <double,8>               current_contacts = *contacts; 
	std::array <std::array<double,8>,5> contacts_store; 
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
	int    c_idx   = 0; 
	int    c_idx_n = 0;
	double Epair   = 0;
	double Epair_n = 0;
	double Em      = 0;
	double Es      = 0; 
	double Em_n    = 0; 
	double Es_n    = 0; 
	double Esys    = *backflow_energy; 
	double Etest   = 0;

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};


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
				contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
			}

			else {
				// std::cout << "monomer swap time." << std::endl;
				// get the initial neighboring energies. 
				Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 	
				Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

				// get the new energies	
				Es_n    = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n    = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

				// get the energy
				energies      [idx_counter] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); // CalculateEnergy (Polymers, Cosolvent, LATTICE, E, contacts, x, y, z); 
				contacts_store[idx_counter] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store[idx_counter][c_idx  ] += 1;
				contacts_store[idx_counter][c_idx_n] -= 1;
				// delete later 
				
				Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
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
			Es    = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em    = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx, x, y, z);

			// prep the swap 
			(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
			(*Polymers)[p_index].chain[m_index-1]->coords                 = ne_list[idx_counter];
			
			// perform the swap (since coords were changed, this swap works)
			(*LATTICE)[ lattice_index (loc_m, y, z) ]                = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
			(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];

			// get the energy 
			Es_n    = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em_n    = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index (loc_m, y, z)], (*LATTICE)[lattice_index (ne_list[idx_counter], y, z)], InteractionMap, &c_idx_n, x, y, z);

			energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
			contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 	
			contacts_store [ idx_counter ][c_idx]   += 1;
			contacts_store [ idx_counter ][c_idx_n] -= 1; 

			// delete later
			
			Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
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

	BackFlowFromTailRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, old_cut, InteractionMap, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

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
	double* sysEnergy, double temperature, int p_index, int x, int y, int z) {

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
		std::cout << "Head regrowth time..." << std::endl;
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

		frontflow_energy = CalculateEnergy (Polymers, Cosolvent, LATTICE, E, &c1_contacts, x, y, z); 

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

		std::cout << "Tail regrowth... " << std::endl;
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

		frontflow_energy = CalculateEnergy (Polymers, Cosolvent, LATTICE, E, &c1_contacts, x, y, z); 
		
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
	bool* IMP_BOOL, int deg_poly, int p_index, int m_index, int x, int y, int z){


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
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	
	std::array <double,8>                current_contacts = *contacts; 
	std::array <double,25>               energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25>               boltzmann; 
	std::array <int,25>                  orientations; 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 
	int original_ori         = (*Polymers)[p_index].chain[m_index-1]->orientation; 

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
	double Epair      = 0;
	double Epair_n    = 0;
	double Esys       = *frontflow_energy; 


	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int c_idx               = 0 ; 
	int c_idx_n             = 0 ; 
	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 
	

	while ( idx_counter < 5 ){

		// std::cout << "ptype = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;

		if ( ne_list[idx_counter] == loc_m ){

			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index (loc_m, y, z), x, y, z); 
			for ( int i{0}; i<5; ++i){

				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
				orientations[5*idx_counter+i] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation; 
			
				Em_n = NeighborEnergetics ( LATTICE, InteractionMap, &cm_n, lattice_index (loc_m, y, z), x, y, z);
				energies       [5*idx_counter+i] = Esys - Em + Em_n; 
				contacts_store [5*idx_counter+i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 

				// go back to the original state
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 

			}
			 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ) {

			// check which region of the polymer the swap is taking place 
			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ){
					self_swap_idx = u; 
					break; 
				}				
			}

			// if it is with other tail units, do the swap 
			// if not, discourage it 
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z); 
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 


			for ( int i{0}; i < 5; ++i ){
				
				if ( self_swap_idx > m_index ){

					orientations[5*idx_counter+i] = original_ori; 
					energies[5*idx_counter+i] = 1e+08; 
					contacts_store[5*idx_counter+i] = {-1, -1, -1, -1, -1, -1, -1, -1}; 
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
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
					Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 

					energies       [5*idx_counter+i] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
					

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
			// std::cout << "performing a solvent swap." << std::endl;
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 

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
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
				Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 
				

				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair ; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
				

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

	TailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void TailRegrowthPlusOrientationFlip_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return;
	}

	// std::cout << "m_index = " << m_index << std::endl; 
	// generate an array for energies 
	
	std::array <double,8>				 current_contacts = *contacts; 
	std::array <double,8>				 copy_contacts    = *contacts; 
	std::array <double,25> 				 energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25> 				 boltzmann; 
	std::array <int,25> 				 orientations; 

	std::array <int,3> loc_m = (*Polymers)[p_index].chain[m_index-1]->coords; 
	int original_ori         = (*Polymers)[p_index].chain[m_index-1]->orientation; 

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
	double Epair      = 0;  
	double Epair_n    = 0; 
	double Etest      = 0;
	double Esys       = *frontflow_energy; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int c_idx               = 0 ; 
	int c_idx_n             = 0 ; 
	int block_counter       = 0 ; 
	int idx_counter         = 0 ; 
	int self_swap_idx       = -1; 

	std::array <int,3> difference = {0,0,0};
	// std::cout << "current_contacts = "; print(current_contacts);
	while ( idx_counter < 5 ){

		// std::cout << "ptype = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;

		if ( ne_list[idx_counter] == loc_m ){

			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index (loc_m, y, z), x, y, z); 
			for ( int i{0}; i<5; ++i){

				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
				orientations[5*idx_counter+i] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation; 
			
				Em_n = NeighborEnergetics ( LATTICE, InteractionMap, &cm_n, lattice_index (loc_m, y, z), x, y, z);
				energies       [5*idx_counter+i] = Esys - Em + Em_n; 
				contacts_store [5*idx_counter+i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 
				// delete later
				Etest = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z);
			    if ( Etest != energies[5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts  ){
			    	std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
					std::cout << "Etest = " << Etest << ", energies[j] = " << energies[5*idx_counter+i] << std::endl;
					std::cout << "contacts_store[j] = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "contacts_test = "; print (copy_contacts); 
					exit(EXIT_FAILURE); 
			    }				
				// delete above later 
				// go back to the original state
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 

			}
			 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ) {

			// check which region of the polymer the swap is taking place 

			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ){
					self_swap_idx = u; 
					break; 
				}				
			}

			// if it is with other tail units, do the swap 
			// if not, discourage it 
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z); 
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
			difference = subtract_arrays (&ne_list[idx_counter], &loc_m); 


			for ( int i{0}; i < 5; ++i ){
				
				if ( self_swap_idx > m_index ){
					orientations[5*idx_counter+i] = original_ori; 
					energies[5*idx_counter+i] = 1e+08; 
					contacts_store[5*idx_counter+i] = {-1, -1, -1, -1, -1, -1, -1, -1}; 
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
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
					Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 

					energies       [5*idx_counter+i] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 	

					// delete later 
					Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
					if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
						std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
						std::cout << "Etest (real) = " << Etest << ", energies [" << 5*idx_counter+i <<"] (calc)= " << energies[5*idx_counter+i] << "." << std::endl;
						std::cout << "contacts_store (calc) = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts (real) = "; print(copy_contacts);
						std::cout << "Shit's fucked." << std::endl;
						exit(EXIT_FAILURE);
					}
					// delete above later 

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
			// std::cout << "performing a solvent swap." << std::endl;
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 

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
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
				Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 
		
				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair ; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
				
				// delete later 
				Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
				if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << 5*idx_counter+i <<"] = " << energies[5*idx_counter+i] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				// delete above later 

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

	(*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords             = ne_list[e_idx/5];
	(*Polymers)[p_index].chain[m_index-1]->orientation        = orientations[e_idx];
	
	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx/5], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
	
	// else { do nothing and maintain structure, because the suggested index is in the the maintain index vector }
	// std::cout << "m_index = " << m_index << std::endl;
	// (*Polymers)[0].printChainCoords(); 

	TailRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index-1, x, y, z);

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of TailRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void BackFlowFromTailRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::vector <int>* old_ori, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	std::array <double,8>                current_contacts = *contacts; 
	
	std::array <double,25>               energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25>               boltzmann; 

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

	// i now have a vector which has the back peddling step at position index 0 
	// these are doubles for energies and other stuff 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0;
	double Epair      = 0; 
	double Epair_n    = 0;  
	double Esys       = *backflow_energy; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int  c_idx                     = 0 ;
	int  c_idx_n			       = 0 ;
	int  idx_counter               = 0 ;
	int  self_swap_idx 		       = -1;

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);
		if ( ne_list[idx_counter] == loc_m ){

			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z);
			for (int i{0}; i < 5; ++i ){
				// std::cout << "monomer self-swap: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index-1]->orientation = test_ori[5*idx_counter+i];
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
				energies[5*idx_counter+i]       = Esys - Em +Em_n ;
				contacts_store[5*idx_counter+i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 

				(*Polymers)[p_index].chain[m_index-1]->orientation = original_ori; 
			}
			
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			for ( int u{0}; u<deg_poly; ++u) {
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list [idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}
			}

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
	
			for (int i{0}; i<5; ++i){
				if ( self_swap_idx > m_index-1 ) {
					
					(*Polymers)[p_index].chain[m_index-1]->orientation = original_ori;
					energies       [5*idx_counter+i]                    = 1e+08; 
					contacts_store [5*idx_counter+i]                    = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				
				}

				else {

					// std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
					// prep the swap 
					(*Polymers)[p_index].chain[m_index-1]->orientation     = test_ori[5*idx_counter+i];
					(*LATTICE) [lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
					
					// perform the swap 
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
					
					// get the energy
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
					Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 


					energies [5*idx_counter+i] = Esys - (Es + Em - Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
					

					// revert back to original structure 
					(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];

					// reset orientation 
					(*Polymers)[p_index].chain[m_index-1]->orientation       = original_ori; 
				}
			}
		}

		else {

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 

			for (int i{0}; i<5; ++i){
				// prep the swap 
				// std::cout << "solventswap: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index-1]->orientation     = test_ori[5*idx_counter+i];
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords  = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index-1];
				

				// get the energy 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
				Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 


				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
				

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];

				// reset orientation 
				(*Polymers)[p_index].chain[m_index-1]->orientation       = original_ori; 
			}
		}

		idx_counter += 1;

	}
	
	double Emin = *std::min_element ( energies.begin(), energies.end() ); 
	
	// std::cout << "Emin = " << Emin << std::endl; 

	for (int i{0}; i<25; ++i){
		boltzmann[i] = std::exp  ( -1/temperature * ( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann [i]; 
	}

	// std::cout << "Boltzmann weights are: "; print(boltzmann); 
	// std::cout << "normalization = " << rboltzmann << std::endl;

	*prob_n_to_o      = (*prob_n_to_o)*boltzmann[0]/rboltzmann;
	*backflow_energy  = energies[0]; 
	*contacts         = contacts_store [0];

	// do the swap again
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords           = ne_list[0];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ]  = (*Polymers)[p_index].chain[m_index-1];
	(*Polymers)[p_index].chain[m_index-1]->orientation = test_ori[0]; 

	// else { do nothing, MAINTAIN, because suggested index is in the maintain index vector }
	BackFlowFromTailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, old_cut, old_ori, InteractionMap, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromTailRegrowthPlusOrientationFlip_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::vector <int>* old_ori, \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == 0){
		*IMP_BOOL = true; 
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 
	
	std::array <double,8>                copy_contacts    = *contacts; 
	std::array <double,8>                current_contacts = *contacts; 
	std::array <double,25>               energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25>               boltzmann; 

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

	// i now have a vector which has the back peddling step at position index 0 
	// these are doubles for energies and other stuff 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0;
	double Epair      = 0; 
	double Epair_n    = 0;  
	double Etest      = 0;
	double Esys       = *backflow_energy; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int  c_idx                     = 0 ;
	int  c_idx_n			       = 0 ;
	int  idx_counter               = 0 ;
	int  self_swap_idx 		       = -1;

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);
		if ( ne_list[idx_counter] == loc_m ){

			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z);
			for (int i{0}; i < 5; ++i ){
				// std::cout << "monomer self-swap: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index-1]->orientation = test_ori[5*idx_counter+i];
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
				energies[5*idx_counter+i]       = Esys - Em +Em_n ;
				contacts_store[5*idx_counter+i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 

				// delete later
				Etest = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
			    if ( Etest != energies[5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
			    	std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
					std::cout << "Etest = " << Etest << ", energies[j] = " << energies[5*idx_counter+i] << std::endl;
					std::cout << "contacts_store[j] = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "contacts_test = "; print (copy_contacts); 
					exit(EXIT_FAILURE); 
			    }
				// delete above later 
				(*Polymers)[p_index].chain[m_index-1]->orientation = original_ori; 
			}
			
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			for ( int u{0}; u<deg_poly; ++u) {
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list [idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}
			}

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
	
			for (int i{0}; i<5; ++i){
				if ( self_swap_idx > m_index-1 ){
					
					// std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
					// std::cout << "self_swap_idx = " << self_swap_idx << std::endl;				
					// std::cout << "Selected position is "; print (ne_list[idx_counter]);	
					(*Polymers)[p_index].chain[m_index-1]->orientation = original_ori ;
					energies [5*idx_counter+i]      = 1e+08; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				
				}

				else {

					// std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
					// prep the swap 
					(*Polymers)[p_index].chain[m_index-1]->orientation     = test_ori[5*idx_counter+i];
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index-1]->coords       = ne_list[idx_counter];
					
					// perform the swap 
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index-1];
					
					// get the energy
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
					Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 


					energies [5*idx_counter+i] = Esys - (Es + Em - Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 

					// delete later
					Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
					if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
						std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
						std::cout << "Etest = " << Etest << ", energies [" << 5*idx_counter+i <<"] = " << energies[5*idx_counter+i] << "." << std::endl;
						std::cout << "contacts_store = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
						std::cout << "Shit's fucked." << std::endl;
						exit(EXIT_FAILURE);
					}
					// delete above later 

					// revert back to original structure 
					(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];

					// reset orientation 
					(*Polymers)[p_index].chain[m_index-1]->orientation       = original_ori; 
				}
			}
		}

		else {
			std::cout << "performing a solvent swap." << std::endl;
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
			
			for (int i{0}; i<5; ++i){
				// prep the swap 
				// std::cout << "solventswap: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index-1]->orientation     = test_ori[5*idx_counter+i];
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index-1]->coords  = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] =  (*Polymers)[p_index].chain[m_index-1];
				

				// get the energy 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
				Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 


				// std::cout << "Epair_n = " << Epair_n << std::endl; 
				// std::cout << "c_idx_n = " << c_idx_n << std::endl;

				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
				

				// delete later 
				Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
				if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << 5*idx_counter+i <<"] = " << energies[5*idx_counter+i] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				// delete above later 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index-1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index-1];

				// reset orientation 
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


	// do the swap again
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index-1]->coords           = ne_list[0];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ]  = (*Polymers)[p_index].chain[m_index-1];
	(*Polymers)[p_index].chain[m_index-1]->orientation = test_ori[0]; 

	// else { do nothing, MAINTAIN, because suggested index is in the maintain index vector }
	BackFlowFromTailRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, old_cut, old_ori, InteractionMap, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index-1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			Start of BackFlowFromTailRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			Start of HeadRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void HeadRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	std::array <double,8>				 current_contacts = *contacts; 
	std::array <double,25> 				 energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25> 				 boltzmann; 
	std::array <int,25> 				 orientations; 

	std::array <int,3> loc_m            = (*Polymers)[p_index].chain[m_index+1]->coords; 
	int                original_ori     = (*Polymers)[p_index].chain[m_index+1]->orientation; 

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
	double Epair      = 0;  
	double Epair_n    = 0; 
	double Esys       = *frontflow_energy; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int c_idx                     = 0 ; 
	int c_idx_n                   = 0 ; 
	int block_counter             = 0 ; 
	int idx_counter               = 0 ; 
	int self_swap_idx             = -1; 

	while ( idx_counter < 5 ){ // 5 locations to test 



		if ( ne_list[idx_counter] == loc_m ){ 
			// std::cout << "no location change, just orientation flip. " << std::endl;
			// this is just a solvent flip 
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z);
			for (int i{0}; i < 5; ++i ){

				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
				orientations[ 5*idx_counter + i ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation;

				Em_n = NeighborEnergetics ( LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				energies [5*idx_counter + i]  = Esys - Em + Em_n; 
				contacts_store [5*idx_counter + i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 

				//go back to the original state 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 

			}
			
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ){

			// check which region of the polymer the swap is taking place 

			for ( int u{0}; u<deg_poly; ++u ){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list[idx_counter] ){
					self_swap_idx = u; 
					break; 
				}
			}

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 

			for (int i{0}; i < 5; ++i ){

				if ( self_swap_idx < m_index ){
					// maintain current state, and sample another state 
					// maintain_idx.push_back(idx_counter); 
					// std::cout << "In the bad zone..." << std::endl;
					energies[5*idx_counter+i] = 1e+08; // very unfavorable state 
					orientations[5*idx_counter+i] = original_ori; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
					block_counter += 1; 
				}
				else {
					// std::cout << "monomer swap. " << std::endl;

					// perturb orientation 
					// (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
					(*Polymers) [p_index].chain[m_index+1]->orientation = rng_uniform(0,25); 
					orientations[5*idx_counter+i] = (*Polymers)[p_index].chain[m_index+1]->orientation; 

					// prep the swap 
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
					

					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

					// get the new energies
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
					Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 

					energies       [5*idx_counter+i] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
					

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
			// std::cout << "performing a solvent swap." << std::endl;
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
			
			

			for (int i{0}; i < 5; ++i ){ // 5 orientations to test 

				(*Polymers)[p_index].chain[m_index+1]->orientation = rng_uniform(0,25); 
				orientations[5*idx_counter+i] = (*Polymers)[p_index].chain[m_index+1]->orientation; 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the new energies
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
				Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 
				
				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair ; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];

				// reset orientation
				(*Polymers)[p_index].chain[m_index+1]->orientation = original_ori; 

			}
		}
		idx_counter += 1; 
	}

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
	
	// do the swap again
	
	(*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)]->coords = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords             = ne_list[e_idx/5];
	(*Polymers)[p_index].chain[m_index+1]->orientation        = orientations[e_idx];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx/5], y, z)];
	(*LATTICE)[ lattice_index (ne_list[e_idx/5], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];

	
	HeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void HeadRegrowthPlusOrientationFlip_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_o_to_n, \
	double* frontflow_energy, double temperature, int deg_poly, int p_index, int m_index, int x, int y, int z){


	if (m_index == deg_poly-1){
		*IMP_BOOL = true;
		return; 
	}

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 

	std::array <double,8>				 current_contacts = *contacts; 
	std::array <double,8>				 copy_contacts    = *contacts; 

	std::array <double,25> 				 energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25> 				 boltzmann; 
	std::array <int,25> 				 orientations; 

	std::array <int,3> loc_m            = (*Polymers)[p_index].chain[m_index+1]->coords; 
	int                original_ori     = (*Polymers)[p_index].chain[m_index+1]->orientation; 

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
	double Epair      = 0;  
	double Epair_n    = 0; 
	double Etest      = 0;
	double Esys       = *frontflow_energy; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int c_idx                     = 0 ; 
	int c_idx_n                   = 0 ; 
	int block_counter             = 0 ; 
	int idx_counter               = 0 ; 
	int self_swap_idx             = -1; 
	std::array <int,3> difference = {0,0,0};

	while ( idx_counter < 5 ){ // 5 locations to test 

		// std::cout << "ptype = " << (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype << std::endl;

		if ( ne_list[idx_counter] == loc_m ){ 
			// std::cout << "no location change, just orientation flip. " << std::endl;
			// this is just a solvent flip 
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z);
			for (int i{0}; i < 5; ++i ){

				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = rng_uniform(0,25); 
				orientations[ 5*idx_counter + i ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation;

				Em_n = NeighborEnergetics ( LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				energies [5*idx_counter + i]  = Esys - Em + Em_n; 
				contacts_store [5*idx_counter + i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 

				// delete later
				Etest = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
			    if ( Etest != energies[5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts  ){
			    	std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
					std::cout << "Etest = " << Etest << ", energies[j] = " << energies[5*idx_counter+i] << std::endl;
					std::cout << "contacts_store[j] = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "contacts_test = "; print (copy_contacts); 
					exit(EXIT_FAILURE); 
			    }
				// delete above later 
				//go back to the original state 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 

			}
			// std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
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

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			difference = subtract_arrays (&ne_list[idx_counter], &loc_m); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 

			// std::cout << "Neighbor bool is " << neighbor_bool << std::endl;

			for (int i{0}; i < 5; ++i ){

				if ( self_swap_idx < m_index ){
					energies[5*idx_counter+i] = 1e+08; // very unfavorable state 
					orientations[5*idx_counter+i] = original_ori; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
					block_counter += 1; 
				}
				else {
			
					(*Polymers) [p_index].chain[m_index+1]->orientation = rng_uniform(0,25); 
					orientations[5*idx_counter+i] = (*Polymers)[p_index].chain[m_index+1]->orientation; 

					// prep the swap 
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
					

					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

					// get the new energies
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 						
					Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 


					energies       [5*idx_counter+i] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
					

					// delete later 
					Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
					if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
						std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
						std::cout << "Etest = " << Etest << ", energies [" << 5*idx_counter+i <<"] = " << energies[5*idx_counter+i] << "." << std::endl;
						std::cout << "contacts_store = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
						std::cout << "Shit's fucked." << std::endl;
						exit(EXIT_FAILURE);
					}
					// delete above later 
					

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
			// std::cout << "performing a solvent swap." << std::endl;
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
			
			
			// std::cout << "Epair = " << Epair << std::endl; 
			// std::cout << "c_idx = " << c_idx << std::endl;
			// std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< std::endl;

			for (int i{0}; i < 5; ++i ){ // 5 orientations to test 

				(*Polymers)[p_index].chain[m_index+1]->orientation = rng_uniform(0,25); 
				orientations[5*idx_counter+i] = (*Polymers)[p_index].chain[m_index+1]->orientation; 

				// prep the swap 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

				// get the new energies
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
				Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 


				// std::cout << "Epair_n = " << Epair_n << std::endl; 
				// std::cout << "c_idx_n = " << c_idx_n << std::endl;
				
				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair ; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
				
				

				// get the energy
				// delete later 
				Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
				if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << 5*idx_counter+i <<"] = " << energies[5*idx_counter+i] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				// delete above later 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];

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
	HeadRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, prob_o_to_n, frontflow_energy, temperature, deg_poly, p_index, m_index+1, x, y, z);

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			end of HeadRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// 			Start of BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::vector <int>* old_ori, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == deg_poly-1){ 
		*IMP_BOOL = true; 
		return; 
	} 

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 

	std::array <double,8>                current_contacts = *contacts; 
	
	std::array <double,25>               energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25>               boltzmann; 

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

	// these are doubles for energies 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0;
	double Epair      = 0; 
	double Epair_n    = 0;  
	double Esys       = *backflow_energy; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int  c_idx                     = 0 ;
	int  c_idx_n			       = 0 ;
	int  idx_counter               = 0 ;
	int  self_swap_idx 		       = -1;

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);

		if ( ne_list[idx_counter] == loc_m ){
			// std::cout << "no location change, just orientation flip. " << std::endl;
			// this is just a solvent flip 
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z);
			for (int i{0}; i < 5; ++i ){
				// std::cout << "self-place: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index+1]->orientation = test_ori[5*idx_counter+i]; 

				Em_n = NeighborEnergetics ( LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				energies [5*idx_counter + i]  = Esys - Em + Em_n; 
				contacts_store [5*idx_counter + i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 

				//go back to the original state 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 	
			
			}

			// std::cout << "current_contacts = "; print(current_contacts); 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ) {

			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list [idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}
			}

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 

	
			for (int i{0}; i<5; ++i){

				if ( self_swap_idx < m_index ){

					energies [5*idx_counter+i]      = 1e+08; 
					contacts_store[5*idx_counter+i] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				
				}
				
				else {
					
					// prep the swap 
					(*Polymers)[p_index].chain[m_index+1]->orientation     = test_ori[5*idx_counter+i];
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
					
					// perform the swap 
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
					

					// get the new energies
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z);
					Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 
					

					energies       [5*idx_counter+i] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
					

					// revert back to original structure 
					(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];

					// reset orientation 
					(*Polymers)[p_index].chain[m_index+1]->orientation       = original_ori; 
				}
			}
		}

		else {

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
			

			for (int i{0}; i<5; ++i){
				// prep the swap 
				// std::cout << "solventswap: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index+1]->orientation            = test_ori[5*idx_counter+i]; 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ]                = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]; 
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1]; 
				

				// get the energy 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 
				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
				

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];
				
				// reset orientation 
				(*Polymers)[p_index].chain[m_index+1]->orientation       = original_ori; 
			}

		}

		idx_counter += 1;

	}

	double Emin = *std::min_element ( energies.begin(), energies.end() ); 


	for (int i{0}; i<25; ++i){
		boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
		rboltzmann  += boltzmann[i]; 
	}


	*prob_n_to_o      = (*prob_n_to_o)*boltzmann[0]/rboltzmann;
	*backflow_energy  = energies[0]; 
	*contacts         = contacts_store [0];



	// do the swap again
	(*LATTICE)[lattice_index(ne_list[0], y, z)]->coords     = loc_m;
	(*Polymers)[p_index].chain[m_index+1]->coords           = ne_list[0];

	// perform the swap (since coords were changed, this swap works)
	(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[0], y, z)];
	(*LATTICE)[ lattice_index (ne_list[0], y, z) ]  = (*Polymers)[p_index].chain[m_index+1];
	(*Polymers)[p_index].chain[m_index+1]->orientation = test_ori[0]; 

	// else { do nothing, MAINTAIN, because suggested index is in the maintain index vector }
	BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, old_cut, old_ori, InteractionMap, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::vector<std::array<int,3>>* old_cut, std::vector <int>* old_ori, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* prob_n_to_o, \
	double* backflow_energy, double temperature, int deg_poly, int p_index, int m_index, int recursion_depth, int x, int y, int z){


	if (m_index == deg_poly-1){ 
		*IMP_BOOL = true; 
		return; 
	} 

	// std::cout << "m_index = " << m_index << std::endl;
	// generate an array for energies 

	std::array <double,8>                copy_contacts    = *contacts; 
	std::array <double,8>                current_contacts = *contacts; 
	
	std::array <double,25>               energies; 
	std::array <std::array<double,8>,25> contacts_store; 
	std::array <double,25>               boltzmann; 

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

	// these are doubles for energies 
	double rboltzmann = 0; // running sum for boltzmann weights 
	double Em         = 0;
	double Es         = 0; 
	double Em_n       = 0; 
	double Es_n       = 0;
	double Epair      = 0; 
	double Epair_n    = 0;  
	double Etest      = 0;
	double Esys       = *backflow_energy; 

	std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};

	// start attempting jumps 
	int  c_idx                     = 0 ;
	int  c_idx_n			       = 0 ;
	int  idx_counter               = 0 ;
	int  self_swap_idx 		       = -1;

	while ( idx_counter < 5 ){
		// std::cout << "Suggested location = "; print(ne_list[idx_counter]);

		if ( ne_list[idx_counter] == loc_m ){
			std::cout << "no location change, just orientation flip. " << std::endl;
			// this is just a solvent flip 
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z);
			for (int i{0}; i < 5; ++i ){
				// std::cout << "self-place: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index+1]->orientation = test_ori[5*idx_counter+i]; 

				Em_n = NeighborEnergetics ( LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 
				energies [5*idx_counter + i]  = Esys - Em + Em_n; 
				contacts_store [5*idx_counter + i] = add_arrays ( subtract_arrays (current_contacts, cm), cm_n); 

				// delete later
				Etest = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
			    if ( Etest != energies[5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
			    	std::cout << "Something is fucked. Either energies do not match or..." << std::endl;
					std::cout << "Etest = " << Etest << ", energies[j] = " << energies[5*idx_counter+i] << std::endl;
					std::cout << "contacts_store[j] = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "contacts_test = "; print (copy_contacts); 
					exit(EXIT_FAILURE); 
			    }
				// delete above later 
				//go back to the original state 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->orientation = original_ori; 	
			
			}

			// std::cout << "current_contacts = "; print(current_contacts); 
		}

		else if ( (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ]->ptype[0] == 'm' ) {

			for ( int u{0}; u<deg_poly; ++u){
				if ( (*Polymers)[p_index].chain[u]->coords == ne_list [idx_counter] ){
					self_swap_idx = u; 
					// std::cout << "selfswap = " << self_swap_idx << ", u = " << u << std::endl;
					break; 
				}
			}

			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
			

			for (int i{0}; i<5; ++i){

				if ( self_swap_idx < m_index ){
					
					(*Polymers)[p_index].chain[m_index+1]->orientation = original_ori;
					energies [5*idx_counter+i]                         = 1e+08; 
					contacts_store[5*idx_counter+i]                    = {-1,-1,-1,-1,-1,-1,-1,-1}; 
				
				}
				
				else {

					std::cout << "monomer swap: idx_counter = " << idx_counter << std::endl;
					
					// prep the swap 
					(*Polymers)[p_index].chain[m_index+1]->orientation     = test_ori[5*idx_counter+i];
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
					(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
					
					// perform the swap 
					(*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
					(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];
					

					// get the new energies
					Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
					Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
					Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 
					

					energies       [5*idx_counter+i] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
					contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
					contacts_store [5*idx_counter+i][c_idx]   += 1; 
					contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
					

					// delete later
					Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
					if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
						std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
						std::cout << "Etest = " << Etest << ", energies [" << 5*idx_counter+i <<"] = " << energies[5*idx_counter+i] << "." << std::endl;
						std::cout << "contacts_store = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
						std::cout << "Shit's fucked." << std::endl;
						exit(EXIT_FAILURE);
					}
					
					// delete above later 
					// std::cout << "current_contacts = "; print(current_contacts); 

					// revert back to original structure 
					(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
					(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
					
					// perform the swap (since coords were changed, this swap works)
					(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
					(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];

					// reset orientation 
					(*Polymers)[p_index].chain[m_index+1]->orientation       = original_ori; 
				}
			}
		}

		else {

			std::cout << "performing a solvent swap." << std::endl;
			Es = NeighborEnergetics (LATTICE, InteractionMap, &cs, lattice_index(ne_list[idx_counter], y, z), x, y, z);
			Em = NeighborEnergetics (LATTICE, InteractionMap, &cm, lattice_index(loc_m, y, z), x, y, z); 
			Epair = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx, x, y, z); 
			

			std::cout << "Epair = " << Epair << std::endl; 
			std::cout << "c_idx = " << c_idx << std::endl;
			std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< std::endl;			

			for (int i{0}; i<5; ++i){
				// prep the swap 
				// std::cout << "solventswap: idx_counter = " << idx_counter << std::endl;
				(*Polymers)[p_index].chain[m_index+1]->orientation       = test_ori[5*idx_counter+i]; 
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
				(*Polymers)[p_index].chain[m_index+1]->coords  = ne_list[idx_counter];
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[ lattice_index (loc_m, y, z) ]                = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]; 
				(*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1]; 	

				// get the energy 
				Es_n = NeighborEnergetics (LATTICE, InteractionMap, &cs_n, lattice_index(ne_list[idx_counter], y, z), x, y, z);
				Em_n = NeighborEnergetics (LATTICE, InteractionMap, &cm_n, lattice_index(loc_m, y, z), x, y, z); 	
				Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[lattice_index(ne_list[idx_counter], y, z)], (*LATTICE)[lattice_index(loc_m, y, z)], InteractionMap, &c_idx_n, x, y, z); 
				

				std::cout << "Epair_n = " << Epair_n << std::endl; 
				std::cout << "c_idx_n = " << c_idx_n << std::endl;

				energies       [5*idx_counter+i] = Esys - (Es+Em) + (Es_n+Em_n) - Epair_n + Epair; 
				contacts_store [5*idx_counter+i] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
				contacts_store [5*idx_counter+i][c_idx]   += 1; 
				contacts_store [5*idx_counter+i][c_idx_n] -= 1; 
				

				// delete later 
				Etest  = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
				if ( Etest != energies [5*idx_counter+i] || contacts_store[5*idx_counter+i] != copy_contacts ){
					std::cout << "Energies are bad, or contacts are not right, or solvation shells are messed up." << std::endl;
					std::cout << "Etest = " << Etest << ", energies [" << 5*idx_counter+i <<"] = " << energies[5*idx_counter+i] << "." << std::endl;
					std::cout << "contacts_store = "; print (contacts_store[5*idx_counter+i], ", "); std::cout << "copy_contacts = "; print(copy_contacts);
					std::cout << "Shit's fucked." << std::endl;
					exit(EXIT_FAILURE);
				}
				// delete above later 

				// revert back to original structure 
				(*LATTICE)[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
				(*Polymers)[p_index].chain[m_index+1]->coords  = loc_m;
				
				// perform the swap (since coords were changed, this swap works)
				(*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]    = (*LATTICE)[lattice_index (loc_m, y, z)];
				(*LATTICE)[lattice_index(loc_m, y, z)]                   = (*Polymers)[p_index].chain[m_index+1];
				
				// reset orientation 
				(*Polymers)[p_index].chain[m_index+1]->orientation       = original_ori; 
			}
			std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< std::endl;

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
	BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, old_cut, old_ori, InteractionMap, E, contacts, IMP_BOOL, prob_n_to_o, backflow_energy, temperature, deg_poly, p_index, m_index+1, recursion_depth+1, x, y, z); 

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

void ChainRegrowthPlusOrientationFlip_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, \
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
		

		for (int i {m_index+1}; i<deg_poly; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
			old_ori.push_back ((*Polymers)[p_index].chain[i]->orientation) ; 
		}


		// regrow the polymer frontwards
		HeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 

		for (int i {m_index+1}; i<deg_poly; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
			new_ori.push_back ((*Polymers)[p_index].chain[i]->orientation) ;
		}

		// std::cout << "new orientation = "; print (new_ori);

		if ( old_cut == new_cut && old_ori == new_ori){
			
			return; 
		}

		if ( !(*IMP_BOOL) ){

			acceptance_after_head_regrowth ( LATTICE, &new_cut, &old_cut, y, z );

			// change orientations of polymer bead to old 
			for (int i{m_index+1}; i<deg_poly; ++i){
				(*Polymers)[p_index].chain[i]->orientation = old_ori[i-m_index-1]; 
			}

			return; 
		}

		c2_contacts     = c1_contacts;
		backflow_energy = frontflow_energy; 
		
		BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, &old_ori, InteractionMap, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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
		
		for ( int i{0}; i<m_index; ++i){
			old_cut.push_back ((*Polymers)[p_index].chain[i]->coords);
			old_ori.push_back ((*Polymers)[p_index].chain[i]->orientation);
		}


		TailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 

		for (int i {0}; i<m_index; ++i){
			new_cut.push_back ((*Polymers)[p_index].chain[i]->coords) ;
			new_ori.push_back ((*Polymers)[p_index].chain[i]->orientation) ;
		}

		// std::cout << "new orientation = "; print (new_ori);

		// std::cout << "created new_cut and new_ori!" << std::endl; 

		if ( old_cut == new_cut && old_ori == new_ori ){
			return; 
		}

		
		if ( !(*IMP_BOOL) ){

			acceptance_after_tail_regrowth ( LATTICE, &new_cut, &old_cut, y, z); 
			
			for (int i{0}; i<m_index; ++i){
				(*Polymers)[p_index].chain[i]->orientation = old_ori[i]; 
			}
			
			return; 
		}

		c2_contacts     = c1_contacts; 
		backflow_energy = frontflow_energy; 

		BackFlowFromTailRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, &old_cut, &old_ori, InteractionMap, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

		// check acceptance criterion 
		double rng_acc = rng_uniform (0.0, 1.0); 
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n ){

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
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void ChainRegrowthPlusOrientationFlip_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, \
	double* sysEnergy, double temperature, int p_index, int x, int y, int z) {


	int deg_poly = (*Polymers)[p_index].deg_poly; 
	int m_index  = rng_uniform (1, deg_poly-1); 

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
		
		std::cout << "Head regrowth... " << std::endl;
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
		HeadRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		std::cout << "Completed regrowth!" << std::endl;
		// exit (EXIT_SUCCESS);
		// CheckStructures (Polymers, Cosolvent, LATTICE, x, y, z); 

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
		c2_contacts     = c1_contacts;
		backflow_energy = frontflow_energy; 
		std::cout << "Begin backflow... " << std::endl;
		
		BackFlowFromHeadRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, &old_cut, &old_ori, InteractionMap, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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
		

		std::cout << "Initiate tail regrowth..." << std::endl;
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

		TailRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, &c1_contacts, IMP_BOOL, &prob_o_to_n, &frontflow_energy, temperature, deg_poly, p_index, m_index, x, y, z); 
		std::cout << "Completed tail regrowth!" << std::endl;
		// exit(EXIT_SUCCESS);
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

		c2_contacts     = c1_contacts; 
		backflow_energy = frontflow_energy; 

		// std::cout << "BEGIN BACK FLOW for Tail! " << std::endl;

		BackFlowFromTailRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, &old_cut, &old_ori, InteractionMap, E, &c2_contacts, IMP_BOOL, &prob_n_to_o, &backflow_energy, temperature, deg_poly, p_index, m_index, recursion_depth, x, y, z); 

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
	
		if ( rng_acc < std::exp(-1/temperature * (frontflow_energy - *sysEnergy)) * prob_n_to_o/prob_o_to_n ){
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

////////////////////////////////////////////////////////////////////////////////////////////
//

void SolventExchange_UNBIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int x, int y, int z) {

	std::set   <int>                   solvation_shell_set; 
	std::array <std::array<int,3>, 26> ne_list, ne_list_; 

	std::cout << "*E[0] = " << (*E)[0] << std::endl;

	// get the first solvation shell 
	for ( Polymer& pmer: *Polymers ){
		for ( Particle*& p: pmer.chain ) {
			ne_list = obtain_ne_list ( p->coords, x, y, z );
			for ( std::array <int,3>& loc: ne_list ){
				if ( (*LATTICE) [ lattice_index (loc, y, z) ]->ptype[0] == 's' ) {
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

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));

	int exc_idx   = -1; 

	double                frontflow_energy   = 0;
	std::array <double,8> frontflow_contacts = *contacts;
	std::array <double,8> copy_contacts      = *contacts;
	std::array <double,8> cs1                = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2                = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs1_n              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_n              = {0,0,0,0,0,0,0,0};

	int c_idx       = 0;
	int c_idx_n     = 0;

	double Esys          = *sysEnergy;
	double Es1           = 0;
	double Es1_n         = 0;
	double Es2           = 0;
	double Es2_n         = 0;
	double Epair         = 0;
	double Epair_n       = 0;

	Particle* tmp_par_ptr {nullptr};
	exc_idx = rng_uniform (0, x*y*z-1);
	int my_idx = 0; 

	if ( exc_idx == solvation_shell_indices[my_idx] ){
		// 
		return;
	}

	if ( (*LATTICE) [exc_idx]->ptype[0] == 's' ) {

		Es1 = NeighborEnergetics (LATTICE, InteractionMap, &cs1, solvation_shell_indices [my_idx], x, y, z);
		Es2 = NeighborEnergetics (LATTICE, InteractionMap, &cs2, exc_idx, x, y, z);
		Epair = IsolatedPairParticleInteraction  ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx, x, y, z); 

		tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[ my_idx ] ];
		
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ] = (*LATTICE)[ exc_idx ];
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ]->coords = location ( (solvation_shell_indices)[ my_idx ], x, y, z);

		(*LATTICE)[ exc_idx ] = tmp_par_ptr; 
		(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); 

		// flip the particle newly added to the solvation shell 
		Es2_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
		Es1_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
		Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z); 
		

		frontflow_energy   = Esys - (Es1 + Es2 - Epair) + (Es1_n + Es2_n - Epair_n);  
		frontflow_contacts = add_arrays ( subtract_arrays (*contacts, add_arrays (cs1, cs2) ), add_arrays (cs1_n, cs2_n) );
		frontflow_contacts [c_idx]   += 1;
		frontflow_contacts [c_idx_n] -= 1;
		

		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		double energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &copy_contacts, x, y, z); 
		if ( energy_n != frontflow_energy || copy_contacts != frontflow_contacts ){
			std::cout << "Either energy of contacts is messed up in id swap..." << std::endl;
			std::cout << "energy_n = " << energy_n << ", frontflow_energy = " << frontflow_energy << "." << std::endl;
			std::cout << "copy_contacts = "; print (copy_contacts, ", "); std::cout << "frontflow_contacts = "; print (frontflow_contacts);
			exit (EXIT_FAILURE);
		}
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}

	else if ( (*LATTICE) [exc_idx]->ptype == "m1" ) {
		(*IMP_BOOL) = false;
		return;
	}
	
	double rng_acc = rng_uniform (0.0, 1.0); 
	if ( rng_acc < std::exp (-1 / temperature * (frontflow_energy - *sysEnergy) ) ) {
		*sysEnergy = frontflow_energy;
		*contacts  = frontflow_contacts;
	}
	else {
		tmp_par_ptr = (*LATTICE) [ solvation_shell_indices [ my_idx ] ];
		(*LATTICE)[ solvation_shell_indices [ my_idx ] ]         = (*LATTICE)[ exc_idx ]; 
		(*LATTICE)[ solvation_shell_indices [ my_idx ] ]->coords = location( solvation_shell_indices[ my_idx ], x, y, z); 

		(*LATTICE)[ exc_idx ]         = tmp_par_ptr;
		(*LATTICE)[ exc_idx ]->coords = location (exc_idx, x, y, z);
		
	}
	CheckStructures (Polymers, Cosolvent, LATTICE, x, y ,z);

	return;

}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


void SolventExchange_UNBIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int x, int y, int z) {

	std::set   <int>                   solvation_shell_set;
	std::array <std::array<int,3>, 26> ne_list, ne_list_  ;

	// get the first solvation shell
	for ( Polymer& pmer: *Polymers ){
		for ( Particle*& p: pmer.chain ) {
			ne_list = obtain_ne_list ( p->coords, x, y, z );
			for ( std::array <int,3>& loc: ne_list ){
				if ( (*LATTICE) [ lattice_index (loc, y, z) ]->ptype[0] == 's' ) {
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

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle ( solvation_shell_indices.begin(), solvation_shell_indices.end(), std::default_random_engine(seed));

	int exc_idx   = -1;

	double                frontflow_energy   = 0;
	std::array <double,8> frontflow_contacts = *contacts;
	std::array <double,8> cs1                = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2                = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs1_n              = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_n              = {0,0,0,0,0,0,0,0};

	int c_idx       = 0;
	int c_idx_n     = 0;

	double Esys           = *sysEnergy;
	double Es1            = 0;
	double Es1_n          = 0;
	double Es2            = 0;
	double Es2_n          = 0;
	double Epair          = 0;
	double Epair_n        = 0;

	Particle* tmp_par_ptr {nullptr};
	exc_idx = rng_uniform (0, x*y*z-1);
	int my_idx = 0; 

	if ( exc_idx == solvation_shell_indices[my_idx] ){
		return;
	}

	if ( (*LATTICE) [exc_idx]->ptype[0] == 's' ) {

		Es1 = NeighborEnergetics (LATTICE, InteractionMap, &cs1, solvation_shell_indices [my_idx], x, y, z);
		Es2 = NeighborEnergetics (LATTICE, InteractionMap, &cs2, exc_idx, x, y, z);
		Epair = IsolatedPairParticleInteraction  ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx, x, y, z); 


		tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[ my_idx ] ];
		
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ] = (*LATTICE)[ exc_idx ];
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ]->coords = location ( (solvation_shell_indices)[ my_idx ], x, y, z);

		(*LATTICE)[ exc_idx ] = tmp_par_ptr;
		(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z );

		// flip the particle newly added to the solvation shell
		Es2_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
		Es1_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
		Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z);


		frontflow_energy   = Esys - (Es1 + Es2-Epair) + (Es1_n + Es2_n-Epair_n);  
		frontflow_contacts = add_arrays ( subtract_arrays (*contacts, add_arrays (cs1, cs2) ), add_arrays (cs1_n, cs2_n) );
		frontflow_contacts [c_idx]   += 1;
		frontflow_contacts [c_idx_n] -= 1;
		
	}

	else if ( (*LATTICE) [exc_idx]->ptype == "m1" ) {
		(*IMP_BOOL) = false;
		return;
	}
	
	double rng_acc = rng_uniform (0.0, 1.0); 
	if ( rng_acc < std::exp (-1 / temperature * (frontflow_energy - *sysEnergy) ) ) {
		*sysEnergy = frontflow_energy;
		*contacts  = frontflow_contacts;
	}
	else {
		tmp_par_ptr = (*LATTICE) [ solvation_shell_indices [ my_idx ] ];
		(*LATTICE)[ solvation_shell_indices [ my_idx ] ]         = (*LATTICE)[ exc_idx ]; 
		(*LATTICE)[ solvation_shell_indices [ my_idx ] ]->coords = location ( solvation_shell_indices[ my_idx ], x, y, z); 

		(*LATTICE)[ exc_idx ]         = tmp_par_ptr;
		(*LATTICE)[ exc_idx ]->coords = location (exc_idx, x, y, z);
	}
	
	return;

}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//             End of acceptance_after_tail_regrowth
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void SolventExchange_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
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

	int                                 ntest              = 5; 
	std::array <double,5>               energies           = {0,0,0,0,0}; 
	std::array <double,5>               boltzmann          = {0,0,0,0,0};
	std::array <int,5>                  orientations       = {0,0,0,0,0}; 
	double                              prob_o_to_n        = 1; 
	double                              prob_n_to_o        = 1; 
	double                              frontflow_energy   = 0; 
	std::array <double,8>               frontflow_contacts = *contacts; 
	std::array<std::array<double,8>,5>  contacts_store     = {*contacts, *contacts, *contacts, *contacts, *contacts}; 
	double                              Emin               = 0; 
	double                              rng_acc            = 0;
	double                              rng                = 0; 
	double                              rsum               = 0; 
	double                              rboltzmann         = 0; 
	int                                 e_idx              = 0; 
    
	std::array <double,8> cs1     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs1_n   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_n   = {0,0,0,0,0,0,0,0};

	int c_idx       = 0;
	int c_idx_n     = 0; 

	double Esys          = *sysEnergy; 
	double Es1           = 0; 
	double Es1_n         = 0; 
	double Es2           = 0;
	double Es2_n         = 0;
	double Epair         = 0; 
	double Epair_n       = 0; 

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
    
    if ( exc_idx == solvation_shell_indices[my_idx] ){
    	// std::cout << "same particle selected." << std::endl;
    	return;
    }

    if ( (*LATTICE)[exc_idx]->ptype[0] == 's' ){
        // std::cout << "Inside loop..." << std::endl;
        // made a copy of the particle in the solvation shell
        // 
		// swap particles
		// initial neighbor energies and contacts
		Es1 = NeighborEnergetics (LATTICE, InteractionMap, &cs1, solvation_shell_indices[my_idx], x, y, z);
		Es2 = NeighborEnergetics (LATTICE, InteractionMap, &cs2, exc_idx, x, y, z);
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx, x, y, z); 

		tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[ my_idx ] ];
		old_ori.push_back ((*LATTICE)[ exc_idx ]->orientation);
		
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ] = (*LATTICE)[ exc_idx ];
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ]->coords = location ( (solvation_shell_indices)[ my_idx ], x, y, z);

		(*LATTICE)[ exc_idx ] = tmp_par_ptr; 
		(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); 

		// flip the particle newly added to the solvation shell 
		for (int j{0}; j < ntest; ++j) { 
		// std::cout << " j = " << j << std::endl;
			(*LATTICE)[ solvation_shell_indices[ my_idx ] ]->orientation = rng_uniform (0, 25); 
			orientations   [ j ] = (*LATTICE)[ solvation_shell_indices[ my_idx ] ]->orientation; 
			Es2_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
			Es1_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z); 

			energies       [ j ] = Esys - (Es1 + Es2-Epair) + (Es1_n + Es2_n-Epair_n);  
			contacts_store [ j ] = add_arrays ( subtract_arrays (*contacts, add_arrays (cs1, cs2) ), add_arrays (cs1_n, cs2_n) );
			contacts_store [j][c_idx]   += 1; 
			contacts_store [j][c_idx_n] -= 1; 
			
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
	frontflow_contacts   = contacts_store [e_idx]; // the chosen contacts

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// reversing the flow 

	// made a copy of the particle in the solvation shell 
	// swap particles 

	tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[my_idx] ];
	
	Es1 = NeighborEnergetics (LATTICE, InteractionMap, &cs1, solvation_shell_indices[my_idx] , x, y, z); 
	Es2 = NeighborEnergetics (LATTICE, InteractionMap, &cs2, exc_idx, x, y, z); 
	Epair = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx, x, y, z); 


	(*LATTICE)[ solvation_shell_indices [my_idx] ]         = (*LATTICE)[ exc_idx ]; 
	(*LATTICE)[ solvation_shell_indices [my_idx] ]->coords = location ( solvation_shell_indices[my_idx], x, y, z); 

	(*LATTICE)[ exc_idx ]         = tmp_par_ptr; 
	(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); // this particle is the particle that was in the solvation shell 

	// flip the particle sent away 
	(*LATTICE) [ exc_idx ]->orientation = old_ori[0]; 

	Es2_n   = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
	Es1_n   = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
	Epair_n = IsolatedPairParticleInteraction  ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z); 

	energies [ 0 ]       = frontflow_energy - (Es1 + Es2 - Epair) + (Es1_n + Es2_n - Epair_n); 
	contacts_store [ 0 ] = add_arrays( subtract_arrays(frontflow_contacts, add_arrays (cs1, cs2)), add_arrays(cs1_n, cs2_n) );
	contacts_store [0][c_idx]   += 1; 
	contacts_store [0][c_idx_n] -= 1; 
	

	for (int j{1}; j < ntest; ++j){

		(*LATTICE)[ exc_idx ]->orientation = rng_uniform (0, 25); 
		Es2_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
		Es1_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z); 

		energies       [ j ] = frontflow_energy - (Es1+Es2-Epair) + (Es1_n+Es2_n-Epair_n); 
		contacts_store [ j ] = add_arrays( subtract_arrays(frontflow_contacts, add_arrays (cs1, cs2)), add_arrays(cs1_n, cs2_n) );
		contacts_store [j][c_idx]   += 1; 
		contacts_store [j][c_idx_n] -= 1; 
		
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
		*contacts  = frontflow_contacts;
	}
	else {
		*IMP_BOOL = false;
		// get the particle back to its original orientation 
		(*LATTICE) [ exc_idx ]->orientation = old_ori[0]; 
	}	


	return; 

}

// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

void SolventExchange_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, 
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, bool* IMP_BOOL, double* sysEnergy, \
	double temperature, int x, int y, int z) {


	std::cout << "*E[0] = " << (*E)[0] << std::endl;
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

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	
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
    std::array <double,8>               frontflow_contacts = *contacts; 
	std::array<std::array<double,8>,5>  contacts_store   = {*contacts, *contacts, *contacts, *contacts, *contacts}; 
    std::array <double,8>               c_contacts1      = *contacts; 
    double                              Emin             = 0; 
    double                              rng_acc          = 0;
    double                              rng              = 0; 
    double                              rsum             = 0; 
    double                              rboltzmann       = 0; 
    int                                 e_idx            = 0; 
    
	std::array <double,8> cs1     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2     = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs1_n   = {0,0,0,0,0,0,0,0};
	std::array <double,8> cs2_n   = {0,0,0,0,0,0,0,0};

	int c_idx       = 0;
	int c_idx_n     = 0; 

    double Esys          = *sysEnergy; 
    double Es1           = 0; 
    double Es1_n         = 0; 
    double Es2           = 0;
    double Es2_n         = 0;
    double Epair         = 0; 
    double Epair_n       = 0; 
    double energy_n      = 0; 

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
    std::cout << "solvation_shell_indices[0] = " << solvation_shell_indices[0] << ", coords = "; 
    print ( (*LATTICE)[solvation_shell_indices[0]]->coords, ", " ); std::cout << "ptype = " <<  (*LATTICE)[solvation_shell_indices[0]]->ptype << ", or = " << (*LATTICE)[solvation_shell_indices[0]]->orientation << std::endl; 
    std::cout << "exc_idx = " << exc_idx << ", coords = "; print ( (*LATTICE)[exc_idx]->coords, ", " ); std::cout << "ptype = " << (*LATTICE)[exc_idx]->ptype << ", or = " << (*LATTICE)[exc_idx]->orientation << std::endl;

    if ( exc_idx == solvation_shell_indices[my_idx] ){
    	std::cout << "same particle selected." << std::endl;
    	return;
    }

    if ( (*LATTICE)[exc_idx]->ptype[0] == 's' ) {

		Es1 = NeighborEnergetics (LATTICE, InteractionMap, &cs1, solvation_shell_indices[my_idx], x, y, z);
		Es2 = NeighborEnergetics (LATTICE, InteractionMap, &cs2, exc_idx, x, y, z);
		Epair = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx, x, y, z); 

		tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[ my_idx ] ];
		old_ori.push_back ((*LATTICE)[ exc_idx ]->orientation);
		
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ] = (*LATTICE)[ exc_idx ];
		(*LATTICE)[ (solvation_shell_indices)[ my_idx ] ]->coords = location ( (solvation_shell_indices)[ my_idx ], x, y, z);

		(*LATTICE)[ exc_idx ] = tmp_par_ptr; 
		(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); 

		// flip the particle newly added to the solvation shell 
		for (int j{0}; j < ntest; ++j) { 
            // std::cout << " j = " << j << std::endl;
			(*LATTICE)[ solvation_shell_indices[ my_idx ] ]->orientation = rng_uniform (0, 25); 
			orientations   [ j ] = (*LATTICE)[ solvation_shell_indices[ my_idx ] ]->orientation; 
			Es2_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
			Es1_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
			Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z); 
			
			energies       [ j ] = Esys - (Es1 + Es2-Epair) + (Es1_n + Es2_n-Epair_n);  
			contacts_store [ j ] = add_arrays ( subtract_arrays (*contacts, add_arrays (cs1, cs2) ), add_arrays (cs1_n, cs2_n) );		
			contacts_store [j][c_idx]   += 1; 
			contacts_store [j][c_idx_n] -= 1; 
			

		
			// DELETE THIS LATER 
            
			energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &c_contacts1, x, y, z);
			if (energies[j] != energy_n || contacts_store[j] != c_contacts1) {
				std::cout << "Either energy or contacts is messed up in id swap in frontflow... " << std::endl;
				std::cout << "energies["<<j<<"] = " << energies[j] << ", energy_n = " << energy_n << ". " << std::endl; 
				std::cout << "contacts_store[" <<j<< "] = "; print (contacts_store[j], ", "); std::cout << "c_contacts = "; print(c_contacts1);
				exit(EXIT_FAILURE);
			}
             
			// DELETE ABOVE LATER 
		
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
	frontflow_contacts   = contacts_store [e_idx]; // the chosen contacts

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// reversing the flow 
	// made a copy of the particle in the solvation shell 
	// swap particles 

	tmp_par_ptr = (*LATTICE)[ solvation_shell_indices[my_idx] ];
	
	Es1 = NeighborEnergetics (LATTICE, InteractionMap, &cs1, solvation_shell_indices[my_idx] , x, y, z);  
	Es2 = NeighborEnergetics (LATTICE, InteractionMap, &cs2, exc_idx, x, y, z); 
	Epair = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx, x, y, z); 
	

	(*LATTICE)[ solvation_shell_indices [my_idx] ]         = (*LATTICE)[ exc_idx ]; 
	(*LATTICE)[ solvation_shell_indices [my_idx] ]->coords = location ( solvation_shell_indices[my_idx], x, y, z); 

	(*LATTICE)[ exc_idx ]         = tmp_par_ptr; 
	(*LATTICE)[ exc_idx ]->coords = location ( exc_idx, x, y, z ); // this particle is the particle that was in the solvation shell 

	// flip the particle sent away 
	(*LATTICE) [ exc_idx ]->orientation = old_ori[0]; 

	Es2_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
	Es1_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
	Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z); 
	

	energies [ 0 ]       = frontflow_energy - (Es1 + Es2 - Epair) + (Es1_n + Es2_n - Epair_n); // CalculateEnergy ( Polymers, Cosolvent, LATTICE, E, &c_contacts2, x, y, z ); 
	contacts_store [ 0 ] = add_arrays( subtract_arrays(frontflow_contacts, add_arrays (cs1, cs2)), add_arrays(cs1_n, cs2_n) );
	contacts_store [0][c_idx]   += 1; 
	contacts_store [0][c_idx_n] -= 1; 
	
	// DELETE THIS LATER
        
	energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &c_contacts1, x, y, z);
	if (energies[0] != energy_n || contacts_store[0] != c_contacts1) {
		std::cout << "Either energy or contacts is messed up in id swap... " << std::endl;
		std::cout << "energies[0] = " << energies[0] << ", energy_n = " << energy_n << ". " << std::endl; 
		std::cout << "contacts_store = "; print (contacts_store[0], ", "); std::cout << "c_contacts1 = "; print(c_contacts1);
		exit(EXIT_FAILURE);
	}
         
	// DELETE ABOVE LATER 

	for (int j{1}; j < ntest; ++j){

		(*LATTICE)[ exc_idx ]->orientation = rng_uniform (0, 25); 
		Es2_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs2_n, solvation_shell_indices[my_idx], x, y, z);
		Es1_n                = NeighborEnergetics  (LATTICE, InteractionMap, &cs1_n, exc_idx, x, y, z);
		Epair_n = IsolatedPairParticleInteraction ((*LATTICE)[solvation_shell_indices[my_idx] ], (*LATTICE)[exc_idx], InteractionMap, &c_idx_n, x, y, z); 
		

		energies       [ j ] = frontflow_energy - (Es1+Es2-Epair) + (Es1_n+Es2_n-Epair_n); 
		contacts_store [ j ] = add_arrays( subtract_arrays(frontflow_contacts, add_arrays (cs1, cs2)), add_arrays(cs1_n, cs2_n) );
		contacts_store [j][c_idx]   += 1; 
		contacts_store [j][c_idx_n] -= 1; 
				

		// DELETE THIS LATER 
        
		energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &c_contacts1, x, y, z);
		if (energies[j] != energy_n || contacts_store[j] != c_contacts1) {
			std::cout << "Either energy or contacts is messed up in id swap in backflow... " << std::endl;
			std::cout << "energies[j] = " << energies[j] << ", energy_n = " << energy_n << ". " << std::endl; 
			std::cout << "contacts_store[j] = "; print (contacts_store[j], ", "); std::cout << "c_contacts = "; print(c_contacts1);
			exit(EXIT_FAILURE);
		}
        
		// DELETE ABOVE LATER 
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

		// DELETE THIS LATER 

		energy_n = CalculateEnergyRevamped (Polymers, Cosolvent, LATTICE, InteractionMap, &c_contacts1, x, y, z);
		if (frontflow_energy != energy_n || frontflow_contacts != c_contacts1) {
			std::cout << "Either energy or contacts is messed up in id swap post backflow... " << std::endl;
			std::cout << "frontflow_energy = " << frontflow_energy << ", energy_n = " << energy_n << ". " << std::endl; 
			std::cout << "frontflow_contacts = "; print (frontflow_contacts, ", "); std::cout << "c_contacts1 = "; print(c_contacts1);
			exit(EXIT_FAILURE);
		}
        
		// DELETE ABOVE LATER 

		*sysEnergy = frontflow_energy;
		*contacts  = frontflow_contacts;
	}
	else {
		*IMP_BOOL = false;
		(*LATTICE) [ exc_idx ]->orientation = old_ori[0]; 
	}	


	return; 

}

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


void PerturbSystem_BIASED (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE,  \
	std::map <std::pair<std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* contacts, std::array <int,9>* attempts, \
	bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, \
	int* move_number, int x, int y, int z) {

	int index = 0; // rng_uniform (0, static_cast<int>((*Polymers).size()-1) ); 
	int r     = rng_uniform (0, 8); 
	// std::array <double,4> c_contacts = *contacts; 

	switch (r) {

		case(0):
			if (v) {
				std::cout << "Performing end rotations..." << std::endl; 
			}
			EndRotation_UNBIASED  (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break; 

		case (1):
			if (v) {
				std::cout << "Performing reptation..." << std::endl; 
			}
			Reptation_UNBIASED    (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (2):
			if (v) {
				std::cout << "Performing biased chain regrowth..." << std::endl;
			}
			ChainRegrowth_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

		case (3):
			if (v) {
				std::cout << "Performing biased chain regrowth with orientation flip..." << std::endl; 
			}
			ChainRegrowthPlusOrientationFlip_BIASED (Polymers, Cosolvent, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

		case (6):
			if (v) {
				std::cout << "Performing solvent flips..." << std::endl;
			}
			SolventFlip_UNBIASED (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break;

		case (7):
			if (v) {
				std::cout << "Performing a biased solvation shell flip..." << std::endl;
			}
			SolvationShellFlip_BIASED_remake2 (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1;
			break;

		case (8):
			if (v) {
				std::cout << "Performing a biased polymer flip..." << std::endl;
			}
			PolymerFlip_BIASED (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (4): 
			if (v) {
				std::cout << "Performing an identity swap... " << std::endl;
			}
			SolventExchange_BIASED (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z);
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;			

		case (5):
			if (v) {
				std::cout << "Performing an unbiased identity swap..." << std::endl;
			}
			SolventExchange_UNBIASED (Polymers, LATTICE, InteractionMap, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z);
			*move_number = r;
			(*attempts)[r] += 1; 
			break; 

	}

	return; 

}


// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
//			End of PerturbSystem_BIASED
// ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~


void PerturbSystem_BIASED_debug (std::vector <Polymer>* Polymers, std::vector <Particle*>* Cosolvent, std::vector <Particle*>* LATTICE, \
	std::map <std::pair <std::string, std::string>, std::tuple <std::string, double, double, int, int>>* InteractionMap, \
	std::array <double,8>* E, std::array <double,8>* contacts, std::array <int,9>* attempts, \
	bool* IMP_BOOL, bool v, double* sysEnergy, double temperature, \
	int* move_number, int x, int y, int z) {

	int index = 0; // rng_uniform (0, static_cast<int>((*Polymers).size()-1) ); 
	int r     = rng_uniform (0, 8); 

	switch (r) {

		case(0):
			if (v) {
				std::cout << "Performing end rotations..." << std::endl; 
			}
			EndRotation_UNBIASED_debug  (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break; 

		case (1):
			if (v) {
				std::cout << "Performing reptation..." << std::endl; 
			}
			Reptation_UNBIASED_debug    (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (2):
			if (v) {
				std::cout << "Performing biased chain regrowth..." << std::endl;
			}
			ChainRegrowth_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

		case (3):
			if (v) {
				std::cout << "Performing biased chain regrowth with orientation flip..." << std::endl; 
			}
			ChainRegrowthPlusOrientationFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1;			
			break; 

		case (4):
			if (v) {
				std::cout << "Performing solvent flips..." << std::endl;
			}
			SolventFlip_UNBIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r; 
			(*attempts)[r] += 1; 
			break;

		case (5):
			if (v) {
				std::cout << "Performing a biased solvation shell flip..." << std::endl;
			}
			SolvationShellFlip_BIASED_remake2_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z); 
			*move_number    = r;
			(*attempts)[r] += 1;
			break;

		case (6):
			if (v) {
				std::cout << "Performing a biased polymer flip..." << std::endl;
			}
			PolymerFlip_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, index, x, y, z);
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (7): 
			if (v) {
				std::cout << "Performing an identity swap... " << std::endl;
			}
			SolventExchange_BIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z);
			*move_number    = r;
			(*attempts)[r] += 1; 
			break;

		case (8):
			if (v) {
				std::cout << "Performing an unbiased identity swap..." << std::endl;
			}
			SolventExchange_UNBIASED_debug (Polymers, Cosolvent, LATTICE, InteractionMap, E, contacts, IMP_BOOL, sysEnergy, temperature, x, y, z); 
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
