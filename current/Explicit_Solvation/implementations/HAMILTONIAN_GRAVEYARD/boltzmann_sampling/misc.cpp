#include <iostream>
#include <vector>
#include <string>
#include <map>
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

/* ==================================================
These are some objects I have defined which I tend to use often. 
Helpful definitions which are employed often in the context of the z=6 lattice I am using.
==================================================== */ 

const std::vector <int> ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1}; 	// unit directions 
const std::vector <std::vector <int>> drns = {ex, nex, ey, ney, ez, nez};  							// vector of unit directions 

std::array <int,3> ax     = {1,0,0}, ay     = {0,1,0} , az     = {0,0,1} , nx     = {-1,0,0} ,  ny     = {0,-1,0}, nz     =  {0,0,-1} ; 
std::array <int,3> axay   = {1,1,0}, axaz   = {1,0,1} , axny   = {1,-1,0}, axnz   = {1,0,-1} , nxay    = {-1,1,0}, nxaz   =  {-1,0,1} , nxny = {-1,-1,0} , nxnz = {-1,0,-1}; 
std::array <int,3> ayaz   = {0,1,1}, aynz   = {0,1,-1}, nyaz   = {0,-1,1}, nynz   = {0,-1,-1};  
std::array <int,3> axayaz = {1,1,1}, axaynz = {1,1,-1}, axnyaz = {1,-1,1}, axnynz = {1,-1,-1},  nxayaz = {-1,1,1}, nxaynz = {-1,1,-1}, nxnyaz = {-1,-1,1}, nxnynz = {-1,-1,-1}; 
std::array <std::array <int,3>, 26> adrns = { ax, ay, az, nx, ny, nz, axay, axaz, axny, axnz, nxay, nxaz, nxny, nxnz, ayaz, aynz, nyaz, nynz, axayaz, axnyaz, axaynz, axnynz, nxayaz, nxaynz, nxnyaz, nxnynz }; 
std::map <int, std::array<double,3>> Or2Dir = { {0, {1.0,0,0}}, {1, {0,1.0,0}}, {2, {0,0,1}}, {3, {-1,0,0}}, {4, {0,-1,0}}, {5, {0,0,-1}}, {6, {1.0/(std::sqrt(2)), 1.0/(std::sqrt(2)), 0}}, {7, {1.0/(std::sqrt(2)), 0, 1.0/(std::sqrt(2))}}, \
{8, {1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}}, {9, {1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}}, {10, {-1.0/(std::sqrt(2)),1.0/(std::sqrt(2)),0}}, {11, {-1.0/(std::sqrt(2)),0,1.0/(std::sqrt(2))}}, \
{12, {-1.0/(std::sqrt(2)),-1.0/(std::sqrt(2)),0}}, {13, {-1.0/(std::sqrt(2)),0,-1.0/(std::sqrt(2))}}, {14, {0,1.0/(std::sqrt(2)),1.0/(std::sqrt(2))}}, {15, {0,1.0/(std::sqrt(2)),-1.0/(std::sqrt(2))}}, {16, {0,-1.0/(std::sqrt(2)), 1.0/(std::sqrt(2))}}, \
{17, {0,-1.0/(std::sqrt(2)), -1.0/(std::sqrt(2))}}, {18, {1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {19, {1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {20, {1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, \
{21, {1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {22, {-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, {23, {-1.0/(std::sqrt(3)),1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}}, {24, {-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),1.0/(std::sqrt(3))}}, \
{25, {-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3)),-1.0/(std::sqrt(3))}} };

//=====================================================
// impose periodic boundary conditions on vector 
//=====================================================

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


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//=====================================================

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

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//=====================================================

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


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

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



//=====================================================
std::array <std::array <int,3>, 3> HingeSwingDirections(std::array <int,3>* HingeToHinge, std::array <int,3>* HingeToKink, int x, int y, int z){
    
    (*HingeToHinge)[0] = modified_modulo((*HingeToHinge)[0], x); 
    (*HingeToHinge)[1] = modified_modulo((*HingeToHinge)[1], y); 
    (*HingeToHinge)[2] = modified_modulo((*HingeToHinge)[2], z); 

    (*HingeToKink)[0] = modified_modulo((*HingeToKink)[0], x); 
    (*HingeToKink)[1] = modified_modulo((*HingeToKink)[1], y); 
    (*HingeToKink)[2] = modified_modulo((*HingeToKink)[2], z); 


    std::array <int,3 > ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1};     // unit directions 
    std::array <std::array <int,3>,6> drns = {ex, nex, ey, ney, ez, nez};                           // vector of unit directions 


    // get rid of HingeToHinge and its negative from drns 
    std::array <int,3> nHingeToHinge = { -(*HingeToHinge).at(0), -(*HingeToHinge).at(1), -(*HingeToHinge).at(2) }; 

    std::array <std::array <int,3>, 3> directions; 
    int i {0}; 
    for (std::array<int,3>& d: drns){
    	
    	if (d == *HingeToHinge || d == nHingeToHinge || d == *HingeToKink){
    		continue; 
    	}
    	else {
    		directions[i] = d;
    		++i; 
    	}
    } 

    // drns.erase(std::remove(drns.begin(), drns.end(), *HingeToHinge), drns.end() ); 
    // drns.erase(std::remove(drns.begin(), drns.end(), nHingeToHinge), drns.end() );  
    // drns.erase(std::remove(drns.begin(), drns.end(), *HingeToKink), drns.end() ); 

    return directions; 

}

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//=====================================================
// function to print out contents of a vector 
//$====================================================

void print(std::vector <int> v){
	for (int i: v){
		std::cout << i << " | ";
	}
	std::cout << std::endl;
}

void print(std::vector <double> v){
	for (double i: v){
		std::cout << i << " | ";
	}
	std::cout << std::endl;
}

void print(std::array <int, 3> v){
	for (int i: v){
		std::cout << i << " | ";
	}
	std::cout << std::endl;
}

void print(std::array <double, 3> v){
	for (double i: v){
		std::cout << i << " | ";
	}
	std::cout << std::endl;
}

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//=====================================================
// function to print out contents of a vector of vectors 
//$====================================================

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


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//======================================================

//=====================================================
// function to add two vectors 
//$====================================================
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


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

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

double distance_between_points (std::array <int,3>* a1, std::array <int,3>* a2, int xlen, int ylen, int zlen){

	std::array <int,3> delta = subtract_arrays (a1, a2); 
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


/*double distance_between_points (std::array <double,3>* a1, std::array <double,3>* a2, int xlen, int ylen, int zlen){

	std::array <double,3> delta = subtract_arrays (a1, a2); 
	impose_pbc ( &delta, xlen, ylen, zlen ); 

	return std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);

}*/


//=====================================================
// function to subtract two vectors 
//$====================================================
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

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//=====================================================
// function to check for avoidance in the walk 
//$====================================================
bool check_avoidance(const std::vector <int> to_check, const std::vector<std::vector <int>> loc_list){
	for (std::vector <int> v: loc_list){
		if (v==to_check){
			return false;
		}
	}
	return true;
}

//=====================================================
// function to get rid of unnecesaary directions for crank shafts 
//=====================================================



//=====================================================
// executing self-avoiding random walk
//$====================================================
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
//=====================================================


//=====================================================
// executing sarw with pbc
//$====================================================
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

//===============================================================================


//=====================================================
// generating lattice points 
//$====================================================
std::vector <std::array <int,3>> create_lattice_pts(int x_len, int y_len, int z_len){
	std::vector<std::array <int,3>> lattice_pts;
	lattice_pts.reserve(x_len*y_len*z_len); 

	for (int i{0}; i < x_len; i++){
		for (int j{0}; j < y_len; j++){
			for (int k{0}; k < z_len; k++){
				lattice_pts.push_back({i,j,k}); 
			}
		}
	} 

	return lattice_pts; 

} 


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


//=====================================================
// generating neighbor lists  
//$====================================================
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


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 



// ===================================================================
// decide on whether polymer configuration is to be accepted or not
// ===================================================================
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


// ===================================================================
// extract information
// ===================================================================

//============================================================
//============================================================
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


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractNumberOfPolymers. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
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


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractContentFromFile. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

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

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

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

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 




bool isSymmetric(std::vector <std::vector <double>> mat){
	int nr = mat.size(); 
	int nc = mat.at(0).size(); 

	if (nr != nc){
		std::cerr << "ERROR: matrix is not a square, or there might be an empty line somewhere in your script." << std::endl;
		exit(EXIT_FAILURE); 
	}

	for (int i{0}; i < nr; i++){
		for (int j{i}; j < nc; j++){

			if (mat.at(i).at(j) != mat.at(j).at(i)){
				std::cerr << "matrix is not symmetric. Bad energy file." << std::endl;
				return false; 
			}
			else {
				continue; 
			}


		}
	}

	return true; 

}


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// ===============================================================


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



/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


// ===============================================================

/*
Polymer makePolymer_ZeroOrientation(std::vector <std::array <int,3> > locations, std::string type_m){
	std::vector <int> pmer_spins; 
    int size_ = locations.size(); 

    std::vector <Particle> ptc_vec; 

    for (int i=0;i<static_cast<int>( size_ ); i++ ){
        Particle p (locations.at(i), type_m, 0); 
        ptc_vec.push_back(p); 
    }

    Polymer pmer (size_, ptc_vec);

    return pmer; 
}
*/

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 



// ===============================================================


void ClusterFlip(std::vector <Particle>* cluster){


	if ((*cluster).at(0).orientation==1){
		// std::cout << "hello, is this being hit? o=1?" << std::endl;

		for (int i=0; i<static_cast<int>((*cluster).size()); i++ ){ 
			(*cluster).at(i).orientation = 0;
			// std::cout << "p.orientation is " << (*cluster).at(i).orientation << std::endl;
			// std::cout <<"--> orientation reported above should be 0" << std::endl;
		}
	}
	else {
		// std::cout << "hello, is this being hit? o=0?" << std::endl;
		for (int i=0; i<static_cast<int>((*cluster).size()); i++ ){ 
			(*cluster).at(i).orientation = 1;
			// std::cout << "p.orientation is " << (*cluster).at(i).orientation << std::endl;
			// std::cout <<"--> orientation reported above should be 1" << std::endl;
		}
		
	}

	return; 
}


/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


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

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

// ===============================================================
// ===============================================================

void StringToFile(std::string filename, std::string to_send){
	std::ofstream send_file(filename, std::ios::app); 
	send_file << to_send;
    send_file.close();

    return; 
}

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 



// ===============================================================
// ===============================================================

void InputParser(int dfreq, int max_iter, bool r,
	std::string positions, std::string topology, std::string dfile, 
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

//============================================================
//============================================================
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



//============================================================
//============================================================
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

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of checkForOverlaps. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
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
        		exit(EXIT_FAILURE);
        	}

        }
    }


    std::cout << "Input file has no overlap between and internally amongst solvent and monomers!" << std::endl;
    return true;

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of checkForOverlaps. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
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
    // std::cout << "Am i INSIDE CHECKCONNECTIVITY here?" << std::endl;
    return true;
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of checkConnectivity. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//============================================================
//============================================================
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
			std::cout << "Bad location is "; print ( location (i, x, y, z) ); 
			std::cout << "LATTICE says "; print ((*LATTICE)[i]->coords);
			std::cerr << "Something is fucked. Pointer does not correspond to position. " << std::endl;
			exit (EXIT_FAILURE);
			return false;
		}

	}

	return true;

}



//============================================================
//============================================================
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


//============================================================
//============================================================
// 
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: ParticleReporter 
//
// PARAMETERS: PolymerVector, x, y, z
//
// WHAT THE FUNCTION DOES: Calculates energy of the current Grid. Critical to correctly evolve system. 
// includes nearest neighbor interactions with and without directional effects.  
//
// DEPENDENCIES: obtain_ne_list 
//
// OPTIMIZATION OPPORTUNITY: I am double counting monomer-monomer interactions. This can possibly be avoided. 
//
// THE CODE: 
/*

*/

//============================================================
//============================================================
// 
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: MonomerReporter
//
// PARAMETERS: PolymerVector, to_check
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

	if ( (*LATTICE)[lattice_index((*check_1), y, z)]->ptype[0] == 'm' || (*LATTICE)[lattice_index((*check_2), y, z)]->ptype [0]== 'm'){
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


//============================================================
//============================================================
// 
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: 
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

double CalculateEnergy (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z, std::array<double,8>* E, std::array <double,8>* contacts) {
    
    double Energy {0.0};
    (*contacts) = {0,0,0,0,0,0,0,0};
    // polymer-polymer interaction energies 
    double delta_p = 0; 
    double delta_o = 0; 
    std::array <double,3> ext1, ext2, scaled_o1, scaled_o2;
    std::array <std::array <int,3>, 26> ne_list;
    std::vector <int> ss1; 

    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){
            ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list 
            
            // std::cout << "Particle loc is "; print (p->coords); 

            for ( std::array <int, 3>& loc: ne_list){

            	if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype[0] == 'm'){
					
            		delta_p   = distance_between_points ( &(p->coords), &((*LATTICE)[lattice_index(loc, y, z)]->coords ), x, y, z );
            		scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[p->orientation] ) );
            		scaled_o2 = scale_arrays (delta_p/2, &(Or2Dir[(*LATTICE)[ lattice_index(loc, y, z )]->orientation ] ) );

            		ext1  = add_arrays  ( &(p->coords), &( scaled_o1 ) );
            		ext2  = add_arrays ( &( (*LATTICE)[ lattice_index(loc, y, z )]->coords ), &( scaled_o2 ) );

            		/*
            		std::cout << "or = " << p->orientation << ", "; print (p->coords); 
            		std::cout << "scaled o1 = "; print (scaled_o1); 
            		std::cout << "ext1 is "; print (ext1); 

            		std::cout << "or = " << (*LATTICE)[lattice_index(loc, y, z)]->orientation << ", "; print ((*LATTICE)[lattice_index(loc, y, z)]->coords); 
            		std::cout << "scaled o2 = "; print (scaled_o2); 
            		std::cout << "ext2 is "; print (ext2);
					*/

            		delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 
            		// std::cout << "delta_o = " << delta_o << std::endl; 

            		if ( delta_o <= delta_p/2 ){
            			Energy          += 0.5* (*E)[0];
            			(*contacts)[0]   += 0.5;
            		}
            		else {
            			Energy          += 0.5* (*E)[1];
            			(*contacts)[1]  += 0.5;
            		}
            	}

            	else if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s1" ) {

            		if ( std::find (ss1.begin(), ss1.end(), lattice_index (loc, y, z)) == ss1.end() ){
            			ss1.push_back ( lattice_index (loc, y, z) ); 
            		}

            		delta_p   = distance_between_points ( &(p->coords), &((*LATTICE)[lattice_index(loc, y, z)]->coords ), x, y, z );
            		scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[p->orientation] ) );
            		scaled_o2 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[ lattice_index(loc, y, z )]->orientation ] ) );

            		ext1  = add_arrays ( &(p->coords), &( scaled_o1 ) );
            		ext2  = add_arrays ( &( (*LATTICE)[ lattice_index(loc, y, z )]->coords ), &( scaled_o2 ) );

            		delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 

            		if ( delta_o <= delta_p/2 ) {
            			Energy         += (*E)[2];
            			(*contacts)[2] += 1;
            		}
            		else {
            			Energy         += (*E)[3]; 
            			(*contacts)[3]  += 1;
            		}
            	}

            	else if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s2" ) {

            		if ( std::find (ss1.begin(), ss1.end(), lattice_index (loc, y, z)) == ss1.end() ){
            			ss1.push_back ( lattice_index (loc, y, z) ); 
            		}

            		delta_p   = distance_between_points ( &(p->coords), &((*LATTICE)[lattice_index(loc, y, z)]->coords ), x, y, z );
            		scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[p->orientation] ) );
            		scaled_o2 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[ lattice_index(loc, y, z )]->orientation ] ) );

            		ext1  = add_arrays ( &(p->coords), &( scaled_o1 ) );
            		ext2  = add_arrays ( &( (*LATTICE)[ lattice_index(loc, y, z )]->coords ), &( scaled_o2 ) );

            		delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 

            		if ( delta_o <= delta_p/2 ) {
            			Energy         += (*E)[4];
            			(*contacts)[4] += 1;
            		}
            		else {
            			Energy         += (*E)[5]; 
            			(*contacts)[5]  += 1;
            		}
            	}
            }
        }
    }

    // std::cout << "Size of solvation shell is " << ss1.size() << "." << std::endl;
    // std::vector <int> ss2; 
    // find the second solvation shell... 
    // ss1 is the first solvation shell... 
    /*
    for (int i: ss1){

    	ne_list = obtain_ne_list( (*LATTICE)[i]->coords, x, y, z); 
    	for ( std::array<int,3>& n: ne_list ) {
    		// only calculate energy if the two particles of different solvent species 

    		if ( (*LATTICE)[ lattice_index (n, y, z)]->ptype[0]=='m' ){
    			// this energy has already been calculated
    			continue;
    		}
    		else if ( (*LATTICE)[lattice_index (n, y, z)]->ptype != (*LATTICE)[i]->ptype ){

    			delta_p   = distance_between_points ( &(n), &((*LATTICE)[i]->coords ), x, y, z );
            	scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[lattice_index(n, y, z)]->orientation ] ) );
            	scaled_o2 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[i]->orientation ] ) );

            	ext1  = add_arrays ( &(n), &( scaled_o1 ) );
            	ext2  = add_arrays ( &( (*LATTICE)[i]->coords ), &( scaled_o2 ) );

            	delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 

            	if ( std::find (ss1.begin(), ss1.end(), lattice_index(n, y, z) ) != ss1.end() ){

	            	if ( delta_o <= delta_p/2 ) {
	            		Energy         += 0.5*(*E)[6];
	            		(*contacts)[6] += 1;
	            	}
	            	else {
	            		Energy          += 0.5*(*E)[7]; 
	            		(*contacts)[7]  += 1;
	            	}
	            }
	            else {
	            	
	            	if ( delta_o <= delta_p/2 ) {
	            		Energy         += (*E)[6];
	            		(*contacts)[6] += 1;
	            	}
	            	else {
	            		Energy          += (*E)[7]; 
	            		(*contacts)[7]  += 1;
	            	}	
	            }
    		}
    	}
    }
	*/ 
    /*
    std::cout << "Printing out solvation shell one..." << std::endl;
    for (int i: ss1){
    	print (location(i, x, y, z));
    }
    std::cout << "Size of solvation shell is " << ss1.size() << std::endl << std::endl; 
    
    std::cout << "Printing out solvation shell two..." << std::endl;
    for (int i: ss2){
    	print (location (i, x,y, z));
    }
    std::cout << "Size of solvation shell is " << ss2.size() << std::endl << std::endl; 
	*/ 

    return Energy; 
}


double CalculateEnergy (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, int x, int y, int z, std::array<double,8>* E) {
    
    double Energy {0.0};
    // (*contacts) = {0,0,0,0,0,0,0,0};
    // polymer-polymer interaction energies 
    double delta_p = 0; 
    double delta_o = 0; 
    std::array <double,3> ext1, ext2, scaled_o1, scaled_o2;
    std::array <std::array <int,3>, 26> ne_list;
    std::vector <int> ss1; 

    for (Polymer& pmer: (*Polymers)) {
        for (Particle*& p: pmer.chain){
            ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list 
            
            // std::cout << "Particle loc is "; print (p->coords); 

            for ( std::array <int, 3>& loc: ne_list){

            	if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype[0] == 'm'){
					
            		delta_p   = distance_between_points ( &(p->coords), &((*LATTICE)[lattice_index(loc, y, z)]->coords ), x, y, z );
            		scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[p->orientation] ) );
            		scaled_o2 = scale_arrays (delta_p/2, &(Or2Dir[(*LATTICE)[ lattice_index(loc, y, z )]->orientation ] ) );

            		ext1  = add_arrays  ( &(p->coords), &( scaled_o1 ) );
            		ext2  = add_arrays ( &( (*LATTICE)[ lattice_index(loc, y, z )]->coords ), &( scaled_o2 ) );

            		delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 

            		if ( delta_o <= delta_p/2 ){
            			Energy          += 0.5* (*E)[0];
            			
            		}
            		else {
            			Energy          += 0.5* (*E)[1];
            			
            		}
            	}

            	else if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s1" ) {

            		if ( std::find (ss1.begin(), ss1.end(), lattice_index (loc, y, z)) == ss1.end() ){
            			ss1.push_back ( lattice_index (loc, y, z) ); 
            		}

            		delta_p   = distance_between_points ( &(p->coords), &((*LATTICE)[lattice_index(loc, y, z)]->coords ), x, y, z );
            		scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[p->orientation] ) );
            		scaled_o2 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[ lattice_index(loc, y, z )]->orientation ] ) );

            		ext1  = add_arrays ( &(p->coords), &( scaled_o1 ) );
            		ext2  = add_arrays ( &( (*LATTICE)[ lattice_index(loc, y, z )]->coords ), &( scaled_o2 ) );

            		delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 

            		if ( delta_o <= delta_p/2 ) {
            			Energy         += (*E)[2];
            			
            		}
            		else {
            			Energy         += (*E)[3]; 
            			
            		}
            	}

            	else if ( (*LATTICE)[ lattice_index(loc, y, z) ]->ptype == "s2" ) {

            		if ( std::find (ss1.begin(), ss1.end(), lattice_index (loc, y, z)) == ss1.end() ){
            			ss1.push_back ( lattice_index (loc, y, z) ); 
            		}

            		delta_p   = distance_between_points ( &(p->coords), &((*LATTICE)[lattice_index(loc, y, z)]->coords ), x, y, z );
            		scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[p->orientation] ) );
            		scaled_o2 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[ lattice_index(loc, y, z )]->orientation ] ) );

            		ext1  = add_arrays ( &(p->coords), &( scaled_o1 ) );
            		ext2  = add_arrays ( &( (*LATTICE)[ lattice_index(loc, y, z )]->coords ), &( scaled_o2 ) );

            		delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 

            		if ( delta_o <= delta_p/2 ) {
            			Energy         += (*E)[4];
            			
            		}
            		else {
            			Energy         += (*E)[5]; 
            			
            		}
            	}
            }
        }
    }

    // std::cout << "Size of solvation shell is " << ss1.size() << "." << std::endl;
    // std::vector <int> ss2; 
    // find the second solvation shell... 
    // ss1 is the first solvation shell... 
    /*
    for (int i: ss1){

    	ne_list = obtain_ne_list( (*LATTICE)[i]->coords, x, y, z); 
    	for ( std::array<int,3>& n: ne_list ) {
    		// only calculate energy if the two particles of different solvent species 

    		if ( (*LATTICE)[ lattice_index (n, y, z)]->ptype[0]=='m' ){
    			// this energy has already been calculated
    			continue;
    		}
    		else if ( (*LATTICE)[lattice_index (n, y, z)]->ptype != (*LATTICE)[i]->ptype ){

    			delta_p   = distance_between_points ( &(n), &((*LATTICE)[i]->coords ), x, y, z );
            	scaled_o1 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[lattice_index(n, y, z)]->orientation ] ) );
            	scaled_o2 = scale_arrays (delta_p/2, &( Or2Dir[(*LATTICE)[i]->orientation ] ) );

            	ext1  = add_arrays ( &(n), &( scaled_o1 ) );
            	ext2  = add_arrays ( &( (*LATTICE)[i]->coords ), &( scaled_o2 ) );

            	delta_o = distance_between_points ( &ext1, &ext2, x, y, z ); 

            	if ( std::find (ss1.begin(), ss1.end(), lattice_index(n, y, z) ) != ss1.end() ){

	            	if ( delta_o <= delta_p/2 ) {
	            		Energy         += 0.5*(*E)[6];
	            		
	            	}
	            	else {
	            		Energy          += 0.5*(*E)[7]; 
	            		
	            	}
	            }
	            else {
	            	
	            	if ( delta_o <= delta_p/2 ) {
	            		Energy         += (*E)[6];
	            		
	            	}
	            	else {
	            		Energy          += (*E)[7]; 
	            		
	            	}	
	            }
    		}
    	}
    }
	*/ 
    /*
    std::cout << "Printing out solvation shell one..." << std::endl;
    for (int i: ss1){
    	print (location(i, x, y, z));
    }
    std::cout << "Size of solvation shell is " << ss1.size() << std::endl << std::endl; 
    
    std::cout << "Printing out solvation shell two..." << std::endl;
    for (int i: ss2){
    	print (location (i, x,y, z));
    }
    std::cout << "Size of solvation shell is " << ss2.size() << std::endl << std::endl; 
	*/ 

    return Energy; 
}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of calculateEnergy. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
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
    // const auto& str = os.str(); 
    // dump_file.write(str.c_str(), static_cast<std::streamsize> (str.size() )); 

    // dump_file.close();

    

    return; 

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of dumpPositionsOfPolymer. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


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


void dumpEnergy (double sysEnergy, int step, std::array <double,8>* contacts , std::string filename){
    std::ofstream dump_file(filename, std::ios::app); 
    // std::ostringstream os; 
    
    dump_file << sysEnergy << " | " << (*contacts)[0]+(*contacts)[1] << " | " << (*contacts)[0] << " | " << (*contacts)[1] << " | " \
            << (*contacts)[2]+(*contacts)[3]+(*contacts)[4]+(*contacts)[5] << " | " << (*contacts)[2] + (*contacts)[3] << " | " << (*contacts)[2] << " | " << (*contacts)[3] << " | " \
            << (*contacts)[4]+(*contacts)[5] << " | " << (*contacts)[4] << " | " << (*contacts)[5] << " | " << (*contacts)[6] + (*contacts)[7] << " | " << (*contacts)[6] << " | " << (*contacts)[7] << " | " << step << "\n";
    
    return; 
}
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

void dumpMoveStatistics (std::array <int,9>* attempts, std::array <int,9>* acceptances, int step, std::string stats_file){
    
    std::ofstream dump_file (stats_file, std::ios::out); 
    dump_file << "For step " << step << ".\n";
    
    dump_file << "End rotations                      - attempts: " << (*attempts)[0] <<", acceptances: " << (*acceptances)[0] << ", acceptance fraction: " << static_cast<double>((*acceptances)[0])/static_cast<double>((*attempts)[0]) << ".\n"; 
    // dump_file << "Bond vibrations                    - attempts: " << (*attempts)[1] <<", acceptances: " << (*acceptances)[1] << ", acceptance fraction: " << static_cast<double>((*acceptances)[1])/static_cast<double>((*attempts)[1]) << ".\n"; 
    // dump_file << "Crank shafts                       - attempts: " << (*attempts)[2] <<", acceptances: " << (*acceptances)[2] << ", acceptance fraction: " << static_cast<double>((*acceptances)[2])/static_cast<double>((*attempts)[2]) << ".\n"; 
    dump_file << "Reptation                          - attempts: " << (*attempts)[3] <<", acceptances: " << (*acceptances)[3] << ", acceptance fraction: " << static_cast<double>((*acceptances)[3])/static_cast<double>((*attempts)[3]) << ".\n"; 
    dump_file << "Chain regrowth                     - attempts: " << (*attempts)[4] <<", acceptances: " << (*acceptances)[4] << ", acceptance fraction: " << static_cast<double>((*acceptances)[4])/static_cast<double>((*attempts)[4]) << ".\n"; 
    dump_file << "Single solvent orientation flips   - attempts: " << (*attempts)[5] <<", acceptances: " << (*acceptances)[5] << ", acceptance fraction: " << static_cast<double>((*acceptances)[5])/static_cast<double>((*attempts)[5]) << ".\n"; 
    dump_file << "Single monomer orientation flips   - attempts: " << (*attempts)[6] <<", acceptances: " << (*acceptances)[6] << ", acceptance fraction: " << static_cast<double>((*acceptances)[6])/static_cast<double>((*attempts)[6]) << ".\n"; 
    dump_file << "Single random site flips           - attempts: " << (*attempts)[7] <<", acceptances: " << (*acceptances)[7] << ", acceptance fraction: " << static_cast<double>((*acceptances)[7])/static_cast<double>((*attempts)[7]) << ".\n"; 
    dump_file << "Single solvent exchange flip       - attempts: " << (*attempts)[8] <<", acceptances: " << (*acceptances)[8] << ", acceptance fraction: " << static_cast<double>((*acceptances)[8])/static_cast<double>((*attempts)[8]) << ".\n"; 
    // dump_file << "Multiple solvent orientation flips - attempts: " << (*attempts)[7] <<", acceptances: " << (*acceptances)[7] << ", acceptance fraction: " << static_cast<double>((*acceptances)[7])/static_cast<double>((*attempts)[7]) << ".\n"; 
    // dump_file << "Multiple monomer orientation flips - attempts: " << (*attempts)[8] <<", acceptances: " << (*acceptances)[8] << ", acceptance fraction: " << static_cast<double>((*acceptances)[8])/static_cast<double>((*attempts)[8]) << ".\n"; 

    return; 

}

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


void dumpLATTICE ( std::vector <Particle*> *LATTICE, int step, int y, int z, std::string filename ){

	std::ofstream dump_file ( filename, std::ios::out ); 
	dump_file << "FINAL STEP: " << step << ".\n"; 
	for ( Particle*& p: (*LATTICE) ){
		dump_file << p->orientation << ", " << p->ptype << ", " << lattice_index(p->coords, y, z) << "\n"; 
	}

	dump_file << "END. \n";
	return; 

}

//============================================================
//============================================================
// 
// NAME OF FUNCTION: TailRotation
//
// PARAMETERS: index of a polymer to perform ZeroIndexRotation, the polymervector, and the dimensions of the box 
// 
// WHAT THE FUNCTION DOES: it will perform a rotation of the monomer at the first index (the tail)
// of the polymer. As it stands, this function only rotates one molecule at the tail.
//
// PLANNED EXTENSION: Multiple molecules at the tail need to be rotated.    
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully.    
//
// THE CODE: 

void TailRotation (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, double* sysEnergy, double temperature, \
	bool* IMP_BOOL, std::array<double,8>* E, std::array<double,8>* contacts ){

    // get the neighborlist of particle at index 1 

	std::array <double,8> c_contacts = *contacts; 

    std::array <int,3> loc_0 = (*Polymers)[index].chain[0]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[1]->coords; 

    double w1 {0}, w2 {0}; // acceptance criteria weights 

    // first check if tail rotation can even be performed in the first place... 
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

	double rweight = 0.0; 

	// find all spots that are available 
	std::vector <std::array <int,3>> idx_v; 
	for (std::array <int,3>& to_rot: ne_list){
		if ( to_rot == loc_1){
			continue; 
		}

		else if ( to_rot == (*Polymers)[index].chain[2]->coords ) {
			continue; 
		}

		// else if ( ! ( MonomerReporter(Polymers, &to_rot) ) ){
		else if ( ! (MonomerReporter(LATTICE, &to_rot, y, z) ) ){
			rweight += 1;
			idx_v.push_back( to_rot );  
		}
	}

    if ( rweight == 0 ){
    	// if nothing can be done, just get out. 
    	*IMP_BOOL = false;
    	return; 
    }

    // std::cout << "Suggested locations are "; print (idx_v);

    // prepare the initial boltzmann sampling
    // get the solvent neighbors in the current location 
    std::array <std::array<int,3>,26> initial_ne_list = obtain_ne_list (loc_0, x, y, z);
    std::array<int,3> solvents_to_flip;
    std::array<int,3> solvent_init_config; 
    int current_idx {0};
    int solvent_count = 0; 
    for (std::array<int,3>& i: initial_ne_list){
    	if ( (*LATTICE)[lattice_index(i, y, z)]->ptype[0] != 'm') {
    		solvent_count += 1; 
    		if (current_idx < 3 ){
	    		solvents_to_flip[current_idx] = lattice_index(i,y,z); 
	    		solvent_init_config [current_idx] = (*LATTICE)[lattice_index(i, y, z)]->orientation; 
	    		current_idx += 1; 
    		}
    	}
    }
    std::cout << "solvent count = " << solvent_count << std::endl; 
    std::array<std::array<int,3>,51> orientation_configs; 
    std::array<double,51> energetic_states; 
    // std::cout << "solvents present = " << solvent_count << std::endl;
    // std::cout << "indices to flip are "; print (solvents_to_flip);

    orientation_configs[0] = solvent_init_config; 
    energetic_states   [0] = *sysEnergy;

    std::cout << "Initial orientations is: "; print (orientation_configs[0]);

    int or_config_c {1}; 
    for (int i{0}; i<5; ++i){
    	(*LATTICE)[ solvents_to_flip[0] ]->orientation = rng_uniform (0, 25); 
    	for (int j{0}; j<5; ++j){
    		(*LATTICE)[ solvents_to_flip[1] ]->orientation = rng_uniform (0, 25); 
    		for (int k{0}; k<2; ++k){
    			(*LATTICE)[ solvents_to_flip[2] ]->orientation = rng_uniform (0, 25); 
    			orientation_configs[or_config_c] = {(*LATTICE)[ solvents_to_flip[0] ]->orientation, (*LATTICE)[ solvents_to_flip[1] ]->orientation, (*LATTICE)[ solvents_to_flip[2] ]->orientation}; 
    			// std::cout << "Flipped orientations is: "; print (orientation_configs[or_config_c]);

    			energetic_states[or_config_c] = CalculateEnergy (Polymers, LATTICE, x, y, z, E); 
    			// std::cout << "energy is = " << energetic_states[or_config_c] << std::endl; 
    			or_config_c += 1; 
    		}
    	}
    }

    for ( double energy: energetic_states){
    	w1 += std::exp (-1/temperature*energy); 
    }
    w1 = w1*solvent_count/26; 

    std::cout << "Orientation of monomer 0 is " << (*Polymers)[0].chain[0]->orientation << std::endl;
    // print (energetic_states); 

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // i have done the boltzmann sampling spins in the original state. 
    // let's move on to the perturbed state - after tailrotation has been performed. 
    // but first, get the LATTICE to its original state 

    // std::cout << "Index: orientation --- ";
    for (int i{0}; i<3; ++i){
    	(*LATTICE)[solvents_to_flip[i]]->orientation = solvent_init_config[i]; 
    	// std::cout << solvents_to_flip[i] << " : " << solvent_init_config[i] << " --- "; 
    }
    std::cout << std::endl;

    // lattice is now in the original state... 
    //************************************************

    // otherwise, let's fuckin go

	int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 

	// made the change to the polymer
	(*Polymers)[index].chain[0]->coords = idx_v[r]; 

	std::cout << "=--------------------------------------------------------------------=" << std::endl;
	std::cout << "Suggested location for monomer to swing to is "; print (idx_v[r]);

	// make the change on the lattice 
	// make the change only to the SOLVENT site... 
	(*LATTICE)[ lattice_index (idx_v[r], y, z) ]->coords = loc_0; 
	std::cout << "Just to check, loc_0 is "; print (loc_0);

	// do the switch 
	// take the pointer of the solvent, and put it where the monomer was on the lattice 
	(*LATTICE)[ lattice_index (loc_0, y, z) ]    = (*LATTICE)[ lattice_index (idx_v[r], y, z)]; 
	(*LATTICE)[ lattice_index (idx_v[r], y, z)] = (*Polymers)[index].chain[0]; 

	// update memory 
	// (*memory).first.push_back  ( loc_0 ); 		// initial location of monomer
	// (*memory).second.push_back ( idx_v[r] );	// final location of monomer 



	// update connmap 
	(*Polymers)[index].ChainToConnectivityMap(); 
	

	std::cout << "Orientation of monomer 0 is " << (*Polymers)[0].chain[0]->orientation << std::endl; 

	// run the boltzmann sampling of the new state... 
	// get the solvent neighbors of the particle at the very end.  
	// reusing the same variable from above for no particular reason apart from efficient resource utilization... 
	std::array <std::array<double,8>, 51> store_contacts; 

	initial_ne_list     = obtain_ne_list (idx_v[r], x, y, z);
	energetic_states[0] = CalculateEnergy (Polymers, LATTICE, x, y, z, E, contacts);
	store_contacts [0] = *contacts; 
	
	current_idx = 0; 

	// obtain the first iteration of the perturbed state 
	solvent_count = 0; 

	std::array <int,3> solvents_to_flip_n; 
	std::array <int,3> solvent_init_config_n; 
	for (std::array<int,3>& i: initial_ne_list){
    	if ( (*LATTICE)[lattice_index(i, y, z)]->ptype[0] != 'm') {
    		solvent_count += 1; 
    		if (current_idx < 3){
    			solvents_to_flip_n[current_idx] = lattice_index(i,y,z); 
    			solvent_init_config_n [current_idx] = (*LATTICE)[lattice_index(i, y, z)]->orientation; 
    			current_idx += 1; 
    		}
    	}
    }
    orientation_configs[0] = {(*LATTICE)[ solvents_to_flip_n[0] ]->orientation, (*LATTICE)[ solvents_to_flip_n[1] ]->orientation, (*LATTICE)[ solvents_to_flip_n[2] ]->orientation}; 
    std::cout << "solvent count = " << solvent_count << std::endl;

    or_config_c = 1; 
    
    std::cout << "Orientation of monomer 0 is " << (*Polymers)[0].chain[0]->orientation << std::endl;

    for (int i{0}; i<5; ++i){
    	(*LATTICE)[ solvents_to_flip_n[0] ]->orientation = rng_uniform (0, 25); 
    	for (int j{0}; j<5; ++j){
    		(*LATTICE)[ solvents_to_flip_n[1] ]->orientation = rng_uniform (0, 25); 
    		for (int k{0}; k<2; ++k){
    			(*LATTICE)[ solvents_to_flip_n[2] ]->orientation = rng_uniform (0, 25); 
    			orientation_configs[or_config_c] = {(*LATTICE)[ solvents_to_flip_n[0] ]->orientation, (*LATTICE)[ solvents_to_flip_n[1] ]->orientation, (*LATTICE)[ solvents_to_flip_n   [2] ]->orientation}; 
    			energetic_states[or_config_c] = CalculateEnergy (Polymers, LATTICE, x, y, z, E, contacts); 
    			store_contacts  [or_config_c] = *contacts; 
    			or_config_c += 1; 
    		}
    	}
    }

    std::cout << "Orientation of monomer 0 is " << (*Polymers)[0].chain[0]->orientation << std::endl;
    // print(energetic_states);
    // calculate the boltzmann weight 
    for ( double energy: energetic_states){
    	w2 += std::exp (-1/temperature*energy); 
    }
    w2 = w2*solvent_count/26; 
    
    // choose an orientation of solvent molecules from the states we have 

    double rng = rng_uniform (0.0, 1.0); 
    std::cout << "rng = " << rng << std::endl; 

    if (rng < w2/w1){
    	std::cout << "w2 = " << w2 <<", w1 = " << w1 << std::endl;
    	std::cout << "We accepted." << std::endl;
    	double state_index = rng_uniform(0.0, 1.0); 
    	double running_boltzmann = 0; 
    	int count = 0; 
    	for (double energy: energetic_states){
    		running_boltzmann += std::exp (-1/temperature*energy)/w2; 
    		if ( state_index <= running_boltzmann ) {
    			break;
    		}
    		count += 1; 
    	}

    	for (int i{0}; i<3; ++i){
    		(*LATTICE)[solvents_to_flip_n[i]]->orientation = orientation_configs[count][i];
    	}
    	*contacts = store_contacts  [count]; 
    	std::cout << "contacts on acceptance is "; print(*contacts);
    	*sysEnergy = energetic_states[count];
    	std::cout << "New energy is " << *sysEnergy << std::endl;

    }

    else {
    	std::cout << "We rejected." << std::endl;
    	// return back to the original state 
    	(*Polymers)[index].chain[0]->coords = loc_0;
    	(*LATTICE)[ lattice_index(loc_0, y, z) ]->coords = idx_v[r];

    	// do the switch 
    	(*LATTICE)[ lattice_index( idx_v[r], y, z) ] = (*LATTICE)[ lattice_index (loc_0, y, z) ];
    	(*LATTICE)[ lattice_index( loc_0, y, z )   ] = (*Polymers)[index].chain[0]; 

    	// now that the solvent and polymer are in their original spatial state, i have to now get the solvents back into
    	// their original orientational state 
    	for (int j {0}; j<3; ++j) {
    		(*LATTICE)[solvents_to_flip[j]]->orientation = solvent_init_config [j];
    	}
    	std::cout << "contacts on rejection is "; print(c_contacts);
    	*contacts = c_contacts;
    }
	
    std::cout << "Orientation of monomer 0 is " << (*Polymers)[0].chain[0]->orientation << std::endl;

	return; 
}



//============================================================
//============================================================
// 
// NAME OF FUNCTION: HeadRotation
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

void HeadRotation (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, double* sysEnergy, double temperature, \
	bool* IMP_BOOL, std::array<double,8>* E, std::array<double,8>* contacts){ 

	// make a copy of the contacts 
	std::array<double,8> c_contacts = *contacts; 

    // get the neighborlist of particle at index 1 
    int dop = (*Polymers)[index].deg_poly; 

    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1]->coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2]->coords;
    
    double w1 {0}, w2 {0}; 

    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

    double rweight = 0; 
	
    // find all spots that are available 
    std::vector <std::array <int,3>> idx_v; 
    for (std::array <int,3>& to_rot: ne_list){

    	if (to_rot == loc_1){
    		continue; 
    	}

    	else if ( to_rot == (*Polymers)[index].chain[dop-3]->coords ) {
			continue; 
		}
    	// else if ( ! (MonomerReporter(Polymers, &to_rot) ) ){
		else if ( ! ( MonomerReporter (LATTICE, &to_rot, y, z ) ) ) {
    		// std::cout << "a free location is "; print(to_rot); 
    		rweight += 1; 
    		idx_v.push_back(to_rot); 
    	}
    }

    if ( rweight == 0){
    	// if rweight is zero, return 
    	*IMP_BOOL = false;
    	return; 
    }

    // prepare the initial boltzmann sampling 
    // get the solvent neighbors in the current location 

    std::array <std::array<int,3>,26> initial_ne_list = obtain_ne_list (loc_0, x, y, z);
    std::array<int,3> solvents_to_flip;
    std::array<int,3> solvent_init_config; 
    int current_idx {0};
    int solvent_count = 0; 
    for (std::array<int,3>& i: initial_ne_list){
    	if ( (*LATTICE)[lattice_index(i, y, z)]->ptype[0] != 'm') {
    		solvent_count += 1; 
    		if (current_idx < 3 ){
	    		solvents_to_flip[current_idx] = lattice_index(i,y,z); 
	    		solvent_init_config [current_idx] = (*LATTICE)[lattice_index(i, y, z)]->orientation; 
	    		current_idx += 1; 
    		}
    	}
    }
    std::cout << "solvent count = " << solvent_count << std::endl; 
    std::array<std::array<int,3>,51> orientation_configs; 
    std::array<double,51> energetic_states; 

    orientation_configs[0] = solvent_init_config; 
    energetic_states   [0] = *sysEnergy;

    std::cout << "Initial orientations is: "; print (orientation_configs[0]);

    int or_config_c {1}; 
    for (int i{0}; i<5; ++i){
    	(*LATTICE)[ solvents_to_flip[0] ]->orientation = rng_uniform (0, 25); 
    	for (int j{0}; j<5; ++j){
    		(*LATTICE)[ solvents_to_flip[1] ]->orientation = rng_uniform (0, 25); 
    		for (int k{0}; k<2; ++k){
    			(*LATTICE)[ solvents_to_flip[2] ]->orientation = rng_uniform (0, 25); 
    			orientation_configs[or_config_c] = {(*LATTICE)[ solvents_to_flip[0] ]->orientation, (*LATTICE)[ solvents_to_flip[1] ]->orientation, (*LATTICE)[ solvents_to_flip[2] ]->orientation}; 
    			// std::cout << "Flipped orientations is: "; print (orientation_configs[or_config_c]);

    			energetic_states[or_config_c] = CalculateEnergy (Polymers, LATTICE, x, y, z, E); 
    			// std::cout << "energy is = " << energetic_states[or_config_c] << std::endl; 
    			or_config_c += 1; 
    		}
    	}
    }

    for ( double energy: energetic_states){
    	w1 += std::exp (-1/temperature*energy); 
    }
    w1 = w1*solvent_count/26; 

    std::cout << "Orientation of monomer 0 is " << (*Polymers)[index].chain[dop-1]->orientation << std::endl;
    // print (energetic_states); 

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // i have done the boltzmann sampling spins in the original state. 
    // let's move on to the perturbed state - after tailrotation has been performed. 
    // but first, get the LATTICE to its original state 

    // std::cout << "Index: orientation --- ";
    for (int i{0}; i<3; ++i){
    	(*LATTICE)[solvents_to_flip[i]]->orientation = solvent_init_config[i]; 
    	// std::cout << solvents_to_flip[i] << " : " << solvent_init_config[i] << " --- "; 
    }
    std::cout << std::endl;

    // lattice is now in the original state... 
    //************************************************

    // let's fuckin go

    int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 
	
	// make the change to the polymer 
	(*Polymers)[index].chain[dop-1]->coords = idx_v[r];

	std::cout << "=--------------------------------------------------------------------=" << std::endl;
	std::cout << "Suggested location for monomer to swing to is "; print (idx_v[r]);

	// make the change on the lattice 
	// make the change only to the SOLVENT site...
	(*LATTICE)[ lattice_index (idx_v[r], y, z)]->coords = loc_0; 

	// do the switch 
	// take the pointer of the solvent, and put it where the monomer was on the lattice 
	(*LATTICE)[ lattice_index (loc_0, y, z)]    = (*LATTICE)[ lattice_index (idx_v[r], y, z) ]; 
	(*LATTICE)[ lattice_index (idx_v[r], y, z)] = (*Polymers)[index].chain[dop-1]; 

	// update memory 
	// (*memory).first.push_back  ( loc_0 ); 
	// (*memory).second.push_back ( idx_v[r] );

	// update connmap 
	(*Polymers)[index].ChainToConnectivityMap(); 
	// (*rweight) = (*rweight)/26.0; 

	std::cout << "Orientation of monomer dop-1 is " << (*Polymers)[index].chain[dop-1]->orientation << std::endl; 

	// run the boltzmann sampling of the new state... 
	// get the solvent neighbors of the particle at the very end.  
	// reusing the same variable from above for no particular reason apart from efficient resource utilization... 
	std::array <std::array<double,8>, 51> store_contacts; 

	initial_ne_list     = obtain_ne_list (idx_v[r], x, y, z);
	energetic_states[0] = CalculateEnergy (Polymers, LATTICE, x, y, z, E, contacts);
	store_contacts [0] = *contacts; 
	
	current_idx = 0; 

	// obtain the first iteration of the perturbed state 
	solvent_count = 0; 

	std::array <int,3> solvents_to_flip_n; 
	std::array <int,3> solvent_init_config_n; 
	for (std::array<int,3>& i: initial_ne_list){
    	if ( (*LATTICE)[lattice_index(i, y, z)]->ptype[0] != 'm') {
    		solvent_count += 1; 
    		if (current_idx < 3){
    			solvents_to_flip_n[current_idx] = lattice_index(i,y,z); 
    			solvent_init_config_n [current_idx] = (*LATTICE)[lattice_index(i, y, z)]->orientation; 
    			current_idx += 1; 
    		}
    	}
    }

	orientation_configs[0] = {(*LATTICE)[ solvents_to_flip_n[0] ]->orientation, (*LATTICE)[ solvents_to_flip_n[1] ]->orientation, (*LATTICE)[ solvents_to_flip_n[2] ]->orientation}; 
    std::cout << "solvent count = " << solvent_count << std::endl;

    or_config_c = 1; 
    
    std::cout << "Orientation of monomer 0 is " << (*Polymers)[0].chain[dop-1]->orientation << std::endl;

    for (int i{0}; i<5; ++i){
    	(*LATTICE)[ solvents_to_flip_n[0] ]->orientation = rng_uniform (0, 25); 
    	for (int j{0}; j<5; ++j){
    		(*LATTICE)[ solvents_to_flip_n[1] ]->orientation = rng_uniform (0, 25); 
    		for (int k{0}; k<2; ++k){
    			(*LATTICE)[ solvents_to_flip_n[2] ]->orientation = rng_uniform (0, 25); 
    			orientation_configs[or_config_c] = {(*LATTICE)[ solvents_to_flip_n[0] ]->orientation, (*LATTICE)[ solvents_to_flip_n[1] ]->orientation, (*LATTICE)[ solvents_to_flip_n   [2] ]->orientation}; 
    			energetic_states[or_config_c] = CalculateEnergy (Polymers, LATTICE, x, y, z, E, contacts); 
    			store_contacts  [or_config_c] = *contacts; 
    			or_config_c += 1; 
    		}
    	}
    }

    std::cout << "Orientation of monomer 0 is " << (*Polymers)[index].chain[dop-1]->orientation << std::endl;
    // print(energetic_states);
    // calculate the boltzmann weight 
    for ( double energy: energetic_states){
    	w2 += std::exp (-1/temperature*energy); 
    }
    w2 = w2*solvent_count/26; 
    
    // choose an orientation of solvent molecules from the states we have 

    double rng = rng_uniform (0.0, 1.0); 
    std::cout << "rng = " << rng << std::endl; 

    if (rng < w2/w1){
    	std::cout << "w2 = " << w2 <<", w1 = " << w1 << std::endl;
    	std::cout << "We accepted." << std::endl;
    	double state_index = rng_uniform(0.0, 1.0); 
    	double running_boltzmann = 0; 
    	int count = 0; 
    	for (double energy: energetic_states){
    		running_boltzmann += std::exp (-1/temperature*energy)/w2; 
    		if ( state_index <= running_boltzmann ) {
    			break;
    		}
    		count += 1; 
    	}

    	for (int i{0}; i<3; ++i){
    		(*LATTICE)[solvents_to_flip_n[i]]->orientation = orientation_configs[count][i];
    	}
    	*contacts = store_contacts  [count]; 
    	std::cout << "contacts on rejection is "; print(*contacts);
    	*sysEnergy = energetic_states[count];
    	std::cout << "New energy is " << *sysEnergy << std::endl;

    }

    else {
    	std::cout << "We rejected." << std::endl;
    	// return back to the original state 
    	(*Polymers)[index].chain[dop-1]->coords = loc_0;
    	(*LATTICE)[ lattice_index(loc_0, y, z) ]->coords = idx_v[r];

    	// do the switch 
    	(*LATTICE)[ lattice_index( idx_v[r], y, z) ] = (*LATTICE)[ lattice_index (loc_0, y, z) ];
    	(*LATTICE)[ lattice_index( loc_0, y, z )   ] = (*Polymers)[index].chain[dop-1]; 

    	// now that the solvent and polymer are in their original spatial state, i have to now get the solvents back into
    	// their original orientational state 
    	for (int j {0}; j<3; ++j) {
    		(*LATTICE)[solvents_to_flip[j]]->orientation = solvent_init_config [j];
    	}
    	std::cout << "contacts on rejection is "; print(c_contacts);
    	*contacts = c_contacts;
    }
	
    std::cout << "Orientation of monomer 0 is " << (*Polymers)[index].chain[dop-1]->orientation << std::endl;

	return; 
}


//============================================================
//============================================================
// 
// NAME OF FUNCTION: EndRotation
//
// PARAMETERS: index of polymer to perform EndRotation on, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a rotation of the monomer at the terminal index (the tail or head)
// of the polymer. As it stands, this function only rotates one molecule at the tail or head. The choice of head 
// or tail rotation comes from a random distribution. 
//
// PLANNED EXTENSION: Multiple molecules at the termini need to be rotated.    
//
// DEPENDENCIES: ZeroIndexRotation, FinalIndexRotation
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void EndRotation (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, double* sysEnergy, double temperature, \
	bool* IMP_BOOL, std::array<double,8>* E, std::array<double,8>* contacts){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        std::cout << "Zero index rotation!" << std::endl;
        TailRotation (Polymers, LATTICE, index, x, y, z, sysEnergy, temperature, IMP_BOOL, E, contacts); 
        return; 
    }
    else {
        std::cout << "Final index rotation!" << std::endl;
        HeadRotation (Polymers, LATTICE, index, x, y, z, sysEnergy, temperature, IMP_BOOL, E, contacts); 
        return; 
    }
   
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of EndRotation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
// 
// NAME OF FUNCTION: KinkJump
//
// PARAMETERS: index of which polymer to perform KinkJump on, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a kink jump when it finds a kink. 
// The kinks are shuffled before any is chosen. 
//
// DEPENDENCIES: findKinks, ChainToConnectivityMap, add_vectors, subtract_vectors
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void KinkJump (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>> >* memory){

    std::vector <int> k_idx = (*Polymers)[index].findKinks(); 

    (*rweight) = 0.0; 

    if (k_idx.size() == 0 ){
        *IMP_BOOL = false;
        // std::cout << "No kinks found in polymer..." << std::endl;
        return;
    }

    std::vector <int> idx_v;
    std::vector <std::array <int,3>> pos_v;  

    for (int idx: k_idx){

        std::array <int,3> d2       = subtract_arrays( &( (*Polymers) [index].chain[idx+2]->coords), &( (*Polymers) [index].chain[idx+1]->coords) ); 

        std::array <int,3> to_check = add_arrays     ( &( (*Polymers) [index].chain[idx]->coords ), &d2); 
        
        impose_pbc(&to_check, x, y, z); 

        // find locations where kink can jump to... 

        // if ( ! ( MonomerReporter (Polymers, &to_check) ) ){
        if ( ! ( MonomerReporter ( LATTICE, &to_check, y, z) ) ) {
        	// std::cout << "kink is at "; print(NewPol[index].chain[idx+1].coords);
        	// std::cout << "a free location is "; print(to_check); 
        	(*rweight) += 1;
        	idx_v.push_back(idx);  
        	pos_v.push_back(to_check); 
        }

    }
    
    if ( (*rweight) == 0){
    	// if nothing can be done, just get out. 
    	*IMP_BOOL = false; 
    	return;
    }



    // otherwise, let's fuckin go
	
	int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 
	
	// store the location in loc_0
	std::array <int,3> loc_0 = (*Polymers)[index].chain[idx_v[r]+1]->coords; 

	// make the change to the coordinates in *Polymers
	(*Polymers)[index].chain[idx_v[r]+1]->coords = pos_v[r];

	// make the change on the lattice
	// make the change only to the SOLVENT site
	(*LATTICE)[ lattice_index (pos_v[r], y, z)]->coords = loc_0;

	// do the switch
	// take the pointer of the solvent, and put it where the monomer was on the lattice
	(*LATTICE)[ lattice_index (loc_0, y, z)]    = (*LATTICE)[ lattice_index (pos_v[r], y, z)]; 
	(*LATTICE)[ lattice_index (pos_v[r], y, z)] = (*Polymers)[index].chain[ idx_v[r] + 1 ];

	// update memory 
	(*memory).first.push_back ( loc_0 );		// initial location of monomer 
	(*memory).second.push_back ( pos_v[r] );	// final location of monomer

	// update connmap
	// (*Polymers)[index].ChainToConnectivityMap(); 
	(*rweight) = (*rweight)/( k_idx.size() ); 

    return; 

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of KinkJump. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// NAME OF FUNCTION: BondVibration
//
// PARAMETERS: 
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a kink jump when it finds a kink. 
// The kinks are shuffled before any is chosen. 
//
// DEPENDENCIES: findKinks, ChainToConnectivityMap, add_vectors, subtract_vectors
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void BondVibration ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>> >* memory ){

	// choose a site to vibrate 
	int m_idx = rng_uniform ( 1, static_cast<int>((*Polymers)[index].chain.size())-2 );
	*rweight = 0; 

	// find d1 
	std::array <int,3> d1 = subtract_arrays ( &((*Polymers)[index].chain[m_idx]->coords)  , &((*Polymers)[index].chain[m_idx-1]->coords) );

	// find d2 
	std::array <int,3> d2 = subtract_arrays ( &((*Polymers)[index].chain[m_idx+1]->coords), &((*Polymers)[index].chain[m_idx]->coords  ) );   

	modified_direction (&d1, x, y, z); 
	modified_direction (&d2, x, y, z); 

	// m_idx-1 is at p, m_idx is at p+d1, m_idx+1 is at p+d1+d2 
	// let p+d1 be displaced in direction dx -> p+d1+dx 
	// so difference between p and p+d1+dx is d1+dx, difference between p+d1+dx and p+d1+d2 is d2-dx
	// so for a valid position d1+dx and d2-dx needs to be in adrns 
	
	std::array <int,26> good_dx = {}; 
	int count = 0, d_idx = 0; 

	std::sort ( adrns.begin(), adrns.end() ); 

	for ( std::array<int,3>& dx: adrns ){
		// impose important condition on dx 
		if ( ( std::binary_search (adrns.begin(), adrns.end(), add_arrays (&d1, &dx) ) ) && ( std::binary_search ( adrns.begin(), adrns.end(), subtract_arrays (&d2, &dx) ) ) ){

			std::array <int,3> temp = add_arrays ( & ( (*Polymers)[index].chain[m_idx]->coords), &dx );
			impose_pbc ( &temp, x, y, z); 

			if ( ! (MonomerReporter (LATTICE, &(temp), y, z ) ) ) {
				good_dx [count] = d_idx;
				count += 1; 
			}
		}
		d_idx += 1; 
	}
	// found all the directions where monomer can be vibrated ...
	// if there are none, 
	if ( count == 0 ){
		*IMP_BOOL = false; 
		return; 
	}

	// if there are some, choose a direction 
	int d_rng = rng_uniform (0 ,count-1); 
	std::array <int,3> new_loc = add_arrays ( & ( (*Polymers)[index].chain[m_idx]->coords), & (adrns[good_dx[d_rng]] ) );
	impose_pbc ( &new_loc, x, y, z); 

	// make the change to the polymer 
	std::array <int,3> loc0 = (*Polymers)[index].chain[m_idx]->coords; 
	(*Polymers)[index].chain[m_idx]->coords = new_loc; 

	// make the change on the lattice
	// make the change only to the solvent site... 

	(*LATTICE) [ lattice_index (new_loc, y, z)]->coords = loc0; 
	(*LATTICE) [ lattice_index (loc0, y, z) ] = (*LATTICE)[ lattice_index (new_loc, y, z)];
	(*LATTICE) [ lattice_index (new_loc, y, z)] = (*Polymers)[index].chain[m_idx]; 
	
	// update memory 
	(*memory).first.push_back  ( loc0 ); 
	(*memory).second.push_back ( new_loc ); 

	(*rweight) = (*rweight)/26;

	return; 

}

//============================================================
//============================================================
// 
// NAME OF FUNCTION: CrankShaft
//
// PARAMETERS: index of a polymer to perform CrankShaft on, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a kink jump when it finds a kink. 
// The kinks are shuffled before any is chosen. 
//
// DEPENDENCIES: findKinks, ChainToConnectivityMap, add_vectors, subtract_vectors
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void CrankShaft (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE,\
 int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, \
 std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory){

    std::vector <int> c_idx = (*Polymers)[index].findCranks(); 

    if ( c_idx.size()==0 ){
        *IMP_BOOL = false;  
        // std::cout << "No cranks..." << std::endl;
        return ; 
    }
    

    // std::shuffle (std::begin (c_idx), std::end(c_idx), std::default_random_engine() ); 
    int idx = c_idx[ rng_uniform (0, static_cast<int> (c_idx.size() - 1 ) ) ]; 
	std::vector <std::array <int,3>> pos_v1;
	std::vector <std::array <int,3>> pos_v2;

	(*rweight) = 0; 

	std::array <int,3> HingeToKink  = subtract_arrays(&( (*Polymers) [index].chain[idx+2]->coords), &( (*Polymers) [index].chain[idx+3]->coords) ); 
    std::array <int,3> HingeToHinge = subtract_arrays(&( (*Polymers) [index].chain[idx+3]->coords ), &( (*Polymers) [index].chain[idx]->coords ) );
    std::array <std::array <int,3>,3> drns = HingeSwingDirections (& (HingeToHinge), &(HingeToKink), x, y, z); 

	// std::array <int,3> d1 = drns[choice];

	for (std::array <int,3>& d: drns ){

        std::array <int,3> to_check_1 = add_arrays ( &( (*Polymers) [index].chain[idx]->coords), &d ); 
        std::array <int,3> to_check_2 = add_arrays ( &( (*Polymers) [index].chain[idx+3]->coords), &d ); 
        impose_pbc(&to_check_1, x, y, z); 
        impose_pbc(&to_check_2, x, y, z);   

		
      	// check if site is unoccupied 
        if ( ! ( MonomerReporter (LATTICE, &to_check_1, &to_check_2, y, z) ) ) {
        	pos_v1.push_back(to_check_1); 
        	pos_v2.push_back(to_check_2); 
        	(*rweight) += 1; 
        }
    }

    // if the site is unoccupied for sure 

    if ( (*rweight) == 0 ){
    	*IMP_BOOL = false; 
    	// std::cout << "Uncrankable." << std::endl;
    	return;
    }
  	

	int r = rng_uniform( 0, pos_v1.size() - 1 ); 

	std::array <int,3> loc_1 = (*Polymers)[index].chain[idx+1]->coords;
	std::array <int,3> loc_2 = (*Polymers)[index].chain[idx+2]->coords;

	// make the changes to the coordinates in *Polymers
	(*Polymers)[index].chain[idx+1]->coords = pos_v1[r]; 
	(*Polymers)[index].chain[idx+2]->coords = pos_v2[r];

	// make the change on the lattice 
	// make the change only to the SOLVENT site 

	(*LATTICE)[ lattice_index (pos_v1[r], y, z)]->coords = loc_1;
	(*LATTICE)[ lattice_index (pos_v2[r], y, z)]->coords = loc_2; 

	// do the switch 
	// take the pointer of the solvent, and put it where the monomer was on the lattice 
	(*LATTICE)[ lattice_index (loc_1, y, z)]     = (*LATTICE)[ lattice_index (pos_v1[r], y, z)]; 
	(*LATTICE)[ lattice_index (pos_v1[r], y, z)] = (*Polymers)[index].chain[idx+1];

	(*LATTICE)[ lattice_index (loc_2, y, z)]     = (*LATTICE)[ lattice_index (pos_v2[r], y, z)]; 
	(*LATTICE)[ lattice_index (pos_v2[r], y, z)] = (*Polymers)[index].chain[idx+2]; 

	// update memory 
	(*memory).first.push_back ( loc_1 );
	(*memory).first.push_back ( loc_2 );
	(*memory).second.push_back ( pos_v1[r] );
	(*memory).second.push_back ( pos_v2[r] );

	// update connmap 
	// (*Polymers) [index].ChainToConnectivityMap();
	(*rweight) = (*rweight)/3.0 ; 
  		
    return ;

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of CrankShaft. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
//
// !!!!!!! 
// NOTE FOR REPTATION CODE: The code is about all updating only the positions 
// of each particle in PolymersInGrid. One of the termini will move to a spot occupied by solvent. 
// If the monomer at index 0 is moving to a new spot, monomer at index 1 will take on coordinates of monomer 0 in the initial config, 
// monomer at index 2 will take coordinates of monomer at index 1 in the initial config, and so on. 
// Reverse the process if monomer at final index is moving to a new spot. 
// !!!!!!!
// 
// NAME OF FUNCTION: ForwardReptation 
//
// PARAMETERS: index of a polymer to reptate forward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: index of polymer on which to perform ForwardReptation, given a Grid, it will perform a forward reptation. 
// which means that the final index (size-1) is moving somewhere in its vicinity. 
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 
////////////////////////////////////////////////////////////

void ForwardReptation (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory){

    int deg_poly = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc0 = (*Polymers)[index].chain[0]->coords; 
    std::array <int,3> locf = (*Polymers)[index].chain[deg_poly-1]->coords;
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( locf, x, y, z ); 
    std::vector <std::array<int,3>> idx_v; 
    idx_v.reserve(26); 

    (*rweight) = 0; 

    for (std::array <int,3>& to_check: ne_list){

    	if ( to_check == (*Polymers) [index].chain[deg_poly-2]->coords ){
    		continue; 
    	}
        else if ( to_check == (*Polymers) [index].chain[0]->coords ) {
            (*rweight) += 1; 
            idx_v.push_back(to_check); 
        }

        else if ( ! (MonomerReporter (LATTICE, &to_check, y, z ) ) ){
			(*rweight) += 1; 
			idx_v.push_back(to_check);
		}
	}

	if ( (*rweight) == 0){
		*IMP_BOOL = false;
		return;
	}

	int r = rng_uniform( 0, idx_v.size()-1 );
	// if everything checks out, do the deed - make it slither forward 
	

	for (int i{0}; i<deg_poly; ++i){

		if ( (*LATTICE)[ lattice_index (idx_v[r], y, z)]->ptype[0] == 's' ){
			if ( i != deg_poly-1 ){
				
				(*memory).first.push_back  ( (*Polymers)[index].chain[i]->coords );
				(*Polymers)[index].chain[i]->coords = (*Polymers)[index].chain[i+1]->coords; 
				(*LATTICE)[ lattice_index( (*Polymers)[index].chain[i]->coords, y, z) ] = (*Polymers)[index].chain[i];
				(*memory).second.push_back ( (*Polymers)[index].chain[i]->coords );

			}
			else {

				(*memory).first.push_back  ( (*Polymers)[index].chain[i]->coords );

				// do the solvent switch 
				(*LATTICE)[ lattice_index (loc0, y, z) ] = (*LATTICE) [ lattice_index (idx_v[r], y, z) ]; 
				(*LATTICE)[ lattice_index (loc0, y, z) ]->coords = loc0;

				// update polymer 
				(*Polymers)[index].chain[i]->coords = idx_v[r]; 
				(*LATTICE)[ lattice_index ( (*Polymers)[index].chain[i]->coords, y, z) ] = (*Polymers)[index].chain[i];	
				(*memory).second.push_back ( (*Polymers)[index].chain[i]->coords );			

			}
		}
		else {

			if ( i != deg_poly-1 ){

				(*memory).first.push_back ( (*Polymers)[index].chain[i]->coords); 
				(*Polymers)[index].chain[i]->coords = (*Polymers)[index].chain[i+1]->coords; 
				(*LATTICE)[ lattice_index ( (*Polymers)[index].chain[i]->coords, y, z ) ] = (*Polymers)[index].chain[i];
				(*memory).second.push_back ( (*Polymers)[index].chain[i]->coords);

			}
			else {

				(*memory).first.push_back ( (*Polymers)[index].chain[i]->coords); 
				(*Polymers)[index].chain[i]->coords = (*memory).first[0];
				(*LATTICE)[ lattice_index ( (*Polymers)[index].chain[i]->coords, y, z ) ] = (*Polymers)[index].chain[i] ;
				(*memory).second.push_back ( (*Polymers)[index].chain[i]->coords); 

			}

		}
	}



	(*rweight) = (*rweight)/26; 

    return ; 

} 






//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ForwardReptation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
//
// !!!!!!! 
// NOTE FOR REPTATION CODE: The code is about all updating only the positions 
// of each particle in PolymersInGrid. One of the termini will move to a spot occupied by solvent. 
// If the monomer at index 0 is moving to a new spot, monomer at index 1 will take on coordinates of monomer 0 in the initial config, 
// monomer at index 2 will take coordinates of monomer at index 1 in the initial config, and so on. 
// Reverse the process if monomer at final index is moving to a new spot. 
// !!!!!!!
// 
// NAME OF FUNCTION: BackwardReptation 
//
// PARAMETERS: index of a polymer to reptate forward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: index of polymer on which to perform ForwardReptation, given a Grid, it will perform a forward reptation. 
// which means that the final index (size-1) is moving somewhere in its vicinity. 
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 
/////////////////////////////////////////////////////////////


void BackwardReptation (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory){

    int deg_poly = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc0 = (*Polymers)[index].chain[deg_poly-1]->coords; 
    std::array <int,3> locf = (*Polymers)[index].chain[0]->coords;
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list( locf, x, y, z ); 
    std::vector <std::array <int,3>> idx_v; 
    idx_v.reserve(26);

    (*rweight) = 0.0; 
     
    for (std::array <int,3>& to_check: ne_list){

    	if ( to_check == (*Polymers) [index].chain[1]->coords){
    		continue; 
    	}

        else if ( to_check == (*Polymers) [index].chain[deg_poly-1]->coords ) {
            (*rweight) += 1;
            idx_v.push_back(to_check); 
        }

        else if ( ! (MonomerReporter (LATTICE, &to_check, y, z) ) ){
    		(*rweight) += 1; 
    		idx_v.push_back(to_check); 
    	}
    }

    if ( (*rweight) == 0 ){
    	*IMP_BOOL = false; 
    	return;
    }

	int r = rng_uniform ( 0, idx_v.size()-1 );
	// std::cout<<"new location is: "; print(idx_v[r]); 
	for (int i{0}; i <deg_poly; ++i){

		if ((*LATTICE)[lattice_index (idx_v[r], y, z)]->ptype[0] == 's'){

			if ( i != deg_poly-1 ){

				(*memory).first.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords );
				(*Polymers) [index].chain[deg_poly-1-i]->coords = (*Polymers)[index].chain[deg_poly-2-i]->coords;
				(*LATTICE)[ lattice_index ((*Polymers)[index].chain[deg_poly-1-i]->coords, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i];	// transfer remaining information 
				(*memory).second.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords ); 

			}
			else {

				(*memory).first.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords );
				
				// do the solvent switch 
				(*LATTICE)[ lattice_index (loc0, y, z) ] = (*LATTICE) [ lattice_index (idx_v[r], y, z) ];
				(*LATTICE)[ lattice_index (loc0, y, z) ]-> coords = loc0;  

				// update polymer 
				(*Polymers) [index].chain[deg_poly-1-i]->coords = idx_v[r]; 
				(*LATTICE)[ lattice_index ( (*Polymers)[index].chain[deg_poly-1-i]->coords, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i];
				(*memory).second.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords );

			}
		}
		else {

			if ( i != deg_poly-1 ){

				(*memory).first.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords );
				(*Polymers)[index].chain[deg_poly-1-i]->coords = (*Polymers)[index].chain[deg_poly-2-i]->coords;
				(*LATTICE)[ lattice_index ( (*Polymers)[index].chain[deg_poly-1-i]->coords, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i];
				(*memory).second.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords);

			}
			else {

				(*memory).first.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords); 
				(*Polymers)[index].chain[deg_poly-1-i]->coords = (*memory).first[0]; 
				(*LATTICE)[ lattice_index ( (*Polymers)[index].chain[deg_poly-1-i]->coords, y, z) ] = (*Polymers)[index].chain[deg_poly-1-i];
				(*memory).second.push_back ( (*Polymers)[index].chain[deg_poly-1-i]->coords ); 

			}

		}
	}
	
	(*rweight) = (*rweight)/26; 

    return ; 

}




//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of BackwardReptation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
//
// !!!!!!! 
// NOTE FOR REPTATION CODE: The code is about all updating only the positions 
// of each particle in PolymersInGrid. One of the termini will move to a spot occupied by solvent. 
// If the monomer at index 0 is moving to a new spot, monomer at index 1 will take on coordinates of monomer 0 in the initial config, 
// monomer at index 2 will take coordinates of monomer at index 1 in the initial config, and so on. 
// Reverse the process if monomer at final index is moving to a new spot. 
// !!!!!!!
// 
// NAME OF FUNCTION: Reptation 
//
// PARAMETERS: index of a polymer to reptate forward or backward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: index of polymer on which to perform Reptation, given a Grid, it will perform a kink jump when it finds a kink. 
// The kinks are shuffled before any is chosen. 
//
// DEPENDENCIES: obtain_ne_list, ChainToConnectivityMap 
// as with every MC move, OccupancyMap, ConnectivityMaps need to be updated very very carefully. 
//
// THE CODE: 

void Reptation (std::vector<Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index, int x, int y, int z, bool* IMP_BOOL, double* rweight, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 

    if (num==0){
        // std::cout << "Backward reptation only!" << std::endl;
        BackwardReptation (Polymers, LATTICE, index, x, y, z, IMP_BOOL, rweight, memory); 
        return; 
    }
    else {
        // std::cout << "Forward reptation!" << std::endl;
        ForwardReptation (Polymers, LATTICE, index, x, y, z, IMP_BOOL, rweight, memory); 
        return; 
    }
    
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of Reptation. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


//============================================================
//============================================================
//
// 
// NAME OF FUNCTION: ChainRegrowth
//
// PARAMETERS: index of a polymer to reptate forward or backward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: It performs a chain regrowth on an existing polymer. It will uniformly choose to perform a tailspin or a headspin. 
// the second half of the code is to make sure the Solvent Vector is up to date after the polymer has been updated. 
//
// DEPENDENCIES: TailSpin, HeadSpin, extract_positions_tail, extract_positions_head 
//
// THE CODE: 
/////////////////////////////////////////////////////////////////////////////////////////////////////////


void ChainRegrowth (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int index_of_polymer, int x, int y, int z, bool* IMP_BOOL, double* rweight, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory, \
	int* index_monomer, int* back_or_front ){

	int deg_of_poly = (*Polymers)[index_of_polymer].deg_poly; 
	 

	// decide which end of the polymer do i want to move around 
    bool first_entry_bool = true; 

	*back_or_front = rng_uniform(0, 1); 

	std::vector <std::array<int,3>> old_cut;
	std::vector <std::array <int,3>> new_cut;	
	

	if ( *back_or_front == 0 ) {
	    // std::cout << "Performing tailspin.\n";	
	    *index_monomer = rng_uniform (1, static_cast<int>(std::floor(deg_of_poly/2)) ); // (1, deg_of_poly-2);
		old_cut.reserve (*index_monomer);
		new_cut.reserve (*index_monomer);
		

		old_cut = extract_positions_tail ( &((*Polymers)[index_of_polymer].chain), *index_monomer );

		TailSpin (Polymers, index_of_polymer, *index_monomer, x, y, z, IMP_BOOL, &first_entry_bool, rweight);  
		

		if ( !(*IMP_BOOL) ){

			// std::cout << "TailSpin was terminated early because there was no way to go forward. " << std::endl;
			for (int i{0}; i<*index_monomer; ++i){

				(*Polymers)[0].chain.at(i)->coords = old_cut.at(i); 

			}

		return; 
		
        }

		for ( int i{0}; i < *index_monomer; ++i ){
			new_cut.push_back ( (*Polymers)[index_of_polymer].chain[i]->coords );
		}
		// update memory 
		(*memory).first  = old_cut ;
		(*memory).second = new_cut ;
		
		std::vector <std::vector<std::array<int,3>>> master_linked_list; 
		std::vector <std::array<int,3>> link; 

		create_linked_list ( (*memory).first, (*memory).second, link, &master_linked_list, 1); 

		for ( std::vector <std::array<int,3>>& linked_list: master_linked_list ){
			
			if ( (*LATTICE)[ lattice_index( linked_list.back(), y, z ) ]->ptype[0] == 's' ){
				
				(*LATTICE)[ lattice_index (linked_list[0], y, z) ] = (*LATTICE)[ lattice_index (linked_list.back(), y, z) ];
				(*LATTICE)[ lattice_index (linked_list[0], y, z) ]->coords = linked_list[0]; 
			}

		}
		// std::cout << "printing out polymers pointer..." << std::endl;
		// print ((*Polymers)[0].chain[0]->coords);
		for ( int i{0}; i<*index_monomer; ++i) {
			// std::cout << "Reassigning LATTICE." << std::endl;			
			(*LATTICE)[ lattice_index ((*Polymers)[0].chain[i]->coords, y, z) ] = (*Polymers)[0].chain[i]; 
		}

	}

	else {

        *index_monomer = rng_uniform (static_cast<int>(std::ceil(deg_of_poly/2)), deg_of_poly-2 );
        // std::cout << "Performing headspin.\n"; 
		old_cut.reserve (deg_of_poly-(*index_monomer));
		new_cut.reserve (deg_of_poly-(*index_monomer));
		old_cut = extract_positions_head ( &((*Polymers)[index_of_polymer].chain), *index_monomer );
		
		HeadSpin (Polymers, index_of_polymer, *index_monomer, deg_of_poly, x, y, z, IMP_BOOL, &first_entry_bool, rweight); 

		if ( !(*IMP_BOOL) ){

			for (int i{*index_monomer+1}; i<deg_of_poly; ++i){
				(*Polymers)[0].chain.at(i)->coords = old_cut[i-(*index_monomer+1)]; 
			}

		return; 

		}


		for (int i{(*index_monomer)+1}; i<deg_of_poly; ++i ){
			new_cut.push_back ( (*Polymers)[index_of_polymer].chain[i]->coords );
		}

		// update memory 
		(*memory).first  = ( old_cut );
		(*memory).second = ( new_cut );

		std::vector <std::vector<std::array<int,3>>> master_linked_list; 
		std::vector <std::array<int,3>> link; 

		create_linked_list ( (*memory).first, (*memory).second, link, &master_linked_list, 1); 

		for ( std::vector <std::array<int,3>>& linked_list: master_linked_list ) {
			// std::cout << "printing out linked list... " << std::endl;
			// print (linked_list);
			
			if ( (*LATTICE)[ lattice_index(linked_list.back(), y, z) ]->ptype[0] == 's' ) {
				
				(*LATTICE)[ lattice_index (linked_list[0], y, z) ] = (*LATTICE)[ lattice_index (linked_list.back(), y, z) ];
				(*LATTICE)[ lattice_index (linked_list[0], y, z) ]->coords = linked_list[0]; 
				
			}
	
		}
		
		for ( int i{(*index_monomer+1)}; i<deg_of_poly; ++i) {
		
			(*LATTICE)[ lattice_index ((*Polymers)[0].chain[i]->coords, y, z) ] = (*Polymers)[0].chain[i]; 
		}
	}

	
	return; 
}




//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ChainRegrowth. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector <std::array <int,3>> extract_positions_tail (std::vector <Particle*>* chain, int pivot_idx){

	std::vector <std::array <int,3>> extracted;
	extracted.reserve(pivot_idx); 
	for (int i{0}; i<pivot_idx; ++i){
		extracted.push_back( (*chain)[i]->coords );
	}
	return extracted; 
}


std::vector <std::array <int,3>> extract_positions_head (std::vector <Particle*>* chain, int pivot_idx){

	int deg_of_poly = static_cast<int>( (*chain).size() ); 
	std::vector <std::array <int,3>> extracted; 
	extracted.reserve(deg_of_poly-pivot_idx-1); 
	for (int i{pivot_idx+1}; i<deg_of_poly; ++i){
		extracted.push_back ( (*chain)[i]->coords ); 
	}
	return extracted;
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of end extractions. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
//
// 
// NAME OF FUNCTION: TailSpin
//
// PARAMETERS: std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* b, bool* IMP_BOOL
// 
// WHAT THE FUNCTION DOES: It performs a tailspin. 
// the second half of the code is to make sure the Solvent Vector is up to date after the polymer has been updated. THIS IS A RECURSIVE FUNCTION. 
//
// DEPENDENCIES: impose_pbc, checkOccupancyTail
//
// THE CODE: 
//////////////////////////////////////////////////////////////

void TailSpin (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, \
                int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight){

	// std::cout << "index of monomer is " << index_of_monomer << std::endl;
    
    if (*first_entry_bool){
    	(*rweight) = 1; 
	    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	    std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
        *first_entry_bool = false; 
    }


	if (index_of_monomer == 0){
		// std::cout << "You have reached the final spot via tail spin!" << std::endl;
		*IMP_BOOL = true; 
		// (*Polymers)[index_of_polymer].ChainToConnectivityMap(); 
		return ; 
	}

    
    // std::cout << "Current pivot point is "; print(  (*PVec) [index_of_polymer].chain[index_of_monomer].coords );
    int rw_tmp = 0; 
    std::vector <std::array <int,3>> ind_v; 
	for (std::array <int,3>& d: adrns){ 
        
		std::array <int, 3> to_check = add_arrays( &( (*Polymers) [index_of_polymer].chain[index_of_monomer]->coords), &d);

		impose_pbc(&to_check, x, y, z); 
        // std::cout << "Growth location is "; print(to_check); 

		if (checkOccupancyTail(&to_check, Polymers, index_of_polymer, index_of_monomer)){ 
			continue; 
		}

		else {
			rw_tmp += 1;
			ind_v.push_back(to_check); 
		}
	}

	if (rw_tmp == 0){
		*IMP_BOOL = false;
	}

	else{
		(*Polymers)[index_of_polymer].chain[index_of_monomer-1]->coords = ind_v[rng_uniform(0, rw_tmp-1)]; 
		(*rweight) = (*rweight) * rw_tmp/26; 
		TailSpin (Polymers, index_of_polymer, index_of_monomer-1, x, y, z, IMP_BOOL, first_entry_bool, rweight); 
	}

	return; 

}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of TailSpin
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
//
// 
// NAME OF FUNCTION: HeadSpin
//
// PARAMETERS: index of a polymer to reptate forward or backward, a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: It performs a chain regrowth on an existing polymer. It will uniformly choose to perform a tailspin or a headspin. 
// the second half of the code is to make sure the Solvent Vector is up to date after the polymer has been updated. 
//
// DEPENDENCIES: impose_pbc, checkOccupancyTail
//
// THE CODE: 
//////////////////////////////////////////////////////////////

void HeadSpin (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, \
                int deg_poly,int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight){
    
     if (*first_entry_bool){
     	(*rweight) = 1; 
	    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	    std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
        *first_entry_bool = false; 
    }

	if (index_of_monomer == deg_poly-1){
		// std::cout << "You have reached the final spot of head spin!" << std::endl;
		*IMP_BOOL = true;
		// (*Polymers)[index_of_polymer].ChainToConnectivityMap(); 
		return ;
	}

	
	int rw_tmp = 0; 
	std::vector <std::array <int,3>> ind_v; 

	for (std::array <int,3>& d: adrns){
        // std::cout << "d is "; print(d); 
		std::array <int, 3> to_check = add_arrays ( &( (*Polymers) [index_of_polymer].chain[index_of_monomer]->coords ), &d);

		impose_pbc (&to_check, x, y, z); 
        // std::cout << "Growth location is "; print(to_check); 

		if (checkOccupancyHead(&to_check, Polymers, index_of_polymer, index_of_monomer)){
			continue; 
		}

		else {
			rw_tmp += 1; 
			ind_v.push_back(to_check); 
		}
	}	
	
	if (rw_tmp == 0){
		*IMP_BOOL = false; 
	}

	else {
		
		(*Polymers)[index_of_polymer].chain[index_of_monomer+1]->coords = ind_v[rng_uniform(0, rw_tmp-1)]; 
		(*rweight) = (*rweight) * rw_tmp/26; 
		HeadSpin (Polymers, index_of_polymer, index_of_monomer+1, deg_poly, x, y, z, IMP_BOOL, first_entry_bool, rweight);

		}

	return ;

}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of TailSpin
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


bool checkOccupancyTail(std::array <int,3>* loc, std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer){
	int Np = static_cast<int>( (*Polymers).size() );  

	for (int pnum = 0; pnum < Np; pnum++){

		if (pnum == index_of_polymer){

			int dpol = static_cast<int> ( (*Polymers)[pnum].chain.size() ); 
			for (int p = index_of_monomer; p < dpol; p++){

				if ( *loc == (*Polymers)[pnum].chain[p]->coords ){
					return true;
				}
				else {
					continue;
				}

			}

		}

		else {

			for (Particle*& p: (*Polymers)[pnum].chain){

				if ( *loc == p->coords ){
					return true;
				}

				else {
					continue;
				}

			}

		}

	}

	return false; 

}



bool checkOccupancyHead(std::array <int,3>* loc, std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer){
	int Np = static_cast<int>( (*Polymers).size() );  

	for (int pnum = 0; pnum < Np; pnum++){

		if (pnum == index_of_polymer){

			for (int p = 0; p < index_of_monomer+1; p++){

				if ( *loc == (*Polymers)[pnum].chain[p]->coords ){
					return true;
				}
				else {
					continue;
				}

			}

		}

		else {

			for (Particle*& p: (*Polymers)[pnum].chain){

				if ( *loc == p->coords ){
					return true;
				}

				else {
					continue;
				}

			}

		}

	}

	return false; 

}


//============================================================
//============================================================
//
// 
// NAME OF FUNCTION: OrientationFlip and PolymerFlip
//
// PARAMETERS: std::vector <Particle>* SolvVect
// 
// WHAT THE FUNCTION DOES: randomly picks a region of space, and perturbs the orientation of the solvent molecule 
// in that region 
// this function is a bit aggressive, imo. 
//
// DEPENDENCIES: impose_pbc
//
// THE CODE: 

void SolventFlip ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int x, int y, int z, double* rweight, int Nsurr, \
	std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory ){

	std::array <std::array<int,3>, 26> ne_list; 
	std::vector <int> solvent_indices; 
	// number of surrounding solvent molecules 
	for ( Polymer& pmer: (*Polymers) ){
		for ( Particle*& p: pmer.chain ){

			ne_list = obtain_ne_list ( p->coords, x, y, z );

			for ( std::array<int,3>& ne: ne_list ){
				if ( (*LATTICE).at(lattice_index (ne, y, z))->ptype[0] =='s' ){
					solvent_indices.push_back ( lattice_index(ne, y, z) ); 
				}
			}
		}
	}

	// get rid of duplicates 
	std::unordered_set<int> s;
	for ( int i: solvent_indices) {
		s.insert(i);
	}
	solvent_indices.assign ( solvent_indices.begin(), solvent_indices.end() ); 
	int Nmer  = static_cast<int>( (*Polymers)[0].chain.size() ); 

    int to_flip = rng_uniform(1, Nsurr); 
    
    *rweight = 1; 

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	
	std::vector <int> rvec (Nsurr);
	std::iota ( std::begin(rvec), std::end(rvec), 0);

	std::shuffle ( rvec.begin(), rvec.end(), std::default_random_engine(seed) );
	rvec = std::vector<int>(rvec.begin(), rvec.begin()+to_flip); 

	// print(rvec);

	int j = 0;
	for (int i: rvec ){

		(*memory).first.push_back ( { solvent_indices.at(i), (*LATTICE).at(solvent_indices.at(i))->orientation });
		LATTICE->at(solvent_indices.at(i))->orientation = rng_uniform (0, 5);  
		(*memory).second.push_back ( { solvent_indices.at(i), (*LATTICE).at(solvent_indices.at(i))->orientation });
		(*rweight) = (*rweight) * static_cast<double>(Nsurr-j)/static_cast<double>(Nmer+Nsurr-j); 
		j = j + 1; 
	}

	// for ( std::array<int,2>& a: (*memory).second){
	// std::cout << "solvent_idx = " << a[0] << ", orientation = " << a[1] << std::endl;
	// }
	
	return; 

}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void SolventFlipSingular ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int x, int y, int z, double* rweight, int Nsurr, \
	std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory ){
	
	// std::cout << "Size of Lattice is " << (*LATTICE).size() << std::endl;

	std::array <std::array<int,3>, 26> ne_list; 
	std::vector <int> solvent_indices; 
	// number of surrounding solvent molecules 
	for ( Polymer& pmer: (*Polymers) ){
		for ( Particle*& p: pmer.chain ){

			ne_list = obtain_ne_list ( p->coords, x, y, z );
			// std::cout << "Index of monomer is: ";
			// print (p->coords);

			for ( std::array<int,3>& ne: ne_list ){
				// std::cout << "neighbor is: "; print(ne);
				// std::cout << "lattice_index(ne , y, z) = " << lattice_index(ne, y, z) << std::endl;
				if ( (*LATTICE)[lattice_index (ne, y, z)]->ptype[0] =='s' ){
					solvent_indices.push_back ( lattice_index(ne, y, z) ); 
				}
			}
		}
	}


	// get rid of duplicates 
	std::unordered_set<int> s;
	for ( int i: solvent_indices) {
		s.insert(i);
	}
	solvent_indices.assign ( solvent_indices.begin(), solvent_indices.end() ); 

	// std::cout << "Got solvent indices." << std::endl;

	int Nmer  = static_cast<int>( (*Polymers)[0].chain.size() ); 
	
    *rweight = static_cast<double>(Nsurr)/static_cast<double>(Nmer+Nsurr); 

	int ridx = rng_uniform (0, Nsurr-1);
	// std::cout << "Nsurr is " << Nsurr << std::endl;
	(*memory).first.push_back( { solvent_indices.at(ridx), (*LATTICE)[solvent_indices.at(ridx)]->orientation } );
	(*LATTICE)[solvent_indices.at(ridx)]->orientation = rng_uniform (0, 25);  
	(*memory).second.push_back( { solvent_indices.at(ridx), (*LATTICE)[solvent_indices.at(ridx)]->orientation } );	

	// std::cout << "Exiting..." << std::endl;

	return; 
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



void SiteFlipSingular ( std::vector <Particle*>* LATTICE, \
	int x, int y, int z, int* nflips, \
	std::pair <std::vector <std::array<int,2>>, std::vector <std::array<int,2>>>* memory ){
	int flip_idx = 0; 

	(*nflips) = rng_uniform(0, static_cast<int>(x*y*z/20));

	for (int i{0}; i<*nflips; ++i){
		flip_idx = rng_uniform (0, x*y*z-1); 
		(*memory).first.push_back ( {flip_idx, (*LATTICE)[flip_idx]->orientation } );
		(*LATTICE)[flip_idx]->orientation = rng_uniform (0, 25);  
		(*memory).second.push_back ( {flip_idx, (*LATTICE)[flip_idx]->orientation } );
	}

	return; 

}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void PolymerFlip ( std::vector <Polymer>* Polymers, \
	double* rweight, int Nsurr, \
	std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory ){
    
	int Nmer  = static_cast<int>( (*Polymers)[0].chain.size() );     
    int to_flip = rng_uniform (1, Nmer );

    *rweight = 1; 

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	
	std::vector <int> rvec (Nmer);
	std::iota ( std::begin(rvec), std::end(rvec), 0); // this is a vector that goes from 0 to Nmer 
	std::shuffle ( rvec.begin(), rvec.end(), std::default_random_engine(seed) );
	rvec = std::vector<int>(rvec.begin(), rvec.begin()+to_flip); 

	int j = 0;
	for (int i: rvec) {
		(*memory).first.push_back( { i, (*Polymers)[0].chain.at(i)->orientation } );
		(*Polymers)[0].chain[i]->orientation = rng_uniform (0, 5); 
		*rweight = (*rweight) * static_cast<double>(Nmer-j)/static_cast<double>(Nmer+Nsurr-j);
		(*memory).second.push_back( { i, (*Polymers)[0].chain.at(i)->orientation } );
		j = j+1; 
		// std::cout << "index of monomer is " << i << ", initially orientation is " << (*memory).first[j-1][1] << std::endl;
		// std::cout << "index of monomer is " << i << ", finally orientation is "   << (*memory).second[j-1][1] << std::endl;
	}	
    
    return; 
}


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


void PolymerFlipSingular ( std::vector <Polymer>* Polymers, int nmonomer){
    
	std::random_device rg;
	std::mt19937 g(rd()); 

	std::vector <int> sample( static_cast<int>( (*Polymers)[0].chain.size() );  ); 
	std::iota (sample.begin(), sample.end(), 0);
	std::random_shuffle ( sample.begin(), sample.end(), g); 

	// this is the vector with the indices of monomers to flip 
	std::vector <int> chosen (sample.begin(), sample.begin()+nmonomer); 
	std::vector <int> initial_orientation; 
	std::vector <double> energetic_states;  
	for (int i{0}; i<chosen.size(); ++i){
		inital_orientation.push_back ((*Polymers)[0].chain[chosen[i]]->orientation);
	}	

	obtain_energetic_states ( Polymers, &energetic_states, &chosen, nmonomer ); 
	

	return; 
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void obtain_energetic_states ( std::vector <Polymer>* Polymers, std::vector<double>* energetic_states, \
	std::vector <int>* config, std:vector<std::vector<int>>* config_coll, \
	int nmonomer){

	for (int i {0}; i<5; ++i){
		(*config).push_back(rng_uniform(0,25)); 
		if ((*config).size() < nmonomer){
			obtain_energetic_states( Polymers, energetic_states, config, config_coll, nmonomer );
		}
		else {
			(*config_coll).push_back (*config); 
			(*config).pop_back(); 
		}
	}
	(*config).pop_back();
	return; 


}




///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


void SingleSolventExchange ( std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int x, int y, int z, bool* IMP_BOOL, double* rweight, \
	std::vector <std::array <int,2>>* s_memory ){

	// find indices of surrouding solvent molecules 
	std::array < std::array<int,3>, 26> ne_list; 
	std::vector <int> surr_solvent_indices; 
	int nmonomer = (*Polymers)[0].chain.size(); 

	for ( Polymer& pmer: (*Polymers) ){
		for ( Particle*& p: pmer.chain ){

			ne_list = obtain_ne_list ( p->coords, x, y, z); 
			for ( std::array<int,3>& ne: ne_list ){
				if ( (*LATTICE).at(lattice_index (ne, y, z))->ptype[0] == 's' ) {
					surr_solvent_indices.push_back ( lattice_index(ne, y, z) ); 
					continue; 
				}
			}
		}
	}

	// std::cout << "Is finding solvent indices was not an issue..." << std::endl;

	// get rid of repeated indices, obtain indices of solvent molecules surrounding polymer 
	std::sort ( surr_solvent_indices.begin(), surr_solvent_indices.end() ); 
	surr_solvent_indices.erase ( std::unique ( surr_solvent_indices.begin(), surr_solvent_indices.end() ), surr_solvent_indices.end() );

	// shuffle surr_solvent_indices 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	std::shuffle (std::begin(surr_solvent_indices), std::end(surr_solvent_indices), std::default_random_engine(seed));

  	int nsurr_solvent = static_cast<int>( surr_solvent_indices.size() );

  	// std::cout << "Number of surrounding solvents is " << nsurr_solvent << std::endl;

	int nexchange = 1 ; 
	int test_idx  = -1; 
	*rweight = 1; 
   
	for ( int j{0}; j < nexchange; ++j ){
		
		test_idx = rng_uniform (0, x*y*z-1); 
		// std::cout << "test_idx = " << test_idx << "." << std::endl;
		if ( (*LATTICE)[test_idx]->ptype != "m1" && test_idx != surr_solvent_indices[j] ){

			*IMP_BOOL = true; 

			// std::cout << "surr_solvent_indices[j] = " << surr_solvent_indices[j] << std::endl;
            // std::cout << "Location of solvent in surr = "; print ( (*LATTICE)[ surr_solvent_indices[j] ]->coords );
			// std::cout << "test_idx = " << test_idx << std::endl;
            // std::cout << "Location of far away solvent = "; print ( (*LATTICE)[ test_idx ]->coords ); 

			(*s_memory).push_back ( {surr_solvent_indices[j], test_idx } );

			// tmp     = (*LATTICE)[ test_idx ];  

			// std::cout << "Made the temporary pointer..." << std::endl;
            
            // if ( (*LATTICE)[test_idx]->ptype != (*LATTICE)[surr_solvent_indices[j]]->ptype ){
            //     std::cout << "IMPORTANT EXCHANGE taking place for indices right above!!!!" << std::endl;
            //     std::cout << "ptype of test idx is " << (*LATTICE)[test_idx]->ptype << std::endl;
            //     std::cout << "ptype of surr solvent is " << (*LATTICE)[surr_solvent_indices[j]]->ptype << std::endl; 
            // }

            std::string t_ptype = (*LATTICE)[test_idx]->ptype; 
            int orientation = (*LATTICE)[test_idx]->orientation; 

			(*LATTICE)[ test_idx ]->ptype        = (*LATTICE)[surr_solvent_indices[j] ]->ptype; 
			(*LATTICE)[ test_idx ]->orientation  = (*LATTICE)[surr_solvent_indices[j] ]->orientation; 
			
			// std::cout << "Transferred the surrounding solvent to another spot!" << std::endl;
            // std::cout << "--------------------------------------------------------------------" << std::endl;
			// change the test_idx
			(*LATTICE)[ surr_solvent_indices[j] ]->ptype            = t_ptype    ; 
			(*LATTICE)[ surr_solvent_indices[j] ]->orientation      = orientation;
			
			// std::cout << "(*LATTICE)[ssi[j]]->coords = "; print ( (*LATTICE)[surr_solvent_indices[j]]->coords );

			*rweight = (*rweight)*static_cast<double>(nsurr_solvent)/static_cast<double>(nsurr_solvent + nmonomer); 

		}
	}
	// std::cout << "ratio of solvent/solvent+monomer = " << static_cast<double>(nsurr_solvent)/static_cast<double>(nsurr_solvent + nmonomer) << std::endl; 
	// std::cout << "Total number of surr solv = " << nsurr_solvent << std::endl;
	// std::cout << "Reached end of Solvent Exchange." << std::endl;

	return; 
}







///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void PerturbSystem (std::vector <Polymer>* Polymers, std::vector <Particle*>* LATTICE, \
	int x, int y, int z, bool v, double* sysEnergy, double temperature, bool* IMP_BOOL, double* rweight, \
	std::array<double,8>* E, std::array<double,8>* contacts, int* nflips, \
	std::array <int,9>* attempts, int* move_number, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory3, \
	std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory2, \
	std::vector <std::array<int,2>>* s_memory, \
	int* monomer_index, int* back_or_front, int Nsurr ){

    int index = rng_uniform(0, static_cast<int>((*Polymers).size())-1); 
    int r = 1; // rng_uniform (1, 7);
 	// std::cout << x << y << z << v << r << index << *IMP_BOOL << rweight << (*attempts)[0] << move_number << std::endl;
 	// LATTICE->begin();

    switch (r) {
        case (1):
            if (v){
               printf("Performing end rotations...\n"); 
            }
            EndRotation		(Polymers, LATTICE, index, x, y, z, sysEnergy, temperature, IMP_BOOL, E, contacts);
            *move_number = 1;
            (*attempts)[0] += 1;
            break;  

        case (2):
            if (v){
               printf("Performing reptation...\n"); 
            }
            Reptation 		(Polymers, LATTICE, index, x, y, z, IMP_BOOL, rweight, memory3); 
            *move_number = 4; 
            (*attempts)[3] += 1;
            break; 
        
        case (3):
        	if (v) {
        		printf("Performing configuration sampling... \n"); 
        		// std::cout << "index of polymer is " << index << std::endl;
        	} 
        	ChainRegrowth 	(Polymers, LATTICE, index, x, y, z, IMP_BOOL, rweight, memory3, monomer_index, back_or_front ); 
            *move_number = 5; 
            (*attempts)[4] += 1;
        	break;
        
        
        case (4): 
        	if (v){
        		printf("Performing single solvent orientation flips... \n");
        	}
        	SolventFlipSingular (Polymers, LATTICE, x, y, z, rweight, Nsurr, memory2); 
            *move_number = 6; 
            (*attempts)[5] += 1;
        	break; 
        
        case (5):
            if (v) {
                printf("Performing single monomer orientation flip... \n");
            }
            PolymerFlipSingular (Polymers, rweight, Nsurr, memory2); 
            *move_number = 7; 
            (*attempts)[6] += 1;
            break;
        
        case (6): 
        	if (v) {
        		printf("Performing a random site flip... \n");
        	}
        	SiteFlipSingular ( LATTICE, x, y, z, nflips, memory2); 
        	*move_number = 8;
        	(*attempts)[7] += 1;
        	break; 

        case (7):
        	if (v) {
        		printf("Performing a solvent exchange... \n"); 
        	}
        	SingleSolventExchange (Polymers, LATTICE, x, y, z, IMP_BOOL,  rweight, s_memory); 
            *move_number = 9; 
            (*attempts)[8] += 1; 
            break;     
        
    }
    
    return;
}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of PerturbSystem
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

void ReversePerturbation (std::vector <Polymer>* Polymers, std::vector<Particle*>* LATTICE, int x, int y, int z, bool v, int move_number, int* nflips, \
	std::pair <std::vector<std::array<int,3>>, std::vector<std::array<int,3>>>* memory3, \
	std::pair <std::vector<std::array<int,2>>, std::vector<std::array<int,2>>>* memory2, \
	std::vector <std::array<int,2>>* s_memory, \
	int monomer_index, int back_or_front){

	switch (move_number){
		case (1):
			if (v) {
				printf("Reversing end rotations...\n");
			}

			// swap pointers 
			// take the location from memory.first and put it back on lattice (which is a solvent)  
			{
				Particle* tmp = (*LATTICE)[lattice_index((*memory3).second[0], y, z)]; 

				(*LATTICE)[ lattice_index((*memory3).second[0], y, z)] = (*LATTICE)[ lattice_index((*memory3).first[0], y, z)];
				(*LATTICE)[ lattice_index((*memory3).second[0], y, z)]->coords = (*memory3).second[0]; 

				(*LATTICE)[ lattice_index((*memory3).first[0], y, z)]  = tmp;  
				(*LATTICE)[ lattice_index((*memory3).first[0], y, z)]->coords  = (*memory3).first[0]; 
			}
			break;

		case (2):
			if (v) {
				printf("Reversing bond vibration...\n");
			}

			// swap pointers 
			// take the location from memory.first and put it back on lattice 
			{
				Particle* tmp = (*LATTICE)[lattice_index((*memory3).second[0], y, z)]; 

				(*LATTICE)[lattice_index((*memory3).second[0], y, z)] = (*LATTICE)[ lattice_index((*memory3).first[0], y, z)]; 
				(*LATTICE)[lattice_index((*memory3).second[0], y, z)]->coords = (*memory3).second[0]; 

				(*LATTICE)[ lattice_index((*memory3).first[0], y, z)] = tmp; 
				(*LATTICE)[ lattice_index((*memory3).first[0], y, z)]->coords = (*memory3).first[0]; 
				
			}
			break;

		case (3):
			if (v) {
				printf("Reversing crank shaft...\n");
			}
			// swap pointers 
			// take the location from memory.first and put it back on the lattice
			{
				Particle* tmp1 = (*LATTICE)[ lattice_index((*memory3).second[0], y, z)];
				Particle* tmp2 = (*LATTICE)[ lattice_index((*memory3).second[1], y, z)];

				(*LATTICE)[lattice_index((*memory3).second[0], y, z)] = (*LATTICE)[lattice_index((*memory3).first[0], y, z)]; 
				(*LATTICE)[lattice_index((*memory3).second[0], y, z)]->coords = (*memory3).second[0]; 
				(*LATTICE)[lattice_index((*memory3).second[1], y, z)] = (*LATTICE)[lattice_index((*memory3).first[1], y, z)]; 
				(*LATTICE)[lattice_index((*memory3).second[1], y, z)]->coords = (*memory3).second[1]; 

				(*LATTICE)[lattice_index((*memory3).first[0], y, z)]  = tmp1;
				(*LATTICE)[lattice_index((*memory3).first[0], y, z)]->coords = (*memory3).first[0]; 
				(*LATTICE)[lattice_index((*memory3).first[1], y, z)]  = tmp2; 
				(*LATTICE)[lattice_index((*memory3).first[1], y, z)]->coords = (*memory3).first[1];

			}
			break;

		case (4):
			if (v) {
				printf("Reversing reptation...\n");
			}
			// swap pointers 
			// take the location from memory.first and put it back on the lattice
			{	
				int deg_poly = static_cast<int>((*memory3).first.size());
				if ( (*memory3).first[0] != (*memory3).second[deg_poly-1]){

					// switch up the solvent molecule
					Particle* tmp = (*LATTICE)[ lattice_index((*memory3).second[deg_poly-1], y, z)]; // this is the monomer at the tail after perturbation 
					(*LATTICE)[ lattice_index((*memory3).second[deg_poly-1], y, z)] = (*LATTICE)[ lattice_index((*memory3).first[0], y, z) ];
					(*LATTICE)[ lattice_index((*memory3).second[deg_poly-1], y, z)]->coords = (*memory3).second[deg_poly-1]; 

					for (int i{0}; i<deg_poly; ++i){
						if ( i != deg_poly-1 ){
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ] = (*LATTICE)[ lattice_index((*memory3).second[i], y, z) ];
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ]->coords = (*memory3).first[i]; 
						}
						else {
							
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ] = tmp;
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ]->coords = (*memory3).first[i]; 
						}
					} 
				}

				else {
					Particle* tmp = (*LATTICE)[ lattice_index((*memory3).second[deg_poly-1], y, z)]; // this is the monomer at the tail after perturbation 
					
					for (int i{0}; i<deg_poly; ++i){
						if ( i != deg_poly-1 ){
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ] = (*LATTICE)[ lattice_index((*memory3).second[i], y, z) ];
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ]->coords = (*memory3).first[i];
						}
						else {
							
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ] = tmp; 
							(*LATTICE)[ lattice_index((*memory3).first[i], y, z) ]->coords = (*memory3).first[i]; 
						}
					}
				}
			}
			break;

		case (5):
			if (v) {
				printf ("Reversing chain regrowth...\n");
			}
			{
				// create a linked list 
				// this will need to be done recursively

				if (back_or_front == 0){
					// change the locations of the pointers on the lattice... 
					int l = static_cast<int>((*memory3).first.size()); 
					for (int i{0}; i<l; ++i){
						(*Polymers)[0].chain[i]->coords = (*memory3).first[i]; 
					}
					
					// print ((*Polymers)[0].chain[0]->coords);
					std::vector <std::vector<std::array<int,3>>> master_linked_list; 
					std::vector <std::array<int,3>> link;
					// std::cout << "Is it crapping out here?" << std::endl;
					create_linked_list ( (*memory3).second, (*memory3).first, link, &master_linked_list, 1); 

					/*
					for ( auto& V: master_linked_list){
						for (auto& v: V){
							for ( int& i: v ){
								std::cout << " | " << i; 
							}
							if ( v != V.back() ){
								std::cout << " -->";
							}
						}
						std::cout << std::endl;
					}
					*/
					// for each link in the list, understand if it is a circular list (where monomer jump between spots) 
					// or a straight list (where a monomer eventually jumps
					// to a spot where a solvent is)
					// std::cout << "master_linked_list.size() = " << master_linked_list.size() << std::endl;
					for ( std::vector <std::array<int,3>>& linked_list: master_linked_list){
						// std::cout << "printing out linked list... " << std::endl;
						// print (linked_list);

						if ( (*LATTICE)[ lattice_index( linked_list.back(), y, z ) ]->ptype[0] == 's' ){
							// std::cout << "linked_list.back() = "; print (linked_list.back() );
							// std::cout << "ptype is " << (*LATTICE)[ lattice_index( linked_list.back(), y, z ) ]->ptype;
							// std::cout << "this is a straight link." << std::endl;
							// transfer the solvent molecule to an appropriate location 
							// print ( (*LATTICE)[ lattice_index (linked_list[0], y, z) ]->coords );
							(*LATTICE)[ lattice_index (linked_list[0], y, z) ] = (*LATTICE)[ lattice_index (linked_list.back(), y, z) ];
							(*LATTICE)[ lattice_index (linked_list[0], y, z) ]->coords = linked_list[0]; 

						}
						else {
							// std::cout << "this is a circular link. do jack shit." << std::endl;
							// std::cout << "linked_list.back() = "; print (linked_list.back() );
							// std::cout << "ptype is " << (*LATTICE)[ lattice_index( linked_list.back(), y, z ) ]->ptype << std::endl;
							// do jack shit 
							// (*LATTICE)

						}

					}
					// std::cout << "printing out polymers pointer..." << std::endl;
					// print ((*Polymers)[0].chain[0]->coords);
					for ( int i{0}; i<monomer_index; ++i) {
						// std::cout << "Reassigning LATTICE." << std::endl;
						
						(*LATTICE)[ lattice_index ((*Polymers)[0].chain[i]->coords, y, z) ] = (*Polymers)[0].chain[i]; 

					}
				}
				
				else {

					int deg_poly = static_cast<int>((*Polymers)[0].chain.size()); 
					for (int i{monomer_index+1}; i<deg_poly; ++i){
						(*Polymers)[0].chain[i]->coords = (*memory3).first[i-(monomer_index+1)]; 
					}
					// print ( (*Polymers)[0].chain[monomer_index+1]->coords); 
					std::vector <std::vector<std::array<int,3>>> master_linked_list; 
					std::vector <std::array<int,3>> link;
					
					// std::cout << "Is it crapping out here?" << std::endl;

					create_linked_list ( (*memory3).second, (*memory3).first, link, &master_linked_list, 1); 
					/*
					for ( auto& V: master_linked_list){
						for (auto& v: V){
							for ( int& i: v ){
								std::cout << " | " << i; 
							}
							if ( v != V.back() ){
								std::cout << " -->";
							}
						}
						std::cout << std::endl;
					}
					*/

					// std::cout << "master_linked_list.size() = " << master_linked_list.size() << std::endl;

					for ( std::vector <std::array<int,3>>& linked_list: master_linked_list){
						// std::cout << "printing out linked list... " << std::endl;
						// print (linked_list);

						if ((*LATTICE)[ lattice_index (linked_list.back(), y, z ) ]->ptype[0] == 's' ){

							// std::cout << "linked_list.back() = "; print (linked_list.back() );
							// std::cout << "ptype is " << (*LATTICE)[ lattice_index( linked_list.back(), y, z ) ]->ptype;
							// std::cout << "this is a straight link." << std::endl;
							// transfer the solvent molecule to an appropriate location 
							// print ( (*LATTICE)[ lattice_index (linked_list[0], y, z) ]->coords );
							(*LATTICE)[ lattice_index (linked_list[0], y, z) ] = (*LATTICE)[ lattice_index (linked_list.back(), y, z) ];
							(*LATTICE)[ lattice_index (linked_list[0], y, z) ]->coords = linked_list[0];
							// std::cout << "type of linked_list[0] is " << (*LATTICE)[ lattice_index (linked_list[0], y, z) ]->ptype << std::endl;
							// print ((*LATTICE)[ lattice_index (linked_list[0], y, z) ]->coords); 
						}
						else {
							//std::cout << "this is a circular link. do jack shit." << std::endl;
							//std::cout << "linked_list.back() = "; print (linked_list.back() );
							// std::cout << "ptype is " << (*LATTICE)[ lattice_index( linked_list.back(), y, z ) ]->ptype << std::endl;
							// do jack shit 
							// (*LATTICE)

						}
					}
					// std::cout << "printing out polymers pointer..." << std::endl;
					// print ((*Polymers)[0].chain[monomer_index+1]->coords);
					for ( int i{monomer_index+1}; i<deg_poly; ++i) {
						// std::cout << "Reassigning LATTICE." << std::endl;
						(*LATTICE)[ lattice_index ((*Polymers)[0].chain[i]->coords, y, z) ] = (*Polymers)[0].chain[i]; 
					}
				}
			}
			break;

		case (6):
			if (v) {
				printf("Reversing a singular solvent flip...");
			}
			(*LATTICE)[ (*memory2).first[0][0] ]->orientation = (*memory2).first[0][1]; 
			break; 

		case (7):
			if (v) {
				printf("Reversing a singular monomer flip...");
			}
			(*LATTICE)[ lattice_index((*Polymers)[0].chain[ (*memory2).first[0][0]]->coords, y, z) ]->orientation = (*memory2).first[0][1]; 
			break;

		case (8):
			if (v) {
				printf("Reversing singular random flip...");
			}
			for (int i{0}; i<*nflips; ++i){
				(*LATTICE)[ (*memory2).first[i][0] ]->orientation = (*memory2).first[i][1]; 
			}
			break; 

		case (9):
			if (v) {
				printf("Reversing a solvent exchange...");
			}
			{
				for ( std::array<int,2>& a: (*s_memory) ) {
				    // std::cout << "First index is "  << a[0] << std::endl; 
				    // std::cout << "Second index is " << a[1] << std::endl;

				    Particle* tmp          = (*LATTICE)[ a[0] ]; 

				    (*LATTICE)[ a[0] ]         = (*LATTICE)[ a[1] ]; 
				    (*LATTICE)[ a[0] ]->coords = location ( a[0], x, y, z); 
				
				    (*LATTICE)[ a[1] ]         = tmp; 
				    (*LATTICE)[ a[1] ]->coords = location ( a[1], x, y, z);    
                }
			}
			break; 

	}

	// (*Polymers)[0].ChainToConnectivityMap (); 

	return; 
}

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
	std::regex characters ("[a-z]");

	int orientation = -1; 
	std::string ptype = "x"; 
	int index = -1; 
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
					// std::cout << "match for index is " << *a << std::endl;
					index = std::stoi ( *a );
					// std::cout << "match for index is " << *a << std::endl;
				}
			}
 

			std::regex_search ( s, match, characters );
			ptype = match[0].str()[0]; 

			// std::cout << "ptype is " << ptype << std::endl; 

			Particle* p_ptr = new Particle (location(index, x, y, z), ptype, orientation); 

			LATTICE.insert ( LATTICE.begin() + index, p_ptr); 

		}

	}

	std::cout << "Created lattice from file!" << std::endl;
	return LATTICE; 

}


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
	            // break;
	            continue; 
	            
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

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractPolymersFromTraj. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// NAME OF FUNCTION: ExtractIndexOfFinalMove
//
// PARAMETERS: std::string filename 
// 
// WHAT THE FUNCTION DOES: it looks at the trajectory file, and extracts the index of the final move. 
// 
// DEPENDENCIES: ExtractContentFromFile, makePolymer, checkValidityOfCoords 
//
// THE CODE: 

int ExtractIndexOfFinalMove(std::string trajectory){

    std::vector <std::string> contents = ExtractContentFromFile(trajectory); 

    int step_number{0}; 
    std::regex stepnum ("Dumping coordinates at step"); 

    for (std::string& s: contents){

        std::stringstream ss(s); 
        std::string temp; 
        int found; 

        if (std::regex_search (s, stepnum) ){
            while (!ss.eof() ){
                ss >> temp; 
                if (std::stringstream(temp) >> found){
                    step_number = found; 
                }
            }
        }
    }

    return step_number; 

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractIndexOfFinalMove 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
// 
// NAME OF FUNCTION: ExtractEnergyOfFinalMove
//
// PARAMETERS: std::string energy_file 
// 
// WHAT THE FUNCTION DOES: it looks at the trajectory file, and extracts the index of the final move. 
// 
// DEPENDENCIES: ExtractContentFromFile, makePolymer, checkValidityOfCoords 
//
// THE CODE: 

double ExtractEnergyOfFinalMove(std::string energy_file){

    std::vector <std::string> contents = ExtractContentFromFile(energy_file); 

    std::array <double,2> energy; 
    // int step_number{0}; 
    std::regex stepnum ("Dumping coordinates at step"); 

    for (std::string& s: contents){

        std::stringstream ss(s); 
        std::string temp; 
        double found; 

        int i{0};
        while (!ss.eof() ){
            ss >> temp; 
            if (std::stringstream(temp) >> found){
                energy[i] = found;
                ++i; 
            }
        }
    }

    return energy[0]; 

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ExtractIndexOfFinalMove 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//============================================================
//============================================================
// 
// NAME OF FUNCTION: CreateSolventVector
//
// PARAMETERS: int x, int y, int z, std::vector <Polymer>* PolymerVector  
// 
// WHAT THE FUNCTION DOES: it looks at the polymer vector and solvates the lattice with solvent molecules 
// 
// DEPENDENCIES: create_lattice_points 
//
// THE CODE: 

void AddSolvent (int x, int y, int z, std::vector <Particle*>* LATTICE){

	// std::vector <std::array<int,3>> lattice_points = create_lattice_pts (x, y, z); 
	// int count = 0;
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


void AddCosolvent (int x, int y, int z, double frac, int Nmonomer, std::vector <Particle*>* LATTICE){

	std::vector <int> used_idx; 

	int nsol2         = std::floor((x*y*z-Nmonomer)*frac); 
	int lattice_index = -1; 

	// generate lattice index 
	int count = 0;

	while ( count<nsol2 ){

		lattice_index = rng_uniform(0, x*y*z-1);
		if ( (*LATTICE)[lattice_index]->ptype[0] == 'm' ){
			used_idx.push_back (lattice_index); 
		}
		else if ( std::find (used_idx.begin(), used_idx.end(), lattice_index) != used_idx.end() ){
			;
		}
		else {
			(*LATTICE)[lattice_index]->ptype = "s2"; 
		}
		count += 1; 

	} 

	return;

}


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of CreateSolventVector
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



//============================================================
//============================================================
// 
// NAME OF FUNCTION: Translation
//
// PARAMETERS: int x, int y, int z, std::vector <Polymer>* PolymerVector  
// 
// WHAT THE FUNCTION DOES: It takes a polymer and translates it in a certain direction 
// 
// DEPENDENCIES: impose_pbc, add_arrays 
//
// THE CODE: 

/*
std::vector <Polymer> Translation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){

	// std::cout << "Index of polymer to be translated is " << index <<"." << std::endl;

	std::vector <Polymer> NewPol = *PolymerVector; 

	int deg_poly = static_cast<int> ((*PolymerVector)[index].chain.size()) ; 
	// bool breakout_completely = false; 
	size_t npoly = (*PolymerVector).size(); 

	// choose a random direction 
	std::array <int,3> rdirection = adrns[ rng_uniform(0,5) ]; 
	// std::cout <<"displacement direction is "; print(rdirection); 

	// choose a displacement length 
	int dl = rng_uniform(1, static_cast<int>(round(x/2)) ) ; 

	// std::cout << "displacement length is " << dl << std::endl;

	// define displacement vector 
	for (int i{0}; i<3; ++i){
		rdirection[i] = rdirection[i]*dl;
	}

	// std::cout << "displacement vector is "; print(rdirection); 

	std::array <int,3> to_check; 

	// displace each particle 
	for (Particle& p: NewPol[index].chain){

		// std::cout << "p.coords is (before assignment)"; print(p.coords); 
		// define a new position for the particle 
		to_check = add_arrays(&(p.coords), &rdirection); 


		// impose periodic boundary conditions 
		impose_pbc(&to_check, x, y, z); 
		// std::cout << "to_check is "; print(to_check); 

		///// start check loop 
		for (size_t i{0}; i < npoly; ++i){

			// do not consider the same polymer 
			if (static_cast<int>(i) == index){
				continue; 
			}

			for (const Particle& p_v: (*PolymerVector)[i].chain){

				if (to_check == p_v.coords){
					(*IMP_BOOL) = false; 
					// std::cout << "Translation process deemed impossible." << std::endl;
					return NewPol;  
				}
			}

		}
		///// end check loop

		p.coords = to_check; 

	}

	// now that every particle has been displaced correctly, make sure that every solvent particle that was originally in the spot of the monomer particle finds another home 
	// the edge case where when a monomer particle occupies a spot that previously occupied by another monomer particle... 

	// find out which particles went from a monomer-occupied region to a previously monomer occupied region 

	std::vector <std::array <int,3>> old_locations = extract_positions_tail(& ((*PolymerVector)[index].chain ), deg_poly );
	std::vector <std::array <int,3>> new_locations = extract_positions_tail(& (NewPol[index].chain), deg_poly );

	// sort both vectors
	std::sort ( old_locations.begin(), old_locations.end() ); 
	std::sort ( new_locations.begin(), new_locations.end() ); 

	// collect common elements 
	std::vector < std::array <int,3>> store; 
	std::set_intersection( old_locations.begin(), old_locations.end(), new_locations.begin(), new_locations.end(), std::back_inserter(store)); 

	// erase items in each vector that match the items in the set 

	for (const std::array <int,3>& a: store ){

		for (size_t j{0}; j < old_locations.size(); ++j){

			if (a == old_locations[j]){
				old_locations.erase(std::remove(old_locations.begin(), old_locations.end(), a), old_locations.end() ); 
				new_locations.erase(std::remove(new_locations.begin(), new_locations.end(), a), new_locations.end() ); 
				break; 
			}

		}
	}

	// now plant all the solvent particles that were originally there before the monomers came in to the spots originally occupied
	// by monomers which are not occupied by monomers anymore. 

	for (size_t i{0}; i<old_locations.size(); ++i){

		for (Particle& p: (*SolvVector)){

			if (p.coords == new_locations[i]){

				p.coords = old_locations[i];
				break; 

			}

		}

	}

	// std::cout << "Relevant polymer coordinates: " << std::endl;
	// NewPol[index].printChainCoords();


	// NewPol[index].ChainToConnectivityMap(); 
	*IMP_BOOL = true; 

	return NewPol;
}
*/


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of Translation
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

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




