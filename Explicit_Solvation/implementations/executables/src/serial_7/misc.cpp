#include "misc.h"
#include "lattice_directions.h"
/* 
=============================================================================================
These are some objects I have defined which I tend to use often. 
Helpful definitions which are employed often in the context of the z=6 lattice I am using.
=============================================================================================
*/ 

bool metropolis_acceptance(double E1, double E2, double kT){

	double dE = E2-E1; 
	double prob = std::exp(-1/kT*dE); 
	double r = rng_uniform(0.0, 1.0); 

	return (r<prob);
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


int heads_or_tails(int m_index, int deg_poly){

	int growth_dir {-1}; 

	if ( deg_poly % 2 == 0 ){
		growth_dir = (0.5 >= (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
	}
	else {
		if ( 0.5 == (m_index+1)/static_cast<double>(deg_poly+1) ){
			growth_dir = rng_uniform (0, 1);
		}
		else {
			growth_dir = (0.5 > (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
		}
	}
	return growth_dir;
}



void input_parser(int dfreq, int lfreq, int max_iter, bool r, \
	std::string positions, std::string topology, std::string dfile, \
	std::string efile, std::string mfile, std::string stats_file, \
	std::string lattice_file_read, std::string lattice_file_write, std::string ssfile){

	
	if (!r) {

		if (dfreq == -1 || lfreq == -1 || max_iter == -1) {
			std::cerr << "ERROR: No value for option f (frequency of dumping) (" << dfreq << ") and/or for option l (ell) (" << lfreq << ") and/or for option M (maximum number of moves to be performed) (" << max_iter << ") was provided. Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
	
		if ( positions== "__blank__" || topology == "__blank__" || dfile== "__blank__" || efile == "__blank__" || \
			mfile == "__blank__" || stats_file == "__blank__" || lattice_file_write == "__blank__" ){
			std::cerr << "polymer coords file is " << positions <<
			",\ntopology is " << topology <<
			",\npolymer coordinate dump file is " << dfile << 
			",\nenergy dump file is " << efile << 
			",\norientation file is " << mfile << 
			",\nmove statistics file is " << stats_file << 
			",\nlattice dump file is " << lattice_file_write << 
			",\nsolvation dump file is " << ssfile << "." << std::endl;

			std::cerr << "ERROR: No value for option p (polymer coordinate file) and/or\nfor option S (solvent coordinate file) and/or\n" <<
			"for option t (energy and geometry file) and/or\nfor option o (name of output dump file) and/or\nfor option e (name of orientation file) and/or\n" <<
			"for option s (name of move stats file) and/or\n for option u (name of energy dump file). Exiting..." << std::endl;
			exit (EXIT_FAILURE);    
		}
		
		// set up these files 
		std::ofstream polymer_dump_file (dfile);
		std::ofstream energy_dump_file (efile);
		std::ofstream orientation_dump_file (mfile); 
		std::ofstream statistics_dump_file (stats_file); 

		polymer_dump_file.close();
		energy_dump_file.close();
		orientation_dump_file.close();
		statistics_dump_file.close();

		if ( lattice_file_read != "__blank__" ){
			std::cerr << "Restart has not been requested. Do not provide a lattice file to read. Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		} 

		if ( ssfile != "__blank__"){
			std::ofstream solvation_dump_file(ssfile);
			solvation_dump_file.close();
		}

	}

	else {

		if ( dfreq == -1 || lfreq == -1 || max_iter == -1 ){
			std::cerr << "ERROR: No value for option f (frequency of dumping) and/or for option l (ell) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
	        exit (EXIT_FAILURE);	
		}

		if ( positions== "__blank__" || topology == "__blank__" || dfile== "__blank__" || efile == "__blank__" || \
			mfile == "__blank__" || stats_file == "__blank__" || lattice_file_read == "__blank__" ) {
			std::cerr << "polymer coords file is " << positions <<
			",\ntopology is " << topology <<
			",\npolymer coordinate dump file is " << dfile << 
			",\nenergy dump file is " << efile << 
			",\norientation file is " << mfile << 
			",\nmove statistics file is " << stats_file << 
			",\nlattice file to read is " << lattice_file_read << 
			",\nsolvation dump file is " << ssfile << "." << std::endl;
			std::cerr << 
			"ERROR: No value for option p (polymer coordinate file) and/or\n" << 
			"for option t (energy and geometry file) and/or\n" << 
			"for option o (name of output dump file) and/or\n" << 
			"for option e (name of orientation file) and/or\n" <<
			"for option s (name of move stats file) and/or\n" << 
			"for option u (name of energy dump file) was provided. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	return;
}


void filter_coords(const std::string& filename, int cutoff_step) {

	std::ifstream input_file(filename);

	if ( !(input_file.good()) ){
		// create the file;
		std::ofstream file(filename);
		file.close();
	}

	else {

		if (!input_file.is_open()) {
			std::cerr << "Error opening file for reading: " << filename << std::endl;
			std::cerr << "Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::string line; 
		std::vector<std::string> contents;
		while (std::getline(input_file, line)){
			contents.push_back(line);
		}
		input_file.close();

		std::ofstream output_file(filename);

		if (!output_file.is_open()) {
			std::cerr << "Error opening file for writing: " << filename << std::endl;
			std::cerr << "Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::regex pattern("^Dumping coordinates at step (\\d+)");
		std::smatch matches;
		int step_number {-1};

		for (std::string& single_line: contents){
			if (std::regex_search(single_line, matches, pattern)) {
				step_number = std::stoi(matches[1].str());
				if (step_number > cutoff_step){
					break;
				}
				else {
					output_file << single_line << std::endl;	
				}
			}
			else {
				output_file << single_line << std::endl;
			}
		}
		output_file.close();
	}
	return;
}

void filter_csv(const std::string& filename, int cutoff_step) {

	std::ifstream input_file(filename);

	if (!(input_file.good())){
		// create the file;
		std::ofstream file(filename);
		file.close();
	}

	else {
		if (!input_file.is_open()) {
			std::cerr << "Error opening file for reading: " << filename << std::endl;
			std::cerr << "Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::string line; 
		std::vector<std::string> contents;
		while (std::getline(input_file, line)){
			contents.push_back(line);
		}
		input_file.close();

		std::ofstream output_file(filename);

		if (!output_file.is_open()) {
			std::cerr << "Error opening file for writing: " << filename << std::endl;
			std::cerr << "Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::regex pattern("(\\d+)(?![^|]*\\|)");
		std::smatch matches;
		int step_number {-1};

		for (std::string& single_line: contents){
			// std::cout << "line = " << single_line << std::endl;
			if (std::regex_search(single_line, matches, pattern)) {
				step_number = std::stoi(matches[0]);
				// std::cout << "step_number = " << step_number << std::endl;
				if (step_number > cutoff_step){
					break;
				}
				else {
					output_file << single_line << std::endl;
				}
			}
			else {
				output_file << single_line << std::endl;
			}
		}
		output_file.close();
	}
	return;
}

void filter_orientations(const std::string& filename, int cutoff_step) {

	std::ifstream input_file(filename);

	if (! (input_file.good())){
		std::ofstream file (filename);
		file.close();
	}

	else {
		if (!input_file.is_open()) {
			std::cerr << "Error opening file for reading: " << filename << std::endl;
			std::cerr << "Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::string line; 
		std::vector<std::string> contents;
		while (std::getline(input_file, line)){
			contents.push_back(line);
		}
		input_file.close();

		std::ofstream output_file(filename);

		if (!output_file.is_open()) {
			std::cerr << "Error opening file for writing: " << filename << std::endl;
			std::cerr << "Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}

		std::regex pattern("^START for Step (\\d+)");
		std::smatch matches;
		int step_number {-1};

		for (std::string& single_line: contents){
			if (std::regex_search(single_line, matches, pattern)) {
				step_number = std::stoi(matches[1].str());
				if (step_number > cutoff_step){
					break;
				}
				else {
					output_file << single_line << std::endl;
				}
			}
			else {
				output_file << single_line << std::endl;
			}
		}
		output_file.close();
	}
	return;

}

void resetting_containers(std::array<double,8>* cs1_i, std::array<double,8>*cs1_f, 
std::array<double,8>* cs2_i, std::array<double,8>*cs2_f,
std::array<double,8>* cpair_i, std::array<double,8>* cpair_f, 
double* Es1_i, double* Es1_f, double* Es2_i, double* Es2_f,
double* Epair_i, double* Epair_f){

	(*cs1_i).fill(0);
	(*cs1_f).fill(0);
	(*cs2_i).fill(0);
	(*cs2_f).fill(0);
	(*cpair_i).fill(0);
	(*cpair_f).fill(0);
	Es1_i = 0;
	Es1_f = 0;
	Es2_i = 0;
	Es2_f = 0;
	Epair_i = 0;
	Epair_f = 0;
	return;

}