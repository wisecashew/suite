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

/* ==================================================
These are some objects I have defined which I tend to use often. 
Helpful definitions which are employed often in the context of the z=6 lattice I am using.
====================================================*/ 

const std::vector <int> ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1}; 	// unit directions 
const std::vector <std::vector <int>> drns = {ex, nex, ey, ney, ez, nez};  							// vector of unit directions 

const std::array <int,3> ax{1,0,0}, nax{-1,0,0}, ay{0,1,0}, nay{0,-1,0}, az{0,0,1}, naz{0,0,-1};		// unit directions
const std::array <std::array <int,3> ,6> adrns = {ax, nax, ay, nay, az, naz}; 

//=====================================================
// impose periodic boundary conditions on vector 
//=====================================================

void impose_pbc(std::vector <int>* vect, int x_len, int y_len, int z_len){
	for (int i{0}; i<3; i++){
		if (i==0){
			(*vect).at(i) = ((((*vect).at(i))%x_len)+x_len)%x_len; 
		}
		else if (i==1){
			(*vect).at(i) = ((((*vect).at(i))%y_len)+y_len)%y_len; 
		}
		else {
			(*vect).at(i) = ((((*vect).at(i))%z_len)+z_len)%z_len; 
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

// modifies array to be a legit direction

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
	for (int i: v){
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
/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


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
std::vector <std::vector <int>> create_lattice_pts(int x_len, int y_len, int z_len){
	std::vector<std::vector <int>> lattice_pts;
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

std::array <std::array <int,3>, 6> obtain_ne_list(std::array <int,3> loc, int x_len, int y_len, int z_len){
	std::array <std::array <int,3>, 6> nl;
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



std::array <double,7> ExtractTopologyFromFile(std::string filename){
    
    std::array <double, 7> info_vec; 
    double info; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(filename); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), Emm_a ("Emm_a"), Emm_n ("Emm_n"), Ems ("Ems"), eof ("END OF FILE"); 
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

    	else if (std::regex_search (s, Ems)){
    		double info = NumberExtractor(s); 
    		info_vec[6] = info; 
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
    int size_ = locations.size(); 
    for (int i=0; i<static_cast<int>(size_); i++){
        unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
        std::mt19937 generator(seed); 
        std::uniform_int_distribution<int> distribution (0,5); 
        pmer_spins.push_back(distribution(generator));
    }

    std::vector <Particle> ptc_vec; 

    for (int i=0;i<static_cast<int>( size_ ); i++ ){
        Particle p (locations.at(i), type_m, pmer_spins.at(i)); 
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

void InputParser(bool a, bool r, int Nacc, int dfreq, int max_iter, 
	std::string positions, std::string topology, std::string dfile, 
	std::string efile, std::string restart_traj){

	// if acceptance criterion is NOT CALLED 

	if (!a) {
		if (Nacc != -1) {
			std::cerr << "ERROR: You do not need to provide a -N option if you are not looking for accepted configuration statistics. Use -a if you want to use -N. Safeguarding against uncontrolled behavior. Exiting..." << std::endl;
            exit (EXIT_FAILURE);
		}

		if (dfreq == -1 || max_iter == -1){
            std::cerr << "ERROR: No value for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
            exit (EXIT_FAILURE);
        }

	}

	// if restart is NOT CALLED
	if (!r){
		if (restart_traj != "blank"){
            std::cerr << "ERROR: You cannot ask for a trajectory coordinate file with -T without the -r flag. Safeguard against uncontrolled behavior. Exiting..." << std::endl;
            exit (EXIT_FAILURE); 
        }
        if (dfreq == -1 || max_iter == -1 ){
            std::cerr << "ERROR: No value for option N (number of accepted MC moves to have) and/or for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
            exit (EXIT_FAILURE);
        }

    }
    
    else {
        if (dfreq == -1 || max_iter == -1 ){
            std::cerr << "ERROR: No value for option N (number of accepted MC moves to have) and/or for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
            exit (EXIT_FAILURE);
        }
        else if (restart_traj == "blank"){
            std::cerr << "ERROR: Need a trajectory file -T if -r option is being called. Exiting..." << std::endl;
            exit(EXIT_FAILURE); 
        }

	}

	// simulation requires an initial "positions" file, a topology file, a coordinate dump, and an energydump file 

    if (positions=="blank" || topology == "blank" || dfile=="blank" || efile == "blank" ){
        std::cerr << "positions is " << positions <<", topology is " << topology <<", coordinate dump file is " << dfile << ", energy dump file is " << efile << std::endl;
        std::cerr << "ERROR: No value for option p (positions file) and/or for option t (energy and geometry file) and/or for option o (name of output dump file)" <<
        " and/or for option u (name of energy dump file) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);    
    }

    // set up outputs 

    if (r){
        std::ofstream dump_file(dfile, std::ios::app); 
        std::ofstream energy_dump_file (efile, std::ios::app); 
    }

    else if (!r){
        std::ofstream dump_file (dfile);
        std::ofstream energy_dump_file (efile);
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


bool checkForOverlaps(std::vector <Polymer> PolymerVector){
    
    std::vector <std::array <int,3>> loc_list; 

    for (Polymer& pmer: PolymerVector){
        for (Particle& p: pmer.chain){
            // check if element exists in vector 
                if (std::find(loc_list.begin(), loc_list.end(), p.coords) != loc_list.end() ){
                    std::cerr << "you have a repeated element." << std::endl;
                    // std::cout << "current element is: " << std::endl;
                    // print(p.coords); 
                    //print(loc_list);
                    return false; 
                    }
            
                else{
                    loc_list.push_back(p.coords);  
                }
            }
        }    
    
    std::cout << "Input file has no overlaps!" << std::endl;
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

bool checkConnectivity(std::vector <Polymer> PolymerVector, int x, int y, int z) {
    std::array <int,3> d1 = {0, 0, 1}; 
    std::array <int,3> dx = {0, 0, x-1}; 
    std::array <int,3> dy = {0, 0, y-1}; 
    std::array <int,3> dz = {0, 0, z-1}; 
    for (Polymer& pmer: PolymerVector){
        size_t length = pmer.chain.size(); 
        std::array <int,3> connection; 
        for (int i{1}; i<static_cast<int>(length); ++i){
            
            connection = subtract_arrays(&(pmer.chain[i].coords), &(pmer.chain[i-1].coords));
            impose_pbc(&connection, x, y, z);
            std::sort(connection.begin(), connection.end()); 
            if (connection == d1 || connection == dx || connection == dy || connection == dz){
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


//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of checkConnectivity. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#





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
// NAME OF FUNCTION: CalculateEnergy
//
// PARAMETERS: some attributes of the Grid Object
// WHAT THE FUNCTION DOES: Calculates energy of the current Grid. Critical to correctly evolve system. 
// includes nearest neighbor interactions with and without directional effects.  
//
// DEPENDENCIES: obtain_ne_list 
//
// OPTIMIZATION OPPORTUNITY: I am double counting monomer-monomer interactions. This can possibly be avoided. 
//
// THE CODE: 

double CalculateEnergy(std::vector <Polymer>* PolymerVector, int x, int y, int z, double Emm_a, double Emm_n, double Ems){
    double Energy {0.0}; 
    bool b {false}; 
    // polymer-polymer interaction energies 
    for (const Polymer& pmer: (*PolymerVector)){
        for (const Particle& p: pmer.chain){
            
            // std::vector <Particle> part_vec = pmer.ConnectivityMap[p];
            
            std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(p.coords, x, y, z); // get neighbor list 

            // consider bonded interactions  
            
            for (std::array <int, 3>& loc: ne_list){
                // std::cout << "neighbor list is: "; print (loc); 
                // get the particle at loc 
            	b = false; 
                for (const Polymer &pmer_: (*PolymerVector)){
                	for (const Particle& p_: pmer_.chain){

                		if (p_.coords == loc){

                			if (p_.orientation == p.orientation){
                				Energy += Emm_a *0.5; 
                				b = true; 
                				break;
                			}
                			else {
                				Energy += Emm_n * 0.5; 
                				b = true; 
                				break; 
                			}

                		}

                		else {
                			continue; 
                		}

                	}

                	if (b){
                		break; 
                	}
                }

                if (!b){
                	Energy += 0.5*Ems;
                	continue;  
                }
                else {
                	continue; 
                }

            }
        }
    }
    
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


void dumpPositionsOfPolymers (std::vector <Polymer>* PolymersInGrid, int step, std::string filename){
    // std::ostringstream os; 
    std::ofstream dump_file(filename, std::ios::app); 
    dump_file <<"Dumping coordinates at step " << step << ".\n";
    
    // dump_file <<"Dumping coordinates at step " << step << ".\n";
    int count = 0; 
    for (const Polymer& pmer: *PolymersInGrid){
        
        dump_file <<"Dumping coordinates of Polymer # " << count << ".\n";
        dump_file<<"START" << "\n";
        for (const Particle& p: pmer.chain){
            for (int i: p.coords){
                dump_file << i << " | "; 
            }
            dump_file << "\n"; 
        }
        ++count ; 
    }
    dump_file <<"END"<<"\n";
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


//============================================================
//============================================================
//
// NAME OF FUNCTION: dumpEnergyOfGrid 
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


void dumpEnergy (double sysEnergy, int step, std::string filename, bool first_call){
    std::ofstream dump_file(filename, std::ios::app); 
    // std::ostringstream os; 
    if (first_call){
        dump_file <<"This file contains energy of the Grid at certain points in the simulation.\n";
        dump_file <<"Energy | Step_Number\n" ;
    }
    
    dump_file << sysEnergy << " | " << step << "\n";
    

    return; 

}



//============================================================
//============================================================
// 
// NAME OF FUNCTION: TailRotation
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

std::vector <Polymer> TailRotation(std::vector <Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){


    // get the neighborlist of particle at index 1 
    std::array <int,3> loc_0 = (*PolymerVector)[index].chain[0].coords; 
    std::array <int,3> loc_1 = (*PolymerVector)[index].chain[1].coords;
    std::vector <Polymer> NewPol {*PolymerVector};
    
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(loc_1, x, y, z) ; 
    bool b {false}; 

	size_t tries = 0; 
    for (std::array <int,3>& to_rot: ne_list){

    	b = false; 

    	if (loc_0 == to_rot){
    		++tries;
    		continue;
    	}

    	for (const Polymer& pmer: (*PolymerVector)){
    		for (const Particle& p: pmer.chain){

    			if (p.coords == to_rot){
    				b = true;
    				break; 
    			}

    		}

    		if (b){
    			break;
    		}

    	}
  		if (!b){
  			// std::cout << "number of tries is " << tries << std::endl;
  			// std::cout << "max number of tries can be " << ne_list.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl;  
			NewPol[index].chain[0].coords = to_rot; 
			break;
  		}
  		else {
  			// std::cout << "Occupied, pick another place to rotate." << std::endl; 
			++tries; 
  		}

    }


	if (tries == ne_list.size()){
		*IMP_BOOL = false; 
	}
	
	NewPol[index].ChainToConnectivityMap(); 
	return NewPol; 
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

std::vector <Polymer> HeadRotation(std::vector <Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){


    // get the neighborlist of particle at index 1 
    int dop = (*PolymerVector)[index].deg_poly; 
    std::array <int,3> loc_0 = (*PolymerVector)[index].chain[dop-1].coords; 
    std::array <int,3> loc_1 = (*PolymerVector)[index].chain[dop-2].coords;
    std::vector <Polymer> NewPol {*PolymerVector};
    
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(loc_1, x, y, z) ; 
    bool b {false}; 

	size_t tries = 0; 
    for (std::array <int,3>& to_rot: ne_list){

    	b = false; 

    	if (loc_0 == to_rot){
    		++tries;
    		continue;
    	}

    	for (const Polymer& pmer: (*PolymerVector)){
    		for (const Particle& p: pmer.chain){

    			if (p.coords == to_rot){
    				b = true;
    				break; 
    			}

    		}

    		if (b){
    			break;
    		}

    	}
  		if (!b){
  			// std::cout << "number of tries is " << tries << std::endl;  
  			// std::cout << "max number of tries can be " << ne_list.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl;  
			NewPol[index].chain[dop-1].coords = to_rot; 
			break;
  		}
  		else {
  			// std::cout << "Occupied, pick another place to rotate." << std::endl; 
			++tries; 
  		}

    }


	if (tries == ne_list.size()){
		*IMP_BOOL = false; 
	}
	
	NewPol[index].ChainToConnectivityMap(); 
	return NewPol; 
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


std::vector <Polymer> EndRotation(std::vector <Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return TailRotation(PolymerVector, index, x, y, z, IMP_BOOL); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return HeadRotation(PolymerVector, index, x, y, z, IMP_BOOL); 

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


std::vector <Polymer> KinkJump(std::vector <Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){

    std::vector <Polymer> NewPol = *PolymerVector; 

    std::vector <int> k_idx = NewPol.at(index).findKinks(); 

    if (k_idx.size() == 0 ){
        *IMP_BOOL = false;
        // std::cout << "No kinks found in polymer..." << std::endl;
        return *PolymerVector;
    }

    // std::shuffle(std::begin(k_idx), std::end(k_idx), std::default_random_engine() ); 
	size_t tries = 0; 

	bool b {false};

    for (int idx: k_idx){

        // std::cout << "idx right before kink spot is " << idx << std::endl; 

        // std::array <int,3> d1 = subtract_arrays(&(NewG.PolymersInGrid[index].chain[idx+1].coords), &(NewG.PolymersInGrid[index].chain[idx].coords) );
        std::array <int,3> d2 = subtract_arrays(&(NewPol[index].chain[idx+2].coords), &(NewPol[index].chain[idx+1].coords) ); 

        std::array <int,3> to_check = add_arrays( &( NewPol[index].chain[idx].coords ), &d2); 
        
        impose_pbc(&to_check, x, y, z); 

        b = false; 

    	for (const Polymer& pmer: (*PolymerVector)){
    		for (const Particle& p: pmer.chain){

    			if (p.coords == to_check){
    				b = true;
    				break; 
    			}

    		}

    		if (b){
    			break;
    		}

    	}
    	// if the site is unoccupied for sure 
  		if (!b){			
  			// std::cout << "number of tries is " << tries << std::endl;  
  			// std::cout << "max number of tries for kjump can be " << k_idx.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl;  
			NewPol[index].chain[idx+1].coords = to_check; 
			break;
  		}
  		else {
  			// std::cout << "Occupied, pick another place to kink jump." << std::endl; 
			++tries; 
  		}

    }

	
	
	if (tries == k_idx.size()) {
        *IMP_BOOL = false; 
		//printf("No spot was available to kink jump! Position will be accepted by default!\n");
	}
	NewPol[index].ChainToConnectivityMap(); 

    return NewPol; 

}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of KinkJump. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



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


std::vector <Polymer> CrankShaft(std::vector <Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){


    // std::array <int,3 > ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1};     // unit directions 
    // std::array <std::array <int,3>,6> ultimate_drns = {ex, nex, ey, ney, ez, nez}; 

    std::vector <Polymer> NewPol {*PolymerVector}; 

    std::vector <int> c_idx = NewPol[index].findCranks(); 

    if ( c_idx.size()==0 ){
        *IMP_BOOL = false; 
        // std::cout << "No cranks in this polymer..." << std::endl; 
        return NewPol; 
    }
    
    bool b {false};  
    

    std::shuffle (std::begin (c_idx), std::end(c_idx), std::default_random_engine() ); 
	size_t tries = 0; 
    for (int idx: c_idx){

        
        
        std::array <int,3> HingeToKink = subtract_arrays(&(NewPol[index].chain[idx+2].coords), &(NewPol[index].chain[idx+3].coords) ); 
        std::array <int,3> HingeToHinge = subtract_arrays(&(NewPol[index].chain[idx+3].coords ), &(NewPol[index].chain[idx].coords ) );

        std::array <std::array <int,3>,3> drns = HingeSwingDirections (& (HingeToHinge), &(HingeToKink), x, y, z); 
        int choice = rng_uniform(0,2); 

        std::array <int,3> d1 = drns[choice];  //subtract_vectors( &(NewG.PolymersInGrid.at(index).chain.at(idx+3).coords), &(NewG.PolymersInGrid.at(index).chain.at(idx+2).coords) );

        std::array <int,3> to_check_1 = add_arrays ( &(NewPol[index].chain[idx].coords), &d1 ); 
        std::array <int,3> to_check_2 = add_arrays ( &(NewPol[index].chain[idx+3].coords), &d1 ); 
        impose_pbc(&to_check_1, x, y, z); 
        impose_pbc(&to_check_2, x, y, z);   

        b = false; 

        
      	for (const Polymer& pmer: (*PolymerVector)) {
      		for (const Particle& p: pmer.chain){

      			if (p.coords == to_check_1 || p.coords == to_check_2){
      				b = true; 
      				break;
      			}

      		}

      		if (b) {
      			break; 
      		}

      	}
      	// if the site is unoccupied for sure 
      	if (!b){
      		// std::cout << "number of tries is " << tries << std::endl;
  			// std::cout << "max number of tries for cshaft can be " << c_idx.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl; 
      		NewPol[index].chain[idx+1].coords = to_check_1;
      		NewPol[index].chain[idx+2].coords = to_check_2;
      		break;
      	}

      	else {
      		// std::cout << "Occupied, pick another place to crank." << std::endl; 
			++tries; 
      	}



    }
    // std::cout << "Loop exit?" << std::endl;
	if (tries == c_idx.size()){
        *IMP_BOOL = false; 
		// printf("No position was available for a crankshaft! position will be accepted by default!\n");
	}

	NewPol[index].ChainToConnectivityMap();

    // std::cout << "Have you reached the end of the Crank?" << std::endl;
    return NewPol;

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

std::vector <Polymer> ForwardReptation(std::vector <Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){

    std::vector <Polymer> NewPol {*PolymerVector}; 
    int deg_poly = (*PolymerVector)[index].deg_poly; 

    // get size of polymer chain 
    // int size = NewPol[index].chain.size(); 
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( NewPol[index].chain[deg_poly-1].coords, x, y, z ); 

    // get rid of second to last monomer position from ne list 
    // ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), NewG.PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() ); 

    bool b {false}; 

	size_t tries = 0; 
    for (std::array <int,3>& to_check: ne_list){

    	b = false; 

    	for (const Polymer& pmer: (*PolymerVector)){
    		for (const Particle& p: pmer.chain){

    			if (p.coords == to_check){
    				b = true; 
    				break; 
    			}

    		}

    		if (b) {
    			break;
    		}

    	}

    	if (!b){

    		for (int i{0}; i<deg_poly; i++){

    			if ( i != deg_poly-1 ){
    				NewPol[index].chain[i].coords = (*PolymerVector)[index].chain[i+1].coords ; 
    			}
    			else {
    				NewPol[index].chain[i].coords = to_check; 
    			}
    		}
    		// std::cout << "number of tries is " << tries << std::endl;
  			// std::cout << "max number of tries for forward reptation can be " << ne_list.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl; 
    		break;

    	}

    	else {
    		// std::cout << "Occupied, pick another place to reptate." << std::endl;
    		++tries; 
    	}

    }
	if (tries == ne_list.size() ){
        *IMP_BOOL = false; 
		// std::cout << "There is no place to slither! Will be accepted by default!" << std::endl;
	}

	NewPol[index].ChainToConnectivityMap(); 

    return NewPol; 

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

std::vector <Polymer> BackwardReptation(std::vector <Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){

    std::vector <Polymer> NewPol {*PolymerVector}; 
    int deg_poly = (*PolymerVector)[index].deg_poly; 

    // get size of polymer chain 
    // int size = NewPol[index].chain.size(); 
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( NewPol[index].chain[0].coords, x, y, z ); 

    // get rid of second to last monomer position from ne list 
    // ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), NewG.PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() ); 

    bool b { false }; 

	size_t tries = 0; 
    for (std::array <int,3>& to_check: ne_list){

    	b = false; 

    	for (const Polymer& pmer: (*PolymerVector)){
    		for (const Particle& p: pmer.chain){

    			if (p.coords == to_check){
    				b = true; 
    				break; 
    			}

    		}

    		if (b) {
    			break;
    		}

    	}

    	if (!b){

    		for (int i{0}; i<deg_poly; ++i){

    			if ( i != deg_poly-1 ){
    				NewPol[index].chain[deg_poly-1-i].coords = (*PolymerVector)[index].chain[deg_poly-2-i].coords ; 
    			}
    			else {
    				NewPol[index].chain[deg_poly-1-i].coords = to_check; 
    			}
    		}
    		// std::cout << "number of tries is " << tries << std::endl;
  			// std::cout << "max number of tries for back reptation can be " << ne_list.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl; 
    		break;
    	}

    	else {
    		// std::cout << "Occupied, pick another place to reptate." << std::endl;
    		++tries; 
    	}

    }
	if (tries == ne_list.size() ){
        *IMP_BOOL = false; 
		// std::cout << "There is no place to slither! Will be accepted by default!" << std::endl;
	}

	NewPol[index].ChainToConnectivityMap(); 

    return NewPol; 

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

std::vector <Polymer> Reptation(std::vector<Polymer>* PolymerVector, int index, int x, int y, int z, bool* IMP_BOOL){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return BackwardReptation(PolymerVector, index, x, y, z, IMP_BOOL); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return ForwardReptation(PolymerVector, index, x, y, z, IMP_BOOL); 

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
// NAME OF FUNCTION: MoveChooser 
//
// PARAMETERS: a well-defined Grid Object ie a Grid which has all its attributes set up (correctly)
// 
// WHAT THE FUNCTION DOES: Given a Grid, it will perform a certain Monte Carlo move on the Grid. 
// The move could be from any of the following: 
// 1. End Rotation
// 2. Kink Jump
// 3. Crank Shaft
// 4. Reptation
// 5. Aggressive End rotation
//
// DEPENDENCIES: rng_uniform, CalculateEnergy, EndRotation, KinkJump, CrankShaft, Reptation, IsingFlip  
//
// THE CODE: 

std::vector <Polymer> MoveChooser(std::vector <Polymer>* PolymerVector, int x, int y, int z, bool v, bool* IMP_BOOL){

    int index = rng_uniform(0, static_cast<int>((*PolymerVector).size())-1); 
    std::vector <Polymer> NewPol;

    Grid G_ ; 
    int r = rng_uniform(1, 4);
    switch (r) {
        case (1):
            if (v){
               printf("Performing end rotations.\n"); 
            }
            // 
            NewPol = EndRotation(PolymerVector, index, x, y, z, IMP_BOOL);
            break;     
        
        case (2):
            if (v){
               printf("Performing crank shaft.\n"); 
            }
            NewPol = CrankShaft(PolymerVector, index, x, y, z, IMP_BOOL);
            break; 

        case (3):
            if (v){
               printf("Performing reptation.\n"); 
            }
            NewPol = Reptation(PolymerVector, index, x, y, z, IMP_BOOL); 
            break; 

        case (4):
            if (v){
               printf("Performing kink jump.\n"); 
            }
            NewPol = KinkJump(PolymerVector, index, x, y, z, IMP_BOOL);
            break; 
    }

    return NewPol;
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of MoveChooser. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

