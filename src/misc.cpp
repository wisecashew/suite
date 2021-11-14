#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <random>
#include <chrono>
#include "misc.h"
#include <sstream> 
#include <fstream> 
#include <regex>
#include <tuple> 

const std::vector <int> ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1}; 	// unit directions 
const std::vector <std::vector <int>> drns = {ex, nex, ey, ney, ez, nez};  							// vector of unit directions 

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


//=====================================================

int modified_modulo(int divident, int divisor){
	double midway = static_cast<double>(divisor/2); 
	int result; 
	if (divident%divisor > midway){
		result = (divident%divisor)-divisor; 
		return result; 
	}
	else {
		return (divident%divisor);
	}

}

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







//=====================================================
std::vector <std::vector <int>> HingeSwingDirections(std::vector <int>* HingeToHinge, std::vector <int>* HingeToKink, int x, int y, int z){
    
    (*HingeToHinge).at(0) = modified_modulo((*HingeToHinge).at(0), x); 
    (*HingeToHinge).at(1) = modified_modulo((*HingeToHinge).at(1), y); 
    (*HingeToHinge).at(2) = modified_modulo((*HingeToHinge).at(2), z); 

    (*HingeToKink).at(0) = modified_modulo((*HingeToKink).at(0), x); 
    (*HingeToKink).at(1) = modified_modulo((*HingeToKink).at(1), y); 
    (*HingeToKink).at(2) = modified_modulo((*HingeToKink).at(2), z); 


    std::vector <int> ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1};     // unit directions 
    std::vector <std::vector <int>> drns = {ex, nex, ey, ney, ez, nez};                           // vector of unit directions 


    // get rid of HingeToHinge and its negative from drns 
    std::vector <int> nHingeToHinge = {-(*HingeToHinge).at(0), -(*HingeToHinge).at(1), -(*HingeToHinge).at(2)}; 

    drns.erase(std::remove(drns.begin(), drns.end(), *HingeToHinge), drns.end() ); 
    drns.erase(std::remove(drns.begin(), drns.end(), nHingeToHinge), drns.end() );  
    drns.erase(std::remove(drns.begin(), drns.end(), *HingeToKink), drns.end() ); 

    return drns; 

}


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


//======================================================

//=====================================================
// function to add two vectors 
//$====================================================
std::vector <int> add_vectors(std::vector <int>* v1, std::vector <int>* v2){
	size_t s = (*v1).size(); 
	std::vector <int> v3 (s,0); 
	for (int i{0}; i < s; i++){
		v3.at(i) = (*v2).at(i) + (*v1).at(i);
	}

	return v3; 
}
//======================================================

//=====================================================
// function to subtract two vectors 
//$====================================================
std::vector <int> subtract_vectors(std::vector <int>* v1, std::vector <int>* v2){
	size_t s = (*v1).size(); 
	std::vector <int> v3 (s,0); 
	for (int i{0}; i < s; i++){
		v3.at(i) = (*v1).at(i) -  (*v2).at(i);
	}

	return v3; 
}




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

	return nl; 

}


// ================================================================
// generating a vector of particles given a list of locations
// ================================================================

std::vector <Particle> loc2part (std::vector <std::vector <int>> loc_list){

	std::vector <Particle> pVec; 
	for (std::vector <int> pos: loc_list){
		Particle p (pos); 
		pVec.push_back(p);
	}

	return pVec; 

}

// ===================================================================
// generating a vector of locations from a vector of particles
// ===================================================================
std::vector <std::vector <int> > part2loc (std::vector <Particle> pVec){

	std::vector <std::vector <int> > loc_list; 
	for (Particle p: pVec){
		loc_list.push_back(p.coords);
	}

	return loc_list; 
}

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


// ===================================================================
// calculate energy of a polymeric chain, with edge lengths  
// ===================================================================

int PolymerEnergySolvation(std::vector <Particle> polymer, int x_len, int y_len, int z_len, int intr_energy, int intr_energymm){

	// get all the neighboring sites for the polymer 
	// this only works for a single polymer 

	std::vector <std::vector <int>> poly_ne_list;
	std::vector <std::vector <int>> nl; 
	std::vector <std::vector <int>> loc_list = part2loc(polymer); 
	
	int count = 0; 
	int z = 6; // coordination number in 3d 
	int mm_intrs = 0; 
	for (Particle p: polymer){

		nl = obtain_ne_list(p.coords, x_len, y_len, z_len); 

		for (std::vector <int> v: loc_list){
			
			nl.erase(std::remove(nl.begin(), nl.end(), v), nl.end());
		}


		if (count==0){
			mm_intrs += (z-nl.size()-1); 
			}
		else if (count == polymer.size()-1){
				mm_intrs += (z-nl.size()-1); 
			}
		else {
				mm_intrs += (z-nl.size()-2);
			}
		
		count += 1; 
		poly_ne_list.insert( poly_ne_list.end(), nl.begin(), nl.end() ); // concatenate all neighbors in a location
	}

	// std::cout << "mm intrx energy is " << intr_energymm << std::endl;
	// std::cout << "# of mm intrxs is " << mm_intrs << std::endl;
	// std::cout <<"0.5*mm_intrs*intr_energymm = " << 0.5*mm_intrs*intr_energymm << std::endl;
	int ms_energy = poly_ne_list.size()*intr_energy; 
	int mm_energy = 0.5*mm_intrs*intr_energymm;

	int net_energy = ms_energy + mm_energy;

	return net_energy;

}



// ===================================================================
// extract information
// ===================================================================


std::vector <double> NumberExtractor(std::string s){

	std::vector <double> info; 
	std::stringstream ss (s); 
	std::string temp; 

	double found; 

	while (!(ss.eof()) ) {
		ss >> temp; 

		if (std::stringstream(temp) >> found){
			info.push_back(found); 
		}

		temp = "" ; //not sure what this is for; 
	}

	return info; 

}


std::tuple <std::vector <double>, std::vector<std::vector<double>> > ExtractTopologyFromFile(std::string filename){
    
    std::vector <double> info_vec; 
    std::vector <std::vector <double>> energy_mat; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(filename); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), mat ("ENERGY INTERACTION MATRIX"); 
    bool out_mat = true; 

    for (std::string s: contents){

    	if (std::regex_search(s, x)){
    		std::vector <double> info = NumberExtractor(s); 
    		info_vec.push_back(info.at(0)); 
    		continue; 
    	}

    	else if (std::regex_search(s, y)){
    		std::vector <double> info = NumberExtractor(s); 
    		info_vec.push_back(info.at(0)); 
    		continue; 
    	}

    	else if (std::regex_search(s, z)){
    		std::vector <double> info = NumberExtractor(s); 
    		info_vec.push_back(info.at(0)); 
    		continue; 
    	}

		else if (std::regex_search(s, kT)){
    		std::vector <double> info = NumberExtractor(s); 
    		info_vec.push_back(info.at(0)); 
    		continue; 
    	}

    	else if (std::regex_search (s, mat)){
    		out_mat = false;
    		continue; 
    	}    	

    	else if (!out_mat){
    		energy_mat.push_back(NumberExtractor(s)); 
    		continue; 
    		}

    	else {
    		std::cerr << "ERROR: There is a nonstandard input provided." << std::endl;
    		exit(EXIT_FAILURE); 
    	}

    }



    return std::make_tuple(info_vec, energy_mat);

};





std::vector <std::string> ExtractContentFromFile(std::string filename){
    std::ifstream myfile (filename); 
    std::string mystring; 
    std::vector <std::string> contents; 

    if (myfile.is_open() ){
        while (myfile.good()) {
            std::getline(myfile, mystring); // pipe file's content into stream 
            contents.push_back(mystring); 
        }
    }

    return contents; 

};





