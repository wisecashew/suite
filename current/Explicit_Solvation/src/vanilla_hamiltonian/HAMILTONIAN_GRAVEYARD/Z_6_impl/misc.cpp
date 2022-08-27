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

/* ==================================================
These are some objects I have defined which I tend to use often. 
Helpful definitions which are employed often in the context of the z=6 lattice I am using.
====================================================*/ 

const std::vector <int> ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1}; 	// unit directions 
const std::vector <std::vector <int>> drns = {ex, nex, ey, ney, ez, nez};  							// vector of unit directions 

std::array <int,3> ax{1,0,0}, nax{-1,0,0}, ay{0,1,0}, nay{0,-1,0}, az{0,0,1}, naz{0,0,-1};		// unit directions
std::array <std::array <int,3> ,6> adrns = {ax, nax, ay, nay, az, naz}; 

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


Polymer makePolymer(std::vector <std::array <int,3> > locations, char type_m){
	std::vector <int> pmer_spins; 
    int size_ = locations.size(); 
    for (int i=0; i<size_; i++){
        unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
        std::mt19937 generator(seed); 
        std::uniform_int_distribution<int> distribution (0,1); 
        pmer_spins.push_back(distribution(generator));
    }

    std::vector <Particle> ptc_vec; 

    for (int i=0;i< size_ ; i++ ){
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

	std::cout << "E1 is " << E1 << std::endl;
	std::cout << "E2 is " << E2 << std::endl;

	std::cout << "Probability of acceptance is " << prob << "." << std::endl;
	std::cout << "RNG is " << r << "." << std::endl;
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

void InputParser(int dfreq, int max_iter, std::string solvent_file,
	std::string positions, std::string topology, std::string dfile, 
	std::string efile, std::string mfile, std::string stats_file){

	// if acceptance criterion is NOT CALLED 


    if (dfreq == -1 || max_iter == -1){
        std::cerr << "ERROR: No value for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);
    }

	// if restart is NOT CALLED
	// if (!r){
	//	if (restart_traj != "blank"){
    //        std::cerr << "ERROR: You cannot ask for a trajectory coordinate file with -T without the -r flag. Safeguard against uncontrolled behavior. Exiting..." << std::endl;
    //        exit (EXIT_FAILURE); 
    //    }
    //    if (dfreq == -1 || max_iter == -1 ){
    //        std::cerr << "ERROR: No value for option N (number of accepted MC moves to have) and/or for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
    //       exit (EXIT_FAILURE);
    //    }
    // }
    
    
    if (dfreq == -1 || max_iter == -1 ){
        std::cerr << "ERROR: No value for option N (number of accepted MC moves to have) and/or for option f (frequency of dumping) and/or for option M (maximum number of moves to be performed) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);
    }
    // else if (restart_traj == "blank"){
    //        std::cerr << "ERROR: Need a trajectory file -T if -r option is being called. Exiting..." << std::endl;
    //        exit(EXIT_FAILURE); 
    // }

	// simulation requires an initial "positions" file, a topology file, a coordinate dump, and an energydump file 

    if ( positions== "blank" || topology == "blank" || dfile== "blank" || efile == "blank" || mfile == "blank" || stats_file == "blank" || solvent_file == "blank" ){
        std::cerr << "polymer coords file is " << positions <<",\n topology is " << topology <<",\n polymer coordinate dump file is " << dfile << ",\n energy dump file is " \
        << efile << ",\n orientation file is " << mfile << ",\n move statistics file is " << stats_file << ",\n solvent coords file is " << stats_file << "." << std::endl;
        std::cerr << "ERROR: No value for option p (polymer coordinate file) and/or\nfor option S (solvent coordinate file) and/or\n" <<
        "for option t (energy and geometry file) and/or\nfor option o (name of output dump file) and/or\nfor option e (name of orientation file) and/or\n" <<
        "for option s (name of move stats file) and/or\n for option u (name of energy dump file) was provided. Exiting..." << std::endl;
        exit (EXIT_FAILURE);    
    }

    // set up outputs 

    // if (r){
    //    std::ofstream dump_file(dfile, std::ios::app); 
    //    std::ofstream energy_dump_file (efile, std::ios::app); 
    // }

    
    std::ofstream polymer_dump_file (dfile);
    std::ofstream energy_dump_file (efile);
    std::ofstream orientation_dump_file (mfile); 
    std::ofstream statistics_dump_file (stats_file); 
    std::ofstream solvent_dump_file (solvent_file); 
    

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


bool checkForSolventMonomerOverlap(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolventVector){
    
    std::vector <std::array <int,3>> loc_list;
    loc_list.reserve( (*SolventVector).size() ); 

    for (const Polymer& pmer: (*PolymerVector)) {
        for (const Particle& p: pmer.chain){
            // check if element exists in vector 
                if (std::find(loc_list.begin(), loc_list.end(), p.coords) != loc_list.end() ){
                    std::cerr << "you have a repeated element." << std::endl;
                    return false; 
                }
            
                else{
                    loc_list.push_back(p.coords);  
                }
        }
    }    
    
    for (const Particle& p: (*SolventVector)){

    	if (std::find(loc_list.begin(), loc_list.end(), p.coords) != loc_list.end() ){
            std::cerr << "you have a repeated element. There is a fuck up." << std::endl;
            return false;    		
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

Particle ParticleReporter (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, std::array <int,3> to_check){

	for (const Polymer& pmer: (*Polymers)) {
		for (const Particle& p: pmer.chain){ 

			if (to_check == p.coords){
				return p;
			}

		}
	}

	for (const Particle& p: *Solvent){

		if (to_check == p.coords){
			return p;
		}

	} 
    // (*PolymerVector)[0].printChainCoords();     
    // (*PolymerVector)[1].printChainCoords(); 
    // std::cout << "Coords of solvent: \n";    

    for (const Particle& p: *Solvent){
        print(p.coords); 
    }
    
	std::cout << "Something is profoundly fucked - reporting from ParticleReporter." << std::endl;
    exit(EXIT_FAILURE);
	// Particle p;
	return (*Solvent)[0]; 

}


void ParticleReporter (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, std::pair <char, int>* properties, std::array <int,3>* to_check){


	for ( const Polymer& pmer: (*Polymers)) {
		for ( const Particle& p: pmer.chain) {

			if ((*to_check) == p.coords){
				(*properties).first  = 'm'; 
				(*properties).second = p.orientation; 
				return;
			}
		}
	}

	for (const Particle& p: *Solvent){

		if ( (*to_check) == p.coords){
			(*properties).first  = 's'; 
			(*properties).second = p.orientation; 
			return;
		}

	} 
    // (*PolymerVector)[0].printChainCoords();     
    // (*PolymerVector)[1].printChainCoords(); 
    // std::cout << "Coords of solvent: \n";    

	std::cout << "location being checked is: "; print(*to_check);
	std::cout << "printing out solvent locations: " << std::endl;
    for ( const Particle& p: *Solvent ){

        print(p.coords); 
    }
    
	std::cout << "Something is profoundly fucked - reporting from void ParticleReporter." << std::endl;
    exit(EXIT_FAILURE);
	return ; 


}

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
		for (Particle& p: pmer.chain){ 

			if ((*to_check) == p.coords){
				return true;
			}
		}
	}

	return false; 

}

bool MonomerReporter (std::vector <Polymer>* Polymers, std::array <int,3>* check_1, std::array <int,3>* check_2){

	for (const Polymer& pmer: (*Polymers)) {
      		for (const Particle& p: pmer.chain){

      			if (p.coords == (*check_1) || p.coords == (*check_2) ){	 
      				return true; 
      			}
      		}
      	}
    return false; 
}


bool MonomerNeighborReporter ( std::vector <Polymer>* Polymers, std::array <int,3>* to_check, int x, int y, int z){

	std::array <int,3> diff = {0,0,0}; 
	for (Polymer& pmer: (*Polymers) ){
		for (Particle& p: pmer.chain) {

			diff = subtract_arrays ( &p.coords, to_check ); 
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

double CalculateEnergy(std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, double Emm_a, double Emm_n, double Ems_a, double Ems_n, double* mm_aligned, double* mm_naligned, int* ms_aligned, int* ms_naligned){
    double Energy {0.0};
    (*mm_aligned)  = 0; 
    (*mm_naligned) = 0; 
    (*ms_aligned)  = 0;  
    (*ms_naligned) = 0;
    // polymer-polymer interaction energies 
    
    std::pair <char, int> properties ( ' ' , -1 );

    for (const Polymer& pmer: (*Polymers)) {
        for (const Particle& p: pmer.chain){
            std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(p.coords, x, y, z); // get neighbor list 
            
            for ( std::array <int, 3>& loc: ne_list){
            	
            	ParticleReporter (Polymers, Solvent, &properties, &loc);

            	if ( properties.first == 'm'){
                    
            		if ( properties.second == p.orientation ){
                        (*mm_aligned) += 0.5; 
            			Energy += 0.5*Emm_a; 
            		}
            		else {
                        (*mm_naligned) += 0.5; 
            			Energy += 0.5*Emm_n; 
            		}
            	}

            	else { // particle is of type solvent 

            		if ( properties.second == p.orientation ){
                        (*ms_aligned)  += 1;
            			Energy += Ems_a;
            		}
            		else {
                        (*ms_naligned) += 1; 
            			Energy += Ems_n; 
            		}
            	}
            }
        }
    }
    // std::cout << "Energy is " << Energy << "." << std::endl;
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
            dump_file << p.orientation << " | ";
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


void dumpPositionOfSolvent(std::vector <Particle>* SolventVector, int step, std::string filename){

	std::ofstream dump_file(filename, std::ios::app); 
    dump_file <<"Dumping coordinates at step " << step << ".\n";
    
    // dump_file <<"Dumping coordinates at step " << step << ".\n";
    int count = 0; 
        
    dump_file <<"Dumping coordinates of solvent # " << count << ".\n";
    dump_file<<"START" << "\n";
    
    for (const Particle& p: *SolventVector){
    	dump_file<<"Orientation: " << p.orientation <<", ";
        for (int i: p.coords){
            dump_file << i << " | "; 
        }
        dump_file << "\n"; 
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


void dumpEnergy (double sysEnergy, int step, double mm_aligned, double mm_naligned, int ms_aligned, int ms_naligned, std::string filename){
    std::ofstream dump_file(filename, std::ios::app); 
    // std::ostringstream os; 
    
    dump_file << sysEnergy << " | " << mm_aligned+mm_naligned << " | " << mm_aligned << " | " << mm_naligned << " | " \
            << ms_aligned+ms_naligned << " | " << ms_aligned << " | " << ms_naligned << " | " << step << "\n";
    
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
    
    std::ofstream dump_file (stats_file, std::ios::app); 
    dump_file << "For step " << step << ".\n";
    

    dump_file << "End rotations                    - attempts: " << (*attempts)[0] <<", acceptances: " << (*acceptances)[0] << ", acceptance fraction: " << static_cast<double>((*acceptances)[0])/static_cast<double>((*attempts)[0]) << ".\n"; 
    dump_file << "Kink jumps                       - attempts: " << (*attempts)[1] <<", acceptances: " << (*acceptances)[1] << ", acceptance fraction: " << static_cast<double>((*acceptances)[1])/static_cast<double>((*attempts)[1]) << ".\n"; 
    dump_file << "Crank shafts                     - attempts: " << (*attempts)[2] <<", acceptances: " << (*acceptances)[2] << ", acceptance fraction: " << static_cast<double>((*acceptances)[2])/static_cast<double>((*attempts)[2]) << ".\n"; 
    dump_file << "Reptation                        - attempts: " << (*attempts)[3] <<", acceptances: " << (*acceptances)[3] << ", acceptance fraction: " << static_cast<double>((*acceptances)[3])/static_cast<double>((*attempts)[3]) << ".\n"; 
    dump_file << "Chain regrowth                   - attempts: " << (*attempts)[4] <<", acceptances: " << (*acceptances)[4] << ", acceptance fraction: " << static_cast<double>((*acceptances)[4])/static_cast<double>((*attempts)[4]) << ".\n"; 
    dump_file << "Solvent orientation flips        - attempts: " << (*attempts)[5] <<", acceptances: " << (*acceptances)[5] << ", acceptance fraction: " << static_cast<double>((*acceptances)[5])/static_cast<double>((*attempts)[5]) << ".\n"; 
    dump_file << "Single solvent orientation flips - attempts: " << (*attempts)[6] <<", acceptances: " << (*acceptances)[6] << ", acceptance fraction: " << static_cast<double>((*acceptances)[6])/static_cast<double>((*attempts)[6]) << ".\n"; 
    dump_file << "Polymer orientation flips        - attempts: " << (*attempts)[7] <<", acceptances: " << (*acceptances)[7] << ", acceptance fraction: " << static_cast<double>((*acceptances)[7])/static_cast<double>((*attempts)[7]) << ".\n"; 
    dump_file << "Local polymer orientation flips  - attempts: " << (*attempts)[8] <<", acceptances: " << (*acceptances)[8] << ", acceptance fraction: " << static_cast<double>((*acceptances)[8])/static_cast<double>((*attempts)[8]) << ".\n"; 

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
void dumpOrientation( std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int step, std::string filename, int x, int y, int z ) {
    // std::cout<< "just inside dumpO..."<<std::endl; 
    std::ofstream dump_file (filename, std::ios::app); 
    dump_file << "START for Step " << step << ".\n";
    std::pair <char, int> properties (' ' , -1); 
    for ( const Polymer& pmer: (*Polymers) ) {
        for ( const Particle& p: pmer.chain ) {
            // std::cout << "particle is at "; print(p.coords);  
            dump_file << p.orientation << " | ";
            std::array <std::array<int,3> ,6> ne_list = obtain_ne_list (p.coords, x, y, z) ;
            // std::cout << "neighbor list: " << std::endl;
            // print(ne_list); 
            for ( std::array <int,3>& ne: ne_list) {
                // std::cout <<"Particle positions: "; 
                // print (ne) ;
                ParticleReporter (Polymers, Solvent, &properties, &ne);
                // std::cout << "Reported~\n"; 
                if (properties.first == 's'){
                    dump_file << properties.second << " | ";  
                } 
            }
            dump_file << "\n"; 
        } 
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

void TailRotation (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){


    // get the neighborlist of particle at index 1 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[0].coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[1].coords;
    // std::vector <Polymer> NewPol {*Polymers};
    
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

	(*rweight) = 0.0; 

	// find all spots that are available 
	std::vector <std::array <int,3>> idx_v; 
	for (std::array <int,3>& to_rot: ne_list){
		if ( to_rot == loc_1){
			continue; 
		}

		else if ( to_rot == (*Polymers)[index].chain[2].coords ) {
			continue; 
		}

		else if ( ! ( MonomerReporter(Polymers, &to_rot) ) ){
			(*rweight) += 1;
			idx_v.push_back( to_rot );  
		}
	}

    if ( (*rweight) == 0 ){
    	// if nothing can be done, just get out. 
    	*IMP_BOOL = false;
    	return; 
    }

    // otherwise, let's fuckin go

	int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 

	for (Particle& p: (*Solvent)){
		if (p.coords == idx_v[ r ]){
			p.coords = loc_0; 
			(*Polymers)[index].chain[0].coords = idx_v[ r ];  
			break;
		}
	}

	// (*Polymers)[index].ChainToConnectivityMap(); 
	(*rweight) = (*rweight)/6.0; 


	// MOVE UPDATING POLYMER SKELETON HAS BEEN PERFORMED. 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// now to perform the *Solvent update. 

	// I have to make sure I am updating Solvent such that the Solvent is still surrounding polymer
	// find all neighbors of the new structure, keeping in mind only one monomer segment was moved.

	// (1) From the initial monomer position, find all the neighbors of the monomer that were solvent. Do the same for the perturbed monomer. 
	// keeping in mind initial monomer position and perturbed monomer position can never be neighbors.

	std::array  <std::array<int,3>, 6> ne_before = obtain_ne_list ( loc_0, x, y, z );                    					  // neighborlist of the monomer before perturbation 
	std::array  <std::array<int,3>, 6> ne_after = obtain_ne_list (idx_v[r], x, y, z);  // neighborlist of the monomer after perturbation 
	std::vector <std::array <int,3>  > solvent_to_delete, solvent_to_add; 
	solvent_to_delete.reserve(5); solvent_to_add.reserve(5); 

	// make sure there is no monomer location in ne_list... 
	for (std::array <int,3>& ne: ne_before) {

		if ( ne == (*Polymers)[index].chain[1].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location of a monomer, continue 
			continue;
		}

		else {
			// if it is a neighbor of another monomer, this will be kept in solvent. 
			if ( MonomerNeighborReporter (Polymers, &ne, x, y, z) ){
				continue; 
			}
			else {
				// if it is a solvent particle but not neighboring any solvent in the new structure, it must NOT BE PRESENT in *Solvent.
				// that is the purpose of solvent_ne!!! Figuring out all the elements that must be cut out from *Solvent. 
				solvent_to_delete.push_back( ne );
			}
		}
	}



	// make sure there is no monomer location in ne_list_new... 
	for (std::array <int,3>& ne: ne_after) {

		if ( ne == (*Polymers)[index].chain[1].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location of a monomer, continue 
			continue;
		}

		else {
			// if it is a solvent particle AND is a neighbor of a monomer in the new structure, IT MUST BE PRESENT in *Solvent. 
			// that is the purpose of solvent_ne_new!!! Figuring out all the elements that NEED TO BE ADDED to *Solvent, in case they are not already there.   
			solvent_to_add.push_back (ne); 
		}
	}


	// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
	// deleting elements that MUST NOT be in *Solvent... 

	
	int s_idx = 0; 
	for ( const std::array<int,3>& s_loc: solvent_to_delete ) {
		s_idx = 0; 
		for ( Particle& p: *Solvent ){
			
			if ( p.coords == s_loc ){
				(*Solvent).erase ( (*Solvent).begin() + s_idx ); 
				break; 
			}
			s_idx += 1; 

		}

	}

	// (3) now that *Solvent has been cleaned out, time to add new things to it 
	bool present = false; 
	for ( const std::array <int,3>& s_loc: solvent_to_add ){
		present = false; 
		for ( Particle& p: *Solvent ){
			if ( p.coords == s_loc ){
				present = true;
				break; 
			}
		}

		if ( present ){
			continue; 
		}

		else {
			Particle temp ( s_loc, 's', 0 );
			(*Solvent).push_back( temp ); 
		}

	}
	
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

void HeadRotation (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    // get the neighborlist of particle at index 1 
    int dop = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc_0 = (*Polymers)[index].chain[dop-1].coords; 
    std::array <int,3> loc_1 = (*Polymers)[index].chain[dop-2].coords;
    
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

    (*rweight) = 0.0; 
	
    // find all spots that are available 
    std::vector <std::array <int,3>> idx_v; 
    for (std::array <int,3>& to_rot: ne_list){

    	if (to_rot == loc_1){
    		continue; 
    	}

    	else if ( to_rot == (*Polymers)[index].chain[dop-3].coords ) {
			continue; 
		}
    	else if ( ! (MonomerReporter(Polymers, &to_rot) ) ){
    		// std::cout << "a free location is "; print(to_rot); 
    		(*rweight) += 1; 
    		idx_v.push_back(to_rot); 
    	}
    }

    if ( (*rweight) == 0){
    	// if rweight is zero, return 
    	*IMP_BOOL = false;
    	return; 
    }

    // if rweight is not zero, lets fuckin go

	int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 
	
	// std::cout << "Suggested rot position is "; print (idx_v [r]); 
	for (Particle& p: (*Solvent)){
		if (p.coords == idx_v[ r ]){
			p.coords = loc_0;
			(*Polymers)[index].chain[dop-1].coords = idx_v[ r ]; 
			break; 
		}
	}

	// (*Polymers)[index].ChainToConnectivityMap(); 
	(*rweight) = (*rweight)/6.0; 


    // MOVE UPDATING POLYMER SKELETON HAS BEEN PERFORMED. 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// now to perform the *Solvent update. 

	// I have to make sure I am updating Solvent such that the Solvent is still surrounding polymer
	// find all neighbors of the new structure, keeping in mind only one monomer segment was moved.

	// (1) From the initial monomer position, find all the neighbors of the monomer that were solvent. Do the same for the perturbed monomer. 
	// keeping in mind initial monomer position and perturbed monomer position can never be neighbors.

	std::array <std::array <int,3>, 6> ne_before = obtain_ne_list ( loc_0, x, y, z); 										 // neighborlist of the monomer before perturbation 
	std::array <std::array <int,3>, 6> ne_after = obtain_ne_list (idx_v[r], x, y, z); 										// neighborlist of the monomer after perturbation

	std::vector <std::array <int,3> > solvent_to_delete, solvent_to_add; 
	solvent_to_delete.reserve(5); solvent_to_add.reserve(5); 

	// make there is no monomer location in ne_list...
	for ( std::array <int,3>& ne: ne_before){

		if ( ne == (*Polymers)[index].chain[dop-2].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location of a monomer, continue 
			continue; 
		}

		else {
			// if it is a neighbor of another monomer, this will be kepy in solvent. 
			if ( MonomerNeighborReporter (Polymers, &ne, x, y, z) ){
				continue; 
			} 
			else {
				// if it is a solvent particle but not neighboring any solvent in the new structure, it must NOT BE PRESENT in *Solvent.
				// that is the purpose of solvent_ne!!! Figuring out all the elements that must be cut out from *Solvent. 
				solvent_to_delete.push_back( ne );
			}
		}
	}


	// make sure there is no monomer location in this... 
	for ( std::array <int,3>& ne: ne_after ) {

		if ( ne == (*Polymers)[index].chain[dop-2].coords ){
			continue;
		}

		else if (MonomerReporter (Polymers, &ne) ){
			// if ne is alocation of a monomer, continue
			continue;
		}

		else {
			// if it is a solvent particle AND is a neighbor of a monomer in the new structure, IT MUST BE PRESENT in *Solvent. 
			// that is the purpose of solvent_ne_new!!! Figuring out all the elements that NEED TO BE ADDED to *Solvent, in case they are not already there.   
			solvent_to_add.push_back (ne); 
		}
	}

	// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
	// deleting elements that MUST NOT be in *Solvent... 

	int s_idx = 0; 
	for ( const std::array <int,3> s_loc: solvent_to_delete ) {
		s_idx = 0; 
		for ( Particle& p: *Solvent ){
			if ( p.coords == s_loc ){
				(*Solvent).erase ( (*Solvent).begin() + s_idx );
				break;
			}
			
			s_idx += 1; 

		}

	}

	// (3) now that *Solvent has been cleaned out, time to add new things to it 
	bool present = false; 
	for ( const std::array <int,3>& s_loc: solvent_to_add ){
		present = false;
		for ( Particle& p: *Solvent ){
			if ( p.coords == s_loc ){
				present = true; 
				break; 
			}
		}

		if ( present ) {
			// if it is already there in solvent, dont do anything. move to the next position. 
			continue; 
		}

		else {
			Particle temp ( s_loc, 's', 0 );
			(*Solvent).push_back ( temp ); 
		}

	}


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

void EndRotation (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        TailRotation (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight); 
        return; 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        HeadRotation (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight); 
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

void KinkJump (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

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

        std::array <int,3> d2       = subtract_arrays( &( (*Polymers) [index].chain[idx+2].coords), &( (*Polymers) [index].chain[idx+1].coords) ); 

        std::array <int,3> to_check = add_arrays     ( &( (*Polymers) [index].chain[idx].coords ), &d2); 
        
        impose_pbc(&to_check, x, y, z); 

        // find locations where kink can jump to... 

        if ( ! ( MonomerReporter (Polymers, &to_check) ) ){
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
	// (*Polymers) [index].chain[idx_v[r]+1].coords = pos_v[r]; 

	// kink jump spot is pos_v[r]
	// initial location of monomer that is jumping is... 
	std::array <int,3> loc_0 = (*Polymers)[index].chain[idx_v[r]+1].coords; 

	for (Particle& p: (*Solvent)){
		if (p.coords == pos_v[r]){
			p.coords = (*Polymers)[index].chain[idx_v[r]+1].coords; 
			(*Polymers) [index].chain[idx_v[r]+1].coords = pos_v[r]; 
			break;
		}
	}

	// (*Polymers)[index].ChainToConnectivityMap(); 
	(*rweight) = (*rweight)/( k_idx.size() ); 


	// MOVE UPDATING POLYMER SKELETON HAS BEEN PERFORMED. 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// now to perform the *Solvent update. 

	// I have to make sure I am updating Solvent such that the Solvent is still surrounding polymer
	// find all neighbors of the new structure, keeping in mind only one monomer segment was moved.

	// (1) From the initial monomer position, find all the neighbors of the monomer that were solvent. Do the same for the perturbed monomer. 
	// keeping in mind initial monomer position and perturbed monomer position can never be neighbors.

	std::array < std::array <int,3>, 6> ne_before     = obtain_ne_list ( loc_0, x, y, z );                    	// neighborlist of the monomer before perturbation
	std::array < std::array <int,3>, 6> ne_after 	  = obtain_ne_list ( pos_v[r], x, y, z); 						// neighborlist of the monomer after perturbation
	std::vector <std::array <int,3> > solvent_to_delete, solvent_to_add; 
	solvent_to_delete.reserve(4); solvent_to_add.reserve(4); 

	// make sure there is no monomer location in ne_list... 
	for (std::array <int,3>& ne: ne_before) {

		if ( ne == (*Polymers)[index].chain[idx_v[r]].coords || ne == (*Polymers)[index].chain[idx_v[r]+2].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location of a monomer, continue 
			continue;
		}

		else {
			// if it is a neighbor of another monomer, this will be kept in solvent. 
			if ( MonomerNeighborReporter (Polymers, &ne, x, y, z) ){
				continue; 
			}
			else {
				// if it is a solvent particle but not neighboring any solvent in the new structure, it must NOT BE PRESENT in *Solvent.
				// that is the purpose of solvent_ne!!! Figuring out all the elements that must be cut out from *Solvent. 
				solvent_to_delete.push_back( ne );
			}
		}
	}

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// make sure there is no monomer location in ne_list_new... 
	for (std::array <int,3>& ne: ne_after) {

		if ( ne == (*Polymers)[index].chain[idx_v[r]].coords || ne == (*Polymers)[index].chain[idx_v[r]+2].coords ){
			continue;
		}

		else if (MonomerReporter (Polymers, &ne) ){
			// if ne is a location of a monomer, continue 
			continue;
		}

		else {
			// if it is a solvent particle AND is a neighbor of a monomer in the new structure, IT MUST BE PRESENT in *Solvent. 
			// that is the purpose of solvent_ne_new!!! Figuring out all the elements that NEED TO BE ADDED to *Solvent, in case they are not already there.   
			solvent_to_add.push_back (ne); 
			
		}
	}

	// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
	// deleting elements that MUST NOT be in *Solvent... 

	int s_idx = 0; 
	for ( const std::array<int,3>& s_loc: solvent_to_delete ) {
		s_idx = 0; 
		for ( Particle& p: *Solvent ){
			
			if ( p.coords == s_loc ){
				(*Solvent).erase ( (*Solvent).begin() + s_idx ); 
				break; 
			}
			s_idx += 1; 
		}
	}

	// (3) now that *Solvent has been cleaned out, time to add new things to it 
	bool present = false; 
	for ( const std::array <int,3>& s_loc: solvent_to_add ){
		present = false; 
		for ( Particle& p: *Solvent ){
			// if it is already present in *Solvent, breakout and continue to next location
			if ( p.coords == s_loc ){
				present = true;
				break; 
			}
		}

		if ( present ){
			continue; 
		}

		else {
			Particle temp ( s_loc, 's', 0 );
			(*Solvent).push_back( temp ); 
		}

	}

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

void CrankShaft (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

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

	std::array <int,3> HingeToKink  = subtract_arrays(&( (*Polymers) [index].chain[idx+2].coords), &( (*Polymers) [index].chain[idx+3].coords) ); 
    std::array <int,3> HingeToHinge = subtract_arrays(&( (*Polymers) [index].chain[idx+3].coords ), &( (*Polymers) [index].chain[idx].coords ) );
    std::array <std::array <int,3>,3> drns = HingeSwingDirections (& (HingeToHinge), &(HingeToKink), x, y, z); 

	// std::array <int,3> d1 = drns[choice];

	for (std::array <int,3>& d: drns ){

        std::array <int,3> to_check_1 = add_arrays ( &( (*Polymers) [index].chain[idx].coords), &d ); 
        std::array <int,3> to_check_2 = add_arrays ( &( (*Polymers) [index].chain[idx+3].coords), &d ); 
        impose_pbc(&to_check_1, x, y, z); 
        impose_pbc(&to_check_2, x, y, z);   

		
      	// check if site is unoccupied 
        if ( ! ( MonomerReporter (Polymers, &to_check_1, &to_check_2) ) ) {
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

	// std::cout << "cranked pos 1 is "; print(pos_v1[r]); 
	// std::cout << "cranked pos 2 is "; print(pos_v2[r]); 
		
	int d = 0; 

	std::array <int,3> loc_1 = (*Polymers)[index].chain[idx+1].coords;
	std::array <int,3> loc_2 = (*Polymers)[index].chain[idx+2].coords;

	for (Particle& p: (*Solvent)){
		if (p.coords == pos_v1[r]){
			++d; 
			p.coords = (*Polymers)[index].chain[idx+1].coords; 
			(*Polymers)[index].chain[idx+1].coords = pos_v1[r];		
			if (d==2){
				break;
			}
		}
		else if (p.coords == pos_v2[r]){
			++d; 
			p.coords = (*Polymers)[index].chain[idx+2].coords; 
			(*Polymers)[index].chain[idx+2].coords = pos_v2[r];
			if (d==2){
				break;
			}
		}
	}
	// (*Polymers) [index].ChainToConnectivityMap();
	(*rweight) = (*rweight)/3.0 ; 
  	

	// MOVE UPDATING POLYMER SKELETON HAS BEEN PERFORMED. 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// now to perform the *Solvent update. 

	// I have to make sure I am updating Solvent such that the Solvent is still surrounding polymer
	// find all neighbors of the new structure, keeping in mind only one monomer segment was moved.

	// (1) From the initial monomer position, find all the neighbors of the monomer that were solvent. Do the same for the perturbed monomer. 
	// you are going to have to repeat this procedure twice because you have to swaps taking place 

    
	std::array < std::array <int,3>, 6> ne_before_1     = obtain_ne_list ( loc_1, x, y, z );                    		// neighborlist of monomer1 before perturbation
	std::array < std::array <int,3>, 6> ne_after_1      = obtain_ne_list ( pos_v1[r], x, y, z); 						// neighborlist of monomer1 after perturbation
	std::array < std::array <int,3>, 6> ne_before_2     = obtain_ne_list ( loc_2, x, y, z );                    		// neighborlist of monomer2 before perturbation
	std::array < std::array <int,3>, 6> ne_after_2      = obtain_ne_list ( pos_v2[r], x, y, z); 						// neighborlist of monomer2 after perturbation
	std::vector <std::array <int,3> > solvent_to_delete, solvent_to_add; 
	solvent_to_delete.reserve(8); solvent_to_add.reserve(8);

	// make sure there is no monomer location in ne_list_1...
	for ( std::array <int,3>& ne: ne_before_1 ) {

		if ( ne == (*Polymers)[index].chain[idx].coords || ne == loc_2 ) {
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ) {
			continue;
		}

		else {
			if ( MonomerNeighborReporter (Polymers, &ne, x, y, z) ){
				continue;
			}
			else{
				solvent_to_delete.push_back( ne );
			}
		}
	}

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// make sure there is no monomer location in ne_list_new_1

	for ( std::array <int,3>& ne: ne_after_1){

		if ( ne == (*Polymers)[index].chain[idx].coords || ne == (*Polymers)[index].chain[idx+2].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ) {
			continue;
		}

		else {
			solvent_to_add.push_back( ne ); 
			continue;
		}

	}

	// make sure there is no monomer location in ne_list_2...
	for ( std::array <int,3>& ne: ne_before_2 ) {

		if ( ne == (*Polymers)[index].chain[idx+3].coords || ne == loc_1 ) {
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ) {
			continue;
		}

		else {
			if ( MonomerNeighborReporter (Polymers, &ne, x, y, z) ){
				continue;
			}
			else{
				solvent_to_delete.push_back( ne );
			}
		}
	}

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// make sure there is no monomer location in ne_list_new_1

	for ( std::array <int,3>& ne: ne_after_2){

		if ( ne == (*Polymers)[index].chain[idx+3].coords || ne == (*Polymers)[index].chain[idx+1].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ) {
			continue;
		}

		else {
				solvent_to_add.push_back( ne ); 
				continue;
		}

	}

	// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
	// deleting elements that MUST NOT be in *Solvent... 

	int s_idx = 0; 
	for ( const std::array<int,3>& s_loc: solvent_to_delete ) {
		s_idx = 0; 
		for ( Particle& p: *Solvent ){
			
			if ( p.coords == s_loc ){
				(*Solvent).erase ( (*Solvent).begin() + s_idx ); 
				break; 
			}
			s_idx += 1; 

		}

	}

	// (3) now that *Solvent has been cleaned out, time to add new things to it 
	bool present = false; 
	for ( const std::array <int,3>& s_loc: solvent_to_add ){
		present = false; 
		for ( Particle& p: *Solvent ){
			// if it is already present in *Solvent, breakout and continue to next location
			if ( p.coords == s_loc ){
				present = true;
				break; 
			}
		}

		if ( present ){
			continue; 
		}

		else {
			Particle temp ( s_loc, 's', 0 );
			(*Solvent).push_back( temp ); 
		}

	}	

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

void ForwardReptation (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    int deg_poly = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc0 = (*Polymers)[index].chain[0].coords; 
    std::array <int,3> locf = (*Polymers)[index].chain[deg_poly-1].coords;
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( locf, x, y, z ); 
    std::vector <std::array<int,3>> idx_v; 
    idx_v.reserve(6); 

    (*rweight) = 0; 

    for (std::array <int,3>& to_check: ne_list){

    	if ( to_check == (*Polymers) [index].chain[deg_poly-2].coords ){
    		continue; 
    	}
        else if ( to_check == (*Polymers) [index].chain[0].coords ) {
            (*rweight) += 1; 
            idx_v.push_back(to_check); 
        }

		else if ( ! (MonomerReporter(Polymers, &to_check) )){
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
	// std::cout << "In acceptance land." << std::endl; 
	// std::cout << "location of new spot is "; print(idx_v[r]); 
	
	for (int i{0}; i<deg_poly; ++i){

		if ( i != deg_poly-1 ){
			(*Polymers)[index].chain[i].coords = (*Polymers)[index].chain[i+1].coords ; 
		}
		else {
			(*Polymers)[index].chain[i].coords = idx_v[r]; 
		}
	}
	
	if ( loc0 != idx_v[r]){
		for (Particle& p: (*Solvent)){
			if (p.coords == idx_v[r]){
				p.coords = loc0; 
				break;
			}
		}
		// (*Polymers)[index].ChainToConnectivityMap();
		(*rweight) = (*rweight)/6; 
	}

	else {
		// (*Polymers)[index].ChainToConnectivityMap();
		(*rweight) = (*rweight)/6; 
		return; 
	}

	
	// MOVE UPDATING POLYMER SKELETON HAS BEEN PERFORMED. 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// now to perform the *Solvent update. 

	// I have to make sure I am updating Solvent such that the Solvent is still surrounding polymer
	// find all neighbors of the new structure, keeping in mind only one monomer segment was moved.

	// (1) From the initial monomer position, find all the neighbors of the monomer that were solvent. 
	// the first and final indexed monomers are the only ones who have done any moving. 
	// their initial and changed environments need most consideration. 

	std::array < std::array <int,3>, 6> ne_before_id0 = obtain_ne_list ( loc0, x, y, z );
	std::array < std::array <int,3>, 6> ne_after_idf  = obtain_ne_list ( (*Polymers)[index].chain[deg_poly-1].coords, x, y, z);
	std::vector < std::array <int,3> > solvent_to_delete, solvent_to_add; 
	solvent_to_delete.reserve(5); solvent_to_add.reserve(5); 

	// checking particles around old position of index0
	for ( std::array <int,3>& ne: ne_before_id0 ) {
		// obviously the current zeroth index position has a monomer on it right now 
		if ( ne == (*Polymers)[index].chain[0].coords ){
			continue; 
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location occupied by a monomer, move on
			continue;
		}

		else {
			// if it is a neighbor of another monomer, then it won't be deleted
			if (MonomerNeighborReporter (Polymers, &ne, x, y, z) ){
				continue;
			}
			else {
				// if it is a solvent particle but not neighboring any solvent in the new structure, it must NOT BE PRESENT in *Solvent.
				// that is the purpose of solvent_ne!!! Figuring out all the elements that must be cut out from *Solvent. 
				solvent_to_delete.push_back(ne);
			}
		}
	}

	// checking particles around new position of indexf
	for ( std::array <int,3>& ne: ne_after_idf ){
		// obviously the (f-1)th index monomer is on the old position of the current fth index 
		if ( ne == (*Polymers)[index].chain[deg_poly-2].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location occupied by a monomer, move on
			continue;
		}

		else {
			// if it is currently a neighbor of another monomer in the present state, IT MUST BE PRESENT in *Solvent.
			// that is the purpose of solvent_ne_new!!! Figuring out all the elements that NEED TO BE ADDED to *Solvent, in case they are not already there.   
			solvent_to_add.push_back(ne); 
		}
	}

	// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
	// deleting elements that MUST NOT be in *Solvent... 

	int s_idx = 0; 
	for ( const std::array<int,3>& s_loc: solvent_to_delete ) {
		s_idx = 0; 
		for ( Particle& p: *Solvent ){
			
			if ( p.coords == s_loc ){
				(*Solvent).erase ( (*Solvent).begin() + s_idx ); 
				break; 
			}
			s_idx += 1; 
		}
	}


	// (3) now that *Solvent has been cleaned out, time to add new things to it 
	bool present = false; 
	for ( const std::array <int,3>& s_loc: solvent_to_add ){
		present = false; 
		for ( Particle& p: *Solvent ){
			// if it is already present in *Solvent, breakout and continue to next location
			if ( p.coords == s_loc ){
				present = true;
				break; 
			}
		}

		if ( present ){
			continue; 
		}

		else {
			Particle temp ( s_loc, 's', 0 );
			(*Solvent).push_back( temp ); 
		}

	}

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


void BackwardReptation (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    int deg_poly = (*Polymers)[index].deg_poly; 
    std::array <int,3> loc0 = (*Polymers)[index].chain[deg_poly-1].coords; 
    std::array <int,3> locf = (*Polymers)[index].chain[0].coords;

    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( locf, x, y, z ); 
    std::vector <std::array <int,3>> idx_v; 

    (*rweight) = 0.0; 

     
    for (std::array <int,3>& to_check: ne_list){

    	if ( to_check == (*Polymers) [index].chain[1].coords){
    		continue; 
    	}

        else if ( to_check == (*Polymers) [index].chain[deg_poly-1].coords ) {
            (*rweight) += 1;
            idx_v.push_back(to_check); 
        }

    	else if ( ! (MonomerReporter(Polymers, &to_check) ) ){
    		(*rweight) += 1; 
    		idx_v.push_back(to_check); 
    	}
    }

    if ( (*rweight) == 0 ){
    	*IMP_BOOL = false; 
    	return;
    }

	int r = rng_uniform( 0, idx_v.size()-1 );

	for (int i{0}; i <deg_poly; ++i){
		if ( i != deg_poly-1 ){
			(*Polymers) [index].chain[deg_poly-1-i].coords = (*Polymers)[index].chain[deg_poly-2-i].coords;
		}
		else {
			(*Polymers) [index].chain[deg_poly-1-i].coords = idx_v[r]; 
		}
	}

	if ( loc0 != idx_v[r] ){
		for (Particle& p: (*Solvent)) {
			if (p.coords == idx_v[r]){
				p.coords = loc0; //(*Polymers)[index].chain[deg_poly-1].coords;
				// print(loc0);
				break;
			}
		}
		// (*Polymers)[index].ChainToConnectivityMap(); 
		(*rweight) = (*rweight)/6; 
	}

	else {
		// (*Polymers)[index].ChainToConnectivityMap(); 
		(*rweight) = (*rweight)/6; 
		return; 
	}

	// MOVE UPDATING POLYMER SKELETON HAS BEEN PERFORMED. 
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// now to perform the *Solvent update. 

	// I have to make sure I am updating Solvent such that the Solvent is still surrounding polymer
	// find all neighbors of the new structure, keeping in mind only one monomer segment was moved.

	// (1) From the initial monomer position, find all the neighbors of the monomer that were solvent. 
	// the first and final indexed monomers are the only ones who have done any moving. 
	// their initial and changed environments need most consideration. 

	std::array < std::array <int,3>, 6> ne_before_idf = obtain_ne_list ( loc0, x, y, z );
	std::array < std::array <int,3>, 6> ne_after_id0  = obtain_ne_list ( (*Polymers)[index].chain[0].coords, x, y, z);
	std::vector < std::array <int,3> > solvent_to_delete, solvent_to_add; 
	solvent_to_delete.reserve(5); solvent_to_add.reserve(5); 

	// checking particles around old position of index0
	for ( std::array <int,3>& ne: ne_before_idf ) {
		// obviously the current final-th index position has a monomer on it right now 
		if ( ne == (*Polymers)[index].chain[deg_poly-1].coords ){
			continue; 
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location occupied by a monomer, move on
			continue;
		}

		else {
			// if it is a neighbor of another monomer, then it won't be deleted
			if (MonomerNeighborReporter (Polymers, &ne, x, y, z) ){
				continue;
			}
			else {
				// if it is a solvent particle but not neighboring any solvent in the new structure, it must NOT BE PRESENT in *Solvent.
				// that is the purpose of solvent_ne!!! Figuring out all the elements that must be cut out from *Solvent. 
				solvent_to_delete.push_back(ne);
			}
		}
	}

	// checking particles around new position of indexf
	for ( std::array <int,3>& ne: ne_after_id0 ){
		// obviously the (f-1)th index monomer is on the old position of the current fth index 
		if ( ne == (*Polymers)[index].chain[1].coords ){
			continue;
		}

		else if ( MonomerReporter (Polymers, &ne) ){
			// if ne is a location occupied by a monomer, move on
			continue;
		}

		else {
			// if it is currently a neighbor of another monomer in the present state, IT MUST BE PRESENT in *Solvent.
			// that is the purpose of solvent_ne_new!!! Figuring out all the elements that NEED TO BE ADDED to *Solvent, in case they are not already there.   
			solvent_to_add.push_back(ne); 
		}
	}

	
	// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
	// deleting elements that MUST NOT be in *Solvent... 

	int s_idx = 0; 
	for ( const std::array<int,3>& s_loc: solvent_to_delete ) {
		s_idx = 0; 
		for ( Particle& p: *Solvent ){
			
			if ( p.coords == s_loc ){
				(*Solvent).erase ( (*Solvent).begin() + s_idx ); 
				break; 
			}
			s_idx += 1; 
		}
	}



	// (3) now that *Solvent has been cleaned out, time to add new things to it 
	bool present = false; 
	for ( const std::array <int,3>& s_loc: solvent_to_add ){
		present = false; 
		for ( Particle& p: *Solvent ){
			// if it is already present in *Solvent, breakout and continue to next location
			if ( p.coords == s_loc ){
				present = true;
				break; 
			}
		}

		if ( present ){
			continue; 
		}

		else {
			Particle temp ( s_loc, 's', 0 );
			(*Solvent).push_back( temp ); 
		}

	}


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

void Reptation (std::vector<Polymer>* Polymers, std::vector <Particle>* Solvent, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = 1; // distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        BackwardReptation (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight); 
        return; 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        ForwardReptation (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight); 
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


void ChainRegrowth (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int index_of_polymer, int x, int y, int z, bool* IMP_BOOL, double* rweight){

	int deg_of_poly = (*Polymers)[index_of_polymer].deg_poly; 
	int index_monomer = rng_uniform (1, deg_of_poly-2); // (1, deg_of_poly-2); 

	// std::cout << "Index of monomer is " << index_monomer << std::endl;

	// decide which end of the polymer do i want to move around 
    bool first_entry_bool = true; 

	int back_or_front = rng_uniform(0, 1); 

	std::vector <std::array<int,3>> old_cut;
	std::vector <std::array <int,3>> new_cut;	
	old_cut.reserve (index_monomer*6);
	new_cut.reserve (index_monomer*6);
	

	if (back_or_front == 0){
		
		old_cut = extract_positions_tail ( &((*Polymers)[index_of_polymer].chain), index_monomer );

		TailSpin (Polymers, index_of_polymer, index_monomer, x, y, z, IMP_BOOL, &first_entry_bool, rweight);  
		
		for ( int i{0}; i < index_monomer; ++i ){
			new_cut.push_back ( (*Polymers)[index_of_polymer].chain[i].coords );
		}

		std::vector <std::array <int,3>> intersection; 

		for ( int i{0}; i < index_monomer; ++i){
			for ( int j{0}; j <index_monomer; ++j){
				if ( old_cut[j] == new_cut [i] ){
					intersection.push_back(old_cut[j]);
					break;
				}
			}				
		}

		for (std::array <int,3>& a: intersection){
			old_cut.erase ( std::remove (old_cut.begin(), old_cut.end(), a), old_cut.end() );
		}

	}

	else {

		old_cut = extract_positions_head ( &((*Polymers)[index_of_polymer].chain), index_monomer );
		HeadSpin (Polymers, index_of_polymer, index_monomer, deg_of_poly, x, y, z, IMP_BOOL, &first_entry_bool, rweight); 

		for (int i{index_monomer+1}; i<deg_of_poly; ++i ){
			new_cut.push_back ( (*Polymers)[index_of_polymer].chain[i].coords );
		}

		std::vector <std::array <int,3>> intersection; 

		for ( int i{index_monomer+1}; i < deg_of_poly; ++i){
			for ( int j{index_monomer+1}; j <deg_of_poly; ++j){
				if ( old_cut[j - (index_monomer+1) ] == new_cut [i - (index_monomer+1) ] ){
					intersection.push_back(old_cut[ j - (index_monomer+1) ]);
					break;
				}
			}				
		}

		for (std::array <int,3>& a: intersection){
			old_cut.erase ( std::remove (old_cut.begin(), old_cut.end(), a), old_cut.end() );
		}


	}

	if (*IMP_BOOL) {
		// initialize some data stores
		std::array <std::array <int,3>, 6> ne_before;
		std::array <std::array <int,3>, 6> ne_after;
		// find neighbors from the old tail to delete 
		std::vector <std::array <int,3>> solvent_to_delete; 
		std::vector <std::array <int,3>> solvent_to_add; 
		solvent_to_delete.reserve (old_cut.size()*6); solvent_to_add.reserve (old_cut.size()*6); 

		if (back_or_front == 0){
			// ################### begin tail corrections... ############################

			// start populating solvent_to_delete with solvent neighbors from old_p
			for ( std::array <int,3>& loc: old_cut ){
				ne_before = obtain_ne_list (loc, x, y, z);

				// make sure there is no monomer in ne_list...
				for ( std::array <int,3>& ne: ne_before ){

					if ( MonomerReporter  (Polymers, &ne) ){
						// if ne is a location of a monomer, continue
						continue;
					}
					else {
						// if it is a neighbor of another monomer, this will be kept in solvent
						if ( MonomerNeighborReporter (Polymers, &ne, x, y, z) ) {
							continue;
						}
						else {
							// if it is a solvent particle but not neighboring any solvent in the new structure, it must NOT BE PRESENT in *Solvent.
							// that is the purpose of solvent_ne!!! Figuring out all the elements that must be cut out from *Solvent. 
							solvent_to_delete.push_back ( ne ); 
						}
					}
				}
			}

			// add on molecules from new_cut to solvent_to_delete
			for ( std::array <int,3>& to_del: new_cut ){
				solvent_to_delete.push_back( to_del );
			}

			// tack on molecules to Solvent that must be in Solvent so that they get appropriately eliminated or not 
			for ( std::array<int,3>& loc: old_cut ) {
				Particle temp ( loc, 's', 0 );
				(*Solvent).push_back ( temp );
			}

			for ( int i{0}; i < index_monomer; ++i){
				ne_after = obtain_ne_list ( (*Polymers)[index_of_polymer].chain[i].coords, x, y, z);
				for ( std::array <int,3>& ne: ne_after ){

					// make sure there is no monomer at the location 
					if ( MonomerReporter (Polymers, &ne) ){
						continue;
					}
					else {
						// if it doesnt have a monomer, it by definition is besides the monomer of interest
						solvent_to_add.push_back (ne);
					}
				}
			}

			// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
			// deleting elements that MUST NOT be in *Solvent... 

			int s_idx = 0; 
			for ( const std::array<int,3>& s_loc: solvent_to_delete ) {
				s_idx = 0; 
				for ( Particle& p: *Solvent ){
			
					if ( p.coords == s_loc ){
						(*Solvent).erase ( (*Solvent).begin() + s_idx ); 
						break; 
					}
					s_idx += 1; 
				}
			}

			// (3) now that *Solvent has been cleaned out, time to add new things to it 
			bool present = false; 
			for ( const std::array <int,3>& s_loc: solvent_to_add ){
				present = false; 
				for ( Particle& p: *Solvent ){
					// if it is already present in *Solvent, breakout and continue to next location
					if ( p.coords == s_loc ){
						present = true;
						break; 
					}
				}
				if ( present ){
					continue; 
				}
				else {
					Particle temp ( s_loc, 's', 0 );
					(*Solvent).push_back( temp ); 
				}
			}

			// ################### end tail corrections... ############################

		}

		else {
			// ################### begin head corrections... ############################

			// start populating solvent_to_delete with solvent neighbors from old_p
			for ( std::array <int,3>& loc: old_cut ){
				ne_before = obtain_ne_list (loc, x, y, z);

				// make sure there is no monomer in ne_list... 
				for ( std::array <int,3>& ne: ne_before ){

					if ( MonomerReporter ( Polymers, &ne) ){
						// if ne is a location of a monomer, continue 
						continue; 
					}
					else {
						// if it is a neighbor of another monomer, this will be kept in solvent 
						if ( MonomerNeighborReporter (Polymers, &ne, x, y, z) ) {
							continue;
						}
						else {
							// if it is a solvent particle but not neighboring any solvent in the new structure, it must NOT BE PRESENT in *Solvent.
							// that is the purpose of solvent_ne!!! Figuring out all the elements that must be cut out from *Solvent. 
							solvent_to_delete.push_back ( ne );
						}
					}
				}
			}

			// add on molecules from new_cut to solvent_to_delete
			for ( std::array <int,3>& to_del: new_cut ){
				solvent_to_delete.push_back( to_del );
			}

			// tack on molecules to Solvent that must be in Solvent so that they get appropriately eliminated or not 
			for ( std::array<int,3>& loc: old_cut ) {
				Particle temp ( loc, 's', 0 );
				(*Solvent).push_back ( temp );
			}			

			for ( int i{index_monomer+1}; i<deg_of_poly; ++i ){
				ne_after = obtain_ne_list ( (*Polymers)[index_of_polymer].chain[i].coords, x, y, z );
				for ( std::array <int,3>& ne: ne_after ){

					// make sure there is no monomer at the location 
					if ( MonomerReporter (Polymers, &ne) ) {
						continue;
					}
					else {
						// if it doesnt have a monomer, it by definition is besides the monomer of interest
						solvent_to_add.push_back (ne); 
					}
				}
			}

			// (2) now that i know which solvent molecules need to be eliminated from *Solvent, and which need to be added, let's do the deed. 
			// deleting elements that MUST NOT be in *Solvent... 

			int s_idx = 0; 
			for ( const std::array<int,3>& s_loc: solvent_to_delete ){
				s_idx = 0;
				for ( Particle& p: *Solvent ){
					if (p.coords == s_loc){
						(*Solvent).erase ( (*Solvent).begin() + s_idx );
						break;
					}
					s_idx += 1;
				}
			}

			// (3) now that *Solvent has been cleaned out, time to add new things to it 
			bool present = false; 
			for ( const std::array <int,3>& s_loc: solvent_to_add ){
				present = false; 
				for ( Particle& p: *Solvent ){
					// if it is already present in *Solvent, breakout and continue to next location
					if ( p.coords == s_loc ){
						present = true;
						break; 
					}
				}
				if ( present ){
					continue; 
				}
				else {
					Particle temp ( s_loc, 's', 0 );
					(*Solvent).push_back( temp ); 
				}
			}
			// ################### end head corrections... ############################
		}

	}

	return; 
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector <std::array <int,3>> extract_positions_tail (std::vector <Particle>* chain, int pivot_idx){

	std::vector <std::array <int,3>> extracted;
	extracted.reserve(pivot_idx); 
	for (int i{0}; i<pivot_idx; ++i){
		extracted.push_back( (*chain)[i].coords );
	}
	return extracted; 
}


std::vector <std::array <int,3>> extract_positions_head (std::vector <Particle>* chain, int pivot_idx){

	int deg_of_poly = static_cast<int>( (*chain).size() ); 
	std::vector <std::array <int,3>> extracted; 
	extracted.reserve(deg_of_poly-pivot_idx-1); 
	for (int i{pivot_idx+1}; i<deg_of_poly; ++i){
		extracted.push_back ( (*chain)[i].coords ); 
	}
	return extracted;
}

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of ChainRegrowth
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

void TailSpin (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight){

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
        
		std::array <int, 3> to_check = add_arrays( &( (*Polymers) [index_of_polymer].chain[index_of_monomer].coords), &d);

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
		(*Polymers)[index_of_polymer].chain[index_of_monomer-1].coords = ind_v[ rng_uniform(0, rw_tmp-1) ]; 
		(*rweight) = (*rweight) * rw_tmp/6; 
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

void HeadSpin (std::vector <Polymer>* Polymers, int index_of_polymer, int index_of_monomer, int deg_poly,int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight){
    
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

	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	// std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
    
    // std::cout << "Current pivot point is "; print(  (*PVec) [index_of_polymer].chain[index_of_monomer].coords );
	int rw_tmp = 0; 
	std::vector <std::array <int,3>> ind_v; 

	for (std::array <int,3>& d: adrns){
        // std::cout << "d is "; print(d); 
		std::array <int, 3> to_check = add_arrays ( &( (*Polymers) [index_of_polymer].chain[index_of_monomer].coords ), &d);

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
		
		(*Polymers)[index_of_polymer].chain[index_of_monomer+1].coords = ind_v[ rng_uniform(0, rw_tmp-1) ]; 
		(*rweight) = (*rweight) * rw_tmp/6; 
		HeadSpin (Polymers, index_of_polymer, index_of_monomer+1, deg_poly, x, y, z, IMP_BOOL, first_entry_bool, rweight);

		}

	return ;

}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of TailSpin
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



bool checkOccupancyTail(std::array <int,3>* loc, std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer){
	int Np = static_cast<int>( (*PVec).size() );  

	for (int pnum = 0; pnum < Np; pnum++){

		if (pnum == index_of_polymer){

			int dpol = static_cast<int> ( (*PVec)[pnum].chain.size() ); 
			for (int p = index_of_monomer; p < dpol; p++){

				if ( *loc == (*PVec)[pnum].chain[p].coords ){
					return true;
				}
				else {
					continue;
				}

			}

		}

		else {

			for (const Particle& p: (*PVec)[pnum].chain){

				if ( *loc == p.coords ){
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



bool checkOccupancyHead(std::array <int,3>* loc, std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer){
	int Np = static_cast<int>( (*PVec).size() );  

	for (int pnum = 0; pnum < Np; pnum++){

		if (pnum == index_of_polymer){

			for (int p = 0; p < index_of_monomer+1; p++){

				if ( *loc == (*PVec)[pnum].chain[p].coords ){
					return true;
				}
				else {
					continue;
				}

			}

		}

		else {

			for (const Particle& p: (*PVec)[pnum].chain){

				if ( *loc == p.coords ){
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

void SolventFlip ( std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, double* rweight){

    // number of surrounding solvent molecules 
    int Nsolv = static_cast<int>( (*Solvent).size() ); 
	int Nmer  = static_cast<int>( (*Polymers)[0].chain.size() ); 

    int to_flip = rng_uniform(1, Nsolv ); 
    
    *rweight = 1; 

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	
	std::vector <int> rvec (to_flip);
	std::iota ( std::begin(rvec), std::end(rvec), 0);

	std::shuffle ( rvec.begin(), rvec.end(), std::default_random_engine(seed) );

	int j = 0;
	for (int i: rvec ){
		(*Solvent)[i].orientation = rng_uniform (0, 5);  
		*rweight = (*rweight) * (Nsolv-j)/(Nmer+Nsolv-j); 
		j = j + 1; 
	}
	
	return; 

}


void SolventFlipSingular ( std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, double* rweight){

	// number of surrounding solvent molecules 
    int Nsolv = static_cast<int>( (*Solvent).size() ); 
	int Nmer  = static_cast<int>( (*Polymers)[0].chain.size() ); 

    *rweight = static_cast<double>(Nsolv)/static_cast<double>(Nmer+Nsolv); 

	int ridx = rng_uniform (0, Nsolv-1);

	(*Solvent)[ridx].orientation = rng_uniform (0, 5);  

	return; 
}


///////////////////////////////////////////////////////////////////////////
void PolymerFlip ( std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, double* rweight){
    
    int Nsolv = static_cast<int> ( (*Solvent).size() ); 
    int Nmer  = static_cast<int> ( (*Polymers)[0].chain.size()  ) ; 
    
    int to_flip = rng_uniform (1, Nmer );

    *rweight = 1; 

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	
	std::vector <int> rvec (to_flip);
	std::iota ( std::begin(rvec), std::end(rvec), 0);

	std::shuffle ( rvec.begin(), rvec.end(), std::default_random_engine(seed) );

	int j = 0;
	for (int i: rvec) {
		(*Polymers)[0].chain[i].orientation = rng_uniform (0, 5); 
		*rweight = (*rweight) * static_cast<double>(Nmer-j)/static_cast<double>(Nmer+Nsolv-j);
		j = j+1; 
	}
    
    return; 
}

void PolymerFlipSingular ( std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, double* rweight){
    
    // obtain the list of solvent particles neighboring the polymer 
    int Nsolv    = static_cast <int> ((*Solvent).size()  ); 
    int Nmer     = static_cast <int> ((*Polymers)[0].chain.size() );
    
    // index of monomer unit 
    int idx = rng_uniform   ( 0, static_cast<int> ((*Polymers)[0].chain.size()-1) ) ;

    (*Polymers)[0].chain[idx].orientation = rng_uniform(0,5); 
    *rweight = static_cast<double>(Nmer)/static_cast<double>(Nmer+Nsolv); 
    
    return; 
}

//////////////////////////////////////////////////////////////

void PerturbSystem (std::vector <Polymer>* Polymers_c, std::vector <Particle>* Solvent_c, int x, int y, int z, bool v, bool* IMP_BOOL, double* rweight, std::array <int,9>* attempts, int* move_number){

    int index = rng_uniform(0, static_cast<int>((*Polymers_c).size())-1); 
    int r = rng_uniform(1, 9);
    switch (r) {
        case (1):
            if (v){
               printf("Performing end rotations.\n"); 
            }
            EndRotation		(Polymers_c, Solvent_c, index, x, y, z, IMP_BOOL, rweight);
            *move_number = 1;
            (*attempts)[0] += 1;
            break;  

        case (2):
            if (v){
               printf("Performing kink jump.\n"); 
            }
            KinkJump		(Polymers_c, Solvent_c, index, x, y, z, IMP_BOOL, rweight);
            *move_number = 2; 
            (*attempts)[1] += 1;
            break;   
        
        case (3):
            if (v){
               printf("Performing crank shaft.\n"); 
            }
            CrankShaft		(Polymers_c, Solvent_c, index, x, y, z, IMP_BOOL, rweight);
            *move_number = 3; 
            (*attempts)[2] += 1;
            break; 
        
        case (4):
            if (v){
               printf("Performing reptation.\n"); 
            }
            Reptation 		(Polymers_c, Solvent_c, index, x, y, z, IMP_BOOL, rweight); 
            *move_number = 4; 
            (*attempts)[3] += 1;
            break; 
        
        case (5):
        	if (v) {
        		printf("Performing configuration sampling. \n"); 
        		std::cout << "index of polymer is " << index << std::endl;
        	} 
        	ChainRegrowth 	(Polymers_c, Solvent_c, index, x, y, z, IMP_BOOL, rweight ); 
            *move_number = 5; 
            (*attempts)[4] += 1;
        	break;
        
        case (6): 
        	if (v){
        		printf("Performing solvent orientation flips. \n");
        	}
        	SolventFlip (Polymers_c, Solvent_c, rweight); 
            *move_number = 6; 
            (*attempts)[5] += 1;
        	break; 
        
        case (7):
            if (v) {
                printf("Performing single solvent orientation flip. \n");
            }
            SolventFlipSingular (Polymers_c, Solvent_c, rweight); 
            *move_number = 7; 
            (*attempts)[6] += 1;
            break;
        
        case (8):
        	if (v){
        		printf("Performing polymer orientation flips. \n");
        	}
        	PolymerFlip ( Polymers_c, Solvent_c, rweight );
            *move_number = 8; 
            (*attempts)[7] += 1;
        	break;
        
        case (9):
            if (v) {
                printf("Performing local polymer orientation flips. \n");
            }
            PolymerFlipSingular ( Polymers_c, Solvent_c, rweight); 
            *move_number = 9; 
            (*attempts)[8] += 1;
            break;
    }
    return;
}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of MoveChooser. 
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



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


std::vector <Polymer> ExtractPolymersFromTraj(std::string trajectory, std::string position, int x, int y, int z){

    int NumberOfPolymers = ExtractNumberOfPolymers(position); 

    std::vector <Polymer> PolymerVector; 
    PolymerVector.reserve(NumberOfPolymers);

    std::vector <std::array <int,3>> locations; 

    std::vector <std::string> contents = ExtractContentFromFile(trajectory); // this extracts every line of the file

    // std::vector <int> step_num_store; 
    std::vector <int> index_store; // need this guy to hold the index of the final set of coordinates. 

    //////////////////////////////////////////////////////////////////////
    // extract the final coordinates from the traj file
    std::regex stepnum ("Dumping coordinates at step"); 
    // int step_num = 0; 
    int j {0}; 
    for (std::string& s: contents){

        std::stringstream ss(s);
        std::string temp; 
        int found; 
        // std::stringstream int_ss; 

        if (std::regex_search(s, stepnum)){
            while(!ss.eof()){
                ss >> temp;
                if (std::stringstream(temp) >> found){
                    // step_num_store.push_back(found);
                    index_store.push_back(j);
                } 
            }
        } 
        ++j; 
    }

    contents.erase(contents.begin(), contents.begin() + (index_store[index_store.size()-1] ) ); 

    //////////////////////////////////////////////////////////////////

    
    bool start_bool {false}, end_bool {false}; 
    std::regex start ("START"), end ("END"), reg_poly ("Dumping coordinates of Polymer"); 
    
    int startCount{0}, endCount{0}; 
    
    for (std::string& s: contents){
         
        std::stringstream ss(s); 
        if (std::regex_search(s, start) ){
            ++startCount;
            start_bool = true; 
            end_bool = false; 
            continue; 
        }

        else if (std::regex_search(s, end) ) {
            ++endCount;
            start_bool = false; 
            end_bool = false; 

            Polymer pmer = makePolymer(locations);
            PolymerVector.push_back(pmer);
            
            locations.clear();
            continue; 
            
        }

        else if (start_bool == end_bool){
            continue;
        }

        else{
            std::array <int,3> loc;
            std::string strr; // temp string 
            int k;            // temp int container 
            
            int j{0};
            while (!ss.eof() ){
                ss >> strr; 
                if (std::stringstream(strr) >> k){
                    loc[j] = k;
                    ++j;
                }
            }  

            if (!checkValidityOfCoords(loc, x, y, z)){
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


std::vector <Particle> CreateSolventVector(int x, int y, int z, std::vector <Polymer>* Polymers){

	std::vector <Particle> Solvent; 
	std::vector <std::array<int,3>> p_locations; 

	for (const Polymer& pmer: *Polymers){
		for (const Particle& p: pmer.chain){

			std::array <std::array <int,3>, 6> ne_list = obtain_ne_list (p.coords, x, y, z); 

			for (std::array <int,3>& ne: ne_list){
				// make sure ne is not occupied by a monomer 
				if ( MonomerReporter (Polymers, &ne) ){
					continue;
				}
				// make sure ne is not already occupied
				else if (std::find(p_locations.begin(), p_locations.end(), ne ) != p_locations.end() ){
					continue; 
				}
				else {
					Particle p = Particle (ne, 's', 0);
					p_locations.push_back (p.coords);
					Solvent.push_back(p);
				}
			}
		}
	}
	return Solvent; 

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



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of Translation
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


