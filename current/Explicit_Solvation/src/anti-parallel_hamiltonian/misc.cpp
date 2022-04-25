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
#include <unordered_set>

/* ==================================================
These are some objects I have defined which I tend to use often. 
Helpful definitions which are employed often in the context of the z=6 lattice I am using.
====================================================*/ 

const std::vector <int> ex{1,0,0}, nex{-1,0,0}, ey{0,1,0}, ney{0,-1,0}, ez{0,0,1}, nez{0,0,-1}; 	// unit directions 
const std::vector <std::vector <int>> drns = {ex, nex, ey, ney, ez, nez};  							// vector of unit directions 

std::array <int,3> ax{1,0,0}, nax{-1,0,0}, ay{0,1,0}, nay{0,-1,0}, az{0,0,1}, naz{0,0,-1};		// unit directions
std::array <std::array <int,3> ,6> adrns = {ax, nax, ay, nay, az, naz}; 

std::array <std::string, 6> orientations = {"+x", "-x", "+y", "-y", "+z", "-z"}; 

std::map <std::array <int,3>, std::string> BondMap { {ax, "+x"}, {nax, "-x"}, {ay, "+y"}, {nay, "-y"}, {az, "+z"}, {naz, "-z"} }; 
std::map <std::string, std::string> Dir2InvMap { {"+x", "-x"}, {"-x", "+x"}, {"+y", "-y"}, {"-y", "+y"}, {"+z", "-z"}, {"-z", "+z"}}; 

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

std::array <int,3> subtract_arrays(const std::array <int,3>* a1, const std::array <int,3>* a2){
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



std::array <double,10> ExtractTopologyFromFile(std::string filename){
    
    std::array <double, 10> info_vec; 
    double info; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(filename); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), Emm_a ("Emm_a"), Emm_n ("Emm_n"), Ems_a ("Ems_a"), Ems_n ("Ems_n"), Eb_f ("Eb_f"), Eb_u ("Eb_u"), eof ("END OF FILE"); 
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

    	else if (std::regex_search (s, Eb_f)) {

    		double info = NumberExtractor (s); 
    		info_vec[8] = info; 
    		continue;
    	}

    	else if (std::regex_search (s, Eb_u)) {

    		double info = NumberExtractor (s); 
    		info_vec[9] = info; 
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
	std::vector <std::string> pmer_spins; 
    int size_ = locations.size(); 
    for (int i=0; i<static_cast<int>(size_); i++){
        unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
        std::mt19937 generator(seed); 
        std::uniform_int_distribution<int> distribution (0,5); 
        pmer_spins.push_back( orientations[distribution(generator)] );
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

/*
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
*/

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
	std::string efile, std::string restart_traj, std::string mfile){

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

    if (positions=="blank" || topology == "blank" || dfile=="blank" || efile == "blank" || mfile == "blank" ){
        std::cerr << "positions is " << positions <<", topology is " << topology <<", coordinate dump file is " << dfile << ", energy dump file is " << efile << ", orientation file is " << mfile << "." << std::endl;
        std::cerr << "ERROR: No value for option p (positions file) and/or for option t (energy and geometry file) and/or for option o (name of output dump file) and/or for option e (name of orientation file)" <<
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
        std::ofstream or_file (mfile); 
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

Particle ParticleReporter (std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVect, std::array <int,3> to_check){

	for (const Polymer& pmer: (*PolymerVector)) {
		for (const Particle& p: pmer.chain){ 

			if (to_check == p.coords){
				return p;
			}

		}
	}

	for (const Particle& p: *SolvVect){

		if (to_check == p.coords){
			return p;
		}

	} 
    // (*PolymerVector)[0].printChainCoords();     
    // (*PolymerVector)[1].printChainCoords(); 
    // std::cout << "Coords of solvent: \n";    

    for (const Particle& p: *SolvVect){
        print(p.coords); 
    }
    
	std::cout << "Something is profoundly fucked - reporting from ParticleReporter." << std::endl;
    exit(EXIT_FAILURE);
	// Particle p;
	return (*SolvVect)[0]; 

}

// used only inside OrientationFlip local
/*int SolventIndexReporter ( std::vector <Particle>* Solvent, std::array <int,3>* to_check){

	int idx {0}; 
	for (const Particle& p: *Solvent){
		idx += 1; 
		if (p.coords == (*to_check)){
			return idx;
		}
	}

	std::cout <<"Something is profoundly fucked (SolventIndexReporter)." << std::endl;
	return -1; 

}*/

/*int IsSolvent( std::vector <Polymer>* Polymers, std::array <int,3>* to_check ){

	for (const Polymer& pmer: (*Polymers) ){
		for (const Particle& p: pmer.chain ){

			if (p.coords == (*to_check)){
				return 0;
			}

		}
	}

	return 1; 

}*/

//============================================================
//============================================================
// 
// !!!!!!!!! CRITICAL OBJECT TO RUN SIMULATION ACCURATELY !!!!!!!!!!!!!!!!!!!
// 
// NAME OF FUNCTION: MonomerCheck
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

bool MonomerReporter (std::vector <Polymer>* PolymerVector, std::array <int,3>* to_check){


	for (Polymer& pmer: (*PolymerVector)) {
		for (Particle& p: pmer.chain){ 

			if ((*to_check) == p.coords){
				return true;
			}
		}
	}

	return false; 

}

bool MonomerReporter (std::vector <Polymer>* PolymerVector, std::array <int,3>* check_1, std::array <int,3>* check_2){

	for (const Polymer& pmer: (*PolymerVector)) {
      		for (const Particle& p: pmer.chain){

      			if (p.coords == (*check_1) || p.coords == (*check_2) ){	 
      				return true; 
      			}
      		}
      	}
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

double CalculateEnergy(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int x, int y, int z, double Emm_a, double Emm_n, double Ems_a, double Ems_n, double Eb_f, double Eb_u, double* m_neighbor, int* a_contacts, int* n_contacts){
    double Energy {0.0};
    (*m_neighbor) = 0; 
    (*a_contacts) = 0; 
    (*n_contacts) = 0;  
    // polymer-polymer interaction energies 
    // (*PolymerVector)[0].printChainCoords();
    int bonded = 0; 
    for (const Polymer& pmer: (*PolymerVector)) {
        for (const Particle& p: pmer.chain){
        	std::cout << "--------\nParticle of interest: "; print(p.coords); 
            std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(p.coords, x, y, z); // get neighbor list 
            
            // identify who is bonded to that atom 
            std::vector <Particle> bonded_particles = pmer.ConnectivityMap.at(p); 

            for (const std::array <int, 3>& loc: ne_list){
                
                // if particle is bonded, check for directionality 
                for ( const Particle& p_: bonded_particles ){
                	bonded = 0; 
                	if ( loc == p_.coords ) {
                		bonded = 1; 
                		std::array <int,3> bond = subtract_arrays ( &(p_.coords), &(p.coords) ); 
                		modified_direction (&bond, x, y, z); 
                		// std::cout << "printing out bond..."; 
                		// print(bond); 
                		// bonded interaction between p and p_ 
                		std::string bond_dirn = BondMap [ bond ] ;
                		// std::cout << "bond dirn is: " << bond_dirn << std::endl;

                		if (bond_dirn == p.orientation){
                			// if bond direction is anti-parallel to p_ orientation 
                			std::cout << "Unfavorable bonded config. + "<< Eb_u << std::endl;
                			print( p_.coords ); 
                			Energy += Eb_u; 
                		}
                		else {
                			print( p_.coords ); 
                			std::cout << "Favorable bonded config. + " << Eb_f << std::endl;
                			Energy += Eb_f; 
                		}
                		
                		break;
                	}
                	else {
                		continue; 
                	}

                }
            	
                if (bonded){
                	continue; 
                }

                else {
                	// THESE ARE NONBONDED INTERACTIONS
                	// PTR IS THE NONBONDED PARTICLE NEXT TO P 
            		Particle ptr {ParticleReporter (PolymerVector, SolvVector, loc) };
            		if ((ptr).ptype == "monomer"){
                    	(*m_neighbor) += 0.5; 
                    	
                    	std::cout << "Monomer of interest: "; print(ptr.coords); 
                    	std::array <int,3> bond = subtract_arrays ( & (ptr.coords), &(p.coords));
                    	modified_direction (&bond, x, y, z); 

                    	// if ANTIPARALLEL with respect to how they are connected 
                    	// if orientation of nonbonded monomer is antiparallel to monomer of interest
                    	
            			if ( ((ptr).orientation == Dir2InvMap[p.orientation]) && (BondMap[bond] == p.orientation) ){

            				std::cout << "Yes, antiparallel. + "<<Emm_a << "." << std::endl;
            				Energy += 0.5*Emm_a; 

            			}
            			// if not antiparallel along direction of bond... 
            			else {
            				std::cout << "Not antiparallel. + "<<Emm_n << "." << std::endl;
            				Energy += 0.5*Emm_n; 
            			}

            		}

            		else { 
            			if ( (ptr).orientation == p.orientation ){
                        	(*a_contacts) += 1;
            				Energy += Ems_a;
            			}
            			else {
                        	(*n_contacts) += 1; 
            				Energy += Ems_n; 
            			}

            		}
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

	std::ofstream dump_file(filename, std::ios::out); 
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


void dumpEnergy (double sysEnergy, int step, double m_contacts, int a_contacts, int n_contacts, std::string filename){
    std::ofstream dump_file(filename, std::ios::app); 
    // std::ostringstream os; 
    
    dump_file << sysEnergy << " | " << m_contacts << " | " << a_contacts << " | " << n_contacts << " | " << step << "\n";
    
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
void dumpOrientation( std::vector <Polymer>* PV, std::vector <Particle>* SV, int step, std::string filename, int x, int y, int z ) {
    // std::cout<< "just inside dumpO..."<<std::endl; 
    std::ofstream dump_file (filename, std::ios::app); 
    dump_file << "START for Step " << step << ".\n";
    Particle t_particle; 
    for ( const Polymer& pmer: (*PV) ) {
        for ( const Particle& p: pmer.chain ) {
            // std::cout << "particle is at "; print(p.coords);  
            dump_file << p.orientation << " | ";
            std::array <std::array<int,3> ,6> ne_list = obtain_ne_list (p.coords, x, y, z) ;
            // std::cout << "neighbor list: " << std::endl;
            // print(ne_list); 
            for ( const std::array <int,3>& ne: ne_list) {
                // std::cout <<"Particle positions: "; 
                // print (ne) ;
                t_particle = ParticleReporter (PV, SV, ne);
                // std::cout << "Reported~\n"; 
                if (t_particle.ptype == "solvent"){
                    dump_file << t_particle.orientation << " | ";  
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

std::vector <Polymer> TailRotation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){


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

			NewPol[index].chain[0].coords = to_rot; 

			for (Particle& p: (*SolvVector)){
				if (p.coords == to_rot){
					p.coords = loc_0; 
					break;
				}

			}

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


std::vector <Polymer> TailRotation_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){


    // get the neighborlist of particle at index 1 
    std::array <int,3> loc_0 = (*PolymerVector)[index].chain[0].coords; 
    std::array <int,3> loc_1 = (*PolymerVector)[index].chain[1].coords;
    std::vector <Polymer> NewPol {*PolymerVector};
    
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

	(*rweight) = 0; 

	// find all spots that are available 
	std::vector <std::array <int,3>> idx_v; 
	for (std::array <int,3>& to_rot: ne_list){
 		
		if ( to_rot == loc_1){
			continue; 
		}

		else if ( to_rot == (*PolymerVector)[index].chain[2].coords ) {
			continue; 
		}

		else if ( ! (MonomerReporter(PolymerVector, &to_rot) ) ){
			// std::cout << "a free location is "; print(to_rot); 
			(*rweight) += 1;
			idx_v.push_back( to_rot );  
		}
		
	}

    if ( (*rweight) == 0 ){
    	*IMP_BOOL = false;
    }


	else {
		int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 
		NewPol[index].chain[0].coords = idx_v[ r ];  

		// std::cout << "Suggested rot position is "; print (idx_v [r]); 

		for (Particle& p: (*SolvVector)){
			if (p.coords == idx_v[ r ]){
				p.coords = loc_0; 
				break;
			}
		}

		NewPol[index].ChainToConnectivityMap(); 
		(*rweight) = (*rweight)/6; 
	}
	
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

std::vector <Polymer> HeadRotation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector,int index, int x, int y, int z, bool* IMP_BOOL){


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

			for (Particle& p: (*SolvVector)){
				if (p.coords == to_rot){
					p.coords = loc_0; 
					break;
				}

			}

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


std::vector <Polymer> HeadRotation_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector,int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){


    // get the neighborlist of particle at index 1 
    int dop = (*PolymerVector)[index].deg_poly; 
    std::array <int,3> loc_0 = (*PolymerVector)[index].chain[dop-1].coords; 
    std::array <int,3> loc_1 = (*PolymerVector)[index].chain[dop-2].coords;
    std::vector <Polymer> NewPol {*PolymerVector};
    
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list(loc_1, x, y, z) ; 

    (*rweight) = 0; 
	
    // find all spots that are available 
    std::vector <std::array <int,3>> idx_v; 

    for (std::array <int,3>& to_rot: ne_list){

    	if (to_rot == loc_1){
    		continue; 
    	}

    	else if ( to_rot == (*PolymerVector)[index].chain[dop-3].coords ) {
			continue; 
		}

    	else if ( ! (MonomerReporter(PolymerVector, &to_rot) ) ){
    		// std::cout << "a free location is "; print(to_rot); 
    		(*rweight) += 1; 
    		idx_v.push_back(to_rot); 
    	}

    }

    if ( (*rweight) == 0){
    	*IMP_BOOL = false;
    }

    else {
    	int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 
    	NewPol[index].chain[dop-1].coords = idx_v[ r ]; 
    	// std::cout << "Suggested rot position is "; print (idx_v [r]); 
    	for (Particle& p: (*SolvVector)){
    		if (p.coords == idx_v[ r ]){
    			p.coords = loc_0;
    			break; 
    		}
    	}

    	NewPol[index].ChainToConnectivityMap(); 
    	(*rweight) = (*rweight)/6; 

    }

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


std::vector <Polymer> EndRotation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return TailRotation(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return HeadRotation(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL); 

    }
    
}


std::vector <Polymer> EndRotation_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return TailRotation_Rosenbluth(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL, rweight); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return HeadRotation_Rosenbluth(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL, rweight); 

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


std::vector <Polymer> KinkJump(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){

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
			for (Particle& p: (*SolvVector)){
				if (p.coords == to_check){
					p.coords = (*PolymerVector)[index].chain[idx+1].coords; 
					break;
				}

			}
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


std::vector <Polymer> KinkJump_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    std::vector <Polymer> NewPol = *PolymerVector; 

    std::vector <int> k_idx = NewPol.at(index).findKinks(); 

    (*rweight) = 0; 

    if (k_idx.size() == 0 ){
        *IMP_BOOL = false;
        // std::cout << "No kinks found in polymer..." << std::endl;
        return *PolymerVector;
    }

    std::vector <int> idx_v;
    std::vector <std::array <int,3>> pos_v;  

    for (int idx: k_idx){

        // std::cout << "idx right before kink spot is " << idx << std::endl; 

        // std::array <int,3> d1 = subtract_arrays(&(NewG.PolymersInGrid[index].chain[idx+1].coords), &(NewG.PolymersInGrid[index].chain[idx].coords) );
        std::array <int,3> d2 = subtract_arrays(&(NewPol[index].chain[idx+2].coords), &(NewPol[index].chain[idx+1].coords) ); 

        std::array <int,3> to_check = add_arrays( &( NewPol[index].chain[idx].coords ), &d2); 
        
        impose_pbc(&to_check, x, y, z); 

        if ( ! ( MonomerReporter (PolymerVector, &to_check) ) ){
        	// std::cout << "kink is at "; print(NewPol[index].chain[idx+1].coords);
        	// std::cout << "a free location is "; print(to_check); 
        	(*rweight) += 1;
        	idx_v.push_back(idx);  
        	pos_v.push_back(to_check); 
        }

    }
    
    if ( (*rweight) == 0){
    	*IMP_BOOL = false; 
    }

	else {
		
		int r = rng_uniform(0, static_cast<int> (idx_v.size() - 1) ); 
		NewPol[index].chain[idx_v[r]+1].coords = pos_v[r]; 

		// std::cout << "kink jump spot is: "; print(pos_v[r]);
		
		for (Particle& p: (*SolvVector)){
			if (p.coords == pos_v[r]){
				p.coords = (*PolymerVector)[index].chain[idx_v[r]+1].coords; 
				break;
			}
		}

		NewPol[index].ChainToConnectivityMap(); 
		(*rweight) = (*rweight)/( k_idx.size() ); 
	}

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


std::vector <Polymer> CrankShaft(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){


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
      		int d = 0; 

      		for (Particle& p: (*SolvVector)){
				if (p.coords == to_check_1){
					++d; 
					p.coords = (*PolymerVector)[index].chain[idx+1].coords; 
					if (d==2){
						break;
					}
				}
				else if (p.coords == to_check_2){
					++d; 
					p.coords = (*PolymerVector)[index].chain[idx+2].coords; 
					if (d==2){
						break;
					}
				}
			}

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




std::vector <Polymer> CrankShaft_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){


    std::vector <Polymer> NewPol {*PolymerVector}; 
    std::vector <int> c_idx = NewPol[index].findCranks(); 

    if ( c_idx.size()==0 ){
        *IMP_BOOL = false;  
        // std::cout << "No cranks..." << std::endl;
        return NewPol; 
    }
    

    // std::shuffle (std::begin (c_idx), std::end(c_idx), std::default_random_engine() ); 
    int idx = c_idx[ rng_uniform (0, static_cast<int> (c_idx.size() - 1 ) ) ]; 
	std::vector <std::array <int,3>> pos_v1;
	std::vector <std::array <int,3>> pos_v2;

	(*rweight) = 0; 

	std::array <int,3> HingeToKink = subtract_arrays(&(NewPol[index].chain[idx+2].coords), &(NewPol[index].chain[idx+3].coords) ); 
    std::array <int,3> HingeToHinge = subtract_arrays(&(NewPol[index].chain[idx+3].coords ), &(NewPol[index].chain[idx].coords ) );
    std::array <std::array <int,3>,3> drns = HingeSwingDirections (& (HingeToHinge), &(HingeToKink), x, y, z); 

	// std::array <int,3> d1 = drns[choice];

	for (std::array <int,3>& d: drns ){

        std::array <int,3> to_check_1 = add_arrays ( &(NewPol[index].chain[idx].coords), &d ); 
        std::array <int,3> to_check_2 = add_arrays ( &(NewPol[index].chain[idx+3].coords), &d ); 
        impose_pbc(&to_check_1, x, y, z); 
        impose_pbc(&to_check_2, x, y, z);   

		
      	// check if site is unoccupied 
        if ( ! ( MonomerReporter (PolymerVector, &to_check_1, &to_check_2) ) ) {
        	pos_v1.push_back(to_check_1); 
        	pos_v2.push_back(to_check_2); 
        	(*rweight) += 1; 
        }
    }

    // if the site is unoccupied for sure 

    if ( (*rweight) == 0 ){
    	*IMP_BOOL = false; 
    	// std::cout << "No allowable crankshafts." << std::endl;
    }

  	else {

  		int r = rng_uniform( 0, pos_v1.size() - 1 ); 

  		// std::cout << "cranked pos 1 is "; print(pos_v1[r]); 
  		// std::cout << "cranked pos 2 is "; print(pos_v2[r]); 

  		NewPol[index].chain[idx+1].coords = pos_v1[r];
  		NewPol[index].chain[idx+2].coords = pos_v2[r];
  		int d = 0; 

  		for (Particle& p: (*SolvVector)){
			if (p.coords == pos_v1[r]){
				++d; 
				p.coords = (*PolymerVector)[index].chain[idx+1].coords; 
				if (d==2){
					break;
				}
			}
			else if (p.coords == pos_v2[r]){
				++d; 
				p.coords = (*PolymerVector)[index].chain[idx+2].coords; 
				if (d==2){
					break;
				}
			}
		}
		NewPol[index].ChainToConnectivityMap();
		(*rweight) = (*rweight)/ 3 ; 
  	}

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

std::vector <Polymer> ForwardReptation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){

    std::vector <Polymer> NewPol {*PolymerVector}; 
    int deg_poly = (*PolymerVector)[index].deg_poly; 

    // get size of polymer chain 
    // int size = NewPol[index].chain.size(); 
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( NewPol[index].chain[deg_poly-1].coords, x, y, z ); 

    // get rid of second to last monomer position from ne list 
    // ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), NewG.PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() ); 

    bool b {false}; 
    bool b_ind1{false};
	size_t tries = 0; 
    for (std::array <int,3>& to_check: ne_list){

    	b = false; 

    	//for (const Polymer& pmer: (*PolymerVector)){
    	for (int pmer_index=0; pmer_index < static_cast<int>((*PolymerVector).size()); pmer_index++) {

    		// for (const Particle& p: pmer.chain){
    		for (int particle_index=0; particle_index < static_cast<int>((*PolymerVector)[pmer_index].chain.size()); particle_index++ ){ 

    			// if there is a match, do the following... 
    			if ((*PolymerVector)[pmer_index].chain[particle_index].coords == to_check){

    				// check if p is the first element of the polymer of the right index! 
    				if (pmer_index == index){

    					if (particle_index == 0){
    						b = false;
    						b_ind1 = true; 
    						break;
    					}
    					else {
    						b = true; 
    						break; 
    					}
    				}
    				else {
    					b = true;
    					break;
    				}
    			}

    		}

    		if (b) {
    			break;
    		}

    		if (b_ind1){
    			break;
    		}

    	}

    	if (!b){
    		// if everything checks out, do the deed - make it slither forward 
    		// std::cout << "In acceptance land." << std::endl; 
    		for (int i{0}; i<deg_poly; i++){

    			if ( i != deg_poly-1 ){
    				NewPol[index].chain[i].coords = (*PolymerVector)[index].chain[i+1].coords ; 
    			}
    			else {
    				NewPol[index].chain[i].coords = to_check; 
    			}
    		}

    		if (b_ind1){
    			break;
    		}
    		// std::cout << "number of tries is " << tries << std::endl;
  			// std::cout << "max number of tries for forward reptation can be " << ne_list.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl; 
  			else {
  				for (Particle& p: (*SolvVector)){
					if (p.coords == to_check){
						p.coords = (*PolymerVector)[index].chain[0].coords; 
						break;
					}
				}
				break;	
  			}
    		

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


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


std::vector <Polymer> ForwardReptation_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    std::vector <Polymer> NewPol {*PolymerVector}; 
    int deg_poly = (*PolymerVector)[index].deg_poly; 

    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( NewPol[index].chain[deg_poly-1].coords, x, y, z ); 
    std::vector <std::array<int,3>> idx_v; 

    (*rweight) = 0; 

    for (std::array <int,3>& to_check: ne_list){

    	if ( to_check == NewPol[index].chain[deg_poly-2].coords ){
    		continue; 
    	}

		else if ( ! (MonomerReporter(PolymerVector, &to_check) )){
			(*rweight) += 1; 
			idx_v.push_back(to_check);
		}
	}

	if ( (*rweight) == 0){
		*IMP_BOOL = false;
	}

	else {
		int r = rng_uniform( 0, idx_v.size()-1 );
		// if everything checks out, do the deed - make it slither forward 
		// std::cout << "In acceptance land." << std::endl; 
		// std::cout << "location of new spot is "; print(idx_v[r]); 
		for (int i{0}; i<deg_poly; ++i){

			if ( i != deg_poly-1 ){
				NewPol[index].chain[i].coords = (*PolymerVector)[index].chain[i+1].coords ; 
			}
			else {
				NewPol[index].chain[i].coords = idx_v[r]; 
			}
		}
		
		for (Particle& p: (*SolvVector)){
			if (p.coords == idx_v[r]){
				p.coords = (*PolymerVector)[index].chain[0].coords; 
				break;
			}
		}	
		
		NewPol[index].ChainToConnectivityMap();
		(*rweight) = (*rweight)/6; 
	}

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

std::vector <Polymer> BackwardReptation(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){

    std::vector <Polymer> NewPol {*PolymerVector}; 
    int deg_poly = (*PolymerVector)[index].deg_poly; 

    // get size of polymer chain 
    // int size = NewPol[index].chain.size(); 
    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( NewPol[index].chain[0].coords, x, y, z ); 

    // get rid of second to last monomer position from ne list 
    // ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), NewG.PolymersInGrid.at(index).chain.at(size-2).coords), ne_list.end() ); 

    bool b { false }; 
    bool b_indf { false }; 
	size_t tries = 0; 
    for (std::array <int,3>& to_check: ne_list){

    	b = false; 

    	for (int pmer_index=0; pmer_index < static_cast<int>((*PolymerVector).size()); ++pmer_index) {

    		// for (const Particle& p: pmer.chain){
    		for (int particle_index=0; particle_index < static_cast<int>((*PolymerVector)[pmer_index].chain.size() ); ++particle_index ){ 

    			// in the event the to_check particle matches with something in the vector
    			// check if that match is the final monomer bead 
    			if ((*PolymerVector)[pmer_index].chain[particle_index].coords == to_check ){ 
    				
    				// std::cout << "hello. there is a hit." << std::endl;
    				// print(to_check); 
    				if (pmer_index == index){

    					if (particle_index == deg_poly-1){
    						b = false; 	
    						b_indf = true;					// false because i do not want the particle to be disregarded
    						break; 
    					}

    					else {
    						b = true;						// true, you have to disregard the to_check position
    						break;
    					}
    				}
					else {
    					b = true;
    					break; 
    				}

    			}
    			

    			// if the current position is occupied by a bead that is not the right bead, move on. 
    		}
    		if (b) {
    			// if the place is occupied and not by a good particle, get out of the loop and choose another to_check
    			break;
    		}

    		if (b_indf){
    			break;
    		}

    	}

    	if (!b){
    		// std::cout << "In acceptance land." << std::endl;
    		// print(to_check);
    		for (int i{0}; i<deg_poly; ++i){

    			if ( i != deg_poly-1 ){
    				NewPol[index].chain[deg_poly-1-i].coords = (*PolymerVector)[index].chain[deg_poly-2-i].coords ; 
    			}
    			else {

    				NewPol[index].chain[deg_poly-1-i].coords = to_check; 
    			}
    		}

    		if (b_indf){
    			break;
    		}

    		else {
  				for (Particle& p: (*SolvVector)){
					if (p.coords == to_check){
						p.coords = (*PolymerVector)[index].chain[deg_poly-1].coords; 
						// std::cout << "Is this break being hit?" << std::endl;
						break;
					}
				}	
				break;
  			}
    		// std::cout << "number of tries is " << tries << std::endl;
  			// std::cout << "max number of tries for back reptation can be " << ne_list.size() << std::endl;
  			// std::cout << "This should be the only tries statement for this particular move." << std::endl; 
    		// break;
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


std::vector <Polymer> BackwardReptation_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    std::vector <Polymer> NewPol {*PolymerVector}; 
    int deg_poly = (*PolymerVector)[index].deg_poly; 

    std::array <std::array <int,3>, 6> ne_list = obtain_ne_list( NewPol[index].chain[0].coords, x, y, z ); 
    std::vector <std::array <int,3>> idx_v; 

    (*rweight) = 0; 

     
    for (std::array <int,3>& to_check: ne_list){

    	if ( to_check == NewPol[index].chain[1].coords){
    		continue; 
    	}

    	else if ( ! (MonomerReporter(PolymerVector, &to_check) ) ){
    		(*rweight) += 1; 
    		idx_v.push_back(to_check); 
    	}
    }

    if ( (*rweight) == 0 ){
    	*IMP_BOOL = false; 
    }

    else {

    	int r = rng_uniform( 0, idx_v.size()-1 );

    	for (int i{0}; i <deg_poly; ++i){
    		if ( i != deg_poly-1 ){
    			NewPol[index].chain[deg_poly-1-i].coords = (*PolymerVector)[index].chain[deg_poly-2-i].coords;
    		}
    		else {
    			NewPol[index].chain[deg_poly-1-i].coords = idx_v[r]; 
    		}
    	}

    	for (Particle& p: (*SolvVector)) {
    		if (p.coords == idx_v[r]){
    			p.coords = (*PolymerVector)[index].chain[deg_poly-1].coords;
    			break;
    		}
    	}

    	NewPol[index].ChainToConnectivityMap(); 
    	(*rweight) = (*rweight)/6; 
    }

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

std::vector <Polymer> Reptation(std::vector<Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return BackwardReptation(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return ForwardReptation(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL); 

    }
    
}


std::vector <Polymer> Reptation_Rosenbluth(std::vector<Polymer>* PolymerVector, std::vector <Particle>* SolvVector, int index, int x, int y, int z, bool* IMP_BOOL, double* rweight){

    unsigned seed = static_cast<unsigned> (std::chrono::system_clock::now().time_since_epoch().count());
    std::mt19937 generator(seed); 
    std::uniform_int_distribution<int> distribution (0,1); 
    int num = distribution(generator); 
    // std::cout << "rng is " << num << std::endl;
    if (num==0){
        // std::cout << "Zero index rotation!" << std::endl;
        return BackwardReptation_Rosenbluth(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL, rweight); 
    }
    else {
        // std::cout << "Final index rotation!" << std::endl;
        return ForwardReptation_Rosenbluth(PolymerVector, SolvVector, index, x, y, z, IMP_BOOL, rweight); 

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

void ChainRegrowth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVec, int index_of_polymer, int x, int y, int z, bool* IMP_BOOL){

	// std::vector <Polymer> copy_pvec {*PolymerVector}; 
	std::vector <std::array<int,3>> old_p;

	int deg_of_poly = (*PolymerVector)[index_of_polymer].deg_poly; 
	int index_monomer = rng_uniform(1, deg_of_poly-2); 

	// std::cout << "Index of monomer is " << index_monomer << std::endl;

	// decide which end of the polymer do i want to move around 
	bool b = false; 
    bool first_entry_bool = true; 

	int back_or_front = rng_uniform(0, 1); 

	if (back_or_front == 0){
		old_p = extract_positions_tail ( &((*PolymerVector)[index_of_polymer].chain), index_monomer );
		// std::cout << "Index of monomer is " << index_monomer << std::endl; 
		// std::cout << "Pivot position is "; 
		// print((*PolymerVector)[index_of_polymer].chain[index_monomer].coords); 
        // std::cout << "Tail spin being performed..." << std::endl;
		TailSpin(PolymerVector, index_of_polymer, index_monomer, x, y, z, &b, IMP_BOOL, &first_entry_bool);  
	}

	else {
		old_p = extract_positions_head ( &((*PolymerVector)[index_of_polymer].chain), index_monomer );
        // std::cout << "Index of monomer is " << index_monomer << std::endl; 
		// std::cout << "Pivot position is "; 
		// print((*PolymerVector)[index_of_polymer].chain[index_monomer].coords); 
        // std::cout << "Head spin being performed..." << std::endl;
		HeadSpin(PolymerVector, index_of_polymer, index_monomer, deg_of_poly, x, y, z, &b, IMP_BOOL, &first_entry_bool); 
	}

	// In the above half of the code, tailspin or headspin will be performed. However, once it has been performed, 
	// you have to make sure the solvent molecules that have been displaced will now have another home 

	if (*IMP_BOOL){

		if (back_or_front == 0 ){
			// check tail end of poly
			// std::cout << "Performed a tail spin..." << std::endl;

			// std::vector <std::array <int,3>> old_p = extract_positions_tail ( &(copy_pvec[index_of_polymer].chain), index_monomer );
			std::vector < std::array <int,3>> new_p = extract_positions_tail ( &((*PolymerVector)[index_of_polymer].chain), index_monomer ) ;  

			// sort both vectors 
			std::sort ( old_p.begin(), old_p.end()); 
			std::sort ( new_p.begin(), new_p.end()); 

			// collect common elements 
			std::vector <std::array <int,3>> store; 
			std::set_intersection( old_p.begin(), old_p.end(), new_p.begin(), new_p.end(), std::back_inserter(store) ); 

			// erase items in each vector that match the items in the set 
			
			for (const std::array <int,3>& a: store){

				for (size_t j{0}; j < old_p.size(); ++j){
					
					if (a == old_p[j]){
						old_p.erase(std::remove(old_p.begin(), old_p.end(), a), old_p.end() ); 
						new_p.erase(std::remove(new_p.begin(), new_p.end(), a), new_p.end() );
						break;
					}
				}
			}

			// make sure that the solvent particles that were originally at the spots NOW occupied by the polymer
			// are set to the spots originally occupied by the polymer 

			for (size_t i{0}; i<old_p.size(); ++i){

				for ( Particle& p: (*SolvVec)){ 

					if (p.coords == new_p[i]){

						p.coords = old_p[i];
						break; 
					}
				}
			}
		}
		else {
			
			// check head end of poly
			// std::vector <std::array <int,3>> old_p = extract_positions_head ( &(copy_pvec[index_of_polymer].chain), index_monomer );
			std::vector <std::array <int,3>> new_p = extract_positions_head ( &((*PolymerVector)[index_of_polymer].chain), index_monomer ) ;  

			// sort both vectors 
			std::sort ( old_p.begin(), old_p.end()); 
			std::sort ( new_p.begin(), new_p.end()); 

			// collect common elements 
			std::vector <std::array <int,3>> store; 
			std::set_intersection( old_p.begin(), old_p.end(), new_p.begin(), new_p.end(), std::back_inserter(store) ); 

			// erase items in each vector that match the items in the set 
			
			for (const std::array <int,3>& a: store){

				for (size_t j{0}; j < old_p.size(); ++j){
					
					if (a == old_p[j]){
						old_p.erase(std::remove(old_p.begin(), old_p.end(), a), old_p.end() ); 
						new_p.erase(std::remove(new_p.begin(), new_p.end(), a), new_p.end() );
						break;
					}
				}
			}

			// make sure that the solvent particles that were originally at the spots NOW occupied by the polymer
			// are set to the spots originally occupied by the polymer 

			for (size_t i{0}; i<old_p.size(); ++i){

				for ( Particle& p: (*SolvVec)){ 

					if (p.coords == new_p[i]){

						p.coords = old_p[i];
						break; 
					}
				}
			}

		}
	}

	return; 

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////


void ChainRegrowth_Rosenbluth(std::vector <Polymer>* PolymerVector, std::vector <Particle>* SolvVec, int index_of_polymer, int x, int y, int z, bool* IMP_BOOL, double* rweight){

	// std::vector <Polymer> copy_pvec {*PolymerVector}; 
	std::vector <std::array<int,3>> old_p;

	int deg_of_poly = (*PolymerVector)[index_of_polymer].deg_poly; 
	int index_monomer = rng_uniform(1, deg_of_poly-2); 

	// std::cout << "Index of monomer is " << index_monomer << std::endl;

	// decide which end of the polymer do i want to move around 
    bool first_entry_bool = true; 

	int back_or_front = 0; //rng_uniform(0, 1); 

	if (back_or_front == 0){
		old_p = extract_positions_tail ( &((*PolymerVector)[index_of_polymer].chain), index_monomer );
		// std::cout << "Index of monomer is " << index_monomer << std::endl; 
		// std::cout << "Pivot position is "; 
		// print((*PolymerVector)[index_of_polymer].chain[index_monomer].coords); 
        // std::cout << "Tail spin being performed..." << std::endl;
		TailSpin_Rosenbluth(PolymerVector, index_of_polymer, index_monomer, x, y, z, IMP_BOOL, &first_entry_bool, rweight);  
	}

	else {
		old_p = extract_positions_head ( &((*PolymerVector)[index_of_polymer].chain), index_monomer );
        // std::cout << "Index of monomer is " << index_monomer << std::endl; 
		// std::cout << "Pivot position is "; 
		// print((*PolymerVector)[index_of_polymer].chain[index_monomer].coords); 
        // std::cout << "Head spin being performed..." << std::endl;
		HeadSpin_Rosenbluth(PolymerVector, index_of_polymer, index_monomer, deg_of_poly, x, y, z, IMP_BOOL, &first_entry_bool, rweight); 
	}

	// In the above half of the code, tailspin or headspin will be performed. However, once it has been performed, 
	// you have to make sure the solvent molecules that have been displaced will now have another home 

	if (*IMP_BOOL){

		if (back_or_front == 0 ){
			// check tail end of poly
			// std::cout << "Performed a tail spin..." << std::endl;

			// std::vector <std::array <int,3>> old_p = extract_positions_tail ( &(copy_pvec[index_of_polymer].chain), index_monomer );
			std::vector < std::array <int,3>> new_p = extract_positions_tail ( &((*PolymerVector)[index_of_polymer].chain), index_monomer ) ;  

			// sort both vectors 
			std::sort ( old_p.begin(), old_p.end()); 
			std::sort ( new_p.begin(), new_p.end()); 

			// collect common elements 
			std::vector <std::array <int,3>> store; 
			std::set_intersection( old_p.begin(), old_p.end(), new_p.begin(), new_p.end(), std::back_inserter(store) ); 

			// erase items in each vector that match the items in the set 
			
			for (const std::array <int,3>& a: store){

				for (size_t j{0}; j < old_p.size(); ++j){
					
					if (a == old_p[j]){
						old_p.erase(std::remove(old_p.begin(), old_p.end(), a), old_p.end() ); 
						new_p.erase(std::remove(new_p.begin(), new_p.end(), a), new_p.end() );
						break;
					}
				}
			}

			// make sure that the solvent particles that were originally at the spots NOW occupied by the polymer
			// are set to the spots originally occupied by the polymer 

			for (size_t i{0}; i<old_p.size(); ++i){

				for ( Particle& p: (*SolvVec)){ 

					if (p.coords == new_p[i]){

						p.coords = old_p[i];
						break; 
					}
				}
			}
		}
		else {
			
			// check head end of poly
			// std::vector <std::array <int,3>> old_p = extract_positions_head ( &(copy_pvec[index_of_polymer].chain), index_monomer );
			std::vector <std::array <int,3>> new_p = extract_positions_head ( &((*PolymerVector)[index_of_polymer].chain), index_monomer ) ;  

			// sort both vectors 
			std::sort ( old_p.begin(), old_p.end()); 
			std::sort ( new_p.begin(), new_p.end()); 

			// collect common elements 
			std::vector <std::array <int,3>> store; 
			std::set_intersection( old_p.begin(), old_p.end(), new_p.begin(), new_p.end(), std::back_inserter(store) ); 

			// erase items in each vector that match the items in the set 
			
			for (const std::array <int,3>& a: store){

				for (size_t j{0}; j < old_p.size(); ++j){
					
					if (a == old_p[j]){
						old_p.erase(std::remove(old_p.begin(), old_p.end(), a), old_p.end() ); 
						new_p.erase(std::remove(new_p.begin(), new_p.end(), a), new_p.end() );
						break;
					}
				}
			}

			// make sure that the solvent particles that were originally at the spots NOW occupied by the polymer
			// are set to the spots originally occupied by the polymer 

			for (size_t i{0}; i<old_p.size(); ++i){

				for ( Particle& p: (*SolvVec)){ 

					if (p.coords == new_p[i]){
						p.coords = old_p[i];
						break; 
					}
				}
			}
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

void TailSpin(std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* b, bool* IMP_BOOL, bool* first_entry_bool){

	// std::cout << "index of monomer is " << index_of_monomer << std::endl;
    
    if (*first_entry_bool){
	    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	    std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
        *first_entry_bool = false; 
    }


	if (index_of_monomer == 0){
		// std::cout << "You have reached the final spot via tail spin!" << std::endl;
		*b = true; 
		*IMP_BOOL = true; 
		(*PVec)[index_of_polymer].ChainToConnectivityMap(); 
		return ; 
	}

    
    // std::cout << "Current pivot point is "; print(  (*PVec) [index_of_polymer].chain[index_of_monomer].coords );
	for (std::array <int,3>& d: adrns){ 
        
        // std::cout << "d is "; print(d); 
		std::array <int, 3> to_check = add_arrays( &( (*PVec) [index_of_polymer].chain[index_of_monomer].coords), &d);

		impose_pbc(&to_check, x, y, z); 
        // std::cout << "Growth location is "; print(to_check); 

		if (checkOccupancyTail(&to_check, PVec, index_of_polymer, index_of_monomer)){
            // std::cout << "This location is occupied! - "; print(to_check); 
            // std::cout << "Moving on to another growth spot..." << std::endl;
			continue; 
		}

		else {

			(*PVec)[index_of_polymer].chain[index_of_monomer-1].coords = to_check; 
            // std::cout << "Moving in deeper...\n\t"; 
			TailSpin (PVec, index_of_polymer, index_of_monomer-1, x, y, z, b, IMP_BOOL, first_entry_bool); 

			if (*b){
				break; 
			}
			else {
                // std::cout << "We went down a rabbithole, but it failed. So moving on..." << std::endl;
				continue; 
			}


		}

	}

	if ( !(*b) ){
		*IMP_BOOL = false; 
	}

	return; 

}



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void TailSpin_Rosenbluth(std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer, int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight){

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
		(*PVec)[index_of_polymer].ChainToConnectivityMap(); 
		return ; 
	}

    
    // std::cout << "Current pivot point is "; print(  (*PVec) [index_of_polymer].chain[index_of_monomer].coords );
    int rw_tmp = 0; 
    std::vector <std::array <int,3>> ind_v; 
	for (std::array <int,3>& d: adrns){ 
        
		std::array <int, 3> to_check = add_arrays( &( (*PVec) [index_of_polymer].chain[index_of_monomer].coords), &d);

		impose_pbc(&to_check, x, y, z); 
        // std::cout << "Growth location is "; print(to_check); 

		if (checkOccupancyTail(&to_check, PVec, index_of_polymer, index_of_monomer)){ 
			continue; 
		}

		else {
			rw_tmp += 1;
			// std::cout <<"a possible position is: "; print(to_check); 
			ind_v.push_back(to_check); 
		}
	}

	if (rw_tmp == 0){
		*IMP_BOOL = false;
	}

	else{
		(*PVec)[index_of_polymer].chain[index_of_monomer-1].coords = ind_v[ rng_uniform(0, rw_tmp-1) ]; 
		(*rweight) = (*rweight) * rw_tmp/6; 
		// std::cout << "rweight is " << *rweight << std::endl;
    	// std::cout << "Moving in deeper...\n\t"; 
		TailSpin_Rosenbluth (PVec, index_of_polymer, index_of_monomer-1, x, y, z, IMP_BOOL, first_entry_bool, rweight); 
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

void HeadSpin(std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer, int deg_poly,int x, int y, int z, bool* b, bool* IMP_BOOL, bool* first_entry_bool){
    
     if (*first_entry_bool){
	    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	    std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
        *first_entry_bool = false; 
    }

	if (index_of_monomer == deg_poly-1){
		// std::cout << "You have reached the final spot of head spin!" << std::endl;
		*b = true; 
		*IMP_BOOL = true;
		(*PVec)[index_of_polymer].ChainToConnectivityMap(); 
		return ;
	}

	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	// std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
    
    // std::cout << "Current pivot point is "; print(  (*PVec) [index_of_polymer].chain[index_of_monomer].coords );

	for (std::array <int,3>& d: adrns){
        // std::cout << "d is "; print(d); 
		std::array <int, 3> to_check = add_arrays ( &( (*PVec) [index_of_polymer].chain[index_of_monomer].coords ), &d);

		impose_pbc (&to_check, x, y, z); 
        // std::cout << "Growth location is "; print(to_check); 

		if (checkOccupancyHead(&to_check, PVec, index_of_polymer, index_of_monomer)){
            // std::cout << "This location is occupied! - "; print(to_check); 
            // std::cout << "Moving on to another growth spot..." << std::endl;
			continue; 
		}

		else {

			(*PVec)[index_of_polymer].chain[index_of_monomer+1].coords = to_check;
            // std::cout << "Moving in deeper...\n\t"; 
			HeadSpin (PVec, index_of_polymer, index_of_monomer+1, deg_poly, x, y, z, b, IMP_BOOL, first_entry_bool);

			if (*b){
				break;
			} 
			else {
                // std::cout << "We went down a rabbithole, but it failed. So moving on..." << std::endl;
				continue; 
			}

		}

	}

	if ( !(*b) ){
		*IMP_BOOL = false; 
	}

	return ;

}



void HeadSpin_Rosenbluth(std::vector <Polymer>* PVec, int index_of_polymer, int index_of_monomer, int deg_poly,int x, int y, int z, bool* IMP_BOOL, bool* first_entry_bool, double* rweight){
    
     if (*first_entry_bool){
     	(*rweight) = 1; 
	    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	    std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
        *first_entry_bool = false; 
    }

	if (index_of_monomer == deg_poly-1){
		// std::cout << "You have reached the final spot of head spin!" << std::endl;
		*IMP_BOOL = true;
		(*PVec)[index_of_polymer].ChainToConnectivityMap(); 
		return ;
	}

	// unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	// std::shuffle (adrns.begin(), adrns.end(), std::default_random_engine(seed));
    
    // std::cout << "Current pivot point is "; print(  (*PVec) [index_of_polymer].chain[index_of_monomer].coords );
	int rw_tmp = 0; 
	std::vector <std::array <int,3>> ind_v; 

	for (std::array <int,3>& d: adrns){
        // std::cout << "d is "; print(d); 
		std::array <int, 3> to_check = add_arrays ( &( (*PVec) [index_of_polymer].chain[index_of_monomer].coords ), &d);

		impose_pbc (&to_check, x, y, z); 
        // std::cout << "Growth location is "; print(to_check); 

		if (checkOccupancyHead(&to_check, PVec, index_of_polymer, index_of_monomer)){
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
		
		(*PVec)[index_of_polymer].chain[index_of_monomer+1].coords = ind_v[ rng_uniform(0, rw_tmp-1) ]; 
		(*rweight) = (*rweight) * rw_tmp/6; 
		HeadSpin_Rosenbluth (PVec, index_of_polymer, index_of_monomer+1, deg_poly, x, y, z, IMP_BOOL, first_entry_bool, rweight);

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

void OrientationFlip (std::vector <Particle>* SolvVect, int x, int y, int z, int size_of_region) {

	// pick a random point on the lattice 
	std::array <int,3> rpoint = {rng_uniform(0,x-1), rng_uniform(0,y-1), rng_uniform (0,z-1)}; 
	std::array <int,3> loc; 

	for (int i{0}; i < size_of_region; ++i){
		for (int j{0}; j < size_of_region; ++j){
			for (int k{0}; k < size_of_region; ++k){

				// std::cout << "i = " << i << ", j = " << j << ", k = " << k <<"." << std::endl;

				loc[0] = (rpoint[0]+i)%x;
				loc[1] = (rpoint[1]+j)%y;
				loc[2] = (rpoint[2]+k)%z;

				for (Particle& p: (*SolvVect)){

					if (p.coords == loc){
						p.orientation = orientations [ rng_uniform(0,5) ]; 
						break;
					}
				}
			}
		}
	}
	return; 
}


void OrientationFlipLocal ( std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z){


	// obtain the list of solvent particles neighboring the polymer 
	std::vector <int> solvent_indices;
    solvent_indices.reserve((*Solvent).size()); 
      
	for (const Polymer& pmer: (*Polymers)) {
		for (const Particle& p: pmer.chain) {

			std::array <std::array <int,3>, 6> ne_list = obtain_ne_list (p.coords, x, y, z); 
            for (const std::array <int,3>& n: ne_list){ 
                
                for (int j{0}; j<static_cast<int>( (*Solvent).size() ); ++j){
                    if (n == (*Solvent)[j].coords){
                        solvent_indices.push_back(j); 
                    }
                }                
            }	
		}
	}

	// get rid of repeated indices
	std::sort ( solvent_indices.begin(), solvent_indices.end() ); 
	solvent_indices.erase ( std::unique ( solvent_indices.begin(), solvent_indices.end() ), solvent_indices.end() );

	// std::cout << "size of solvent_indices is " << solvent_indices.size() << std::endl; 

	int stopping_idx = rng_uniform(0, static_cast<int>(solvent_indices.size() - 1) ); 

	// std::cout << "stopping_idx is " << stopping_idx << std::endl;
	int counter = 0; 
	for (int idx: solvent_indices){

		if (counter == stopping_idx){
			break;
		}
		// std::cout << "orientation is " << (*Solvent)[idx].orientation <<", location is: "; 
		// print( (*Solvent)[idx].coords ); 
		(*Solvent)[idx].orientation = orientations [ rng_uniform (0, 5) ]; 
		++counter; 

	}

	return; 

}


void PolymerFlip ( std::vector <Polymer>* PolVec ){

    for (Polymer& pmer: (*PolVec) ){
        for (Particle& p: pmer.chain){
            p.orientation = orientations [ rng_uniform(0,5) ]; 
        }
    }
    // std::cout << "Performed flip." << std::endl;
    return; 
}

void PolymerFlipLocal ( std::vector <Polymer>* PolVec ){

    for (Polymer& pmer: (*PolVec) ){
        
        int idx = rng_uniform (0, static_cast<int> (pmer.chain.size()-1) ) ;
        pmer.chain[idx].orientation = orientations [ rng_uniform(0,5) ] ; 
        break; 
        
    }
    // std::cout << "Performed flip." << std::endl;
    return; 
}


//============================================================
//============================================================
// 
// NAME OF FUNCTION:                  
//
// PARAMETERS: *PolymerVector, *SolvVector, x, y, z, v, *BOOL
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
/*
std::vector <Polymer> MoveChooser(std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, bool v, bool* IMP_BOOL){

    int index = rng_uniform(0, static_cast<int>((*Polymers).size())-1); 
    std::vector <Polymer> NewPol = (*Polymers); 
    int r = rng_uniform(1, 8);
    switch (r) {
        case (1):
            if (v){
               printf("Performing end rotations.\n"); 
            }
            // 
            
            NewPol = EndRotation(Polymers, Solvent, index, x, y, z, IMP_BOOL);
            break;     
        
        case (2):
            if (v){
               printf("Performing crank shaft.\n"); 
            }
            
            NewPol = CrankShaft(Polymers, Solvent, index, x, y, z, IMP_BOOL);
            break; 

        case (3):
            if (v){
               printf("Performing reptation.\n"); 
            }
            
            NewPol = Reptation(Polymers, Solvent, index, x, y, z, IMP_BOOL); 
            break; 

        case (4):
            if (v){
               printf("Performing kink jump.\n"); 
            }

            NewPol = KinkJump(Polymers, Solvent, index, x, y, z, IMP_BOOL);
            break; 

        case (5):
        	if (v) {
        		printf("Performing configuration sampling. \n"); 
        		std::cout << "index of polymer is " << index << std::endl;
        	}
        	NewPol = *Polymers; 
        	ChainRegrowth(&NewPol, Solvent, index, x, y, z, IMP_BOOL ); 
        	break;

        case (6): 
        	if (v){
        		printf("Performing solvent orientation flips. \n");
        	}
        	OrientationFlip(Solvent, x, y, z, 4); 
        	// OrientationFlipLocal (Polymers, Solvent, x, y, z); 
            // std::cout << "performed flip local." << std::endl;
        	break; 

        case (7):
        	if (v){
        		printf("Performing translation. \n");
        	}
        	NewPol = Translation(Polymers, Solvent, index, x, y, z, IMP_BOOL);
        	break;
        
        case (8):
            if (v){
                printf("Performing a polymer orientation flip.\n"); 
            }
            NewPol = *Polymers;
            PolymerFlip(&NewPol); 
            break; 
    }

    return NewPol;
}
*/
//////////////////////////////////////////////////////////////

std::vector <Polymer> MoveChooser_Rosenbluth (std::vector <Polymer>* Polymers, std::vector <Particle>* Solvent, int x, int y, int z, bool v, bool* IMP_BOOL, double* rweight){

    int index = rng_uniform(0, static_cast<int>((*Polymers).size())-1); 
    std::vector <Polymer> NewPol = (*Polymers); 
    int r = rng_uniform(1, 9);
    switch (r) {
        case (1):
            if (v){
               printf("Performing end rotations.\n"); 
            }
            // 
            
            NewPol = EndRotation_Rosenbluth (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight);
            break;     
        
        case (2):
            if (v){
               printf("Performing crank shaft.\n"); 
            }
            
            NewPol = CrankShaft_Rosenbluth (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight);
            break; 

        case (3):
            if (v){
               printf("Performing reptation.\n"); 
            }
            
            NewPol = Reptation_Rosenbluth (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight); 
            break; 

        case (4):
            if (v){
               printf("Performing kink jump.\n"); 
            }

            NewPol = KinkJump_Rosenbluth (Polymers, Solvent, index, x, y, z, IMP_BOOL, rweight);
            break; 

        case (5):
        	if (v) {
        		printf("Performing configuration sampling. \n"); 
        		std::cout << "index of polymer is " << index << std::endl;
        	}
        	NewPol = *Polymers; 
        	ChainRegrowth_Rosenbluth (&NewPol, Solvent, index, x, y, z, IMP_BOOL, rweight ); 
        	break;

        case (6): 
        	if (v){
        		printf("Performing local solvent orientation flips. \n");
        	}
        	// OrientationFlip(Solvent, x, y, z, 4); 
        	OrientationFlipLocal (Polymers, Solvent, x, y, z); 
        	(*rweight) = 1; 
            // std::cout << "performed flip local." << std::endl;
        	break; 
        
        case (7):
            if (v){
                printf("Performing solvent orientation flips. \n"); 
            }
            OrientationFlip(Solvent, x, y, z, 4);
            (*rweight) = 1; 
            break; 

        case (8):
        	if (v){
        		printf("Performing translation. \n");
        	}
        	NewPol = Translation(Polymers, Solvent, index, x, y, z, IMP_BOOL);
        	break;

        case (9):
        	if (v){
        		printf("Performing polymer orientation flips. \n");
        	}
        	NewPol = *Polymers; 
        	PolymerFlipLocal (&NewPol);
        	(*rweight) = 1; 
        	break; 
    }

    return NewPol;
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


std::vector <Particle> CreateSolventVector(int x, int y, int z, std::vector <Polymer>* PolymerVector){

	// create the lattice 
	std::vector <std::array <int,3>> lattice = create_lattice_pts(x,y,z); 

	// remove all the points where monomer segments are present 
	int monomer_count = 0; 

	for (Polymer pmer: *PolymerVector){
		for (Particle p: pmer.chain){

			++monomer_count; 
			lattice.erase(std::find(lattice.begin(), lattice.end(), p.coords)); 

		}
	}

	int nsolpart = x*y*z - monomer_count; 

	std::vector <Particle> SolvPartVector; 
	SolvPartVector.reserve(nsolpart); 

	for (int i{0}; i<nsolpart; ++i){

		Particle p = Particle (lattice[i], "solvent", orientations[ rng_uniform (0,5) ]);
		SolvPartVector.push_back(p);  

	}

	return SolvPartVector; 

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


	NewPol[index].ChainToConnectivityMap(); 
	*IMP_BOOL = true; 

	return NewPol;
}



//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//             End of Translation
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


