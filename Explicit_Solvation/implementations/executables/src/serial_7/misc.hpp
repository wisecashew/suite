#ifndef _MISC_H_
#define _MISC_H_
#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <fstream>
#include <string>
#include <array>
#include <utility>
#include <functional>
#include <array>
#include <random>
#include <chrono>
#include <set>
#include <cmath>
#include <chrono>
#include <sstream>
#include <regex>
#include <tuple>
#include <iterator>
#include <unordered_set>
#include <algorithm>
#include <filesystem>
#include <limits>
#include <type_traits>
#include <cassert>
#include <getopt.h> 
#include <stdlib.h> 

constexpr int CONTACT_SIZE = 8;

// this is a bunch of miscellaneous functions I use 
// check if the new location added to random walk has been visited before 
bool check_avoidance             (std::vector <int> to_check, std::vector<std::vector <int>> loc_list); 

// run a sarw without checking for pbc 
void sarw                        (std::vector<std::vector<int>>* loc_list, int dop); 

// run a sarw with periodic boundary conditions 
void sarw_pbc                    (std::vector<std::vector<int>>* loc_list, int dop, int x_len, int y_len, int z_len); 

// get the arccos of an input
double branchless_acos(double input);

// get the norm of an array
double norm(std::array<double,3> arr);

// given a list of locations for a particular type of Particle, create a vector of Particles 
// std::vector <Particle> loc2part           (std::vector <std::vector <int>> loc_list, std::string s); 

// given a vector of Particles, obtain a vector of locations 
// std::vector <std::vector <int> > part2loc (std::vector <Particle> pVec); 

// run an acceptance criterion for two polymers 
bool acceptance(int dE, double kT); 
bool metropolis_acceptance(double E1, double E2, double kT);

// sending a string to a file
void StringToFile(std::string filename, std::string to_send);

// splitting up a string into components separated by a delimiter
std::vector<std::string> split (const std::string &s, char delim);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// imposing modified modulo methods 
void   modified_direction (std::array<int,3>* a, int x, int y, int z);
int    modified_modulo    (int divident, int divisor);
double modified_modulo    (double divident, int divisor);
//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// imposing periodic boundary conditions 
void impose_pbc (std::vector <int>* vect , int x_len, int y_len, int z_len); 
void impose_pbc (std::array  <int,3>* arr, int x_len, int y_len, int z_len);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// obtain neighbor list 
std::vector <std::vector <int>>      obtain_ne_list (std::vector <int>   loc, int x_len, int y_len, int z_len); 
std::array  <std::array <int,3>, 26> obtain_ne_list (std::array  <int,3> loc, int x_len, int y_len, int z_len);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
void create_linked_list ( std::vector<std::array<int,3>> v1, std::vector<std::array<int,3>> v2, std::vector <std::array<int,3>> link, std::vector <std::vector <std::array<int,3>>>* master_linked_list, int beginning);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// random number generation 
int    rng_uniform (int start, int end);
double rng_uniform (double start, double end);


// extract the numbers from a string
double NumberExtractor(std::string s);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
int                  lattice_index ( std::array<int,3> location, int y, int z );
std::array <int,3>   location      ( int lattice_index, int x, int y, int z);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// std::vector and std::array arithmetic
// adding the two 
/*
std::vector <int>      add_containers (std::vector <int>* v1, std::vector <int>* v2); 
std::array  <int,3>    add_containers  (std::array <int,3>* a1, std::array <int,3>* a2);
std::array  <double,3> add_containers  (std::array <int,3>* a1, std::array <double,3>* a2);
std::array  <double,3> add_containers  (std::array <double,3>* a1, std::array <double,3>* a2);
std::array  <int,CONTACT_SIZE>    add_containers  (std::array <int,CONTACT_SIZE>*     a1, std::array <int,CONTACT_SIZE>*     a2);
std::array  <double,CONTACT_SIZE> add_containers  (std::array <double,CONTACT_SIZE>*  a1, std::array <double,CONTACT_SIZE>*  a2);
std::array  <double,CONTACT_SIZE> add_containers  (std::array <double,CONTACT_SIZE>   a1, std::array <double,CONTACT_SIZE>   a2);


// subtracting the two 
std::vector <int>      subtract_vectors (std::vector <int>* v1,     std::vector <int>* v2); 
std::array  <int,3>    subtract_containers  (std::array <int,3>* a1,    std::array <int,3>* a2);
std::array  <double,3> subtract_containers  (std::array <double,3>* a1, std::array <double,3>* a2);
std::array  <int,CONTACT_SIZE>    subtract_containers (std::array <int,CONTACT_SIZE>*    a1, std::array <int,CONTACT_SIZE>*    a2);
std::array  <double,CONTACT_SIZE> subtract_containers (std::array <double,CONTACT_SIZE>* a1, std::array <double,CONTACT_SIZE>* a2);
std::array  <double,CONTACT_SIZE> subtract_containers (std::array <double,CONTACT_SIZE>  a1, std::array <double,CONTACT_SIZE>  a2); 
*/

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// getting distance between two points 
double              distance_between_points (std::array <int,3>*    a1, std::array <int,3>* a2   , int xlen, int ylen, int zlen);
double              distance_between_points (std::array <double,3>* a1, std::array <double,3>* a2, int xlen, int ylen, int zlen);
double              take_dot_product        (int o1, int o2); 
double              take_dot_product        (const std::array <double,3>& o1, const std::array <double,3>& o2); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// check input of main driver code 
void input_parser (int dfreq, 
	int lfreq, 
	int max_iter, 
	bool r, 
	std::string positions, 
	std::string topology, 
	std::string dfile, 
	std::string efile, 
	std::string mfile, 
	std::string stats_file, 
	std::string lattice_file_read,
	std::string lattice_file_write,
	std::string ssfile); 

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
void filter_coords       (const std::string& filename, int cutoffStep);
void filter_csv          (const std::string& filename, int cutoffStep);
void filter_orientations (const std::string& filename, int cutoffStep);

//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
// print methods  
// print out a vector 
void print (std::vector <int> v, std::string c="\n"); 
void print (std::vector <double> v, std::string c="\n");

// print out an array 
void print (std::array <int,3> a, std::string c="\n"); 
void print (std::array <int,6> a, std::string c="\n"); 
void print (std::array <double,3> a, std::string c="\n");
void print (std::array <double,4> a, std::string c="\n"); 
void print (std::array <double, 6> v, std::string c="\n");
void print (std::array <std::array<int,3>,6> aa );

// print out a list of vectors
void print (std::vector <std::vector <int>> v); 
void print (std::vector <std::vector <double>> v); 
void print (std::vector <std::vector <std::array<int,3>>> v);
void print (std::vector <std::vector <int>>vv);
void print (std::vector <std::vector <double>>vv);
void print (std::vector <std::string>vv);
void print ( std::array <std::array<int,3>,6> aa );
void print ( std::vector <std::array<int,3>> aa );
void print ( std::vector <std::vector <std::array<int,3>>> V);

template <typename T>
void print (T vec, std::string end="\n"){

	for (int i{0}; i < static_cast<int>(vec.size()); ++i){
		std::cout << vec[i] << " | "; 
	}
	std::cout << end; 
	return;
}

template <typename T>
void reset(T &x){
	x = T();
}

// template <typename ContainerType>
// ContainerType add_containers(const ContainerType& c1, const ContainerType& c2);

// template <typename ContainerType>
// ContainerType subtract_containers(const ContainerType& c1, const ContainerType& c2);

template <typename ContainerType>
ContainerType add_containers(const ContainerType& c1, const ContainerType& c2) {
	// Ensure both containers have the same size
	assert(c1.size() == c2.size());

	// Deduce the result type using the common type between the two containers' elements

	// Create a container of the same size for the result
	ContainerType result = c1;          // Can be array or vector; assumes both containers are the same type
	// std::fill(result.begin(), result.end(), 0); // result.fill(CommonType(0));    // For arrays

	// Perform element-wise addition
	for (std::size_t i = 0; i < c1.size(); ++i) {
		result[i] = (c1[i]) + (c2[i]);
	}

	return result;
}

template <typename ContainerType>
ContainerType subtract_containers(const ContainerType& c1, const ContainerType& c2) {
	// Ensure both containers have the same size
	assert(c1.size() == c2.size());

	// Deduce the result type using the common type between the two containers' elements

	// Create a container of the same size for the result
	ContainerType result = c1;          // Can be array or vector; assumes both containers are the same type
	// std::fill(result.begin(), result.end(), 0); // result.fill(CommonType(0));    // For arrays

	// Perform element-wise addition
	for (std::size_t i = 0; i < c1.size(); ++i) {
		result[i] = (c1[i]) - (c2[i]);
	}

	return result;
}


template <typename Scalar, typename T>
std::vector<decltype(std::declval<Scalar>() * std::declval<T>())> scale_vectors(Scalar scalar, const std::vector<T>& container) {
    using ResultType = decltype(scalar * std::declval<T>());
    std::vector<ResultType> result;
    result.reserve(container.size());
    
    for (const T& value : container) {
        result.push_back(scalar * value);
    }
    return result;
}

template <typename Scalar, typename T, std::size_t N>
std::array<decltype(std::declval<Scalar>() * std::declval<T>()), N> 
scale_arrays (Scalar scalar, const std::array<T, N>& container) {
    using ResultType = decltype(scalar * std::declval<T>());
    std::array<ResultType, N> result{};
    
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = scalar * container[i];
    }
    
    return result;
}


void resetting_containers(std::array<double,CONTACT_SIZE>* cs1_i, std::array<double,CONTACT_SIZE>*cs1_f,
std::array<double,CONTACT_SIZE>* cs2_i, std::array<double,CONTACT_SIZE>*cs2_f,
std::array<double,CONTACT_SIZE>* cpair_i, std::array<double,CONTACT_SIZE>* cpair_f,
double* Es1_i, double* Es1_f, double* Es2_i, double* Es2_f,
double* Epair_i, double* Epair_f);

inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
{
	s.erase(s.find_last_not_of(t) + 1);
	return s;
}

// trim from beginning of string (left)
inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
{
	s.erase(0, s.find_first_not_of(t));
	return s;
}

// trim from both ends of string (right then left)
inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
{
	return ltrim(rtrim(s, t), t);
}


#endif 
