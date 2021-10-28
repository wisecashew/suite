#ifndef _MISC_H_
#define _MISC_H_
#include "classes.h"

// this is a bunch of miscellaneous functions I use 
// to manipulate std::vectors. Furthermore, this will be used 
// as a skeleton for building functions for grids 

// add two vectors 
std::vector <int> add_vectors(std::vector <int>* v1, std::vector <int>* v2); 

// subtract two vectors 
std::vector <int> subtract_vectors(std::vector <int>* v1, std::vector <int>* v2); 

// check if the new location added to random walk has been visited before 
bool check_avoidance(std::vector <int> to_check, std::vector<std::vector <int>> loc_list); 

// print out a vector 
void print(std::vector <int> v); 

// print out a list of vectors
void print(std::vector <std::vector <int>> v); 

// impose periodic boundary conditions on a vector 
void impose_pbc(std::vector <int>* vect, int x_len, int y_len, int z_len); 


// run a sarw without checking for pbc 
void sarw(std::vector<std::vector<int>>* loc_list, int dop); 

// run a sarw with periodic boundary conditions 
void sarw_pbc(std::vector<std::vector<int>>* loc_list, int dop, int x_len, int y_len, int z_len); 


// create a lattice 
std::vector <std::vector <int>> create_lattice_pts(int x_len, int y_len, int z_len); 


// obtain the neighbor vector (list) 
// implicitly assumed in 3d 

std::vector <std::vector <int>> obtain_ne_list(std::vector <int> loc, int x_len, int y_len, int z_len); 

// given a list of locations for a particular type of particle, create a vector of particles 
std::vector <Particle> loc2part (std::vector <std::vector <int>> loc_list, std::string s); 

// given a vector of particles, obtain a vector of locations 
std::vector <std::vector <int> > part2loc (std::vector <Particle> pVec); 

#endif 
