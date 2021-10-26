#include <iostream>
#include <vector> 
#include "misc.h"
#include "grid.h"


void Grid::get_loclist(){
	std::vector <std::vector <int>> loclist; 
	for (Particle p: this->GRID){
		loclist.push_back(p.loc);
	}
	this->loc_list = loclist;
	return; 
};

void Grid::print_loclist(){
	for(std::vector <int> v: this->loc_list){
		print(v); 
	}
	return ;
};



std::vector <int> seed {0,0,0}; 


// assuming seed is a valid location 
void Grid::polymer_insertion(int dop, std::vector <int> seed){
	// this is a problematic thing I am doing 
	this->loc_list = { seed }; 

	sarw_pbc(&(this->loc_list), dop, this->x_len, this->y_len, this->z_len );
	// the assumption of this method is that the seed is an available spot 
	
	return; 

};
