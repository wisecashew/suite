#include <iostream>
#include <vector> 
#include <algorithm> 
#include "classes.h"
#include "misc.h" 


// methods for class Particle: 
// ####################################


void Particle::print_loc(){
    print(this->loc);
    return;
}


// END OF METHODS FOR CLASS PARTICLE
// #####################################

// methods for class Polymer :




// END OF CLASS POLYMER 
// ######################################


// methods for class Grid 
void Grid::get_loclist(){
    std::vector <std::vector <int>> loclist; 
    for (Particle p: this->GRID){
        loclist.push_back(p.loc);
    }
    this->loc_list = loclist;
    return; 
};


void Grid::print_loclist(){
    std::cout << "Printing out locations of particles currently in the Grid:" << std::endl;
    for(std::vector <int> v: this->loc_list){
        print(v); 
    }
    return ;
};


// std::vector <int> seed {0,0,0}; 


// assuming seed is a valid location 
void Grid::polymer_insertion(int dop, std::vector <int> seed = {0,0,0}){
    // this is a problematic thing I am doing 
    this->loc_list = { seed }; 

    sarw_pbc(&(this->loc_list), dop, this->x_len, this->y_len, this->z_len );
    // the assumption of this method is that the seed is an available spot 
    
    return; 

};

std::vector <std::vector<int>> Grid::solvate(){
    std::vector <std::vector <int>> ulattice = create_lattice_pts(this->x_len, this->y_len, this->z_len); 
    // ulattice = create_lattice_pts(this->x_len, this->y_len, this->z_len); 
    std::vector <std::vector <int>> solv_pts; 
    solv_pts.reserve(ulattice.size()-this->loc_list.size()); 

    for (std::vector <int> v: loc_list){
        ulattice.erase(std::find(ulattice.begin(), ulattice.end(), v)); 
    }

    return ulattice; 

};

