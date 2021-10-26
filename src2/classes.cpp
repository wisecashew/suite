#include <iostream>
#include <vector> 
#include <string> 
#include <algorithm> 
#include "classes.h"
#include "misc.h" 


// methods for class Particle: 
// ####################################


void Particle::print_loc(){
    print(this->loc);
    return;
}

std::vector <Particle> loc2part(std::vector<std::vector<int>> loc_list, std::string s ){
    std::vector <Particle> pvec; 
    for (std::vector <int> loc: loc_list){
        Particle p(loc, s); 
        pvec.push_back(p);
    }; 

    return pvec; 
}


// END OF METHODS FOR CLASS PARTICLE
// #####################################

// methods for class Polymer :

void Polymer::print_loc(){
    for (Particle p: this->chain){
        p.print_loc(); 
        std::cout <<"particle type is " << p.ptype << std::endl;
    }
    return; 
}


// END OF CLASS POLYMER 
// ######################################


// methods for class Solvent:

void Solvent::print_loc(){
    for (Particle p: this->solvation){
        p.print_loc();
        std::cout <<"particle type is " << p.ptype << std::endl;
    }
    return; 
}


// END OF CLASS SOLVENT
// #####################################


// methods for class Grid 
void Grid::print_polymer(){
    std::cout << "Printing out locations of monomeric units of polymer on the lattice:" << std::endl; 
    this->polymer.print_loc(); 
    
    return; 
};


void Grid::print_solvent(){
    std::cout << "Printing out locations of solvent particles currently on the lattice: " << std::endl;
    this->solvent.print_loc(); 

    return;
}; 

void Grid::print_occupied(){
    std::cout << "Printing out locations of currently occupied sites on the lattice: " << std::endl;
    for (Particle p: this->occupied){
        print(p.loc);
    }
}

// std::vector <int> seed {0,0,0}; 


// assuming seed is a valid location 
void Grid::polymer_insertion(int dop, std::vector <int> seed ){

    std::vector <std::vector <int>> loc_list = { seed }; 

    sarw_pbc(&(loc_list), dop, this->x_len, this->y_len, this->z_len );
    // the assumption of this method is that the seed is an available spot 
    
    // this->occupied = loc2part(loc_list);
    this->polymer = loc2part(loc_list, "polymer"); 

    for (Particle p: this->polymer.chain){
        this->occupied.push_back(p); 
    }

    return; 

};




void Grid::solvate(){
    std::vector <std::vector <int>> ulattice = create_lattice_pts(this->x_len, this->y_len, this->z_len); 
    
    std::vector <std::vector <int>> solv_pts; 
    
    solv_pts.reserve(ulattice.size()-this->polymer.chain.size()); 

    for (Particle p: polymer.chain){
        ulattice.erase(std::find(ulattice.begin(), ulattice.end(), p.loc)); 
    };

    this->solvent = loc2part(ulattice, "solvent"); //loc2part(ulattice); 

    for (Particle p: this->solvent.solvation){
        this->occupied.push_back(p);  
    }

    return; 

};



