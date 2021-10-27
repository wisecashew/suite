#include <iostream>
#include <vector> 
#include <string> 
#include <map>
#include <algorithm> 
#include <random>
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


std::vector <std::vector <int> > Particle::nlist(int x_len, int y_len, int z_len){
    return obtain_ne_list(this->loc, x_len, y_len, z_len); 
}


// END OF METHODS FOR CLASS PARTICLE
// #####################################

// methods for class Polymer :

void Polymer::print_loc(){
    for (Particle p: this->chain){
        p.print_loc(); 
    }
    return; 
}


// only works for linear polymers 
void Polymer::obtain_connectivity() {
    int size = this->chain.size(); 
    // std::cout << "size of chain is " << size << std::endl; 
    for (int i{0}; i<size; i++){
        if (i==0){
            
            this->conn[this->chain.at(i)] = { (this->chain.at(i+1)) };
        }
        else if (i==(size-1)){
            // std::vector <Particle> pvec = { (this->chain.at(i-1)) }; 
            this->conn[this->chain.at(i)] = { (this->chain.at(i-1)) }; 
        }
        else {
            std::vector <Particle> pvec = {this->chain.at(i-1), this->chain.at(i+1)};
            this->conn[this->chain.at(i)] = {this->chain.at(i-1), this->chain.at(i+1)}; 
        }
    }
    return;
}; 

// get the plocs attribute 
void Polymer::get_plocs(){
    this->p_locs.reserve(this->chain.size());

    for (Particle p: this->chain ){
        this->p_locs.push_back(p.loc);
    }
    return ;
}


// find if there are any kinks in the polymer structure 
std::vector <int> Polymer::find_kinks(){

    std::vector <int> kink_indices; 
    
    // obtain all location where kinks exist in the polymer 
    for (int i{0}; i<this->dop-2; i++){
        std::vector <int> v1 = subtract_vectors(&(this->p_locs.at(i+1)), &(this->p_locs.at(i)));
        std::vector <int> v2 = subtract_vectors(&(this->p_locs.at(i+2)), &(this->p_locs.at(i+1)));

        if (v1==v2){            
        }
        else {
            kink_indices.push_back(i);  
        }

    }

    return kink_indices; 

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


// assuming seed is a valid location, insert a polymer into the grid 
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



// solvate every unoccupied site with a solvent particle 
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







// BEGIN MONTE CARLO MOVES 

// START THE END ROTATION 

void Grid::ZeroIndexRotation(){

    Particle pNextToEdge = this->polymer.chain.at(1);         // .at(.size()-2) for the other edge 

    // get this particles neighborlist 
    std::vector <std::vector <int> > ne_list = pNextToEdge.nlist(this->x_len, this->y_len, this->z_len); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(0).loc ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(2).loc ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the self-avoidance criterion from misc.cpp 
    

    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 
    for (std::vector <int> to_rot:  ne_list ) {
        if (check_avoidance( to_rot, this->polymer.p_locs) ){
            this->polymer.chain.at(0).loc = to_rot;
            this->polymer.p_locs.at(0) = to_rot;
            this->polymer.conn[pNextToEdge].at(0) = polymer.chain.at(0);
            break;
        }
        else{
            std::cout << "occupied..." << std::endl; 
        }
    }

    return; 
}


void Grid::FinalIndexRotation(){
    int dop = polymer.chain.size(); 
    Particle pNextToEdge = this->polymer.chain.at(dop-2);         // .at(.size()-2) for the other edge 
    // get this particles neighborlist 
    std::vector <std::vector <int> > ne_list = pNextToEdge.nlist(this->x_len, this->y_len, this->z_len); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(dop-1).loc ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(dop-3).loc ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the self-avoidance criterion from misc.cpp 
    

    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 
    for (std::vector <int> to_rot:  ne_list ) {
        if (check_avoidance( to_rot, this->polymer.p_locs) ){
            this->polymer.chain.at(dop-1).loc = to_rot;
            this->polymer.p_locs.at(dop-1) = to_rot;
            this->polymer.conn[pNextToEdge].at(0) = polymer.chain.at(dop-1);
            break;
        }
        else{
            std::cout << "occupied..." << std::endl; 
        }
    }



    return ; 
}


void Grid::end_rotation() {
    // consider the polymer in the grid 

    // consider the particle right BESIDE the edge 
    int r = (rand()%2); 
    std::cout <<"random number r is " << r <<"!"<< std::endl;


    if (r==0){
        std::cout << "Zero Index is being rotated!" << std::endl;
        this->ZeroIndexRotation(); 
    }    

    else{
        std::cout << "Final index is being rotated!" << std::endl;
        FinalIndexRotation();
    }

    return; 
}

// END OF END ROTATION 

// START OF KINK JUMP 




