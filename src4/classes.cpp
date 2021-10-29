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











// END OF METHODS FOR CLASS PARTICLE 
//////////////////////////////////////////////////////////////////////////////////////////














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


// find if there are any kinks in the polymer structure 
std::vector <int> Polymer::find_kinks(){
    int dop = this->chain.size(); 
    std::vector <int> kink_indices; 
    
    // obtain all location where kinks exist in the polymer 
    for (int i{0}; i< dop-2; i++){
        std::vector <int> v1 = subtract_vectors(&(this->chain.at(i+1).loc), &(this->chain.at(i).loc) );
        std::vector <int> v2 = subtract_vectors(&(this->chain.at(i+2).loc), &(this->chain.at(i+1).loc) );

        if (v1==v2){
            continue;            
        }
        else {
            kink_indices.push_back(i);  
        }

    }

    if (kink_indices.size() == 0){
        std::cout << "no kinks in current structure..." << std::endl;
    }

    return kink_indices; 

}


// find if there any cranks in the polymer structure 
std::vector <int> Polymer::find_cranks(){
    int dop = this->chain.size(); 
    std::vector <int> crank_indices;
    // obtain all locations where kinks exist in the polymer 
    for (int i{0}; i< dop-3; i++){
        std::vector <int> v1 = subtract_vectors(&(this->chain.at(i+1).loc), &(this->chain.at(i).loc) ); 
        std::vector <int> v2 = subtract_vectors(&(this->chain.at(i+2).loc), &(this->chain.at(i+3).loc) );

        //

        if (v1==v2){
            crank_indices.push_back(i); 
        } 
        else {
            continue;
        }
    }

    if (crank_indices.size() == 0) {
        std::cout << "there are no cranks in the current structure..." << std::endl;
    } 

    return crank_indices;
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
    this->polymer.chain = loc2part(loc_list, "polymer"); 

    for (Particle p: this->polymer.chain){
        this->occupied.push_back(p); 
    }
    this->polymer.obtain_connectivity();
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


    // get this particles neighborlist 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(this->polymer.chain.at(1).loc, this->x_len, this->y_len, this->z_len); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(0).loc ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(2).loc ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the self-avoidance criterion from misc.cpp 
    

    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 
    for (std::vector <int> to_rot:  ne_list ) {
        if (check_avoidance( to_rot, part2loc(this->polymer.chain)) ){

            this->polymer.chain.at(0).loc = to_rot;
            this->polymer.obtain_connectivity(); 
        }
        else{
            std::cout << "occupied..." << std::endl; 
        }
    }


    return; 
}


void Grid::FinalIndexRotation(){

    int dop = polymer.chain.size(); 
    
    // get this particles neighborlist 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(this->polymer.chain.at(dop-2).loc, this->x_len, this->y_len, this->z_len); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(dop-1).loc ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(dop-3).loc ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the self-avoidance criterion from misc.cpp 
    

    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 
    for (std::vector <int> to_rot:  ne_list ) {
        if (check_avoidance( to_rot, part2loc(this->polymer.chain) ) ){

            this->polymer.chain.at(dop-1).loc = to_rot;
            this->polymer.obtain_connectivity(); 

        }
        else{
            std::cout << "occupied..." << std::endl; 
        }
    }



    return ; 
}


void Grid::end_rotation() {

    int r = (rand()%2); 
    std::cout <<"random number r is " << r <<"!"<< std::endl;


    if (r==0){
        std::cout << "Zero Index is being rotated!" << std::endl;
        this->ZeroIndexRotation(); 
    }    

    else{
        std::cout << "Final index is being rotated!" << std::endl;
        this->FinalIndexRotation();
    }

    return; 
}

// END OF END ROTATION 

// START OF KINK JUMP 
void Grid::kink_jump() {

    std::vector <int> k_idx = this->polymer.find_kinks();

    std::shuffle(std::begin(k_idx), std::end(k_idx), std::default_random_engine() ); 

    for (int idx: k_idx){
        std::cout << "idx right before kink spot is " << idx << std::endl;
        std::vector <int> d1 = subtract_vectors(&(this->polymer.chain.at(idx+1).loc), &(this->polymer.chain.at(idx).loc));
        std::vector <int> d2 = subtract_vectors(&(this->polymer.chain.at(idx+2).loc), &(this->polymer.chain.at(idx+1).loc));

        std::vector <int> to_check = add_vectors(&(this->polymer.chain.at(idx).loc), &d2); 
        if (check_avoidance(to_check, part2loc(this->polymer.chain)) ) {


            this->polymer.chain.at(idx+1).loc = to_check; 
            this->polymer.obtain_connectivity();
            return; 
        } 
        else {
            std::cout << "occupied..." << std::endl;
        }
    }

    return;

}

// END OF KINK JUMP 

// START OF CRANK SHAFT 
void Grid::crank_shaft() {
    std::vector <int> c_idx = this->polymer.find_cranks(); 

    std::shuffle(std::begin(c_idx), std::end(c_idx), std::default_random_engine() ); 

    for (int idx: c_idx){
        std::vector <int> d1 = subtract_vectors(&(this->polymer.chain.at(idx+3).loc), &(this->polymer.chain.at(idx+2).loc)); 
        impose_pbc(&d1, this->x_len, this->y_len, this->z_len); 

        std::vector <int> to_check_1 = add_vectors(&(this->polymer.chain.at(idx).loc), &d1);
        std::vector <int> to_check_2 = add_vectors(&(this->polymer.chain.at(idx+3).loc), &d1); 

        impose_pbc(&to_check_1, this->x_len, this->y_len, this->z_len);
        impose_pbc(&to_check_2, this->x_len, this->y_len, this->z_len);

        if ((check_avoidance(to_check_1, part2loc(this->polymer.chain)) ) && (check_avoidance(to_check_2, part2loc(this->polymer.chain)) ) ) {

            // reupdate polymer now that I have changed it 
            this->polymer.chain.at(idx+1).loc = to_check_1;
            this->polymer.chain.at(idx+2).loc = to_check_2; 
            this->polymer.obtain_connectivity(); 
            return;
        }
        else{
            std::cout << "occupied..." << std::endl;
        }

    }    
    return ;
}

// END OF CRANK SHAFT 

// START OF REPTATION MOVE 

void Grid::FinalToZero(){
    
    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->polymer.chain.at(0).loc, this->x_len, this->y_len, this->z_len); 

    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(1).loc ), ne_list.end() );

    this->polymer.chain.pop_back();

    std::vector <std::vector <int>> popFinal = part2loc(this->polymer.chain); 

    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 

    for (std::vector <int> to_check: ne_list){
        if (check_avoidance(to_check, popFinal)){
            popFinal.insert(popFinal.begin(), to_check); 

            // since I have changed the position of my polymer, I am required to update it 

            this->polymer.chain = loc2part(popFinal, "polymer");  
            this->polymer.obtain_connectivity(); 
            return; 
        }
    }

    return; 
}

void Grid::ZeroToFinal(){

    int size = this->polymer.chain.size(); 
    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->polymer.chain.at(size-1).loc, this->x_len, this->y_len, this->z_len); 

    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(size-1).loc ), ne_list.end() ); 

    this->polymer.chain.erase(this->polymer.chain.begin() ); 

    std::vector <std::vector <int>> popZero = part2loc(this->polymer.chain) ; // erase first element of plocs

    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 

    for (std::vector <int> to_check: ne_list){
        if (check_avoidance(to_check, popZero)){
            popZero.push_back(to_check); 

            // since I have changed the position of my polymer, I am required to update it 
            this->polymer.chain = loc2part(popZero, "polymer"); 
            this->polymer.obtain_connectivity(); 
            return; 
        }
    }

    return; 

}

void Grid::reptation() {

    int r = (rand()%2); 
    std::cout <<"random number r is " << r <<"!"<< std::endl;


    if (r==0){
        std::cout << "Zero index is being sent to final!" << std::endl;
        this->ZeroToFinal(); 
    }    

    else{
        std::cout << "Final index is being send to zero!" << std::endl;
        this->FinalToZero();
    }

    return;

}

int Grid::CalcEn(){

    return PolymerEnergySolvation(this->polymer.chain, this->x_len, this->y_len, this->z_len, -1); 

}





// INVOLVE MONTE CARLO 

void Grid::ZeroIndexRotation_MC(){

    std::vector <Particle> cPolymer = this->polymer.chain;         // current configuration in cPolymer
    int E1 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 
    std::cout << "E1 is " << E1 << std::endl;
    // get this particles neighborlist 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(this->polymer.chain.at(1).loc, this->x_len, this->y_len, this->z_len); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(0).loc ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(2).loc ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the self-avoidance criterion from misc.cpp 
    
    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 
    for (std::vector <int> to_rot:  ne_list ) {
        if (check_avoidance( to_rot, part2loc(cPolymer)) ){
            
            cPolymer.at(0).loc = to_rot;
            int E2 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); // -1 is interaction energy

            if ( acceptance( E2-E1, 1)){
                this->polymer.chain.at(0).loc = to_rot;
                this->polymer.obtain_connectivity(); 
                return;
            }
            else {
                std::cout << "Monte Carlo said no..." << std::endl;
                // return;
            }

        }
        else {
            std::cout << "occupied..." << std::endl; 
        }
    }

    return; 
}



// 


void Grid::FinalIndexRotation_MC(){

    int dop = polymer.chain.size(); 

    std::vector <Particle> cPolymer = this->polymer.chain;         // current configuration of polymer 
    int E1 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 
    std::cout << "E1 is " << E1 << std::endl;
    // get this particles neighborlist 
    std::vector <std::vector <int> > ne_list = obtain_ne_list(this->polymer.chain.at(dop-2).loc, this->x_len, this->y_len, this->z_len); 

    // get the locations of particles it is connected to, and then erase them from the neighbor list 
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(dop-1).loc ), ne_list.end() );
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(dop-3).loc ), ne_list.end() );    

    // find a location that is unoccupied 
    // use the self-avoidance criterion from misc.cpp 
    

    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 
    for (std::vector <int> to_rot:  ne_list ) {
        if (check_avoidance( to_rot, part2loc(cPolymer)) ){

            cPolymer.at(dop-1).loc = to_rot;
            int E2 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 

            if (acceptance (E2-E1, 1)){

                this->polymer.chain.at(dop-1).loc = to_rot;
                this->polymer.obtain_connectivity(); 
                return; 
            }
            else {
                std::cout << "Monte Carlo said no..." << std::endl; 
                // 
            }
        }
        else{
            std::cout << "occupied..." << std::endl; 
        }
    
    }
    return ; 
}


void Grid::end_rotation_MC(){
    // shitty rng being used 

    int r = (rand()%2); 
    std::cout <<"random number r is " << r <<"!"<< std::endl;


    if (r==0){
        std::cout << "Zero Index is being rotated!" << std::endl;
        this->ZeroIndexRotation_MC(); 
    }    

    else{
        std::cout << "Final index is being rotated!" << std::endl;
        this->FinalIndexRotation_MC();
    }

    return; 
}
//

//

//


void Grid::FinalToZero_MC(){
    
    std::vector <Particle> cPolymer = this->polymer.chain;         // current configuration in cPolymer

    // this is the critical energy calculation move 

    int E1 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 
    std::cout << "E1 is " << E1 << std::endl;

    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->polymer.chain.at(0).loc, this->x_len, this->y_len, this->z_len); 

    // neighbor list of everything around particle at location 0, except particle at location 1
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(1).loc ), ne_list.end() );

    cPolymer.pop_back(); // pop the final element in the vector 

    std::vector <std::vector <int>> popFinal = part2loc(cPolymer); 
    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 

    for (std::vector <int> to_check: ne_list){
        if (check_avoidance(to_check, popFinal)){
            popFinal.insert(popFinal.begin(), to_check); 
            cPolymer = loc2part(popFinal, "polymer"); 
            int E2 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 

            if (acceptance((E2-E1), 1)){
            this->polymer.chain = loc2part(popFinal, "polymer");  
            this->polymer.obtain_connectivity(); 
            return;                
            }
            else{
                std::cout << "Monte Carlo said no..." << std::endl;
            }
            popFinal.erase(popFinal.begin()); // get rid of thay element you added in the beginning

           
        }
        else {
            std::cout << "occupied..." << std::endl;
        }
    }

    return; 
}


void Grid::ZeroToFinal_MC(){

    std::vector <Particle> cPolymer = this->polymer.chain; 
    int size = this->polymer.chain.size(); 

    // critical energy calculation move 

    int E1 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 
    std::cout << "E1 is " << E1 << std::endl;

    std::vector <std::vector <int>> ne_list = obtain_ne_list(this->polymer.chain.at(size-1).loc, this->x_len, this->y_len, this->z_len); 

    // neighbor list of everything around particle at location (final) except particle at location (final-1)
    ne_list.erase(std::remove(ne_list.begin(), ne_list.end(), this->polymer.chain.at(size-1).loc ), ne_list.end() ); 

    cPolymer.erase(cPolymer.begin() ); // pop the first element in the vector 

    std::vector <std::vector <int>> popZero = part2loc(cPolymer) ; // erase first element of plocs
    std::shuffle(std::begin(ne_list), std::end(ne_list), std::default_random_engine() ); 

    for (std::vector <int> to_check: ne_list){
        if (check_avoidance(to_check, popZero)){
            popZero.push_back(to_check); 
            cPolymer = loc2part(popZero, "polymer"); 
            int E2 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 

            if (acceptance((E2-E1), 1)){
            // since I have changed the position of my polymer, I am required to update it 
            this->polymer.chain = loc2part(popZero, "polymer"); 
            this->polymer.obtain_connectivity(); 
            return; 
            }
            else {
                std::cout << "Monte Carlo said no..." << std::endl;
            }
            popZero.pop_back(); // pop the final element you just added


        }
        else {
            std::cout << "occupied..." << std::endl;
        }
    }

    return; 

}


void Grid::reptation_MC() {

    int r = (rand()%2); 
    std::cout <<"random number r is " << r <<"!"<< std::endl;


    if (r==0){
        std::cout << "Zero index is being sent to final!" << std::endl;
        this->ZeroToFinal_MC(); 
    }    

    else{
        std::cout << "Final index is being send to zero!" << std::endl;
        this->FinalToZero_MC();
    }

    return;

}



// START OF KINK JUMP 
void Grid::kink_jump_MC() {

    std::vector <Particle> cPolymer = this->polymer.chain; 
    int E1 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 
    std::cout << "E1 is " << E1 << std::endl;    


    std::vector <int> k_idx = this->polymer.find_kinks();

    std::shuffle(std::begin(k_idx), std::end(k_idx), std::default_random_engine() ); 

    for (int idx: k_idx){
        std::cout << "idx right before kink spot is " << idx << std::endl;
        std::vector <int> d1 = subtract_vectors(&(this->polymer.chain.at(idx+1).loc), &(this->polymer.chain.at(idx).loc));
        std::vector <int> d2 = subtract_vectors(&(this->polymer.chain.at(idx+2).loc), &(this->polymer.chain.at(idx+1).loc));
        std::vector <int> to_check = add_vectors(&(this->polymer.chain.at(idx).loc), &d2); 
        impose_pbc(&to_check, this->x_len, this->y_len, this->z_len);
        if (check_avoidance(to_check, part2loc(this->polymer.chain)) ) {

            cPolymer.at(idx+1).loc = to_check; 
            int E2 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 

            if (acceptance(E2-E1, 1)){
            this->polymer.chain.at(idx+1).loc = to_check; 
            this->polymer.obtain_connectivity();
            return; 
            }
            else {
                std::cout << "Monte Carlo said no..." << std::endl;
            }
            cPolymer = this->polymer.chain; // reset cPolymer 
        } 
        else {
            std::cout << "occupied..." << std::endl;
        }
    }

    return;

}



// start of crank shaft 
void Grid::crank_shaft_MC() {

    std::vector <Particle> cPolymer = this->polymer.chain; 
    int E1 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 
    std::cout << "E1 is " << E1 << std::endl;

    std::vector <int> c_idx = this->polymer.find_cranks(); 

    if (c_idx.size()==0){
        std::cout << "no cranks in structure..." << std::endl;
        return;
    }

    std::shuffle(std::begin(c_idx), std::end(c_idx), std::default_random_engine() ); 

    for (int idx: c_idx){
        std::vector <int> d1 = subtract_vectors(&(this->polymer.chain.at(idx+3).loc), &(this->polymer.chain.at(idx+2).loc)); 
        impose_pbc(&d1, this->x_len, this->y_len, this->z_len); 

        std::vector <int> to_check_1 = add_vectors(&(this->polymer.chain.at(idx).loc), &d1);
        std::vector <int> to_check_2 = add_vectors(&(this->polymer.chain.at(idx+3).loc), &d1); 

        impose_pbc(&to_check_1, this->x_len, this->y_len, this->z_len);
        impose_pbc(&to_check_2, this->x_len, this->y_len, this->z_len);

        if ((check_avoidance(to_check_1, part2loc(this->polymer.chain)) ) && (check_avoidance(to_check_2, part2loc(this->polymer.chain)) ) ) {

            // reupdate polymer now that I have changed it 
            cPolymer.at(idx+1).loc = to_check_1; 
            cPolymer.at(idx+2).loc = to_check_1; 

            int E2 = PolymerEnergySolvation(cPolymer, this->x_len, this->y_len, this->z_len, -1); 

            if (acceptance(E2-E1, 1)){
            this->polymer.chain.at(idx+1).loc = to_check_1;
            this->polymer.chain.at(idx+2).loc = to_check_2; 
            this->polymer.obtain_connectivity(); 
            return;
            }
            else {
                std::cout << "Monte Carlo said no..." << std::endl;
            }
            cPolymer = this->polymer.chain; 
        }
        else{
            std::cout << "occupied..." << std::endl;
        }

    }    
    return ;
}


