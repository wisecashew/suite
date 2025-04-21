#include "lattice_directions.hpp"
#include "Simulation.hpp"

//////////////////////////////////////////////////////////
//   
//                Verify structures
//
//////////////////////////////////////////////////////////

bool Simulation::check_validity_of_coords(std::array<int,3> v){
	if (v.at(0)>this->x || v.at(0) < 0 || v.at(1)>this->y || v.at(1)<0 || v.at(2)>this->z || v.at(2)<0){
		return false;
	}
	else {
		return true;
	}
}

bool Simulation::check_for_overlaps_within_polymers_raw(){

	std::vector <std::array <int,3>> loc_list; 
	for (Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
			// print(p->coords);
			// check if element exists in vector 
			if (std::find(loc_list.begin(), loc_list.end(), p->coords) != loc_list.end() ){
				std::cerr << "you have a repeated element." << std::endl;
				print(p->coords); 
				print(loc_list);
				return false; 
			}
			else{
				loc_list.push_back(p->coords);  
			}
		}
	}
	std::cout << "Input file has no overlaps! ";
	return true;
}

bool Simulation::check_for_overlaps_within_polymers(){

	std::vector <std::array <int,3>> loc_list; 

	for (Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
		// check if element exists in vector 
			if (std::find(loc_list.begin(), loc_list.end(), p->coords) != loc_list.end() ){
				std::cerr << "you have a repeated element." << std::endl;
				// std::cout << "current element is: " << std::endl;
				print(p->coords); 
				print(loc_list);
				return false; 
			}
			else{
				loc_list.push_back(p->coords);  
			}
		}
	}
	std::cout << "Input file has no overlaps! ";
	return true;
}

bool Simulation::check_for_solvent_monomer_overlap(){

	std::vector <std::array <int,3>> loc_list;

	for (Polymer& pmer: this->Polymers) {
		for (Particle*& p: pmer.chain){
			// check if element exists in vector 
				if (std::find(loc_list.begin(), loc_list.end(), p->coords) != loc_list.end() ){
					std::cerr << "you have a repeated element." << std::endl;
					return false; 
				}
			
				else{
					loc_list.push_back(p->coords);  
				}
		}
	}    
	
	for (Polymer& pmer: this->Polymers) {
		for (Particle*& p: pmer.chain){

			if ( this->Lattice[ lattice_index (p->coords, y, z)]->coords == p->coords ) {

				if ( this->Lattice[ lattice_index(p->coords, y, z)]->ptype[0] == 's' ){
					std::cerr << "Some kind of bad solvent-monomer overlap that has taken a place. A monomer is being represented by a solvent. Something's fucked." << std::endl;
					std::cerr << "Location is: "; print (p->coords); 
					std::cerr << "Type is: " << (this->Lattice[ lattice_index (p->coords, y, z)]->ptype) << std::endl; 
					std::cerr << "Type is: " << p->ptype << std::endl; 
					std::cerr << "Location on lattice is: "; print (this->Lattice[ lattice_index (p->coords, y, z)]->coords); 
					exit (EXIT_FAILURE); 
				}
				else {
					continue;
				}

			}

			else {
				std::cerr << "Something is wrong with the Lattice map. Monomer unit lost. Something is fucked. " << std::endl;
				std::cerr << "Output from this->Polymers: "; print (p->coords);
				std::cerr << "Output from this->Lattice: "; print (this->Lattice[ lattice_index (p->coords, y, z)]->coords);
				exit(EXIT_FAILURE);
			}

		}
	}

	std::cout << "Input file has no overlap between and internally amongst solvent and monomers!" << std::endl;
	return true;
}

bool Simulation::check_for_overlaps_on_lattice(){

	bool pmer_loop_flag = false; 
	bool particle_found_flag = false; 
	for ( int i{0}; i< static_cast<int>(this->Lattice.size()); ++i ){
		
		pmer_loop_flag = false; 
		particle_found_flag = false; 

		if ( this->Lattice[i]->ptype[0] == 'm' ){

			for ( Polymer& pmer: this->Polymers ){
				for (Particle*& p: pmer.chain){
					if (this->Lattice[i]->coords == p->coords){
						pmer_loop_flag = true; 
						break;
					}
				}
				if (pmer_loop_flag) {
					particle_found_flag = true; 
					break;
				}
			}

			if (particle_found_flag){
				continue;
			}
			else {
				std::cout << "Something is fucked. There is a random monomer floating around." << std::endl;
				std::cout << "Location = "; print (this->Lattice[i]->coords);
				exit (EXIT_FAILURE); 
			}
		}

	}

	return true;
}

bool Simulation::check_connectivity_raw(){

	for (Polymer& pmer: this->Polymers){
		size_t length = pmer.chain.size(); 
		std::array <int,3> connection = {0,0,0}; 
		std::sort (adrns.begin(), adrns.end() ); 

		for (int i{1}; i<static_cast<int>(length); ++i){
			
			connection = subtract_containers((pmer.chain[i]->coords), (pmer.chain[i-1]->coords));
			impose_pbc(&connection, this->x, this->y, this->z);
			modified_direction ( &connection, this->x, this->y, this->z); 

			if ( binary_search ( adrns.begin(), adrns.end(), connection) ) {
				continue;
			}
			else {
				std::cerr << "Shit, you have bad connectivity inside one (or maybe more) polymers. Check input file." << std::endl;
				return false; 
			}

		}
	}

	std::cout << "Input polymers are well-connected! \n";
	return true;
}

bool Simulation::check_connectivity(){

	for (Polymer& pmer: this->Polymers){
		size_t length = pmer.chain.size(); 
		std::array <int,3> connection = {0,0,0}; 
		std::sort (adrns.begin(), adrns.end() ); 

		for (int i{1}; i<static_cast<int>(length); ++i){
			
			connection = subtract_containers((pmer.chain[i]->coords), (pmer.chain[i-1]->coords));
			impose_pbc(&connection, this->x, this->y, this->z);
			modified_direction ( &connection, this->x, this->y, this->z); 

			if ( binary_search ( adrns.begin(), adrns.end(), connection) ) {
				continue;
			}
			else {
				std::cerr << "Shit, you have bad connectivity inside one (or maybe more) polymers. Check input file." << std::endl;
				return false; 
			}

		}
	}

	std::cout << "Input polymers are well-connected! \n";
	return true;
}

bool Simulation::check_pointers_on_lattice(){
	for ( int i{0}; i<this->x*this->y*this->z; ++i ) {

		if ( location(i, this->x, this->y, this->z) == (this->Lattice)[i]->coords ) {
			continue;
		}
		else {
			std::cout << "Problem with pointers." << std::endl;
			std::cout << "Bad location is "; print ( location (i, x, y, z) ); 
			std::cout << "Lattice says "; print ((this->Lattice)[i]->coords);
			std::cerr << "Something is fucked. Pointer does not correspond to position. " << std::endl;
			exit (EXIT_FAILURE);
			return false;
		}

	}

	return true;
}

void Simulation::check_structures(){

	// std::cout << "Checking validity of coords...";
	// std::cout << std::boolalpha << "checkForOverlaps says: " << this->check_for_overlaps_within_polymers() << "." << std::endl; 
	if (!this->check_for_overlaps_within_polymers()){
		std::cout << "Something is fucked up overlaps-wise in the polymer itself." << std::endl;
		exit(EXIT_FAILURE);
	}
	// std::cout << "Looks good." << std::endl;

	// std::cout << "Checking if no stray monomers are on the Lattice...";
	if (!this->check_for_overlaps_on_lattice()){
		std::cout << "Random monomer floating in the lattice! Breaking out..." << std::endl; 
		exit (EXIT_FAILURE);
	}
	// std::cout << "No random monomers floating around!" << std::endl;

	if (!this->check_for_solvent_monomer_overlap()){
		std::cout << "Something is fucked up solvent-monomer overlaps-wise. " << std::endl; 
		exit(EXIT_FAILURE);
	}

	// std::cout << std::boolalpha << "checkConnectivity says: " << this->check_connectivity() << "." << std::endl;
	if (!this->check_connectivity()){
		std::cout << "Something is fucked up connectivity-wise." << std::endl; 
		this->Polymers[0].print_chain();
		exit(EXIT_FAILURE);
	}

	for (Particle*& p: this->Cosolvent){

		if ( this->Lattice.at( lattice_index(p->coords, this->y, this->z) )->ptype != "s2" && this->Lattice.at( lattice_index(p->coords, this->y, this->z) )->ptype != p->ptype && (this->Lattice).at(lattice_index(p->coords, this->y, this->z))->orientation != p->orientation && \
		this->Lattice.at( lattice_index(p->coords, this->y, this->z) )->coords != p->coords ) {
			std::cerr << "Cosolvent and Lattice do not align. Something's fucked." << std::endl;
			exit (EXIT_FAILURE);
		}

	}

	// std::cout << "Cosolvent placement looks good!" << std::endl;

	for (Particle*& p: this->Solvent){

		if ( this->Lattice.at( lattice_index(p->coords, this->y, this->z) )->ptype != "s1" && this->Lattice.at( lattice_index(p->coords, this->y, this->z) )->ptype != p->ptype && (this->Lattice).at(lattice_index(p->coords, this->y, this->z))->orientation != p->orientation && \
		this->Lattice.at( lattice_index(p->coords, this->y, this->z) )->coords != p->coords ) {
			std::cerr << "Cosolvent and Lattice do not align. Something's fucked." << std::endl;
			exit (EXIT_FAILURE);
		}
	}

	// std::cout << "Solvent placement looks good!" << std::endl;

	if (this->check_pointers_on_lattice()){
		// std::cout << "Okay. Lattice is in good shape." << std::endl; 
	}
	else {
		std::cerr <<"Something is fucked with pointers on Lattice." << std::endl; 
		exit (EXIT_FAILURE);
	}

	return; 
}
