#include "Simulation.h"

//////////////////////////////////////////////////////////
//
//                Set-up the Lattice
//
//////////////////////////////////////////////////////////

void Simulation::set_up_energy_calculator(){

	if (this->potts) {
		this->calculate_energy_ptr = &Simulation::accelerate_calculate_energy_potts;
	}
	else {
		if (this->frac_c > 0.5){
			std::cout << "We are running the energy computation by looping over solvent particles." << std::endl;
			this->calculate_energy_ptr = &Simulation::accelerate_calculate_energy_solvent; 
		}
		else {
			std::cout << "We are running the energy computation by looping over cosolvent particles." << std::endl;
			this->calculate_energy_ptr = &Simulation::accelerate_calculate_energy_cosolvent; 
		}
	}
	return;
}

void Simulation::set_up_local_dump(){

	if (this->potts){
		this->dump_local_ptr = &Simulation::dump_potts;
	}
	else if (this->SSfile == "__blank__"){
		this->dump_local_ptr = &Simulation::dump_local_no_ss;
	}
	else {
		this->dump_local_ptr = &Simulation::dump_local_all;
	}
	return;
}

void Simulation::set_up_style_of_run(){
	if (this->potts){
		this->run_ptr = &Simulation::run_potts;
	}

	else {
		if (this->v){
			this->run_ptr = &Simulation::run_debug;
		}
		else if (this->isotropic){
			this->run_ptr = &Simulation::run_isotropic;
		}
		else {
			this->run_ptr = &Simulation::run_straight;
		}
	}
	return;
}

void Simulation::set_up_system(){

	if (this->r){
		if (this->potts){
			this->set_up_lattice_for_restart();
			this->Polymers.clear   ();
			this->Polymers.reserve (0);
			this->Solvent.clear    ();
			this->Solvent.reserve  (0);
			this->Cosolvent.clear  ();
			this->Cosolvent.reserve(0);
			filter_csv(this->efile, this->step_number);
		}
		else {
			this->set_up_lattice_for_restart();
			this->set_up_polymers_for_restart();
			this->set_up_files_for_restart();
			this->enhanced_flipper.setup(1, 5);
			this->enhanced_swing.setup(1, 5);
		}
	}
	else {
		if (this->potts){
			std::cout << "Setting up a Potts simulation." << std::endl;
			this->set_up_from_scratch_potts();
		}
		else {
			std::cout << "Setting up an FHP simulation." << std::endl;
			this->set_up_from_scratch();
		}
	}

	return;
}

//////////////////////////////////////////////

void Simulation::set_up_add_solvent(){

	int c_idx {-1};
	std::array <int,3> loc = {0,0,0}; 
	for (int k{0}; k<z; ++k){
		for (int j{0}; j<y; ++j){
			for (int i{0}; i<x; ++i){
				loc = {i, j, k};
				if (lattice_index(loc,y,z)-c_idx != 1){
					std::cerr << "Something is fucked in Lattice creation." << std::endl;
					exit (EXIT_FAILURE);
				}
				c_idx = lattice_index (loc, y, z); 
				Particle* p_ptr = new Particle ( loc, "s1", rng_uniform (0, 25) ); 
				this->Lattice.insert( this->Lattice.begin() + lattice_index(loc, y, z), p_ptr) ;
			}
		}
	}
	return; 
}

void Simulation::set_up_add_cosolvent(){

	std::cout << "\n--------------------------------------------------------------------\n" << std::endl;
	// std::cout << "Begin adding cosolvent... \n"; 
	int Nmonomer = 0;
	for (Polymer& pmer: this->Polymers){
		Nmonomer += pmer.chain.size();
	}
	int nsol2         = std::floor ((this->x*this->y*this->z-Nmonomer)*this->frac_c); 
	std::cout << "Number of particles of cosolvent being added is " << nsol2 << "." << std::endl;
	std::vector <int> indices (this->x*this->y*this->z);
	std::iota (indices.begin(), indices.end(), 0); 
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	// get indices to loop through 
	std::shuffle ( indices.begin(), indices.end(), std::default_random_engine(seed) );
	int count = 0;
	int i = 0; 

	// find solvation shell indices 
	if (this->s) {
		std::set <int> solvation_shell_set;
		std::array <std::array<int,3>, 26> ne_list; 
		for ( Polymer& pmer: this->Polymers){
			for ( Particle*& p: pmer.chain ){
				ne_list = obtain_ne_list ( p->coords, this->x, this->y, this->z); 
				for ( std::array <int,3>& loc: ne_list ){
					if ( (this->Lattice).at(lattice_index(loc, y, z))->ptype[0] == 's' ){
						solvation_shell_set.insert(lattice_index (loc, this->y, this->z)); 
					}
				}
			}
		}
		std::vector <int> solvation_shell_indices(solvation_shell_set.begin(), solvation_shell_set.end());
		
		while ( (count < nsol2) && (count < static_cast<int>(solvation_shell_indices.size()))) {

			this->Lattice.at(solvation_shell_indices.at(i))->ptype = "s2"; 
			this->Cosolvent.push_back ( this->Lattice.at(solvation_shell_indices.at(i)));
			count += 1;
			i     += 1;
		}
	}

	i = 0; 
	while ( count < nsol2 ) {

		if ( this->Lattice.at( indices[i] )->ptype[0] == 'm' ){
			; 
		}
		else if ( this->Lattice.at(indices[i])->ptype == "s2" ) {
			;
		}
		else {
			this->Lattice.at( indices[i] )->ptype = "s2"; 
			this->Cosolvent.push_back((this->Lattice).at( indices[i] )); 
			count += 1; 
		}
		i += 1; 
	}

	// now that we have all the cosolvent on the lattice, it is now time to fill up the Solvent vector
	int nsol1 = 0;
	for (int i{0}; i < static_cast<int>(this->Lattice.size()); ++i){
		if ((this->Lattice[i])->ptype == "s1"){
			this->Solvent.push_back(this->Lattice[i]);
			nsol1 += 1;
		}
	}

	std::cout << "Number of solvent particles in the system is " << nsol1 << "." << std::endl;

	return; 
}

void Simulation::set_up_align_lattice(){
	for (Particle*& p: this->Lattice) {
		p->orientation = 0;
	}

	return; 
}

void Simulation::set_up_align_solvation_shell(){

	std::array  <std::array<int,3>, 26> ne_list;
	std::vector <int> solvent_indices; 
	
	for ( Polymer& pmer: this->Polymers ){
		for ( Particle*& p: pmer.chain ) {
			
			ne_list = obtain_ne_list (p->coords, this->x, y, z);
			for ( std::array<int,3>& ne: ne_list ){
				if ( this->Lattice.at(lattice_index(ne, y, z))->ptype[0] == 's' ) {
					solvent_indices.push_back (lattice_index(ne, y, z));
				}
			}
		}
	}  

	// get rid of duplicates 
	std::unordered_set <int> s (solvent_indices.begin(), solvent_indices.end() );
	solvent_indices.assign (s.begin(), s.end() ); 

	for ( Polymer& pmer: this->Polymers ) {
		for ( Particle*& p: pmer.chain ) {
			p->orientation = 0;
		}
	}

	for (int i: solvent_indices){
		this->Lattice[i]->orientation = 0;
	}

	return;
}

//////////////////////////////////////////////
//
// RESTART SET UP
//
//////////////////////////////////////////////

void Simulation::set_up_lattice_for_restart(){

	std::vector <Particle*> Lattice;
	std::vector <Particle*> Solvent;
	std::vector <Particle*> Cosolvent;

	Lattice.reserve(this->x*this->y*this->z);

	std::vector<std::string> contents = this->extract_content_from_file(this->lattice_file_read);
	std::regex start ("FINAL STEP: "), end ("END");
	std::regex numbers ("[0-9]+");
	std::regex characters ("[a-z][0-9]+");
	std::regex pattern (R"([^,]+)");
	std::smatch mat;

	int final_step = 0;
	bool start_recording = false;

	std::string part1;
	std::string part2;
	std::string part3;
	

	std::cout << "About to dive into the lattice file..." << std::endl;
	Particle* p_ptr {nullptr};

	// find the last step dumped in lattice file
	for ( std::string& s: contents) {
		if ( std::regex_search (s, start) ) {
			std::regex_search (s, mat, numbers);
			final_step = std::stoi (mat[0].str());
		}
	}

	this->step_number = final_step;

	for ( std::string& s: contents) {
		// std::cout << "line = " << s << std::endl;
		if ( std::regex_search (s, start) ) {
			std::regex_search ( s, mat, numbers ); 
			if ( final_step == std::stoi(mat[0].str()) ) {
				start_recording = true;
			}
		}

		else if ( std::regex_search (s, end) && start_recording ){
			break; 
		}

		else if (start_recording) {
			// std::cout << s << std::endl;
			std::smatch matches;
			std::string::const_iterator search_start (s.begin());
			std::vector <std::string> parts; 
			while (std::regex_search(search_start, s.cend(), matches, pattern)) {
				for (auto& match : matches) {
					parts.push_back(match.str());
				}
				search_start = matches.suffix().first;
			}
			parts[1] = trim(parts[1]);
			p_ptr = new Particle (location(std::stoi(parts[2]), this->x, this->y, this->z), parts[1], std::stoi(parts[0])); 

			Lattice.insert (Lattice.begin() + std::stoi(parts[2]), p_ptr);
			if (parts[1] == "s1"){
				Solvent.push_back(p_ptr);
			}
			else if (parts[1] == "s2"){
				Cosolvent.push_back(p_ptr);
			}
			else {
				;
				// std::cout << "Ptype is " << parts[1] << "." << std::endl;
			}
		}
	}

	if (!start_recording) {
		std::cout <<"Something is wrong with the restart!" << std::endl;
		std::cerr << "Exiting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	else {
		std::cout << "Created lattice from file!" << std::endl;
		this->Lattice   = Lattice;
		this->Solvent   = Solvent;
		this->Cosolvent = Cosolvent;
	}
	return;
}

void Simulation::set_up_polymers_for_restart(){

	std::vector <Polymer> polymer_container; 
	std::vector <std::array<int,3>> locations; 
	std::vector <int> spins;

	std::vector <std::string> contents = this->extract_content_from_file(this->dfile);

	int  final_step_num = this->step_number;
	bool step_bool {false}, start_bool {false}, end_bool {false}; 

	std::regex start ("START"), end ("END"), step ("step " + std::to_string(final_step_num) );
	std::regex step_generic ("step"); 
	std::regex reg_poly ("Dumping coordinates of Polymer"); 
	std::regex numbers ("[0-9]+"); 

	int startCount{0}, endCount{0}; 
	std::array <int,3> loc;

	std::smatch match; 

	for (std::string& s: contents){
		
		if ( std::regex_search (s, step_generic) ){

			std::regex_search ( s, match, numbers ); 
			std::regex_token_iterator<std::string::iterator> rend; 
			std::regex_token_iterator<std::string::iterator> a (s.begin(), s.end(), numbers);

			if ( std::stoi(*a) > final_step_num ){
				std::cerr << "\n\nYou coordinates file and restart file are not in sync. \nIt is probably because the restart file you chose is not for that simulation.\nPlease check them out. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
		}

		if ( std::regex_search (s, step) ) {
			step_bool = true; 
			continue; 
		}

		if (step_bool){
			if ( std::regex_search(s, start) ){
				++startCount;
				start_bool = true; 
				end_bool   = false; 
				continue; 
			}

			else if (std::regex_search(s, end) ) {
				++endCount;
				start_bool = false; 
				end_bool   = false; 
				step_bool  = false; 

				// std::cout << "Final line = " << s << std::endl;
				Polymer pmer = Polymer (locations, spins);
				for (int i=0; i<pmer.deg_poly; ++i){
					pmer.chain[i] = this->Lattice[lattice_index(pmer.chain[i]->coords, this->y, this->z)];
				}
				polymer_container.push_back(pmer);
				
				locations.clear();
				break;
				
			}

			else if (start_bool == end_bool){
				continue;
			}

			else{
				// std::cout << s << std::endl;
				std::regex_search ( s, match, numbers ); 
				std::regex_token_iterator<std::string::iterator> rend; 
				std::regex_token_iterator<std::string::iterator> a ( s.begin(), s.end(), numbers );

				for ( int i=0; i<4; ++i ){
					if ( i ==0 ) {
						// std::cout << "x-coords is " << *a << std::endl; 
						loc[i] = std::stoi ( *a );
						*a++; 
					}
					else if ( i==1 ){
						// std::cout << "y-coord is " << *a << std::endl; 
						loc[i] = std::stoi ( *a );
						*a++; 
					}
					else if ( i==2 ){
						// std::cout << "z-coord is " << *a << std::endl; 
						loc[i] = std::stoi ( *a );
						*a++;
					}
					else if ( i==3 ){
						spins.push_back ( std::stoi(*a) ); 
					}
				}  
				
				if (!this->check_validity_of_coords(loc)){
					std::cerr << "Coordinates are out of bounds. Bad input file." << std::endl;
					exit(EXIT_FAILURE); 
				}
				// std::cout << "location = "; print(loc);
				locations.push_back(loc); 
				
			}
		}
	}

	this->Polymers = polymer_container;
	return;
}

void Simulation::set_up_files_for_restart(){

	filter_coords(this->dfile,  this->step_number);
	filter_csv   (this->efile,  this->step_number);
	filter_csv   (this->SSfile, this->step_number);
	filter_orientations(this->mfile, this->step_number);

	return;
}

void Simulation::set_up_for_restart(){

	if (this->potts){
		this->set_up_lattice_for_restart();
		this->Polymers.clear   ();
		this->Polymers.reserve (0);
		this->Solvent.clear    ();
		this->Solvent.reserve  (0);
		this->Cosolvent.clear  ();
		this->Cosolvent.reserve(0);
		filter_csv(this->efile, this->step_number);
	}
	else {
		this->set_up_lattice_for_restart();
		this->set_up_polymers_for_restart();
		this->set_up_files_for_restart();
		this->enhanced_flipper.setup(1, 5);
		this->enhanced_swing.setup(1, 5);
	}
	return;
}

//////////////////////////////////////////////

void Simulation::set_up_from_scratch(){

	// set up the system from scratch
	this->step_number = 0;

	// extract the number of polymers from the file
	this->extract_number_of_polymers();

	// initialize custom data structures
	// this data structure will hold the coordinates of the polymer
	std::vector <Polymer> Polymers;
	Polymers.reserve(this->Npoly);

	// this data structure will hold the coordinates of the cosolvent
	std::vector <Particle*> Cosolvent;
	std::vector <Particle*> Solvent;

	// this data structure will hold the lattice
	std::vector <Particle*> Lattice; 
	Lattice.reserve((this->x)*(this->y)*(this->z));

	this->Lattice   =  Lattice;
	this->Polymers  =  Polymers;
	this->Solvent   =  Solvent;
	this->Cosolvent =  Cosolvent;

	this->extract_polymers_from_file();

	// just raw addition of solvent particles
	// these particles will be overwritten if frac_c > 0
	this->set_up_add_solvent();

	// populate the lattice
	for (Polymer& pmer: (this->Polymers)) {
		for (Particle*& p: pmer.chain){
			// now that I have my polymer coordinates, time to create the grand lattice 
			(this->Lattice).at(lattice_index (p->coords, y, z) ) = p; 
		}
	}

	// throw the cosolvents in
	this->set_up_add_cosolvent (); 

	if (this->A){
		this->set_up_align_lattice();
	}
	else if (this->S){
		this->set_up_align_solvation_shell();
	}
	
	this->check_structures();
	this->enhanced_flipper.setup(1, 5);
	this->enhanced_swing.setup(1, 5);
	return;
}

//////////////////////////////////////////////

void Simulation::set_up_from_scratch_potts(){

	this->Npoly = 0;
	std::vector <Particle*> Cosolvent;
	std::vector <Particle*> Solvent  ;
	std::vector <Polymer>   Polymers ;
	Polymers.reserve (0);
	Cosolvent.reserve(0);
	Solvent.reserve  (0);

	// create the Lattice data structure
	std::vector <Particle*> Lattice;
	Lattice.reserve((this->x)*(this->y)*(this->z));

	// update the object
	this->Lattice   = Lattice;
	this->Solvent   = Solvent;
	this->Cosolvent = Cosolvent;
	this->Polymers  = Polymers; 

	// populate the lattice 
	int c_idx {-1}; 
	std::array <int,3> loc = {0, 0, 0};
	for (int k{0}; k<this->z; ++k){
		for (int j{0}; j<this->y; ++j){
			for (int i{0}; i<this->x; ++i){
				loc = {i, j, k};
				if (lattice_index(loc, this->y, this->z) - c_idx != 1){
					std::cerr << "Something is fucked in Lattice creation." << std::endl;
					exit (EXIT_FAILURE);
				}
				c_idx = lattice_index(loc, this->y, this->z);
				Particle* p_ptr = new Particle (loc, "sp", rng_uniform(0, 25));
				this->Lattice.insert(this->Lattice.begin() + lattice_index(loc, this->y, this->z), p_ptr);
			}
		}
	}

	return;

}

//////////////////////////////////////////////

void Simulation::set_up_FHP(){

	this->extract_topology_from_file();                                  // I have the geometry and energies.
	this->set_up_system();                                               // now that i have all the info, i can set up the simulation lattice
	this->set_up_local_dump();                                           // sets up the dump function
	this->set_up_energy_calculator();                                    // sets up the energy function
	this->set_up_style_of_run();                                         // sets up the run function
	this->initialize_pairwise_function_map();                            // initialize the pairwise function
	this->initialize_neighbor_function_map();                            // initialize the neighbor function
	this->accelerate_calculate_energy();                                 // get the energy of the system
	this->dump_local();                                                  // dump out the conditions at step number 0
	this->debug_checks_energy_contacts(this->sysEnergy, this->contacts); // run the final debugging check
	this->check_structures();                                            // run another structure check

	return;

}

void Simulation::set_up_Potts(){

	this->extract_topology_for_potts();       // get the geometry, energies, and composition of the system.
	this->set_up_system();                    // now that i have the topology, I can now set up the system.
	this->set_up_local_dump();                // sets up the dump function
	this->set_up_energy_calculator();         // set up the energy calculator
	this->set_up_style_of_run();              // set up the run function
	this->initialize_pairwise_function_map(); // initialize the pairwise function
	this->initialize_neighbor_function_map(); // initialize the neighbor function
	this->accelerate_calculate_energy();      // get the energy of the system
	this->dump_local();                       // dump out the conditions at step number 0

	return;

}
