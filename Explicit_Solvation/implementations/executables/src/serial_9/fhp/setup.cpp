#include "FHP.hpp"

void FHP::setup_parser(
	int dfreq,    int lfreq,         \
	int max_iter, bool restart_bool, \
	std::string topology,             std::string energy_dump_file,      \
	std::string coords_dump_file,     std::string orientation_dump_file, \
	std::string solvation_shell_file, std::string stats_file,            \
	std::string lattice_file_write,   std::string lattice_file_read){

	if (!restart_bool){
		if (dfreq == -1 || lfreq == -1 || max_iter == -1){
			std::cerr << "ERROR: No value for option l (frequency of dumping) (" << lfreq << ")" << " and/or for option M (maximum number of moves to be performed) (" << max_iter << ") was provided. Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
		else if (topology     == "__blank__" || \
		energy_dump_file      == "__blank__" || \
		coords_dump_file      == "__blank__" || \
		solvation_shell_file  == "__blank__" || \
		stats_file            == "__blank__" || \
		lattice_file_write    == "__blank__" || \
		orientation_dump_file == "__blank__") {
			std::cerr << 
			"Topology is in " << topology <<
			",\nenergy dump file is "      << energy_dump_file      << 
			",\ncoords dump file is "      << coords_dump_file      << 
			",\nsolvation shell file is "  << solvation_shell_file  << 
			",\nmove statistics file is "  << stats_file            << 
			",\nlattice dump file is "     << lattice_file_write    << 
			",\norientation dump file is " << orientation_dump_file << "." << std::endl;
			// std::cerr << "ERROR: No value for " <<
			// "for option t (energy and geometry file) and/or\n" <<
			// "for option s (name of move stats file) and/or\n for option u (name of energy dump file). Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
		else {
			std::ofstream out_lattice_dump_file(lattice_file_write,    std::ios::out);
			std::ofstream out_coords_file      (coords_dump_file,      std::ios::out);
			std::ofstream out_orientation_file (orientation_dump_file, std::ios::out);
			std::ofstream out_energy_dump_file (energy_dump_file,      std::ios::out);
			std::ofstream out_statistics_dump_file (stats_file,        std::ios::out);

			out_lattice_dump_file.close();
			out_coords_file.close();
			out_orientation_file.close();
			out_energy_dump_file.close();
			out_statistics_dump_file.close();
		}
	}
	else {
		if (dfreq == -1 || lfreq == -1 || max_iter == -1){
			std::cerr << "ERROR: No value for option l (frequency of dumping) (" << lfreq << ")" << " and/or for option M (maximum number of moves to be performed) (" << max_iter << ") was provided. Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
		else if (topology     == "__blank__" || \
		energy_dump_file      == "__blank__" || \
		coords_dump_file      == "__blank__" || \
		solvation_shell_file  == "__blank__" || \
		stats_file            == "__blank__" || \
		lattice_file_write    == "__blank__" || \
		lattice_file_read     == "__blank__" || \
		orientation_dump_file == "__blank__") {
			std::cerr << 
			"Topology is in " << topology <<
			",\nenergy dump file is "      << energy_dump_file      << 
			",\ncoords dump file is "      << coords_dump_file      << 
			",\nsolvation shell file is "  << solvation_shell_file  << 
			",\nmove statistics file is "  << stats_file            << 
			",\nlattice dump file is "     << lattice_file_write    << 
			",\nlattice read file is "     << lattice_file_read     << 
			",\norientation dump file is " << orientation_dump_file << "." << std::endl;
		}
	}

	return;
}

void FHP::setup_add_solvent(){

	int c_idx {-1};
	std::array <int,3> loc = {0, 0, 0};
	for (int k{0}; k<this->z; ++k){
		for (int j{0}; j<this->y; ++j){
			for (int i{0}; i<this->x; ++i){
				loc = {i, j, k};
				if (lattice_index(loc, this->y, this->z)-c_idx != 1){
					std::cerr << "Something is fucked in Lattice creation." << std::endl;
					exit(EXIT_FAILURE);
				}
				else {
					c_idx = lattice_index(loc, this->y, this->z);
					Particle* p_ptr = new Particle (rng_uniform(0, 25), "s1", loc);
					this->Lattice.insert(this->Lattice.begin() + lattice_index(loc, this->y, this->z), p_ptr);
				}
			}
		}
	}
	return;
}

void FHP::setup_add_cosolvent(){

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

void FHP::setup_align_lattice(){
	for (Particle*& p: this->Lattice) {
		p->orientation = 0;
	}
	return;
}

void FHP::setup_align_solvation_shell(){

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

void FHP::setup_from_scratch(){

	// set up the system from scratch
	this->step_number = 0;

	// extract the number of polymers from the file
	// std::cout << "Extract number of polymers." << std::endl;
	this->extract_number_of_polymers();

	// initialize the custom data structures
	// this data structure will hold the coordinates of the polymer
	std::vector <Polymer> Polymers;
	Polymers.reserve(this->N_poly);

	// this data structure will hold the coordinates of the cosolvent
	std::vector <Particle*> Cosolvent;
	std::vector <Particle*> Solvent;

	// this data structure will hold the lattice
	std::vector <Particle*> Lattice;
	Lattice.reserve(this->x * this->y * this->z);

	this->Lattice   = Lattice;
	this->Polymers  = Polymers;
	this->Solvent   = Solvent;
	this->Cosolvent = Cosolvent;

	this->extract_polymers_from_file();

	// just raw addition of solvent particles
	// these particles will be overwritten if frac_c > 0
	this->setup_add_solvent();

	// populate the lattice
	for (Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
			(this->Lattice).at(lattice_index(p->coords, this->y, this->z)) = p;
		}
	}

	// throw the cosolvent in
	this->setup_add_cosolvent();

	if (this->A){
		this->setup_align_lattice();
	}
	else if (this->S){
		this->setup_align_solvation_shell();
	}
	else {
		// do nothing
	}

	this->check_structures();
	this->enhanced_flipper.setup(1, 5);
	this->enhanced_swing.setup(1, 5);
	return;

}

void FHP::setup_lattice_from_restart(){

	// set up the lattice object    
	std::vector <Particle*> Lattice;
	std::vector <Particle*> Solvent;
	std::vector <Particle*> Cosolvent;

	Lattice.reserve((this->x)*(this->y)*(this->z));

	// get the content from the file
	std::vector<std::string> contents = extract_content_from_file(this->inp_lattice_file_read);

	// set up the regexs so that you can start jumping in 
	std::regex start ("FINAL STEP: "), end ("END");
	std::regex numbers ("[0-9]+");
	std::regex characters ("[a-z][0-9]+");
	std::regex pattern (R"([^,]+)");
	std::smatch mat;

	// some more variable instantiations...
	int final_step       = 0;
	bool start_recording = false;
	std::string part1;
	std::string part2;
	std::string part3;
	
	// begin
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
			p_ptr = new Particle (std::stoi(parts[0]), parts[1], location(std::stoi(parts[2]), this->x, this->y, this->z)); 

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

void FHP::setup_polymers_for_restart(){

	std::vector <Polymer> polymer_container;
	std::vector <std::array<int,3>> locations;
	std::vector <int> spins;

	std::vector <std::string> contents = this->extract_content_from_file(this->out_coord_dump);

	int final_step_num = this->step_number;
	bool step_bool {false}, start_bool {false}, end_bool{false};

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
				locations.push_back(loc);

			}
		}
	}

	this->Polymers = polymer_container;
	return;

}

void FHP::setup_lattice(){

	// to restart simulation
	if (this->r){
		this->setup_lattice_from_restart();
		this->setup_polymers_for_restart();
		filter_coords(this->out_coord_dump, this->step_number);
		filter_csv(this->out_solvation_shell_dump, this->step_number);
		filter_csv(this->out_energy_dump, this->step_number);
		filter_orientations(this->out_orientation_dump, this->step_number);
		this->enhanced_flipper.setup(1, 5);
		this->enhanced_swing.setup(1, 5);
	}
	// set up the lattice
	else {
		this->setup_from_scratch();
	}
	return;
}

void FHP::setup_energy_calculator(){
	
	if (this->dry && this->field){
		std::cout << "We are running a dry simulation with field." << std::endl;
		this->calculate_energy_ptr = &FHP::calculate_energy_dryfield;
	}
	else if (this->frac_c > 0.5){
		std::cout << "We are running the energy computation by looping over solvent particles." << std::endl;
		this->calculate_energy_ptr = &FHP::calculate_energy_solvent; 
	}
	else {
		std::cout << "We are running the energy computation by looping over cosolvent particles." << std::endl;
		this->calculate_energy_ptr = &FHP::calculate_energy_cosolvent; 
	}
	return;
}

void FHP::setup_energetics(){
	// std::cout << "Initialize initial function map." << std::endl;
	this->initialize_interaction_function_map();
	// std::cout << "Initialize neighbor function map." << std::endl;
	this->initialize_neighbor_function_map();
	// std::cout << "Set up the energy calculator." << std::endl;
	this->setup_energy_calculator(); 
	// std::cout << "Compute energy." << std::endl;
	this->energy_compute();
	// std::cout << "Ran the compute." << std::endl;
	return;
}

void FHP::setup_dump_stats(){
	if (this->isotropic){
		this->dump_stats_ptr = &FHP::dump_stats_isotropic;
	}
	else{
		this->dump_stats_ptr = &FHP::dump_stats_all;
	}
	return;
}

void FHP::setup_dump(){

	this->setup_dump_stats();

	if (this->out_solvation_shell_dump == "__blank__"){
		this->dump_local_ptr = &FHP::dump_local_no_ss;
	}
	else{
		this->dump_local_ptr = &FHP::dump_local_all;
	}
	return;
}

void FHP::setup_run(){
	if (this->v){
		this->run_ptr = &FHP::run_debug;
	}
	else if (this->dry){
		this->run_ptr = &FHP::run_dry;
	}
	else if (this->isotropic){
		this->run_ptr = &FHP::run_isotropic;
	}
	else if (this->polymer){
		this->run_ptr = &FHP::run_simple;
	}
	else{
		this->run_ptr = &FHP::run_straight;
	}
	return;
}

void FHP::setup(){
	// std::cout << "Extracting topology..." << std::endl;
	this->extract_topology_from_file();		// extract the topology
	// std::cout << "Setting up the lattice..." << std::endl;
	this->setup_lattice();					// handles the restart
	// std::cout << "Setting up the energetics..." << std::endl;
	this->setup_energetics();				// set up the energetics
	// std::cout << "Setting up the dumps..." << std::endl;
	this->setup_dump();						// set up the dump functions
	// std::cout << "Setting up the run..." << std::endl;
	this->setup_run();                      // set up the run

	std::cout << "Please be careful about what your energetics are. If you want a polymer simulation, probably  does not make physical sense to have monomer-solvent/cosolvent interactions." << std::endl;
	std::cout << "If you want an isotropic simulation, having aligned and misaligned interactions also does not make sense." << std::endl;

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if ((this->dry || this->polymer) && this->frac_c > 0){
		std::cout << "This is a dry simulation and you have cosolvent. There is a problem. Exiting..." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (this->isotropic){
		for (auto& pair: this->InteractionMap){
			if (std::get<0>(pair.second) != "isotropic"){
				std::cout << "Some interaction is not isotropic, yet you have asked for an isotropic simulation. What's up? Exiting..." << std::endl;
				exit(EXIT_FAILURE);
			} 
		}
	}

	if (this->polymer){
		for (auto& pair: this->InteractionMap){
			if (pair.first.first=="m1" && pair.first.second=="m1"){
				continue;
			}
			else if (std::get<0>(pair.second) != "isotropic" ){
				std::cout << "If you are running a polymer simulation, why do you have non-isotropic monomer-solvent interactions? This will only mess your simulation. Exiting..." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	return;
}
