#include "Potts.hpp"

void Potts::setup_parser(int lfreq, int max_iter,       \
bool restart_bool,              std::string topology,   \
std::string energy_dump_file,   std::string stats_file, \
std::string lattice_file_write, std::string lattice_file_read) {

	if (!restart_bool){
		if (lfreq == -1 || max_iter == -1){
			std::cerr << "ERROR: No value for option l (frequency of dumping) (" << lfreq << ")" << " and/or for option M (maximum number of moves to be performed) (" << max_iter << ") was provided. Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
		else if (topology == "__blank__" || energy_dump_file == "__blank__" || stats_file == "__blank__" || lattice_file_write == "__blank__"){
			std::cerr << "Topology is in " << topology <<
			",\nenergy dump file is " << energy_dump_file << 
			",\nmove statistics file is " << stats_file << 
			",\nlattice dump file is " << lattice_file_write << std::endl;
			std::cerr << "ERROR: No value for " <<
			"for option t (energy and geometry file) and/or\n" <<
			"for option s (name of move stats file) and/or\n for option u (name of energy dump file). Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
		else {
			std::ofstream out_lattice_dump_file(lattice_file_write, std::ios::out);
			std::ofstream out_energy_dump_file (energy_dump_file, std::ios::out);
			std::ofstream out_statistics_dump_file (stats_file, std::ios::out);

			out_lattice_dump_file.close();
			out_energy_dump_file.close();
			out_statistics_dump_file.close();
		}
	}
	else {
		if (lfreq == -1 || max_iter == -1){
			std::cerr << "ERROR: No value for option l (frequency of dumping) (" << lfreq << ")" << " and/or for option M (maximum number of moves to be performed) (" << max_iter << ") was provided. Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
		else if (topology == "__blank__" || energy_dump_file == "__blank__" || stats_file == "__blank__" || lattice_file_write == "__blank__" || lattice_file_read == "__blank__"){
			std::cerr << "Topology is in " << topology <<
			",\nenergy dump file is " << energy_dump_file << 
			",\nmove statistics file is " << stats_file << 
			",\nlattice dump file is " << lattice_file_write << 
			",\nlattice read file is " << lattice_file_read << "." << std::endl;
			std::cerr << "ERROR: No value for " <<
			"for option t (energy and geometry file) and/or\n" <<
			"for option s (name of move stats file) and/or\n for option u (name of energy dump file). Exiting..." << std::endl;
			exit (EXIT_FAILURE);
		}
	}

	return;
}

void Potts::setup_lattice_from_scratch(){

	std::vector <Particle*> Lattice;
	Lattice.reserve((this->x)*(this->y)*(this->z));

	// update the object
	this->Lattice = Lattice;

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
				Particle* p_ptr = new Particle (rng_uniform(0, 25), "sp", loc);
				this->Lattice.insert(this->Lattice.begin() + lattice_index(loc, this->y, this->z), p_ptr);
			}
		}
	}
	return;
}

void Potts::setup_lattice_from_restart(){

	// set up the lattice object    
	std::vector <Particle*> Lattice;
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
	}
	return;
}

void Potts::setup_lattice(){

	// to restart simulation
	if (this->r){
		this->setup_lattice_from_restart();
		filter_csv(this->out_energy_dump, this->step_number);
	}
	else {
		this->setup_lattice_from_scratch();
	}
	return;
}

void Potts::setup_energetics(){
	this->initial_energetics_map();
	this->neighbor_energetic_map();
	this->energy_compute();
	return;
}

void Potts::setup(){

	this->extract_topology();
	this->setup_lattice();
	this->setup_energetics();
	return;
}