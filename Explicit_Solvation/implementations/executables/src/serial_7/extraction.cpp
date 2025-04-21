#include "Simulation.hpp"

//////////////////////////////////////////////////////////
//
//         Extract information from inputs
//
//////////////////////////////////////////////////////////

std::vector <std::string> Simulation::extract_content_from_file(std::string filename){

	std::ifstream myfile (filename); 

	if ( !myfile ){
		std::cerr << "File named " << filename << " could not be opened!" << std::endl; 
		exit(EXIT_FAILURE);
	}

	std::string mystring; 
	std::vector <std::string> contents; 

	if (myfile.is_open() ){
		while ( std::getline(myfile, mystring) ) {
			// pipe file's content into stream 
			contents.push_back(mystring); 
		}
	}

	return contents; 
}

void Simulation::extract_number_of_polymers(){

	std::regex start("START"); 
	std::regex end ("END"); 
	std::ifstream position_file (this->positions); 

	if ( !position_file ){
		std::cerr << "File named " << this->positions << " could not be opened!" << std::endl; 
		exit(EXIT_FAILURE);
	}

	std::string my_string; 

	int number_of_starts {0}, number_of_ends {0}; 

	if (position_file.is_open()) {
		while ( std::getline(position_file, my_string) ) {
			
			// std::cout << myString << std::endl; 
			// std::getline(myfile, myString); // pipe file's content into stream 

			if (std::regex_search(my_string, start)){
				++number_of_starts; 
			}
			else if (std::regex_search(my_string, end)){
				++number_of_ends; 
			}
			else if ( my_string.empty()){
				std::cerr << "ERROR: Empty line found. Bad positions file. " << std::endl;
				exit (EXIT_FAILURE);
			}
		}
	}

	if (number_of_starts==number_of_ends){
		this->Npoly = number_of_starts;
	}
	else {
		std::cerr << "Number of starts is not the same as number of ends. Bad input file." << std::endl;
		exit (EXIT_FAILURE);
	}

	return;
}

void Simulation::extract_polymers_from_file(){

	std::vector <std::array <int,3>> locations; 

	std::vector <std::string> contents = this->extract_content_from_file(this->positions); // this extracts every line of the file
	// std::cout << "Extracted contents from file " << positions << "." << std::endl; 

	std::regex start ("START"), end ("END"); 

	int startCount{0}, endCount {0}; 
	
	for (std::string s: contents){
			
		std::stringstream ss(s); 
		if (std::regex_search(s, start) ){
			++startCount;
			continue; 
		}

		else if (std::regex_search(s, end) ) {
			++endCount;
			Polymer pmer = Polymer(locations);
			this->Polymers.push_back(pmer);
			locations.clear();
			
		}

		else{
			std::array <int,3> loc;
			int j{0};  
			for (int i=0; ss>>i; ){
				loc[j] = i;
				++j;
			}

			if (!this->check_validity_of_coords(loc)){
				std::cerr << "Coordinates are out of bounds. Bad input file." << std::endl;
				exit(EXIT_FAILURE); 
			}
			locations.push_back(loc); 
		}
	}

	// std::cout << "Right behind check_for_overlaps." << std::endl;

	if(!(this->check_for_overlaps_within_polymers_raw())){
		std::cerr << "ERROR: There is a problem with the input file for positions. Overlap detected." << std::endl;
		exit(EXIT_FAILURE);  
	}

	// std::cout << "Passed through all contents." << std::endl;

	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	// throw in a check for connectivity of polymer chains 
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

	if (!(this->check_connectivity_raw())){
		std::cerr << "ERROR: There is a problem with the input file for positions. Monomer units are not adjacent to one another on the lattice." << std::endl;
		exit(EXIT_FAILURE); 
	}

	return;
}

void Simulation::extract_topology_from_file(){

	// now, we will get geometry and energy surface. 
	double info; 
	std::array <int,17>     input_hit_check = {0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0};
	std::pair  <std::string, std::string> p; 
	std::tuple <std::string, double, double, int, int> t; 
	std::string interaction_type; 
	std::string mystring; 
	std::vector <std::string> contents = this->extract_content_from_file(this->topology); 

	// defining all the regexes... 
	std::regex x ("^x"),            y       ("^y"),       \
	           z ("^z"),            kT      ("^kT"),      \
	           Emm_a  ("^Emm_a"),   Emm_n   ("^Emm_n"),   \
	           Ems1_a ("^Ems1_a"),  Ems1_n  ("^Ems1_n"),  \
	           Ems2_a ("^Ems2_a"),  Ems2_n  ("^Ems2_n"),  \
	           Es1s2_a("^Es1s2_a"), Es1s2_n ("^Es1s2_n"), \
	           m1_m1  ("^m1-m1"),   m1_s1   ("^m1-s1"),   \
	           m1_s2 ("^m1-s2"),    s1_s2  ("^s1-s2"),    \
	           frac    ("^frac"),    eof ("END OF FILE");

	for (std::string s: contents){

		if (std::regex_search(s, x)){
			info = NumberExtractor(s); 
			input_hit_check[0] += 1;
			if (input_hit_check [0] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->x = info;
			continue; 
		}

		else if (std::regex_search(s, y)){
			info = NumberExtractor(s); 
			input_hit_check[1] += 1;
			if (input_hit_check [1] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->y = info;
			continue; 
		}

		else if (std::regex_search(s, z)){
			info = NumberExtractor(s); 
			input_hit_check[2] += 1;
			if (input_hit_check [2] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->z = info;
			continue; 
		}

		else if (std::regex_search(s, kT)){
			info = NumberExtractor(s); 
			input_hit_check[3] += 1;
			if (input_hit_check [3] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->T = info;
			continue; 
		}

		else if (std::regex_search(s, m1_m1)){
			interaction_type = trim (split(s, ':')[1]);
			if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
				p = std::make_pair ("m1", "m1");
				t = std::make_tuple(interaction_type, 0, 0, 0, 1);
				(this->InteractionMap)[p] = t;
			}
			else {
				std::cout << "\"" << s << "\"" << std::endl;
				std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
				exit(EXIT_FAILURE);     			
			}
			input_hit_check[4] += 1;
			if (input_hit_check [4] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			continue;
		}

		else if (std::regex_search(s, m1_s1)){

			interaction_type = trim (split(s, ':')[1]);
			if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
				p = std::make_pair ("m1", "s1");
				(this->InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 2, 3);
				p = std::make_pair ("s1", "m1") ;
				(this->InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 2, 3);
			}
			else {
				std::cout << "\"" << s << "\"" << std::endl;
				std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
				exit(EXIT_FAILURE); 
			}
			input_hit_check[5] += 1;
			if (input_hit_check [5] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			continue;
		}

		else if (std::regex_search(s, m1_s2)){
			interaction_type = trim (split(s, ':')[1]);
			if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
				p = std::make_pair ("m1", "s2");
				(this->InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 4, 5);
				p = std::make_pair ("s2", "m1");
				(this->InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 4, 5);	
			}
			else {
				std::cout << "\"" << s << "\"" << std::endl;
				std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
				exit(EXIT_FAILURE); 	    		
			}
			input_hit_check[6] += 1;
			if (input_hit_check [6] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			continue;
		}

		else if (std::regex_search(s, s1_s2)){
			interaction_type = trim (split(s, ':')[1]);
			if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
				p = std::make_pair ("s1", "s2");
				(this->InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 6, 7);
				p = std::make_pair ("s2", "s1");
				(this->InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 6, 7);	
			}
			else {
				std::cout << "\"" << s << "\"" << std::endl;
				std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
				exit(EXIT_FAILURE); 	    
			}
			input_hit_check[7] += 1;
			if (input_hit_check [7] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			continue;
		}

		else if (std::regex_search (s, Emm_a)){
			info = NumberExtractor(s); 
			input_hit_check[16] += 1;
			if (input_hit_check [16] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[0] = info;
			continue; 
		}

		else if (std::regex_search (s, Emm_n)){
			info = NumberExtractor(s); 
			input_hit_check[8] += 1;
			if (input_hit_check [8] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[1] = info;
			continue; 
		}

		else if (std::regex_search (s, Ems1_a)){
			info = NumberExtractor(s); 
			input_hit_check[9] += 1;
			if (input_hit_check [9] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[2] = info;
			continue; 
		}

		else if (std::regex_search (s, Ems1_n)){
			info = NumberExtractor(s);
			input_hit_check[10] += 1;
			if (input_hit_check [10] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[3] = info;
			continue;
		}
		else if (std::regex_search (s, Ems2_a)){
			info = NumberExtractor(s);
			input_hit_check[11] += 1;
			if (input_hit_check [11] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[4] = info;
			continue;
		}
		else if (std::regex_search (s, Ems2_n)){
			info = NumberExtractor(s);
			input_hit_check[12] += 1;
			if (input_hit_check [12] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[5] = info;
			continue;
		}
		else if (std::regex_search (s, Es1s2_a)){
			info = NumberExtractor(s);
			input_hit_check[13] += 1;
			if (input_hit_check [13] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[6] = info;
			continue;
		}
		else if (std::regex_search (s, Es1s2_n)){
			info = NumberExtractor(s);
			input_hit_check[14] += 1;
			if (input_hit_check [14] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[7] = info;
			continue;
		}
		else if (std::regex_search (s, frac)) {
			info = NumberExtractor(s);
			input_hit_check[15] += 1;
			if (input_hit_check [15] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->frac_c = info; 
			// std::cout << "From extractor, frac_c = " << info << std::endl;
			continue; 
		}

		else if (std::regex_search(s, eof)){
			int sum = std::accumulate(input_hit_check.begin(), input_hit_check.end(), 0);
			if (sum != 17){
				std::cerr << "Some input was not provided." << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		}

		else {
			std::cout << s << std::endl;
			std::cerr << "ERROR: There is a nonstandard input provided in topology for an FHP simulation file." << std::endl;
			exit(EXIT_FAILURE); 
		}

	}

	p = std::make_pair ("m1", "m1");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[0];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[1];

	p = std::make_pair ("m1", "s1");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[2];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[3];

	p = std::make_pair ("s1", "m1");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[2];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[3];

	p = std::make_pair ("m1", "s2");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[4];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[5];

	p = std::make_pair ("s2", "m1");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[4];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[5];

	p = std::make_pair ("s2", "s1");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[6];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[7];

	p = std::make_pair ("s1", "s2");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[6];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[7];

	return;
}

