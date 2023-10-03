#include "Simulation.h"
#include "lattice_directions.h"
#include "History.h"

//////////////////////////////////////////////////////////
//
//                Extract information from inputs
//
//////////////////////////////////////////////////////////

const char* ws = " \t\n\r\f\v";
// trim from end of string (right)
inline std::string& rtrim(std::string& s, const char* t = ws)
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from beginning of string (left)
inline std::string& ltrim(std::string& s, const char* t = ws)
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from both ends of string (right then left)
inline std::string& trim(std::string& s, const char* t = ws)
{
    return ltrim(rtrim(s, t), t);
}


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

    std::vector <Polymer> PolymerVector; 
    PolymerVector.reserve(this->Npoly);

    std::vector <std::array <int,3>> locations; 

    std::vector <std::string> contents = this->extract_content_from_file(this->positions); // this extracts every line of the file

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
            PolymerVector.push_back(pmer);
            
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

    if(!(this->check_for_overlaps_within_polymers_raw(&PolymerVector))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Overlap detected." << std::endl;
        exit(EXIT_FAILURE);  
    }

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // throw in a check for connectivity of polymer chains 
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    if (!(this->check_connectivity_raw(&PolymerVector))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Monomer units are not adjacent to one another on the lattice." << std::endl;
        exit(EXIT_FAILURE); 
    }

    this->PSETS.Polymers = PolymerVector;

    return;
}

void Simulation::extract_topology_from_file(){

    double info; 
    std::array <double, 13> info_vec; 
    std::array <int,17> input_hit_check = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::pair  <std::string, std::string> p; 
    std::tuple <std::string, double, double, int, int> t; 
    std::string interaction_type; 
    std::string mystring; 
    std::vector <std::string> contents = this->extract_content_from_file(this->topology); 
    std::regex x ("x"), y ("y"), z ("z"), kT ("kT"), Emm_a ("Emm_a"),\
    Emm_n ("Emm_n"), Ems1_a ("Ems1_a"), Ems1_n ("Ems1_n"), Ems2_a ("Ems2_a"), \
    Ems2_n ("Ems2_n"), Es1s2_a("Es1s2_a"), Es1s2_n ("Es1s2_n"), \
    m1_m1 ("m1-m1"), m1_s1 ("m1-s1"), m1_s2 ("m1-s2"), s1_s2 ("s1-s2"),\
    frac ("frac"), eof ("END OF FILE"); 

    for (std::string s: contents){

    	if (std::regex_search(s, x)){
    		info = NumberExtractor(s); 
    		info_vec[0]=info; 
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
    		info_vec[1] = info; 
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
    		info_vec[2] = info; 
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
    		info_vec[3] = info ; 
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
    		info_vec[4] = info; 
    		input_hit_check[16] += 1;
    		if (input_hit_check [16] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[0] = info;
    		continue; 
    	}

    	else if (std::regex_search (s, Emm_n)){
    		info = NumberExtractor(s); 
    		info_vec[5] = info; 
    		input_hit_check[8] += 1;
    		if (input_hit_check [8] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[1] = info;
    		continue; 
    	}

    	else if (std::regex_search (s, Ems1_a)){
    		info = NumberExtractor(s); 
    		info_vec[6] = info; 
    		input_hit_check[9] += 1;
    		if (input_hit_check [9] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[2] = info;
    		continue; 
    	}

    	else if (std::regex_search (s, Ems1_n)){

    		info = NumberExtractor(s);
    		info_vec[7] = info; 
    		input_hit_check[10] += 1;
    		if (input_hit_check [10] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[3] = info;
    		continue;
    	}
    	else if (std::regex_search (s, Ems2_a)){

    		info = NumberExtractor(s);
    		info_vec[8] = info; 
    		input_hit_check[11] += 1;
    		if (input_hit_check [11] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[4] = info;
    		continue;
    	}
    	else if (std::regex_search (s, Ems2_n)){

    		info = NumberExtractor(s);
    		info_vec[9] = info; 
    		input_hit_check[12] += 1;
    		if (input_hit_check [12] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[5] = info;
    		continue;
    	}
    	else if (std::regex_search (s, Es1s2_a)){

    		info = NumberExtractor(s);
    		info_vec[10] = info; 
    		input_hit_check[13] += 1;
    		if (input_hit_check [13] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[6] = info;
    		continue;
    	}
    	else if (std::regex_search (s, Es1s2_n)){

    		info = NumberExtractor(s);
    		info_vec[11] = info; 
    		input_hit_check[14] += 1;
    		if (input_hit_check [14] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->E[7] = info;
    		continue;
    	}
    	else if (std::regex_search (s, frac)) {
    		info = NumberExtractor(s);
    		info_vec[12] = info;
    		input_hit_check[15] += 1;
    		if (input_hit_check [15] > 1){
    			std::cout << "line = " << s << "." << std::endl;
    			std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
    			exit (EXIT_FAILURE);
    		}
    		this->frac = info; 
    		continue; 
    	}

    	else if (std::regex_search(s, eof)){
    		// std::cout << "End of topology file." << std::endl;
    		break;
    	}

    	else {
    		std::cout << s << std::endl;
    		std::cerr << "ERROR: There is a nonstandard input provided in topology file." << std::endl;
    		exit(EXIT_FAILURE); 
    	}

    }

    p = std::make_pair ("m1", "m1");
    std::get<1>((this->InteractionMap)[p]) = info_vec[4];
    std::get<2>((this->InteractionMap)[p]) = info_vec[5];

    p = std::make_pair ("m1", "s1");
	std::get<1>((this->InteractionMap)[p]) = info_vec[6];
    std::get<2>((this->InteractionMap)[p]) = info_vec[7];

    p = std::make_pair ("s1", "m1");
    std::get<1>((this->InteractionMap)[p]) = info_vec[6];
    std::get<2>((this->InteractionMap)[p]) = info_vec[7];

    p = std::make_pair ("m1", "s2");
    std::get<1>((this->InteractionMap)[p]) = info_vec[8];
    std::get<2>((this->InteractionMap)[p]) = info_vec[9];

	p = std::make_pair ("s2", "m1");
    std::get<1>((this->InteractionMap)[p]) = info_vec[8];
    std::get<2>((this->InteractionMap)[p]) = info_vec[9];

	p = std::make_pair ("s2", "s1");
    std::get<1>((this->InteractionMap)[p]) = info_vec[10];
    std::get<2>((this->InteractionMap)[p]) = info_vec[11];

	p = std::make_pair ("s1", "s2");
    std::get<1>((this->InteractionMap)[p]) = info_vec[10];
    std::get<2>((this->InteractionMap)[p]) = info_vec[11];

    return;
}

//////////////////////////////////////////////////////////
//
//                Create the Lattice
//
//////////////////////////////////////////////////////////


void Simulation::add_solvent(){

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
				this->PSETS.Lattice.insert( this->PSETS.Lattice.begin() + lattice_index(loc, y, z), p_ptr) ;
			}
		}
	}
	return; 
}

void Simulation::add_cosolvent(){

	std::cout << "\n--------------------------------------------------------------------\n" << std::endl;
	std::cout << "Begin adding cosolvent... \n"; 
	int Nmonomer = 0;
	for (Polymer& pmer: this->PSETS.Polymers){
		Nmonomer += pmer.chain.size();
	}
	int nsol2         = std::floor ((this->x*this->y*this->z-Nmonomer)*frac); 
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
		for ( Polymer& pmer: this->PSETS.Polymers){
			for ( Particle*& p: pmer.chain ){
				ne_list = obtain_ne_list ( p->coords, x, y, z); 
				for ( std::array <int,3>& loc: ne_list ){
					if ( (this->PSETS.Lattice).at(lattice_index(loc, y, z))->ptype[0] == 's' ){
						solvation_shell_set.insert(lattice_index (loc, y, z)); 
					}
				}
			}
		}
		std::vector <int> solvation_shell_indices(solvation_shell_set.begin(), solvation_shell_set.end());
		
		while ( (count < nsol2) && (count < static_cast<int>(solvation_shell_indices.size()))) {

			this->PSETS.Lattice.at(solvation_shell_indices.at(i))->ptype = "s2"; 
			this->PSETS.Cosolvent.push_back ( this->PSETS.Lattice.at(solvation_shell_indices.at(i)));
			count += 1;
			i     += 1;
		}
	}
	// std::cout << "count = " << count << std::endl;
	i = 0; 
	while ( count < nsol2 ) {

		if ( this->PSETS.Lattice.at( indices[i] )->ptype[0] == 'm' ){
			; 
		}
		else if ( this->PSETS.Lattice.at(indices[i])->ptype == "s2" ) {
			;
		}
		else {
			this->PSETS.Lattice.at( indices[i] )->ptype = "s2"; 
			this->PSETS.Cosolvent.push_back((this->PSETS.Lattice).at( indices[i] )); 
			count += 1; 
		}
		i += 1; 
	}
	return; 
}

void Simulation::align_lattice(){
	for (Particle*& p: this->PSETS.Lattice) {
		p->orientation = 0;
	}

	return; 
}

void Simulation::align_solvation_shell(){

    std::array  <std::array<int,3>, 26> ne_list;
    std::vector <int> solvent_indices; 
    
    for ( Polymer& pmer: this->PSETS.Polymers ){
        for ( Particle*& p: pmer.chain ) {
            
            ne_list = obtain_ne_list (p->coords, x, y, z);
            for ( std::array<int,3>& ne: ne_list ){
                if ( this->PSETS.Lattice.at(lattice_index(ne, y, z))->ptype[0] == 's' ) {
                    solvent_indices.push_back (lattice_index(ne, y, z));
                }
            }
        }
    }  

    // get rid of duplicates 
    std::unordered_set <int> s (solvent_indices.begin(), solvent_indices.end() );
    solvent_indices.assign (s.begin(), s.end() ); 

    for ( Polymer& pmer: this->PSETS.Polymers ) {
        for ( Particle*& p: pmer.chain ) {
            p->orientation = 0;
        }
    }

    for (int i: solvent_indices){
        this->PSETS.Lattice[i]->orientation = 0;
    }

    return;
}

void Simulation::extract_polymers_from_restart(){


	std::vector <Polymer> PolymerVector;
    std::vector <std::array<int,3>> locations; 
    std::vector <int> spins;

    std::vector <std::string> contents = this->extract_content_from_file(this->dfile); // this extracts every line of the file

    // std::vector <int> step_num_store; 
    std::vector <int> index_store; // need this guy to hold the index of the final set of coordinates. 

    // std::cout << "Is content being extracted?" << std::endl;

    int final_step_num = this->step_number;
    bool step_bool {false}, start_bool {false}, end_bool {false}; 

    std::regex start ("START"), end ("END"), step ("step " + std::to_string(final_step_num) );
    std::regex step_generic ("step"); 
    std::regex reg_poly ("Dumping coordinates of Polymer"); 
    std::regex numbers ("[0-9]+"); 

    int startCount{0}, endCount{0}; 
    std::array <int,3> loc;
    // std::stringstream ss; 
    std::smatch match; 

    // std::cout << "final step number is " << final_step_num << std::endl;

    for (std::string& s: contents){
    	
    	if ( std::regex_search (s, step_generic) ){

    		std::regex_search ( s, match, numbers ); 
			std::regex_token_iterator<std::string::iterator> rend; 
			std::regex_token_iterator<std::string::iterator> a ( s.begin(), s.end(), numbers );

			// std::cout << "*a is " << std::stoi(*a) << std::endl;

			if ( std::stoi(*a) > final_step_num ){
				// std::cout << "*a is " << std::stoi(*a) << std::endl;
				// std::cout << "final_step_num = " << final_step_num << std::endl;				
				std::cerr << "\n\nYou coordinates file and restart file are not in sync. \nIt is probably because the restart file you chose is not for that simulation.\nPlease check them out. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}

    	}

    	if ( std::regex_search ( s, step) ) {
    		// std::cout << s << std::endl;
    		step_bool = true; 
    		continue; 
    	}

        // sending to the stringstream
        // ss (s); 

        if ( step_bool ){

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

	            Polymer pmer = Polymer (locations, spins);
	            PolymerVector.push_back(pmer);
	            
	            locations.clear();
	            break;
	            // continue; 
	            
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
				
				// std::cout << "Location is "; print (loc);
	            
	            if (!this->check_validity_of_coords(loc)){
	            	std::cerr << "Coordinates are out of bounds. Bad input file." << std::endl;
	            	exit(EXIT_FAILURE); 
	            }
	        
	            locations.push_back(loc); 
	            
	        }
    	}
    }

    if(!(this->check_for_overlaps_within_polymers_raw(&PolymerVector))){
        std::cerr << "ERROR: There is a problem with the input file for positions. Overlap detected." << std::endl;
        exit(EXIT_FAILURE);  
    }

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
    // throw in a check for connectivity of polymer chains 
    //~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

    if (!(this->check_connectivity())){
        std::cerr << "ERROR: There is a problem with the input file for positions. Monomer units are not adjacent to one another on the lattice." << std::endl;
        exit(EXIT_FAILURE); 
    }
    this->PSETS.Polymers = PolymerVector;
    
    return; 
}

void Simulation::extract_lattice_from_restart(){
	std::vector <std::string> contents = this->extract_content_from_file(this->lattice_file_read);
	std::regex start ("FINAL STEP: "), end ("END"); 
	std::regex numbers ("[0-9]+"); 
	std::regex characters ("[a-z][0-9]+");

	int orientation   = -1; 
	std::string ptype = "x"; 
	int index         = -1; 
	std::smatch match; 

	int  final_step = 0;
	bool start_recording = false; 

	// find the last step dumped in lattice file
	for ( std::string& s: contents) {

		if ( std::regex_search (s, start) ) {
			std::regex_search (s, match, numbers);
			final_step = std::stoi (match[0].str());
		}
	}

	this->step_number = final_step;

	for ( std::string& s: contents) {

		// send content into stringstream 

		if ( std::regex_search (s, start) ) {
			std::regex_search ( s, match, numbers ); 
			if ( final_step == std::stoi(match[0].str()) ) {
				start_recording = true;
			}
		}

		else if ( std::regex_search (s, end) && start_recording ){
			break; 
		}

		else if (start_recording) {
			std::regex_search ( s, match, numbers );
			std::regex_token_iterator<std::string::iterator> rend; 
			std::regex_token_iterator<std::string::iterator> a ( s.begin(), s.end(), numbers );

			for ( int i=0; i<2; ++i ){
				if ( i == 0 ){
					orientation = std::stoi ( *a );
					*a++;  
				}
				else {
					*a++; 
					index = std::stoi ( *a );
				}
			}
 

			std::regex_search ( s, match, characters );
			ptype = match[0].str(); 

			// std::cout << "ptype is " << ptype << std::endl; 

			Particle* p_ptr = new Particle (location(index, x, y, z), ptype, orientation); 

			this->PSETS.Lattice.insert ( this->PSETS.Lattice.begin() + index, p_ptr); 

		}

	}

	if (!start_recording) {
		std::cout <<"Something is wrong with the restart!" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::cout << "Created lattice from file!" << std::endl;
	return;
}

void Simulation::set_up_lattice_from_restart(){

	this->extract_lattice_from_restart();
	this->extract_polymers_from_restart();

	for (Polymer& pmer: this->PSETS.Polymers){
		for (Particle*& p: pmer.chain){

			if ( !( this->PSETS.Lattice[ lattice_index(p->coords, y, z) ]->ptype == p->ptype && this->PSETS.Lattice[ lattice_index(p->coords, y, z) ]->coords == p->coords && (this->PSETS.Lattice) [ lattice_index(p->coords, y, z) ]->orientation == p->orientation) ) {
				std::cerr << "There is a problem with extraction for restart..." << std::endl;
				exit (EXIT_FAILURE); 
			}
			this->PSETS.Lattice[lattice_index(p->coords, y, z)] = p;

		}
	}

	for (int i{0}; i < x*y*z; ++i){
		if ( this->PSETS.Lattice.at(i)->ptype == "s2" ){
			this->PSETS.Cosolvent.push_back((this->PSETS.Lattice).at(i));
		}
	}

	return;
}

void Simulation::set_up_lattice(){

	// initialize custom data structures
	// this data structure will hold the coordinates of the polymer
	std::vector <Polymer> Polymers;
	Polymers.reserve(this->Npoly);

	// this data structure will hold the coordinates of the cosolvent
	std::vector <Particle*> Cosolvent;

	// this data structure will hold the lattice
	std::vector <Particle*> Lattice; 
	Lattice.reserve(this->x*this->y*this->z);

	this->PSETS.Lattice   = Lattice;
	this->PSETS.Polymers  = Polymers;
	this->PSETS.Cosolvent = Cosolvent;
	

	if (!this->r){

		this->extract_polymers_from_file ();
		this->add_solvent();

		// populate the lattice
		for (Polymer& pmer: (this->PSETS.Polymers)) {
			for (Particle*& p: pmer.chain){
				// now that I have my polymer coordinates, time to create the grand lattice 
				(this->PSETS.Lattice).at(lattice_index (p->coords, y, z) ) = p; 
			}
		}

		this->add_cosolvent (); 

		if (this->A){
			this->align_lattice();
		}
		else if (this->S){
			this->align_solvation_shell();
		}

	}

	else {
		this->set_up_lattice_from_restart();
	}

	this->check_structures();
	return;
}

//////////////////////////////////////////////////////////
//   
//                Check Structures
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

bool Simulation::check_for_overlaps_within_polymers_raw(std::vector <Polymer>* Polymers){

	std::vector <std::array <int,3>> loc_list; 

	for (Polymer& pmer: *Polymers){
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

bool Simulation::check_for_overlaps_within_polymers(){

	std::vector <std::array <int,3>> loc_list; 

	for (Polymer& pmer: this->PSETS.Polymers){
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

    for (Polymer& pmer: this->PSETS.Polymers) {
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
    
    for (Polymer& pmer: this->PSETS.Polymers) {
        for (Particle*& p: pmer.chain){

        	if ( this->PSETS.Lattice[ lattice_index (p->coords, y, z)]->coords == p->coords ) {

        		if ( this->PSETS.Lattice[ lattice_index(p->coords, y, z)]->ptype[0] == 's' ){
        			std::cerr << "Some kind of bad solvent-monomer overlap that has taken a place. A monomer is being represented by a solvent. Something's fucked." << std::endl;
        			std::cerr << "Location is: "; print (p->coords); 
        			std::cerr << "Type is: " << (this->PSETS.Lattice[ lattice_index (p->coords, y, z)]->ptype) << std::endl; 
        			std::cerr << "Type is: " << p->ptype << std::endl; 
        			std::cerr << "Location on lattice is: "; print (this->PSETS.Lattice[ lattice_index (p->coords, y, z)]->coords); 
        			exit (EXIT_FAILURE); 
        		}
    			else {
    				continue;
    			}

        	}

        	else {
        		std::cerr << "Something is wrong with the Lattice map. Monomer unit lost. Something is fucked. " << std::endl;
        		std::cerr << "Output from *Polymers: "; print (p->coords);
        		std::cerr << "Output from this->PSETS.Lattice: "; print (this->PSETS.Lattice[ lattice_index (p->coords, y, z)]->coords);
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
	for ( int i{0}; i< static_cast<int>(this->PSETS.Lattice.size()); ++i ){
		
		pmer_loop_flag = false; 
		particle_found_flag = false; 

		if ( this->PSETS.Lattice[i]->ptype[0] == 'm' ){

			for ( Polymer& pmer: this->PSETS.Polymers ){
				for (Particle*& p: pmer.chain){
					if (this->PSETS.Lattice[i]->coords == p->coords){
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
				std::cout << "Location = "; print (this->PSETS.Lattice[i]->coords);
				exit (EXIT_FAILURE); 
			}
		}

	}

	return true;
}

bool Simulation::check_connectivity_raw(std::vector<Polymer>* Polymers){

    for (Polymer& pmer: *Polymers){
        size_t length = pmer.chain.size(); 
        std::array <int,3> connection = {0,0,0}; 
        std::sort (adrns.begin(), adrns.end() ); 

        for (int i{1}; i<static_cast<int>(length); ++i){
            
            connection = subtract_arrays(&(pmer.chain[i]->coords), &(pmer.chain[i-1]->coords));
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

    for (Polymer& pmer: this->PSETS.Polymers){
        size_t length = pmer.chain.size(); 
        std::array <int,3> connection = {0,0,0}; 
        std::sort (adrns.begin(), adrns.end() ); 

        for (int i{1}; i<static_cast<int>(length); ++i){
            
            connection = subtract_arrays(&(pmer.chain[i]->coords), &(pmer.chain[i-1]->coords));
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

		if ( location(i, this->x, this->y, this->z) == (this->PSETS.Lattice)[i]->coords ) {
			continue;
		}
		else {
			std::cout << "Problem with pointers." << std::endl;
			std::cout << "Bad location is "; print ( location (i, x, y, z) ); 
			std::cout << "Lattice says "; print ((this->PSETS.Lattice)[i]->coords);
			std::cerr << "Something is fucked. Pointer does not correspond to position. " << std::endl;
			exit (EXIT_FAILURE);
			return false;
		}

	}

	return true;
}

void Simulation::check_structures(){

	std::cout << "Checking validity of coords..." << std::endl;
    std::cout << std::boolalpha << "checkForOverlaps says: " << this->check_for_overlaps_within_polymers() << "." << std::endl; 
    if (!this->check_for_overlaps_within_polymers()){
        std::cout << "Something is fucked up overlaps-wise in the polymer itself." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (!this->check_for_overlaps_on_lattice()){
    	std::cout << "Random monomer floating in the lattice! Breaking out..." << std::endl; 
        exit (EXIT_FAILURE);
    }
    std::cout << "No random monomers floating around!!" << std::endl;

    if (!this->check_for_solvent_monomer_overlap()){
    	std::cout << "Something is fucked up solvent-monomer overlaps-wise. " << std::endl; 
        exit(EXIT_FAILURE);
    }

    std::cout << std::boolalpha << "checkConnectivity says: " << this->check_connectivity() << "." << std::endl;
    if (!this->check_connectivity()){
    	std::cout << "Something is fucked up connectivity-wise." << std::endl; 
        exit(EXIT_FAILURE);
    }

    if (this->check_pointers_on_lattice()){
        std::cout << "We good. Lattice is in good shape." << std::endl; 
    }
    else {
        std::cerr <<"Something is fucked with pointers on Lattice." << std::endl; 
        exit (EXIT_FAILURE);
    }


    for (Particle*& p: this->PSETS.Cosolvent){

    	if ( this->PSETS.Lattice.at( lattice_index(p->coords, this->y, this->z) )->ptype != "s2" && this->PSETS.Lattice.at( lattice_index(p->coords, this->y, this->z) )->ptype != p->ptype && (this->PSETS.Lattice).at(lattice_index(p->coords, this->y, this->z))->orientation != p->orientation && \
    	this->PSETS.Lattice.at( lattice_index(p->coords, this->y, this->z) )->coords != p->coords ) {
    		std::cerr << "Cosolvent and Lattice do not align. Something's fucked." << std::endl;
    		exit (EXIT_FAILURE);
    	}

    }
    std::cout << "Cosolvent vector looks good!" << std::endl;

    return; 
}

//////////////////////////////////////////////////////////
//   
//                Run the energetic calculation
//
//////////////////////////////////////////////////////////

void Simulation::pair_interaction(Particle* p1, Particle* p2, double* energy){

	double dot_product = 0;
	double theta_1     = 0;
	double theta_2     = 0;
	double magnitude   = 0;

	std::array <int,3> connvec = {0, 0, 0};

	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 

	std::string interaction_type = std::get<0>((this->InteractionMap)[particle_pair]);
	switch (interaction_type[0]){

	case 'i':
		if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){

			(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 0.5;
			*energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 
		}
		else {
			// std::cout << "particle_pair = " << particle_pair.first << ", " << particle_pair.second << std::endl;
			(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1; 	
			*energy += std::get<1>((this->InteractionMap)[particle_pair]); 
		}
		break;

	case 'p':
		dot_product = take_dot_product ( p1->orientation, p2->orientation );
		if (dot_product > 0.54){
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 0.5; 
				
				*energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 

			}
			else {
				(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1;
				*energy += std::get<1>((this->InteractionMap)[particle_pair]); 
			}	
		}
		else {
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				// std::cout << "misaligned parallel m1-m1..." << std::endl;
				(this->contacts)[std::get<4>((this->InteractionMap)[particle_pair])] += 0.5; 
				*energy += std::get<2>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(this->contacts)[std::get<4>((this->InteractionMap)[particle_pair])] += 1; 
				*energy += std::get<2>((this->InteractionMap)[particle_pair]); 
			}
		}
		break;

	case 'a':
		connvec = subtract_arrays ( &(p2->coords), &(p1->coords) );
		modified_direction (&connvec, x, y, z); 
		magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
		theta_1   = std::acos (take_dot_product (scale_arrays (1/magnitude , &connvec), Or2Dir[p1->orientation] ) );
		theta_2   = std::acos (take_dot_product (scale_arrays (-1/magnitude, &connvec), Or2Dir[p2->orientation] ) );

		if ( theta_1 + theta_2 > M_PI/2 ){
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ) {
				// std::cout << "misaligned antiparallel m1-m1..." << std::endl;
				(this->contacts)[std::get<4>((this->InteractionMap)[particle_pair])] += 0.5; 
				*energy += std::get<2>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(this->contacts)[std::get<4>((this->InteractionMap)[particle_pair])] += 1; 	
				*energy += std::get<2>((this->InteractionMap)[particle_pair]);
			}
		}
		else {
			if ( particle_pair.first == "m1" && particle_pair.second == "m1" ){
				// std::cout << "aligned antiparallel m1-m1..." << std::endl;
				(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 0.5; 
				*energy += std::get<1>((this->InteractionMap)[particle_pair])*0.5; 
			}
			else {
				(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1; 
				*energy += std::get<1>((this->InteractionMap)[particle_pair]); 
			}
		}
		break;
	}

	return;
}

void Simulation::modified_pair_interaction(Particle* p1, Particle* p2, double* energy){

	double dot_product = 0;
	double theta_1     = 0;
	double theta_2     = 0;
	double magnitude   = 0;

	std::array <int,3> connvec = {0, 0, 0};

	std::pair <std::string, std::string> particle_pair = std::make_pair (p1->ptype, p2->ptype); 

	switch ( std::get<0>((this->InteractionMap)[particle_pair])[0] ) {

		case 'i':
			(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1; 
			*energy += std::get<1>((this->InteractionMap)[particle_pair]);
			break;

		case 'p':
				dot_product = take_dot_product ( p1->orientation, p2->orientation );
				if (dot_product > 0.54){
					(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1; 
					*energy += std::get<1>((this->InteractionMap)[particle_pair]); 
				}
				else {
					(this->contacts)[std::get<4>((this->InteractionMap)[particle_pair])] += 1; 
					*energy += std::get<2>((this->InteractionMap)[particle_pair]); 
				}	
				break;

		case 'a':
			connvec = subtract_arrays ( &(p2->coords), &(p1->coords) );
			modified_direction ( &connvec, x, y, z); 
			magnitude = std::sqrt (connvec[0]*connvec[0] + connvec[1]*connvec[1] + connvec[2]*connvec[2]);
			theta_1   = std::acos (take_dot_product (scale_arrays (1/magnitude , &connvec), Or2Dir[p1->orientation] ) );
			theta_2   = std::acos (take_dot_product (scale_arrays (-1/magnitude, &connvec), Or2Dir[p2->orientation] ) );

			if ( theta_1 + theta_2 > M_PI/2 ){

				(this->contacts)[std::get<4>((this->InteractionMap)[particle_pair])] += 1; 
				*energy += std::get<2>((this->InteractionMap)[particle_pair]); 
			}
			else {
				(this->contacts)[std::get<3>((this->InteractionMap)[particle_pair])] += 1; 
				*energy += std::get<1>((this->InteractionMap)[particle_pair]); 					
			}
			break;

	}

	return; 
}

void Simulation::calculate_energetics(){

	double Energy {0.0};
	this->contacts = {0,0,0,0,0,0,0,0}; 

	std::array <std::array <int,3>, 26> ne_list; 

	// run energy computations for every monomer bead 
	// auto start = std::chrono::high_resolution_clock::now();

	for (Polymer& pmer: this->PSETS.Polymers) {
		for (Particle*& p: pmer.chain){
			ne_list = obtain_ne_list(p->coords, x, y, z); // get neighbor list 
			for ( std::array <int,3>& loc: ne_list) {
				this->pair_interaction(p, this->PSETS.Lattice[lattice_index(loc, y, z)], &Energy);
				// PairInteractionForRevampedEnergyCalc (p, (this->PSETS.Lattice)[ lattice_index(loc, y, z) ], InteractionMap, &Energy, contacts, x, y, z);
			}
		}
	}

	for (Particle*& p: this->PSETS.Cosolvent){
		ne_list = obtain_ne_list ( p->coords, x, y, z );
		for ( std::array <int,3>& loc: ne_list ){
			if ( this->PSETS.Lattice[lattice_index(loc, y, z) ]->ptype == "m1" || this->PSETS.Lattice[ lattice_index(loc, y, z) ]->ptype == "s2"){
				continue; 
			}
			this->modified_pair_interaction(p, this->PSETS.Lattice[lattice_index(loc, y, z)], &Energy);
			// ParticlePairEnergyContribution (p, (this->PSETS.Lattice)[ lattice_index(loc, y, z) ], InteractionMap, &Energy, contacts, x, y, z);
		}
	}

	this->sysEnergy = Energy;
	return; 
}

//////////////////////////////////////////////////////////
//   
//                Dump functions
//
//////////////////////////////////////////////////////////

void Simulation::dump_energy(int step_num){

	std::ofstream dump_file(this->efile, std::ios::app); 

	dump_file << sysEnergy << " | " \
			<< (this->contacts)[0]+(this->contacts)[1] << " | " << (this->contacts)[0] << " | " << (this->contacts)[1] << " | " \
			<< (this->contacts)[2]+(this->contacts)[3] << " | " << (this->contacts)[2] << " | " << (this->contacts)[3] << " | " \
			<< (this->contacts)[4]+(this->contacts)[5] << " | " << (this->contacts)[4] << " | " << (this->contacts)[5] << " | " \
			<< (this->contacts)[6]+(this->contacts)[7] << " | " << (this->contacts)[6] << " | " << (this->contacts)[7] << " | " << step_num << "\n";

	return; 
}

void Simulation::dump_polymers(int step_num){
    std::ofstream dump_file(this->dfile, std::ios::app); 
    dump_file <<"Dumping coordinates at step " << step_num << ".\n";
    
    int count = 0; 
    for (Polymer& pmer: this->PSETS.Polymers){
        
        dump_file <<"Dumping coordinates of Polymer # " << count << ".\n";
        dump_file<<"START" << "\n";
        for (Particle*& p: pmer.chain){
            for (int i: p->coords){
                dump_file << i << " | "; 
            }
            dump_file << p->orientation << " | ";
            dump_file << "\n"; 
        }
        ++count ; 
        dump_file <<"END"<<"\n";
    }
    
    dump_file << "~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#\n";
    return; 
}

void Simulation::dump_solvation_shell_orientations(int step_num){
    
    std::ofstream dump_file (this->mfile, std::ios::app); 
    
    dump_file << "START for Step " << step_num << ".\n";
    std::vector <int> solvent_indices; 
    
    for ( Polymer& pmer: this->PSETS.Polymers ) {
        for ( Particle*& p: pmer.chain ) {
            
            dump_file << p->orientation << " | ";
            std::array <std::array<int,3>, 26> ne_list = obtain_ne_list (p->coords, x, y, z) ;
            
            for ( std::array <int,3>& ne: ne_list) {
                
                // std::cout << "Reported~\n"; 
                
                if (this->PSETS.Lattice[lattice_index(ne, y, z)]->ptype[0] == 's' && std::find(solvent_indices.begin(), solvent_indices.end(), lattice_index(ne, y, z)) == solvent_indices.end()){
                    dump_file << (this->PSETS.Lattice[lattice_index(ne, y, z)])->orientation << " | ";  
                } 
            }
            dump_file << "\n"; 
        } 
    }
    
    dump_file << "END. \n";
    return;
}

void Simulation::dump_statistics(int step_num){

	std::ofstream dump_file (this->stats_file, std::ios::out); 
	dump_file << "At step " << step_num << "...\n";

	dump_file << "End rotations without bias         - attempts: " << (this->attempts)[0] <<", acceptances: " << (this->acceptances)[0] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[0])/static_cast<double>((this->attempts)[0]) << std::endl; 
	dump_file << "Reptation without bias             - attempts: " << (this->attempts)[1] <<", acceptances: " << (this->acceptances)[1] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[1])/static_cast<double>((this->attempts)[1]) << std::endl; 
	dump_file << "Chain regrowth with overlap bias   - attempts: " << (this->attempts)[2] <<", acceptances: " << (this->acceptances)[2] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[2])/static_cast<double>((this->attempts)[2]) << std::endl; 
	dump_file << "Chain regrowth with ori flip       - attempts: " << (this->attempts)[3] <<", acceptances: " << (this->acceptances)[3] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[3])/static_cast<double>((this->attempts)[3]) << std::endl; 
	dump_file << "Solvent flips without bias         - attempts: " << (this->attempts)[4] <<", acceptances: " << (this->acceptances)[4] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[4])/static_cast<double>((this->attempts)[4]) << std::endl;
	dump_file << "Solvation shell flip with bias     - attempts: " << (this->attempts)[5] <<", acceptances: " << (this->acceptances)[5] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[5])/static_cast<double>((this->attempts)[5]) << std::endl;
	dump_file << "Polymer flips                      - attempts: " << (this->attempts)[6] <<", acceptances: " << (this->acceptances)[6] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[6])/static_cast<double>((this->attempts)[6]) << std::endl;
	dump_file << "Solvent exchange with bias         - attempts: " << (this->attempts)[7] <<", acceptances: " << (this->acceptances)[7] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[7])/static_cast<double>((this->attempts)[7]) << std::endl;
	dump_file << "Solvent exchange without bias      - attempts: " << (this->attempts)[8] <<", acceptances: " << (this->acceptances)[8] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[8])/static_cast<double>((this->attempts)[8]) << std::endl;

	return;
}

void Simulation::dump_lattice(int step_num){
	std::ofstream dump_file ( this->lattice_file_write, std::ios::out ); 
	dump_file << "FINAL STEP: " << step_num << ".\n"; 
	for ( Particle*& p: this->PSETS.Lattice ){
		dump_file << p->orientation << ", " << p->ptype << ", " << lattice_index(p->coords, y, z) << "\n"; 
	}

	dump_file << "END. \n";
	return; 
}

//////////////////////////////////////////////////////////
//   
//                Perturb the system
//
//////////////////////////////////////////////////////////


void Simulation::swing_monomer(int m, std::vector<History>* history_store, double* energy_forw, std::array <int,8>* contacts, double* prob_forw){

    if (m == this->PSETS.Polymers[0].size() - 1){
        this->IMP_BOOL = true;
        return;
    }

    History move;

    std::array <double,5> energies;
    std::array <double,8>               current_contacts = *contacts; 
    std::array <std::array<double,8>,5> contacts_store; 
    std::array <double,5> boltzmann;
    std::array <int,3>    loc_m = this->PSETS.Polymers[0].chain[m+1]->coords;

    // generate possible locations to jump to
    std::array <std::array <int,3>, 26> ne_list = obtain_ne_list (this->PSETS.Polymers[0].chain[m]->coords, this->x, this->y, this->z);

    // randomly select five of them
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle (ne_list.begin(), ne_list.end(), std::default_random_engine(seed)); 

    // doubles for energy transfer 
    double rboltzmann = 0; // running sum for boltzmann weights 
    double Em         = 0;
    double Es         = 0; 
    double Em_n       = 0; 
    double Es_n       = 0; 
    double Esys       = *energy_forw; 
    double Epair      = 0;
    double Epair_n    = 0;

    std::array <double,8> cm   = {0,0,0,0,0,0,0,0};
    std::array <double,8> cs   = {0,0,0,0,0,0,0,0};
    std::array <double,8> cm_n = {0,0,0,0,0,0,0,0};
    std::array <double,8> cs_n = {0,0,0,0,0,0,0,0};    

    // start attempting jumps 

    int block_counter       = 0 ; 
    int idx_counter         = 0 ; 
    int self_swap_idx       = -1; 
    int c_idx               = 0;
    int c_idx_n             = 0;

    while (idx_counter<5) {
        if (ne_list[idx_counter] == loc_m){
            // current position
            energies[idx_counter]       = Esys;
            contacts_store[idx_counter] = current_contacts;
            block_counter += 1;
        }

        else if ( this->PSETS.Lattice[lattice_index(ne_list[idx_counter], this->y, this->z)]->ptype[0] == 'm' ){

            for (int u{0}; u<deg_poly; ++u){
                if ( this->PSETS.Polymers[0].chain[u]->coords == ne_list[idx_counter] ){
                    self_swap_idx = u;
                    break;
                }
            }

            // if it is with other head units, do the swap 
            // if not, discourage it 
            // std::cout << "self_swap_idx = " << self_swap_idx << std::endl; 

            if (self_swap_idx < m){
                // maintain current state, and sample another state 
                energies[idx_counter] = 1e+08; // very unfavorable state 
                contacts_store[idx_counter] = {-1,-1,-1,-1,-1,-1,-1,-1}; 
                block_counter += 1;   
            }
            else {
                // get the initial neighboring energies. 
                Es = this->neighbor_energetics (&cs, lattice_index(ne_list[idx_counter], this->y, this->z));
                Em = this->neighbor_energetics (&cm, lattice_index(loc_m, this->y, this->z)); 
                Epair = isolated_pair_particle_interaction (this->PSETS.Lattice[lattice_index (ne_list[idx_counter], this->y, this->z)], this->PSETS.Lattice[lattice_index (loc_m, this->y, this->z)], &c_idx);

                // prep the swap 
                (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
                (*Polymers)[p_index].chain[m_index+1]->coords       = ne_list[idx_counter];
                
                // perform the swap (since coords were changed, this swap works)
                (*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
                (*LATTICE)[ lattice_index (ne_list[idx_counter], y, z) ] = (*Polymers)[p_index].chain[m_index+1];

                // get the new energies
                Es_n = neighbor_energetics(&cs_n, lattice_index(ne_list[idx_counter], this->y, this->z));
                Em_n = neighbor_energetics(&cm_n, lattice_index(loc_m, this->y, this->z)); 
                Epair_n = isolated_pair_particle_interaction(this->PSETS.Lattice[lattice_index (ne_list[idx_counter], this->y, this->z)], this->PSETS.Lattice[lattice_index (loc_m, this->y, this->z)], &c_idx_n);

                energies       [idx_counter] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
                contacts_store [idx_counter] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) ); 
                contacts_store [idx_counter][c_idx] += 1;
                contacts_store [idx_counter][c_idx_n] -= 1; 

                // revert back to original structure 
                this->PSETS.Lattice[lattice_index(loc_m, this->y, this->z)]->coords = ne_list[idx_counter];
                this->PSETS.Polymers[0].chain[m+1]->coords  = loc_m;
                
                // perform the swap (since coords were changed, this swap works)
                this->PSETS.Lattice[lattice_index(ne_list[idx_counter], y, z)]    = this->PSETS.Lattice[lattice_index (loc_m, y, z)];
                this->PSETS.Lattice[lattice_index(loc_m, y, z)]                   = this->PSETS.Polymers[0].chain[m+1];

            }

        }

        else {
            // std::cout << "solvent swap time." << std::endl;
            Es = neighbor_energetics (&cs, lattice_index(ne_list[idx_counter], this->y, this->z));
            Em = neighbor_energetics (&cm, lattice_index(loc_m, this->y, this->z)); 
            Epair = IsolatedPairParticleInteraction (this->PSETS.Lattice[lattice_index (ne_list[idx_counter], this->y, this->z)], (*LATTICE)[lattice_index (loc_m, this->y, this->z)], &c_idx);

            // prep the swap 
            this->PSETS.Lattice[lattice_index(ne_list[idx_counter], y, z)]->coords = loc_m;
            this->PSETS.Polymers[0].chain[m_index+1]->coords = ne_list[idx_counter];
            
            // perform the swap (since coords were changed, this swap works)
            this->PSETS.Lattice[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[idx_counter], y, z)];
            this->PSETS.Lattice[ lattice_index (ne_list[idx_counter], y, z) ] = this->PSETS.Polymers[p_index].chain[m_index+1];

            // get the new energies
            Es_n = neighbor_energetics (&cs_n, lattice_index(ne_list[idx_counter], this->y, this->z));
            Em_n = neighbor_energetics (&cm_n, lattice_index(loc_m, this->y, this->z));
            Epair_n = isolated_pair_particle_interaction (this->PSETS.Lattice[lattice_index (ne_list[idx_counter], y, z)], this->PSETS.Lattice[lattice_index (loc_m, this->y, this->z)], &c_idx_n);

            energies       [ idx_counter ] = Esys - (Es+Em-Epair) + (Es_n+Em_n-Epair_n); 
            contacts_store [ idx_counter ] = add_arrays (subtract_arrays (current_contacts, add_arrays (cs, cm)), add_arrays (cs_n, cm_n) );    
            contacts_store [ idx_counter ][ c_idx ]   += 1;
            contacts_store [ idx_counter ][ c_idx_n ] -= 1;

            // revert back to original structure 
            this->PSETS.Lattice[lattice_index(loc_m, y, z)]->coords = ne_list[idx_counter];
            this->PSETS.Polymers[0].chain[m_index+1]->coords  = loc_m;
            
            // perform the swap (since coords were changed, this swap works)
            this->PSETS.Lattice[lattice_index(ne_list[idx_counter], y, z)]    = this->PSETS.Lattice[lattice_index (loc_m, y, z)];
            this->PSETS.Lattice[lattice_index(loc_m, y, z)]         = this->PSETS.Polymers[0].chain[m_index+1];

        }
        idx_counter += 1;
    }

    if ( block_counter == 5 ){
        this->IMP_BOOL = false; 
        return; 
    }

    // now that i have all the energies, and boltzmann weights, i can choose a configuration 
    // std::cout << "Energies are "; print(energies); 

    double Emin = *std::min_element ( energies.begin(), energies.end() ); 
    // std::cout << "Emin = " << Emin << std::endl; 

    for (int i{0}; i<5; ++i){
        boltzmann[i] = std::exp(-1/temperature*( energies[i] - Emin ) ); 
        rboltzmann  += boltzmann[i]; 
    }

    // std::cout << "Boltzmann weights are: "; print(boltzmann); 
    double rng_acc = rng_uniform (0.0, 1.0); 
    double rsum    = 0; 
    int    e_idx   = 0; 
    
    // std::cout << "normalization = " << rboltzmann << std::endl;
    // std::cout << "rng = " << rng_acc << std::endl;

    for (int j{0}; j<5; ++j){
        rsum += boltzmann[j]/rboltzmann; 
        if ( rng_acc < rsum ){
            e_idx = j; 
            break; 
        }   
    }

    // now that I have chosen a configuration, go with it
    *prob_forw     *= boltzmann[e_idx]/rboltzmann; 
    *energy_forw    = energies[e_idx];
    *contacts_forw  = contacts_store [e_idx]; 

    // update the history 
    
    // do the swap again
    (*LATTICE)[lattice_index(ne_list[e_idx], y, z)]->coords = loc_m;
    (*Polymers)[0].chain[m_index+1]->coords           = ne_list[e_idx];

    // perform the swap (since coords were changed, this swap works)
    (*LATTICE)[ lattice_index (loc_m, y, z) ] = (*LATTICE)[lattice_index(ne_list[e_idx], y, z)];
    (*LATTICE)[ lattice_index (ne_list[e_idx], y, z) ]  = (*Polymers)[0].chain[m_index+1];
    move.LocH[(*Polymers)[0].chain[m_index+1]].push_back(lattice_index((*Polymers)[0].chain[m_index+1]->coords, this->y, this->z));

    return;

}






void Simulation::iced_regrowth(){

	// assumption that only one polymer is in Polymers
	int deg_poly = this->PSETS.Polymers[0].chain.size();
	int m_index  = deg_poly/2 + 1; // rng_uniform(1, deg_poly-2);

	double prob_forw {1}; 
	double prob_back {1}; 

    double energy_curr {this->sysEnergy};
    double energy_back {this->sysEnergy};
	double energy_forw {this->sysEnergy};

    std::array <int,8> contacts_curr {this->contacts};
    std::array <int,8> contacts_forw {this->contacts};
    std::array <int,8> contacts_back {this->contacts};


    History history_store;

	// decide which half of the polymer is going to regrow
	int growth_dir {-1}; 

	if ( deg_poly % 2 == 0 ){
		growth_dir = (0.5 >= (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
	}
	else {
		if ( 0.5 == (m_index+1)/static_cast<double>(deg_poly+1) ){
			growth_dir = rng_uniform (0, 1);
		}
		else {
			growth_dir = (0.5 > (m_index+1)/static_cast<double>(deg_poly)) ? 0 : 1; 
		}
	}

	if (growth_dir) {
        
        for (int m{m_index+1}; m<deg_poly; ++m){
            history_store.LocH [m] = lattice_index(this->PSETS.Polymers[0].chain[m]->coords, this->y, this->z);
            history_store.SpinH[this->PSETS.Polymers[0].chain[m]].push_back(this->PSETS.Polymers[0].chain[m]->orientation);
        }

		// start the head regrowth...
        for (int m {m_index}; m<deg_poly; ++m){
            this->swing_monomer(m, &history_store, &energy_forw, &contacts_forw, &prob_forw);
            this->kick_orientation(m, &history_store, &energy_forw, &contacts_forw, &prob_forw);
        }
        this->check_structures();

        this->calculate_energetics();
        if (this->sysEnergy != energy_forw || this->contacts != contacts_forw){
            std::cout << "energy_forw = " << energy_forw << " while sysEnergy = " << this->sysEnergy << "." << std::endl;
            std::cout << "contacts_forw = "; print(contacts_forw, ", "); std::cout << "while contacts = "; print(contacts);
            std::cout << "We have a problem." << std::endl;
            exit(EXIT_FAILURE);
        }

        contacts_back = contacts_forw;
        energy_back   = energy_forw;

        // flow backwards...
        for (int m{m_index}; m<deg_poly; ++m){
            this->unswing_monomer(m, &history_store, &energy_back, &contacts_back, &prob_back);
            this->restore_orientation(m, &history_store, &energy_back, &contacts_back, &prob_back);
        }
        this->check_structures(); 

        if (energy_back != energy_curr || contacts_back != contacts_curr){
            std::cout << "Energy after back flow is " << energy_back <<" as opposed to energy_curr = " << energy_curr << std::endl;
            std::cout << "contacts_back = "; print(contacts_back, ", "); std::cout << "while contacts_curr = "; print(contacts_curr);
            std::cout << "We have a problem." << std::endl; 
            exit(EXIT_FAILURE);
        }


        // check for acceptance
        double rng_acc = rng_uniform (0.0, 1.0);
        if ( rng_acc < std::exp(-1/temperature * (energy_forw - energy_curr)) * prob_forw/prob_back){

            this->adopt_new_structure (&history_store); 
            this->sysEnergy = energy_forw; 
            this->contacts  = contacts_forw;

        }
	}

    return;

}

//////////////////////////////////////////////////////////
//   
//                Run the simulation
//
//////////////////////////////////////////////////////////

void Simulation::run() {

	for(int i{this->step_number+1}; i < (this->step_number+this->max_iter+1); ++i){

		// perturb the system!
		// this->perturb_system(); 

		if (this->IMP_BOOL){
			this->acceptances[this->move_number] += 1;
		}

		if ((i % this->dfreq == 0)){
			this->dump_polymers(i);
			this->dump_energy(i); 
			if (i % this->lfreq == 0){
				this->dump_lattice(i); 
			}
		}
		this->IMP_BOOL = true; 

	}

	// dump out the stats, lattice
	this->dump_statistics(this->step_number+this->max_iter);
	this->dump_lattice(this->step_number+this->max_iter);

	return;
}

void Simulation::run_VERBOSE() {

	for(int i{this->step_number+1}; i < (this->step_number+this->max_iter+1); ++i){
		if(this->v && (i%this->dfreq==0)){
            for (int j{0}; j< int((this->PSETS.Polymers[0]).chain.size()); ++j){
                print((this->PSETS.Polymers[0]).chain[j]->coords, ", "); std::cout << "o = " << (this->PSETS.Polymers[0]).chain[j]->orientation << std::endl;
            }
            std::cout << "Contacts = "; print (this->contacts);
            std::cout << "-------------------------- " << std::endl;
            std::cout << "Step number: " << i << "." << std::endl; 
            std::cout << "Executing..." << std::endl << std::endl;
		}

		// perturb the system!
		// this->perturb_system_V(); 

		if (this->IMP_BOOL){
			this->acceptances[this->move_number] += 1;
		}

		if ((i % this->dfreq == 0)){
			this->dump_polymers(i);
			this->dump_energy(i); 
			if (i % this->lfreq == 0){
				this->dump_lattice(i); 
			}
		}
		this->IMP_BOOL = true; 

	}

	// dump out the stats, lattice
	this->dump_statistics(this->step_number+this->max_iter);
	this->dump_lattice(this->step_number+this->max_iter);

	return;
}

