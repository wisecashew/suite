#ifndef _CREATION_H_
#define _CREATION_H_
#include "particlesets.h"

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 

//                 DEFINITIONS FOR CLASS Simulation

/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 


class Simulation {
public:
	// to first order, let's define certain properties 
	const int x;
	const int y;
	const int z;
	const int max_iter;
	const int dfreq;
	const int lfreq;
	const int Npoly;
	const double T;
	const double frac;
	const bool IMP_BOOL;
	const bool r;
	const bool s;
	const bool v;
	const bool A; 
	const bool S;
	const std::array<double,8> E;
	const std::string positions;          // name of file with initial coords of polymer 
	const std::string topology;           // name of file with topology of system 
	const std::string dfile;              // name of coordinate dump file 
	const std::string efile;              // name of energy dump file 
	const std::string mfile;              // name of orientation dump file 
	const std::string stats_file;         // name of file with move statisitics 
	const std::string lattice_file_write; // name of file where lattice will be dumped to 
	const std::string lattice_file_read ; // name of file from which lattice will be read 

    int    step_number;
    int    move_number;  
    double sysEnergy;

    ParticleSets PSETS; 

    std::array <int,9>    attempts    = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::array <int,9>    acceptances = {0, 0, 0, 0, 0, 0, 0, 0, 0}; 
    std::array <double,8> contacts;

    std::map <std::pair<std::string, std::string>, std::tuple<std::string, double, double, int, int>> InteractionMap; 

    // constructor
	Simulation(int max_iter, int dfreq, int lfreq, 
           bool r, bool s, bool v, bool A, bool S, const std::string& positions,
           const std::string& topology, const std::string& dfile, const std::string& efile,
           const std::string& mfile, const std::string& stats_file, const std::string& lattice_file_write,
           const std::string& lattice_file_read)
    : max_iter(max_iter), dfreq(dfreq), lfreq(lfreq), 
      r(r), s(s), v(v), positions(positions), topology(topology), dfile(dfile), efile(efile),
      mfile(mfile), stats_file(stats_file), lattice_file_write(lattice_file_write),
      lattice_file_read(lattice_file_read) {

      	// this->Npoly = ExtractNumberOfPolymers(this->positions);
      	this->extract_number_of_polymers();
      	this->extract_topology_from_file(); 
      	this->set_up_lattice(); 
      	this->calculate_energy();

      };

    // destructor 
    ~Simulation(){};

    // extraction methods
    void extract_number_of_polymers(); 
    void extract_topology_from_file();

    // lattice set up methods
    void set_up_lattice();

    // calculate energy
    void pair_interaction(Particle*);
    void particle_pair_energy_contribution(Particle*);
    void calculate_energy();

    // verification methods
    void check_structures();

    // dump positions
    void dump_energy  (int);
    void dump_polymers(int);
    void dump_lattice (int);

    // perturbation methods


    // execution methods
    void run(); 
    void run_VERBOSE();
    void run_DEBUG();

};


void Simulation::extract_topology_from_file(){

    double info; 
    std::array <double, 13> info_vec; 
    std::array <int,17> input_hit_check = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::pair  <std::string, std::string> p; 
    std::tuple <std::string, double, double, int, int> t; 
    std::string interaction_type; 
    std::string mystring; 
    std::vector <std::string> contents = ExtractContentFromFile(this->topology); 
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
    		this->x = x;
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
    		this->y = y;
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
    		this->z = z;
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
    		this->T = T;
    		continue; 
    	}

    	else if (std::regex_search(s, m1_m1)){
    		interaction_type = trim (split(s, ':')[1]);
    		if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel"){
    			p = std::make_pair ("m1", "m1");
    			t = std::make_tuple(interaction_type, 0, 0, 0, 1);
    			(*InteractionMap)[p] = t;
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
    			(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 2, 3);
    			p = std::make_pair ("s1", "m1") ;
    			(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 2, 3);
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
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 4, 5);
	    		p = std::make_pair ("s2", "m1");
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 4, 5);	
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
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 6, 7);
	    		p = std::make_pair ("s2", "s1");
	    		(*InteractionMap)[p] = std::make_tuple(interaction_type, 0, 0, 6, 7);	
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

	if (!this->r){

		Polymers = ExtractPolymersFromFile (this->positions, this->x, this->y, this->z);
		AddSolvent(&Lattice, this->x, this->y, this->z);

		// populate the lattice
		for (Polymer& pmer: (Polymers)) {
			for (Particle*& p: pmer.chain){
				// now that I have my polymer coordinates, time to create the grand lattice 
				(LATTICE).at(lattice_index (p->coords, y, z) ) = p; 
			}
		}

		AddCosolvent (&Polymers, &Cosolvent, &Lattice, this->s, this->frac, this->x, this->y, this->z); 

		if (this->A){
			AlignTheLattice(&Lattice);
		}
		else if (this->S){
			AlignTheSolvationShell(&Polymers, &Lattice, this->x, this->y, this->z);
		}

	}

	else {
		SetUpLatticeFromRestart(&Polymers, &Cosolvent, &Lattice, this->step_number, this->lattice_file_read, this->dfile, this->positions, this->x, this->y, this->z);
	}

	this->PSETS.Lattice   = Lattice;
	this->PSETS.Polymers  = Polymers;
	this->PSETS.Cosolvent = Cosolvent;
	this->check_structures();

}


void Simulation::run() {

	for(int i{this->step_number+1}; i < (this->step_number+this->max_iter+1); ++i){

		// perturb the system!
		this->PerturbSystem(); 

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

	// dump out the lattice
	this->dump_lattice(this->step_number+this->max_iter);
	return;
}


void Simulation::run_VERBOSE() {

	for(int i{this->step_number+1}; i < (this->step_number+this->max_iter+1); ++i){
		if(this->v && (i%this->dfreq==0)){
            for (int j{0}; j< int((Polymers)[0].chain.size()); ++j){
                print((this->LATTICE[0])[0].chain[j]->coords, ", "); std::cout << "o = " << (this->LATTICE[0])[0].chain[j]->orientation << std::endl;
            }
            std::cout << "Contacts = "; print (this->contacts);
            std::cout << "-------------------------- " << std::endl;
            std::cout << "Step number: " << i << "." << std::endl; 
            std::cout << "Executing..." << std::endl << std::endl;
		}

		// perturb the system!
		this->PerturbSystem(); 

		if (this->IMP_BOOL){
			this->acceptances[this->move_number] += 1;
		}

		if ((i % this->dfreq == 0)){
			this->dump_positions_of_polymers(i);
			this->dump_energy(i); 
			if (i % this->lfreq == 0){
				this->dump_lattice(i); 
			}
		}
		this->IMP_BOOL = true; 

	}

	// dump out the lattice
	this->dump_lattice(this->step_number+this->max_iter);

	return;
}

#endif



