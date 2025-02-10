#include "FHP.hpp"

//////////////////////////////////////////////////////////
//
//                Dump functions
//
//////////////////////////////////////////////////////////

void FHP::dump_energy(){

	std::ofstream dump_file(this->out_energy_dump, std::ios::app); 

	dump_file << this->Energy << " | " \
			<< (this->contacts)[0]+(this->contacts)[1] << " | " << (this->contacts)[0] << " | " << (this->contacts)[1] << " | " \
			<< (this->contacts)[2]+(this->contacts)[3] << " | " << (this->contacts)[2] << " | " << (this->contacts)[3] << " | " \
			<< (this->contacts)[4]+(this->contacts)[5] << " | " << (this->contacts)[4] << " | " << (this->contacts)[5] << " | " \
			<< (this->contacts)[6]+(this->contacts)[7] << " | " << (this->contacts)[6] << " | " << (this->contacts)[7] << " | " << this->step_number << "\n";

	return; 
}

void FHP::dump_polymers(){
	std::ofstream dump_file(this->out_coord_dump, std::ios::app); 
	dump_file <<"Dumping coordinates at step " << this->step_number << ".\n";
	
	int count = 0; 
	for (Polymer& pmer: this->Polymers){
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

void FHP::dump_solvation_shell_orientations(){
	
	std::ofstream dump_file (this->out_orientation_dump, std::ios::app); 
	
	dump_file << "START for Step " << this->step_number << ".\n";
	std::vector <int> solvent_indices; 
	
	for ( Polymer& pmer: this->Polymers ) {
		for ( Particle*& p: pmer.chain ) {
			
			dump_file << p->orientation << " | ";
			std::array <std::array<int,3>, 26> ne_list = obtain_ne_list (p->coords, this->x, this->y, this->z);
			
			for ( std::array <int,3>& ne: ne_list) {
				if (this->Lattice[lattice_index(ne, this->y, this->z)]->ptype[0] == 's' && std::find(solvent_indices.begin(), solvent_indices.end(), lattice_index(ne, this->y, this->z)) == solvent_indices.end()){
					dump_file << (this->Lattice[lattice_index(ne, this->y, this->z)])->orientation << " | ";  
				} 
			}
			dump_file << "\n"; 
		} 
	}
	
	dump_file << "END. \n";
	return;
}

void FHP::dump_solvation_shell(){

	std::ofstream dump_file(this->out_solvation_shell_dump, std::ios::app);
	double energy{0};

	// define set for solvent and cosolvent indices
	std::set  <int>      solvent_indices;
	std::set  <int>      cosolvent_indices;

	// define some stores
	std::array<double,2> stats    = {0,0};
	std::array<double,CONTACT_SIZE_FHP> contacts = {0,0,0,0,0,0,0,0};

	// define a store neighboring particles
	std::array<std::array<int,3>,26> ne_list;

	// get the first solvation shell 
	for ( Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
		ne_list = obtain_ne_list (p->coords, this->x, this->y, this->z); 
			for ( std::array <int,3>& loc: ne_list ){
				if ( this->Lattice[lattice_index (loc, y, z)]->ptype == "s1" ){
					solvent_indices.insert (lattice_index (loc, y, z));
				}
				else if ( this->Lattice[lattice_index(loc, this->y, this->z)]->ptype == "s2" ){
					cosolvent_indices.insert (lattice_index (loc, y, z));
				}
			}
		}
	}

	// get the total solvent and cosolvent particles
	int total_solvent_particles   {static_cast<int>(  solvent_indices.size())};
	int total_cosolvent_particles {static_cast<int>(cosolvent_indices.size())};
	std::pair <std::string, std::string> particle_pair = std::make_pair ("s1", "s2");

	auto func = this->InitialInteractionFunctionMap.find(particle_pair);

	for (int s_idx: solvent_indices) {
		// get the neighbors
		ne_list = obtain_ne_list ((this->Lattice)[s_idx]->coords, this->x, this->y, this->z);
		for (std::array<int,3>& ne: ne_list){
			if ((this->Lattice)[lattice_index(ne, this->y, this->z)]->ptype == "s2"){
				func->second(this, Lattice[s_idx], Lattice[lattice_index(ne, this->y, this->z)], contacts, energy);
				auto it = solvent_indices.find(lattice_index(ne, this->y, this->z));
				if (it == solvent_indices.end()){
					stats[0]   += contacts[6]; // outside solvation shell, aligned
					stats[1]   += contacts[7]; // outside solvation shell, misaligned
					contacts[6] = 0;
					contacts[7] = 0;
				}
				else{
					stats[0]   += 0.5 * contacts[6]; // inside solvation shell, aligned
					stats[1]   += 0.5 * contacts[7]; // inside solvation shell, misaligned
					contacts[6] = 0;
					contacts[7] = 0;
				}
			}
		}
	}

	for (int c_idx: cosolvent_indices){
		// get the neighbors
		ne_list = obtain_ne_list(this->Lattice[c_idx]->coords, this->x, this->y, this->z);
		for (std::array<int,3>& ne: ne_list){
			auto it = solvent_indices.find(lattice_index(ne, y, z)); // check if the neighbor is within the solvation shell
			// check if the particle 
			if(this->Lattice[lattice_index(ne, this->y, this->z)]->ptype == "s1" && it == solvent_indices.end()){
				stats[0]   += contacts[6];
				stats[1]   += contacts[7];
				contacts[6] = 0;
				contacts[7] = 0;
			}
			else{
				contacts[6] = 0;
				contacts[7] = 0;
			}
		}
	}

	dump_file << total_solvent_particles + total_cosolvent_particles << " | " << total_solvent_particles \
	<< " | " << total_cosolvent_particles << " | " << stats[0] << " | " << stats[1] << " | " << this->step_number << "\n";

	return; 

}


void FHP::dump_stats_all(){

	std::ofstream dump_file (this->out_stats_dump, std::ios::out); 
	dump_file << "At step " << this->step_number << "...\n";

	dump_file << "End rotations without bias         - attempts: " << (this->attempts)[0] <<", acceptances: " << (this->acceptances)[0] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[0])/static_cast<double>((this->attempts)[0]) << std::endl; 
	dump_file << "Reptation without bias             - attempts: " << (this->attempts)[1] <<", acceptances: " << (this->acceptances)[1] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[1])/static_cast<double>((this->attempts)[1]) << std::endl; 
	dump_file << "Regrowth                           - attempts: " << (this->attempts)[2] <<", acceptances: " << (this->acceptances)[2] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[2])/static_cast<double>((this->attempts)[2]) << std::endl;
	dump_file << "Polymer flips                      - attempts: " << (this->attempts)[3] <<", acceptances: " << (this->acceptances)[3] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[3])/static_cast<double>((this->attempts)[3]) << std::endl;
	dump_file << "Lattice flip                       - attempts: " << (this->attempts)[4] <<", acceptances: " << (this->acceptances)[4] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[4])/static_cast<double>((this->attempts)[4]) << std::endl;
	dump_file << "Solvation shell flip with bias     - attempts: " << (this->attempts)[5] <<", acceptances: " << (this->acceptances)[5] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[5])/static_cast<double>((this->attempts)[5]) << std::endl;
	dump_file << "Solvation shell exchange           - attempts: " << (this->attempts)[6] <<", acceptances: " << (this->acceptances)[6] << "." << std::endl;
	dump_file << "Solvent exchange                   - attempts: " << (this->attempts)[7] <<", acceptances: " << (this->acceptances)[7] << "." << std::endl;

	return;
}

void FHP::dump_stats_isotropic(){

	std::ofstream dump_file (this->out_stats_dump, std::ios::out); 
	dump_file << "At step " << this->step_number << "...\n";

	dump_file << "End rotations without bias         - attempts: " << (this->attempts)[0] <<", acceptances: " << (this->acceptances)[0] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[0])/static_cast<double>((this->attempts)[0]) << std::endl; 
	dump_file << "Reptation without bias             - attempts: " << (this->attempts)[1] <<", acceptances: " << (this->acceptances)[1] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[1])/static_cast<double>((this->attempts)[1]) << std::endl; 
	dump_file << "Regrowth                           - attempts: " << (this->attempts)[2] <<", acceptances: " << (this->acceptances)[2] << ", acceptance fraction: " << static_cast<double>((this->acceptances)[2])/static_cast<double>((this->attempts)[2]) << std::endl;
	dump_file << "Solvation shell exchange           - attempts: " << (this->attempts)[3] <<", acceptances: " << (this->acceptances)[3] << "." << std::endl;
	dump_file << "Solvent exchange                   - attempts: " << (this->attempts)[4] <<", acceptances: " << (this->acceptances)[4] << "." << std::endl;

	return;

}

void FHP::dump_stats(){
	(this->*dump_stats_ptr)();
	return;
}

void FHP::dump_local_all(){
	if ((this->step_number % this->dump_freq) == 0){
		this->dump_energy();
		this->dump_polymers();
		this->dump_solvation_shell();
		this->dump_solvation_shell_orientations();
		this->dump_stats();
	}

	if ((this->step_number % this->lattice_dfreq) == 0){
		this->dump_lattice();
	}

	return;
}

void FHP::dump_local_no_ss(){
	// std::cout << "step_number % this->dfreq = " << step_number % this->dfreq << std::endl;
	if ((this->step_number % this->dump_freq) == 0){
		this->dump_energy();
		this->dump_polymers();
		this->dump_solvation_shell_orientations();
		this->dump_stats();
	}

	if ((this->step_number % this->lattice_dfreq) == 0){
		this->dump_lattice();
	}

	return;
}

void FHP::dump_lattice(){

	if ((this->step_number % this->lattice_dfreq) == 0){
		std::ofstream dump_file (this->out_lattice_file_write, std::ios::out); 
		dump_file << "FINAL STEP: " << this->step_number << ".\n"; 
		for ( Particle*& p: this->Lattice ){
			dump_file << p->orientation << ", " << p->ptype << ", " << lattice_index(p->coords, y, z) << "\n"; 
		}
		dump_file << "END. \n";
	}
	return; 
}

void FHP::dump_lattice_end(){

	std::ofstream dump_file (this->out_lattice_file_write, std::ios::out); 
	dump_file << "FINAL STEP: " << this->step_number << ".\n"; 
	for ( Particle*& p: this->Lattice ){
		dump_file << p->orientation << ", " << p->ptype << ", " << lattice_index(p->coords, y, z) << "\n"; 
	}
	dump_file << "END. \n";

	return; 
}

void FHP::dump_opening_tiles(){

	std::pair  <std::string, std::string> mm_pair    = std::make_pair("m1", "m1");
	std::pair  <std::string, std::string> m1s1_pair  = std::make_pair("m1", "s1");
	std::pair  <std::string, std::string> m1s2_pair  = std::make_pair("m1", "s2");
	std::pair  <std::string, std::string> s1s2_pair  = std::make_pair("s1", "s2");

	std::cout << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cout << "This is a Flory-Huggins-Potts simulation." << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
	std::cout << "Preparing for take-off." << std::endl << std::endl;
	std::cout << "Geometric information about simulation cell: " << std::endl;
	std::cout << "x = " << this->x <<", y = " << this->y << ", z = "<< this->z << "." << std::endl << std::endl;
	std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
	std::cout << "Temperature = " << this->T << ", and cosolvent fraction x_c = " << this->frac_c << "."  << std::endl;
	std::cout << "Particle-particle energetics: \n(" << std::get<0>(this->InteractionMap[mm_pair])  << ") Emm_a = "  
	<< this->energy_surface[0] <<", Emm_n = "  << this->energy_surface[1] << std::endl; // "." << std::endl << std::endl;
	std::cout << "(" << std::get<0>(this->InteractionMap[m1s1_pair])  << ") Ems1_a = "  
	<< this->energy_surface[2] <<", Ems1_n = "  << this->energy_surface[3] << std::endl;
	std::cout << "(" << std::get<0>(this->InteractionMap[m1s2_pair])  << ") Ems2_a = "  
	<< this->energy_surface[4] <<", Ems2_n = "  << this->energy_surface[5] << std::endl;
	std::cout << "(" << std::get<0>(this->InteractionMap[s1s2_pair])  << ") Es1s2_a = "  
	<< this->energy_surface[6] <<", Es1s2_n = "  << this->energy_surface[7] << std::endl;
	std::cout << "Initial energy of the system is " << this->Energy << "."  << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
	return;
}

std::set <int> FHP::get_solvation_shell(){
	std::set   <int> solvation_shell_set; 
	std::array <std::array<int,3>,26> ne_list; 
	// int dop = static_cast<int>((*Polymers)[0].chain.size() ); 

	// get the first solvation shell 
	// auto start = std::chrono::high_resolution_clock::now(); 
	for (Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
			ne_list = obtain_ne_list ( p->coords, this->x, this->y, this->z); 
			for (std::array <int,3>& loc: ne_list){
				if (this->Lattice[lattice_index (loc, this->y, this->z)]->ptype[0] == 's'){
					solvation_shell_set.insert (lattice_index (loc, this->y, this->z)); 
				}
			}
		}
	}
	return solvation_shell_set;
}
