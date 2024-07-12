#include "Simulation.h"

//////////////////////////////////////////////////////////
//
//                Dump functions
//
//////////////////////////////////////////////////////////

void Simulation::dump_energy(){

	std::ofstream dump_file(this->efile, std::ios::app); 

	dump_file << sysEnergy << " | " \
			<< (this->contacts)[0]+(this->contacts)[1] << " | " << (this->contacts)[0] << " | " << (this->contacts)[1] << " | " \
			<< (this->contacts)[2]+(this->contacts)[3] << " | " << (this->contacts)[2] << " | " << (this->contacts)[3] << " | " \
			<< (this->contacts)[4]+(this->contacts)[5] << " | " << (this->contacts)[4] << " | " << (this->contacts)[5] << " | " \
			<< (this->contacts)[6]+(this->contacts)[7] << " | " << (this->contacts)[6] << " | " << (this->contacts)[7] << " | " << this->step_number << "\n";

	return; 
}

void Simulation::dump_polymers(){
	std::ofstream dump_file(this->dfile, std::ios::app); 
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

void Simulation::dump_solvation_shell_orientations(){
	
	std::ofstream dump_file (this->mfile, std::ios::app); 
	
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

void Simulation::dump_solvation_shell(){

	std::ofstream dump_file(this->SSfile, std::ios::app);
	std::array<std::array<int,3>,26> ne_list;
	std::array<double,8> contacts = {0,0,0,0,0,0,0,0};
	std::array<int,3> stats = {0, 0, 0}; // {aligned, misaligned, total number of particles}
	std::set  <int> solvation_shell_set; 

	// get the first solvation shell 
	// auto start = std::chrono::high_resolution_clock::now(); 
	for ( Polymer& pmer: this->Polymers){
		for (Particle*& p: pmer.chain){
		ne_list = obtain_ne_list (p->coords, this->x, this->y, this->z); 
			for ( std::array <int,3>& loc: ne_list ){
				if ( this->Lattice[lattice_index (loc, y, z)]->ptype[0] == 's' ){
					solvation_shell_set.insert (lattice_index (loc, y, z)); 
				}
			}
		}
	}

	stats[2] = solvation_shell_set.size(); 
	double energy{0}; 

	for (const int loc: solvation_shell_set){
		ne_list = obtain_ne_list(location(loc, this->x, this->y, this->z), this->x, this->y, this->z);
		for (std::array <int,3>& ne: ne_list){
			if (this->Lattice[lattice_index(ne, this->y, this->z)]->ptype[0] != 'm' && this->Lattice[lattice_index(ne, this->y, this->z)]->ptype != this->Lattice[loc]->ptype){
				auto func = this->PairwiseFunctionMap.find({"s1", "s2"});
				func->second(this, Lattice[loc], Lattice[lattice_index(ne, this->y, this->z)], &contacts, &energy);
				auto it   = solvation_shell_set.find(lattice_index(ne, this->y, this->z));
				if (it != solvation_shell_set.end()){
					stats[0] += contacts[6]; // outside solvation shell
					stats[1] += contacts[7]; // outside solvation shell
					contacts[6] = 0;
					contacts[7] = 0;
				}
				else {
					stats[0] += 0.5 * contacts[6]; // inside solvation shell
					stats[1] += 0.5 * contacts[7]; // inside solvation shell
					contacts[6] = 0;
					contacts[7] = 0;
				}
			}
		}
	}
	dump_file << stats[0] << " | " << stats[1] << " | " << stats[2] << " | " << this->step_number << "\n";
	return; 

}

void Simulation::dump_statistics(){

	std::ofstream dump_file (this->stats_file, std::ios::out); 
	dump_file << "At step " << this->step_number << "...\n";

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

void Simulation::dump_local_all(){
	if ((this->step_number % this->dfreq) == 0){
		this->dump_energy();
		this->dump_polymers();
		this->dump_solvation_shell();
		this->dump_solvation_shell_orientations();
	}
	return;
}

void Simulation::dump_local_no_ss(){
	// std::cout << "step_number % this->dfreq = " << step_number % this->dfreq << std::endl;
	if ((this->step_number % this->dfreq) == 0){
		this->dump_energy();
		this->dump_polymers();
		this->dump_solvation_shell_orientations();
	}
	return;
}

void Simulation::dump_lattice(){

	if ((this->step_number % this->lfreq) == 0){
		std::cout << "dump lattice @ step number = " << this->step_number << ", lfreq = " << this->lfreq << std::endl;
		std::ofstream dump_file (this->lattice_file_write, std::ios::out); 
		dump_file << "FINAL STEP: " << this->step_number << ".\n"; 
		for ( Particle*& p: this->Lattice ){
			dump_file << p->orientation << ", " << p->ptype << ", " << lattice_index(p->coords, y, z) << "\n"; 
		}
		dump_file << "END. \n";
	}
	return; 
}
