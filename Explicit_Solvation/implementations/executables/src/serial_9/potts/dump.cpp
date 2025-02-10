#include "Potts.hpp"

void Potts::dump_stats(){

	std::ofstream dump_file (this->out_stats_dump, std::ios::out); 
	dump_file << "At step " << this->step_number << "...\n";
	dump_file << "Lattice flip - attempts: " << this->n_attempts << ", acceptances: " << this->n_acceptances << ", acc fraction = " << static_cast<double>(this->n_acceptances)/static_cast<double>(this->n_attempts) << "." << std::endl;
	dump_file.close();
	return;

}

void Potts::dump_energy(){
	
	std::ofstream edump_file(this->out_energy_dump, std::ios::app);
		edump_file << this->Energy << " | " \
		<< (this->contacts)[0]+(this->contacts)[1] << " | " << (this->contacts)[0] << " | " << (this->contacts)[1] << " | " << this->step_number << "\n";
		edump_file.close();
	
	return;
}

void Potts::dump_lattice(){

	std::ofstream dump_file (this->out_lattice_file_write, std::ios::out); 
	dump_file << "STEP: " << this->step_number << ".\n"; 
	for ( Particle*& p: this->Lattice ){
		dump_file << p->orientation << ", " << p->ptype << ", " << lattice_index(p->coords, y, z) << "\n"; 
	}
	dump_file << "END. \n";
	dump_file.close();
	return;

}

void Potts::dump(){

	if ((this->step_number % this->dump_freq) == 0){
		this->dump_stats();
		this->dump_energy();
	}

	if ((this->step_number % this->lattice_dfreq) == 0) {
		this->dump_lattice();
	}

	return; 

}
