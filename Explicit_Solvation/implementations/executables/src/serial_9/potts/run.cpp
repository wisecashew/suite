#include "Potts.hpp"

void Potts::run(){

	// set up the simulation
	this->setup();

	// print out the opening tiles
	this->print_opening_tiles();

	// run the simulation!
	this->simulate();

	// final dumps of lattice and move statistics
	this->dump_lattice();
	this->dump_stats();

	return;
}
