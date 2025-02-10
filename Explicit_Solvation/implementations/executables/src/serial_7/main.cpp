#include "Simulation.hpp"

int main (int argc, char** argv) {

	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	auto start = std::chrono::high_resolution_clock::now();
	
	// set up the simulation object 
	Simulation mySim(argc, argv);

	// right now, i have only instantiated certain files and certain pathways. 
	if (mySim.potts){
		// set up the potts simulation
		mySim.set_up_Potts();
	}
	else {
		// set up the fhp simulation
		mySim.set_up_FHP();
	}
	
	// print out certain properties of the simulation
	std::cout << "Positions: "             << mySim.positions          << std::endl;
	std::cout << "Topology: "              << mySim.topology           << std::endl;
	std::cout << "Coord dump file: "       << mySim.dfile              << std::endl;
	std::cout << "Energy dump file: "      << mySim.efile              << std::endl;
	std::cout << "Orientation dump file: " << mySim.mfile              << std::endl;
	std::cout << "Stats file: "            << mySim.stats_file         << std::endl;
	std::cout << "Lattice file to write: " << mySim.lattice_file_write << std::endl;
	std::cout << "Lattice file to read: "  << mySim.lattice_file_read  << std::endl;
	std::cout << "Solvation shell file: "  << mySim.SSfile             << std::endl;

	// print out the opening tiles
	mySim.opening_tiles();

	// run the simulation!
	mySim.run();

	// final dumps of lattice and move statistics
	mySim.dump_lattice_end();
	mySim.dump_statistics();
	
	auto stop     = std::chrono::high_resolution_clock::now(); 
	auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 

	std::cout << "Time required for computation is " << duration.count()/1e+6 << " seconds." << std::endl;
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

	return 0;

}
