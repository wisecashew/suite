#include "potts/Potts.hpp"
#include "fhp/FHP.hpp"

int main (int argc, char** argv) {

	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	auto start = std::chrono::high_resolution_clock::now();

	// read the first input
	if (argc > 1){
		std::string sim_style = argv[1];
		if (sim_style == "potts") {
			Potts my_sim (argc, argv);
			my_sim.run();
		}
		else if (sim_style == "fhp"){
			FHP my_sim (argc, argv);
			my_sim.run();
		}
		else {
			std::cout << "The first argument has to be one of \"potts\", \"fhp\"." << std::endl;
			std::cout << "Enter \'-h\' after an appropriate first argument is entered." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	std::cout << "Outside loop." << std::endl;
	auto stop     = std::chrono::high_resolution_clock::now(); 
	auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start); 

	std::cout << "Time required for computation is " << duration.count()/1e+6 << " seconds." << std::endl;
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

	return 0;

}
