#include "Potts.hpp"


void Potts::print_opening_tiles(){

	std::pair  <std::string, std::string> spsp_pair  = std::make_pair("sp", "sp");
	std::cout << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl;
	std::cout << "This is a Potts model simulation." << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;
	std::cout << "Preparing for take-off." << std::endl << std::endl;
	std::cout << "Geometric information about simulation cell: " << std::endl;
	std::cout << "x = " << this->x <<", y = " << this->y << ", z = "<< this->z << "." << std::endl << std::endl;
	std::cout << "Thermodynamic and energetic information about simulation: " << std::endl; 
	std::cout << "Temperature = " << this->T << "." << std::endl;
	std::cout << "Particle-particle energetics: ("   << std::get<0>(this->InteractionMap[spsp_pair])  << ") Espsp_a = "  << this->energy_surface[0] <<", Espsp_n = "  << this->energy_surface[1] << "." << std::endl << std::endl;
	std::cout << "--------------------------------------------------------------------" << std::endl << std::endl;

	return;
}
