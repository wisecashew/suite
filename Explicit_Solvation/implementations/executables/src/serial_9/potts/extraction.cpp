#include "Potts.hpp"

void Potts::extract_topology(){

	double info;
	std::array <int,7> input_hit_check  = {0, 0, 0, 0, 0, 0, 0};
	std::pair  <std::string, std::string> p;
	std::tuple <std::string, double, double, int, int> t;
	std::string mystring;
	std::string interaction_type;
	std::vector <std::string> contents = extract_content_from_file(this->inp_topology);

	// defining all the regexes... 
	std::regex x ("^x"),             y       ("^y"),       \
	           z ("^z"),             kT      ("^T"),      \
	           Espsp_a ("^Espsp_a"), Espsp_n ("^Espsp_n"), \
	           sp_sp ("^sp-sp"),     eof     ("^END OF FILE");

	for (std::string s: contents){
		// std::cout << "input_hit_check = "; print(input_hit_check);
		if (std::regex_search(s, x)){
			info = NumberExtractor(s); 
			input_hit_check[0] += 1;
			if (input_hit_check [0] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->x = static_cast<int>(info);
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
			this->y = static_cast<int>(info);
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
			this->z = static_cast<int>(info);
			continue;
		}

		else if (std::regex_search(s, kT)){
			info = NumberExtractor(s);
			input_hit_check[3] += 1;
			if (input_hit_check [3] > 1) {
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->T = info;
			continue; 
		}

		else if (std::regex_search(s, sp_sp)){
			interaction_type = trim (split(s, ':')[1]);
			if (interaction_type == "isotropic" || interaction_type == "parallel" || interaction_type == "antiparallel" || interaction_type == "symmetric"){
				p = std::make_pair ("sp", "sp");
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

		else if (std::regex_search (s, Espsp_a)){
			info = NumberExtractor(s);
			input_hit_check[5] += 1;
			if (input_hit_check [5] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[0] = info;
			continue; 
		}

		else if (std::regex_search (s, Espsp_n)){
			info = NumberExtractor(s);
			input_hit_check[6] += 1;
			if (input_hit_check [6] > 1){
				std::cout << "line = " << s << "." << std::endl;
				std::cerr << "You have not provided all of the topology or provided redundant topology, or your inputs might be incorrectly written. Please check. Exiting..." << std::endl;
				exit (EXIT_FAILURE);
			}
			this->energy_surface[1] = info;
			continue; 
		}

		else if (std::regex_search(s, eof)){
			int sum = std::accumulate(input_hit_check.begin(), input_hit_check.end(), 0);
			if (sum != 7){
				std::cerr << "Some input was not provided." << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		}

		else {
			std::cout << s << std::endl;
			std::cerr << "ERROR: There is a nonstandard input provided in topology for a potts simulation file." << std::endl;
			exit(EXIT_FAILURE); 
		}
	}

	p = std::make_pair ("sp", "sp");
	std::get<1>((this->InteractionMap)[p]) = this->energy_surface[0];
	std::get<2>((this->InteractionMap)[p]) = this->energy_surface[1];

    return;
}