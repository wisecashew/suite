#include "Partition.h"

int main (int argc, char** argv) {

	// INSTANTIATE USER INPUT VARIABLES 
	int opt;            // storage variable to hold options from getopt() [line 36] 
	int Nm        {-1}; // degree of polymerization
	int Nmm_max   {-1}; // maximum number of contacts
	double pw     {-1}; // p_omega for contact alignment 
	double T      {-1}; // temperature of simulation
	double E_MM_A {-1}; // aligned monomer-monomer energy
	double E_MM_N {-1}; // boolean for restarts (default: no restarts) 
	double E_MS_A {-1}; // boolean for biased solvation shell 
	double E_MS_N {-1}; // boolean for aligned lattice 
	double alpha  {-1}; // alpha parameter for beta distribution
	double beta   {-1}; // beta parameter for beta distribution
	std::string dumpfile {"__blank__"};

	// define the long options
	static struct option long_options[] {
		{"degree-of-polymerization", required_argument, 0, 'D'},
		{"pw",                       required_argument, 0, 'p'},
		{"temperature",              required_argument, 0, 'T'},
		{"max-contacts",             required_argument, 0, 'N'},
		{"dumpfile",                 required_argument, 0, 'd'},
		{"emma",                     required_argument, 0, 'M'},
		{"emmn",                     required_argument, 0, 'm'},
		{"emsa",                     required_argument, 0, 'S'},
		{"emsn",                     required_argument, 0, 's'},
		{"alpha",                    required_argument, 0, 'a'},
		{"beta",                     required_argument, 0, 'b'},
		{0, 0, 0, 0} // End of options
	};

	int option_index = 0;
	while ((opt = getopt_long(argc, argv, ":D:d:T:N:p:M:m:S:s:a:b:h", long_options, &option_index)) != -1){
		switch(opt) {
			case 'D':
				Nm = std::stoi(optarg);
				break;

			case 'T':
				T = std::stod(optarg);
				break;

			case 'N':
				Nmm_max = std::stoi(optarg);
				break;

			case 'p':
				pw = std::stod(optarg);
				break;

			case 'M':
				E_MM_A = std::stod(optarg);
				break;

			case 'm':
				E_MM_N = std::stod(optarg);
				break;

			case 'S':
				E_MS_A = std::stod(optarg);
				break;

			case 's':
				E_MS_N = std::stod(optarg);
				break;

			case 'a':
				alpha = std::stod(optarg);
				break;

			case 'b':
				beta = std::stod(optarg);
				break;

			case 'd':
				dumpfile = optarg; 
				break;

			case '?':
				std::cout << "ERROR: Unknown option " << static_cast<char>(optopt) << " was provided." << std::endl;
				exit(EXIT_FAILURE);
				break;

			case ':':
				std::cout << "ERROR: Missing arg for " << static_cast <char> (optopt) << "." << std::endl;
				exit(EXIT_FAILURE);
				break;

			case 'h':
				std::cout << "\n" <<
				"This is a partition function calculator for a single chain.\n" <<
				"Author: satyend@princeton.edu\n" <<
				"\n" <<
				"----------------------------------------------------------------------------------------------------------------------------------\n" <<
				"These are all the inputs the software accepts currently. \n" <<
				"help                       [-h]                             (NO ARG)      (NOT REQUIRED) Prints out this message.\n" <<
				"Degree of polymerization   [-D, --degree-of-polymerization] (INTEGER ARG) (REQUIRED)     Enter number of repeating units in your polymer. \n" <<
				"Maximum number of contacts [-N, --max-contacts]             (INTEGER ARG) (REQUIRED)     Enter number of repeating units in your polymer. \n" <<
				"Temperature                [-T, --temperature]              (DOUBLE ARG)  (REQUIRED)     Enter temperature for computation.\n" <<
				"Alignment probability      [-p, --pw]                       (DOUBLE ARG)  (REQUIRED)     Enter alignment probability for interactions.\n" <<
				"Aligned MM interactions    [-M, --emma]                     (DOUBLE ARG)  (REQUIRED)     Enter strength of aligned monomer-monomer interactions.\n" <<
				"Misaligned MM interactions [-m, --emmn]                     (DOUBLE ARG)  (REQUIRED)     Enter strength of misaligned monomer-monomer interactions.\n" <<
				"Aligned MS interactions    [-S, --emsa]                     (DOUBLE ARG)  (REQUIRED)     Enter strength of aligned monomer-solvent interactions.\n" <<
				"Misaligned MS interactions [-s, --emsa]                     (DOUBLE ARG)  (REQUIRED)     Enter strength of misaligned monomer-solvent interactions.\n" <<
				"Alpha parameter            [-a, --alpha]                    (DOUBLE ARG)  (REQUIRED)     Enter alpha parameter for beta distribution.\n" <<
				"Beta parameter             [-b, --beta]                     (DOUBLE ARG)  (REQUIRED)     Enter beta parameter for beta distribution.\n" <<
				"Dump file                  [-d, --dumpfile]                 (STRING ARG)  (REQUIRED)     Enter name of dump file for single chain.\n\n";
				exit(EXIT_SUCCESS);
				break;
		}
	}

	std::cout << "----------------------------------------------------\n";
	std::cout << "Thermodynamic information:\n" << 
	"\tDegree of polymerization is " << Nm << ".\n" << 
	"\tMaximum number of contacts is " << Nmm_max << ".\n" <<
	"\tTemperature of computation is " << T <<".\n" <<
	"\tAlignment probability is " << pw << ".\n" << 
	"\tAligned MM interaction is " << E_MM_A << ".\n" << 
	"\tMisaligned MM interaction is " << E_MM_N << ".\n" << 
	"\tAligned MS interaction is " << E_MS_A << ".\n" << 
	"\tMisaligned MS interaction is " << E_MS_N << ".\n" <<
	"\tAlpha is " << alpha << ".\n" << 
	"\tBeta is " << beta << ".\n" << 
	"----------------------------------------------------"<< std::endl;

	auto start = std::chrono::high_resolution_clock::now();

	// Use a loop to generate values from start to end in increments of step
	Partition my_partition (Nm, Nmm_max, T, pw, E_MM_A, E_MM_N, E_MS_A, E_MS_N, alpha, beta, dumpfile);
	my_partition.get_partition_weights();
	my_partition.get_partition     ();
	my_partition.get_average_energy(); 
	my_partition.get_average_Nmm   ();
	my_partition.get_average_Nmm_a ();
	my_partition.get_average_Nms_a ();
	my_partition.get_Cv();
	my_partition.write ();

	// this is a parallelizer over T
	/*
	#pragma omp parallel for
	for (std::size_t i = 0; i < T_range.size(); ++i){
		Partition my_partition (Nm, T_range[i], pw, E_MM_A, E_MM_N, E_MS_A, E_MS_N, alpha, beta, dumpfile);
		std::cout << "T = " << T_range[i] << "." << std::endl;
		my_partition.get_partition_weights();
		my_partition.get_partition     ();
		my_partition.get_average_energy(); 
		my_partition.get_average_Nmm   ();
		my_partition.get_average_Nmm_a ();
		my_partition.get_average_Nms_a ();
		my_partition.get_Cv();
		my_partition.write ();
	}
	*/
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
	std::cout << "Time for computation is " << duration.count()/1e+6 << " seconds." << std::endl;

	return 0;
}
