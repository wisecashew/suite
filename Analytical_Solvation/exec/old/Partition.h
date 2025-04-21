#ifndef _PARTITION_H_
#define _PARTITION_H_
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <numeric>
#include <fstream>
#include <chrono>
#include <type_traits>
#include <getopt.h>
#include <stdlib.h>
#include <omp.h>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/gamma.hpp>

// Typedef for high-precision floating-point type
typedef boost::multiprecision::cpp_dec_float_50 big_float;

class Partition {
public:
	// properties to be plugged in
	int    Nm;
	int    Nmm_max;
	double T;
	double pw;
	double EMM_A;
	double EMM_N;
	double EMS_A;
	double EMS_N;
	
	// properties that will be evaluated
	int    z;
	int    contact_constraint_constant;
	double Cv;
	double ave_E;
	double ave_Nmm;
	double ave_Nms;
	double ave_Nmm_a;
	double ave_Nms_a;

	// get the beta distribution parameters
	double alpha; // = 6.004256058771041;
	double beta ; // = 27.319457320842343;
	double loc  ; // = 31;
	double scale; // = ((loc+1) * 12 + 1) - (loc);

	// set up the BIG FLOATS!
	big_float Z;
	std::vector<big_float> factorial_cache;
	std::vector<big_float> pw_power_cache;
	std::vector<big_float> ipw_power_cache;


	// define the containers
	std::vector <big_float> weights_store;
	std::vector <big_float> energy_store;
	std::vector <big_float> energy2_store;
	std::vector <big_float> Nmm_store;
	std::vector <big_float> Nmm_a_store;
	std::vector <big_float> Nms_a_store;

	// get the dump file
	std::string dumpfile;

	// constructor
	Partition(int Nm, int Nmm_max,
	double T,
	double pw,
	double EMM_A, double EMM_N,
	double EMS_A, double EMS_N,
	double alpha, double beta,
	std::string dfile):
	Nm (Nm), Nmm_max (Nmm_max), T(T), pw (pw),
	EMM_A (EMM_A), EMM_N (EMM_N),
	EMS_A (EMS_A), EMS_N (EMS_N),
	alpha (alpha), beta (beta),
	dumpfile(dfile) {
		this->z  = 26;
		this->Cv = -1;
		this->ave_Nmm   = -1;
		this->ave_Nmm_a = -1;
		this->ave_Nms_a = -1;
		this->contact_constraint_constant = 2*this->Nm + (z-2) * (this->Nm-2) + (z-1) * (this->Nm);
		this->loc       = this->Nm - 1;
		this->scale     = this->Nmm_max - this->loc;
		std::ofstream f (this->dumpfile, std::ios::out);
		f.close();
		this->init_factorial_cache();
		this->init_pw_power_cache();
		this->init_ipw_power_cache();

		std::cout << "pw_power_cache" << std::endl;
		for (auto&  pw_power: this->pw_power_cache){
			std::cout << pw_power << ", ";
		}
		std::cout << "ipw_power_cache" << std::endl;
		for (auto&  pw_power: this->ipw_power_cache){
			std::cout << pw_power << ", ";
		}


	};

	// destructor
	~Partition(){};

	// methods
	int Nms_tot(int Nmm){
		return this->contact_constraint_constant - 2 * (Nmm - (this->Nm - 1));
	}

	// get the beta distribution
	double p_chain(int Nmm) {
		// Normalize Nmm to the range [0, 1] as required for the Beta distribution
		double x        = (Nmm - this->loc) / this->scale;
		double beta_pdf = 0;
		if (x <= 0 || x >= 1.0) {
			beta_pdf = 0.0;
		}
		else {

		// Compute in log space for numerical stability
		beta_pdf = (alpha - 1.0) * std::log(x) +
						(beta - 1.0) * std::log1p(-x) -
						std::log(boost::math::beta(alpha, beta)) -
						std::log(scale);
		beta_pdf = std::exp(beta_pdf);
		}

		return beta_pdf;
	}

	void init_factorial_cache() {
		std::cout << "Initializing factorial cache." << std::endl;
		int max_n = this->contact_constraint_constant - 2 * (this->Nm - 1);
		this->factorial_cache.resize(max_n+1);
		this->factorial_cache[0] = 1;
		for(int i = 1; i <= max_n; i++) {
			this->factorial_cache[i] = this->factorial_cache[i-1] * i;
		}
		return;
	}

	void init_pw_power_cache(){
		std::cout << "Initializing pw powers cache." << std::endl;
		int max_n = this->contact_constraint_constant - 2 * (this->Nm - 1);
		this->pw_power_cache.resize(max_n+1);
		this->pw_power_cache[0] = 1;
		for(int i = 1; i <= max_n; i++) {
			this->pw_power_cache[i] = this->pw_power_cache[i-1] * (this->pw);
		}
		return;
	}

	void init_ipw_power_cache(){
		std::cout << "Initializing ipw powers cache." << std::endl;
		int max_n = this->contact_constraint_constant - 2 * (this->Nm - 1);
		this->ipw_power_cache.resize(max_n+1);
		this->ipw_power_cache[0] = 1;
		for(int i = 1; i <= max_n; i++) {
			this->ipw_power_cache[i] = this->ipw_power_cache[i-1] * (1-this->pw);
		}
		return;
	}

	// Function to calculate combinations (comb_mm)
	big_float comb_mm(int Nmm, int Nmm_a) {
		big_float factorial_tot   = this->factorial_cache[Nmm]/(this->factorial_cache[Nmm_a] * this->factorial_cache[Nmm-Nmm_a]) * this->pw_power_cache[Nmm_a] * this->ipw_power_cache[Nmm-Nmm_a]; // (factorial_Nmm / (factorial_Nmm_a * factorial_diff)) * boost::multiprecision::pow(this->pw, bf_Nmm_a) * boost::multiprecision::pow(1 - this->pw, bf_Nmm - bf_Nmm_a);
		
		return factorial_tot < 0 ? big_float(0) : factorial_tot;
	}

	// Function to calculate combinations (comb_mm)
	big_float comb_ms(int Nmm, int Nms_a) {
		int Nms = this->Nms_tot(Nmm);
		big_float factorial_tot   = this->factorial_cache[Nms]/(this->factorial_cache[Nms_a] * this->factorial_cache[Nms-Nms_a]) * this->pw_power_cache[Nms_a] * this->ipw_power_cache[Nms-Nms_a]; // (factorial_Nmm / (factorial_Nmm_a * factorial_diff)) * boost::multiprecision::pow(this->pw, bf_Nmm_a) * boost::multiprecision::pow(1 - this->pw, bf_Nmm - bf_Nmm_a);
		std::cout << this->factorial_cache[Nms] << " | " << this->factorial_cache[Nms_a] << " | " << this->factorial_cache[Nms-Nms_a] << " | " << factorial_tot << std::endl;
		// assert(factorial_tot >= 0 && "comb_ms returned negative value");
		return factorial_tot < 0 ? big_float(0) : factorial_tot;
	}

	// get the energy
	double Energy(int Nmm, int Nmm_a, int Nms_a) {
		double energy = Nmm_a * this->EMM_A + (Nmm - Nmm_a) * this->EMM_N + Nms_a * EMS_A + (Nms_tot(Nmm) - Nms_a) * EMS_N;
		return energy;
	}

	// get the boltzmann factor
	big_float boltzmann(int Nmm, int Nmm_a, int Nms_a){
		big_float entropic_factor = this->p_chain(Nmm) * this->comb_mm(Nmm, Nmm_a) * this->comb_ms(Nmm, Nms_a);
		big_float bf_energy       = big_float(-1/this->T * Energy(Nmm, Nmm_a, Nms_a));
		big_float energy_factor   = boost::multiprecision::exp(bf_energy);
		big_float boltz           = entropic_factor * energy_factor;
		
		return boltz < 0 ? big_float(0) : boltz;
	}

	// run the for loop
	void get_partition_weights(){
		
		double    S_chain = 0;
		big_float S_Nmm   = 0;
		big_float S_Nms   = 0;
		big_float boltzmann_factor = 0;
		big_float energy_factor    = 0;
		auto start = std::chrono::high_resolution_clock::now();
		for (int Nmm{this->Nm-1}; Nmm < this->Nmm_max; ++Nmm){
			std::cout << "@ Nmm = " << Nmm << "...";
			S_chain = this->p_chain(Nmm);
			for (int Nmm_a{0}; Nmm_a < Nmm+1; ++Nmm_a){
				S_Nmm = this->comb_mm(Nmm, Nmm_a);
				for (int Nms_a {0}; Nms_a < this->Nms_tot(Nmm)+1; ++Nms_a){
					S_Nms            = this->comb_ms(Nmm, Nms_a);
					energy_factor    = big_float(this->Energy(Nmm, Nmm_a, Nms_a));
					boltzmann_factor = S_chain * S_Nmm * S_Nms * boost::multiprecision::exp(-1/this->T * energy_factor);
					std::cout << "Boltzmann factor = " << boltzmann_factor << std::endl;
					this->weights_store.push_back(boltzmann_factor);
					this->energy_store.push_back (energy_factor * boltzmann_factor);
					this->energy2_store.push_back(energy_factor * energy_factor * boltzmann_factor);
					this->Nmm_store.push_back  (Nmm   * boltzmann_factor);
					this->Nmm_a_store.push_back(Nmm_a * boltzmann_factor);
					this->Nms_a_store.push_back(Nms_a * boltzmann_factor);
				}
			}
			auto stop     = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds> (stop-start);
			std::cout << " done! Time spent: " << static_cast<double>(duration.count())/1e+6 << " seconds." << std::endl;
		}
		return;
	}
	
	// get the partition function
	void get_partition(){
		this->Z = -this->T*boost::multiprecision::log(std::accumulate(this->weights_store.begin(), this->weights_store.end(), big_float(0.0)));
		return;
	}

	// get the average energy
	void get_average_energy(){
		this->ave_E = (std::accumulate(this->energy_store.begin(), this->energy_store.end(), big_float(0.0))/this->Z).convert_to<double>();
		return;
	}

	// get average Nmm contacts
	void get_average_Nmm(){
		this->ave_Nmm = (std::accumulate(this->Nmm_store.begin(), this->Nmm_store.end(), big_float(0.0))/this->Z).convert_to<double>();
		return;
	}

	// get average aligned Nmm contacts
	void get_average_Nmm_a(){
		this->ave_Nmm_a = (std::accumulate(this->Nmm_a_store.begin(), this->Nmm_a_store.end(), big_float(0.0))/this->Z).convert_to<double>();
		return;
	}

	// get averages aligned Nms contacts
	void get_average_Nms_a(){
		this->ave_Nms_a = (std::accumulate(this->Nms_a_store.begin(), this->Nms_a_store.end(), big_float(0.0))/this->Z).convert_to<double>();
		return;
	}

	// get the heat capacity
	void get_Cv(){
		double ave_E2   = (std::accumulate(this->energy2_store.begin(), this->energy2_store.end(), big_float(0.0))/this->Z).convert_to<double>();
		this->Cv = ave_E2 - ave_E * ave_E;
		this->Cv = this->Cv / (this->T*this->T);
		return;
	}
	
	// write down the thing
	void write(){
		std::ofstream outfile(this->dumpfile, std::ios::app);
		outfile << this->T << " | " << this->Z << " | " << this->ave_E << " | " << this->Cv << " | " << this->ave_Nmm << " | " << this->ave_Nmm_a << " | " << this->ave_Nms_a << "\n";
		outfile.close();
		return;
	}

};

#endif
