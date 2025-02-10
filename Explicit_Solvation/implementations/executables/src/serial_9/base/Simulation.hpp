#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include "Polymer.hpp"

class Simulation {
public:

	/*~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
	//
	//			Properties defintions block
	//
	//~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#*/ 
	// to first order, let's define certain properties 
	// class properties are organized according to their data type
	// as in, all the `int`s come first, then `double`s, then `bool`s, and so on...
	int x;                // length along x-axis of box. 
	int y;                // length along y-axis of box. 
	int z;                // length along z-axis of box. 
	int max_iter;         // number of moves to perform in the simulation.
	int dump_freq;        // frequency of dumps to energydump and coords file. 
	int lattice_dfreq;    // frequency of dumping out the entire lattice.
	int step_number {0};  // current step number of the simulation.
	int move_number {-1}; // number of last selected move.
	// end of `int` properties

    // global, thermodynamic properties of the system
	double T         {0}; // temperature of the system
	double Energy    {0}; // energy of the system
	// end of `double` properties

	bool v;          // if true, run the debugger
	bool r;          // if true, the simulation is being restarted
	bool A;          // if true, the entire lattice has orientation 0.
	bool acceptance; // if true, the suggested perturbation has been accepted.  
	// end of `bool` properties

	std::string inp_topology;           // input; name of file with topology (energetics + geometry) of system 
    std::string inp_lattice_file_read;  // input; name of file from which lattice will be read 
	std::string out_energy_dump;        // output; name of energy dump file 
	std::string out_stats_dump;         // output; name of file with move statisitics 
	std::string out_lattice_file_write; // output; name of file where lattice will be dumped to 
	// end of `std::string` properties

	// define the entire set of particles that we are interested in 
	std::vector <Particle*> Lattice;

    // virtual destructor for proper cleanup in derived classes
    virtual ~Simulation() = default;	

};

#endif
