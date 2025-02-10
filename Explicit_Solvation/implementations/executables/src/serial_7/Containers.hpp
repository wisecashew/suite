#ifndef _CONTAINERS_H_
#define _CONTAINERS_H_
#include "misc.hpp"

class Container {
public:

	// properties
	double                E1_initial;
	double                E1_final;
	double                E2_initial;
	double                E2_final;
	double                Epair_initial;
	double                Epair_final;
	std::array <double,CONTACT_SIZE> c1_initial;
	std::array <double,CONTACT_SIZE> c1_final;
	std::array <double,CONTACT_SIZE> c2_initial;
	std::array <double,CONTACT_SIZE> c2_final;
	std::array <double,CONTACT_SIZE> cpair_initial;
	std::array <double,CONTACT_SIZE> cpair_final;

	// constructor
	Container(){
		this->E1_initial    = 0;
		this->E1_final      = 0;
		this->E2_initial    = 0;
		this->E2_final      = 0;
		this->Epair_initial = 0;
		this->Epair_final   = 0;
		this->c1_initial.fill(0);
		this->c1_final.fill(0);
		this->c2_initial.fill(0);
		this->c2_final.fill(0);
		this->cpair_initial.fill(0);
		this->cpair_final.fill(0);
	};

	// destructor
	~Container(){};

	void reset(){
		this->E1_initial    = 0;
		this->E1_final      = 0;
		this->E2_initial    = 0;
		this->E2_final      = 0;
		this->Epair_initial = 0;
		this->Epair_final   = 0;
		this->c1_initial.fill(0);
		this->c1_final.fill(0);
		this->c2_initial.fill(0);
		this->c2_final.fill(0);
		this->cpair_initial.fill(0);
		this->cpair_final.fill(0);
		return;
	}

};

#endif
