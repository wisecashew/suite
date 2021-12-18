#ifndef _ISING_CLASSES_H_
#define _ISING_CLASSES_H_



class Cluster{
public:
    std::string cluster_name;
    std::vector <ClusterParticle> cparticles; 

    // constructors 
    Cluster(std::string name): cluster_name(name){};
};




class ClusterParticle{
public:
    Particle p; 
    std::string cluster_name = "None"; 


    // constructors 
    ClusterParticle(Particle p_, std::string name_): p(p_), cluster_name(name_) {};
    ClusterParticle(Particle p_): p(p_), cluster_name("None") {}; 
};





#endif 