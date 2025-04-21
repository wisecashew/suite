#include <iostream>
#include <cmath>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/factorials.hpp>

// Parameters for the beta distribution
const double a = 6.004256058771041;
const double b = 27.319457320842343;
const double loc = 31;
const double scale = loc * (loc-1) / 2;

// count Nms 
int Nms_tot(int Nmm, int Nm){
	return 770 - 2 * (Nmm - (Nm - 1));
}

// get the beta distribution
double p_chain(double Nmm) {
    // Normalize Nmm to the range [0, 1] as required for the Beta distribution
    double x = (Nmm - loc) / scale;

    if (x <= 0.0 || x >= 1.0) {
        return 0.0;  // PDF is 0 outside the range of the Beta distribution
    }

    // Calculate the Beta PDF using Boost's beta function
    double beta_pdf = std::pow(x, a - 1) * std::pow(1 - x, b - 1) / boost::math::beta(a, b);

    // Adjust for the scale of the distribution
    return beta_pdf / scale;
}

// Function to calculate combinations (comb_mm)
double comb_mm(int Nmm, int Nmm_a, double pw) {
    double factorial_Nmm = boost::math::factorial<double>(Nmm);
    double factorial_Nmm_a = boost::math::factorial<double>(Nmm_a);
    double factorial_diff = boost::math::factorial<double>(Nmm - Nmm_a);
    
    return (factorial_Nmm / (factorial_Nmm_a * factorial_diff)) * std::pow(pw, Nmm_a) * std::pow(1 - pw, Nmm - Nmm_a);
}

// Function to calculate combinations (comb_ms)
double comb_ms(int Nmm, int Nms_a, int Nm, double pw) {
    int Nms_tot_value = Nms_tot(Nmm, Nm);
    double factorial_Nms_tot = boost::math::factorial<double>(Nms_tot_value);
    double factorial_Nms_a = boost::math::factorial<double>(Nms_a);
    double factorial_diff = boost::math::factorial<double>(Nms_tot_value - Nms_a);
    
    return (factorial_Nms_tot / (factorial_Nms_a * factorial_diff)) * std::pow(pw, Nms_a) * std::pow(1 - pw, Nms_tot_value - Nms_a);
}



// get the partition function
double partition_function(double T){



}


int main (int argc, char** argv) {


	return 0;
}