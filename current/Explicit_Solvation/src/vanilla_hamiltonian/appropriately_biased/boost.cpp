#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <iostream>
#include <cmath>

int main()
{
   // using namespace boost::multiprecision;

    double u = std::exp (900.0);

   boost::multiprecision::cpp_dec_float_50 v = 2.71;

   // loop for e^900
   for(unsigned i = 0; i < 900; ++i){
      v *= 2.71;
   }
   
   boost::multiprecision::cpp_dec_float_50 x = boost::math::expm1<boost::multiprecision::cpp_dec_float_50>(900.0);

   std::cout << "u = " << u << std::endl; 
   std::cout << "v = " << v << std::endl; 
   std::cout << "x = " << x << std::endl;  

   return 0;
}
