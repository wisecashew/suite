#include <boost/multiprecision/cpp_dec_float.hpp>
#include <iostream>
#include <array>

int main()
{
   // using namespace boost::multiprecision;

   boost::multiprecision::cpp_dec_float_50 u = 10, v = 11;
   for(unsigned i = 1; i <= 1000; ++i){
      u *= 10;
   }
   for(unsigned i = 1; i <= 1000; ++i){
      v *= 11;
   }
   
        

   // prints 93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000 (i.e. 100!)
   std::cout << u << std::endl;
   std::cout << v << std::endl;

   std::cout << u/v << std::endl;

   std::array <boost::multiprecision::cpp_dec_float_50,3> arr = {u, v, u/v};
   std::cout << arr[0] << ", " << arr[1] << ", " << arr[2] << "." << std::endl;


   return 0;
}
