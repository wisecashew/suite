#include <iostream>
#include <vector>
#include <array> 
#include "template.h"

int main() {

    std::array <int,4> a = {1,2,3,4};
    std::vector <double> v = {0.0, 4.0, 2.0 ,2.2}; 
    std::array <double,3> ad = {-1.2, -2, 3.14};

    print (a); 
    print (v); 
    print (ad);

}