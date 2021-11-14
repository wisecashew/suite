#include <iostream> 
#include <vector> 
#include <string>
#include <algorithm> 
#include <map>
#include <random>
#include <chrono>
#include <regex> 
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <stdlib.h>

  int getIntFromString(std::string s){
    std::stringstream ss; 
    ss << s;  // convert the string s into stringstream 
    std::string temp_s; 
    int i; 
    while (!ss.eof()) {
      ss >> temp_s; 
      if (std::stringstream(temp_s) >> i){
        return i; 
      }
    }

    return -1;
  }


int main(int argc, char* argv[])
{

  std::cout << getIntFromString("check check=33") << std::endl;

  std::regex x ("x ="); 
  std::regex y ("y ="); 
  std::regex z ("z ="); 
  std::regex kT ("kT ="); 
  std::regex mat("ENERGY INTERACTION MATRIX");
  std::string myString; 
  std::ifstream file ("energy.txt"); 
  if (file.is_open()) {
    while (file.good()) {
      std::getline(file, myString); 

      if (std::regex_search(myString, x)){
        int x_val = getIntFromString(myString);
        std::cout << "xval is " << x_val << std::endl;
      }
      else if (std::regex_search(myString, y)){
        int y_val = getIntFromString(myString); 
        std::cout << "yval is " << y_val << std::endl;
      }
      else if (std::regex_search(myString, z)){
        int z_val = getIntFromString(myString); 
        std::cout << "zval is " << z_val << std::endl;
      }
      else if (std::regex_search(myString, kT)){
        int kT = getIntFromString(myString); 
        std::cout << "kT is " << kT << std::endl;
        break;
      }


  }

}


  
  return 0;

}