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



int main(int argc, char* argv[]){

    std::ifstream myfile ("energy.txt"); 
    std::string mystring; 
    std::vector <std::string> contents; 

    if (myfile.is_open() ){
        while (myfile.good()) {
            std::getline(myfile, mystring); // pipe file's content into stream 
            contents.push_back(mystring); 
        }
    }

for (auto input: contents){
  std::stringstream ss;
  ss << input;
  int found;
  std::string temp;

  while(std::getline(ss, temp,' ')) {
    if(std::stringstream(temp)>>found)
      {
        std::cout<<found<<std::endl;
      }
    }
  }
return 0;
}