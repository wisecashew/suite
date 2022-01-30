#include <iostream> 
#include <vector>
#include <string> 
#include <iterator>
#include <map>
#include <algorithm> 
#include <regex>
#include <fstream>
#include <sstream>
#include <regex>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <random>
#include <array> 


int main() {

    int numberofpoly = 1; 
    std::vector <int> step_num_store; 
    std::vector <int> index_store; 
    std::ifstream myfile("coords.txt"); 

    std::string mystring; 
    std::vector <std::string> contents; 

    if (myfile.is_open() ) {
        while (std::getline(myfile, mystring) ) {
            contents.push_back(mystring); 
        }
    }

    std::regex stepnum ("Dumping coordinates at step"); 
    // int step_num = 0; 
    int j {0}; 
    for (std::string& s: contents){

        std::stringstream ss(s);
        std::string temp; 
        int found; 
        std::stringstream int_ss; 

        if (std::regex_search(s, stepnum)){
            while(!ss.eof()){
                ss >> temp;
                if (std::stringstream(temp) >> found){
                    step_num_store.push_back(found);
                    index_store.push_back(j);
                } 
            }
        } 
        ++j; 
    }

    //for (int f: step_num_store){
    //    std::cout << f << " | "; 
    // }
    // std::cout << std::endl;

    std::regex final_step ("Dumping coordinates at step " + std::to_string(step_num_store[step_num_store.size()-1]) );

    contents.erase(contents.begin(), contents.begin() + (index_store[index_store.size()-1] ) ); 
    std::cout << "Edited contents are: " << std::endl;
    // for (auto& s: contents){
    //     std::cout << s << std::endl;
    // }

    // I have the critical contents in here. 
    std::vector <std::array <int,3> > locations; 
    std::regex start_coords ("START"), end_coords ("END"); 
    bool start_bool {false}, end_bool {false}; 
    for (auto& s: contents){

        std::stringstream ss(s); 
        
        if (std::regex_search(s, start_coords)){
            start_bool = true; 
            end_bool = false; 
            continue; 
        }

        else if (std::regex_search(s, end_coords)){
            start_bool = false; 
            end_bool = false; 
            continue; 
        }

        else if (start_bool == end_bool) {
            continue; 
        }

        else{ 
            std::cout <<"inside assignment loop" << std::endl;
            std::cout << s << std::endl;
            std::array <int,3> loc; 
            std::string strr; // temp string 
            int k;    // temp container 
            int y{0};     // index 
            while (!ss.eof()){
                ss >> strr;  
                
                if (std::stringstream(strr) >> k){
                    loc[y] = k; 
                    ++y;
                }
            }

            locations.push_back(loc); 

        }

    }

    // 
    std::cout << "Locations are: " << std::endl;
    for (auto& lll: locations){
        for (auto& p: lll){
            std::cout << p << " | ";
        }
        std::cout << std::endl;
    }


    return 0; 
}