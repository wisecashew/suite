#include <iostream> 
#include <sstream> 
#include <fstream>
#include <vector>
#include <string> 
#include <regex> 

int main() {


	// parse through the content file, and get total number of polymers in the box 
	std::regex start ("START"); 
	std::ifstream myfile ("testfile.txt"); 
	std::regex end ("END");
	std::string mystring; 
	std::vector <std::string> myStartContents; 
	std::vector <std::string> myEndContents; 
	 if (myfile.is_open() ) { // check if file is open 
		while (myfile.good() ) {
			std::getline(myfile, mystring); // pipe file's content into stream 
			if (std::regex_search(mystring, start)){
				std::cout << mystring << std::endl; // pipe stream's content to standard output
				myStartContents.push_back(mystring); 
			}
			else if (std::regex_search(mystring, end)){
				std::cout << mystring << std::endl; 
				myEndContents.push_back( mystring );
			}
		}
	}
	std::cout << "Number of polymers thrown in (start metric): " << myStartContents.size() << std::endl;
	std::cout << "Number of polymers thrown in (end metric): " << myEndContents.size() << std::endl;
	
	/*
	std::string teststr = "this is a test str"; 

	std::regex reg ("test"); 
	bool tval = std::regex_search(teststr, reg); 
	std::cout << "truth value is " << tval << std::endl;

	std::cout << "========" << std::endl; 
	std::ifstream myfile ("testfile.txt");  // declare the file you want to read from 
	std::string mystring; 
	std::vector <std::string> mycontents;  

	if (myfile.is_open() ) { // check if file is open 
		while (myfile.good() ) {
			std::getline(myfile, mystring); // pipe file's content into stream 
			std::cout << mystring << std::endl; // pipe stream's content to standard output
			mycontents.push_back(mystring); 
		}
	}
	std::cout << "================" << std::endl; 
	std::cout << "let's check if mycontents has everything..." << std::endl; 

	for (auto s: mycontents){
		std::cout << s << std::endl;
	}


	
	std::cout << "================" << std::endl;
	std::vector <std::vector <int>> loc_list; 
	
	// create a list of all coordinates 
	
	std::regex start("START"); 
	std::regex end ("END"); 
	for (auto s: mycontents){
		std::vector <int> loc; 
		std::stringstream ss(s);

		if (std::regex_search(s, start)){
			continue; 
		}
		else if (std::regex_search(s, end)){
			continue; 
		}

		bool b = false; 
		for (int i =0; ss >> i;){
			b = true; 
			loc.push_back(i);
			std::cout << i << ", ";
		}

	std::cout << std::endl;
	if (b){
		loc_list.push_back(loc); 
		}
	
	}
	for (auto l: loc_list){
		for (auto i: l){
			std::cout << i << " "; 
		}
		std::cout << std::endl;
	}
	*/ 
	return 0;
}
