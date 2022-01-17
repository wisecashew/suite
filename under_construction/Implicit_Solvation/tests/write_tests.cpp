#include <iostream> 
#include <chrono> 
#include <string> 
#include <fstream> 
#include <sstream> 

int main() {
	std::ofstream dump_file("oss.txt");
    std::string s{ "this is a silly string."}; 
    std::ostringstream os; 

    auto start = std::chrono::high_resolution_clock::now(); 

    std::ofstream dump_file_1("oss.txt", std::ios::app); 

    for (int i{0}; i< 1000000; ++i){
        os << s <<"\n";
    }
    const auto& str = os.str() ; 
    dump_file_1.write(str.c_str(), static_cast<std::streamsize> (str.size() )); 

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

    printf("Time required by oss is %lld milliseconds.", duration.count() ); 

    start = std::chrono::high_resolution_clock::now(); 

    std::ofstream dump_file_2("throw.txt", std::ios::app); 

    for (int i{0}; i< 1000000; ++i){
        dump_file_2 << s << "\n";
    }
    

    stop = std::chrono::high_resolution_clock::now(); 
    duration = std::chrono::duration_cast<std::chrono::milliseconds> (stop-start); 

    printf("Time required by throw operator is %lld milliseconds.", duration.count() ); 

    return 0; 
}
