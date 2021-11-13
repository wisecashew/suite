#include <iostream> 
#include <vector> 
#include <string>
#include <algorithm> 
#include <map>
#include <random>
#include <chrono>
#include <getopt.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  int opt;

  int a; 

  while ((opt = getopt(argc, argv, ":a:b:X")) != -1) 
  {
     switch (opt) 
     {
      case 'a':
      	std::cout << "Option a has arg: "<<optarg << std::endl;
      	a = atoi(optarg);
        break;
      case 'b':
      	std::cout << "Option b has arg: "<<optarg << std::endl;
        break;
      case 'X':
      	std::cout << "Option X was provided." << std::endl;
        break;
      case '?':
        std::cout << "Unknown option " << optarg << " was provided." << std::endl;
        break;
      case ':':
      	std::cout << "Missing arg for " << static_cast<char>(optopt) << std::endl;
        break;
     }
  }

    /* Get all of the non-option arguments */
  std::cout << "a has value " << a << std::endl;
  if (optind < argc) 
  {
    printf("Non-option args: ");
    while (optind < argc)
      printf("%s ", argv[optind++]);
    printf("\n");
  }
  
  return 0;
}