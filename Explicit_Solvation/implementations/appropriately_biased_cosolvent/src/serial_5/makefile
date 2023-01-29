CC = g++
OMPCC = g++-12

# use this set of flags if on Linux 
LINUXFLAGS = -std=c++14 -pedantic -O3 -march=native -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

# use this set of flags if on MacOS M1
MACFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10 # -fopenmp
OMPFLAGS = -fopenmp -O3 -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

McLattE: main.cpp misc.cpp classes.cpp
	@echo "Compiling main for simulation..."
	@$(CC) $(LINUXFLAGS) main.cpp misc.cpp classes.cpp -o McLattE
	@echo "Compiled!" 

dMcLattE: debug.cpp misc.cpp classes.cpp
	@echo "Compiling the debugging script..."
	@$(CC) $(LINUXFLAGS) debug.cpp misc.cpp classes.cpp -o dMcLattE
	@echo "Compiled!"

mac_MLE: main.cpp misc.cpp classes.cpp
	@echo "Compiling main for simulation..."
	@$(CC) $(MACFLAGS) main.cpp misc.cpp classes.cpp -o McLattE
	@echo "Compiled!" 

mac_dMLE: debug.cpp misc.cpp classes.cpp
	@echo "Compiling the debugging script..."
	@$(CC) $(MACFLAGS) debug.cpp misc.cpp classes.cpp -o dMcLattE
	@echo "Compiled!"
