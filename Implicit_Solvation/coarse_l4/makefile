CC = g++
OMPCC = g++-12
ICC = icpx
# use this set of flags if on Linux 
LINUXFLAGS = -std=c++14 -pedantic -O3 -march=native -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type
INTELFLAFS = -O3 -std=c++14 -lto -ipa -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type


# use this set of flags if on MacOS M1
MACFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10
OMPFLAGS = -fopenmp -O3 -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

lin_CG: main.cpp misc.cpp classes.cpp 
	@echo "Compiling the coarse-grained simulation..."
	@$(CC) $(LINUXFLAGS) main.cpp misc.cpp classes.cpp -o CoarseGround_L1
	@echo "Compiled!"

lin_dCG: debug.cpp misc.cpp classes.cpp
	@echo "Compiling the debugging script..."
	@$(CC) $(LINUXFLAGS) debug.cpp misc.cpp classes.cpp -o dCoarseGround_L1
	@echo "Compiled!" 

lin_iCGd: main.cpp misc.cpp classes.cpp 
	@$(ICC) $(INTELFLAGS) main.cpp misc.cpp classes.cpp -o iCoarseGround_L1

mac_CG: main.cpp misc.cpp classes.cpp
	@echo "Compiling the coarse-grained simulation..."
	@$(CC) $(MACFLAGS) main.cpp misc.cpp classes.cpp -o CoarseGround_L1
	@echo "Compiled!" 

mac_dCG: debug.cpp misc.cpp classes.cpp 
	@echo "Compiling the debugging script..."
	@$(CC) $(MACFLAGS) debug.cpp misc.cpp classes.cpp -o dCoarseGround_L1
	@echo "Compiled!"
