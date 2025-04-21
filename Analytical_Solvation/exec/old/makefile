CC = g++
OMPCC = g++-12

# use this set of flags if on Linux 
FFLAGS = -std=c++14 -pedantic -Ofast -march=native -Wall -Wextra \
-Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference \
-Wno-return-local-addr -Wreturn-type -Wuninitialized -Wmaybe-uninitialized \
-Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type \
-Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function \
-Wunused-but-set-variable -Wunused-variable -Wwrite-strings
PFLAGS = -std=c++14 -pedantic -Ofast -march=native -Wall -Wextra \
-Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference \
-Wno-return-local-addr -Wreturn-type -Wuninitialized -Wmaybe-uninitialized \
-Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type \
-Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function \
-Wunused-but-set-variable -Wunused-variable -Wwrite-strings -fopenmp

# use this set of flags if on MacOS M1
MACFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10 # -fopenmp
OMPFLAGS = -fopenmp -O3 -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

main: main.cpp
	@echo -n "Compiling main for simulation... "
	@$(CC) $(FFLAGS) main.cpp -o sct
	@echo "Compiled!"

parallel: main.cpp
	@echo -n "Compiling main for simulation... "
	@$(CC) $(PFLAGS) main.cpp -o p_sct
	@echo "Compiled!"
