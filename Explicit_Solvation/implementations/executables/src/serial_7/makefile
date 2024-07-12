CC = g++
OMPCC = g++-12
ICC = icpx

# use this set of flags if on Linux 
FFLAGS = -std=c++14 -pedantic -O3 -march=native -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -g
WFLAGS = -std=c++14 -pedantic -O3 -march=native -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings
IFLAGS = -std=c++14 -O3 -flto -ipo -march=native -pedantic -Wall -Wextra -Wall -Warray-bounds -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings 

# use this set of flags if on MacOS M1
MACFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10 # -fopenmp
OMPFLAGS = -fopenmp -O3 -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

main: main.cpp misc.cpp Particle.cpp Polymer.cpp setup.cpp energetics.cpp verify.cpp dump.cpp lattice_directions.cpp extraction.cpp run.cpp debug.cpp perturb.cpp run.cpp
	@echo "Compiling main for simulation..."
	@$(CC) $(FFLAGS) main.cpp misc.cpp Particle.cpp Polymer.cpp setup.cpp extraction.cpp energetics.cpp verify.cpp dump.cpp debug.cpp lattice_directions.cpp perturb.cpp run.cpp -o flogotts
	@echo "Compiled!" 

wmain: main.cpp misc.cpp Particle.cpp Polymer.cpp setup.cpp energetics.cpp verify.cpp dump.cpp lattice_directions.cpp extraction.cpp run.cpp
	@echo "Compiling main for simulation..."
	@$(CC) $(WFLAGS) main.cpp misc.cpp Particle.cpp Polymer.cpp setup.cpp extraction.cpp energetics.cpp verify.cpp dump.cpp lattice_directions.cpp -o flogotts
	@echo "Compiled!" 
