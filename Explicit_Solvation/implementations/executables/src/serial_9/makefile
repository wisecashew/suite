# get the appropriate compiler
CC := g++

# use this set of flags if on Linux 
FLAGS := -std=c++14 -pedantic -Ofast -march=native -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -Wunused-parameter -Wno-error=unused-parameter # -g

# use this set of flags if on MacOS M1
MACFLAGS := -std=c++14 -pedantic -Ofast -mcpu=apple-a10 -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-stack-address -Wreturn-type -Wuninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-const-variable -Wunused-variable -Wwrite-strings

# include directories
INCLUDES := -I/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Explicit_Solvation/implementations/executables/src/serial_9/base \
            -I/scratch/gpfs/satyend/MC_POLYMER/polymer_lattice/lattice_md/Explicit_Solvation/implementations/executables/src/serial_9/potts

# source files
SRC := main.cpp $(wildcard base/*.cpp) $(wildcard potts/*.cpp) $(wildcard fhp/*.cpp) $(wildcard cg/*.cpp)

main: $(SRC)
	@echo -n "Compiling the main for simulation... "
	@$(CC) $(FLAGS) $(SRC) -o McLattE.exe
	@echo "Compiled!"

