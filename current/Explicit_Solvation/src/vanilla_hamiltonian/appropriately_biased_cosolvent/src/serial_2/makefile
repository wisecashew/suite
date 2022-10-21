CC = g++
OMPCC = g++-12

# use this set of flags if on Linux 
LINUXFLAGS = -std=c++14 -pedantic -O3 -march=native -Wall -Wextra -Wdisabled-optimization -Werror -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

# use this set of flags if on MacOS M1
MACFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10 -fopenmp
OMPFLAGS = -fopenmp -O3 -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

all: clean linux unbiased mac
	@echo "All simulation objects should be compiled."

isotropic: main.cpp misc.cpp classes.cpp 
	@echo "Compiling an isotropic simulation..."
	@$(CC) $(LINUXFLAGS) main.cpp misc.cpp classes.cpp -o isotropic
	@echo "Compiled!"

clinux: main.cpp misc.cpp classes.cpp
	@echo "Compiling main for cosolvent_regime simulation..."
	@$(CC) $(LINUXFLAGS) main.cpp misc.cpp classes.cpp -o cosolv_ani
	@echo "Compiled!" 

flinux: main.cpp fh_misc.cpp classes.cpp
	@echo "Compiling main for flory-huggins simulation..."
	@$(CC) $(LINUXFLAGS) main.cpp fh_misc.cpp classes.cpp -o flory-huggins
	@echo "Compiled!" 

ldebug: debug.cpp misc.cpp classes.cpp 
	@echo "Compiling debugging simulation on linux..."
	@$(CC) $(LINUXFLAGS) debug.cpp misc.cpp classes.cpp -o debug
	@echo "Compiled!" 

home: main.cpp misc.cpp classes.cpp 
	@echo "Compiling unbiased simulation..."
	@$(CC) $(CFLAGS) main.cpp misc.cpp classes.cpp -o cosolv_ani
	@echo "Compiled!"

mac: main.cpp misc.cpp classes.cpp
	@echo "Compiling simulation on a mac..."
	@$(OMPCC) $(OMPFLAGS) main.cpp misc.cpp classes.cpp -o cosolv_ani
	@echo "Compiled!"

debug: debug.cpp misc.cpp classes.cpp
	@echo "Compiling simulation in debug mode..."
	@$(OMPCC) $(OMPFLAGS) debug.cpp misc.cpp classes.cpp -o debug
	@echo "Compiled!"

clean: 
	@rm -f cosolv_ani || true 
	@echo "Removed executable for biased simulation."
	@rm -f ub_cosolv_ani || true 
	@echo "Removed executable for unbiased simulation."
