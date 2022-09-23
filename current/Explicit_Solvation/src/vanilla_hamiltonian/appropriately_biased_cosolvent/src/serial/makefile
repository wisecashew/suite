CC = g++
OMPCC = g++-12
# use this set of flags if on Linux 
LINUXFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -march=native -fopenmp

# use this set of flags if on MacOS M1
MACFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10 -fopenmp
OMPFLAGS = -fopenmp -O3 -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Wnonnull -Wnull-dereference -Wno-return-local-addr -Wreturn-type

all: clean linux unbiased mac
	@echo "All simulation objects should be compiled."

linux: main.cpp misc.cpp classes.cpp
	@echo "Compiling main simulation..."
	@$(CC) $(LINUXFLAGS) main.cpp misc.cpp classes.cpp -o cosolv_ani
	@echo "Compiled!" 

home: main.cpp misc.cpp classes.cpp 
	@echo "Compiling unbiased simulation..."
	@$(CC) $(CFLAGS) main.cpp misc.cpp classes.cpp -o cosolv_ani
	@echo "Compiled!"

mac: main.cpp misc.cpp classes.cpp
	@echo "Compiling parallel simulation..."
	@$(OMPCC) $(OMPFLAGS) main.cpp misc.cpp classes.cpp -o cosolv_ani
	@echo "Compiled!"

debug: debug2.cpp misc.cpp classes.cpp
	@echo "Compiling parallel simulation..."
	@$(OMPCC) $(OMPFLAGS) debug2.cpp misc.cpp classes.cpp -o debug
	@echo "Compiled!"

clean: 
	@rm -f cosolv_ani || true 
	@echo "Removed executable for biased simulation."
	@rm -f ub_cosolv_ani || true 
	@echo "Removed executable for unbiased simulation."
