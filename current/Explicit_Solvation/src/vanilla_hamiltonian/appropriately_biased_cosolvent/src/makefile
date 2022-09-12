CC = g++

# use this set of flags if on Linux 
CFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -march=native

# use this set of flags if on MacOS M1
# CFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10

all: clean main unbiased
	@echo "All simulation objects should be compiled."

main: main.cpp misc.cpp classes.cpp
	@echo "Compiling biased simulation..."
	@$(CC) $(CFLAGS) main.cpp misc.cpp classes.cpp -o cosolv_ani
	@echo "Compiled!" 

unbiased: unbiased.cpp misc.cpp classes.cpp 
	@echo "Compiling unbiased simulation..."
	@$(CC) $(CFLAGS) unbiased.cpp misc.cpp classes.cpp -o ub_cosolv_ani 
	@echo "Compiled!"

clean: 
	@rm -f cosolv_ani || true 
	@echo "Removed executable for biased simulation."
	@rm -f ub_cosolv_ani || true 
	@echo "Removed executable for unbiased simulation."
