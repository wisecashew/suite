CC = g++

# use this set of flags if on Linux 
# CFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -march=native

# use this set of flags if on MacOS M1
CFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -mcpu=apple-a10

all: clean biased unbiased
	@echo "All simulation objects should be compiled."

biased: biased.cpp misc.cpp classes.cpp
	@echo "Compiling biased simulation..."
	@$(CC) $(CFLAGS) biased.cpp misc.cpp classes.cpp -o biased
	@echo "Compiled!" 

unbiased: unbiased.cpp misc.cpp classes.cpp 
	@echo "Compiling unbiased simulation..."
	@$(CC) $(CFLAGS) unbiased.cpp misc.cpp classes.cpp -o unbiased 
	@echo "Compiled!"

clean: 
	@rm -f biased || true 
	@echo "Removed executable for biased simulation."
	@rm -f unbiased || true 
	@echo "Removed executable for unbiased simulation."
