CC = g++
CFLAGS = -std=c++14 -Wall -Wextra -Wdisabled-optimization -Werror -pedantic -Ofast -march=native

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
	@rm biased || true > /dev/null
	@echo "Removed executable for biased simulation."
	@rm unbiased || true 
	@echo "Removed executable for unbiased simulation."
