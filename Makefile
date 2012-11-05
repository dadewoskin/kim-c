CC = gcc
LIBS = -lm
#PLL = 
PLL = _omp
FLAGS = -O2 -fopenmp -std=c99 -pedantic
DEPS = rhs.h
OBJ = rhs$(PLL).o kimFE$(PLL).o

kimFE$(PLL): $(OBJ)
	gcc -o kimFE$(PLL) $(OBJ) $(LIBS) $(FLAGS)

#rhs$(PLL).o: $(DEPS)
#kimFE$(PLL).o: $(DEPS)
$(OBJ).o: $(DEPS)

clean:
	rm kimFE$(PLL) $(OBJ)
