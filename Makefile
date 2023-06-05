all: CG.cpp main.cpp
	g++ CG.cpp main.cpp -fopenmp -lgomp

clean:
	rm -f *.out