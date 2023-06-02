all: CG_method.cpp main.cpp
	g++  CG_method.cpp main.cpp -o main.out

clean:
	rm -f *.out