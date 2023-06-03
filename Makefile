all: CG_method.cpp main.cpp
	g++ CG_method.cpp main.cpp -o main.out -lgomp

clean:
	rm -f *.out