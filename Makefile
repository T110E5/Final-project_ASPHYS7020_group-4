all: CG_method.cpp SOR_method.cpp main.cpp
	g++ CG_method.cpp SOR_method.cpp main.cpp -o main.out

clean:
	rm -f *.out