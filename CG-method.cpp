#include <cstdio>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>


double const L = 1.0,
			 k = 2*M_PI/L;

double dotproduct(double* A, double * B, int N) {
	double sum = 0.0;
	for (int i = 0; i < N; i++) {
		sum += (A[i] * B[i]);
	}
	return sum;
}

double* vectorAdd(double* des, double* A, double* B, int N) {
	for (int i = 0; i < N; i++) {
		des[i] = A[i] + B[i];
	}
	return des;
}

double* vectorSubtract(double* des, double* A, double* B, int N) {
	for (int i = 0; i < N; i++) {
		des[i] = A[i] - B[i];
	}
	return des;
}

double* vectorZoom(double* A, double scaler, int N) {
	for (int i = 0; i < N; i++) {
		A[i] = scaler * A[i];
	}
	return A;
}

double* matrix_Times_vector(double* des, double* M, double* V, int N) {
	for (int i = 0; i < N; i++) {
		des[i] = 0;
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			des[i] += M[i * N + j] * V[j];
		}
	}
	return des;
}

double* ref_func(double* x, double* y, int N) {
	int grid = N * N;
	double* ref = (double*)calloc(grid, sizeof(double));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			ref[i * N + j] = sin(k * x[i]) * sin(k * y[j]);
		}
	}
	return ref;
	free(ref);
}

double* density(double* x, double* y,int N) {
	int grid = N * N;
	double* d = (double*)calloc(grid, sizeof(double));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			d[i * N + j] = -2 * k * k * sin(k * x[i]) * sin(k * y[j]);
		}
	}
	return d;
	free(d);
}

double* linspace(double start, double stop, int step) {
	double* list = (double*)calloc(step, sizeof(double));
	double dx = (stop - start) / (step-1);
	for (int i = 0; i < step; i++) {
		list[i] = start + i * dx;
	}
	return list;
	free(list);
}

int main(int argc, char* argv[]){
	int side = 4;
	int N = side * side;
	double h = L / (side - 1);
	double* A = (double*)calloc(N*N, sizeof(double));
	for (int i = 0; i < N; i++) {
		if ((i / side != 0) && (i / side != side - 1) && (i % side != 0) && (i % side != side - 1)) {
			int num = i * N + i;
			A[num] = -4/(h*h);
			A[num - 1] = A[num + 1] = 1/(h*h);
			A[num - side] = A[num + side] = 1/(h*h);
		}
		else {
			int num = i * N + i;
			A[num] = 1;
		}
		
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf(" %.2f ", A[i * N + j]);
		}
		printf(" \n");
	}

	double	* x = linspace(0, L, side),
			* y = linspace(0, L, side);
	
	for (int i = 0; i < side; i++) {
		printf("%.2f", x[i]);
	}
	printf("\n");

	double	* rho = density(x, y, side),
			* ref_phi = ref_func(x,y,side);
	/*
	for (int i = 0; i < side; i++) {
		printf("%.5f ", ref_phi[i]);
	}
	*/
	printf(" \n");
	double alpha, beta, d1, d2;
	double	* phi		= (double*)calloc(N, sizeof(double)),
			* r			= (double*)calloc(N, sizeof(double)),
			* r_temp	= (double*)calloc(N, sizeof(double)),
			* p			= (double*)calloc(N, sizeof(double)),
			* Ap		= (double*)calloc(N, sizeof(double)),
			* tempv		= (double*)calloc(N, sizeof(double));
	int iteration = 0;
	double tolerance = 1.e-10;
	/*
	for (int i = 0; i < N; i++) {
		phi[i] = 1;
	}
	*/
	/*
	for (int i = 0; i < side; i++) {
		for (int j = 0; j < side; j++) {
			if ((i == 0) || (i == side - 1) || (j == 0) || (j == side - 1)) {
				int num = i * side + j;
				rho[num] = 1;
			}
		}
		
	}
	*/
	for (int i = 0; i < side; i++) {
		for (int j = 0; j < side; j++) {
			printf(" %.2f ", rho[i * side + j]);
		}
		printf(" \n");
	}
	printf(" \n");
	tempv = matrix_Times_vector(tempv, A, phi, N);
	r = vectorSubtract(r,rho,tempv,N);
	p = r;
	printf("%.16f\n", dotproduct(r, r, N));
	while (dotproduct(r,r,N)>tolerance) {
		
		// A * p
		Ap = matrix_Times_vector(Ap, A, p, N);

		// r*r
		d1 = dotproduct(r, r, N);
		// p^-1 * A * p
		d2 = dotproduct(p, Ap, N);
		alpha = d1 / d2;
		// phi = phi + alpha * p 
		tempv = vectorZoom(p, alpha, N);
		phi = vectorAdd(phi, phi, tempv, N);
		// rnew = r - alpha * A * p
		tempv = vectorZoom(Ap, alpha, N);
		r_temp = vectorSubtract(r_temp, r, tempv, N);
		// beta = (rnew * rnew) / (r * r)  
		beta = dotproduct(r_temp, r_temp, N) / d1;
		// p = r + beta * p
		tempv = vectorZoom(p, beta, N);
		p = vectorAdd(p, r, tempv, N);
		//printf("beta %.2f \n", beta);
		r = r_temp;
		iteration++;
		printf("alpha %f\nbeta %f\nd1 %f\nd2 %f\n", alpha, beta, d1, d2);
		if (iteration > 1) {
			//break;
		}
		
		
	}
	
	for (int i = 0; i < side; i++) {
		for (int j = 0; j < side; j++) {
			printf(" %f ", phi[i * side + j]);
		}
		printf(" \n");
	}
	printf(" \n");
	for (int i = 0; i < side; i++) {
		for (int j = 0; j < side; j++) {
			printf(" %.2f ", ref_phi[i * side + j]);
		}
		printf(" \n");
	}
	//printf("test%d", -2 * 1);
	return EXIT_SUCCESS;
}
