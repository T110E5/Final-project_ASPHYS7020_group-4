#include <cstdlib>
#include <math.h>
#include <time.h>
#include <omp.h>
double dotproduct(double* A, double * B, int N) {
	double sum = 0.0;
	#pragma omp parallel for reduction(+ : sum)
		for (int i = 0; i < N; i++)
		{
			sum += A[i] * B[i];
			#pragma omp barrier
		}
		return sum;
}


double* vectorAdd(double* A, double* B, int N) {
	double* temp = (double*)calloc(N, sizeof(double));
	for (int i = 0; i < N; i++) {
		temp[i] = A[i] + B[i];
	}
	return temp;
	free(temp);
}

double* vectorSubtract(double* A, double* B, int N) {
	double* temp = (double*)calloc(N, sizeof(double));
	for (int i = 0; i < N; i++) {
		temp[i] = A[i] - B[i];
	}
	return temp;
	free(temp);
}

double* vectorZoom(double* A, double scaler, int N) {
	for (int i = 0; i < N; i++) {
		A[i] = scaler * A[i];
	}
	return A;
}

double* matrix_Times_vector(double* M, double* V, int N) {
	double* temp = (double*)calloc(N, sizeof(double));
	# pragma omp parllel for reduction (+:sum)
	# pragma omp for collapse(2)
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++) {
			temp[i] += M[i * N + j] * V[j];
			# pragma omp barrier
		}
	return temp;
	free(temp);
}
/*
int main(int argc, char* argv[]){
	double a[6] = { 0 };
	double b[6] = { 0 };


	for (int i = 0; i < 6; i++) {
		a[i] = double(i);
		b[i] = double(i);
	}

	
	for (int i = 0; i < 6; i++) {
		printf("%f\n",c[i]);
	}
	//printf("hello world");
	
	return EXIT_SUCCESS;
}
*/