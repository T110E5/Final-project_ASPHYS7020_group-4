#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include "def.h"


void CG_method(double *x0, double tol){ // tol: tolerance
	int i, j, k;

	double **A;
	double *x, *b, *d, *tmp, *res, *res_new;
	double alpha, beta;

    A = (double **) malloc(ROW*COL * sizeof(double));
    for (i=0;i<ROW*COL;i++)
    {
      A[i] = (double *) malloc(ROW*COL * sizeof(double));
    }
    x      = (double *) malloc(ROW*COL * sizeof(double));
    b      = (double *) malloc(ROW*COL * sizeof(double));
    d      = (double *) malloc(ROW*COL * sizeof(double));
	tmp  = (double *) malloc(ROW*COL * sizeof(double));
    res    = (double *) malloc(ROW*COL * sizeof(double));
    res_new= (double *) malloc(ROW*COL * sizeof(double));
	// Initialize	when k=0 => xk = 0, rk = b, dk = rk; definition of matrix A, alpha and beta

	make_A(A);
	dot_MV(A, x, tmp); // A*d_k-1

	for (int i=0; i<ROW; i++){
		for (int j=0; j<COL; j++){
			res[COL*i+j] = b[COL*i+j] - tmp[COL*i+j];
			d[COL*i+j] = res[COL*i+j];
		}
	}

	// Next x, r
	for (k=0; k<iter_max ; k++){

		dot_MV(A, d, tmp);
		alpha = dot_VV(res, res)/ dot_VV(d ,tmp);
		for (int i=0; i<ROW; i++){
			for (int j=0; j<COL; j++){
				x[COL*i+j] = x[COL*i+j] + alpha*d[COL*i+j];
				res_new[COL*i+j] = res[COL*i+j] - alpha*tmp[COL*i+j];	
			}
		}

	// Convergence check |Ax -b|^2 =< tol
		if (res_square(res_new) < tol){

			free(A);
			free(x);
			free(b);
			free(d);
			free(res);
			free(res_new);

			break;
		}

	// Next d
		beta = dot_VV(res_new, res_new)/dot_VV(res, res);
		for (int i=0; i<ROW; i++){
			for (int j=0; j<COL; j++){
				d[COL*i+j] = res_new[COL*i+j] - beta*d[COL*i+j];
				res[COL*i+j] = res_new[COL*i+j];
			}
		}
	// Iterate

	}
	printf("Iteration: %f", k);
}

//------------------
//	Functions
//------------------
double res_square(double *res_new){
	double sum = 0;

	for (int i=0; i<ROW*COL; i++){
		sum += res_new[i];
	}

	return sum;
}

double dot_MV(double **A, double *x, double *b){

    for (int i=0; i<ROW*COL; i++){
            b[i] = 0;
    }

    for (int i=0; i<ROW*COL; i++){
        for (int j=0; j<ROW*COL; j++){
            b[i] += A[i][j] * x[j];
        }
    }
}

double dot_VV(double *A, double *B) {
	double sum = 0.0;
	for (int i = 0; i < ROW*COL; i++) {
		sum += A[i] * B[i];
	}
	return sum;
}

double make_A(double **A){ // generate the matrix A

}
/*
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
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			temp[i] += M[i * N + j] * V[j];
		}
	}
	return temp;
	free(temp);
}



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