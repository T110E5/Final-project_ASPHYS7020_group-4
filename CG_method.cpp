#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include "def.h"
#include "CG_method.h"


void CG_method(double **u, double **rho, double dx, double dy){
	double **A;
	double *x, *b, *d, *tmp, *res, *res_new;
	double alpha, beta;

	int iter;

	// start-time

    A = (double **) malloc(ROW*COL * sizeof(double));
    for (int i=0;i<ROW*COL;i++){
      A[i] = (double *) malloc(ROW*COL * sizeof(double));
    }
    x		= (double *) malloc(ROW*COL * sizeof(double));
    b      	= (double *) malloc(ROW*COL * sizeof(double));
    d      	= (double *) malloc(ROW*COL * sizeof(double));
	tmp		= (double *) malloc(ROW*COL * sizeof(double));
    res    	= (double *) malloc(ROW*COL * sizeof(double));
    res_new	= (double *) malloc(ROW*COL * sizeof(double));

	// Initialize	when k=0 => xk = 0, rk = b, dk = rk; definition of matrix A, alpha and beta
	make_MA(A);
	make_Vx(u, x);
	make_Vb(b, rho, dx);
	dot_MV(A, x, tmp); // A*d_k-1

	for (int i=0; i<ROW; i++){
		for (int j=0; j<COL; j++){
			res[COL*i+j] = b[COL*i+j] - tmp[COL*i+j];
			d[COL*i+j] = res[COL*i+j];
		}
	}

	// Loop
	// Next x, r
	for (int k=0; k<iter_max ; k++){

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

			// New array u
			for (int i=0; i<ROW; i++){
				for (int j=0; j<COL; j++){
					u[i][j] = x[ROW*i+j];
				}
			}

			free(A);
			free(x);
			free(b);
			free(d);
			free(res);
			free(res_new);

			iter = k;
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
	// end-time

	printf("Iteration: %d", iter);
}

//-----------------------
//	Calculate Functions
//-----------------------
double res_square(double *res_new){
	double sum = 0;

	for (int i=0; i<ROW*COL; i++){
		sum += pow(res_new[i], 2);
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

//----------------------------------------
//	Make matrix and vectors of CG method
//----------------------------------------
void make_MA(double **A){ // generate the matrix A
	for (int i=0; i<ROW*ROW; i++){
        for (int j=0; j<ROW*ROW; j++){
			A[i][j] = 0;
            if (i==j){
                A[i][j] = 4;
            }
            else if (i==j-1 || i==j+1){
                A[i][j]=-1;
            }
            else if (i==j-ROW || i==j+ROW){
                A[i][j]=-1;
            }            
		}
	}
}

double make_Vx(double **u, double *x){ // generate the vector x
	for (int i=0; i<ROW; i++){
        for (int j=0; j<COL; j++){
        	x[ROW*i+j] = u[i][j];
        }
    }
}

double make_Vb(double *b, double **rho, double dx){ // generate the vector b
	for (int i=0; i<ROW; i++){
        for (int j=0; j<COL; j++){
			if (i==0 || i==ROW-1 || j==0 || j==COL-1){
				b[ROW*i+j] = 0; // Boundary
			}
			else{
				b[ROW*i+j] = (dx * dx) * rho[i][j];
			}
        }
    }
}