#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include "def.h"
#include "CG.h"
#include <omp.h>

void CG_method(double **u, double **rho, double dx, double dy){
	double **A;
	double *x, *b, *d, *tmp, *res, *res_new;
	double alpha, beta;

	int iter;

	// double start=omp_get_wtime()

	A = (double **) malloc(ROW*COL * sizeof(double));
    	for (int i=0;i<ROW*COL;i++){
      		A[i] = (double *) malloc(ROW*COL * sizeof(double));
    	}
    	x	= (double *) malloc(ROW*COL * sizeof(double));
    	b      	= (double *) malloc(ROW*COL * sizeof(double));
    	d      	= (double *) malloc(ROW*COL * sizeof(double));
	tmp	= (double *) malloc(ROW*COL * sizeof(double));
    	res    	= (double *) malloc(ROW*COL * sizeof(double));
    	res_new	= (double *) malloc(ROW*COL * sizeof(double));
	
	// Initialize	when k=0 => xk = 0, rk = b, dk = rk; definition of matrix A, alpha and beta
	make_MA(A);
	make_Vx(u, x);
	make_Vb(b, rho, dx);
	dot_MV(A, x, tmp); // A*d_k-1
	int i,j;
	int tid;
	int NUM_THREADS=6;
	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel
	{
		tid=omp_get_thread_num();
		int nt=omp_get_num_threads();
	#pragma omp for
		
		for (int i=0; i<ROW; i++){
			for (int j=0; j<COL; j++){
				res[COL*i+j] = b[COL*i+j] - tmp[COL*i+j];
				d[COL*i+j] = res[COL*i+j];
				printf("%d/%d\n",tid,nt);
			}
		}
	}

	// Loop
	// Next x, r
	int k;
	for (int k=0; k<iter_max ; k++){
		dot_MV(A, d, tmp);
		alpha = dot_VV(res, res)/ dot_VV(d ,tmp);
		#pragma omp parallel for shared(x,res_new,res,d,tmp) privated(i,j)
		for (int i=0; i<ROW; i++){
			for (int j=0; j<COL; j++){
				x[COL*i+j] = x[COL*i+j] + alpha*d[COL*i+j];
				res_new[COL*i+j] = res[COL*i+j] - alpha*tmp[COL*i+j];
			}
		}

	// Convergence check |Ax -b|^2 =< tol
		if (res_square(res_new) < tol){

			// New array u
			int i,j;
			#pragma omp parallel for shared(u,x) private(i,j)
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
		#pragma omp parallel for shared(d,res,res_new) private(i,j)
		for (int i=0; i<ROW; i++){
			for (int j=0; j<COL; j++){
				d[COL*i+j] = res_new[COL*i+j] - beta*d[COL*i+j];
				res[COL*i+j] = res_new[COL*i+j];
			}
		}
	// Iterate

	}
	// double end =omp_get_wtime();
	//printf("Work time %f seconds",end-start);

	printf("Iteration: %d\n", iter);
}

//-----------------------
//	Calculate Functions
//-----------------------
double res_square(double *res_new){
	double sum = 0;
	int i;
	#pragma omp parallel for shared(sum, res_new) private(i) reduction(+:sum)
	for (int i=0; i<ROW*COL; i++){
		sum += pow(res_new[i], 2);
	}

	return sqrt(sum);
}

void dot_MV(double **A, double *x, double *b){

    for (int i=0; i<ROW*COL; i++){
            b[i] = 0;
    }
	//first define a zero matrix to add up
    int i,j;
    #pragma omp parallel for shared(b,A,x) private(i,j) reduction(+:b)
    for (int i=0; i<ROW*COL; i++){
        for (int j=0; j<ROW*COL; j++){
            b[i] += A[i][j] * x[j];
        }
    }
}

double dot_VV(double *A, double *B) {
	double sum = 0.0;
	int i;
	#pragma omp parallel for shared(A,B,sum) private(i) reduction(+:sum)
	for (int i = 0; i < ROW*COL; i++) {
		sum += A[i] * B[i];
	}
	return sum;
}

//----------------------------------------
//	Make matrix and vectors of CG method
//----------------------------------------
void make_MA(double **A){ // generate the matrix A
    	int i,j,k,l;
	for (k=0;k<ROW;k++){
        for (l=0;l<COL;l++){
            if (k==l){
                if (k==0 || k==ROW-1){
                  for (i=0;i<ROW;i++){
                      A[COL*k+i][ROW*l+i]   = 1;
                  }
                }
                else{
                  for (i=0;i<ROW;i++){
                      if (i == 0){
                          A[COL*k+i][ROW*l+i]   = -1;
                          A[COL*k+i+1][ROW*l+i] = 1;
                      }
                      else if (i == ROW-1){
                          A[COL*k+i][ROW*l+i]   = -1;
                          A[COL*k+i-1][ROW*l+i] = 1;
                      }
                      else {
                          A[COL*k+i][ROW*l+i]   = -4;
                          A[COL*k+i-1][ROW*l+i] = 1;
                          A[COL*k+i+1][ROW*l+i] = 1;
                      }
                  }
                }
            }
            else if ( abs(k-l) == 1 && k!=0 && k!=ROW-1){
                for (i=0;i<ROW;i++){
                  if (i==0 || i==ROW-1)
                    A[COL*k+i][ROW*l+i] = 0;
                  else
                    A[COL*k+i][ROW*l+i] = 1;
                }
            }
            else{
                for (i=0;i<ROW;i++){
                    for (j=0;j<COL;j++){
                        A[COL*k+i][ROW*l+j] = 0;
                    }
                }
            }
		}
	}
}
void make_Vx(double **u, double *x){ // generate the vector x
	int i,j;
	#pragma omp parallel for shared(x,u) private(i,j)
	for (int i=0; i<ROW; i++){
        for (int j=0; j<COL; j++){
        	x[ROW*i+j] = u[i][j];
        }
    }
}

void make_Vb(double *b, double **rho, double dx){ // generate the vector b
	#pragma omp parallel for collapse(2) shared(b,rho) private(i,j)
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