#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <iostream>
#define pi M_PI
#define NUM_THREADS 8


// initial condition
double error=0.0;
int i,j;
const double L=1.0;
const double cfl=1.0; 
int N=1000;
int Row=N;
int Col=N;
void ref_func(double **M, double dx, double dy){
    for (int i=0; i<Row; i++){
        for(int j=0; j<Col; j++){
            M[i][j]=sin(pi*i*dx)*cos(pi*j*dy);
        }
    }
}

void init_rho(double **rho, double dx, double dy){
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for collapse(2) shared(rho) private(i,j)
    for(int i=0; i<Row; i++){
        for(int j=0; j<Col; j++){
            rho[i][j]=-( 2.0*pi*pi*sin(pi*i*dx)*cos(pi*j*dy));
        }
    }
}
//Initialize the potential to one
void init_u(double **u){
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for collapse(2) shared(u) private(i,j)
    for (int i=0; i<Row; i++){
        for(int j=0; j<Col; j++){
            u[i][j]=0;
        }
    }
}
int main (void){
    omp_set_num_threads(NUM_THREADS);
    int i, j;
    double dx=L/(N-1);
    double dy=L/(N-1);
    int iter;
    double tol=1.e-7;
    
    double **u     = (double **) malloc(Row * sizeof(double));
    double **u_in  = (double **) malloc(Row * sizeof(double));
    double **u_ref = (double **) malloc(Row * sizeof(double));
    double **rho   = (double **) malloc(Row * sizeof(double));
    for (int i=0;i<Row;i++){
        u[i]      = (double *) malloc(Col * sizeof(double));
        u_in[i]   = (double *) malloc(Col * sizeof(double));
        u_ref[i]  = (double *) malloc(Col * sizeof(double));
        rho[i]    = (double *) malloc(Col * sizeof(double));        
    }
    init_u(u);
    init_rho(rho,dx,dy);
    //Compute reference function
    ref_func(u_ref, dx, dy);
//SOR method______________________________________________________________
    double w = 1.8;

    double start, end;
    double time_used;
    //backup the data
    u_in=u;
    start = omp_get_wtime();
    for (int iter=1;iter<10000000;iter++){
    	double err = 0.0;
        double res;
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for shared(res,u) private(i,j) reduction(+:err)
        for (i=1; i<N-1; i++){
        	for (j=1;j<N-1;j++){
        		int node = (i+j)%2;
        		if (node ==0){
        			res = (u[i+1][j]+u[i-1][j]+u[i][j-1]+u[i][j+1]-4 * u_in[i][j]-dx * dx * rho[i][j])/(dx * dx);
            			u[i][j]=u_in[i][j]+w*dx*dx*res/4;
            			err += dx * dx * abs(res/u[i][j])/(N * N);
            		}
        		else if (node==1){
        			res = (u[i+1][j]+u[i-1][j]+u[i][j-1]+u[i][j+1]-4 * u_in[i][j]-dx * dx * rho[i][j])/(dx * dx);
            			u[i][j]=u_in[i][j]+w * dx * dx * res/4;
            			err += dx * dx * abs(res/u[i][j])/(N * N);
        		}
        	}
        }
    // Boundary conditions
    	for (int j=0; j<Col; j++){
    		u[0][j] = 0;
        	u[Row-1][j] = 0;
    	}

    	for (int i=0; i<Row; i++) {
        	u[i][0] = u[i][1];
        	u[i][Col-1] = u[i][Col-2];
    	}
    	//printf("err %f \n",err);
    	
    	iter+=1;
    	end = omp_get_wtime();
    	if (err<tol){
        	time_used = end-start;
        	#pragma omp parallel for collapse(2) shared(u,u_ref) private(i,j) reduction(+:error)
        	for (i=0; i<Row; i++){
        		for (j=0; j<Col; j++){
        			#pragma omp atomic
        			error += abs(u[i][j]-u_ref[i][j])/(N * N); 
                	}
            	}
            	printf("error is %f \n Spending %f\n iterations %d \n", error, time_used, iter);
            	break;
        }
    }
//End SOR________________________________________________________________
    // Free memory
    free(u);
    free(rho);
    free(u_ref);
    return 0;
}
		
