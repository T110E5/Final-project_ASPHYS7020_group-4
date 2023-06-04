#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "def.h"
#include "CG_method.h"
#include "SOR_method.h"

int main (){
    char *fileroot;
//-----------------------------
//  2D Poisson Equation
//-----------------------------
    double **u, **u_ref, **rho;

    // Derived constants
    double dx = L / (N - 1); // spatial resolution
    double dy = L / (N - 1); // spatial resolution

    u     = (double **) malloc(ROW*COL * sizeof(double));
    u_ref = (double **) malloc(ROW*COL * sizeof(double));
    for (int i=0;i<ROW*COL;i++){
        u[i]      = (double *) malloc(ROW*COL * sizeof(double));
        u_ref[i]  = (double *) malloc(ROW*COL * sizeof(double));
    }    

    // Initialize rho and u values
    init_rho(rho, dx, dy);
    init_u(u);

    // Compute u_ref values
    ref_func(u, dx, dy);

//-----------------------------
//  Poisson Solver
//-----------------------------
    switch (method){
        case 0:
            /* CG method code */
            //CG_method();
            // save array u
            fileroot = "./CG_result.txt";
            save_M(u, fileroot);
            break;

        case 1:
            /* SOR method code */
            //SOR_method();
            // save array u
            fileroot = "./SOR_result.txt";
            save_M(u, fileroot);            
            break;
    }

    free(u);
    free(u_ref);

    return 0;    
}

//---------------------
//  Setting
//---------------------

// Define a reference analytical solution
double ref_func(double **M, double dx, double dy) {
    double k = 2.0 * pi / L;                      // wavenumber
    for (int i=0; i<ROW; i++){
        for (int j=0; j<COL; j++){
            M[i][j] = (amp * amp) * sin(k * i * dx) * sin(k * j * dy);
        }
    }
}

// Define an initial density distribution
double init_rho(double **rho, double dx, double dy) {
    double k = 2.0 * M_PI / L;
    for (int i=0; i<ROW; i++){
        for (int j=0; j<COL; j++){
            rho[i][j] = -2.0 * (amp * amp) * (k * k) * sin(k * i * dx) * sin(k * j * dy);
        }
    }
}

// Initialize u values to 0
double init_u(double **u) {
    for (int i=0; i<ROW; i++){
        for (int j=0; j<COL; j++){
            u[i][j] = 0;
        }
    }
}

// Save the array to the text file
void save_M(double **M, const std::string& filename){
    FILE* file = fopen(filename.c_str(), "w");

    if (file != nullptr){
        for (int i = 0; i < ROW; i++){
            for (int j = 0; j < COL; j++){
                fprintf(file, "%f", M[i][j]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
        printf("Array data saved to file \n");
    }
    else{
        printf("Unable to open file");
    }
}