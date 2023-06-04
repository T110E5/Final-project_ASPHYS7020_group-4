#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "def.h"
#include "CG_method.h"
#include "SOR_method.h"

int main (){
    // 2D Poisson Equation


    // Poisson Solver
    switch (method){
        case 0:
            /* CG method code */
            //CG_method();
            break;

        case 1:
            /* SOR method code */
            //SOR_method();
            break;
    }

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