// Initial Setting

#define N  128  // number of computing cells
#define ROW N   // it must be equal to N (grids)
#define COL N   // it must be equal to N (grids)
#include <cmath>
#define pi  M_PI
#define iter_max 1000
#define NUM_THREADS 12
const double tol    = 1e-6; // tol: tolerance

//------------------------
// case 0: CG method
// case 1: SOR method
//------------------------
#define method 0

// 2D Poisson Equation: Constants
const double L      = 1.0;  // 1-D computational domain size
const double u0     = 0.0;  // background density
const double amp    = 1.0;  // sinusoidal amplitude
const double cfl    = 1.0;  // Courant condition factor