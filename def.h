#include <cmath>

#define pi M_PI
#define ROW 5
#define COL 5
#define iter_max 1000

#define method 0 // case 0: CG method ; case 1: SOR method

// 2D Poisson Equation: Constants
const double L = 1.0;   // 1-D computational domain size
const int N = 16;       // number of computing cells
const double u0 = 1.0;  // background density
const double amp = 0.5; // sinusoidal amplitude
const double cfl = 1.0; // Courant condition factor