// Function declarations

double res_square(double *res_new);
double dot_MV(double **A, double *x, double *b);
double dot_VV(double *A, double *B);
void make_MA(double **A);
double make_Vx(double **u, double *x);
double make_Vb(double *b, double **rho, double dx);

void CG_method(double **u, double **rho, double dx, double dy);