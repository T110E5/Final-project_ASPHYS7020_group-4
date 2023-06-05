// Function declarations

double res_square(double *res_new);
void dot_MV(double **A, double *x, double *b);
double dot_VV(double *A, double *B);
void make_MA(double **A);
void make_Vx(double **u, double *x);
void make_Vb(double *b, double **rho, double dx);

void CG_method(double **u, double **rho, double dx, double dy);