// Function declarations

double res_square(double *res_new);
double dot_MV(double **A, double *x, double *b);
double dot_VV(double *A, double *B);
double make_MA(double **A);
double make_Vx(double *x);
double make_Vb(double *b);

void CG_method(double *x0, double tol);