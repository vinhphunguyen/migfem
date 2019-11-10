	
int      FindSpan(int n, int p, double u, double U[]);
void     BasisFuns( int i, double u, int p, double U[], double* N);
void     dersBasisFuns(int i, double u, int p, int order, double knot[], double **ders);
double   OneBasisFun(int p, int m, double U[], int i, double u);
void     dersOneBasisFuns(int p, int m, double U[], int i, double u, int n, double* ders);
double** init2DArray(int x, int y);
void     free2Darray(double **array, int x);

