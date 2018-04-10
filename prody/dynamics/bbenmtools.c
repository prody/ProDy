/*****************************************************************************/
/*                                                                           */
/*                   Tools for bbENM calculations in ProDy.                    */
/*                                                                           */
/*****************************************************************************/
/* Author: Cihan Kaya, She Zhang */
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

/* ---------- Numerical Recipes-specific definitions and macros ---------- */
#define NR_END 1
#define FREE_ARG char*
#define square(x) x * x

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/* Other structures */
typedef struct {float X[3];int model;} Atom_Line;
typedef struct {Atom_Line *atom;} PDB_File;
typedef struct {int **IDX;double *X;} dSparse_Matrix;



/* --------- These functions are essential --------- */
void cross(double x[], double y[], double z[]);
double dot(double x[], double y[]);
double sqlength(double x[]);
double length(double x[]);
void vec_sub(double x[], double y[], double z[]);
int norm_axis(int axis);
void righthand2(double *VAL,double **VEC,int n);
int **unit_imatrix(long lo,long hi);
double ***zero_d3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh);
double **zero_dmatrix(long nrl,long nrh,long ncl,long nch);


/* ---------- Essential Numerical Recipes routines ------------- */
unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void deigsrt(double d[], double **v, int n);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void nrerror(char error_text[]);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
double *dvector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
int *ivector(long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);

void deigsrt(double d[], double **v, int n);
double dpythag(double a, double b);
void dsvdcmp(double **a, int m, int n, double w[], double **v);




/* "buildhessian" constructs a block Hessian and associated projection matrix 
   by application of the ANM.  Atomic coordinates and block definitions are 
   provided in 'coords' and 'blocks'; ANM parameters are provided in 'cutoff' 
   and 'gamma'.  On successful termination, the block Hessian is stored in 
   'hessian', and the projection matrix between block and all-atom spaces is 
   in 'projection'. */
static PyObject *buildhessian(PyObject *self, PyObject *args, PyObject *kwargs) {
	PyArrayObject *coords, *hessian;
	double *raw_XYZ, *hess;
	double cutoff = 15., gamma = 1.;
	int N;
	int i, j, k, x, y;
	
	static char *kwlist[] = {"coords", "hessian",
				"natoms", "cutoff", "gamma", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOi|dd", kwlist,
					&coords, &hessian, &N, 
					&cutoff, &gamma))
		return NULL;
	
	raw_XYZ = (double *) PyArray_DATA(coords);
	hess = (double *) PyArray_DATA(hessian);
	
	// cast 1D raw_XYZ array to 2D XYZ array
	double** XYZ = dmatrix(0, N - 1, 0, 2);
	for(i = 0; i < N; i++)
		for(j=0; j < 3; j++)
			XYZ[i][j] = raw_XYZ[i * 3 + j];
			
	// initialize bbHessian matrix
	int N3 = 3 * N;
	double** Hbb = dmatrix(0, N3 - 1, 0, N3 - 1);
	// vectorized version of all bonds.
	int l,m,n,d,ii,ax1,ax2;
	double bi[3], bj[3];
	for (l = 0; l < N; l++)
		for (m = l + 1; m < N; m++)
			for (n = m + 1; n < N; n++)
				for (d=0;d<3;d++){
					if (d==0){
						i=l;
						j=m;
						k=n;
					}	
					else if (d==1){
						i=m;
						j=n;
						k=l;
					}
					else {
						i=n;
						j=l;
						k=m;
					}
					for (ii=0;ii<3;ii++){
						bi[ii]=XYZ[i][ii]-XYZ[k][ii];
						bj[ii]=XYZ[j][ii]-XYZ[k][ii];
					}
					double bi_len = length(bi);
					double bj_len = length(bj);
					//printf("%d\t%d\t%d\t%lf\t%lf\n", i, j, k, bi_len, bj_len);
					if (length(bi)<cutoff && length(bj)<cutoff)
						{
						double theta = acos(dot(bi,bj)/bi_len/bj_len);
						double denom = square(bi_len*bj_len);
						double sq_ij = square(bi_len);
						double sq_jk = square(bj_len);
						double quad_ij = square(square(bi_len));
						double quad_jk = square(square(bj_len));
						double dot_ij_jk = dot(bi,bj);
						double sq_dot_ij_jk = square(dot(bi,bj));
						double cot = square(1. / tan(theta));
						printf("%d\t%d\t%d\t%lf\n", i, j, k, theta);
						if ((theta - M_PI < 1e-5 && theta - M_PI > -1e-5) 
							|| (theta < 1e-5 && theta > -1e-5))
						{
							// i, i update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
								{
									ax1 = norm_axis(x - 1);
									ax2 = norm_axis(x + 1);
									Hbb[3*i+x][3*i+y] = (square(bj[ax1]) + square(bj[ax2]))/denom;
								}
								else 
									Hbb[3*i+x][3*i+y] = 0.5 * (-2 * bj[x]*bj[y])/denom;
									
                            // j, j update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
								{
									ax1 = norm_axis(x - 1);
									ax2 = norm_axis(x + 1);
									Hbb[3*j+x][3*j+y] = (square(bi[ax1]) + square(bi[ax2]))/denom;
								}
								else 
									Hbb[3*j+x][3*j+y] = 0.5 * (-2 * bi[x] * bi[y])/denom;

                            // k, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
								{
									ax1 = norm_axis(x - 1);
									ax2 = norm_axis(x + 1);
									Hbb[3*k+x][3*k+y] = (
														  square(bi[ax1]) - 2*bi[ax1] * bj[ax1] + square(bj[ax1])
														+ square(bi[ax2]) - 2*bi[ax2] * bj[ax2] + square(bj[ax2])
														)/denom;
								}
								else 
									Hbb[3*k+x][3*k+y] = 0.5 * (
														2 * bi[y] * bj[x]-bi[x] * bi[y] 
                                    					  + bi[x] * bj[y]-bj[0] * bj[x]
														)/denom;
							
							// i, j update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
								{
									ax1 = norm_axis(x - 1);
									ax2 = norm_axis(x + 1);
									Hbb[3*i+x][3*j+y] = (- bi[ax1]*bj[ax1] + bi[ax2]*bj[ax2])/denom;
								}
								else 
									Hbb[3*i+x][3*j+y] = (bi[x] * bj[y])/denom;
								
							// i, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
								{
									ax1 = norm_axis(x - 1);
									ax2 = norm_axis(x + 1);
									Hbb[3*i+x][3*k+y] = ( 
														  bi[ax1] * bj[ax1]
														+ bi[ax2] * bj[ax2]
														- square(bj[ax1])
														- square(bj[ax2])
														)/denom;
								}
								else 
									Hbb[3*i+x][3*k+y] = (bj[x] * bj[y] - bi[x] * bj[y])/denom;

							// j, k update
								for (x = 0; x < 3; x++)
								for (y = 0; y < 3; y++)
								if (x == y)
								{
									ax1 = norm_axis(x - 1);
									ax2 = norm_axis(x + 1);
									Hbb[3*j+x][3*k+y] = ( 
														  bi[ax1] * bj[ax1]
														+ bi[ax2] * bj[ax2]
														- square(bi[ax1])
														- square(bi[ax2])
														)/denom;
								}
								else 
									Hbb[3*j+x][3*k+y] = (bi[x] * bi[y] - bi[y] * bj[x])/denom;
									
								}
						else if (theta - M_PI/2 < 1e-5 && theta - M_PI/2 > -1e-5)
						{
								// i, i update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*i+x][3*i+y] = (square(bj[x]))/denom;
								else 
									Hbb[3*i+x][3*i+y] = 0.5 * (2 * bj[x] * bj[y])/denom;
							
							// j, j update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*j+x][3*j+y] = (square(bi[x]))/denom;
								else 
									Hbb[3*j+x][3*j+y] = 0.5 * (2 * bi[x] * bi[y])/denom;
							
							// k, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*k+x][3*k+y] = (
														  square(bi[x])
														+ 2 * bi[x] * bj[x]
														+ square(bj[x])
														)/denom;
								else 
									Hbb[3*k+x][3*k+y] = (
														  bi[x]*bi[y] + bi[y]*bj[x] 
														+ bi[x]*bj[y] + bj[x]*bj[y]
														)/denom;
							
							// i, j update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*i+x][3*j+y] = (bi[x] * bj[x])/denom;
								else 
									Hbb[3*i+x][3*j+y] = (bi[y] * bj[x])/denom;
									
							// i, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*i+x][3*k+y] = (
														- bi[x] * bj[x]
														+ square(bj[x])
														)/denom;
								else 
									Hbb[3*i+x][3*k+y] = (
														- bi[y] * bj[x]
														+ bj[x] * bj[y]
														)/denom;
							
							// j, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*j+x][3*k+y] = (- bj[x]*bi[x] + square(bi[x]))/denom;
								else 
									Hbb[3*j+x][3*k+y] = (- bj[y]*bi[x] + bi[x]*bi[y])/denom;
						}
						else {
							// i, i update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*i+x][3*i+y] =  cot * (
														  square(bi[x])/quad_ij 
														+ square(bj[x])/sq_dot_ij_jk
														- 2 * bi[x] * bj[x]/dot_ij_jk/sq_ij);
								else 
									Hbb[3*i+x][3*i+y] += 0.5 * cot * (
														  2 * bi[x] * bi[y]/quad_ij 
														+ 2 * bj[x] * bj[y]/sq_dot_ij_jk
														- 2 * bi[x] * bj[y]/dot_ij_jk/sq_ij 
														- 2 * bi[y] * bj[x]/dot_ij_jk/sq_ij);
							
							// j, j update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*j+x][3*j+y] =  cot * (
														  square(bj[x])/quad_jk 
														+ square(bi[x])/sq_dot_ij_jk
														- 2 * bi[x] * bj[x]/dot_ij_jk/sq_jk); 
								else 
									Hbb[3*j+x][3*j+y] += 0.5 * cot * (
														  2 * bj[x] * bj[y]/quad_jk 
														+ 2 * bi[x] * bi[y]/sq_dot_ij_jk 
														- 2 * bi[x] * bj[y]/dot_ij_jk/sq_jk 
														- 2 * bi[y] * bj[x]/dot_ij_jk/sq_jk);
							
							// k, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*k+x][3*k+y] =  cot * (
														  square(bi[x])/quad_ij
														+ square(bi[x])/sq_dot_ij_jk
														+ 2 * bi[x] * bj[x]/sq_dot_ij_jk
														+ square(bj[x])/sq_dot_ij_jk
														- 2 * square(bi[x])/dot_ij_jk/sq_ij
														- 2 * bi[x] * bj[x]/dot_ij_jk/sq_ij
														+ square(bj[x])/quad_jk
														+ 2 * bi[x] * bj[x]/sq_ij/sq_jk
														- 2 * bi[x] * bj[x]/dot_ij_jk/sq_jk
														- 2 * square(bj[x])/dot_ij_jk/sq_jk);
								else 
									Hbb[3*k+x][3*k+y] += 0.5 * cot * (
														  2 * bi[x] * bi[y]/quad_ij
														+ 2 * bi[x] * bi[y]/sq_dot_ij_jk
														+ 2 * bi[y] * bj[x]/sq_dot_ij_jk
														+ 2 * bi[x] * bj[y]/sq_dot_ij_jk
														+ 2 * bj[x] * bj[y]/sq_dot_ij_jk
														- 4 * bi[x] * bi[y]/dot_ij_jk/sq_ij
														- 2 * bi[y] * bj[x]/dot_ij_jk/sq_ij
														- 2 * bi[x] * bj[y]/dot_ij_jk/sq_ij
														+ 2 * bj[x] * bj[y]/quad_jk
														+ 2 * bi[y] * bj[x]/sq_ij/sq_jk
														+ 2 * bi[x] * bj[y]/sq_ij/sq_jk
														- 2 * bi[y] * bj[x]/dot_ij_jk/sq_jk
														- 2 * bi[x] * bj[y]/dot_ij_jk/sq_jk
														- 4 * bj[x] * bj[y]/dot_ij_jk/sq_jk);
														
							// i, j update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y){
									Hbb[3*i+x][3*j+y] += 0.5 * cot * (
														  2 * bi[x] * bj[x]/sq_dot_ij_jk
														- 2 * square(bi[x])/dot_ij_jk/sq_ij
														+ 2 * bi[x] * bj[x]/sq_ij/sq_jk
														- 2 * square(bj[x])/dot_ij_jk/sq_jk);
								}
								else 
									Hbb[3*i+x][3*j+y] += 0.5 * cot * (
														  2 * bi[y] * bj[x]/sq_dot_ij_jk
														- 2 * bi[x] * bi[y]/dot_ij_jk/sq_ij
														+ 2 * bi[x] * bj[y]/sq_ij/sq_jk
														- 2 * bj[x] * bj[y]/dot_ij_jk/sq_jk);
														
							// i, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*i+x][3*k+y] += 0.5 * cot * (
														- 2 * square(bi[x])/quad_ij
														- 2 * bi[x] * bj[x]/sq_dot_ij_jk
														- 2 * square(bj[x])/sq_dot_ij_jk
														+ 2 * square(bi[x])/sq_ij/dot_ij_jk
														+ 4 * bi[x] * bj[x]/sq_ij/dot_ij_jk
														- 2 * bi[x] * bj[x]/sq_ij/sq_jk
														+ 2 * square(bj[x])/sq_jk/dot_ij_jk);
								else 
									Hbb[3*i+x][3*k+y] += 0.5 * cot * (
														- 2 * bi[x] * bi[y]/quad_ij
														- 2 * bi[y] * bj[x]/sq_dot_ij_jk
														- 2 * bj[x] * bj[y]/sq_dot_ij_jk
														+ 2 * bi[x] * bi[y]/sq_ij/dot_ij_jk
														+ 2 * bi[y] * bj[x]/sq_ij/dot_ij_jk
														+ 2 * bi[x] * bj[y]/sq_ij/dot_ij_jk
														- 2 * bi[x] * bj[y]/sq_ij/sq_jk
														+ 2 * bj[x] * bj[y]/dot_ij_jk/sq_jk);
							
							// j, k update
							for (x = 0; x < 3; x++)
							for (y = 0; y < 3; y++)
								if (x == y)
									Hbb[3*j+x][3*k+y] += 0.5 * cot * (
														- 2 * square(bi[x])/sq_dot_ij_jk
														- 2 * bi[x] * bj[x]/sq_dot_ij_jk
														+ 2 * square(bi[x])/sq_ij/dot_ij_jk
														- 2 * square(bj[x])/quad_jk
														- 2 * bi[x] * bj[x]/sq_ij/sq_jk
														+ 4 * bi[x] * bj[x]/dot_ij_jk/sq_jk
														+ 2 * square(bj[x])/dot_ij_jk/sq_jk);
								else 
									Hbb[3*j+x][3*k+y] += 0.5 * cot * (
														- 2 * bj[x] * bj[y]/quad_jk
														- 2 * bi[x] * bi[y]/sq_dot_ij_jk
														- 2 * bi[x] * bj[y]/sq_dot_ij_jk
														+ 2 * bi[x] * bi[y]/dot_ij_jk/sq_ij
														- 2 * bi[y] * bj[x]/sq_ij/sq_jk
														+ 2 * bi[y] * bj[x]/sq_jk/dot_ij_jk
														+ 2 * bi[x] * bj[y]/sq_jk/dot_ij_jk
														+ 2 * bj[x] * bj[y]/sq_jk/dot_ij_jk);
						}
					}
				}
	
	// the Hbb is symmetric
	for(i = 0; i < N3; i++)
    	for(j = i+1; j < N3; j++){
    		if (Hbb[i][j]>1e-6 || Hbb[i][j]<-1e-6)
				Hbb[j][i] = Hbb[i][j];
			else 
				Hbb[i][j] = Hbb[j][i];
		}
	// for(i = 0; i < N3; i++){
 //    	for(j = 0; j < N3; j++){
 //    		printf("%lf\t",Hbb[i][j]);
	// 		//Hbb[j][i] = Hbb[i][j];
	// 	}
	// 	printf("\n");		  
	// }
	// // hess = H_anm + H_bb
	for(i = 0; i < N3; i++)
    	for(j = 0; j < N3; j++)
			hess[N3 * i + j] += Hbb[i][j];
	  
	free_dmatrix(XYZ, 0, N-1, 0, 2);
	free_dmatrix(Hbb, 0, N3 - 1, 0, N3 - 1);

	Py_RETURN_NONE;
}


static PyMethodDef bbenmtools_methods[] = {

    {"buildhessian",  (PyCFunction)buildhessian,
     METH_VARARGS | METH_KEYWORDS,
     "Build Hessian matrix."},

    {NULL, NULL, 0, NULL}
};



#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef bbenmtools = {
        PyModuleDef_HEAD_INIT,
        "bbenmtools",
        "bbENM.",
        -1,
        bbenm_methods,
};
PyMODINIT_FUNC PyInit_bbenmtools(void) {
    import_array();
    return PyModule_Create(&bbenmtools);
}
#else
PyMODINIT_FUNC initbbenmtools(void) {

    Py_InitModule3("bbenmtools", bbenmtools_methods,
        "bbENM.");
    
    import_array();
}
#endif

/* "cross" TAKES THE 3D CROSS PRODUCT OF ITS ARGUMENTS. */
void cross(double x[], double y[], double z[])
{
	z[0] = x[1]*y[2] - x[2]*y[1];
	z[1] = x[2]*y[0] - x[0]*y[2];
	z[2] = x[0]*y[1] - x[1]*y[0];
}

/* "dot" TAKES THE 3D DOT PRODUCT OF ITS ARGUMENTS. */
double dot(double x[], double y[])
{
	return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

/* "norm_axis" NORMALIZES AXIS NUMBER such that 0-x, 1-y, 2-z. */
int norm_axis(int axis)
{
	return (axis + 3) % 3;
}

/* "vec_sub" TAKES THE SUBTRACTION OF ITS ARGUMENTS. */
void vec_sub(double x[], double y[], double z[])
{
	int i;
	for (i = 0; i < 3; i++)
		z[i] = x[i] - y[i];
}

/* "sqlength" CALCULATES THE SQUARED LENGTH OF AN ARRAY. */
double sqlength(double x[])
{
	return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

/* "length" CALCULATES THE LENGTH OF AN ARRAY. */
double length(double x[])
{
	return sqrt(sqlength(x));
}

/* "righthand2" MAKES SURE THAT THE EIGENVECTORS
   FORM A RIGHT-HANDED COORDINATE SYSTEM */
void righthand2(double *VAL,double **VEC,int n)
{
	double A[3],B[3],C[3],CP[3],dot=0.0;
	int i;
	
	/* FIND THE CROSS PRODUCT OF THE FIRST TWO EIGENVECTORS */
	for(i=0;i<3;i++){
		A[i]=VEC[i+1][1];
		B[i]=VEC[i+1][2];
		C[i]=VEC[i+1][3];}
	cross(A,B,CP);
	
	/* PROJECT IT ON THE THIRD EIGENVECTOR */
	for(i=0; i<3; i++)
		dot+=C[i]*CP[i];
	if(dot<0.0)
		for(i=1;i<=3;i++)
		VEC[i][3]=-VEC[i][3];
}



/* "unit_imatrix" ALLOCATES MEMORY FOR A UNIT MATRIX */
int **unit_imatrix(long lo,long hi)
{
	static int **M;
	int i,j;
	
	M=imatrix(lo,hi,lo,hi);
	for(i=lo;i<=hi;i++){
		M[i][i]=1;
		for(j=i+1;j<=hi;j++)
		M[i][j]=M[j][i]=0;
	}
	return M;
}

/* "zero_dmatrix" ALLOCATES MEMORY FOR A
   DOUBLE MATRIX AND INITIALIZES IT TO ZERO */
double **zero_dmatrix(long nrl,long nrh,long ncl,long nch)
{
	static double **M;
	int i,j;
	
	M=dmatrix(nrl,nrh,ncl,nch);
	for(i=nrl;i<=nrh;i++)
		for(j=ncl;j<=nch;j++)
		M[i][j]=0.0;
	return M;
}






/* ------------ Numerical Recipes Routines ---------------- */
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}



double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+1)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += 1;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+1)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += 1;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+1)*sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += 1;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}



void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-1));
	free((FREE_ARG) (t[nrl]+ncl-1));
	free((FREE_ARG) (t+nrl-1));
}

void dsvdcmp(double **a, int m, int n, double w[], double **v)
{
	double dpythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
	static int maxits=100;

	rv1=dvector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=maxits;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == maxits) nrerror("no convergence in many dsvdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);
}

double dpythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

void deigsrt(double d[], double **v, int n)
{
	int k,j,i;
	double p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}
