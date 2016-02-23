#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "math.h"
#include "stdio.h"

#define NR_END 1
#define FREE_ARG char*
#define square(x) x * x

double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
double **zero_dmatrix(long nrl,long nrh,long ncl,long nch);
void nrerror(char error_text[]);

static PyObject *calcSM(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject *coords, *sm, *eigvecs, *eigvals;
  int numCA, i, j, k, nmodes;
  double *XYZ, *SM, *lambda, *U, kbt=1.;
  double **stiff_matrix;
  double r_ij, x_ij, y_ij, z_ij;
  static char *kwlist[] = {"coords", "sm", "eigvecs", "eigvals",
          "natoms","n_modes",
          "k_B_x_T",NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOii|d", kwlist, 
          &coords, &sm, &eigvecs, &eigvals,
          &numCA, &nmodes,
          &kbt))
    return NULL;

  XYZ = (double *) PyArray_DATA(coords);
  SM = (double *) PyArray_DATA(sm);
  U = (double *) PyArray_DATA(eigvecs);
  lambda = (double *) PyArray_DATA(eigvals);

  stiff_matrix=dmatrix(0,numCA-1,0,numCA-1);

  //printf("%lf,%lf,%lf\n",U[0],U[1],U[2]);

  for (i=0; i<numCA; i++){
    for (j=i+1; j<numCA; j++){
      r_ij = sqrt((XYZ[j*3]-XYZ[i*3])*(XYZ[j*3]-XYZ[i*3])+\
        (XYZ[j*3+1]-XYZ[i*3+1])*(XYZ[j*3+1]-XYZ[i*3+1])+\
        (XYZ[j*3+2]-XYZ[i*3+2])*(XYZ[j*3+2]-XYZ[i*3+2]));
      x_ij = (XYZ[j*3]-XYZ[i*3])/r_ij;
      y_ij = (XYZ[j*3+1]-XYZ[i*3+1])/r_ij;
      z_ij = (XYZ[j*3+2]-XYZ[i*3+2])/r_ij;
      double u_ij_sup_k[3]={0.0, 0.0, 0.0};
      double d_ij_sup_k=0.0;
      double sum1=0.0;
      double sum2=0.0;
      double cos_alpha_ij=0.0;
      for(k=6; k<nmodes; k++){
      //      u_ij_sup_k[0]=(eigvecs[k][ind_3j  ]-eigvecs[k][ind_3i  ]);
        printf("%d\n",k);
        u_ij_sup_k[0]=(U[(k)*3*numCA+j*3]-U[(k)*3*numCA+i*3]);

      //      u_ij_sup_k[1]=(eigvecs[k][ind_3j+1]-eigvecs[k][ind_3i+1]);
        u_ij_sup_k[1]=(U[(k)*3*numCA+j*3+1]-U[(k)*3*numCA+i*3+1]);

      //      u_ij_sup_k[2]=(eigvecs[k][ind_3j+2]-eigvecs[k][ind_3i+2]);
        u_ij_sup_k[2]=(U[(k)*3*numCA+j*3+2]-U[(k)*3*numCA+i*3+2]);

        cos_alpha_ij=(  (x_ij*u_ij_sup_k[0]) +\
          (y_ij*u_ij_sup_k[1]) +\
          (z_ij*u_ij_sup_k[2])  );

        d_ij_sup_k=sqrt(kbt/lambda[k])*cos_alpha_ij;

       // if (i == 0 && j==1 && k==6)
         // printf("%d,%d,%lf,%lf,%lf,%lf,%lf\n",i,j,lambda[k],d_ij_sup_k,u_ij_sup_k[0],u_ij_sup_k[1],u_ij_sup_k[2]);
      
        sum1+=fabs(lambda[k]*d_ij_sup_k);
        sum2+=fabs(d_ij_sup_k);
       
      }
      stiff_matrix[i][j]=sum1/sum2;
      stiff_matrix[j][i]=stiff_matrix[i][j];
     
    }
  }
  for (i=0;i<numCA;i++)
    for (j=0;j<numCA;j++)
      SM[i*numCA+j]=stiff_matrix[i][j];

  free_dmatrix(stiff_matrix,0,numCA-1,0,numCA-1);

  Py_RETURN_NONE;
}

static PyMethodDef smtools_methods[] = {

    {"calcSM",  (PyCFunction)calcSM,
     METH_VARARGS | METH_KEYWORDS,
     "Build stiffness matrix."},

    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef smtools = {
        PyModuleDef_HEAD_INIT,
        "smtools",
        "SM tools.",
        -1,
        smtools_methods,
};
PyMODINIT_FUNC PyInit_smtools(void) {
    import_array();
    return PyModule_Create(&smtools);
}
#else
PyMODINIT_FUNC initsmtools(void) {

    Py_InitModule3("smtools", smtools_methods,
        "SM tools.");

    import_array();
}
#endif

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

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

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