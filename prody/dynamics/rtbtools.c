/*****************************************************************************/
/*                                                                           */
/*                   Tools for RTB calculations in ProDy.                    */
/*                                                                           */
/*****************************************************************************/
/* Author: Tim Lezon, Ahmet Bakan */
#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"


/* ---------- Numerical Recipes-specific definitions and macros ---------- */
#define NR_END 1
#define FREE_ARG char*

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
int bless_from_tensor(double **HB,double ***HT,int **CT,int nblx);
int calc_blessian_mem(PDB_File *PDB,dSparse_Matrix *PP1,int nres,int nblx,
		      int elm,double **HB,double cut,double gam,double scl,
		      double mlo,double mhi);
void copy_dsparse(dSparse_Matrix *A,dSparse_Matrix *B,int lo,int hi);
void copy_prj_ofst(dSparse_Matrix *PP,double *proj,int elm,int bdim);
void cross(double x[], double y[], double z[]);
int dblock_projections2(dSparse_Matrix *PP,PDB_File *PDB,
			int nres,int nblx,int bmx);
void dsort_PP2(dSparse_Matrix *MM,int n,int idx);
int find_contacts1(int **CT,PDB_File *PDB,int nres,int nblx,double cut);
void hess_superrow_mem(double **HR,int **CT,PDB_File *PDB,int nres,
		       int who,double cut,double gam,double mscl,double mlo,
		       double mhi);
void init_bst(int *BST,dSparse_Matrix *PP,int elm,int n,int idx);
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
  PDB_File PDB;
  dSparse_Matrix PP,HH;
  PyArrayObject *coords, *blocks, *hessian, *projection;
  double *XYZ,*hess,*proj;
  long *BLK;
  double **HB;
  double cutoff = 15., gamma = 1., scl=1., mlo=1., mhi=-1.;
  int natm, nblx, bmx;
  int hsize,elm,bdim,i,j;

  static char *kwlist[] = {"coords", "blocks", "hessian", "projection",
			   "natoms", "nblocks", "maxsize", "cutoff",
			   "gamma", "scale", "memlo", "memhi", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOiii|ddddd", kwlist,
				   &coords, &blocks, &hessian, &projection,
				   &natm, &nblx, &bmx, &cutoff, &gamma, &scl,
				   &mlo, &mhi))
    return NULL;

  XYZ = (double *) PyArray_DATA(coords);
  BLK = (long *) PyArray_DATA(blocks);
  hess = (double *) PyArray_DATA(hessian);
  proj = (double *) PyArray_DATA(projection);



  /* First allocate a PDB_File object to hold the coordinates and block
     indices of the atoms.  This wastes a bit of memory, but it prevents
     the need to re-write all of the RTB functions that are used in
     standalone C code. */
  PDB.atom=malloc((size_t)((natm+2)*sizeof(Atom_Line)));
  if(!PDB.atom) return PyErr_NoMemory();
  for(i=1;i<=natm;i++){
    PDB.atom[i].model=BLK[i-1];
    for(j=0;j<3;j++)
      PDB.atom[i].X[j]=XYZ[j*natm+i-1];
  }



  /* Find the projection matrix */
  hsize = 18*bmx*nblx > 12*natm ? 12*natm : 18*bmx*nblx;
  HH.IDX=imatrix(1,hsize,1,2);
  HH.X=dvector(1,hsize);
  elm=dblock_projections2(&HH,&PDB,natm,nblx,bmx);
  PP.IDX=imatrix(1,elm,1,2);
  PP.X=dvector(1,elm);
  for(i=1;i<=elm;i++){
    PP.IDX[i][1]=HH.IDX[i][1];
    PP.IDX[i][2]=HH.IDX[i][2];
    PP.X[i]=HH.X[i];
  }
  free_imatrix(HH.IDX,1,hsize,1,2);
  free_dvector(HH.X,1,hsize);
  dsort_PP2(&PP,elm,1);


  /* Calculate the block Hessian */
  HB=dmatrix(1,6*nblx,1,6*nblx);
  bdim=calc_blessian_mem(&PDB,&PP,natm,nblx,elm,HB,cutoff,gamma,scl,mlo,mhi);


  /* Cast the block Hessian and projection matrix into 1D arrays. */
  copy_prj_ofst(&PP,proj,elm,bdim);
  for(i=1;i<=bdim;i++)
    for(j=1;j<=bdim;j++)
      hess[bdim*(i-1)+j-1]=HB[i][j];


  free(PDB.atom);
  free_imatrix(PP.IDX,1,elm,1,2);
  free_dvector(PP.X,1,elm);
  free_dmatrix(HB,1,6*nblx,1,6*nblx);


  Py_RETURN_NONE;
}


static PyMethodDef rtbtools_methods[] = {

    {"buildhessian",  (PyCFunction)buildhessian,
     METH_VARARGS | METH_KEYWORDS,
     "Build Hessian matrix and projections."},

    {NULL, NULL, 0, NULL}
};



#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef rtbtools = {
        PyModuleDef_HEAD_INIT,
        "rtbtools",
        "RTB tools.",
        -1,
        rtbtools_methods,
};
PyMODINIT_FUNC PyInit_rtbtools(void) {
    import_array();
    return PyModule_Create(&rtbtools);
}
#else
PyMODINIT_FUNC initrtbtools(void) {

    Py_InitModule3("rtbtools", rtbtools_methods,
        "RTB tools.");

    import_array();
}
#endif



/* "bless_from_tensor" transfers the block Hessian
   from  the tensor HT into the array HB */
int bless_from_tensor(double **HB,double ***HT,int **CT,int nblx)
{
  int *I1,*I2,i,j,p,sb,ii,jj,max,a,b,imx;

  max=6*nblx;
  I1=ivector(1,max);
  I2=ivector(1,max);

  /* I1[i]==i iff there is a non-zero element in column i
     (removes zeroes that are caused by single-node blocks) */
  for(i=1;i<=max;i++){
    I1[i]=0;
    for(j=i;j<=max;j++)
      HB[i][j]=HB[j][i]=0.0;
  }
  for(ii=1;ii<=nblx;ii++){
    for(i=1;i<=6;i++){
      for(jj=ii;jj<=nblx;jj++){
	sb=CT[ii][jj];
	if(sb!=0){
	  p = jj==ii ? i : 1;
	  for(j=p;j<=6;j++)
	    if(HT[sb][i][j]!=0)
	      I1[6*(jj-1)+j]=6*(jj-1)+j;
	}
      }
    }
  }

  /* If I1[i]!=0, then I2[i] is a sequential index */
  imx=0;
  for(i=1;i<=max;i++){
    if(I1[i]!=0) imx++;
    I2[i]=imx;
  }

  for(ii=1;ii<=nblx;ii++){
    for(i=1;i<=6;i++){
      for(jj=ii;jj<=nblx;jj++){
	sb=CT[ii][jj];
	if(sb!=0){
	  p = jj==ii ? i : 1;
	  for(j=p;j<=6;j++)
	    if(HT[sb][i][j]!=0){
	      a=I2[6*(ii-1)+i];
	      b=I2[6*(jj-1)+j];
	      HB[a][b]=HB[b][a]=HT[sb][i][j];
	    }
	}
      }
    }
  }
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
  return imx;
}


/* "calc_blessian_mem" calculates the block Hessian. */
int calc_blessian_mem(PDB_File *PDB,dSparse_Matrix *PP1,int nres,int nblx,
		      int elm,double **HB,double cut,double gam,double scl,
		      double mlo,double mhi)
{
  dSparse_Matrix *PP2;
  double **HR,***HT;
  int **CT,*BST1,*BST2;
  int ii,i,j,k,p,q,q1,q2,ti,tj,bi,bj,sb,nc,out;


  /* ------------------- INITIALIZE LOCAL VARIABLES ------------------- */

  /* HR holde three rows (corresponding to 1 residue) of the full Hessian */
  HR=zero_dmatrix(1,3*nres,1,3);

  /* CT is an array of contacts between blocks */
  CT=unit_imatrix(0,nblx);

  /* Copy PP1 to PP2 and sort by second element */
  PP2=(dSparse_Matrix *)malloc((size_t)sizeof(dSparse_Matrix));
  PP2->IDX=imatrix(1,elm,1,2);
  PP2->X=dvector(1,elm);
  copy_dsparse(PP1,PP2,1,elm);
  dsort_PP2(PP2,elm,2);

  /* BST1: for all j: BST1[i]<=j<BST[i+1], PP1->IDX[j][1]=i */
  /* BST2: for all j: BST2[i]<=j<BST2[i+1], PP2->IDX[j][2]=i */
  BST1=ivector(1,3*nres+1);
  BST2=ivector(1,6*nblx+1);
  init_bst(BST1,PP1,elm,3*nres+1,1);
  init_bst(BST2,PP2,elm,6*nblx+1,2);
  /* ------------------- LOCAL VARIABLES INITIALIZED ------------------ */



  /* ------------- FIND WHICH BLOCKS ARE IN CONTACT --------------- */
  nc=find_contacts1(CT,PDB,nres,nblx,cut);


  /* Allocate a tensor for the block Hessian */
  HT=zero_d3tensor(1,nc,1,6,1,6);


  /* Calculate each super-row of the full Hessian */
  for(ii=1;ii<=nres;ii++){

    if(PDB->atom[ii].model!=0){

      /* ----------------- FIND SUPER-ROW OF FULL HESSIAN --------------- */
      hess_superrow_mem(HR,CT,PDB,nres,ii,cut,gam,scl,mlo,mhi);


      /* Update elements of block hessian */
      q1=BST1[3*(ii-1)+2];
      q2=BST1[3*(ii-1)+3];
      /* Sum over elements of projection matrix corresponding to residue ii:
	 for each k in the following loop, PP1->IDX[k][1]==3*ii + 0,1,2 */
      for(k=BST1[3*ii-2];k<BST1[3*ii+1];k++){
	if(k<q1) q=1;
	else if(k<q2) q=2;
	else q=3;
	i=PP1->IDX[k][2];
	bi=(i-1)/6+1;
	ti=i-6*(bi-1);
	/* Sum over all elements of projection matrix with column j>=i */
	for(p=BST2[i];p<=elm;p++){
	  j=PP2->IDX[p][2];
	  bj=(j-1)/6+1;
	  sb=CT[bi][bj];
	  if(i<=j && sb!=0){  /* the first condition should ALWAYS hold */
	    tj=j-6*(bj-1);
	    HT[sb][ti][tj]+=(PP1->X[k]*PP2->X[p]*HR[PP2->IDX[p][1]][q]);
	  }
	}
      }
    }
  }


  /* Print the block Hessian in sparse format */
  out=bless_from_tensor(HB,HT,CT,nblx);

  /* Free up memory */
  free_dmatrix(HR,1,3*nres,1,3);
  free_d3tensor(HT,1,nc,1,6,1,6);
  free_imatrix(CT,0,nblx,0,nblx);
  free_ivector(BST1,1,3*nres+1);
  free_ivector(BST2,1,6*nblx+1);
  free_imatrix(PP2->IDX,1,elm,1,2);
  free_dvector(PP2->X,1,elm);
  return out;
}



/* "copy_prj_ofst" copies the projection matrix to the Python object, 
   offsetting the column elements if there are blocks with fewer than 
   six degrees of freedom. 
   PP[i].X = proj[6*nblx*(PP[i].IDX[1]-1)+PP[i].IDX[2]-1
*/
void copy_prj_ofst(dSparse_Matrix *PP,double *proj,int elm,int bdim)
{
  int *I1,*I2,max=0,i,j=0;

  for(i=1;i<=elm;i++)
    if(PP->IDX[i][2]>max)
      max=PP->IDX[i][2];
  I1=ivector(1,max);
  I2=ivector(1,max);
  for(i=1;i<=max;i++) I1[i]=0;
  for(i=1;i<=elm;i++)
    I1[PP->IDX[i][2]]=PP->IDX[i][2];
  for(i=1;i<=max;i++){
    if(I1[i]!=0) j++;
    I2[i]=j;
  }
  for(i=1;i<=elm;i++)
    if(PP->X[i]!=0.0)
      proj[bdim*(PP->IDX[i][1]-1) + I2[PP->IDX[i][2]] - 1] = PP->X[i];
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}


/* "copy_dsparse" COPIES ELEMENTS lo THROUGH hi
   OF SPARSE MATRIX 'A' TO SPARSE MATRIX 'B' */
void copy_dsparse(dSparse_Matrix *A,dSparse_Matrix *B,int lo,int hi)
{
  int i;

  for(i=lo;i<=hi;i++){
    B->IDX[i][1]=A->IDX[i][1];
    B->IDX[i][2]=A->IDX[i][2];
    B->X[i]=A->X[i];
  }
}


/* "cross" TAKES THE 3D CROSS PRODUCT OF ITS ARGUMENTS. */
void cross(double x[], double y[], double z[])
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}


/* "dblock_projections2" CALCULATES THE PROJECTION
   FROM FULL RESIDUE SPACE TO RIGID BLOCK SPACE */
int dblock_projections2(dSparse_Matrix *PP,PDB_File *PDB,
			int nres,int nblx,int bmx)
{
  double **X,**I,**IC,*CM,*W,**A,**ISQT;
  double x,tr,dd,df;
  int *IDX,nbp,b,i,j,k,ii,jj,aa,bb,elm;


  /* INITIALIZE BLOCK ARRAYS */
  elm=0;
  X=dmatrix(1,bmx,1,3);
  IDX=ivector(1,bmx);
  CM=dvector(1,3);
  I=dmatrix(1,3,1,3);
  IC=dmatrix(1,3,1,3);
  W=dvector(1,3);
  A=dmatrix(1,3,1,3);
  ISQT=dmatrix(1,3,1,3);

  /* CYCLE THROUGH BLOCKS */
  for(b=1;b<=nblx;b++){

    /* CLEAR MATRICES */
    for(j=1;j<=3;j++){
      CM[j]=0.0;
      for(i=1;i<=3;i++) I[i][j]=0.0;
      for(i=1;i<=bmx;i++) X[i][j]=0.0;
    }

    /* STORE VALUES FOR CURRENT BLOCK */
    nbp=0;
    for(i=1;i<=nres;i++){
      if(PDB->atom[i].model==b){
	IDX[++nbp]=i;
	for(j=1;j<=3;j++){
	  x=(double)PDB->atom[i].X[j-1];
	  X[nbp][j]=x;
	  CM[j]+=x;
	}
      }
    }

    /* TRANSLATE BLOCK CENTER OF MASS TO ORIGIN */
    for(j=1;j<=3;j++) CM[j]/=(double)nbp;
    for(i=1;i<=nbp;i++)
      for(j=1;j<=3;j++)
	X[i][j]-=CM[j];

    /* CALCULATE INERTIA TENSOR */
    for(k=1;k<=nbp;k++){
      dd=0.0;
      for(j=1;j<=3;j++){
	df=X[k][j];
	dd+=df*df;
      }
      for(i=1;i<=3;i++){
	I[i][i]+=(dd-X[k][i]*X[k][i]);
	for(j=i+1;j<=3;j++){
	  I[i][j]-=X[k][i]*X[k][j];
	  I[j][i]=I[i][j];
	}
      }
    }

    /* DIAGONALIZE INERTIA TENSOR */
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++)
	IC[i][j]=I[i][j];
    dsvdcmp(IC,3,3,W,A);
    deigsrt(W,A,3);
    righthand2(W,A,3);

    /* FIND ITS SQUARE ROOT */
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++){
	dd=0.0;
	for(k=1;k<=3;k++)
	  dd+=A[i][k]*A[j][k]/sqrt(W[k]);
	ISQT[i][j]=dd;
      }

    /* UPDATE PP WITH THE RIGID MOTIONS OF THE BLOCK */
    tr=1.0/sqrt((double)nbp);
    for(i=1;i<=nbp;i++){

      /* TRANSLATIONS: 3*(IDX[i]-1)+1 = x-COORDINATE OF RESIDUE IDX[i];
	 6*(b-1)+1 = x-COORDINATE OF BLOCK b */
      for(j=1;j<=3;j++){
	elm++;
	PP->IDX[elm][1] = 3*(IDX[i]-1)+j;
	PP->IDX[elm][2] = 6*(b-1)+j;
	PP->X[elm] = tr;
      }

      /* ROTATIONS */
      if(nbp>1){
	for(ii=1;ii<=3;ii++){
	  for(jj=1;jj<=3;jj++){
	    if(jj==1) {aa=2; bb=3;}
	    else if(jj==2) {aa=3; bb=1;}
	    else {aa=1; bb=2;}
	    dd=ISQT[ii][aa]*X[i][bb]-ISQT[ii][bb]*X[i][aa];
	    elm++;
	    PP->IDX[elm][1] = 3*(IDX[i]-1)+jj;
	    PP->IDX[elm][2] = 6*(b-1)+3+ii;
	    PP->X[elm] = dd;
	  }
	}
      }
    }
  }
  free_dmatrix(X,1,bmx,1,3);
  free_ivector(IDX,1,bmx);
  free_dvector(CM,1,3);
  free_dmatrix(I,1,3,1,3);
  free_dmatrix(IC,1,3,1,3);
  free_dvector(W,1,3);
  free_dmatrix(A,1,3,1,3);
  free_dmatrix(ISQT,1,3,1,3);

  return elm;
}


/* "dsort_PP2" SORTS THE PROJECTION MATRIX IN ASCENDING ORDER OF THE
   INDEX 'idx'.  ADAPTED FROM THE NUMERICAL RECIPES 'HEAPSORT' ROUTINE. */
void dsort_PP2(dSparse_Matrix *MM,int n,int idx)
{
  double x;
  int i,ir,j,l,hi,i1,i2,ndx;
  unsigned long rra,*ra;

  if(n<2) return;
  ndx = idx==1 ? 2 : 1;

  /* CREATE A VECTOR TO INDEX THE ELEMENTS OF MM */
  hi=0;
  for(i=1;i<=n;i++)
    if(MM->IDX[i][ndx]>hi)
      hi=MM->IDX[i][ndx];
  ra=lvector(1,n);
  for(i=1;i<=n;i++)
    ra[i]=(long)hi*(MM->IDX[i][idx]-1)+MM->IDX[i][ndx];


  /* SORT */
  l=(n >> 1)+1;
  ir=n;
  for(;;){
    if(l > 1){
      rra=ra[--l];
      i1=MM->IDX[l][idx];
      i2=MM->IDX[l][ndx];
      x=MM->X[l];
    }
    else {
      rra=ra[ir];
      i1=MM->IDX[ir][idx];
      i2=MM->IDX[ir][ndx];
      x=MM->X[ir];
      ra[ir]=ra[1];
      MM->IDX[ir][idx]=MM->IDX[1][idx];
      MM->IDX[ir][ndx]=MM->IDX[1][ndx];
      MM->X[ir]=MM->X[1];
      if (--ir == 1) {
	ra[1]=rra;
	MM->IDX[1][idx]=i1;
	MM->IDX[1][ndx]=i2;
	MM->X[1]=x;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	MM->IDX[i][idx]=MM->IDX[j][idx];
	MM->IDX[i][ndx]=MM->IDX[j][ndx];
	MM->X[i]=MM->X[j];
	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i]=rra;
    MM->IDX[i][idx]=i1;
    MM->IDX[i][ndx]=i2;
    MM->X[i]=x;
  }
  free_lvector(ra,1,n);
}



/* "find_contacts1" FINDS WHICH BLOCKS ARE IN CONTACT, AND ASSIGNS EACH
   PAIR OF CONTACTING BLOCKS A UNIQUE INDEX.  IT RETURNS THE TOTAL NUMBER
   OF CONTACTS BETWEEN BLOCKS. */
int find_contacts1(int **CT,PDB_File *PDB,int nres,int nblx,double cut)
{
  int nc,i,j,k,ii,jj;
  double csq=cut*cut,df,dd;

  for(i=1;i<=nres;i++){
    ii=PDB->atom[i].model;
    for(j=i+1;j<=nres;j++){
      jj=PDB->atom[j].model;

      if(ii!=jj && ii!=0 && jj!=0 && CT[ii][jj]==0){
	dd=0.0;
	for(k=0;k<3;k++){
	  df=(double)PDB->atom[i].X[k]-PDB->atom[j].X[k];
	  dd+=df*df;
	}
	if(dd<csq)
	  CT[ii][jj]=CT[jj][ii]=1;
      }

    }
  }
  nc=0;
  for(i=1;i<=nblx;i++)
    for(j=i;j<=nblx;j++)
      if(CT[i][j]!=0){
	nc++;
	CT[i][j]=CT[j][i]=nc;
      }
  return nc;
}


/* "hess_superrow_mem" calculates the 'who'-th super-row
   of the Hessian, using 'cut' as the cutoff and 'gam' as the
   spring constant for all interactions. */
void hess_superrow_mem(double **HR,int **CT,PDB_File *PDB,int nres,
		       int who,double cut,double gam,double mscl,
		       double mlo,double mhi)
{
  int i,j,k,jj;
  double DX[3],csq=cut*cut,dsq,df;
  double s0,scl;

  s0=pow(mscl,0.25);

  /* Clear the diagonal super-element */
  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
      HR[3*(who-1)+i][j]=0.0;

  /* Calculate the submatrices */
  for(jj=1;jj<=nres;jj++){

    if(jj!=who && PDB->atom[jj].model!=0 &&
       CT[PDB->atom[who].model][PDB->atom[jj].model]!=0){

      dsq=0.0;
      for(k=0;k<3;k++){
	DX[k] = (double)PDB->atom[who].X[k] - PDB->atom[jj].X[k];
	dsq+=(DX[k]*DX[k]);
      }


      if(dsq<csq){

	/* --------- Membrane scaling -------- */
	scl=1.0;
	if(mhi<mlo || (PDB->atom[who].X[2] < mhi && PDB->atom[who].X[2] > mlo)) scl*=s0;
	if(mhi<mlo || (PDB->atom[jj].X[2] < mhi && PDB->atom[jj].X[2] > mlo)) scl*=s0;


	for(i=1;i<=3;i++){
	  for(j=i;j<=3;j++){

	    df=gam*DX[i-1]*DX[j-1]/dsq;


	    /* Strong backbone bonds:
	       NOTE:  *Not currently available! May be implemented later 
	    if((int)fabs(PDB->atom[who].resnum-PDB->atom[jj].resnum)==1 &&
	       PDB->atom[who].chain==PDB->atom[jj].chain)
	      df*=100.0;
	    */


	    /* -------- MEMBRANE RULES -------- */
	    /* Scale lateral components */
	    if(i!=3) df*=scl;
	    if(j!=3) df*=scl;


	    /* Off-diagonal super-elements */
	    HR[3*(jj-1)+i][j]=HR[3*(jj-1)+j][i]=-df;

	    /* Diagonal super-elements */
	    HR[3*(who-1)+i][j]+=df;
	    if(i!=j)
	      HR[3*(who-1)+j][i]+=df;
	  }
	}
      } /* <----- if(dsq<csq) */
      else
	for(i=1;i<=3;i++)
	  for(j=1;j<=3;j++)
	    HR[3*(jj-1)+i][j]=HR[3*(jj-1)+j][i]=0.0;
    } /* <---- if(jj!=who &&...) */
  }
}



/* "init_bst" INITIALIZES THE n-COMPONENT VECTOR 'BST': GIVEN THE 'elm'
   ELEMENT SPARSE MATRIX 'PP', SORTED BY INDEX 'idx', INITIALIZES 'BST'
   SUCH THAT FOR ALL j: BST[i]<=j<BST[i+1], PP->IDX[j][idx]=i */
void init_bst(int *BST,dSparse_Matrix *PP,int elm,int n,int idx)
{
  int i;

  for(i=1;i<n;i++) BST[i]=0;
  for(i=elm;i>0;i--) BST[PP->IDX[i][idx]]=i;
  BST[n]=elm+1;
  for(i=n-1;i>0;i--)
    if(BST[i]==0)
      BST[i]=BST[i+1];
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



/* "zero_d3tensor" ALLOCATES MEMORY FOR A DOUBLE
   3-TENSOR AND INITIALIZES IT TO ZERO */
double ***zero_d3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  static double ***T;
  int i,j,k;

  T=d3tensor(nrl,nrh,ncl,nch,ndl,ndh);
  for(i=nrl;i<=nrh;i++)
    for(j=ncl;j<=nch;j++)
      for(k=ndl;k<=ndh;k++)
	T[i][j][k]=0.0;
  return T;
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





