/* Author: Tim Lezon, Ahmet Bakan */

#include "Python.h"
#include "numpy/arrayobject.h"
//#include "membranmutil.h"
//#include "nrutil.h"
//#include "nr.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define NR_END 1

#define DEFGAM 1.0   /* Default force constant */
#define DEFCUT 11.0  /* Default cutoff distance in angstroms */
#define DEFMSCL 16.0  /* Default membrane spring scale factor */
#define DEFMLO -13.4 /* Default lower membrane boundary */
#define DEFMHI 13.4  /* Default upper membrane boundary */


#define FREE_ARG char*

#define PDB_MAX_LINE 90

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))





/* PDB file-related structures */
typedef struct {int num;char chain;} Resid;
typedef struct {char HEAD[7];int atmnum;char ATOM[5];char RES[4];char chain;int resnum;float X[3];float beta;char ELEMENT[3];int model;} Atom_Line;
typedef struct {char **HEADER;Atom_Line *atom;} PDB_File;

/* Rigid block-related structures */
typedef struct {char RES[4];char chain;int resnum;float disp;} Disp_Iso;
typedef struct {char RES[4];char chain;int resnum;float X[3];} Disp_Aniso;
typedef struct {int blknum;char LORES[4];char lochain;int lonum;char HIRES[4];char hichain;int hinum;} Rigid_Block;

typedef struct {int **IDX;double *X;} dSparse_Matrix;

int assign_rigid_blocks(PDB_File *,Rigid_Block *,int,int,int *);
int bless_from_tensor(double **,double ***,int **,int);
void copy_dsparse(dSparse_Matrix *,dSparse_Matrix *,int,int);
char **cmatrix(long nrl,long nrh,long ncl,long nch);
void cross(double [],double [],double []);
void dblock_hessian5(PDB_File *,dSparse_Matrix *,int,int,int,double);
void dsort_PP2(dSparse_Matrix *,int,int);
int find_contacts1(int **,PDB_File *,int,int,double);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
void free_PDB(PDB_File *PDB,int nhed,int nres);
char *get_param(char *file,char *param);
void hess_superrow1(double **,int **,PDB_File *,int,int,double,double);
void init_bst(int *,dSparse_Matrix *,int,int,int);
void pdb_hmr(char *file,int *nhed,int *nmod,int *nca);
void pdb_init(PDB_File *PDB,int nhed,int nres);
void print_hess_tensor(double ***,int **,int);
void print_prj_ofst(dSparse_Matrix *,int);
void read_blockfile(char *,char **,Rigid_Block **,int *);
int read_pdb3(char *file,PDB_File *PDB,int nhed,int nres);
void righthand2(double *,double **,int);
int **unit_imatrix(long,long);
double ***zero_d3tensor(long,long,long,long,long,long);
double **zero_dmatrix(long,long,long,long);

void read_command_line(int,char *[],int *,double *,double *,
		       double *,double *,int *,double *);
int calc_blessian_mem(PDB_File *,dSparse_Matrix *,int,int,int,double **,double);
void hess_superrow_mem(double **,int **,PDB_File *,int,int,double,double);
void bless_to_hess(double **,double **,dSparse_Matrix *,int,int);
double scalefunc0(double,double);
void backproject(double **,dSparse_Matrix *,int,int,int,int);
int dblock_projections2(dSparse_Matrix *,PDB_File *,int,int,int);

/* Numerical Recipes */
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void deigsrt(double d[], double **v, int n);
double dpythag(double a, double b);
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void deigsrt(double d[], double **v, int n);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void nrerror(char error_text[]);
void free_lvector(unsigned long *v, long nl, long nh);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
double *dvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
int *ivector(long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);

double MHI,MLO,CUT,MSCL;


static PyObject *buildhessian(PyObject *self, PyObject *args, PyObject *kwargs) {
  /* C-specific variables for RTB routines */
  //Rigid_Block *BLX;
  PDB_File PDB;
  dSparse_Matrix PP,HH;
  //char *pdbfile;
  double **HB;
  //double gam;
  int elm,bdim,i,j;
  //int nres,nrg,nh1,nm1,all,prm,nz=0;




  PyArrayObject *coords, *blocks, *hessian, *projection;
  double cutoff = 15., gamma = 1.;
  int natm, nblx, bmx;

  static char *kwlist[] = {"coords", "blocks", "hessian", "projection",
			   "natoms", "nblocks", "maxsize", "cutoff", 
			   "gamma", NULL};



  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOiii|dd", kwlist,
				   &coords, &blocks, &hessian, &projection,
				   &natm, &nblx, &bmx, &cutoff, &gamma))
    return NULL;

  double *XYZ = (double *) PyArray_DATA(coords);
  long *BLK = (long *) PyArray_DATA(blocks);
  double *hess = (double *) PyArray_DATA(hessian);
  double *proj = (double *) PyArray_DATA(projection);

  /*
  for(i=0;i<natm;i++)
    fprintf(stderr,"%d\t%d\n",i,BLK[i]);
  Py_RETURN_NONE;
  */


  /* ---------- This is where the body of the function goes. -------- */

  /* First allocate a PDB_File object to hold the coordinates and block 
     indices of the atoms.  This wastes a bit of memory, but it prevents 
     the need to re-write all of the RTB functions that are used in 
     standalone C code. */
  pdb_init(&PDB,0,natm);
  /*
  if(!PDB)
    return PyErr_NoMemory();
  */
  for(i=1;i<=natm;i++){
    PDB.atom[i].model=BLK[i-1];
    for(j=0;j<3;j++)
      PDB.atom[i].X[j]=XYZ[j*natm+i-1];
  }



  /*------------------------------------ */
  /* 
  read_command_line(argc,argv,&prm,&MLO,&MHI,&CUT,&MSCL,&all,&gam);
  read_blockfile(argv[1],&pdbfile,&BLX,&nrg);
  pdb_hmr(pdbfile,&nh1,&nm1,&nres);
  pdb_init(&PDB1,nh1,nres);
  read_pdb3(pdbfile,&PDB1,nh1,nres);
  if(all==0){
    nblx=assign_rigid_blocks(&PDB1,BLX,nres,nrg,&bmx);
  }
  else{
    nblx=nres;
    bmx=1;
    for(i=1;i<=nres;i++)
      PDB1.atom[i].model=i;
  }
  fprintf(stderr,"\n%s: %d block ranges specified\n",argv[1],nrg);
  fprintf(stderr,"%s: %d residues, %d asymmetric units\n",pdbfile,nres,nm1);
  fprintf(stderr,"%d total rigid blocks in structure\n",nblx);
  fprintf(stderr,"largest block contains %d residues\n",bmx);
  free(pdbfile);
  free(BLX);
  */


  /* Find the projection matrix */
  HH.IDX=imatrix(1,18*bmx*nblx,1,2);
  HH.X=dvector(1,18*bmx*nblx);
  elm=dblock_projections2(&HH,&PDB,natm,nblx,bmx);
  PP.IDX=imatrix(1,elm,1,2);
  PP.X=dvector(1,elm);
  for(i=1;i<=elm;i++){
    PP.IDX[i][1]=HH.IDX[i][1];
    PP.IDX[i][2]=HH.IDX[i][2];
    PP.X[i]=HH.X[i];
  }
  free_imatrix(HH.IDX,1,18*bmx*nblx,1,2);
  free_dvector(HH.X,1,18*bmx*nblx);
  dsort_PP2(&PP,elm,1);
  

  /* Calculate the block Hessian */
  HB=dmatrix(1,6*nblx,1,6*nblx);
  bdim=calc_blessian_mem(&PDB,&PP,natm,nblx,elm,HB,gamma);
  
  /*------------------------------------ */







  /* 
     Put the block Hessian and projection matrix into Python objects:
     HH[i][j] = hess[6*nblx*(i-1)+j-1]
     PP[i].X = proj[6*nblx*(PP[i].IDX[1]-1)+PP[i].IDX[2]-1
  */
  
  
  for(i=1;i<=bdim;i++)
    for(j=1;j<=bdim;j++)
      hess[bdim*(i-1)+j-1]=HB[i][j];
  for(i=1;i<=elm;i++)
    proj[bdim*(PP.IDX[i][1]-1) + PP.IDX[i][2]-1] = PP.X[i];
  


  free_PDB(&PDB,0,natm);
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


/* "read_command_line" checks for the input arguments. */
void read_command_line(int argc,char *argv[],int *prm,double *mlo,double *mhi,
		       double *cut,double *mscl,int *all,double *gam)
{
  char *param;
  double x;
  int i,ok=0,nopt=0;

  /* Assign default values */
  *prm=-1;
  *cut=DEFCUT;
  *mscl=DEFMSCL;
  *mlo=DEFMLO;
  *mhi=DEFMHI;
  *all=0;
  *gam=DEFGAM;
  

  /* Allowed flags:
     -p    Print the projection matrix
     -a    Use all residues instead of the specified blocks
     -s    Specify membrane scaling factor
     -c    Specify cutoff distance 
     -b    Specify membrane boundaries
     -k    Specify force constant
  */

  if(argc>1 && argc<=13){


    /* Assign values from blockfile, if any */
    param=(char *)malloc((size_t) (99*sizeof(char)));
    if((param=get_param(argv[1],"mlo"))!=NULL)
      sscanf(param,"%lf",mlo);
    if((param=get_param(argv[1],"mhi"))!=NULL)
      sscanf(param,"%lf",mhi);
    if((param=get_param(argv[1],"cut"))!=NULL)
      sscanf(param,"%lf",cut);
    if((param=get_param(argv[1],"mscl"))!=NULL)
      sscanf(param,"%lf",mscl);

    ok=1;
    i=2;
    while(i<argc){
      if(strcmp(argv[i],"-p")==0){       /* Print mode */
	*prm=1;
	nopt++;
      }
      else if(strcmp(argv[i],"-a")==0){  /* Forget blocks; use all residues */
	*all=1;
	nopt++;
      }
      else if(strcmp(argv[i],"-s")==0){  /* Membrane scaling factor */
	sscanf(argv[++i],"%lf",mscl);
	nopt++;
      }
      else if(strcmp(argv[i],"-c")==0){  /* Cutoff distance */
	sscanf(argv[++i],"%lf",cut);
	nopt++;
      }
      else if(strcmp(argv[i],"-b")==0){  /* Membrane boundaries */
	sscanf(argv[++i],"%lf",mlo);
	sscanf(argv[++i],"%lf",mhi);
	if(*mlo > *mhi){
	  x=*mhi;
	  *mhi=*mlo;
	  *mlo=x;
	}
	nopt++;
      }
      else if(strcmp(argv[i],"-k")==0){  /* Force constant */
	sscanf(argv[++i],"%lf",gam);
	nopt++;
      }
      else{
	fprintf(stderr,"\n%s: Unknown argument: %s\n\n",argv[0],argv[i]);
	exit(1);}
      i++;
    }
    if(nopt>6) ok=0;
  }
  if(ok==0){
    fprintf(stderr,"\nUsage:\n%s blockfile.blk [OPTIONS]\n\n",argv[0]);
    fprintf(stderr,"OPTIONS:\n%s%s%s%s%s%s",
	    "\t-p\tPrint the projection matrix instead of Hessian\n",
	    "\t-a\tDisregard blocks and use all residues (standard ANM)\n",
	    "\t-s MSCL\tSpecify z-scaling factor in membrane\n",
	    "\t-c CUT\tSpecify cut-off distance\n",
	    "\t-k GAM\tSpecify force constant\n",
	    "\t-b MLO MHI\tSpecify membrane boundaries\n\n");
    fprintf(stderr,
	    "Output:\nPrints Hessian (default) or RTB projection matrix\n\n");
    exit(1);
  }
  fprintf(stderr,"%s%d\n%s%d\n%s%f\n%s%f\n%s%f\n%s%f\n%s%d\n%s%f\n",
	  "Print: ",*prm,"All: ",*all,"mscl: ",*mscl,"mlo: ",*mlo,
	  "mhi: ",*mhi,"cut: ",*cut,"nopt: ",nopt,"gam: ",*gam);
  return;
}
  

/* "calc_blessian_mem" is the membrane version of dblock_hessian5: 
   in it, 'hess_superrow1' is replaced with 'hess_superrow_mem'. */
int calc_blessian_mem(PDB_File *PDB,dSparse_Matrix *PP1,
		      int nres,int nblx,int elm,double **HB,double gam)
{
  dSparse_Matrix *PP2;
  double **HR,***HT;
  int **CT,*BST1,*BST2;
  int ii,i,j,k,p,q,q1,q2,ti,tj,bi,bj,sb,nc,out;

  /* TESTING */
  double dd;
  int ok,jj;



  /* ------------------- INITIALIZE LOCAL VARIABLES ------------------- */
  fprintf(stderr,"calc_blessian_mem: initializing local variables...\n");

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
  fprintf(stderr,"calc_blessian_mem: finding block contacts...\n");
  nc=find_contacts1(CT,PDB,nres,nblx,CUT);


  /* Allocate a tensor for the block Hessian */
  HT=zero_d3tensor(1,nc,1,6,1,6);


  /* Calculate each super-row of the full Hessian */
  fprintf(stderr,"calc_blessian_mem: calculating full hessian...\n");
  for(ii=1;ii<=nres;ii++){

    if(PDB->atom[ii].model!=0){

      /* ----------------- FIND SUPER-ROW OF FULL HESSIAN --------------- */
      hess_superrow_mem(HR,CT,PDB,nres,ii,CUT,gam);
	  

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
  fprintf(stderr,"calc_blessian_mem: projecting into block space...\n");
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


/* "hess_superrow_mem" calculates the 'who'-th super-row 
   of the Hessian, using 'cut' as the cutoff and 'gam' as the 
   spring constant for all interactions. */
void hess_superrow_mem(double **HR,int **CT,PDB_File *PDB,int nres,
		    int who,double cut,double gam)
{
  int i,j,k,jj;
  double DX[3],csq=cut*cut,dsq,df,dd,scl=1.0;

  /* Clear the diagonal super-element */
  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
      HR[3*(who-1)+i][j]=0.0;

  /* Calculate the submatrices */
  for(jj=1;jj<=nres;jj++){

    /* scalefunc0: Uniform z-scaling inside membrane */
    scl=scalefunc0(PDB->atom[who].X[2],PDB->atom[jj].X[2]);


    if(jj!=who && PDB->atom[jj].model!=0 && 
       CT[PDB->atom[who].model][PDB->atom[jj].model]!=0){
      dsq=0.0;
      for(k=0;k<3;k++){
	DX[k] = (double)PDB->atom[who].X[k] - PDB->atom[jj].X[k];
	dsq+=(DX[k]*DX[k]);
      }


      if(dsq<csq){
	for(i=1;i<=3;i++){
	  for(j=i;j<=3;j++){

	    df=gam*DX[i-1]*DX[j-1]/dsq;


	    /* Strong backbone bonds */
	    if((int)fabs(PDB->atom[who].resnum-PDB->atom[jj].resnum)==1 && 
	       PDB->atom[who].chain==PDB->atom[jj].chain)
	      df*=100.0;


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


    /* **** THIS IS PROBABLY NOT NECESSARY, AS THESE ELEMENTS ARE LIKELY
       NOT USED IN PROJECTING INTO THE BLOCK SPACE **** */
    else if(jj!=who && PDB->atom[jj].model!=0)
      for(i=1;i<=3;i++)
	for(j=1;j<=3;j++)
	  HR[3*(jj-1)+i][j]=HR[3*(jj-1)+j][i]=0.0;
  }
}


/* "bless_to_hess" CONVERTS A BLOCK HESSIAN BACK INTO 
   A FULL HESSIAN USING THE PROJECTION MATRIX PP */
void bless_to_hess(double **HB,double **HF,dSparse_Matrix *PP,int nres,int elm)
{
  int *I1,*I2,ii,jj,i,j,a,b,max=0;

  /* Get a list of indices */
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


  /* Zero the matrix */
  for(i=1;i<=3*nres;i++)
    for(j=i;j<=3*nres;j++)
      HF[i][j]=HF[j][i]=0.0;


  /* ----- Calculate full Hessian ----- */
  for(ii=1;ii<=elm;ii++)
    for(jj=1;jj<=elm;jj++){

      i=PP->IDX[ii][1];
      a=I2[PP->IDX[ii][2]];
      j=PP->IDX[jj][1];
      b=I2[PP->IDX[jj][2]];
      HF[i][j]+=(PP->X[ii]*PP->X[jj]*HB[a][b]);
      HF[j][i]=HF[i][j];
    }
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}


/* "scalefunc0" returns the factor by which the x- and y- rows of 
   the Hessian superelement are to be scaled.  It returns 1 if neither
   residue is in the membrane, sqrt(MSCL) if both residues are in the 
   membrane, and MSCL^(1/4) if only one residue is in the membrane.
 */
double scalefunc0(double z1,double z2)
{
  double scl=1.0;
  static double s0;
  static int first=1;

  /* CALCULATE THE CONSTANTS ON FIRST CALL */
  if(first==1){
    first=0;
    s0=sqrt(MSCL);
    s0=sqrt(s0);
  }

  if(z1<MHI && z1>MLO) scl*=s0;
  if(z2<MHI && z2>MLO) scl*=s0;
  
  return scl;
}


/* "backproject" projects block modes back into the full-residue space.  
   'VEC' is a 'dim' by 'nev' matrix that holds the block vectors in its 
   first 'imx' rows.  The matrix 'PP' has 'elm' elements and projects 
   between block and full spaces. */ 
void backproject(double **VEC,dSparse_Matrix *PP,
		 int dim,int nev,int imx,int elm)
{
  double **VT,dd;
  int *I1,*I2,max,i,j,k,q;
  int imax=0,kmax=0;


  if(imx>dim){
    fprintf(stderr,
	    "\nbackproject: block dimension %d exceeds full dimension %d\n\n",
	    imx,dim);
    exit(1);}


  /* I1[i] = 0 if column 'i' contains only zeros; otherwise I1[i] = i
     I2[i] = rank of 'i' among non-zero columns */
  max=0;
  for(i=1;i<=elm;i++)
    if(PP->IDX[i][2]>max)
      max=PP->IDX[i][2];
  I1=ivector(1,max);
  I2=ivector(1,max);
  for(i=1;i<=max;i++) I1[i]=0;
  for(i=1;i<=elm;i++) I1[PP->IDX[i][2]]=PP->IDX[i][2];
  j=0;
  for(i=1;i<=max;i++){
    if(I1[i]!=0) j++;
    I2[i]=j;
  }


  /* Copy the first 'imx' rows of 'VEC' to a temporary array 
     and clear all of the elements of 'VEC' */
  VT=dmatrix(1,imx,1,nev);
  for(j=1;j<=nev;j++){
    for(i=1;i<=imx;i++){
      VT[i][j]=VEC[i][j];
      VEC[i][j]=0.0;
    }
    for(i=imx+1;i<=dim;i++) VEC[i][j]=0.0;
  }
  

  /* Calculate: VEC_{ij} = \sum_k PP_{ik}*VT_{kj} */
  for(q=1;q<=elm;q++){
    i=PP->IDX[q][1];
    k=I2[PP->IDX[q][2]];
    if(i>imax) imax=i;
    if(k>kmax) kmax=k;
    dd=PP->X[q];
    for(j=1;j<=nev;j++) VEC[i][j]+=dd*VT[k][j];
  }

  free_dmatrix(VT,1,imx,1,nev);
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}



/* "dblock_projections2" CALCULATES THE PROJECTION 
   FROM FULL RESIDUE SPACE TO RIGID BLOCK SPACE */
int dblock_projections2(dSparse_Matrix *PP,PDB_File *PDB,
			int nres,int nblx,int bmx)
{
  double **X,**I,**IC,*CM,*W,**A,**ISQT;
  double alpha,beta,gamma,x,tr,dd,df;
  int *IDX,nbp,b,i,j,k,p,ii,jj,kk,pp,aa,bb,cc,elm;


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
    /* must this be a right-handed coordinate system? */

    /* FIND ITS SQUARE ROOT */
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++){
	dd=0.0;
	for(k=1;k<=3;k++)
	  dd+=A[i][k]*A[j][k]/sqrt(W[k]);
	ISQT[i][j]=dd;
      }

    /* UPDATE PP WITH THE RIGID MOTIONS OF THE BLOCK */
    tr=1.0/sqrt((float)nbp);
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

/* "assign_rigid_blocks" ASSIGNS EACH RESIDUE IN THE PDB STRUCTURE 
   TO A RIGID BLOCK AND RETURNS THE TOTAL NUMBER OF RIGID BLOCKS */
int assign_rigid_blocks(PDB_File *PDB,Rigid_Block *BLX,
			int nres,int nrg,int *bmx)
{
  char ch;
  int *UB,**BMP,nb,nbt,ok,kwn,mold,rr,i,j,k;

  /* DETERIMINE THE NUMBER OF UNIQUE BLOCKS PER CHAIN */
  UB=ivector(1,nrg);
  for(i=1;i<=nrg;i++) UB[i]=0;
  nb=0;
  for(i=1;i<=nrg;i++){
    ok=1;
    for(j=1;j<=nb;j++)
      if(BLX[i].blknum==UB[j]){
	ok=0;
	break;
      }
    if(ok==1)
      UB[++nb]=BLX[i].blknum;
  }


  /* ASSIGN EACH RESIDUE TO A BLOCK: BMP[i][1] CONTAINS THE NUMBER OF 
     BLOCK i AS PROVIDED BY THE .blk FILE.  BMP[i][2] CONTAINS THE 
     NUMBER OF BLOCK i THAT WILL BE USED IN THE RTB CALCULATION. */
  BMP=imatrix(1,nb,1,2);
  nbt=kwn=0;
  for(i=1;i<=nb;i++)
    BMP[i][1]=BMP[i][2]=0;
  mold=PDB->atom[1].model;

  for(i=1;i<=nres;i++){
    rr=PDB->atom[i].resnum;
    ch=PDB->atom[i].chain;
    if(PDB->atom[i].model!=mold){
      for(j=1;j<=nb;j++)
	BMP[j][1]=BMP[j][2]=0;
      mold=PDB->atom[i].model;
      kwn=0;
    }
    ok=0;
    for(j=1;j<=nrg;j++){
      if((BLX[j].lochain<ch && ch<BLX[j].hichain) ||
	 (BLX[j].lochain==ch && ch<BLX[j].hichain && rr>=BLX[j].lonum) ||
	 (BLX[j].lochain<ch && ch==BLX[j].hichain && rr<=BLX[j].hinum) ||
	 (BLX[j].lochain==ch && BLX[j].hichain==ch && 
	  rr>=BLX[j].lonum && rr<=BLX[j].hinum)){
	for(k=1;k<=kwn;k++)
	  if(BLX[j].blknum==BMP[k][1]){
	    ok=1;
	    PDB->atom[i].model=BMP[k][2];
	    break;
	  }
	if(ok==0){
	  BMP[++kwn][1]=BLX[j].blknum;
	  BMP[kwn][2]=PDB->atom[i].model=++nbt;
	  ok=1;
	}
	break;
      }
    }
    if(ok==0)
      PDB->atom[i].model=0;
  }
  free_ivector(UB,1,nrg);
  free_imatrix(BMP,1,nb,1,nb);


  /* FIND THE SIZE OF THE LARGEST BLOCK */
  UB=ivector(1,nbt);
  for(i=1;i<=nbt;i++) UB[i]=0;
  for(i=1;i<=nres;i++) 
    if(PDB->atom[i].model!=0)
      UB[PDB->atom[i].model]++;
  (*bmx)=0;
  for(i=1;i<=nbt;i++){
    if(UB[i]>(*bmx))
      (*bmx)=UB[i];
  }
  free_ivector(UB,1,nbt);

  return nbt;
}



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


/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
char **cmatrix(long nrl,long nrh,long ncl,long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  static char **m;

  /* allocate pointers to rows */
  m=(char **) malloc((size_t)((nrow+1)*sizeof(char*)));
  if(!m){
    fprintf(stderr,"\nallocation failure 1 in cmatrix\n\n");
    exit(1);}
  m++;
  m-=nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(char *) malloc((size_t)((nrow*ncol+1)*sizeof(char)));
  if(!m[nrl]){
    fprintf(stderr,"\nallocation failure 2 in cmatrix\n\n");
    exit(1);}
  m[nrl]++;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
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


/* "dsort_PP2" SORTS THE PROJECTION MATRIX IN ASCENDING ORDER OF THE  
   INDEX 'idx'.  ADAPTED FROM THE NUMERICAL RECIPES 'HEAPSORT' ROUTINE. */
void dsort_PP2(dSparse_Matrix *MM,int n,int idx)
{
  double x;
  int i,ir,j,l,hi,i1,i2,ndx;
  long rra,*ra;

  if(n<2) return;
  if(idx<1 || idx>2){
    fprintf(stderr,"dsort_PP2: bad index value %d\n\n",idx);
    exit(1);}
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


  /* CHECK FOR PROPER SORTING */
  for(i=1;i<n;i++){
    if(MM->IDX[i+1][idx]<MM->IDX[i][idx] || 
       (MM->IDX[i+1][idx]==MM->IDX[i][idx] && 
	MM->IDX[i+1][ndx]<MM->IDX[i][ndx])){
      fprintf(stderr,"\n\nImproper sort in dsort_PP2:\n");
      fprintf(stderr,"%d:\t%d\t%d\n",i,MM->IDX[i][idx],MM->IDX[i][ndx]);
      fprintf(stderr,"%d:\t%d\t%d\n\n",
	      i+1,MM->IDX[i+1][idx],MM->IDX[i+1][ndx]);
      exit(1);
    }
  }
}



/* "find_contacts" FINDS WHICH BLOCKS ARE IN CONTACT, AND ASSIGNS EACH 
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


/* free a char matrix allocated by matrix() */
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-1));
	free((FREE_ARG) (m+nrl-1));
}


/* "free_PDB" FREES MEMORY FROM A PDB STRUCTURE */
void free_PDB(PDB_File *PDB,int nhed,int nres)
{
  int i;

  free_cmatrix(PDB->HEADER,1,nhed,0,PDB_MAX_LINE-1);
  free(PDB->atom);
}


/* "get_param" RETURNS THE PARAMETER OF THE 
   SPECIFIED NAME IN THE SPECIFIED FILE */
char *get_param(char *file,char *param)
{
  FILE *data;
  char LINE[200],*s;
  static char *garp;

  garp=(char *)calloc(99,sizeof(char));
  strcpy(garp,param);
  strcat(garp,"=");
  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nget_param: unable to open %s\n\n",file);
    exit(1);}
  while(!feof(data)){
    fgets(LINE,200,data);
    if(strstr(LINE,garp)!=NULL && LINE[0]!='#'){
      s=strpbrk(LINE,"=")+1;
      strcpy(garp,s);
      garp[strlen(garp)-1]='\0';
      return garp;
    }
  }
  return NULL;
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


/* "pdb_hmr" ESTIMATES THE NUMBER OF HEADER 
   LINES, MODELS AND RESIDUES IN A PDB FILE */
void pdb_hmr(char *file,int *nhed,int *nmod,int *nca)
{
  FILE *data;
  char HED[7],ATM[5],calt,c;
  int res,i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\npdb_hmr: unable to open %s\n\n",file);
    exit(1);}
  (*nca)=(*nhed)=(*nmod)=0;
  for(;;){

    /* Cols 1-6 are field name */
    for(i=0;i<6;i++)
      HED[i]=getc(data);

    /* Check for ATOM lines */
    if(!strncmp(HED,"ATOM  ",6)){
      for(i=7;i<=12;i++)            // Skip atom numbers
	c=getc(data);
      for(i=13;i<=16;i++)           // Find atom name
	ATM[i-13]=getc(data);
      calt=getc(data);              // Find alternative location
      if(strstr(ATM,"CA")!=NULL && (calt=='A' || calt==' '))
	(*nca)++;
    }

    /* Check for lines to be included in the header */
    else if(!strncmp(HED,"HEADER",6) || !strncmp(HED,"TITLE ",6) || 
	    !strncmp(HED,"COMPND",6) || !strncmp(HED,"SEQRES",6) || 
	    !strncmp(HED,"HELIX ",6) || !strncmp(HED,"SHEET ",6) || 
	    !strncmp(HED,"TURN  ",6))
      (*nhed)++;

    /* Check for rotation matrices */
    else if(!strncmp(HED,"MODEL ",6))
    //else if(!strncmp(HED,"MTRIX1",6))
      (*nmod)++;

    else if(!strncmp(HED,"END   ",6) || feof(data))  // Identify END of file
      break;
    do{
      c=getc(data);
    }while(c!='\n');
  }
  fclose(data);
  if(*nmod==0)
    *nmod=1;
}



/* "pdb_init" INITIALIZES A PDB_File STRUCTURE: 'HEADER' IS A CHARACTER
   ARRAY WITH ROWS (1,nhed) AND COLUMNS (0,PDB_MAX_LINE-1); Calpha IS
   AN ARRAY OF Atom_Line STRUCTURES WITH ELEMENTS (1,nres).  */
void pdb_init(PDB_File *PDB,int nhed,int nres)
{

  /* Allocate memory for PDB.HEADER */
  PDB->HEADER=cmatrix(1,nhed,0,PDB_MAX_LINE-1);

  /* Allocate memory for PDB.atom */
  PDB->atom=malloc((size_t)((nres+2)*sizeof(Atom_Line)));
  if(!PDB->atom){
    fprintf(stderr,"\npdb_init: fail to allocate atom\n\n");
    exit(1);}
}



/* "print_prj_ofst" PRINTS THE PROJECTION MATRIX, OFFSETTING THE COLUMN 
   ELEMENTS IF THERE ARE BLOCKS WITH FEWER THAN SIX DEGREES OF FREEDOM */
void print_prj_ofst(dSparse_Matrix *PP,int elm)
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
  fprintf(stderr,"Printing projection matrix...\n");
  for(i=1;i<=elm;i++)
    if(PP->X[i]!=0.0)
      printf("%8d%8d% 25.15e\n",PP->IDX[i][1],I2[PP->IDX[i][2]],PP->X[i]);
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}


/* "read_blockfile" READS THE INPUT blk FILE FOR 
   INFORMATION ON THE SOURCE PDB AND THE BLOCK STRUCTURE */
void read_blockfile(char *file,char **pdbfile,Rigid_Block **BLK,int *nrg)
{
  FILE *data;
  Rigid_Block *Btmp;
  char *ptmp,HED[7],BNK[8],c;
  int blk=0,rng,nn,i;
  long int blkst;


  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_blockfile: unable to open %s\n\n",file);
    exit(1);}
  for(;;){

    /* GET THE FIRST CHARACTER OF EACH LINE */
    c=getc(data);


    /* COMMENT OR WHITESPACE: SKIP OVER REST OF LINE */
    if(c=='#' || isspace(c))
      while(c!='\n' && c!=EOF) 
	c=getc(data);


    /* OTHERWISE CHECK FOR READ ROW HEADING */
    else{
      i=0;
      while(!isspace(c) && i<6){
	HED[i++]=c;
	c=getc(data);
      }
      HED[i]='\0';


      /* PDB FILENAME */
      if(!strncmp(HED,"PDB",3)){
	i=0;
	do{
	  c=getc(data);
	  i++;
	}while(c!='\n' && c!=EOF);

	/* Allocate space for filename, then go back and read it */
	ptmp=(char *)malloc((size_t)(i+1)*sizeof(char));
	fseek(data,-i,SEEK_CUR);
	i=0;
	do{
	  c=getc(data);
	  if(!isspace(c))
	    ptmp[i++]=c;
	}while(c!='\n' && c!=EOF);
	ptmp[i]='\0';
	*pdbfile=ptmp;
      }


      /* BLOCK DATA */
      else if(!strncmp(HED,"BLOCK",5)  && blk==0){
	blk=1;
	nn=0;
	blkst=ftell(data)-strlen(HED)-1;

	/* Count the number of ranges */
	while(!strncmp(HED,"BLOCK",5)){
	  while(c!='\n' && c!=EOF)
	    c=getc(data);
	  c=getc(data);
	  if(c!='#' && !isspace(c)){
	    i=0;
	    while(c!='\n' && c!=EOF && i<6){
	      HED[i++]=c;
	      c=getc(data);
	    }
	    HED[i]='\0';
	    nn++;
	  }
	}


	/* Allocate space for ranges, then go back and read the data */
	Btmp=(Rigid_Block *)malloc((size_t)(nn+1)*sizeof(Rigid_Block));
	fseek(data,blkst,SEEK_SET);


	/* --------------- READ ALL THE BLOCKS --------------- */
	rng=1;
	do{
	  c=getc(data);
	  if(c!='#' && !isspace(c) && c!=EOF){
	    i=0;
	    while(!isspace(c) && i<6){
	      HED[i++]=c;
	      c=getc(data);
	    }
	    HED[i]='\0';
	    if(!strncmp(HED,"BLOCK",5)){
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      sscanf(BNK,"%d",&Btmp[rng].blknum);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      strncpy(Btmp[rng].LORES,BNK,3);
	      while(isspace(c)) c=getc(data);
	      Btmp[rng].lochain=c;
	      c=getc(data);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      sscanf(BNK,"%d",&Btmp[rng].lonum);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      strncpy(Btmp[rng].HIRES,BNK,3);
	      while(isspace(c)) c=getc(data);
	      Btmp[rng].hichain=c;
	      c=getc(data);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      sscanf(BNK,"%d",&Btmp[rng].hinum);
	      while(c!='\n' && c!=EOF) c=getc(data);
	      rng++;
	    }
	  }
	  else{
	    while(c!='\n' && c!=EOF) c=getc(data);
	    if(c==EOF){
	      fclose(data);
	      *BLK=Btmp;
	      *nrg=nn;
	      return;
	    }
	  }
	}while(!strncmp(HED,"BLOCK",5));
      }

      /* END OF FILE */
      if(!strncmp(HED,"END",3)){
	fclose(data);
	*BLK=Btmp;
	*nrg=nn;
	return;
      }

      /* EAT UP REST OF LINE */
      while(c!='\n' && c!=EOF) c=getc(data);
      if(c==EOF){
	fclose(data);
	*BLK=Btmp;
	*nrg=nn;
	return;
      }
    } // <------- END of else{...
  } // <--------- END of for(;;){...
  fclose(data);
  free(ptmp);
}


/* "read_pdb3" READS HEADER AND Calpha INFORMATION FROM A PDB FILE */
int read_pdb3(char *file,PDB_File *PDB,int nhed,int nres)
{
  FILE *data;
  char LINE[PDB_MAX_LINE],HED[7],ATM[5],HOLD[9],c,calt;
  int ca,hd,mdl,i,j;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_pdb3: can't open %s\n\n",file);
    exit(1);}

  ca=hd=mdl=1;
  for(;;){
    
    i=0;
    do{                                  /* Skip remaining columns */
      c=getc(data);
      LINE[i++]=c;
    }while(c!='\n' && !feof(data));
    LINE[i]='\0';
    for(i=0;i<=5;i++)                       /* Cols 1-6 are field name */
      HED[i]=LINE[i];
    HED[6]='\0';



    /* ---------- Check for ATOM lines ---------- */
    if(!strncmp(HED,"ATOM  ",6)){// || !strncmp(HED,"HETATM",6)){

      //fprintf(stderr,"ca = %d\n",ca);

      for(i=0;i<=5;i++) PDB->atom[ca].HEAD[i]=LINE[i];
      PDB->atom[ca].HEAD[6]='\0';

      for(i=12;i<=15;i++)                     /* Cols 13-16 are atom name*/
	ATM[i-12]=LINE[i];
      ATM[4]='\0';
      calt=LINE[16];                          /* Col 17 is alt. location */

      /* Keep atom if it is alpha carbon */
      if(strstr(ATM,"CA")!=NULL && (calt=='A' || calt==' ')){

	for(i=6;i<=10;i++)                    /* Cols 7-11 are atom # */
	  HOLD[i-6]=LINE[i];
	HOLD[5]='\0';
	sscanf(HOLD,"%d",&PDB->atom[ca].atmnum);

	for(i=12;i<=15;i++)                     /* Cols 13-16 are atom name*/
	  PDB->atom[ca].ATOM[i-12]=LINE[i];
	PDB->atom[ca].ATOM[4]='\0';

	for(i=17;i<=19;i++)                   /* Cols 18-20 are res. name */
	  PDB->atom[ca].RES[i-17]=LINE[i];
	PDB->atom[ca].RES[3]='\0';
	PDB->atom[ca].chain=LINE[21];       /* Col 22 is chain name */

	/* Assign a model number */
	PDB->atom[ca].model=mdl;

	for(i=22;i<=25;i++)                   /* Cols 23-26 are residue # */
	  HOLD[i-22]=LINE[i];
	HOLD[4]='\0';
	sscanf(HOLD,"%d",&PDB->atom[ca].resnum);
	for(i=30;i<=37;i++)                   /* Cols 31-38 are x-coordinate */
	  HOLD[i-30]=LINE[i];
	HOLD[8]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[0]);
	for(i=38;i<=45;i++)                   /* Cols 39-46 are y-coordinate */
	  HOLD[i-38]=LINE[i];
	HOLD[8]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[1]);
	for(i=46;i<=53;i++)                   /* Cols 47-54 are z-coordinate */
	  HOLD[i-46]=LINE[i];
	HOLD[8]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[2]);
	for(i=60;i<=65;i++)              /* Columns 61-66 are beta */
	  HOLD[i-60]=LINE[i];
	HOLD[6]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].beta);
	for(i=76;i<=77;i++)              /* Columns 77-78 are element */
	  PDB->atom[ca].ELEMENT[i-76] = LINE[i]=='\n' ? '\0' : LINE[i];
	PDB->atom[ca].ELEMENT[2]='\0';
	if(++ca>nres){
	  fclose(data);
	  return;
	}
      }
    }/* <---- End of 'if(!strncmp(HED,"ATOM  ",6)){... */


    /* ---------- Check for lines to be included in the header ---------- */
    else if(!strncmp(HED,"HEADER",6) || !strncmp(HED,"TITLE ",6) || 
	    !strncmp(HED,"COMPND",6) || !strncmp(HED,"SEQRES",6) ||  
	    !strncmp(HED,"HELIX ",6) || !strncmp(HED,"SHEET ",6) || 
	    !strncmp(HED,"TURN  ",6)){
      sprintf(PDB->HEADER[hd],"%s",LINE);
      hd++;
    }

    /* ---------- Check for MODEL lines ----------- */
    else if(!strncmp(HED,"ENDMDL",6)) 
      mdl++;

    /* ---------- Check for end of file ---------- */
    else if(!strncmp(HED,"END",3) || feof(data)){
      fclose(data);
      return;
    }
  }
}



/* "righthand2" MAKES SURE THAT THE EIGENVECTORS 
   FORM A RIGHT-HANDED COORDINATE SYSTEM */
void righthand2(double *VAL,double **VEC,int n)
{
  double A[3],B[3],C[3],CP[3],dot=0.0;
  int rev=0,i,j;

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

void dsvdcmp(double **a, int m, int n, double w[], double **v)
{
  //double dpythag(double a, double b);
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





