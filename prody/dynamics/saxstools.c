#include <Python.h>
#include <math.h>
#include "saxsConstants.h"
//#include "saxs.h"
#ifdef NPY_1_7_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NPE_PY_ARRAY_OBJECT PyArrayObject
#else
//TODO Remove this as soon as support for Numpy version before 1.7 is dropped
#define NPE_PY_ARRAY_OBJECT PyObject
#endif
#include <numpy/arrayobject.h>


//##########################FUNCTION PROTOTYPES##################################
//#################Functions to calculate structure factors######################
float DS(float x);
float RS(float x);
float ADE(float x);
float CYT(float x);
float GUA(float x);
float THY(float x);
float URA(float x);
float SOL(float x);
float GLY(float x);
float ALA(float x);
float VAL(float x);
float LEU(float x);
float ILE(float x);
float MET(float x);
float PHE(float x);
float TRP(float x);
float PRO(float x);
float SER(float x);
float THR(float x);
float CYS(float x);
float TYR(float x);
float ASN(float x);
float GLN(float x);
float ASP(float x);
float GLU(float x);
float LYS(float x);
float ARG(float x);
float HIS(float x);
float get_f(const char resname[],    float w, float q);
float get_f_numeric(int resname_num, float w, float q);
float min_dist_to_mol_v2(double *X, double *Y, double *Z, \
			 double  x, double  y, double  z, \
			 int atom_count, int *mol_type, int *type);

/* float min_dist_to_mol_v3(double *X, double *Y, double *Z, \ */
/* 			 double x, double y, double z,	  \ */
/* 			 int atom_count, int *mol_type); */

int cgSolvate(FILE *fpdb, double *X, double *Y, double *Z, \
	      int *total_atom, char *cgatom[], double *W,  \
	      float delta, float wDNA, float wRNA, float wPROT, int pdb_flag, \
	      float thickness, float closest_dist, int solvent_flag);
int cgSolvateNumeric_v2(char *fpdb_file, double *X, double *Y, double *Z, double *W, \
			float wDNA, float wRNA, float wPROT,		\
			float thickness, float closest_dist,		\
			int pdb_flag, int solvent_flag);

void calcSAXS(float *I, float *X,float *Y, float *Z, \
	      int total_atom, const char *cgatom[], float *W, float *Q, int nexp);
void calcSAXSNumeric(double *I, double *X,double *Y, double *Z, int total_atom, \
		     int *cgatom_num, double *W, double *Q, int nexp);
//###############################################################################



//###########################PYTHON BINDING SECTION##############################
static char module_docstring[] =
  "This module provides an interface for calculating SAXS profiles of proteins or nucleic acids in C.";
static char calcSAXSNumeric_docstring[] =
  "Calculate SAXS profile of a given protein or nucleic acid.";
static char cgSolvateNumeric_docstring[] =
  "Build a coarse-grained solvation shell around a biological macromolecule.";

static PyObject *saxstools_calcSAXSNumeric(PyObject *self, PyObject *args);
static PyObject *saxstools_cgSolvateNumeric(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"calcSAXSNumeric", saxstools_calcSAXSNumeric, METH_VARARGS, calcSAXSNumeric_docstring},
  {"cgSolvateNumeric", saxstools_cgSolvateNumeric, METH_VARARGS, cgSolvateNumeric_docstring},
  {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef saxstools = {
        PyModuleDef_HEAD_INIT,
        "saxstools",
        "SAXS tools.",
        -1,
        module_methods,
};
PyMODINIT_FUNC PyInit_saxstools(void) {
    import_array();
    return PyModule_Create(&saxstools);
}
#else
PyMODINIT_FUNC initsaxstools(void) {

    Py_InitModule3("saxstools", module_methods,
        module_docstring);
    import_array();
}
#endif

/*
PyMODINIT_FUNC initsaxstools(void)
{
  PyObject *m = Py_InitModule3("saxstools", module_methods, module_docstring);
  // PyObject *m = PyModule_Create("saxstools", module_methods, module_docstring);
  if (m == NULL)
    return 1;

  import_array();
}
*/

static PyObject *saxstools_calcSAXSNumeric(PyObject *self, PyObject *args)
{
  int total_atom, nexp;
  //  const char *cgatom[total_atom];
  PyObject *I_obj, *X_obj, *Y_obj, *Z_obj, *cgatom_num_obj, *W_obj, *Q_obj;
  
  /* Parse the input tuple */
  if (!PyArg_ParseTuple(args, "OOOOiOOOi", &I_obj, &X_obj, &Y_obj, &Z_obj, &total_atom, &cgatom_num_obj, &W_obj, &Q_obj, &nexp))
    return NULL;

  /* Interpret the input objects as numpy arrays. */
  PyObject *I_array = PyArray_FROM_OTF(I_obj, NPY_DOUBLE, NPY_ARRAY_INOUT_ARRAY);
  //  PyObject *I_array = PyArray_FROM_OTF(I_obj, NPY_FLOAT32, NPY_ARRAY_INOUT_ARRAY);
  PyObject *X_array = PyArray_FROM_OTF(X_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  //  PyObject *X_array = PyArray_FROM_OTF(X_obj, NPY_FLOAT32, NPY_ARRAY_IN_ARRAY);
  PyObject *Y_array = PyArray_FROM_OTF(Y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  PyObject *Z_array = PyArray_FROM_OTF(Z_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  PyObject *cgatom_num_array = PyArray_FROM_OTF(cgatom_num_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);
  PyObject *W_array = PyArray_FROM_OTF(W_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  PyObject *Q_array = PyArray_FROM_OTF(Q_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  /* If that didn't work, throw an exception. */
  if (I_array == NULL || X_array == NULL || Y_array == NULL || Z_array==NULL || W_array==NULL || Q_array==NULL || cgatom_num_array==NULL) 
    {
      Py_XDECREF(I_array);
      Py_XDECREF(X_array);
      Py_XDECREF(Y_array);
      Py_XDECREF(Z_array);
      Py_XDECREF(cgatom_num_array);
      Py_XDECREF(W_array);
      Py_XDECREF(Q_array);
      return NULL;
    }

  nexp = (int)PyArray_DIM(I_array, 0);
  
  /* Get pointers to the data as C-types. */
  double *I = (double*)PyArray_DATA(I_array);
  //  float *I = (float*)PyArray_DATA(I_array);
  double *X = (double*)PyArray_DATA(X_array);
  //  float *X = (float*)PyArray_DATA(X_array);
  double *Y = (double*)PyArray_DATA(Y_array);
  double *Z = (double*)PyArray_DATA(Z_array);
  int *cgatom_num = (int*)PyArray_DATA(cgatom_num_array);
  double *W = (double*)PyArray_DATA(W_array);
  double *Q = (double*)PyArray_DATA(Q_array);

  /* Call the external C function to compute the SAXS profile of a given
     solvated biological macromolecular coordinates.*/
  calcSAXSNumeric(I, X, Y, Z, total_atom, cgatom_num, W, Q, nexp);

  /* Clean up. */
  Py_DECREF(I_array);
  Py_DECREF(X_array);
  Py_DECREF(Y_array);
  Py_DECREF(Z_array);
  Py_DECREF(cgatom_num_array);
  Py_DECREF(W_array);
  Py_DECREF(Q_array);

  //  if (value < 0.0) {
  //    PyErr_SetString(PyExc_RuntimeError,
  //                    "Chi-squared returned an impossible value.");
  //    return NULL;
  //  }

  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("");
  return ret;
}

static PyObject *saxstools_cgSolvateNumeric(PyObject *self, PyObject *args)
{
  int total_atom=0, pdb_flag=0, solvent_flag=1;
  //int MAX_ATOM=50000;
  //  const char *cgatom[total_atom];
  float wDNA=0.07f;
  float wRNA=0.126f;
  float wPROT=0.04f;
  float thickness=3.0f;
  float closest_dist=3.5f;

  
  char *fpdb_file;
  
  PyObject *X_obj, *Y_obj, *Z_obj, *W_obj;

  /* Parse the input tuple */
  if (!PyArg_ParseTuple(args, "sOOOOfffffiii", &fpdb_file, &X_obj, &Y_obj, &Z_obj, &W_obj, &wDNA, &wRNA, &wPROT, &thickness, &closest_dist, &pdb_flag, &solvent_flag, &MAX_ATOM))
    return NULL;

  /* Interpret the input objects as numpy arrays. */
  PyObject *X_array = PyArray_FROM_OTF(X_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  PyObject *Y_array = PyArray_FROM_OTF(Y_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  PyObject *Z_array = PyArray_FROM_OTF(Z_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  /* PyObject *cgatom_num_array = PyArray_FROM_OTF(cgatom_num_obj, NPY_INT64, NPY_ARRAY_IN_ARRAY); */
  PyObject *W_array = PyArray_FROM_OTF(W_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (X_array == NULL || Y_array == NULL || Z_array==NULL || W_array==NULL)
    {
      Py_XDECREF(X_array);
      Py_XDECREF(Y_array);
      Py_XDECREF(Z_array);
      /* Py_XDECREF(cgatom_num_array); */
      Py_XDECREF(W_array);
      return NULL;
    }
  //  printf("Here I am \n");  
  /* Get pointers to the data as C-types. */
  double *X = (double*)PyArray_DATA(X_array);
  double *Y = (double*)PyArray_DATA(Y_array);
  double *Z = (double*)PyArray_DATA(Z_array);
  /* int *cgatom_num = (int*)PyArray_DATA(cgatom_num_array); */
  double *W = (double*)PyArray_DATA(W_array);
  
  /* Call the external C function to compute the chi-squared. */
  total_atom=cgSolvateNumeric_v2(fpdb_file, X, Y, Z, W, wDNA, wRNA, wPROT, thickness, closest_dist, pdb_flag, solvent_flag);
  
  /* /\* Clean up. *\/ */
  Py_DECREF(X_array);
  Py_DECREF(Y_array);
  Py_DECREF(Z_array);
  /* Py_DECREF(cgatom_num_array); */
  Py_DECREF(W_array);

  /* //  if (value < 0.0) { */
  /* //    PyErr_SetString(PyExc_RuntimeError, */
  /* //                    "Chi-squared returned an impossible value."); */
  /* //    return NULL; */
  /* //  } */

  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("i", total_atom);
  return ret;
}

//##################################FUNCTIONS####################################
float DS(float x)
{
  float a,b,c,d, temp;

  a = 3.5301846863022959E+01f;
  b = 8.5378345850638038E-01f;
  c = -1.7507051229758559E+01f;
  d = 1.7010560191799904E+00f;

  temp = a + b*x + c*x*x + d*x*x*x;
  return temp;
}

float RS(float x)
{
  float a,b,c,d, temp;

  a = 4.0230401562412005E+01f;
  b = 1.6287412345143366E+00f;
  c = -4.3118717065025962E+01f;
  d = 2.4712585246295806E+01f;

  temp = a + b*x + c*x*x + d*x*x*x;
  return temp;
}

float ADE(float x)
{
  float a,b,c,d, temp;

  a = 3.0473077532191034E+01f;
  b = 2.7267438703396873E-02f;
  c = -9.4174795997194121E+00f;
  d = -1.8632787145940051E+00f;

  temp = a + b*x + c*x*x + d*x*x*x;
  return temp;
}

float CYT(float x)
{
  float a,b,c,d, temp;

  a = 2.2594403303998384E+01f;
  b = -1.3979255030930837E-01f;
  c = 2.0602701684525548E-01f;
  d = -4.8964867374677956E+00f;

  temp = a + b*x + c*x*x + d*x*x*x;
  return temp;
}

float GUA(float x)
{
  float a,b,c,d, temp;

  a = 3.5413612820294524E+01f;
  b = 3.1135141139744515E-01f;
  c = -1.8995234746163160E+01f;
  d = 2.7546664066782145E+00f;

  temp = a + b*x + c*x*x + d*x*x*x;
  return temp;
}

float THY(float x)
{
  float a,b,c,d, temp;

  a = 2.1171045387193391E+01f;
  b = -2.5406818426818600E-01f;
  c = 8.5747955915722898E+00f;
  d = -1.2569339716672346E+01f;

  temp = a + b*x + c*x*x + d*x*x*x;
  return temp;
}

float URA(float x)
{
  float a,b,c,d, temp;

  a = 2.2101478120069594E+01f;
  b = -1.3059135164353608E-01f;
  c = 1.4058971332981493E-01f;
  d = -4.4641771309272844E+00f;

  temp = a + b*x + c*x*x + d*x*x*x;
  return temp;
}

float SOL(float x)
{
  float a,b,c,d, temp, x_sqrd=x*x;
  //  float temp, x_sqrd=x*x;
  a = 9.9989107859081461E+00f;
  b = 7.1847528938221994E-03f;
  c = -1.1436678741857631E+00f;
  d = 1.3146929900084206E-01f;
  //  x_sqrd=x*x;

  //  temp = a + b*x + c*x*x + d*x*x*x;
  temp = a + b*x + c*x_sqrd + d*x*x_sqrd;
  return temp;
}

float GLY(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 9.9691190155764158E+00f;
    b = -1.7528320201296514E-02f;
    c = 1.6450018157708988E+00f;
    d = -7.0378838370050911E-01f;
    f = -6.8495906226686154E-01f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float ALA(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 9.0385953871608269E+00f;
    b = -9.3522348098657587E-02f;
    c = 7.1811578286517612E+00f;
    d = -4.1100094434439107E+00f;
    f = -1.5204739994705183E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float VAL(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 7.1816423864658629E+00f;
    b = -4.7643636342215978E-01f;
    c = 2.6878626217077645E+01f;
    d = -3.0036144646732087E+01f;
    f = 7.6041571666551446E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float LEU(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 6.2519392408994952E+00f;
    b = -7.0384372758014790E-01f;
    c = 4.4699858222716180E+01f;
    d = -6.7435461752575918E+01f;
    f = 2.8904017698517240E+01f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float ILE(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 6.2512251193666968E+00f;
    b = -6.5569054939799687E-01f;
    c = 4.2523603838224325E+01f;
    d = -6.1110285327663291E+01f;
    f = 2.5111019303628737E+01f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}


float MET(float x)
 {
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.6545461441433780E+01f;
    b = -3.0728785802722175E-01f;
    c = 3.3057451408162395E+00f;
    d = -1.3719988781068070E+01f;
    f = 8.9621624126553510E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float PHE(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 9.2253730534361544E+00f;
    b = -1.0235816923538334E+00f;
    c = 4.0056964609302334E+01f;
    d = -6.6733866971096148E+01f;
    f = 3.1299821660208639E+01f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float TRP(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.5684405852404351E+01f;
    b = -9.4431662966303997E-01f;
    c = 2.3675802486236719E+01f;
    d = -4.8238875823804648E+01f;
    f = 2.4504594828865912E+01f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float PRO(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 8.6193900746569110E+00f;
    b = -3.2220720095487299E-01f;
    c = 1.8997202349961341E+01f;
    d = -1.6823509153368370E+01f;
    f = 1.4683423688962776E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float SER(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.3987709443384963E+01f;
    b = -5.3597112650067072E-02f;
    c = 1.6532170902823629E+00f;
    d = -2.1698155818409752E+00f;
    f = -1.0079162195407658E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}
float THR(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.3058035944130433E+01f;
    b = -1.6805767593043786E-01f;
    c = 8.4026246743192505E+00f;
    d = -7.3396417884658378E+00f;
    f = -1.4446678521071559E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float CYS(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.9124654812324472E+01f;
    b = -5.9927863522006228E-02f;
    c = -3.3836190674229298E+00f;
    d = -2.4665166405439520E+00f;
    f = 1.9122022270707959E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float TYR(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.4149376249271270E+01f;
    b = 2.3351533285852905E-01f;
    c = -6.6303117968389449E+00f;
    d = 7.7120460565621229E+00f;
    f = 1.2754732796902342E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float ASN(float x)
{
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.9938673457674692E+01f;
    b = -5.4516265978723853E-02f;
    c = -6.1462879847152454E+00f;
    d = -2.2543374567856151E+00f;
    f = 2.6914441182441431E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float GLN(float x)
 {
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 1.9005567570863825E+01f;
    b = -2.0319122385487057E-02f;
    c = -1.2083581566255768E+01f;
    d = -1.4429406093358430E+00f;
    f = 9.7256610585859082E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
 }

float ASP(float x)
 {
    float a,b,c,d,f,temp, x_sqrd=x*x;

    a = 2.0166013924521497E+01f;
    b = -4.7853403267681283E-02f;
    c = -6.9784835381725507E+00f;
    d = -1.9994797588755993E+00f;
    f = 3.2626175074094688E+00f;

    temp = a + b*x + c*x_sqrd + d*x_sqrd*x + f*x_sqrd*x_sqrd;
    return temp;
}

float GLU(float x)
 {
   float a,b,c,d,f,temp,x_sqrd=x*x;

    a = 1.9232515360604562E+01f;
    b = 3.5086173204513208E-03f;
    c = -1.2596769305364237E+01f;
    d = -4.9237630472478966E-01f;
    f = 9.3597483117941245E+00f;

    //    temp = a + b*x + c*x*x + d*x*x*x + f*x*x*x*x;
    temp = a + b*x + c*x_sqrd + d*x*x_sqrd + f*x_sqrd*x_sqrd;
    return temp;
 }

float LYS(float x)
 {
   float a,b,c,d,f,temp,x_sqrd=x*x;

    a = 1.0962754376635287E+01f;
    b = -3.1372344486261888E-02f;
    c = 7.3756682779136329E-01f;
    d = -3.0693869455133744E+00f;
    f = 5.0544680881861872E+00f;

    //    temp = a + b*x + c*x*x + d*x*x*x + f*x*x*x*x;
    temp = a + b*x + c*x_sqrd + d*x*x_sqrd + f*x_sqrd*x_sqrd;
    return temp;
}

float ARG(float x)
{
  float a,b,c,d,f,temp,x_sqrd=x*x;

    a = 2.3271735452541201E+01f;
    b = 6.7314944953355738E-01f;
    c = -3.6116270352466870E+01f;
    d = 2.1362171853330082E+01f;
    f = 1.6290629563851400E+01f;

    //    temp = a + b*x + c*x*x + d*x*x*x + f*x*x*x*x;
    temp = a + b*x + c*x_sqrd + d*x*x_sqrd + f*x_sqrd*x_sqrd;
    return temp;
}

float HIS(float x)
{
  float a,b,c,d,f,temp,x_sqrd=x*x;

    a = 2.1448835902382307E+01f;
    b = -9.8889323535974627E-02f;
    c = -7.0305612380447311E+00f;
    d = -4.4298193773300865E+00f;
    f = 6.2697213262776561E+00f;

    //    temp = a + b*x + c*x*x + d*x*x*x + f*x*x*x*x;
    temp = a + b*x + c*x_sqrd + d*x*x_sqrd + f*x_sqrd*x_sqrd;
    return temp;
}

float get_f(const char resname[], float w, float q)
{
  float f=0.0001f;

  if(strcmp(resname,"RNA")==0){f = RS(q);}
  if(strcmp(resname,"DNA")==0){f = DS(q);}
  if(strcmp(resname,"ADE")==0){f = ADE(q);}
  if(strcmp(resname,"GUA")==0){f = GUA(q);}
  if(strcmp(resname,"THY")==0){f = THY(q);}
  if(strcmp(resname,"URA")==0){f = URA(q);}
  if(strcmp(resname,"CYT")==0){f = CYT(q);}

  if(strcmp(resname,"GLY")==0){f = GLY(q);}
  if(strcmp(resname,"ALA")==0){f = ALA(q);}
  if(strcmp(resname,"VAL")==0){f = VAL(q);}
  if(strcmp(resname,"ILE")==0){f = ILE(q);}
  if(strcmp(resname,"LEU")==0){f = LEU(q);}
  if(strcmp(resname,"MET")==0){f = MET(q);}
  if(strcmp(resname,"PHE")==0){f = PHE(q);}
  if(strcmp(resname,"TRP")==0){f = TRP(q);}
  if(strcmp(resname,"PRO")==0){f = PRO(q);}
  if(strcmp(resname,"SER")==0){f = SER(q);}
  if(strcmp(resname,"THR")==0){f = THR(q);}
  if(strcmp(resname,"CYS")==0){f = CYS(q);}
  if(strcmp(resname,"TYR")==0){f = TYR(q);}
  if(strcmp(resname,"ASN")==0){f = ASN(q);}
  if(strcmp(resname,"GLN")==0){f = GLN(q);}
  if(strcmp(resname,"ASP")==0){f = ASP(q);}
  if(strcmp(resname,"GLU")==0){f = GLU(q);}
  if(strcmp(resname,"LYS")==0){f = LYS(q);}
  if(strcmp(resname,"ARG")==0){f = ARG(q);}
  if(strcmp(resname,"HIS")==0){f = HIS(q);}
  if(strcmp(resname,"HSD")==0){f = HIS(q);}

  if(strcmp(resname,"SOL")==0){f = w*SOL(q);}
  if(strcmp(resname,"TIP")==0){f = w*SOL(q);}
  if(strcmp(resname,"SPC")==0){f = w*SOL(q);}
  //printf("%s %f %f %f\n", resname,w,q,f);
  return f;
}

float get_f_numeric(int resname_num, float w, float q)
{
  float f=0.0001f;

  if(resname_num==1){f = RS(q);}
  if(resname_num==2){f = DS(q);}
  if(resname_num==3){f = ADE(q);}
  if(resname_num==4){f = GUA(q);}
  if(resname_num==5){f = THY(q);}
  if(resname_num==6){f = URA(q);}
  if(resname_num==7){f = CYT(q);}

  if(resname_num==8){f = GLY(q);}
  if(resname_num==9){f = ALA(q);}
  if(resname_num==10){f = VAL(q);}
  if(resname_num==11){f = ILE(q);}
  if(resname_num==12){f = LEU(q);}
  if(resname_num==13){f = MET(q);}
  if(resname_num==14){f = PHE(q);}
  if(resname_num==15){f = TRP(q);}
  if(resname_num==16){f = PRO(q);}
  if(resname_num==17){f = SER(q);}
  if(resname_num==18){f = THR(q);}
  if(resname_num==19){f = CYS(q);}
  if(resname_num==20){f = TYR(q);}
  if(resname_num==21){f = ASN(q);}
  if(resname_num==22){f = GLN(q);}
  if(resname_num==23){f = ASP(q);}
  if(resname_num==24){f = GLU(q);}
  if(resname_num==25){f = LYS(q);}
  if(resname_num==26){f = ARG(q);}
  if(resname_num==27){f = HIS(q);}
  if(resname_num==28){f = HIS(q);}

  if(resname_num==29){f = w*SOL(q);}
  if(resname_num==30){f = w*SOL(q);}
  if(resname_num==31){f = w*SOL(q);}
  //  printf("Here: %d %f %f %f\n", resname_num,w,q,f);
  return f;
}


/* char *resnameNum_to_resname(int resnameNum, float w) */
/* { */
/*   char f[4]={'\0', '\0', '\0', '\0'}; */

/*   if(resname_num==1){f = RS;} */
/*   if(resname_num==2){f = DS;} */
/*   if(resname_num==3){f = ADE;} */
/*   if(resname_num==4){f = GUA;} */
/*   if(resname_num==5){f = THY;} */
/*   if(resname_num==6){f = URA;} */
/*   if(resname_num==7){f = CYT;} */

/*   if(resname_num==8){f = GLY;} */
/*   if(resname_num==9){f = ALA;} */
/*   if(resname_num==10){f = VAL;} */
/*   if(resname_num==11){f = ILE;} */
/*   if(resname_num==12){f = LEU;} */
/*   if(resname_num==13){f = MET;} */
/*   if(resname_num==14){f = PHE;} */
/*   if(resname_num==15){f = TRP;} */
/*   if(resname_num==16){f = PRO;} */
/*   if(resname_num==17){f = SER;} */
/*   if(resname_num==18){f = THR;} */
/*   if(resname_num==19){f = CYS;} */
/*   if(resname_num==20){f = TYR;} */
/*   if(resname_num==21){f = ASN;} */
/*   if(resname_num==22){f = GLN;} */
/*   if(resname_num==23){f = ASP;} */
/*   if(resname_num==24){f = GLU;} */
/*   if(resname_num==25){f = LYS;} */
/*   if(resname_num==26){f = ARG;} */
/*   if(resname_num==27){f = HIS;} */
/*   if(resname_num==28){f = HIS;} */

/*   if(resname_num==29){f = SOL;} */
/*   if(resname_num==30){f = SOL;} */
/*   if(resname_num==31){f = SOL;} */
/*   //printf("%s %f %f %f\n", resname,w,q,f); */
/*   return f; */
/* } */

// function to find minimum distance betwen water and
// molecule.
float min_dist_to_mol_v2(double *X, double *Y, double *Z, \
			 double x, double y, double z,		\
			 int atom_count, int *mol_type, int *type)
{
  int i=0;
  float dist = 0.0f;
  float dist_x=0.0f, dist_y=0.0f, dist_z=0.0f;
  //  float min_dist = 10000f;
  float min_dist = 100000000.0f;
  for(i=0;i<atom_count;i++)
    {
      //      dist  = (X[i]-x)*(X[i]-x);
      //      dist += (Y[i]-y)*(Y[i]-y);
      //dist += (Z[i]-z)*(Z[i]-z);
      dist_x= (float) (X[i]-x);
      dist_y= (float) (Y[i]-y);
      dist_z= (float) (Z[i]-z);
      dist=((dist_x*dist_x)+(dist_y*dist_y)+(dist_z*dist_z));
      //      dist = sqrt(dist);
      //      if (dist < 3.0)
      if (dist < 9.0)
	{
	  min_dist = 0;
	  break;
	}
      if (dist < min_dist)
	{
	  min_dist = dist;
	  *type = mol_type[i];
	}
    }
  //  return sqrt(min_dist);
  return min_dist;
}


/* float min_dist_to_mol_v3(double *X, double *Y, double *Z, \ */
/* 			 double x, double y, double z,	  \ */
/* 			 int atom_count, int *mol_type) */
/* { */
/*   int i; */
/*   float dist = 0.0; */
/*   float dist_x=0.0, dist_y=0.0, dist_z=0.0; */
/*   //  float min_dist = 10000; */
/*   float min_dist = 100000000; */
/*   for(i=0;i<atom_count;i++) */
/*     { */
/*       //      dist  = (X[i]-x)*(X[i]-x); */
/*       //      dist += (Y[i]-y)*(Y[i]-y); */
/*       //dist += (Z[i]-z)*(Z[i]-z); */
/*       dist_x= (X[i]-x); */
/*       dist_y= (Y[i]-y); */
/*       dist_z= (Z[i]-z); */
/*       dist=((dist_x*dist_x)+(dist_y*dist_y)+(dist_z*dist_z)); */
/*       //      dist = sqrt(dist); */
/*       //      if (dist < 3.0) */
/*       if (dist < 9.0)  */
/* 	{ */
/* 	  min_dist = 0; */
/* 	  break; */
/* 	} */
/*       if (dist < min_dist) */
/* 	{ */
/* 	  min_dist = dist; */
/* 	  *type = mol_type[i]; */
/* 	} */
/*     } */
/*   //  return sqrt(min_dist); */
/*   return min_dist; */
/* } */


// Function to Coarse-grain the input PDB and solvate
int cgSolvate(FILE *fpdb, double *X, double *Y, double *Z, int *total_atom, char *cgatom[], double *W, \
	       float delta, float wDNA, float wRNA, float wPROT,  int pdb_flag, \
	       float thickness, float closest_dist, int solvent_flag)
{

  char *pdb_out_file="cg_solvated.pdb";
  float x,y,z;
  char resname[4]="   \0", atomname[5]="    \0";
  char cx[9]="        \0", cy[9]="        \0", cz[9]="        \0";
  char pdb_line[80];
  int mol_type_flag=0, flag=0;
  float min_dist=0;
  float maxx=-10000, maxy=-10000, maxz=-10000;
  float minx=10000, miny=10000, minz=10000;
  int atom_count=0;
  char *atom[MAX_ATOM];
  int RNA_count=0, DNA_count=0, PROT_count=0;
  int mol_type[MAX_ATOM];
  int w,i,j,k,mtype=0;
  int water_count;
  FILE *fout_pdb=NULL;

  float closest_dist_plus_thickness=(closest_dist+thickness);
  float closest_dist_plus_thickness_sqrd=closest_dist_plus_thickness*closest_dist_plus_thickness;
  float closest_dist_sqrd=closest_dist*closest_dist;


  // Read PDB file and coarge-grain the file
  while(fgets(pdb_line, 80, fpdb) != NULL)
    {
      if(strncmp(pdb_line,"ATOM",4)==0)
	{
	  strncpy(resname,pdb_line+17, 3);
	  strncpy(atomname,pdb_line+12, 4);

	  //mol_type_flag:  0-DNA, 1-RNA, 2-Protein
	  if(strcmp(atomname," O2'")==0 || strcmp(atomname," O2*")==0 || strcmp(resname,"RNA")==0 )
	    {
	      mol_type_flag=1;
	      if(RNA_count==0 && strcmp(resname,"RNA")!=0 && DNA_count==1)
		{
		  cgatom[atom_count-1]="RNA\0";
		  mol_type[atom_count-1]=1;
		  W[atom_count-1]=wRNA;
		  RNA_count++;
		  DNA_count--;
		}
	    }

	  if(strcmp(atomname," O5'")==0 || strcmp(atomname," O5*")==0)
	    {
	      if(mol_type_flag==1)
		{
		  cgatom[atom_count]="RNA\0";
		  RNA_count++;
		}
	      else
		{
		  cgatom[atom_count]="DNA\0";
		  DNA_count++;
		  mol_type_flag=0;
		}
	      atom[atom_count]=" O5'\0";
	      flag=1;
	    }
          if((strcmp(resname,"  A")==0 || strcmp(resname," DA")==0 || strcmp(resname,"DA ")==0 \
	      || strcmp(resname," RA")==0 || strcmp(resname,"RA ")==0 || strcmp(resname,"ADE")==0 ) \
	     && (strcmp(atomname," C5 ")==0 || strcmp(atomname,"  C5")==0))
	    {cgatom[atom_count]="ADE\0";flag=1;atom[atom_count]=" C5 \0";}

	  if((strcmp(resname,"  C")==0 || strcmp(resname," DC")==0 || strcmp(resname,"DC ")==0 \
	      || strcmp(resname," RC")==0 || strcmp(resname,"RC ")==0 || strcmp(resname,"CYT")==0 ) \
	     && (strcmp(atomname," N3 ")==0 || strcmp(atomname,"  N3")==0))
	    {cgatom[atom_count]="CYT\0";flag=1;atom[atom_count]=" N3 \0";}

          if((strcmp(resname,"  G")==0 || strcmp(resname," DG")==0 || strcmp(resname,"DG ")==0 \
	      || strcmp(resname," RG")==0 || strcmp(resname,"RG ")==0 || strcmp(resname,"GUA")==0 ) \
	     && (strcmp(atomname," C4 ")==0 || strcmp(atomname,"  C4")==0))
	    {cgatom[atom_count]="GUA\0";flag=1;atom[atom_count]=" C4 \0";}

          if((strcmp(resname,"  U")==0 || strcmp(resname," DU")==0 || strcmp(resname,"DU ")==0 \
	      || strcmp(resname," RU")==0 || strcmp(resname,"RU ")==0 || strcmp(resname,"URA")==0 ) \
	     && (strcmp(atomname," N3 ")==0 || strcmp(atomname,"  N3")==0))
	    {cgatom[atom_count]="URA\0";flag=1;atom[atom_count]=" N3 \0";}

          if((strcmp(resname,"  T")==0 || strcmp(resname," DT")==0 || strcmp(resname,"DT ")==0 \
	      || strcmp(resname," RT")==0 || strcmp(resname,"RT ")==0 || strcmp(resname,"THY")==0 ) \
	     && (strcmp(atomname," N3 ")==0 || strcmp(atomname,"  N3")==0))
	    {cgatom[atom_count]="THY\0";flag=1;atom[atom_count]=" N3 \0";}

	  //Check if atom belongs to protein
         if(strcmp(atomname," CA ")==0 || strcmp(atomname,"  CA")==0 || strcmp(atomname,"CA  ")==0 )
	   {
	     cgatom[atom_count]="   \0";
	     if(strcmp(resname,"GLY")==0) cgatom[atom_count]="GLY\0";
	     if(strcmp(resname,"ALA")==0) cgatom[atom_count]="ALA\0";
	     if(strcmp(resname,"VAL")==0) cgatom[atom_count]="VAL\0";
	     if(strcmp(resname,"LEU")==0) cgatom[atom_count]="LEU\0";
	     if(strcmp(resname,"ILE")==0) cgatom[atom_count]="ILE\0";
	     if(strcmp(resname,"MET")==0) cgatom[atom_count]="MET\0";
	     if(strcmp(resname,"PHE")==0) cgatom[atom_count]="PHE\0";
	     if(strcmp(resname,"TRP")==0) cgatom[atom_count]="TRP\0";
	     if(strcmp(resname,"PRO")==0) cgatom[atom_count]="PRO\0";
	     if(strcmp(resname,"SER")==0) cgatom[atom_count]="SER\0";
	     if(strcmp(resname,"THR")==0) cgatom[atom_count]="THR\0";
	     if(strcmp(resname,"CYS")==0) cgatom[atom_count]="CYS\0";
	     if(strcmp(resname,"TYR")==0) cgatom[atom_count]="TYR\0";
	     if(strcmp(resname,"ASN")==0) cgatom[atom_count]="ASN\0";
	     if(strcmp(resname,"GLN")==0) cgatom[atom_count]="GLN\0";
	     if(strcmp(resname,"ASP")==0) cgatom[atom_count]="ASP\0";
	     if(strcmp(resname,"GLU")==0) cgatom[atom_count]="GLU\0";
	     if(strcmp(resname,"LYS")==0) cgatom[atom_count]="LYS\0";
	     if(strcmp(resname,"ARG")==0) cgatom[atom_count]="ARG\0";
	     if(strcmp(resname,"HIS")==0) cgatom[atom_count]="HIS\0";
	     if(strcmp(resname,"HSD")==0) cgatom[atom_count]="HIS\0";
	     if(strcmp(resname,"HYP")==0) cgatom[atom_count]="PRO\0";
	      PROT_count++;
	      mol_type_flag=2;
	      flag=1;
	      atom[atom_count]=" CA \0";
	   }

	  if(flag==1)
	    {
	      strncpy(cx,pdb_line+30, 8);
	      strncpy(cy,pdb_line+38, 8);
	      strncpy(cz,pdb_line+46, 8);

        x = (float) atof(cx);
        y = (float) atof(cy);
        z = (float) atof(cz);

	      if(x > maxx) maxx = x;
	      if(y > maxy) maxy = y;
	      if(z > maxz) maxz = z;

	      if(x < minx) minx = x;
	      if(y < miny) miny = y;
	      if(z < minz) minz = z;

	      X[atom_count] = x;
	      Y[atom_count] = y;
	      Z[atom_count] = z;

	      mol_type[atom_count]=mol_type_flag;
	      if (mol_type_flag==0) {W[atom_count] = wDNA;}
	      if (mol_type_flag==1) {W[atom_count] = wRNA;}
	      if (mol_type_flag==2) {W[atom_count] = wPROT;}

	      atom_count++;
	      flag=0;

	    }
	}
    }

  // Output pdb if required
  if(pdb_flag ==1)
    {
      fout_pdb = fopen(pdb_out_file, "w");
      for(i=0;i<atom_count;i++)
	{
	  fprintf(fout_pdb,"ATOM %6d %4s %3s %1d%4d    %8.3f%8.3f%8.3f\n", \
	       i+1,atom[i],cgatom[i],mol_type[i],i+1,X[i],Y[i],Z[i]);
	}

    }


  // add water box if solvent_flag==1
  // This routine can add water box of any size
  water_count = 1;
  if(solvent_flag == 1)
    {
      for(i=0; i< (maxx-minx)/WATER_BOX_SIZE+1; i++)
	{
	  for(j=0; j< (maxy-miny)/WATER_BOX_SIZE+1; j++)
	    {
	      for(k=0; k< (maxz-minz)/WATER_BOX_SIZE+1; k++)
		{
		  for(w=0;w<NUM_WATER_ATOMS;w++)
		    {

		      x = water[w][0]+minx + (i*WATER_BOX_SIZE);
		      y = water[w][1]+miny + (j*WATER_BOX_SIZE);
		      z = water[w][2]+minz + (k*WATER_BOX_SIZE);
		      min_dist = min_dist_to_mol_v2(X,Y,Z,x,y,z,atom_count,mol_type,&mtype);

		      if (mtype==0) {W[atom_count+water_count-1] = wDNA;}
		      if (mtype==1) {W[atom_count+water_count-1] = wRNA;}
		      if (mtype==2) {W[atom_count+water_count-1] = wPROT;}
		      cgatom[atom_count+water_count-1] = "SOL\0";
		      X[atom_count+water_count-1] = x;
		      Y[atom_count+water_count-1] = y;
		      Z[atom_count+water_count-1] = z;

		      if(min_dist < closest_dist_sqrd) {continue;}
		      if(min_dist > closest_dist_plus_thickness_sqrd) {continue;}

		      if(pdb_flag ==1)
			{
			  fprintf(fout_pdb,"ATOM %6d  OW  %3s %1d%4d    %8.3f%8.3f%8.3f\n", \
				  water_count,cgatom[atom_count+water_count-1],mtype,water_count,x,y,z);
			}
		      water_count++;

		      if(water_count==10000)
			water_count =1;
		    }
		}
	    }
	}
    }

  // close pdb file
  if(pdb_flag ==1) fclose(fout_pdb);

  /* fprintf(stderr,"\n@> ******* Summary **************************\n\n"); */
  /* fprintf(stderr,"@> Total number of RNA Nucleotides = %d\n",RNA_count); */
  /* fprintf(stderr,"@> Total number of DNA Nucleotides = %d\n",DNA_count); */
  /* fprintf(stderr,"@> Total number of Protein Residues = %d\n",PROT_count); */
  /* fprintf(stderr,"@> Total non-water CG atoms in the pdb file = %d\n",atom_count); */
  /* fprintf(stderr,"@> Total number of CG water atoms added = %d\n",water_count-1); */
  /* fprintf(stderr,"@> Total number of CG atoms in the pdb file = %d\n",atom_count+water_count-1); */
  /* fprintf(stderr,"\n@> ******************************************\n"); */

  if((atom_count-PROT_count)%2 != 0)
    {
      fprintf(stderr,"WARNING: Something wrong with the input PDB file\n");
      fprintf(stderr,"Perhaps some nucleotides have missing atoms?\n\n");
      fprintf(stderr,"Check your input PDB file. \n\n");
      fprintf(stderr,"Read documentation: -H option for more details. \n\n");
    }

  *total_atom = atom_count + water_count - 1;
  return 0;
}

int cgSolvateNumeric_v2(char *fpdb_file,\
			double *X,\
			double *Y,\
			double *Z,\
			double *W,\
			float wDNA,\
			float wRNA,\
			float wPROT,\
			float thickness,\
			float closest_dist,\
			int pdb_flag,\
			int solvent_flag)
{
  //  fprintf(stdout, "Here I am and opening string is %s\n", fpdb_file);
  FILE *fpdb=fopen(fpdb_file, "r");
  if(fpdb==NULL)
    {
      fprintf(stderr, "ERROR! No such file: %s", fpdb_file);
      exit(EXIT_FAILURE);
    }

  int total_atom=0;
  char *pdb_out_file="cg_solvated.pdb";
  float x,y,z;
  char resname[4]="   \0", atomname[5]="    \0";
  char cx[9]="        \0", cy[9]="        \0", cz[9]="        \0";
  char pdb_line[80];
  int mol_type_flag=0, flag=0;
  float min_dist=0;
  float maxx=-10000.0f, maxy=-10000.0f, maxz=-10000.0f;
  float minx=10000.0f, miny=10000.0f, minz=10000.0f;
  int atom_count=0;
  char *atom[MAX_ATOM];
  char *cgatom[MAX_ATOM];

  int RNA_count=0, DNA_count=0, PROT_count=0;
  int mol_type[MAX_ATOM];
  int w,i,j,k,mtype=0;
  int water_count;
  FILE *fout_pdb=NULL;

  float closest_dist_plus_thickness=(closest_dist+thickness);
  float closest_dist_plus_thickness_sqrd=closest_dist_plus_thickness*closest_dist_plus_thickness;
  float closest_dist_sqrd=closest_dist*closest_dist;

  //  fprintf(stdout, "Here I am 2\n");
  // Read PDB file and coarge-grain the file
  while(fgets(pdb_line, 80, fpdb) != NULL)
    {
      if(strncmp(pdb_line,"ATOM",4)==0)
	{
	  strncpy(resname,pdb_line+17, 3);
	  strncpy(atomname,pdb_line+12, 4);

	  //mol_type_flag:  0-DNA, 1-RNA, 2-Protein
	  if(strcmp(atomname," O2'")==0 || strcmp(atomname," O2*")==0 || strcmp(resname,"RNA")==0 )
	    {
	      mol_type_flag=1;
	      if(RNA_count==0 && strcmp(resname,"RNA")!=0 && DNA_count==1)
		{
		  cgatom[atom_count-1]="RNA\0";
		  mol_type[atom_count-1]=1;
		  W[atom_count-1]=wRNA;
		  RNA_count++;
		  DNA_count--;
		}
	    }

	  if(strcmp(atomname," O5'")==0 || strcmp(atomname," O5*")==0)
	    {
	      if(mol_type_flag==1)
		{
		  cgatom[atom_count]="RNA\0";
		  RNA_count++;
		}
	      else
		{
		  cgatom[atom_count]="DNA\0";
		  DNA_count++;
		  mol_type_flag=0;
		}
	      atom[atom_count]=" O5'\0";
	      flag=1;
	    }
          if((strcmp(resname,"  A")==0 || strcmp(resname," DA")==0 || strcmp(resname,"DA ")==0 \
	      || strcmp(resname," RA")==0 || strcmp(resname,"RA ")==0 || strcmp(resname,"ADE")==0 ) \
	     && (strcmp(atomname," C5 ")==0 || strcmp(atomname,"  C5")==0))
	    {cgatom[atom_count]="ADE\0";flag=1;atom[atom_count]=" C5 \0";}

	  if((strcmp(resname,"  C")==0 || strcmp(resname," DC")==0 || strcmp(resname,"DC ")==0 \
	      || strcmp(resname," RC")==0 || strcmp(resname,"RC ")==0 || strcmp(resname,"CYT")==0 ) \
	     && (strcmp(atomname," N3 ")==0 || strcmp(atomname,"  N3")==0))
	    {cgatom[atom_count]="CYT\0";flag=1;atom[atom_count]=" N3 \0";}

          if((strcmp(resname,"  G")==0 || strcmp(resname," DG")==0 || strcmp(resname,"DG ")==0 \
	      || strcmp(resname," RG")==0 || strcmp(resname,"RG ")==0 || strcmp(resname,"GUA")==0 ) \
	     && (strcmp(atomname," C4 ")==0 || strcmp(atomname,"  C4")==0))
	    {cgatom[atom_count]="GUA\0";flag=1;atom[atom_count]=" C4 \0";}

          if((strcmp(resname,"  U")==0 || strcmp(resname," DU")==0 || strcmp(resname,"DU ")==0 \
	      || strcmp(resname," RU")==0 || strcmp(resname,"RU ")==0 || strcmp(resname,"URA")==0 ) \
	     && (strcmp(atomname," N3 ")==0 || strcmp(atomname,"  N3")==0))
	    {cgatom[atom_count]="URA\0";flag=1;atom[atom_count]=" N3 \0";}

          if((strcmp(resname,"  T")==0 || strcmp(resname," DT")==0 || strcmp(resname,"DT ")==0 \
	      || strcmp(resname," RT")==0 || strcmp(resname,"RT ")==0 || strcmp(resname,"THY")==0 ) \
	     && (strcmp(atomname," N3 ")==0 || strcmp(atomname,"  N3")==0))
	    {cgatom[atom_count]="THY\0";flag=1;atom[atom_count]=" N3 \0";}

	  //Check if atom belongs to protein
         if(strcmp(atomname," CA ")==0 || strcmp(atomname,"  CA")==0 || strcmp(atomname,"CA  ")==0 )
	   {
	     cgatom[atom_count]="   \0";
	     if(strcmp(resname,"GLY")==0) cgatom[atom_count]="GLY\0";
	     if(strcmp(resname,"ALA")==0) cgatom[atom_count]="ALA\0";
	     if(strcmp(resname,"VAL")==0) cgatom[atom_count]="VAL\0";
	     if(strcmp(resname,"LEU")==0) cgatom[atom_count]="LEU\0";
	     if(strcmp(resname,"ILE")==0) cgatom[atom_count]="ILE\0";
	     if(strcmp(resname,"MET")==0) cgatom[atom_count]="MET\0";
	     if(strcmp(resname,"PHE")==0) cgatom[atom_count]="PHE\0";
	     if(strcmp(resname,"TRP")==0) cgatom[atom_count]="TRP\0";
	     if(strcmp(resname,"PRO")==0) cgatom[atom_count]="PRO\0";
	     if(strcmp(resname,"SER")==0) cgatom[atom_count]="SER\0";
	     if(strcmp(resname,"THR")==0) cgatom[atom_count]="THR\0";
	     if(strcmp(resname,"CYS")==0) cgatom[atom_count]="CYS\0";
	     if(strcmp(resname,"TYR")==0) cgatom[atom_count]="TYR\0";
	     if(strcmp(resname,"ASN")==0) cgatom[atom_count]="ASN\0";
	     if(strcmp(resname,"GLN")==0) cgatom[atom_count]="GLN\0";
	     if(strcmp(resname,"ASP")==0) cgatom[atom_count]="ASP\0";
	     if(strcmp(resname,"GLU")==0) cgatom[atom_count]="GLU\0";
	     if(strcmp(resname,"LYS")==0) cgatom[atom_count]="LYS\0";
	     if(strcmp(resname,"ARG")==0) cgatom[atom_count]="ARG\0";
	     if(strcmp(resname,"HIS")==0) cgatom[atom_count]="HIS\0";
	     if(strcmp(resname,"HSD")==0) cgatom[atom_count]="HIS\0";
	     if(strcmp(resname,"HYP")==0) cgatom[atom_count]="PRO\0";
	      PROT_count++;
	      mol_type_flag=2;
	      flag=1;
	      atom[atom_count]=" CA \0";
	   }

	 if(flag==1)
	   {
	     strncpy(cx,pdb_line+30, 8);
	     strncpy(cy,pdb_line+38, 8);
	     strncpy(cz,pdb_line+46, 8);

       x = (float) atof(cx);
       y = (float) atof(cy);
       z = (float) atof(cz);

	     if(x > maxx) maxx = x;
	     if(y > maxy) maxy = y;
	     if(z > maxz) maxz = z;

	     if(x < minx) minx = x;
	     if(y < miny) miny = y;
	     if(z < minz) minz = z;

	     X[atom_count] = x;
	     Y[atom_count] = y;
	     Z[atom_count] = z;

	     mol_type[atom_count]=mol_type_flag;
	     if (mol_type_flag==0) {W[atom_count] = wDNA;}
	     if (mol_type_flag==1) {W[atom_count] = wRNA;}
	     if (mol_type_flag==2) {W[atom_count] = wPROT;}

	     atom_count++;
	     flag=0;

	   }
	}
    }
  //  fprintf(stdout, "Here I am 3\n");
  // Output pdb if required
  if(pdb_flag==1)
    {
      fout_pdb = fopen(pdb_out_file, "w");
      for(i=0;i<atom_count;i++)
	{
	  fprintf(fout_pdb,"ATOM %6d %4s %3s %1d%4d    %8.3f%8.3f%8.3f\n", \
	       i+1,atom[i],cgatom[i],mol_type[i],i+1,X[i],Y[i],Z[i]);
	}

    }
  //  fprintf(stdout, "Here I am 4\n");

  // add water box if solvent_flag==1
  // This routine can add water box of any size
  water_count = 1;
  if(solvent_flag == 1)
    {
      //      printf("Aslinda buradayim ulen!\n");
      for(i=0; i< (maxx-minx)/WATER_BOX_SIZE+1; i++)
	{
	  for(j=0; j< (maxy-miny)/WATER_BOX_SIZE+1; j++)
	    {
	      for(k=0; k< (maxz-minz)/WATER_BOX_SIZE+1; k++)
		{
		  for(w=0;w<NUM_WATER_ATOMS;w++)
		    {

		      x = water[w][0]+minx + (i*WATER_BOX_SIZE);
		      y = water[w][1]+miny + (j*WATER_BOX_SIZE);
		      z = water[w][2]+minz + (k*WATER_BOX_SIZE);
		      min_dist = min_dist_to_mol_v2(X,Y,Z,x,y,z,atom_count,mol_type,&mtype);

		      if (mtype==0) {W[atom_count+water_count-1] = wDNA;}
		      if (mtype==1) {W[atom_count+water_count-1] = wRNA;}
		      if (mtype==2) {W[atom_count+water_count-1] = wPROT;}
		      cgatom[atom_count+water_count-1] = "SOL\0";
		      X[atom_count+water_count-1] = x;
		      Y[atom_count+water_count-1] = y;
		      Z[atom_count+water_count-1] = z;

		      if(min_dist < closest_dist_sqrd) {continue;}
		      if(min_dist > closest_dist_plus_thickness_sqrd) {continue;}

		      if(pdb_flag ==1)
			{
			  fprintf(fout_pdb,"ATOM %6d  OW  %3s %1d%4d    %8.3f%8.3f%8.3f\n", \
				  water_count,cgatom[atom_count+water_count-1],mtype,water_count,x,y,z);
			}
		      water_count++;

		      if(water_count==10000)
			water_count =1;
		    }
		}
	    }
	}
    }
  //  fprintf(stdout, "Here I am 5\n");
  // close pdb file
   if(pdb_flag ==1) fclose(fout_pdb);

  /* fprintf(stderr,"\n@> ******* Summary **************************\n\n"); */
  /* fprintf(stderr,"@> !Total number of RNA Nucleotides = %d\n",RNA_count); */
  /* fprintf(stderr,"@> !Total number of DNA Nucleotides = %d\n",DNA_count); */
  /* fprintf(stderr,"@> !Total number of Protein Residues = %d\n",PROT_count); */
  /* fprintf(stderr,"@> !Total non-water CG atoms in the pdb file = %d\n",atom_count); */
  /* fprintf(stderr,"@> !Total number of CG water atoms added = %d\n",water_count-1); */
  /* fprintf(stderr,"@> !Total number of CG atoms in the pdb file = %d\n",atom_count+water_count-1); */
  /* fprintf(stderr,"\n@> ******************************************\n"); */

  if((atom_count-PROT_count)%2 != 0)
    {
      fprintf(stderr,"WARNING: Something wrong with the input PDB file\n");
      fprintf(stderr,"Perhaps some nucleotides have missing atoms?\n\n");
      fprintf(stderr,"Check your input PDB file. \n\n");
      fprintf(stderr,"Read documentation: -H option for more details. \n\n");
    }

  total_atom = atom_count + water_count - 1;
  fclose(fpdb);
  return (total_atom);
}


//  Function to Calculate SAXS Intensity
void calcSAXS(float *I, float *X,float *Y, float *Z, int total_atom, const char *cgatom[], float *W, float *Q, int nexp)
{
  int i,j,n,k = 0;
  //  float q, fi, fj, r;
  float q, fi;
  float qx,qy,qz, qdotr;
  float A,B,pi,GR;
  int nstart, nstop;
  float Two_jbyn=0.0;
  float Two_pi_jbyGR=0.0;
  float cosasin=0.0;
  n = 201;
  nstart = -(n-1)/2;
  nstop = (n-1)/2;

  pi= 3.1415f;
  GR = (float) ((1 + sqrt(5.0))/2.0);

  for(k=0;k<nexp;k++)
    {
      q = Q[k];
      I[k] = 0.0;

      for(j=nstart;j<=nstop;j++)
	{

	  Two_jbyn = 2.0f*j/n;
	  Two_pi_jbyGR = 2*pi*j/GR;
    cosasin = (float) (q*cos( asin(Two_jbyn) ));
	  //	  qx = q * cos( asin(2.0*j/n) ) * cos(2*pi*j/GR);
    qx = (float) (cosasin * cos(Two_pi_jbyGR));
	  //	  qy = q * cos( asin(2.0*j/n) ) * sin(2*pi*j/GR);
    qy = (float) (cosasin * sin(Two_pi_jbyGR));
	  qz = q * Two_jbyn;
	  A = 0.0;
	  B = 0.0;

	  for(i=0;i<total_atom;i++)
	    {
	      //if( !(strcmp(cgatom[j],"SOL")==0 && strcmp(cgatom[i],"SOL")==0) )
		{
		  fi = get_f(cgatom[i],W[i],q);
		  qdotr = qx*X[i] + qy*Y[i] + qz*Z[i];
      A = A + (float) (fi*cos(qdotr));
      B = B + (float) (fi*sin(qdotr));
		}
	    }
	  A = A*A;
	  B = B*B;
	  I[k] = I[k] + A+B;
	}
      I[k] = I[k]/n;

      I[k] = (float) log10(I[k]);
    }
}

//  Function to Calculate SAXS intensities by using Fast-SAXS method. 
void calcSAXSNumeric(double *I, double *X,double *Y, double *Z, int total_atom, int *cgatom_num, double *W, double *Q, int nexp)
{
  int i,j,n,k = 0;
  int printDataFile=0;
  //  float q, fi, fj, r;
  float q, fi;
  float qx,qy,qz, qdotr;
  float A,B,pi,GR;
  int nstart, nstop;
  float Two_jbyn=0.0;
  float Two_pi_jbyGR=0.0;


  float cosasin=0.0;

  n = 201;
  nstart = -(n-1)/2;
  nstop = (n-1)/2;

  pi= 3.1415f;
  GR = (float) ((1 + sqrt(5.0))/2.0);
  float Two_byn=2.0f/n;
  float Two_pi_byGR=2*pi/GR;
  FILE *SAXS_PROFILE=NULL;
  if(printDataFile)
    {
      SAXS_PROFILE=fopen("saxs_profile.dat", "w");
      if(SAXS_PROFILE==NULL)
	{
	  fprintf(stderr, "ERROR: Could not open %s file for writing\n", "saxs_profile.dat");
	  exit(EXIT_FAILURE);
	}
    }
  //  printf("Total number of atoms %d\n", total_atom);
  //  fprintf(stdout, "%.3f\t%.3f\t%.3f\n", Q[0], Q[1], Q[2]);

  for(k=0;k<nexp;k++)
    //  for(k=0;k<1;k++)
    {
      q = (float) Q[k];
      I[k] = 0.0;

      for(j=nstart;j<=nstop;j++)
	{

	  //	  Two_jbyn=2.0*j/n;
	  Two_jbyn=Two_byn*j;
	  //	  Two_pi_jbyGR=2*pi*j/GR;
	  Two_pi_jbyGR=	Two_pi_byGR*j;
    cosasin= (float) (q*cos( asin(Two_jbyn) ));
	  //	  qx = q * cos( asin(2.0*j/n) ) * cos(2*pi*j/GR);
    qx = (float) (cosasin * cos(Two_pi_jbyGR));
	  //	  qy = q * cos( asin(2.0*j/n) ) * sin(2*pi*j/GR);
    qy = (float) (cosasin * sin(Two_pi_jbyGR));
	  qz = q * Two_jbyn;
	  A = 0.0;
	  B = 0.0;
	  //	  printf("%.3f\t%.3f\t%.3f\n", qx, qy, qz);
	  //	  for(i=0;i<total_atom;i++)
	  for(i=0;i<total_atom;i++)
	    {
	      //if( !(strcmp(cgatom[j],"SOL")==0 && strcmp(cgatom[i],"SOL")==0) )
	      //		{
		  //		  printf("%.3lf\n", W[i]);
	      //		  printf("%d\n", cgatom_num[i]);
      fi = get_f_numeric(cgatom_num[i], (float) W[i], q);

      qdotr = (float) (qx*X[i] + qy*Y[i] + qz*Z[i]);
		  //	  printf("%.3f\t%.3f\n", fi, qdotr);
      A = A + (float) (fi*cos(qdotr));
      B = B + (float) (fi*sin(qdotr));
		  //		}
	    }
	  A = A*A;
	  B = B*B;

	  I[k] = I[k] + A+B;
	  //	  printf("%.3f\t%.3f\t%.3f\n", I[k], A, B);
	}
      I[k] = I[k]/n;

      I[k] = log10(I[k]);
      //     fprintf(SAXS_PROFILE, "   %.5lf    %.5lf\n", q, I[k]);
    }

  if(printDataFile)
    {
      for(k=0;k<nexp;k++)
	{
	  fprintf(SAXS_PROFILE, "   %.5lf    %.5lf\n", Q[k], I[k]);
	}
      fclose(SAXS_PROFILE);
    }
}
//#############################END OF FUNCTIONS##################################
