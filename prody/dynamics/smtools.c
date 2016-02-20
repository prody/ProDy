#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "math.h"
double calcStiffnessMatrixElement(int N, int n_modes, int ind_3i, int ind_3j, double *R_ij_normalized_vec, double *inv_sqrt_eigvals, double *eigvals, double *eigvecs_flat);

static PyObject *calcSM(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject *R_ij_normalized_vec_array, *inv_sqrt_eigvals_array, *eigvals_array, *eigvecs_flat_array, *value_array;
  int numCA, nmodes, ind3i, ind3j;
  double *R_ij_normalized_vec;
  double *inv_sqrt_eigvals;
  double *eigvals;
  double *eigvecs_flat;  
  double *value;
  double m_element;
  
  static char *kwlist[] = {"numCalphas","n_modes","ind_3i","ind_3j",
          "R_ij_sup_0_normalized_vec", "inv_sqrt_eigvals", "eigvals", "eigvecs_flat", "value", NULL};
  /* Parse the input tuple */
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiiiOOOOO", kwlist, 
          &numCA, &nmodes, &ind3i, &ind3j, 
          &R_ij_normalized_vec_array, &inv_sqrt_eigvals_array, &eigvals_array, &eigvecs_flat_array, &value_array))
    return NULL;

  /* Interpret the input objects as numpy arrays. */
  //PyObject *R_ij_normalized_vec_array = PyArray_FROM_OTF(R_ij_normalized_vec_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  //PyObject *inv_sqrt_eigvals_array = PyArray_FROM_OTF(inv_sqrt_eigvals_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  //PyObject *eigvals_array = PyArray_FROM_OTF(eigvals_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  //PyObject *eigvecs_flat_array = PyArray_FROM_OTF(eigvecs_flat_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  // /* If that didn't work, throw an exception. */
  // if (R_ij_normalized_vec_array == NULL || inv_sqrt_eigvals_array == NULL || eigvals_array == NULL || eigvecs_flat_array==NULL) 
  //   {
  //     Py_XDECREF(R_ij_normalized_vec_array);
  //     Py_XDECREF(inv_sqrt_eigvals_array);
  //     Py_XDECREF(eigvals_array);
  //     Py_XDECREF(eigvecs_flat_array);
  //     return NULL;
  //   }

  /* Get pointers to the data as C-types. */
  R_ij_normalized_vec = (double*)PyArray_DATA(R_ij_normalized_vec_array);
  inv_sqrt_eigvals = (double*)PyArray_DATA(inv_sqrt_eigvals_array);
  eigvals = (double*)PyArray_DATA(eigvals_array);
  eigvecs_flat = (double*)PyArray_DATA(eigvecs_flat_array);
  value = (double*)PyArray_DATA(value_array);

  /* Call the external C function to compute the chi-squared. */
  m_element=calcStiffnessMatrixElement(numCA, nmodes, ind3i, ind3j, R_ij_normalized_vec, inv_sqrt_eigvals, eigvals, eigvecs_flat);

  value = &m_element;
  /* Clean up. */
  // Py_DECREF(R_ij_normalized_vec_array);
  // Py_DECREF(inv_sqrt_eigvals_array);
  // Py_DECREF(eigvals_array);
  // Py_DECREF(eigvecs_flat_array);

  //  if (value < 0.0) {
  //    PyErr_SetString(PyExc_RuntimeError,
  //                    "Chi-squared returned an impossible value.");
  //    return NULL;
  //  }

  /* Build the output tuple */
  // PyObject *ret = Py_BuildValue("d", value);
  // return ret;
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

double calcStiffnessMatrixElement(int N, int n_modes, int ind_3i, int ind_3j, double *R_ij_normalized_vec, double *inv_sqrt_eigvals, double *eigvals, double *eigvecs_flat)
{
  int k=0;
  double u_ij_sup_k[3]={0.0, 0.0, 0.0};
  double d_ij_sup_k=0.0;
  double sum1=0.0;
  double sum2=0.0;
  double cos_alpha_ij=0.0;
  for(k=6; k<n_modes; k++)
    {
      //      u_ij_sup_k[0]=(eigvecs[k][ind_3j  ]-eigvecs[k][ind_3i  ]);
      u_ij_sup_k[0]=(eigvecs_flat[k*3*N+ind_3j]-eigvecs_flat[k*3*N+ind_3i  ]);

      //      u_ij_sup_k[1]=(eigvecs[k][ind_3j+1]-eigvecs[k][ind_3i+1]);
      u_ij_sup_k[1]=(eigvecs_flat[k*3*N+ind_3j+1]-eigvecs_flat[k*3*N+ind_3i+1]);

      //      u_ij_sup_k[2]=(eigvecs[k][ind_3j+2]-eigvecs[k][ind_3i+2]);
      u_ij_sup_k[2]=(eigvecs_flat[k*3*N+ind_3j+2]-eigvecs_flat[k*3*N+ind_3i+2]);

      cos_alpha_ij=(  (R_ij_normalized_vec[0]*u_ij_sup_k[0]) +\
          (R_ij_normalized_vec[1]*u_ij_sup_k[1]) +\
          (R_ij_normalized_vec[2]*u_ij_sup_k[2])  );

      d_ij_sup_k=inv_sqrt_eigvals[k]*cos_alpha_ij;
      sum1+=fabs(eigvals[k]*d_ij_sup_k);
      sum2+=fabs(d_ij_sup_k);
    }
  return (sum1/sum2);
}

