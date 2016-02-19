#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>
#include "smtools.h"
static char module_docstring[] =
  "This module provides an interface for calculating elements of stiffness matrix in C.";
static char calcStiffnessMatrixElement_docstring[] =
  "Calculate elements of stiffness matrix given a mode.";

static PyObject *calcStiffnessMatrixElement_calcStiffnessMatrixElement(PyObject *self, PyObject *args);


static PyMethodDef module_methods[] = {
  {"calcStiffnessMatrixElement", calcStiffnessMatrixElement_calcStiffnessMatrixElement, METH_VARARGS, calcStiffnessMatrixElement_docstring},
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC init_calcStiffnessMatrixElement(void)
{
  PyObject *m = Py_InitModule3("_calcStiffnessMatrixElement", module_methods, module_docstring);
  if (m == NULL)
    return;

  /* Load `numpy` functionality. */
  import_array();
}

static PyObject *calcStiffnessMatrixElement_calcStiffnessMatrixElement(PyObject *self, PyObject *args)
{
  int N, n_modes, ind_3i, ind_3j;
  PyObject *R_ij_normalized_vec_obj, *inv_sqrt_eigvals_obj, *eigvals_obj, *eigvecs_flat_obj;

  /* Parse the input tuple */
  if (!PyArg_ParseTuple(args, "iiiiOOOO", &N, &n_modes, &ind_3i, &ind_3j, &R_ij_normalized_vec_obj, &inv_sqrt_eigvals_obj, &eigvals_obj, &eigvecs_flat_obj))
    return NULL;

  /* Interpret the input objects as numpy arrays. */
  PyObject *R_ij_normalized_vec_array = PyArray_FROM_OTF(R_ij_normalized_vec_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *inv_sqrt_eigvals_array = PyArray_FROM_OTF(inv_sqrt_eigvals_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *eigvals_array = PyArray_FROM_OTF(eigvals_obj, NPY_DOUBLE, NPY_IN_ARRAY);
  PyObject *eigvecs_flat_array = PyArray_FROM_OTF(eigvecs_flat_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  /* If that didn't work, throw an exception. */
  if (R_ij_normalized_vec_array == NULL || inv_sqrt_eigvals_array == NULL || eigvals_array == NULL || eigvecs_flat_array==NULL) 
    {
      Py_XDECREF(R_ij_normalized_vec_array);
      Py_XDECREF(inv_sqrt_eigvals_array);
      Py_XDECREF(eigvals_array);
      Py_XDECREF(eigvecs_flat_array);
      return NULL;
    }

  /* Get pointers to the data as C-types. */
  double *R_ij_normalized_vec    = (double*)PyArray_DATA(R_ij_normalized_vec_array);
  double *inv_sqrt_eigvals       = (double*)PyArray_DATA(inv_sqrt_eigvals_array);
  double *eigvals                = (double*)PyArray_DATA(eigvals_array);
  double *eigvecs_flat                = (double*)PyArray_DATA(eigvecs_flat_array);

  /* Call the external C function to compute the chi-squared. */
  double value=calcStiffnessMatrixElement(N, n_modes, ind_3i, ind_3j, R_ij_normalized_vec, inv_sqrt_eigvals, eigvals, eigvecs_flat);

  /* Clean up. */
  Py_DECREF(R_ij_normalized_vec_array);
  Py_DECREF(inv_sqrt_eigvals_array);
  Py_DECREF(eigvals_array);
  Py_DECREF(eigvecs_flat_array);

  //  if (value < 0.0) {
  //    PyErr_SetString(PyExc_RuntimeError,
  //                    "Chi-squared returned an impossible value.");
  //    return NULL;
  //  }

  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("d", value);
  return ret;
}
