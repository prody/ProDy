/* Author: Tim Lezon, Ahmet Bakan */

#include "Python.h"
#include "numpy/arrayobject.h"


static PyObject *buildhessian(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *coords, *blocks, *hessian, *projection;
    double cutoff = 15., gamma = 1.;
    int natoms, nblocks;

    static char *kwlist[] = {"coords", "blocks", "hessian", "projection",
                             "natoms", "nblocks", "cutoff", "gamma", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOO|iidd", kwlist,
                                     &coords, &blocks, &hessian, &projection,
                                     &natoms, &nblocks, &cutoff, &gamma))
        return NULL;


    double *xyz = (double *) PyArray_DATA(coords);
    long *blk = (long *) PyArray_DATA(blocks);
    double *hess = (double *) PyArray_DATA(hessian);
    double *proj = (double *) PyArray_DATA(projection);


    if (!1)
        return PyErr_NoMemory();


    Py_RETURN_NONE;
}


static PyMethodDef rtbtools_methods[] = {

    {"buildhessian",  (PyCFunction)buildhessian,
     METH_VARARGS | METH_KEYWORDS,
     "Return Hessian matrix."},

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