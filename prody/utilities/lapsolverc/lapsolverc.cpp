#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cstdint>

#include "dense_wrap.hpp"

namespace py = pybind11;

const char *doc_dense = R"pbdoc(
Min cost bipartite matching via shortest augmenting paths for dense matrices

This is an O(n^3) implementation of a shortest augmenting path
algorithm for finding min cost perfect matchings in dense
graphs.  In practice, it solves 1000x1000 problems in around 1
second.

    rids, cids = solve_dense(costs)
    total_cost = costs[rids, cids].sum()

Params
------
costs : MxN array
    Array containing costs.

Returns
-------
rids : array
    Array of row ids of matching pairings
cids : array
    Array of column ids of matching pairings

)pbdoc";


PYBIND11_MODULE(lapsolverc, m) {
    m.doc() = R"pbdoc(
        Linear assignment problem solvers using native c-extensions.
    )pbdoc";

    // pybind11 first tries each overload (in order) without conversion.
    // If no match is found, it tries again with conversion, unless disallowed.
    // This conversion will cast e.g. double to int.
    // https://pybind11.readthedocs.io/en/stable/advanced/functions.html#overload-resolution-order
    m.def("solve_dense", solve_dense_wrap<int32_t, py::array::c_style>, py::arg().noconvert());
    m.def("solve_dense", solve_dense_wrap<int64_t, py::array::c_style>, py::arg().noconvert());
    m.def("solve_dense", solve_dense_wrap<float, py::array::c_style>, py::arg().noconvert());
    m.def("solve_dense", solve_dense_wrap<double, py::array::c_style>);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}

/////////////////////////////////////////////////////////////////////////////
//
// Python specific setup
//
/////////////////////////////////////////////////////////////////////////////
// static PyMethodDef LAPMethods[] = {
// 	{"ccealign", ccealign_ccealign, METH_VARARGS, "Fast linear assignment problem (LAP) solvers for Python based on c-extensions."},
// 	{NULL, NULL, 0, NULL}
// };

// #if PY_MAJOR_VERSION >= 3

// static struct PyModuleDef lapsolverc = {
//         PyModuleDef_HEAD_INIT,
//         "lapsolverc",
//         "Fast linear assignment problem (LAP) solvers.",
//         -1,
//         LAPMethods,
// };
// PyMODINIT_FUNC PyInit_lapsolverc(void) {
//     return PyModule_Create(&lapsolverc);
// }
// #else
// PyMODINIT_FUNC initlapsolverc(void)
// {
// 	(void) Py_InitModule3("lapsolverc", LAPMethods, "Fast linear assignment problem (LAP) solvers.");
// }
// #endif

