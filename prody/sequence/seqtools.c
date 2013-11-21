/* ProDy: A Python Package for Protein Dynamics Analysis
 *
 * Copyright (C) 2010-2012 Ahmet Bakan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * Author: Ahmet Bakan
 * Copyright (C) 2010-2012 Ahmet Bakan
 */

#include "Python.h"
#include "numpy/arrayobject.h"
#define NUMCHARS 27
const int twenty[20] = {1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13,
                        14, 16, 17, 18, 19, 20, 22, 23, 25};
const int unambiguous[23] = {0, 1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14,
                             15, 16, 17, 18, 19, 20, 21, 22, 23, 25};

static PyObject *msaeye(PyObject *self, PyObject *args,
                                   PyObject *kwargs) {

    PyArrayObject *msa, *array;
    double unique = 0;
    int turbo = 1;

    static char *kwlist[] = {"msa", "array", "unique", "turbo", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|di", kwlist,
                                     &msa, &array, &unique, &turbo))
        return NULL;

    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);

    /* get dimensions */
    long number = msa->dimensions[0], length = msa->dimensions[1];

    /* get pointers to data */
    char *iraw, *jraw, *raw = (char *) PyArray_DATA(msa);

    long i, j;
    /* allocate memory */
    double *jrow, *sim =  (double *) raw;
    _Bool *unq = (_Bool *) raw; /* to avoid uninitialized warnings*/
    if (unique)
        unq = (_Bool *) PyArray_DATA(array);
    else
        sim = (double *) PyArray_DATA(array);

    /* arrays to store refined sequences*/
    unsigned char *iseq = malloc(length * sizeof(unsigned char));
    if (!iseq)
        return PyErr_NoMemory();

    unsigned char **seq = malloc(number * sizeof(unsigned char *));
    if (!seq) {
        turbo = 0;
    }

    if (turbo) {
        /* allocate rows that will store columns of MSA */
        seq[0] = iseq;
        for (i = 1; i < number; i++) {
            seq[i] = malloc(length * sizeof(unsigned char));
            if (!seq[i]) {
                for (j = 1; j < i; j++)
                    free(seq[j]);
                free(seq);
                turbo = 0;
            }
        }
    }

    /* initialize jseq, so that we don't get uninitialized warning */
    unsigned char *jseq = iseq;

    unsigned char a, b;
    long k, diff;

    /* zero sim array */
    if (unique) {
        for (i = 0; i < number; i++)
            unq[i] = 1;
    } else {
        for (i = 0; i < number; i++) {
            jrow = sim + i * number;
            for (j = 0; j < number; j++)
                jrow[j] = 0;
            jrow[i] = 1;
        }
    }

    double ncols, score, seqid;

    /* START calculation */
    /* calculate first row of MI matrix and all column probabilities */
    i = 0;
    iraw = raw;
    for (j = 1; j < number; j++) {
        ncols = score = 0.;
        jraw = raw + length * j;
        diff = j - 1;
        if (turbo) /* in turbo mode, there is a row for refined sequences */
            jseq = seq[j];
        for (k = 0; k < length; k++) {
            if (diff) {
                a = iseq[k];
            } else {
                a = (unsigned char) iraw[k];
                if (a > 90)
                    a -= 96;
                else
                    a -= 64;
                if (a < 1 || a > 26)
                    a = 0; /* gap character */
                iseq[k] = a;
            }

            b = (unsigned char) jraw[k];
            if (b > 90)
                b -= 96;
            else
                b -= 64;
            if (b < 1 || b > 26)
                b = 0; /* gap character */
            if (turbo)  /* we keep the refined chars for all sequences*/
                jseq[k] = b;

            if (a || b) {
                ncols++;
                if (a == b)
                    score++;
            }
        }

        seqid = score / ncols;
        if (unique) {
            if (seqid >= unique)
                unq[j] = 0;
        } else
            if (ncols)
                sim[j] = sim[number * j] = seqid;
    }

    if (turbo)
        free(iseq);

    /* calculate rest of identities */
    for (i = 1; i < number; i++) {

        if (unique && !unq[i])
            continue;

        if (turbo)
            iseq = seq[i];
        else
            iraw = raw + length * i;

        for (j = i + 1; j < number; j++) {
            ncols = score = 0.;

            if (turbo) {
                jseq = seq[j];
                for (k = 0; k < length; k++) {
                    a = iseq[k];
                    b = jseq[k];
                    if (a || b) {
                        ncols++;
                        if (a == b)
                            score++;
                    }
                }
            } else {
                jraw = raw + length * j;
                diff = j - i - 1;
                for (k = 0; k < length; k++) {
                    if (diff) {
                        a = iseq[k];
                    } else {
                        a = (unsigned char) iraw[k];
                        if (a > 90)
                            a -= 96;
                        else
                            a -= 64;
                        if (a < 1 || a > 26)
                            a = 0; /* gap character */
                        iseq[k] = a;
                    }

                    b = (unsigned char) jraw[k];
                    if (b > 90)
                        b -= 96;
                    else
                        b -= 64;
                    if (b < 1 || b > 26)
                        b = 0; /* gap character */

                    if (a || b) {
                        ncols++;
                        if (a == b)
                            score++;
                    }
                }
            }

            seqid = score / ncols;
            if (unique) {
                if (seqid >= unique)
                    unq[j] = 0;
            } else
                if (ncols)
                    sim[i * number + j] = sim[i + number * j] = seqid;
         }
    }

    /* free memory */
    if (turbo)
        for (j = 1; j < number; j++)
            free(seq[j]);
    free(seq);

    return Py_BuildValue("O", array);
}

static PyMethodDef seqtools_methods[] = {

    {"msaeye",  (PyCFunction)msaeye,
     METH_VARARGS | METH_KEYWORDS,
     "Return sequence identity matrix calculated for given character \n"
     "array that contains an MSA."},

    {NULL, NULL, 0, NULL}
};



#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef seqtools = {
        PyModuleDef_HEAD_INIT,
        "seqtools",
        "Sequence similarity/identity analysis tools.",
        -1,
        seqtools_methods,
};
PyMODINIT_FUNC PyInit_seqtools(void) {
    import_array();
    return PyModule_Create(&seqtools);
}
#else
PyMODINIT_FUNC initseqtools(void) {

    Py_InitModule3("seqtools", seqtools_methods,
        "Sequence similarity/identity analysis tools.");

    import_array();
}
#endif
