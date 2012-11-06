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


static PyObject *parseFasta(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. have same number of lines at equal
       lengths. */

	char *filename;
	PyObject *arrobj;
	PyArrayObject *msa;
	
	if (!PyArg_ParseTuple(args, "sO", &filename, &arrobj))
		return NULL;
	
    msa = (PyArrayObject *) 
        PyArray_ContiguousFromObject(arrobj, PyArray_CHAR, 2, 2);
    if (msa == NULL)
        return NULL;

    long i = 0, lenseq = msa->dimensions[1];
    long lenline = 0, lenlast = 0, numlines = 0; 
    long size = lenseq + 100, iline = 0;
    char line[size];

    FILE *file = fopen(filename, "rb");
    while (fgets(line, size, file) != NULL) {
        if (line[0] == '>')
            continue;

        for (i = 0; i < strlen(line); i++)
            if (line[i] == ' ' || line[i] == '\n')
                break;
        lenline = i;
        lenlast = lenseq % lenline;
        numlines = (lenseq - lenlast) / lenline;
        break;
    }
    
    
    fseek(file, 0, SEEK_SET);


    int slash = 0, dash = 0, j = 0;
    long index = 0, ccount = 0;
    char *data = (char *)PyArray_DATA(msa);
    char clabel[size], ckey[size];
	PyObject *labels, *dict, *plabel, *pkey, *pcount;
	labels = PyList_New(0);
    dict = PyDict_New();

    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] != '>')
            continue;
        for (i = 1; i < size; i++)
            if (line[i] != ' ')
                break;
        strcpy(line, line + i);

        /* parse label */
        slash = 0;
        dash = 0;
        for (i = 0; i < size; i++)
            if (line[i] == '\n' || line[i] == ' ') 
                break;
            else if (line[i] == '/' && slash == 0 &&  dash == 0)
                slash = i;
            else if (line[i] == '-' && slash > 0 && dash == 0)
                dash = i;
        
        if (slash > 0 && dash > slash) {
            strncpy(ckey, line, slash);
            strncpy(clabel, line, i);
            
            clabel[i] = '\0';
            ckey[slash] = '\0';
            pkey = PyString_FromString(ckey);
            plabel = PyString_FromString(clabel);
            pcount = PyInt_FromLong(ccount);
            if (plabel == NULL || pcount == NULL ||
                PyList_Append(labels, plabel) < 0 ||
                PyDict_SetItem(dict, pkey, pcount)) {
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                Py_DECREF(arrobj);
                Py_XDECREF(pcount);
                Py_XDECREF(plabel);
                Py_XDECREF(pkey);
                return NULL;
            }
            Py_DECREF(pkey);
            Py_DECREF(plabel);
            Py_DECREF(pcount); 
        } else {
            strncpy(clabel, line, i);
            clabel[i] = '\0';
            plabel = PyString_FromString(clabel);
            pcount = PyInt_FromLong(ccount);
            if (plabel == NULL || pcount == NULL ||
                PyList_Append(labels, plabel) < 0 ||
                PyDict_SetItem(dict, plabel, pcount)) {
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                Py_DECREF(arrobj);
                Py_XDECREF(pcount);
                Py_XDECREF(plabel);
                return NULL;
            }
            Py_DECREF(plabel);
            Py_DECREF(pcount);
         }

        
        /* parse sequence */
        
        for (i = 0; i < numlines; i++) {
            if (fgets(line, size, file) == NULL) {
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                Py_DECREF(arrobj);
                return NULL;
            }
            for (j = 0; j < lenline; j++)
                data[index++] = line[j];
        }
        
        if (lenlast) {
            if (fgets(line, size, file) == NULL) {
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                Py_DECREF(arrobj);
                return NULL;
            }
            for (j = 0; j < lenlast; j++)
                data[index++] = line[j];
        }
        ccount++;
    }

    fclose(file);
    Py_XDECREF(arrobj);
	return Py_BuildValue("(OO)", labels, dict);
}


static PyObject *parseSelex(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. start and end at the same column. */

	char *filename;
	PyObject *arrobj;
	PyArrayObject *msa;
	
	if (!PyArg_ParseTuple(args, "sO", &filename, &arrobj))
		return NULL;

    msa = (PyArrayObject *) 
        PyArray_ContiguousFromObject(arrobj, PyArray_CHAR, 2, 2);
    if (msa == NULL)
        return NULL;

    long i = 0, beg = 0, end = 0, lenseq = msa->dimensions[1]; 
    long size = lenseq + 100, iline = 0;
    char line[size];

    FILE *file = fopen(filename, "rb");
    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] == '#' || line[0] == '/' || line[0] == '%')
            continue;

        for (i = 0; i < size; i++)
            if (line[i] == ' ')
                break;
        for (; i < size; i++)
            if (line[i] != ' ')
                break;
        beg = i;
        end = beg + lenseq;
        break;
    }
    iline--;
    fseek(file, - strlen(line), SEEK_CUR);

    int slash = 0, dash = 0;
    long index = 0, ccount = 0;
    char *data = (char *)PyArray_DATA(msa);
    char clabel[beg], ckey[beg];
	PyObject *labels, *dict, *plabel, *pkey, *pcount;
	labels = PyList_New(0);
    dict = PyDict_New();

    int space = beg - 1; /* index of space character before sequence */
    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] == '#' || line[0] == '/' || line[0] == '%')
            continue;
            
        if (line[space] != ' ') {
            PyErr_SetString(PyExc_IOError, 
                            "failed to parse msa, at line");
            return NULL;
        } 

        /* parse label */
        
        slash = 0;
        dash = 0;
        for (i = 0; i < size; i++)
            if (line[i] == ' ') 
                break;
            else if (line[i] == '/' && slash == 0 &&  dash == 0)
                slash = i;
            else if (line[i] == '-' && slash > 0 && dash == 0)
                dash = i;
        if (slash > 0 && dash > slash) {
            strncpy(ckey, line, slash);
            strncpy(clabel, line, i);
            clabel[i] = '\0';
            ckey[slash] = '\0';
            pkey = PyString_FromString(ckey);
            plabel = PyString_FromString(clabel);
            pcount = PyInt_FromLong(ccount);
            if (plabel == NULL || pcount == NULL ||
                PyList_Append(labels, plabel) < 0 ||
                PyDict_SetItem(dict, pkey, pcount)) {
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                Py_DECREF(arrobj);
                Py_XDECREF(pcount);
                Py_XDECREF(plabel);
                Py_XDECREF(pkey);
                return NULL;
            }
            Py_DECREF(pkey);
            Py_DECREF(plabel);
            Py_DECREF(pcount);            
        } else {
            strncpy(clabel, line, i);
            clabel[i] = '\0';
            plabel = PyString_FromString(clabel);
            pcount = PyInt_FromLong(ccount);
            if (plabel == NULL || pcount == NULL ||
                PyList_Append(labels, plabel) < 0 ||
                PyDict_SetItem(dict, plabel, pcount)) {
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                Py_DECREF(arrobj);
                Py_XDECREF(pcount);
                Py_XDECREF(plabel);
                return NULL;
            }
            Py_DECREF(plabel);
            Py_DECREF(pcount);
         }
        
        /* parse sequence */
        for (i = beg; i < end; i++)
            data[index++] = line[i];
        ccount++;
    }
    fclose(file);
    Py_XDECREF(arrobj);
	return Py_BuildValue("(OO)", labels, dict);
}


static PyObject *calcShannonEntropy(PyObject *self, PyObject *args,
                                    PyObject *kwargs) {

	PyObject *arrobj, *result;
	PyArrayObject *msa, *entropy;
	int dividend = 0;
	
    static char *kwlist[] = {"msa", "entropy", "dividend", NULL};
		
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|i", kwlist,
	                                 &arrobj, &result, &dividend))
		return NULL;
    
    msa = (PyArrayObject *) 
        PyArray_ContiguousFromObject(arrobj, PyArray_CHAR, 2, 2);
    if (msa == NULL)
        return NULL;
    
    entropy = (PyArrayObject *) 
        PyArray_ContiguousFromObject(result, PyArray_DOUBLE, 1, 1);
    if (entropy == NULL)
        return NULL;

    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
   
    if (entropy->dimensions[0] != lenseq) {
        Py_XDECREF(arrobj);
        Py_XDECREF(result);
        PyErr_SetString(PyExc_IOError, 
                        "msa and entropy array shapes do not match");
        return NULL;
    }

    char *seq = (char *)PyArray_DATA(msa);
    double *ent = (double *)PyArray_DATA(entropy);

    /* start here */
    long size = numseq * lenseq; 
    double count[128]; /* number of ASCII characters*/
    double shannon = 0, probability = 0, numgap = 0, denom = numseq;
    long i = 0, j = 0;
    
    double ambiguous = 0;
    int twenty[20] = {65, 67, 68, 69, 70, 71, 72, 73, 75, 76, 
                      77, 78, 80, 81, 82, 83, 84, 86, 87, 89};
    for (i = 0; i < lenseq; i++) {

        /* zero counters */
        for (j = 65; j < 91; j++)
            count[j] = 0;
        for (j = 97; j < 123; j++)
            count[j] = 0;
        
        /* count characters in a column*/
        for (j = i; j < size; j += lenseq)
            count[(int) seq[j]]++;
        for (j = 65; j < 91; j++)
            count[j] += count[j + 32];
        
        /* handle ambiguous amino acids */
        if (count[66]) {
            ambiguous = count[66] / 2.; /* B */
            count[66] = 0;
            count[68] += ambiguous; /* D */
            count[78] += ambiguous; /* N */
        }
        if (count[90]) {
            ambiguous = count[90] / 2.; /* Z */
            count[90] = 0;
            count[69] += ambiguous; /* E */
            count[81] += ambiguous; /* Q */
        }
        if (count[74]) {
            ambiguous = count[74] / 2.; /* J */
            count[74] = 0;
            count[73] += ambiguous; /* I */
            count[76] += ambiguous; /* L */
        }
        if (count[88]) {
            ambiguous = count[88] / 20.; /* X */
            count[88] = 0;
            for (j = 0; j < 20; j++)
                count[twenty[j]] += ambiguous;
        }
        
        /* non-gap counts */
        numgap = numseq;
        for (j = 65; j < 91; j++)
            numgap -= count[j];
        
        shannon = 0;
        denom = numseq;
        if (dividend)
            denom = numseq - numgap;
        else if (numgap > 0) {
            probability = numgap / numseq;
            shannon += probability * log(probability);
        }

        for (j = 65; j < 91; j++) {
            if (count[j] > 0) {
                probability = count[j] / denom;
                shannon += probability * log(probability);
            }
        }
        ent[i] = -shannon;
    }
    /* end here */
    Py_XDECREF(arrobj);
    Py_XDECREF(result);
    return Py_BuildValue("");

}


static PyObject *calcMutualInfo(PyObject *self, PyObject *args,
                                PyObject *kwargs) {

	PyObject *arrobj, *result;
	PyArrayObject *msa, *mutinfo;
	
    static char *kwlist[] = {"msa", "mutinfo", NULL};
		
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwlist,
	                                 &arrobj, &result))
		return NULL;
    
    msa = (PyArrayObject *) 
        PyArray_ContiguousFromObject(arrobj, PyArray_CHAR, 2, 2);
    if (msa == NULL)
        return NULL;
    
    mutinfo = (PyArrayObject *) 
        PyArray_ContiguousFromObject(result, PyArray_DOUBLE, 1, 1);
    if (mutinfo == NULL)
        return NULL;

    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
   
    if (mutinfo->dimensions[0] != lenseq || mutinfo->dimensions[1] != lenseq) {
        Py_XDECREF(arrobj);
        Py_XDECREF(result);
        PyErr_SetString(PyExc_IOError, 
                        "msa and mutinfo array shapes do not match");
        return NULL;
    }

    char *seq = (char *)PyArray_DATA(msa); /*size: numseq x lenseq */
    double *mut = (double *)PyArray_DATA(mutinfo); /*size: lenseq x lenseq */

    /* start here */
    long i = 0, j = 0;
    
    double **count;
    count = (double **) malloc((size_t) numseq * sizeof(double*));
    if (!count) {
    free((void *) count);
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return 0;
    }


    for (i = 0; i < numseq; i++)  {
        count[i] = (double *) malloc((size_t) 27 * sizeof(double));  
        if (!count[i]) {
            for (j = 0; j <= i; j++) {
                free((void *) count[j]);
            }
            free((void *) count);
            PyErr_SetString(PyExc_MemoryError, "out of memory");
            return NULL;
        }
        for (j = 0; j <= 27; j++)
            count[i][0] = 0;
    } 
    double joint[90][90];


    for (i = 0; i < lenseq; i++) {
        for (j = i + 1; j < lenseq; j++) {
            
        }
    }

    for (i = 0; i < numseq; i++){  
        free((void *) count[i]);
    }  
    free((void *) count);

    /* end here */
    Py_XDECREF(arrobj);
    Py_XDECREF(result);
    return Py_BuildValue("");

}


static PyMethodDef msatools_methods[] = {

	{"parseSelex",  (PyCFunction)parseSelex, METH_VARARGS, 
	 "Return list of labels and a dictionary mapping labels to sequences \n"
	 "after parsing the sequences into empty numpy character array."},

	{"parseFasta",  (PyCFunction)parseFasta, METH_VARARGS, 
	 "Return list of labels and a dictionary mapping labels to sequences \n"
	 "after parsing the sequences into empty numpy character array."},

	{"calcShannonEntropy",  (PyCFunction)calcShannonEntropy, 
     METH_VARARGS | METH_KEYWORDS, 
	 "Calculate information entropy for given character array into given \n"
     "double array."},

	{"calcMutualInfo",  (PyCFunction)calcMutualInfo, 
     METH_VARARGS | METH_KEYWORDS, 
	 "Calculate mutual information for given character array into given \n"
     "2D double array."},
     
	{NULL, NULL, 0, NULL}
	
};


PyMODINIT_FUNC initmsatools(void) {

	Py_InitModule3("msatools", msatools_methods,
	    "Multiple sequence alignment IO and analysis tools.");
	    
    import_array();
}
