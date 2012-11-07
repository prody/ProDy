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
    PyObject *result = Py_BuildValue("(OO)", labels, dict);
    Py_DECREF(labels);
    Py_DECREF(dict);
	return result;
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
            if (!plabel || !pcount || 
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
            if (!plabel || !pcount ||
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
    PyObject *result = Py_BuildValue("(OO)", labels, dict);
    Py_DECREF(labels);
    Py_DECREF(dict);
	return result;
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
    Py_RETURN_NONE;

}


static PyObject *calcMutualInfo(PyObject *self, PyObject *args,
                                PyObject *kwargs) {

	PyArrayObject *msa, *mutinfo;
	int debug = 0;
	
    static char *kwlist[] = {"msa", "mutinfo", "debug", NULL};
		
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|i", kwlist,
	                                 &msa, &mutinfo, &debug))
		return NULL;
    /* check dimensions */
    
    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
   
    if (mutinfo->dimensions[0] != lenseq || mutinfo->dimensions[1] != lenseq) {
        PyErr_SetString(PyExc_IOError, 
                        "msa and mutinfo array shapes do not match");
        return NULL;
    }
    
    /* get pointers to data */
    
    char *seq = (char *)PyArray_DATA(msa); /*size: numseq x lenseq */
    double *mut = (double *)PyArray_DATA(mutinfo); /*size: lenseq x lenseq */

    /* start here */
    long i = 0, j = 0, k = 0, l = 0, diff = 0, offset = 0;
    
    int *iseq = malloc((size_t) numseq * sizeof(double));
    if (!iseq) {
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
    }
    
    double **probs;
    double **joint;
    
    probs = (double **) malloc((size_t) lenseq * sizeof(double*));
    if (!probs) {
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
    }
    joint = (double **) malloc((size_t) 27 * sizeof(double*));
    if (!joint) {
        free((void *) probs);
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
    }
    
    for (i = 0; i < lenseq; i++)  {
        probs[i] = (double *) malloc((size_t) 27 * sizeof(double));  
        if (!probs[i]) {
            for (j = 0; j <= i; j++) {
                free((void *) probs[j]);
            }
            free((void *) probs);
            PyErr_SetString(PyExc_MemoryError, "out of memory");
            return NULL;
        }
        for (j = 0; j < 27; j++)
            probs[i][j] = 0;
    }

    for (i = 0; i < 27; i++)  {
        joint[i] = (double *) malloc((size_t) 27 * sizeof(double));  
        if (!joint[i]) {
            for (j = 0; j <= i; j++) {
                free((void *) joint[j]);
            }
            for (j = 0; j <= lenseq; j++) {
                free((void *) probs[j]);
            }
            free((void *) probs);
            free((void *) joint);
            PyErr_SetString(PyExc_MemoryError, "out of memory");
            return NULL;
        }
    }
    
    double *jrow;    
    
    double p_incr = 1. / numseq;
    double p_half = p_incr / 2.;
    double p_twth = p_incr / 20.;
    int twenty[] = {1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 
                    14, 16, 17, 18, 19, 20, 22, 23, 25};

    double mi = 0, jp = 0;
    int a, b;
    /*printf("\n");*/
    for (i = 0; i < lenseq; i++) {
        for (j = i + 1; j < lenseq; j++) {

            for (k = 0; k < 27; k++)
                for (l = 0; l < 27; l++)
                    joint[k][l] = 0;
        
            diff = j - i - 1;
            for (k = 0; k < numseq; k++) {
                offset = k * lenseq;
                if (diff) {
                    a = iseq[k];
                } else {
                    a = (int) seq[offset + i];
                    if (a > 90)
                        a -= 96;
                    else
                        a -= 64;
                    if (a < 1 || a > 26)
                        a = 0; /* gap character */
                    iseq[k] = a;
                }
                
                b = (int) seq[offset + j];
                if (b > 90)
                    b -= 96;
                else
                    b -= 64;
                if (b < 1 || b > 26)
                    b = 0; /* gap character */
                joint[a][b] += p_incr;
                /*
                printf("%li %li %c %c %f \n", 
                       i, j, (char) a + 64, (char) b + 64, joint[a][b]);
                */
                
                if (!i) {
                    if (b == 2) { /* B */
                        probs[j][4] += p_half; /* D */
                        probs[j][14] += p_half; /* N */
                    } else if (b == 26) { /* Z */
                        probs[j][5] += p_half; /* E */
                        probs[j][17] += p_half; /* Q */
                    } else if (b == 10) { /* J */
                        probs[j][9] += p_half; /* I */
                        probs[j][12] += p_half; /* L */
                    } else if (b == 24) { /* X */
                        for (l = 0; l < 20; l++)
                            probs[j][twenty[l]] += p_twth;
                    } else {
                        probs[j][b] += p_incr;
                    }        
                    
                    if (!diff) {
                        if (a == 2) { /* B */
                            probs[i][4] += p_half; /* D */
                            probs[i][14] += p_half; /* N */
                        } else if (a == 26) { /* Z */
                            probs[i][5] += p_half; /* E */
                            probs[i][17] += p_half; /* Q */
                        } else if (a == 10) { /* J */
                            probs[i][9] += p_half; /* I */
                            probs[i][12] += p_half; /* L */
                        } else if (a == 24) { /* X */
                            for (l = 0; l < 20; l++)
                                probs[i][twenty[l]] += p_twth;
                        } else {
                            probs[i][a] += p_incr;
                        }        
                    }
                }
            }
            
            if (debug) {
                printf("\nJoint probability\n");
                double sum = 0;
                for (k = 0; k < 27; k++) {
                    for (l = 0; l < 27; l++) {
                        printf("%.2f ", joint[k][l]*10);
                        sum += joint[k][l];
                    }
                    printf("\n");
                }
                printf("sum %.2f\n", sum);
            }
            
            for (k = 0; k < 27; k++) {
                jrow = joint[k]; 
                /* B */
                jp = jrow[2];  
                if (jp > 0) {
                    jrow[4] = jrow[14] = jp / 2;
                    jrow[2] = 0;
                }
                jp = joint[2][k]; 
                if (jp > 0) {
                    joint[4][k] = joint[14][k] = jp / 2;
                    joint[2][k] = 0;
                }
                /* Z */
                jp = jrow[26]; 
                if (jp > 0) {
                    jrow[5] = jrow[17] = jp / 2;
                    jrow[26] = 0;
                }
                jp = joint[26][k]; 
                if (jp > 0) {
                    joint[5][k] = joint[17][k] = jp / 2;
                    joint[26][k] = 0;
                }
                /* J */
                jp = jrow[10]; 
                if (jp > 0) {
                    jrow[9] = jrow[12] = jp / 2;
                    jrow[10] = 0;
                }
                jp = joint[10][k]; 
                if (jp > 0) {
                    joint[9][k] = joint[12][k] = jp / 2;
                    joint[10][k] = 0;
                }
                /* X */
                jp = jrow[24]; 
                if (jp > 0) {
                    jp = jp / 20.;
                    for (l = 0; l < 20; l++)
                        jrow[twenty[l]] = jp;    
                    jrow[24] = 0;
                }
                jp = joint[24][k]; 
                if (jp > 0) {
                    jp = jp / 20.;
                    for (l = 0; l < 20; l++)
                        joint[twenty[l]][k] = jp;    
                    joint[24][l] = 0;
                }
            }
            
            if (debug) {
                printf("\nJoint probability\n");
                double sum = 0;
                for (k = 0; k < 27; k++) {
                    for (l = 0; l < 27; l++) {
                        printf("%.2f ", joint[k][l]*10);
                        sum += joint[k][l];
                    }
                    printf("\n");
                }
                printf("sum %.2f\n", sum);
            }
            
            mi = 0;
            for (k = 0; k < 27; k++) {
                jrow = joint[k];
                for (l = 0; l < 27; l++) {
                    jp = jrow[l];
                    if (jp > 0) {
                        if (debug)
                            printf("%c (%.3f) %c (%.3f) - (%.3f)\n",
                                    (char) k + 64, probs[i][k],  
                                    (char) l + 64, probs[j][l], jp);
                                
                        mi += jp * log(jp / probs[i][k] / probs[j][l]);
                    }
                }
            }        
            mut[i * lenseq + j] = mi;
            mut[i + lenseq * j] = mi;
            /*printf("%li %li %f\n", i, j, mi);*/
        }
    }
    if (debug) {
        printf("Probability table\n");
        for (i=0; i<lenseq; i++) { 
            for (j=0; j<27; j++) {
                printf("%.2f ", probs[i][j]);  
            }
            printf("\n");
        }
    }
    for (i = 0; i < lenseq; i++){  
        free((void *) probs[i]);
    }  
    free((void *) probs);
    
    for (i = 0; i < 27; i++){  
        free((void *) joint[i]);
    }  
    free((void *) joint);
    free((void *) iseq);
    /* end here */
    Py_RETURN_NONE;
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
