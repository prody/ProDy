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
#define LENLABEL 100


static char *intcat(char *msg, int line) {
   
    /* Concatenate integer to a string. */
   
    char lnum[10];
    snprintf(lnum, 10, "%i", line);
    strcat(msg, lnum);
    return msg;
}
    

static int parseLabel(PyObject *labels, PyObject *mapping, char line[],
                      char clabel[], char ckey[], long ccount, int size) {
    
    /* Parse protein database identifier from sequence label, and map the 
       sequence. */
    
    int i, slash = 0, dash = 0;
    
    for (i = 0; i < size; i++)
        if (line[i] == '\n' || line[i] == ' ') 
            break;
        else if (line[i] == '/' && slash == 0 && dash == 0)
            slash = i;
        else if (line[i] == '-' && slash > 0 && dash == 0)
            dash = i;
    
    PyObject *plabel, *pcount;
    if (slash > 0 && dash > slash) {
        strncpy(ckey, line, slash);
        strncpy(clabel, line, i);
        clabel[i] = '\0';
        ckey[slash] = '\0';
        
        PyObject *pkey = PyString_FromString(ckey);
        plabel = PyString_FromString(clabel);
        pcount = PyInt_FromLong(ccount);
        if (!plabel || !pcount || PyList_Append(labels, plabel) < 0 ||
            PyDict_SetItem(mapping, pkey, pcount)) {
            Py_XDECREF(pcount);
            Py_XDECREF(plabel);
            Py_XDECREF(pkey);
            return 0;
        }
        Py_DECREF(pkey);
        Py_DECREF(plabel);
        Py_DECREF(pcount); 
    } else {
        strncpy(clabel, line, i);
        clabel[i] = '\0';

        plabel = PyString_FromString(clabel);
        pcount = PyInt_FromLong(ccount);
        if (!plabel || !pcount || PyList_Append(labels, plabel) < 0 ||
            PyDict_SetItem(mapping, plabel, pcount)) {
            Py_XDECREF(pcount);
            Py_XDECREF(plabel);
            return 0;
        }
        Py_DECREF(plabel);
        Py_DECREF(pcount);
     }
    return 1;
}


static PyObject *parseFasta(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. have same number of lines at equal
       lengths. */

    char *filename;
    long lenseq, numseq; /* seq length and expected max num of sequences */
    
    if (!PyArg_ParseTuple(args, "sii", &filename, &lenseq, &numseq))
        return NULL;
    
    long i = 0, lenline = 0, lenlast = 0, numlines = 0; 
    long size = lenseq + LENLABEL, iline = 0;
    char errmsg[LENLABEL] = "failed to parse fasta file at line ";

    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping)
        return PyErr_NoMemory();

    char *line = malloc(size * sizeof(char));
    if (!line) 
        return PyErr_NoMemory();
    char *data = malloc(lenseq * numseq * sizeof(char));
    if (!data) {
        free(line);
        return PyErr_NoMemory();
    }

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

    int j = 0;
    long index = 0, ccount = 0;
    char clabel[LENLABEL], ckey[LENLABEL];

    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] != '>')
            continue;
        for (i = 1; i < size; i++)
            if (line[i] != ' ')
                break;
        strcpy(line, line + i);

        /* parse label */
        if (!parseLabel(labels, mapping, line, clabel, ckey, ccount, size)) {
            free(line);
            free(data);
            PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
            return NULL;
        }

        /* parse sequence */
        for (i = 0; i < numlines; i++) {
            if (fgets(line, size, file) == NULL) {
                free(line);
                free(data);
                PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
                return NULL;
            }
            for (j = 0; j < lenline; j++)
                data[index++] = line[j];
        }
        
        if (lenlast) {
            if (fgets(line, size, file) == NULL) {
                free(line);
                free(data);
                PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
                return NULL;
            }
            for (j = 0; j < lenlast; j++)
                data[index++] = line[j];
        }

        ccount++;
    }

    fclose(file);
    free(line);
    data = realloc(data, lenseq * ccount * sizeof(char));
    npy_intp dims[2] = {ccount, lenseq};
    PyObject *msa = PyArray_SimpleNewFromData(2, dims, PyArray_CHAR, data);
    PyObject *result = Py_BuildValue("(OOO)", msa, labels, mapping);
    Py_DECREF(msa);
    Py_DECREF(labels);
    Py_DECREF(mapping);

    return result;
}


static PyObject *writeFasta(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. have same number of lines at equal
       lengths. */

    char *filename;
    PyObject *labels;
    PyArrayObject *msa;
    
    if (!PyArg_ParseTuple(args, "sOOi", &filename, &labels, &msa))
        return NULL;
    
    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
    
    if (numseq != PyList_Size(labels)) {
        PyErr_SetString(PyExc_ValueError,
            "size of labels and msa array does not match");
        return NULL;
    }
    
    FILE *file = fopen(filename, "wb");
    
    printf("%li", lenseq);
    int nlines = lenseq / 60;
    int remainder = lenseq - 60 * nlines;
    long i, j;
    for (i = 0; i <= numseq; i++) {
        
        for (j = 0; j <= nlines; j++) {
        
        }
        
        if (remainder) {
            
        }
        
    }
    
    fclose(file);
    
    return Py_BuildValue("s", filename);
}

static PyObject *parseSelex(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. start and end at the same column. */

    char *filename;
    long lenseq, numseq; /* seq length and expected max num of sequences */
    
    if (!PyArg_ParseTuple(args, "sii", &filename, &lenseq, &numseq))
        return NULL;

    long i = 0, beg = 0, end = 0; 
    long size = lenseq + LENLABEL, iline = 0;
    char errmsg[LENLABEL] = "failed to parse selex/stockholm file at line ";

    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping)
        return PyErr_NoMemory();
    char *line = malloc(size * sizeof(char));
    if (!line)
        return PyErr_NoMemory();
    char *data = malloc(lenseq * numseq * sizeof(char));
    if (!data) {
        free(line);
        return PyErr_NoMemory();
    }

    /* figure out where the sequence starts and ends in a line*/
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

    long index = 0, ccount = 0;
    char clabel[LENLABEL], ckey[LENLABEL];

    int space = beg - 1; /* index of space character before sequence */
    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] == '#' || line[0] == '/' || line[0] == '%')
            continue;
            
        if (line[space] != ' ') {
            free(line);
            free(data);
            PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
            return NULL;
        } 

        /* parse label */
        if (!parseLabel(labels, mapping, line, clabel, ckey, ccount, size)) {
            free(line);
            free(data);
            PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
            return NULL;
        }
        
        /* parse sequence */
        for (i = beg; i < end; i++)
            data[index++] = line[i];
        ccount++;
    }
    fclose(file);
    free(line);
    
    data = realloc(data, lenseq * ccount * sizeof(char));
    npy_intp dims[2] = {ccount, lenseq};
    PyObject *msa = PyArray_SimpleNewFromData(2, dims, PyArray_CHAR, data);
    PyObject *result = Py_BuildValue("(OOO)", msa, labels, mapping);
    Py_DECREF(msa);
    Py_DECREF(labels);
    Py_DECREF(mapping);
    
    return result;
}


static PyMethodDef msaio_methods[] = {

    {"parseFasta",  (PyCFunction)parseFasta, METH_VARARGS, 
     "Return list of labels and a dictionary mapping labels to sequences \n"
     "after parsing the sequences into empty numpy character array."},

    {"writeFasta",  (PyCFunction)writeFasta, METH_VARARGS, 
     "Return filename after writing MSA in FASTA format."},

    {"parseSelex",  (PyCFunction)parseSelex, METH_VARARGS, 
     "Return list of labels and a dictionary mapping labels to sequences \n"
     "after parsing the sequences into empty numpy character array."},

    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initmsaio(void) {

    Py_InitModule3("msaio", msaio_methods,
        "Multiple sequence alignment IO tools.");
        
    import_array();
}
