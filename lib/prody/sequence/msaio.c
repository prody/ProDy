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
#define FASTALINELEN 1000
#define SELEXLINELEN 10000

static char *intcat(char *msg, int line) {
   
    /* Concatenate integer to a string. */
   
    char lnum[10];
    snprintf(lnum, 10, "%i", line);
    strcat(msg, lnum);
    return msg;
}
    

static int parseLabel(PyObject *labels, PyObject *mapping, char line[], 
                      int length) {
    
    /* Append label to *labels*, extract identifier, and index label 
       position in the list. Return 1 when successful, 0 on failure. */
    
    int i, ch, slash = 0, dash = 0;//, ipipe = 0, pipes[4] = {0, 0, 0, 0};
    
    for (i = 0; i < length; i++) {
        ch = line[i];
        if (ch < 32 && ch != 20)
            break;
        else if (ch == '/' && slash == 0 && dash == 0)
            slash = i;
        else if (ch == '-' && slash > 0 && dash == 0)
            dash = i;
        //else if (line[i] == '|' && ipipe < 4)
        //    pipes[ipipe++] = i;
    }
    
    PyObject *label, *index;
    #if PY_MAJOR_VERSION >= 3
    label = PyUnicode_FromStringAndSize(line, i);
    index = PyLong_FromSsize_t(PyList_Size(labels));
    #else
    label = PyString_FromStringAndSize(line, i);
    index = PyInt_FromSsize_t(PyList_Size(labels));
    #endif

    if (!label || !index || PyList_Append(labels, label) < 0) {
        PyObject *none = Py_None;
        PyList_Append(labels, none);
        Py_DECREF(none);

        Py_XDECREF(index);
        Py_XDECREF(label);
        return 0;
    }
    
    PyObject *key = label;
        
    if (slash > 0 && dash > slash) {
        Py_DECREF(label);
        #if PY_MAJOR_VERSION >= 3
        key = PyUnicode_FromStringAndSize(line, slash);
        #else
        key = PyString_FromStringAndSize(line, slash);
        #endif
    }

    if (PyDict_Contains(mapping, key)) {
        PyObject *item = PyDict_GetItem(mapping, key); /* borrowed */
        if (PyList_Check(item)) {
            PyList_Append(item, index);
            Py_DECREF(index);
        } else {
            PyObject *list = PyList_New(2); /* new reference */
            PyList_SetItem(list, 0, item);
            PyList_SetItem(list, 1, index); /* steals reference, no DECREF */
            PyDict_SetItem(mapping, key, list);
            Py_DECREF(list);
        }
    } else {
        PyDict_SetItem(mapping, key, index);
        Py_DECREF(index);
    }     

    Py_DECREF(key);
    return 1;
}


static PyObject *parseFasta(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. have same number of lines at equal
       lengths. */

    char *filename;
    long filesize;
    int aligned;
    
    if (!PyArg_ParseTuple(args, "sii", &filename, &filesize, &aligned))
        return NULL;
    
    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping)
        return PyErr_NoMemory();
        
    char *line = malloc((FASTALINELEN + 1) * sizeof(char));
    if (!line) 
        return PyErr_NoMemory();
    char *data = malloc(filesize * sizeof(char));
    if (!data) {
        free(line);
        return PyErr_NoMemory();
    }
        
    long iline = 0, i, seqlen = 0, curlen = 0;
    char errmsg[LENLABEL] = "failed to parse FASTA file at line ";

    char ch;
    long index = 0, ccount = -1, clabel = 0;

    FILE *file = fopen(filename, "rb");
    while (fgets(line, FASTALINELEN, file) != NULL) {
        iline++;
        if (line[0] == '>') {
            if (seqlen != curlen) {
                if (seqlen) {
                    free(line);
                    free(data);
                    fclose(file);
                    PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
                    return NULL;
                } else
                    seqlen = curlen;
                ccount++;
            }
            // `line + 1` is to omit `>` character
            clabel += parseLabel(labels, mapping, line + 1, FASTALINELEN);
            curlen = 0;
        } else {
            for (i = 0; i < FASTALINELEN; i++) {
                ch = line[i];
                if (ch < 32)
                    break;
                else {
                    data[index++] = ch;
                    curlen++;
                }
            }
        }
    }
    fclose(file);
    
    if (seqlen != curlen) {
        free(line);
        free(data);
        PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
        return NULL;
    }

    free(line);
    data = realloc(data, index * sizeof(char));
    npy_intp dims[2] = {index / seqlen, seqlen};
    PyObject *msa = PyArray_SimpleNewFromData(2, dims, PyArray_CHAR, data);
    PyObject *result = Py_BuildValue("(OOOi)", msa, labels, mapping, clabel);
    Py_DECREF(msa);
    Py_DECREF(labels);
    Py_DECREF(mapping);
    return result;
}


static PyObject *writeFasta(PyObject *self, PyObject *args, PyObject *kwargs) {

    /* Write MSA where inputs are: labels in the form of Python lists 
    and sequences in the form of Python numpy array and write them in
    FASTA format in the specified filename.*/
    
    char *filename;
    int line_length = 60;
    PyObject *labels;
    PyArrayObject *msa;
    
    static char *kwlist[] = {"filename", "labels", "msa", "line_length", NULL};
    
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sOO|i", kwlist, 
                                     &filename, &labels, &msa, &line_length))
        return NULL;
    
    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa); 

    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
    
    if (numseq != PyList_Size(labels)) {
        PyErr_SetString(PyExc_ValueError,
            "size of labels and msa array does not match");
        return NULL;
    }
    
    FILE *file = fopen(filename, "wb");
    
    int nlines = lenseq / line_length;
    int remainder = lenseq - line_length * nlines;
    int i, j, k;
    int count = 0;
    char *seq = msa->data;
    int lenmsa = strlen(seq);
    #if PY_MAJOR_VERSION >= 3
    PyObject *plabel;
    #endif
    for (i = 0; i < numseq; i++) {
        #if PY_MAJOR_VERSION >= 3
        plabel = PyUnicode_AsEncodedString(
                PyList_GetItem(labels, (Py_ssize_t) i), "utf-8", 
                               "label encoding");
        char *label =  PyBytes_AsString(plabel);
        Py_DECREF(plabel);
        #else
        char *label =  PyString_AsString(PyList_GetItem(labels, 
                                                        (Py_ssize_t) i));
        #endif
        fprintf(file, ">%s\n", label);

        for (j = 0; j < nlines; j++) {
            for (k = 0; k < 60; k++)
                if (count < lenmsa)
                    fprintf(file, "%c", seq[count++]);
            fprintf(file, "\n");
        }
        if (remainder)
            for (k = 0; k < remainder; k++)
                if (count < lenmsa)
                    fprintf(file, "%c", seq[count++]);

        fprintf(file, "\n");
        
    }
    fclose(file);
    return Py_BuildValue("s", filename);
}

static PyObject *parseSelex(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. start and end at the same column. */

    char *filename;
    long filesize;
    int aligned;
    
    if (!PyArg_ParseTuple(args, "sii", &filename, &filesize, &aligned))
        return NULL;

    long i = 0, beg = 0, end = 0;
    long size = SELEXLINELEN + 1, iline = 0, seqlen = 0;
    char errmsg[LENLABEL] = "failed to parse SELEX/Stockholm file at line ";

    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping)
        return PyErr_NoMemory();
    char *line = malloc(size * sizeof(char));
    if (!line)
        return PyErr_NoMemory();
    char *data = malloc(filesize * sizeof(char));
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
        for (; i < size; i++)
            if (line[i] < 32)
                break;
        end = i;
        seqlen = end - beg;
        break;
    }
    iline--;
    fseek(file, - strlen(line), SEEK_CUR);

    long index = 0, ccount = 0, clabel = 0;

    int space = beg - 1; /* index of space character before sequence */
    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] == '#' || line[0] == '/' || line[0] == '%')
            continue;
            
        if (line[space] != ' ') {
            free(line);
            free(data);
            fclose(file);
            PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
            return NULL;
        } 

        clabel += parseLabel(labels, mapping, line, space);
        
        for (i = beg; i < end; i++)
            data[index++] = line[i];
        ccount++;
    }
    fclose(file);
    free(line);

    data = realloc(data, index * sizeof(char));
    npy_intp dims[2] = {index / seqlen, seqlen};
    PyObject *msa = PyArray_SimpleNewFromData(2, dims, PyArray_CHAR, data);
    PyObject *result = Py_BuildValue("(OOOi)", msa, labels, mapping, clabel);
    Py_DECREF(msa);
    Py_DECREF(labels);
    Py_DECREF(mapping);
    
    return result;
}


static PyObject *writeSelex(PyObject *self, PyObject *args, PyObject *kwargs) {
    
    /* Write MSA where inputs are: labels in the form of Python lists 
    and sequences in the form of Python numpy array and write them in
    SELEX (default) or Stockholm format in the specified filename.*/
    
    char *filename;
    PyObject *labels;
    PyArrayObject *msa;
    int stockholm; 
    int label_length = 31;
    
    static char *kwlist[] = {"filename", "labels", "msa", "stockholm", 
                             "label_length", NULL};
    
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sOO|ii", kwlist, &filename,
                                     &labels, &msa, &stockholm, &label_length))
        return NULL;
        
    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa); 

    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
    
    if (numseq != PyList_Size(labels)) {
        PyErr_SetString(PyExc_ValueError,
                        "size of labels and msa array does not match");
        return NULL;
    }
    
    FILE *file = fopen(filename, "wb");
    
    int i, j;
    int pos = 0;
    char *seq = msa->data;
    if (stockholm)
        fprintf(file, "# STOCKHOLM 1.0\n");
    
    char *outline = (char *) malloc((label_length + lenseq + 2) * 
                                    sizeof(char));

    outline[label_length + lenseq] = '\n'; 
    outline[label_length + lenseq + 1] = '\0';
    
    #if PY_MAJOR_VERSION >= 3
    PyObject *plabel;
    #endif
    for (i = 0; i < numseq; i++) {
        #if PY_MAJOR_VERSION >= 3
        plabel = PyUnicode_AsEncodedString(
                PyList_GetItem(labels, (Py_ssize_t) i), "utf-8", 
                               "label encoding");
        char *label =  PyBytes_AsString(plabel);
        Py_DECREF(plabel);
        #else
        char *label = PyString_AsString(PyList_GetItem(labels, (Py_ssize_t)i));
        #endif
        int labelbuffer = label_length - strlen(label);

        strcpy(outline, label);

        if (labelbuffer > 0)
            for(j = strlen(label); j < label_length; j++)
                outline[j] = ' ';
        
        for (j = label_length; j < (lenseq + label_length); j++)
            outline[j] = seq[pos++];

        fprintf(file, "%s", outline);
    }
    
    if (stockholm)
        fprintf(file, "//\n");

    free(outline);
    fclose(file);
    return Py_BuildValue("s", filename);
}


static PyMethodDef msaio_methods[] = {

    {"parseFasta",  (PyCFunction)parseFasta, METH_VARARGS, 
     "Return list of labels and a dictionary mapping labels to sequences \n"
     "after parsing the sequences into empty numpy character array."},

    {"writeFasta",  (PyCFunction)writeFasta, METH_VARARGS | METH_KEYWORDS, 
     "Return filename after writing MSA in FASTA format."},

    {"parseSelex",  (PyCFunction)parseSelex, METH_VARARGS, 
     "Return list of labels and a dictionary mapping labels to sequences \n"
     "after parsing the sequences into empty numpy character array."},

    {"writeSelex",  (PyCFunction)writeSelex, METH_VARARGS | METH_KEYWORDS, 
    "Return filename after writing MSA in SELEX or Stockholm format."},

    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef msaiomodule = {
        PyModuleDef_HEAD_INIT,
        "msaio",
        "Multiple sequence alignment IO tools.",
        -1,
        msaio_methods
};
PyMODINIT_FUNC PyInit_msaio(void) {
    import_array();
    return PyModule_Create(&msaiomodule);
}
#else
PyMODINIT_FUNC initmsaio(void) {

    (void) Py_InitModule3("msaio", msaio_methods,
                          "Multiple sequence alignment IO tools.");
        
    import_array();
}
#endif


