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

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static char *intcat(char *msg, int line) {
   
    /* Concatenate integer to a string. */
   
    char lnum[10];
    snprintf(lnum, 10, "%i", line);
    strcat(msg, lnum);
    return msg;
}
    

static int parseLabel(PyObject *labels, PyObject *mapping, char line[],
                      long ccount, int size) {
    
    /* Parse protein database identifier from sequence label, and map the 
       sequence. */
    
    int i, slash = 0, dash = 0, pipes[4] = {0, 0, 0, 0};
    
    for (i = 0; i < size; i++)
        if (line[i] == '\n' || line[i] == ' ') 
            break;
        else if (line[i] == '/' && slash == 0 && dash == 0)
            slash = i;
        else if (line[i] == '-' && slash > 0 && dash == 0)
            dash = i;
        else if (line[i] == '|') {
            pipes[0]++;
            pipes[pipes[0]] = i;
        }
    
    PyObject *plabel, *pcount;
    if (slash > 0 && dash > slash) {
        #if PY_MAJOR_VERSION >= 3
        PyObject *pkey = PyUnicode_FromStringAndSize(line, slash);
        plabel = PyUnicode_FromStringAndSize(line, i);
        pcount = PyLong_FromLong(ccount);
        #else
        PyObject *pkey = PyString_FromStringAndSize(line, slash);
        plabel = PyString_FromStringAndSize(line, i);
        pcount = PyInt_FromLong(ccount);
        #endif
        if (!plabel || !pcount || PyList_Append(labels, plabel) < 0 ||
            PyDict_SetItem(mapping, pkey, pcount)) {
            Py_XDECREF(pcount);
            Py_XDECREF(plabel);
            Py_XDECREF(pkey);
            return 0;
        }
        Py_DECREF(pkey);
    } else {
        #if PY_MAJOR_VERSION >= 3
        plabel = PyUnicode_FromStringAndSize(line, i);
        pcount = PyLong_FromLong(ccount);
        #else
        plabel = PyString_FromStringAndSize(line, i);
        pcount = PyInt_FromLong(ccount);
        #endif
        if (!plabel || !pcount || PyList_Append(labels, plabel) < 0 ||
            PyDict_SetItem(mapping, plabel, pcount)) {
            Py_XDECREF(pcount);
            Py_XDECREF(plabel);
            return 0;
        }
     }
    Py_DECREF(plabel);
    Py_DECREF(pcount);
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

    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] != '>')
            continue;
        for (i = 1; i < size; i++)
            if (line[i] != ' ')
                break;

        /* parse label */
        if (!parseLabel(labels, mapping, line + i, ccount, size)) {
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
        if (!parseLabel(labels, mapping, line, ccount, size)) {
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
    
    for (i = 0; i < numseq; i++) {
        #if PY_MAJOR_VERSION >= 3
        char *label = PyBytes_AsString(PyList_GetItem(labels, (Py_ssize_t)i));
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


/*PyMODINIT_FUNC initmsaio(void) {

    (void) Py_InitModule3("msaio", msaio_methods,
                          "Multiple sequence alignment IO tools.");
        
    import_array();
}*/


#if PY_MAJOR_VERSION >= 3

static int msaio_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int msaio_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "msaio",
        NULL,
        sizeof(struct module_state),
        msaio_methods,
        NULL,
        msaio_traverse,
        msaio_clear,
        NULL
};

#define INITERROR return NULL

PyObject *
PyInit_msaio(void)

#else
#define INITERROR return

void
initmsaio(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule3("msaio", msaio_methods,
                          "Multiple sequence alignment IO tools.");
#endif

    import_array();
    
    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("myextension.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

