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

    FILE *file = fopen(filename, "r");
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

    FILE *file = fopen(filename, "r");
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


static PyMethodDef msatools_methods[] = {

	{"parseSelex",  (PyCFunction)parseSelex, METH_VARARGS | METH_KEYWORDS, 
	 "Return list of labels and a dictionary mapping labels to sequences \n"
	 "after parsing the sequences into empty numpy character array."},

	{"parseFasta",  (PyCFunction)parseFasta, METH_VARARGS | METH_KEYWORDS, 
	 "Return list of labels and a dictionary mapping labels to sequences \n"
	 "after parsing the sequences into empty numpy character array."},

	{NULL, NULL, 0, NULL}
	
};


PyMODINIT_FUNC initmsatools(void) {

	Py_InitModule3("msatools", msatools_methods,
	    "Multiple sequence alignment tools.");
	    
    import_array();
}
