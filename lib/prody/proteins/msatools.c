#include "Python.h"
#include "numpy/arrayobject.h"

static PyObject *parseSelex(PyObject *self, PyObject *args) {

	char *filename;
	PyObject *arrobj;
	PyArrayObject *msa;
	
	if (!PyArg_ParseTuple(args, "sO", &filename, &arrobj))
		return NULL;

    msa = (PyArrayObject *) 
        PyArray_ContiguousFromObject(arrobj, PyArray_CHAR, 2, 2);
    if (msa == NULL)
        return NULL;

    FILE *file = fopen(filename, "r");
    
    long i = 0, beg = 0, end = 0, lenseq = msa->dimensions[1]; 
    long size = lenseq + 100;
 
    char line[size];
    
    char *data = (char *)PyArray_DATA(msa);
    
    char title_chars[40];
	PyObject *titles;
	titles = PyList_New(0);
	PyObject *title;

    long index = 0, iline = 0;
    while (fgets(line, size, file) != NULL) {
        iline++;
        
        if (line[0] == '#' || line[0] == '/' || line[0] == '%')
            continue;
            
        /* parse title */
        for (i = 0; i < size; i++)
            if (line[i] == ' ')
                break;
        strncpy(title_chars, line, i);
        title_chars[i] = '\0';
        title = PyString_FromString(title_chars);
        if (title == NULL || PyList_Append(titles, title) < 0) {
            PyErr_SetString(PyExc_IOError, 
                            "failed to parse msa, at line");
            return NULL;
        }
        Py_DECREF(title);

        for (; i < size; i++)
            if (line[i] != ' ')
                break;
        beg = i;
        end = beg + lenseq;
        
        /* parse sequence */
        for (i = beg; i < end; i++) {
            data[index] = line[i];
            index++;
        }
        break;
    }
    
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

        /* parse title */
        for (i = 0; i < size; i++)
            if (line[i] == ' ')
                break;
        strncpy(title_chars, line, i);
        title_chars[i] = '\0';
        title = PyString_FromString(title_chars);
        if (title == NULL || PyList_Append(titles, title) < 0) {
            PyErr_SetString(PyExc_IOError, 
                            "failed to parse msa, at line");
            return NULL;
        }
        Py_DECREF(title);
        
        /* parse sequence */
        for (i = beg; i < end; i++) {
            data[index] = line[i];
            index++;
        }
    }
    fclose(file);
    Py_XDECREF(arrobj);
	return titles;
}


static PyMethodDef msatools_methods[] = {
	{"parseSelex",  (PyCFunction)parseSelex, METH_VARARGS | METH_KEYWORDS, 
	 "Return list of titles, and fill in the sequences into empty numpy \n"
	 "character array, which may have more rows than number of sequences."},
	{NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initmsatools(void) {
	Py_InitModule3("msatools", msatools_methods,
	    "Multiple sequence alignment tools.");
    import_array();
}
