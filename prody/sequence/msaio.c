#include "Python.h"
#include <ctype.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#define LENLABEL 100
#define FASTALINELEN 10000
#define SELEXLINELEN 10000

static char *intcat(char *msg, int line) {

    /* Concatenate integer to a string. */

    char lnum[10];
    snprintf(lnum, 10, "%i", line);
    strcat(msg, lnum);
    return msg;
}


static int labelLength(char *line, int length, int tabs) {

    /* Return the length of the label at the start of *line*, which ends at
       the first control character, or the first one other than a tab when
       *tabs* is true.  */

    int i, ch;
    for (i = 0; i < length; i++) {
        ch = (unsigned char) line[i];
        if (ch < 32 && !(tabs && ch == '\t'))
            break;
    }
    return i;
}


static PyObject *decodeLabel(char *line, int length) {

    /* Return a new reference to *line* decoded as a label.  */

    #if PY_MAJOR_VERSION >= 3
    return PyUnicode_DecodeUTF8(line, length, "replace");
    #else
    return PyString_FromStringAndSize(line, length);
    #endif
}


static int parseLabel(PyObject *labels, PyObject *mapping, char *line,
                      int length, int tabs) {

    /* Append label to *labels*, extract identifier, and index label
       position in the list. Return 1 when successful, 0 on failure. */

    int i, ch, slash = 0, dash = 0;//, ipipe = 0, pipes[4] = {0, 0, 0, 0};

    length = labelLength(line, length, tabs);
    for (i = 0; i < length; i++) {
        ch = line[i];
        if (ch == '/' && slash == 0 && dash == 0)
            slash = i;
        else if (ch == '-' && slash > 0 && dash == 0)
            dash = i;
        //else if (line[i] == '|' && ipipe < 4)
        //    pipes[ipipe++] = i;
    }

    PyObject *label, *index;
    label = decodeLabel(line, i);
    #if PY_MAJOR_VERSION >= 3
    index = PyLong_FromSsize_t(PyList_Size(labels));
    #else
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

    if (slash > 0 && dash > slash) {
        Py_DECREF(label);
        label = decodeLabel(line, slash);
    }

    if (PyDict_Contains(mapping, label)) {
        PyObject *item = PyDict_GetItem(mapping, label); /* borrowed */
        if (PyList_Check(item)) {
            PyList_Append(item, index);
            Py_DECREF(index);
        } else {
            PyObject *list = PyList_New(2); /* new reference */
            PyList_SetItem(list, 0, item);
            Py_INCREF(item);
            PyList_SetItem(list, 1, index); /* steals reference, no DECREF */
            PyDict_SetItem(mapping, label, list);
            Py_DECREF(list);
        }
    } else {
        PyDict_SetItem(mapping, label, index);
        Py_DECREF(index);
    }

    Py_DECREF(label);
    return 1;
}


static PyObject *parseFasta(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the memory pointed by the
       Numpy array passed as Python object. */

    char *filename;
    PyArrayObject *msa;

    if (!PyArg_ParseTuple(args, "sO", &filename, &msa))
        return NULL;

    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping)
        return PyErr_NoMemory();

    char *line = malloc((FASTALINELEN) * sizeof(char));
    if (!line)
        return PyErr_NoMemory();

    char *data = (char *) PyArray_DATA(msa);

    int aligned = 1;
    char ch, errmsg[LENLABEL] = "failed to parse FASTA file at line ";
    long index = 0, count = 0;
    long iline = 0, i, seqlen = 0, curlen = 0;

    FILE *file = fopen(filename, "rb");
    while (fgets(line, FASTALINELEN, file) != NULL) {
        iline++;
        if (line[0] == '>') {
            if (seqlen != curlen) {
                if (seqlen) {
                    aligned = 0;
                    free(line);
                    free(data);
                    fclose(file);
                    PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
                    return NULL;
                } else
                    seqlen = curlen;
            }
            // `line + 1` is to omit `>` character
            count += parseLabel(labels, mapping, line + 1, FASTALINELEN, 0);
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

    free(line);
    if (aligned && seqlen != curlen) {
        PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
        return NULL;
    }

    npy_intp dims[2] = {index / seqlen, seqlen};
    PyArray_Dims arr_dims;
    arr_dims.ptr = dims;
    arr_dims.len = 2;
    PyArray_Resize(msa, &arr_dims, 0, NPY_CORDER);
    PyObject *result = Py_BuildValue("(OOOi)", msa, labels, mapping, count);
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

    long numseq = PyArray_DIMS(msa)[0], lenseq = PyArray_DIMS(msa)[1];

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
    char *seq = PyArray_DATA(msa);
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

static int splitSelexLine(char *line, long *lbeg, long *lend,
                          long *sbeg, long *send) {

    /* Find the label and the sequence in a SELEX/Stockholm *line*.  The
       sequence is the last whitespace separated item, and the label is the
       text before it.  Return 0 for markup lines, -1 for blank lines, 1
       for sequence lines, and 2 for the // line that ends the alignment. */

    long i;

    if (line[0] == '/' && line[1] == '/') {
        for (i = 2; line[i] && isspace((unsigned char) line[i]); i++);
        if (!line[i])
            return 2;
    }
    if (line[0] == '#')
        return 0;

    i = strlen(line);
    while (i > 0 && isspace((unsigned char) line[i - 1]))
        i--;
    if (i == 0)
        return -1;
    *send = i;
    while (i > 0 && !isspace((unsigned char) line[i - 1]))
        i--;
    *sbeg = i;
    while (i > 0 && isspace((unsigned char) line[i - 1]))
        i--;
    *lend = i;
    for (i = 0; i < *lend && isspace((unsigned char) line[i]); i++);
    *lbeg = i;
    return 1;
}


static int readLine(FILE *file, char **line, long *size) {

    /* Read a whole line from *file* into *line*, doubling the size of the
       buffer until the line fits.  Return 1 when a line is read, 0 at the
       end of file, and -1 when memory runs out. */

    long length;
    char *longer;

    if (fgets(*line, *size, file) == NULL)
        return 0;
    length = strlen(*line);
    while (length == *size - 1 && (*line)[length - 1] != '\n') {
        longer = realloc(*line, 2 * *size * sizeof(char));
        if (!longer)
            return -1;
        *line = longer;
        *size *= 2;
        if (fgets(*line + length, *size - length, file) == NULL)
            break;
        length += strlen(*line + length);
    }
    return 1;
}


static int sameLabel(PyObject *labels, long index, char *label, long length) {

    /* Return 1 if item *index* of *labels* is *label*, 0 otherwise. */

    int same = 0;
    PyObject *item = PyList_GetItem(labels, (Py_ssize_t) index); /* borrowed */
    PyObject *other = decodeLabel(label, labelLength(label, length, 1));
    #if PY_MAJOR_VERSION >= 3
    if (item && other && PyUnicode_Check(item))
        same = PyUnicode_Compare(item, other) == 0;
    #else
    if (item && other && PyString_Check(item))
        same = strcmp(PyString_AsString(item), PyString_AsString(other)) == 0;
    #endif
    Py_XDECREF(other);
    PyErr_Clear();
    return same;
}


static PyObject *parseSelex(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the the memory pointed by the
       Numpy array passed as Python object.  An alignment may be split into
       blocks separated by blank lines, and each block must repeat the labels
       of the first block in the same order.  */

    char *filename;
    PyArrayObject *msa;

    if (!PyArg_ParseTuple(args, "sO", &filename, &msa))
        return NULL;

    long lbeg = 0, lend = 0, sbeg = 0, send = 0;
    long size = SELEXLINELEN + 1, iline = 0;
    long numseq = 0, seqlen = 0, row, column, width = 0, count = 0;
    int kind, pass, inblock, status;
    char errmsg[LENLABEL] = "failed to parse SELEX/Stockholm file at line ";

    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping)
        return PyErr_NoMemory();
    char *line = malloc(size * sizeof(char));
    if (!line)
        return PyErr_NoMemory();
    char *data = (char *) PyArray_DATA(msa);

    FILE *file = fopen(filename, "rb");
    if (!file) {
        free(line);
        Py_DECREF(labels);
        Py_DECREF(mapping);
        return PyErr_SetFromErrnoWithFilename(PyExc_IOError, filename);
    }

    /* the first pass counts sequences in the first block and adds up the
       widths of blocks, the second pass copies each block into place */
    for (pass = 0; pass < 2; pass++) {
        rewind(file);
        iline = 0;
        row = 0;
        column = 0;
        inblock = 0;
        do {
            status = readLine(file, &line, &size);
            if (status < 0) {
                fclose(file);
                free(line);
                Py_DECREF(labels);
                Py_DECREF(mapping);
                return PyErr_NoMemory();
            }
            if (!status)
                kind = -1;
            else {
                iline++;
                kind = splitSelexLine(line, &lbeg, &lend, &sbeg, &send);
            }
            if (kind == 2) {
                /* nothing after the end of the alignment is read */
                kind = -1;
                status = 0;
            }
            if (kind == 0)
                continue;
            if (kind == -1) {
                /* a blank line or the end of file ends a block */
                if (inblock) {
                    if (!pass && !column)
                        numseq = row;
                    else if (row != numseq)
                        goto fail;
                    column += width;
                    row = 0;
                    inblock = 0;
                }
                continue;
            }
            if (lend == lbeg || (column && row == numseq))
                goto fail;
            if (!inblock) {
                width = send - sbeg;
                inblock = 1;
            } else if (send - sbeg != width)
                goto fail;
            if (pass) {
                if (!column)
                    count += parseLabel(labels, mapping, line + lbeg,
                                        lend - lbeg, 1);
                else if (!sameLabel(labels, row, line + lbeg, lend - lbeg))
                    goto fail;
                memcpy(data + row * seqlen + column, line + sbeg, width);
            }
            row++;
        } while (status);
        seqlen = column;
    }
    fclose(file);
    free(line);

    if (!numseq || !seqlen) {
        Py_DECREF(labels);
        Py_DECREF(mapping);
        PyErr_SetString(PyExc_IOError,
                        "no sequences found in SELEX/Stockholm file");
        return NULL;
    }

    npy_intp dims[2] = {numseq, seqlen};
    PyArray_Dims arr_dims;
    arr_dims.ptr = dims;
    arr_dims.len = 2;
    PyArray_Resize(msa, &arr_dims, 0, NPY_CORDER);
    PyObject *result = Py_BuildValue("(OOOi)", msa, labels, mapping, count);
    Py_DECREF(labels);
    Py_DECREF(mapping);

    return result;

  fail:
    fclose(file);
    free(line);
    Py_DECREF(labels);
    Py_DECREF(mapping);
    PyErr_SetString(PyExc_IOError, intcat(errmsg, iline));
    return NULL;
}


static PyObject *writeSelex(PyObject *self, PyObject *args, PyObject *kwargs) {

    /* Write MSA where inputs are: labels in the form of Python lists
    and sequences in the form of Python numpy array and write them in
    SELEX (default) or Stockholm format in the specified filename.  Labels
    are padded to *label_length*, and at least one space separates a label
    from its sequence.*/

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

    long numseq = PyArray_DIMS(msa)[0], lenseq = PyArray_DIMS(msa)[1];

    if (numseq != PyList_Size(labels)) {
        PyErr_SetString(PyExc_ValueError,
                        "size of labels and msa array does not match");
        return NULL;
    }

    FILE *file = fopen(filename, "wb");
    if (!file)
        return PyErr_SetFromErrnoWithFilename(PyExc_IOError, filename);

    int i, j;
    long pos = 0;
    char *seq = PyArray_DATA(msa);
    if (stockholm)
        fprintf(file, "# STOCKHOLM 1.0\n");

    #if PY_MAJOR_VERSION >= 3
    PyObject *plabel;
    #endif
    for (i = 0; i < numseq; i++) {
        #if PY_MAJOR_VERSION >= 3
        plabel = PyUnicode_AsEncodedString(
                PyList_GetItem(labels, (Py_ssize_t) i), "utf-8",
                               "label encoding");
        char *label =  PyBytes_AsString(plabel);
        #else
        char *label = PyString_AsString(PyList_GetItem(labels, (Py_ssize_t)i));
        #endif

        fputs(label, file);
        j = strlen(label);
        do
            fputc(' ', file);
        while (++j < label_length);
        fwrite(seq + pos, sizeof(char), lenseq, file);
        fputc('\n', file);
        pos += lenseq;

        #if PY_MAJOR_VERSION >= 3
        Py_DECREF(plabel);
        #endif
    }

    if (stockholm)
        fprintf(file, "//\n");

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


