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
#define NUMCHARS 27
const int twenty[20] = {1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 
                        14, 16, 17, 18, 19, 20, 22, 23, 25};
const int unambiguous[23] = {0, 1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 
                             15, 16, 17, 18, 19, 20, 21, 22, 23, 25};


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
    PyArrayObject *msa;
    
    if (!PyArg_ParseTuple(args, "sO", &filename, &msa))
        return NULL;
    
    long i = 0, lenseq = msa->dimensions[1];
    long lenline = 0, lenlast = 0, numlines = 0; 
    long size = lenseq + LENLABEL, iline = 0;
    char *line = malloc(size * sizeof(char));
    if (!line) {
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
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
    char *data = (char *) PyArray_DATA(msa);
    char clabel[LENLABEL], ckey[LENLABEL];
    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping) {
        free(line);
        PyErr_SetString(PyExc_MemoryError, 
                        "failed to create a list or dictionary object");
        return NULL;
    }

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
            PyErr_SetString(PyExc_IOError, 
                            "failed to parse msa, at line");
            return NULL;
        }

        /* parse sequence */
        for (i = 0; i < numlines; i++) {
            if (fgets(line, size, file) == NULL) {
                free(line);
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                return NULL;
            }
            for (j = 0; j < lenline; j++)
                data[index++] = line[j];
        }
        
        if (lenlast) {
            if (fgets(line, size, file) == NULL) {
                free(line);
                PyErr_SetString(PyExc_IOError, 
                                "failed to parse msa, at line");
                return NULL;
            }
            for (j = 0; j < lenlast; j++)
                data[index++] = line[j];
        }

        ccount++;
    }

    fclose(file);
    free(line);

    PyObject *result = Py_BuildValue("(OO)", labels, mapping);
    Py_DECREF(labels);
    Py_DECREF(mapping);

    return result;
}


static PyObject *parseSelex(PyObject *self, PyObject *args) {

    /* Parse sequences from *filename* into the the memory pointed by the
       Numpy array passed as Python object.  This function assumes that
       the sequences are aligned, i.e. start and end at the same column. */

    char *filename;
    PyArrayObject *msa;
    
    if (!PyArg_ParseTuple(args, "sO", &filename, &msa))
        return NULL;

    long i = 0, beg = 0, end = 0, lenseq = msa->dimensions[1]; 
    long size = lenseq + LENLABEL, iline = 0;
    char *line = malloc(size * sizeof(char));
    if (!line) {
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
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
    char *data = (char *) PyArray_DATA(msa);
    char clabel[LENLABEL], ckey[LENLABEL];
    PyObject *labels = PyList_New(0), *mapping = PyDict_New();
    if (!labels || !mapping) {
        free(line);
        PyErr_SetString(PyExc_MemoryError, 
                        "failed to create a list or dictionary object");
        return NULL;
    }

    int space = beg - 1; /* index of space character before sequence */
    while (fgets(line, size, file) != NULL) {
        iline++;
        if (line[0] == '#' || line[0] == '/' || line[0] == '%')
            continue;
            
        if (line[space] != ' ') {
            free(line);
            PyErr_SetString(PyExc_IOError, 
                            "failed to parse msa, at line");
            return NULL;
        } 

        /* parse label */
        if (!parseLabel(labels, mapping, line, clabel, ckey, ccount, size)) {
            free(line);
            PyErr_SetString(PyExc_IOError, 
                            "failed to parse msa, at line");
            return NULL;
        }
        
        /* parse sequence */
        for (i = beg; i < end; i++)
            data[index++] = line[i];
        ccount++;
    }
    fclose(file);
    free(line);
    
    PyObject *result = Py_BuildValue("(OO)", labels, mapping);
    Py_DECREF(labels);
    Py_DECREF(mapping);
    
    return result;
}


static PyObject *calcShannonEntropy(PyObject *self, PyObject *args,
                                    PyObject *kwargs) {

    PyArrayObject *msa, *entropy;
    int ambiquity = 1, dividend = 0;
    
    static char *kwlist[] = {"msa", "entropy", "ambiquity", "dividend", NULL};
        
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ii", kwlist,
                                     &msa, &entropy, &ambiquity, &dividend))
        return NULL;
    
    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
   
    if (entropy->dimensions[0] != lenseq) {
        PyErr_SetString(PyExc_IOError, 
                        "msa and entropy array shapes do not match");
        return NULL;
    }

    char *seq = (char *) PyArray_DATA(msa);
    double *ent = (double *) PyArray_DATA(entropy);

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
        if (ambiquity) {
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

    Py_RETURN_NONE;
}

static void sortJoint(double *joint[]) {
    
    int k, l, t;
    double *jrow, jp, *krow;
    
    /* X */
    jrow = joint[24]; 
    jp = jrow[24]; /* XX */
    if (jp > 0) {
        jp = jp / 20;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[t][t] = jp;
        }    
        jrow[24] = 0;
    }
    /* XB */ 
    jp = jrow[2]; 
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[4] = krow[14] = jp;
        }    
        jrow[2] = 0;
    }
    /* XJ */ 
    jp = jrow[10]; 
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[9] = krow[12] = jp;
        }    
        jrow[10] = 0;
    }
    /* XZ */ 
    jp = jrow[26]; 
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[5] = krow[17] = jp;
        }    
        jrow[26] = 0;
    }

    
    /* B */
    jrow = joint[2]; 
    jp = jrow[2]; /* BB */
    if (jp > 0) {
        joint[4][4] = joint[14][14] = jp / 2;
        jrow[2] = 0;
    }    
    /* BX */ 
    jp = jrow[24]; 
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[4][t] = joint[14][t] = jp;
        }    
        jrow[24] = 0;
    }
    /* BJ */  
    jp = jrow[10]; 
    if (jp > 0) {
        joint[4][9] = joint[14][9] = joint[4][12] = joint[14][12] = jp / 4;
        jrow[10] = 0;
    }    
    /* BZ */
    jp = jrow[26]; 
    if (jp > 0) {
        joint[4][5] = joint[14][5] = joint[4][17] = joint[14][17] = jp / 4;
        jrow[26] = 0;
    }  
    
    /* Z */
    jrow = joint[26]; 
    jp = jrow[26]; /* ZZ */
    if (jp > 0) {
        joint[5][5] = joint[17][17] = jp / 2;
        jrow[26] = 0;
    }
    /* ZX */ 
    jp = jrow[24]; 
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[5][t] = joint[17][t] = jp;
        }    
        jrow[24] = 0;
    }
    /* ZJ */  
    jp = jrow[10]; 
    if (jp > 0) {
        joint[5][9] = joint[17][9] = joint[5][12] = joint[17][12] = jp / 4;
        jrow[10] = 0;
    }    
    /* ZB */
    jp = jrow[2]; 
    if (jp > 0) {
        joint[5][4] = joint[17][4] = joint[5][14] = joint[17][14] = jp / 4;
        jrow[2] = 0;
    }  
    
    /* J */
    jrow = joint[10];
    jp = jrow[10]; /* JJ */
    if (jp > 0) {
        joint[9][9] = joint[12][12] = jp / 2;
        joint[10][10] = 0;
    }
    /* JX */ 
    jp = jrow[24]; 
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[9][t] = joint[12][t] = jp;
        }    
        jrow[24] = 0;
    }
    /* JB */
    jp = jrow[2]; 
    if (jp > 0) {
        joint[9][4] = joint[12][4] = joint[9][14] = joint[12][14] = jp / 4;
        jrow[2] = 0;
    }
    /* BZ */
    jp = jrow[26]; 
    if (jp > 0) {
        joint[9][5] = joint[12][5] = joint[9][17] = joint[12][17] = jp / 4;
        jrow[26] = 0;
    }  
    
            
    /*for (k = 0; k < NUMCHARS; k++) {*/
    for (t = 0; t < 23; t++) {
        k = unambiguous[t];
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
            joint[24][k] = 0;
        }
    }
    
}

static void zeroJoint(double *joint[]) {

    int k, l;
    double *jrow;        
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        for (l = 0; l < NUMCHARS; l++)
            jrow[l] = 0;
    }
}

static double calcMI(double *joint[], double *iprb, double *jprb, int debug) {

    int k, l;
    double *jrow, jp, mi = 0;
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        for (l = 0; l < NUMCHARS; l++) {
            jp = jrow[l];
            if (jp > 0)
                mi += jp * log(jp / iprb[k] / jprb[l]);
            if (debug && jp > 0)
                printf("%c%c %.4f / %.4f / %.4f\n", (char)k+64, (char)l+64, 
                        jp, iprb[k], jprb[l]);           

        }
    }
    return mi;
}

static void printJoint(double *joint[]) {
    int i, j;
    double csum[NUMCHARS], rsum, sum = 0, *row;
    printf("\nJoint probability matrix\n");    
    printf("  ");
    for (i = 0; i < NUMCHARS; i++) {
        printf("%c_%-2i ", i + 64, i);
        csum[i] = 0;
    }
    printf("\n");
    for (i = 0; i < NUMCHARS; i++) {
        rsum = 0;
        printf("%c ", i + 64);
        row = joint[i];
        for (j = 0; j < NUMCHARS; j++) {
            printf("%.2f ", row[j]*10);
            rsum += row[j];
            csum[j] += row[j];
            sum += row[j];          
        }
        printf("%.2f\n", rsum * 10);
    }
    printf("+ ");
    for (i = 0; i < NUMCHARS; i++)
        printf("%.2f ", csum[i] * 10);
    printf("%.2f\n", sum);
}

static void printProbs(double *probs[], long lenseq) {
    int i, j;
    double sum;
    double *row;
    printf("\nProbability matrix\n");    
    for (i = 0; i < NUMCHARS; i++)
        printf("%c_%-2i ", i + 64, i);
    printf("SUM\n");
    for (i = 0; i < lenseq; i++) {
        sum = 0;
        row = probs[i];
        for (j = 0; j < NUMCHARS; j++) {
            printf("%.2f ", row[j] * 10);
            sum += row[j];
        }
        printf("%.2f\n", sum);
    }
}

static PyObject *calcMutualInfo(PyObject *self, PyObject *args,
                                PyObject *kwargs) {

    PyArrayObject *msa, *mutinfo;
    int ambiquity = 1, turbo = 1, debug = 0;
    
    static char *kwlist[] = {"msa", "mutinfo", "ambiquity", "turbo", 
                             "debug", NULL};
        
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iii", kwlist, &msa, 
                                     &mutinfo, &ambiquity, &turbo, &debug))
        return NULL;

    /* check dimensions */
    long numseq = msa->dimensions[0], lenseq = msa->dimensions[1];
   
    if (mutinfo->dimensions[0] != lenseq || 
        mutinfo->dimensions[1] != lenseq) {
        PyErr_SetString(PyExc_IOError, 
                        "msa and mutinfo array shapes do not match");
        return NULL;
    }
    
    /* get pointers to data */
    
    char *seq = (char *) PyArray_DATA(msa); /*size: numseq x lenseq */
    double *mut = (double *) PyArray_DATA(mutinfo); /*size: lenseq x lenseq */

    long i, j;
    /* allocate memory */
    unsigned char *iseq = malloc(numseq * sizeof(unsigned char));
    if (!iseq) {
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
    }
    
    unsigned char **trans = malloc(lenseq * sizeof(unsigned char *));
    if (!trans) {
        turbo = 0;
    }
    
    if (turbo) {
        /* allocate rows that will store columns of MSA */
        trans[0] = iseq;
        for (i = 1; i < lenseq; i++) {
            trans[i] = malloc(numseq * sizeof(unsigned char));
            if (!trans[i]) {
                for (j = 1; j <= i; j++)
                    free(trans[j]);
                free(trans);
                turbo = 0;
            }
        }
    }
    unsigned char *jseq = iseq; /* so that we don't get uninitialized warning*/
    
    /* lenseq*27, a row for each column in the MSA */
    double **probs = malloc(lenseq * sizeof(double *));
    if (!probs) {
        if (turbo)
            for (j = 1; j < lenseq; j++)
                free(trans[j]);
        free(trans);
        free(iseq);
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
    }

    /* 27x27, alphabet characters and a gap*/
    double **joint = malloc(NUMCHARS * sizeof(double *));
    if (!joint) {
        if (turbo)
            for (j = 1; j < lenseq; j++)
                free(trans[j]);
        free(trans);
        free(iseq);
        free(probs);
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        return NULL;
    }
    
    for (i = 0; i < lenseq; i++) {
        probs[i] = malloc(NUMCHARS * sizeof(double));
        if (!probs[i]) {
            for (j = 0; j <= i; j++)
                free(probs[j]);
            free(probs);
            if (turbo)
                for (j = 1; j < lenseq; j++)
                    free(trans[j]);
            free(trans);
            free(iseq);
            PyErr_SetString(PyExc_MemoryError, "out of memory");
            return NULL;
        }
        for (j = 0; j < NUMCHARS; j++)
            probs[i][j] = 0;
    }

    for (i = 0; i < NUMCHARS; i++)  {
        joint[i] = malloc(NUMCHARS * sizeof(double));  
        if (!joint[i]) {
            for (j = 0; j <= i; j++)
                free(joint[j]);
            for (j = 0; j <= lenseq; j++)
                free(probs[j]);
            free(probs);
            free(joint);
            if (turbo)
                for (j = 1; j < lenseq; j++)
                    free(trans[j]);
            free(trans);
            free(iseq);
            PyErr_SetString(PyExc_MemoryError, "out of memory");
            return NULL;
        }
    }

    int a, b;
    long k, l, diff, offset;
    double *prow, *jrow, p_incr = 1. / numseq, prb = 0;
    
    if (debug)
        printProbs(probs, lenseq);
    
    i = 0;
    prow = probs[0];
    /* calculate first row of MI matrix, while calculating probabilities */
    for (j = 1; j < lenseq; j++) {
        jrow = probs[j];
        zeroJoint(joint);
        diff = j - 1;
        if (turbo)
            jseq = trans[j]; 
        for (k = 0; k < numseq; k++) {
            offset = k * lenseq;
            if (diff) {
                a = iseq[k];
            } else {
                a = (unsigned char) seq[offset + i];
                if (a > 90)
                    a -= 96;
                else
                    a -= 64;
                if (a < 1 || a > 26)
                    a = 0; /* gap character */
                iseq[k] = a;
            }
            
            b = (unsigned char) seq[offset + j];
            if (b > 90)
                b -= 96;
            else
                b -= 64;
            if (b < 1 || b > 26)
                b = 0; /* gap character */
            if (turbo)
                jseq[k] = b;
            joint[a][b] += p_incr;
            if (!diff)
                prow[a] += p_incr;
            jrow[b] += p_incr;
        }
        
        if (ambiquity) {
            for (k = 0; k < lenseq; k++) {
                prow = probs[k];
                prb = prow[2];
                if (prb > 0) { /* B -> D, N  */
                    prow[4] = prow[14] = prb / 2.;
                    prow[2] = 0;
                }
                prb = prow[10];
                if (prb > 0) { /* J -> I, L  */
                    prow[9] = prow[12] = prb / 2.;
                    prow[10] = 0;
                }
                prb = prow[26]; 
                if (prb > 0) { /* Z -> E, Q  */
                    prow[5] = prow[17] = prb / 2.;
                    prow[26] = 0;
                }
                if (prow[24] > 0) { /* X -> 20 AA */
                    prb = prow[24] / 20.; 
                    for (l = 0; l < 20; l++)
                        prow[twenty[l]] += prb;
                    prow[24] = 0;
                }
            }
            if (debug)
                printJoint(joint);
            sortJoint(joint);
            if (debug)
                printJoint(joint);
        }
        if (debug)
            printProbs(probs, lenseq);
        mut[i * lenseq + j] = mut[i + lenseq * j] = 
            calcMI(joint, probs[0], jrow, debug);
    }
    
    if (turbo)
        free(iseq);
        
    if (debug)
        printProbs(probs, lenseq);

    
    /* calculate rest of MI matrix */
    for (i = 1; i < lenseq; i++) {
        if (turbo)
            iseq = trans[i];
        for (j = i + 1; j < lenseq; j++) {
            zeroJoint(joint);

            if (turbo) {
                iseq = trans[j];
                for (k = 0; k < numseq; k++)
                    joint[iseq[k]][jseq[k]] += p_incr;
            } else {         
                diff = j - i - 1;
                for (k = 0; k < numseq; k++) {
                    offset = k * lenseq;
                    if (diff) {
                        a = iseq[k];
                    } else {
                        a = (unsigned char) seq[offset + i];
                        if (a > 90)
                            a -= 96;
                        else
                            a -= 64;
                        if (a < 1 || a > 26)
                            a = 0; /* gap character */
                        iseq[k] = a;
                    }
                    
                    b = (unsigned char) seq[offset + j];
                    if (b > 90)
                        b -= 96;
                    else
                        b -= 64;
                    if (b < 1 || b > 26)
                        b = 0; /* gap character */
                    joint[a][b] += p_incr;
                }
            }
            if (ambiquity)
                sortJoint(joint);
            mut[i * lenseq + j] = mut[i + lenseq * j] = 
                calcMI(joint, probs[i], probs[j], debug);
        }
    }

    /* free memory */
    for (i = 0; i < lenseq; i++){  
        free(probs[i]);
    }  
    free(probs);
    for (i = 0; i < NUMCHARS; i++){  
        free(joint[i]);
    }  
    free(joint);
    if (turbo)
        for (j = 1; j < lenseq; j++)
            free(trans[j]);
    free(trans);

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
