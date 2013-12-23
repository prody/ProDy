/* Author: Ahmet Bakan, Wenzhi Mao */

#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#define NUMCHARS 27
const int twenty[20] = {1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13,
                        14, 16, 17, 18, 19, 20, 22, 23, 25};
const int unambiguous[23] = {0, 1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14,
                             15, 16, 17, 18, 19, 20, 21, 22, 23, 25};


static PyObject *msaentropy(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *entropy;
    int ambiguity = 1, omitgaps = 0;

    static char *kwlist[] = {"msa", "entropy", "ambiguity", "omitgaps", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ii", kwlist,
                                     &msa, &entropy, &ambiguity, &omitgaps))
        return NULL;


    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);

    long number = PyArray_DIMS(msa)[0], length = PyArray_DIMS(msa)[1];

    char *seq = (char *) PyArray_DATA(msa);
    double *ent = (double *) PyArray_DATA(entropy);

    /* start here */
    long size = number * length;
    double count[128]; /* number of ASCII characters*/
    double shannon = 0, probability = 0, numgap = 0, denom = number;
    long i = 0, j = 0;

    double ambiguous = 0;
    int twenty[20] = {65, 67, 68, 69, 70, 71, 72, 73, 75, 76,
                      77, 78, 80, 81, 82, 83, 84, 86, 87, 89};
    for (i = 0; i < length; i++) {

        /* zero counters */
        for (j = 65; j < 91; j++)
            count[j] = 0;
        for (j = 97; j < 123; j++)
            count[j] = 0;

        /* count characters in a column*/
        for (j = i; j < size; j += length)
            count[(int) seq[j]]++;
        for (j = 65; j < 91; j++)
            count[j] += count[j + 32];

        /* handle ambiguous amino acids */
        if (ambiguity) {
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
        numgap = number;
        for (j = 65; j < 91; j++)
            numgap -= count[j];

        shannon = 0;
        denom = number;
        if (omitgaps)
            denom = number - numgap;
        else if (numgap > 0) {
            probability = numgap / number;
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
    return Py_BuildValue("O", entropy);
}


static void sortJoint(double **joint) {

    /* Sort probability of ambiguous amino acids. */

    int k, l, t;
    double *jrow, jp, *krow;

    /* X */
    jrow = joint[24];
    /* XX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 400;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            krow = joint[t];
            for (l = 0; l < 20; l++)
                joint[t][twenty[l]] += jp;
        }
        jrow[24] = 0;
    }
    /* XB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[4] += jp;
            krow[14] += jp;
        }
        jrow[2] = 0;
    }
    /* XJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[9] += jp;
            krow[12] += jp;
        }
        jrow[10] = 0;
    }
    /* XZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[5] += jp;
            krow[17] += jp;
        }
        jrow[26] = 0;
    }

    /* B */
    jrow = joint[2];
    /* BB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 4;
        joint[4][4] += jp;
        joint[4][14] += jp;
        joint[14][4] += jp;
        joint[14][14] += jp;
        jrow[2] = 0;
    }
    /* BX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[4][t] += jp;
            joint[14][t] += jp;
        }
        jrow[24] = 0;
    }
    /* BJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 4;
        joint[4][9] += jp;
        joint[4][12] += jp;
        joint[14][9] += jp;
        joint[14][12] += jp;
        jrow[10] = 0;
    }
    /* BZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 4;
        joint[4][5] += jp;
        joint[4][17] += jp;
        joint[14][5] += jp;
        joint[14][17] += jp;
        jrow[26] = 0;
    }

    /* Z */
    jrow = joint[26];
    /* ZZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 4;
        joint[5][5] += jp;
        joint[5][17] += jp;
        joint[17][5] += jp;
        joint[17][17] += jp;
        jrow[26] = 0;
    }
    /* ZX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[5][t] += jp;
            joint[17][t] += jp;
        }
        jrow[24] = 0;
    }
    /* ZJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 4;
        joint[5][9] += jp;
        joint[5][12] += jp;
        joint[17][9] += jp;
        joint[17][12] += jp;
        jrow[10] = 0;
    }
    /* ZB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 4;
        joint[5][4] += jp;
        joint[5][14] += jp;
        joint[17][4] += jp;
        joint[17][14] += jp;
        jrow[2] = 0;
    }

    /* J */
    jrow = joint[10];
    /* JJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 4;
        joint[9][9] += jp;
        joint[9][12] += jp;
        joint[12][9] += jp;
        joint[12][12] += jp;
        joint[10][10] = 0;
    }
    /* JX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[9][t] += jp;
            joint[12][t] += jp;
        }
        jrow[24] = 0;
    }
    /* JB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 4;
        joint[9][4] += jp;
        joint[9][14] += jp;
        joint[12][4] += jp;
        joint[12][14] += jp;
        jrow[2] = 0;
    }
    /* BZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 4;
        joint[9][5] += jp;
        joint[9][17] += jp;
        joint[12][5] += jp;
        joint[12][17] += jp;
        jrow[26] = 0;
    }

    /*for (k = 0; k < NUMCHARS; k++) {*/
    for (t = 0; t < 23; t++) {
        k = unambiguous[t];
        jrow = joint[k];

        /* B */
        jp = jrow[2];
        if (jp > 0) {
            jp = jp / 2;
            jrow[4] += jp;
            jrow[14] += jp;
            jrow[2] = 0;
        }
        jp = joint[2][k];
        if (jp > 0) {
            jp  = jp / 2;
            joint[4][k] += jp;
            joint[14][k] += jp;
            joint[2][k] = 0;
        }

        /* J */
        jp = jrow[10];
        if (jp > 0) {
            jp  = jp / 2;
            jrow[9] += jp;
            jrow[12] += jp;
            jrow[10] = 0;
        }
        jp = joint[10][k];
        if (jp > 0) {
            jp = jp / 2;
            joint[9][k] += jp;
            joint[12][k] += jp;
            joint[10][k] = 0;
        }

        /* Z */
        jp = jrow[26];
        if (jp > 0) {
            jp = jp / 2;
            jrow[5] += jp;
            jrow[17] += jp;
            jrow[26] = 0;
        }
        jp = joint[26][k];
        if (jp > 0) {
            jp = jp / 2;
            joint[5][k] += jp;
            joint[17][k] += jp;
            joint[26][k] = 0;
        }

        /* X */
        jp = jrow[24];
        if (jp > 0) {
            jp = jp / 20.;
            for (l = 0; l < 20; l++)
                jrow[twenty[l]] += jp;
            jrow[24] = 0;
        }
        jp = joint[24][k];
        if (jp > 0) {
            jp = jp / 20.;
            for (l = 0; l < 20; l++)
                joint[twenty[l]][k] += jp;
            joint[24][k] = 0;
        }
    }

}


static void zeroJoint(double **joint) {

    /* Fill NUMCHARSxNUMCHARS joint array with zeros. */

    int k, l;
    double *jrow;
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        for (l = 0; l < NUMCHARS; l++)
            jrow[l] = 0;
    }
}


static double calcMI(double **joint, double **probs, long i, long j, int dbg) {

    /* Calculate mutual information for a pair of columns in MSA. */

    int k, l;
    double *jrow, *iprb = probs[i], *jprb = probs[j], jp, mi = 0, inside;
    /*double isum = 0, jsum = 0, sum = 0;*/
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        /*isum += iprb[k];
        jsum += jprb[k];*/
        for (l = 0; l < NUMCHARS; l++) {
            jp = jrow[l];
            /*sum += jp;*/
            if (jp > 0) {
                inside = jp / iprb[k] / jprb[l];
                if (inside != 1)
                    mi += jp * log(inside);
            }
            /*if (dbg && jp > 0)
                printf("(%li,%li) %c%c %.4f / %.4f / %.4f\n",
                        i, j, (char)k+64, (char)l+64, jp, iprb[k], jprb[l]);*/
        }
    }
    /*if (dbg)
        if (sum != 1.00000 || isum != 1.00000 || jsum != 1.00000)
            printf("(%li,%li) %f/%f/%f\n", i, j, sum, isum, jsum);*/
    return mi;
}


static double jointEntropy(double **joint){

    double ent = 0.0, prob;
    double *row;
    int i, j;
    for (i = 0; i < NUMCHARS; i++) {
        row = joint[i];
        for (j = 0; j < NUMCHARS; j++) {
            prob = row[j];
            if (prob)
                ent -= prob * log(prob);
        }
    }
    return ent;
}


static void printJoint(double **joint, long k, long l) {

    /* Print joint probability matrix for debugging purposes. */

    int i, j;
    double csum[NUMCHARS], rsum, sum = 0, *row;
    printf("\nJoint probability matrix (%li,%li)\n", k, l);
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


static void printProbs(double **probs, long length) {

    /* Print probability matrix for debugging purposes. */

    int i, j;
    double sum;
    double *row;
    printf("\nProbability matrix\n");
    for (i = 0; i < NUMCHARS; i++)
        printf("%c_%-2i ", i + 64, i);
    printf("SUM\n");
    for (i = 0; i < length; i++) {
        sum = 0;
        row = probs[i];
        for (j = 0; j < NUMCHARS; j++) {
            printf("%.2f ", row[j] * 10);
            sum += row[j];
        }
        printf("%.2f\n", sum);
    }
}


static PyObject *msamutinfo(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *mutinfo;
    int ambiguity = 1, turbo = 1, debug = 0, norm = 0;

    static char *kwlist[] = {"msa", "mutinfo",
                             "ambiguity", "turbo", "norm", "debug", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiii", kwlist,
                                     &msa, &mutinfo,
                                     &ambiguity, &turbo, &norm, &debug))
        return NULL;

    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);

    /* check dimensions */
    long number = PyArray_DIMS(msa)[0], length = PyArray_DIMS(msa)[1];

    /* get pointers to data */
    char *seq = (char *) PyArray_DATA(msa); /*size: number x length */
    double *mut = (double *) PyArray_DATA(mutinfo);


    long i, j;
    /* allocate memory */
    unsigned char *iseq = malloc(number * sizeof(unsigned char));
    if (!iseq)
        return PyErr_NoMemory();

    /* hold transpose of the sorted character array */
    unsigned char **trans = malloc(length * sizeof(unsigned char *));
    if (!trans) {
        turbo = 0;
    }

    if (turbo) {
        /* allocate rows that will store columns of MSA */
        trans[0] = iseq;
        for (i = 1; i < length; i++) {
            trans[i] = malloc(number * sizeof(unsigned char));
            if (!trans[i]) {
                for (j = 1; j < i; j++)
                    free(trans[j]);
                free(trans);
                turbo = 0;
            }
        }
    }
    unsigned char *jseq = iseq; /* so that we don't get uninitialized warning*/

    /* length*27, a row for each column in the MSA */
    double **probs = malloc(length * sizeof(double *)), *prow;
    if (!probs) {
        if (turbo)
            for (j = 1; j < length; j++)
                free(trans[j]);
        free(trans);
        free(iseq);
        return PyErr_NoMemory();
    }

    /* 27x27, alphabet characters and a gap*/
    double **joint = malloc(NUMCHARS * sizeof(double *)), *jrow;
    if (!joint) {
        if (turbo)
            for (j = 1; j < length; j++)
                free(trans[j]);
        free(trans);
        free(iseq);
        free(probs);
        return PyErr_NoMemory();
    }

    for (i = 0; i < length; i++) {
        prow = malloc(NUMCHARS * sizeof(double));
        if (!prow) {
            for (j = 0; j < i; j++)
                free(probs[j]);
            free(probs);
            free(joint);
            if (turbo)
                for (j = 1; j < length; j++)
                    free(trans[j]);
            free(trans);
            free(iseq);
            return PyErr_NoMemory();
        }
        probs[i] = prow;
        for (j = 0; j < NUMCHARS; j++)
            prow[j] = 0;
    }

    for (i = 0; i < NUMCHARS; i++)  {
        joint[i] = malloc(NUMCHARS * sizeof(double));
        if (!joint[i]) {
            for (j = 0; j < i; j++)
                free(joint[j]);
            free(joint);
            for (j = 0; j < length; j++)
                free(probs[j]);
            free(probs);
            if (turbo)
                for (j = 1; j < length; j++)
                    free(trans[j]);
            free(trans);
            free(iseq);
            return PyErr_NoMemory();
        }
    }

    if (debug)
        printProbs(probs, length);

    unsigned char a, b;
    long k, l, diff, offset;
    double p_incr = 1. / number, prb = 0;
    prow = probs[0];


    /* START mutinfo calculation */
    /* calculate first row of MI matrix and all column probabilities */
    i = 0;
    mut[0] = 0;
    for (j = 1; j < length; j++) {
        mut[j * length + j] = 0; /* using empty, so needed for diagonal */
        jrow = probs[j];
        zeroJoint(joint);
        diff = j - 1;
        if (turbo) /* in turbo mode, there is a row for refined sequences */
            jseq = trans[j];
        for (k = 0; k < number; k++) {
            offset = k * length;
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
                prow[a] += p_incr;
            }

            b = (unsigned char) seq[offset + j];
            if (b > 90)
                b -= 96;
            else
                b -= 64;
            if (b < 1 || b > 26)
                b = 0; /* gap character */
            if (turbo)  /* we keep the refined chars for all sequences*/
                jseq[k] = b;
            joint[a][b] += p_incr;
            jrow[b] += p_incr;
        }

        if (ambiguity) {

            if (debug)
                printProbs(probs, length);
            if (diff)
                k = j;
            else
                k = 0;
            for (; k <= j; k++) {
                prow = probs[k];
                prb = prow[2];
                if (prb > 0) { /* B -> D, N  */
                    prb = prb / 2.;
                    prow[4] += prb;
                    prow[14] += prb;
                    prow[2] = 0;
                }
                prb = prow[10];
                if (prb > 0) { /* J -> I, L  */
                    prb = prb / 2.;
                    prow[9] += prb;
                    prow[12] += prb;
                    prow[10] = 0;
                }
                prb = prow[26];
                if (prb > 0) { /* Z -> E, Q  */
                    prb = prb / 2.;
                    prow[5] += prb;
                    prow[17] += prb;
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
                printProbs(probs, length);
            if (debug)
                printJoint(joint, i, j);
            sortJoint(joint);
            if (debug)
                printJoint(joint, i, j);
        }
        if (norm)
            mut[j] = mut[length * j] = calcMI(joint, probs, i, j, debug) /
                                      jointEntropy(joint);
        else
            mut[j] = mut[length * j] = calcMI(joint, probs, i, j, debug);
    }
    if (debug)
        printProbs(probs, length);
    if (turbo)
        free(iseq);


    /* calculate rest of MI matrix */
    long ioffset;
    for (i = 1; i < length; i++) {
        ioffset = i * length;
        if (turbo)
            iseq = trans[i];

        for (j = i + 1; j < length; j++) {
            zeroJoint(joint);

            if (turbo) {
                jseq = trans[j];
                for (k = 0; k < number; k++)
                    joint[iseq[k]][jseq[k]] += p_incr;

            } else {
                diff = j - i - 1;
                for (k = 0; k < number; k++) {
                    offset = k * length;
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
            if (ambiguity)
                sortJoint(joint);
        if (norm)
            mut[ioffset + j] = mut[i + length * j] =
                calcMI(joint, probs, i, j, debug) / jointEntropy(joint);
        else
            mut[ioffset + j] = mut[i + length * j] =
                calcMI(joint, probs, i, j, debug);
        }
    }

    /* free memory */
    for (i = 0; i < length; i++){
        free(probs[i]);
    }
    free(probs);
    for (i = 0; i < NUMCHARS; i++){
        free(joint[i]);
    }
    free(joint);
    if (turbo)
        for (j = 1; j < length; j++)
            free(trans[j]);
    free(trans);

    return Py_BuildValue("O", mutinfo);
}


static PyObject *msaocc(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *occ;
    int dim;
    int count = 0;

    static char *kwlist[] = {"msa", "occ", "dim", "count", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOi|i", kwlist,
                                     &msa, &occ, &dim, &count))
        return NULL;

    long number = PyArray_DIMS(msa)[0], length = PyArray_DIMS(msa)[1];
    char *seq = (char *) PyArray_DATA(msa), *row, ch;
    double *cnt = (double *) PyArray_DATA(occ);

    long i, j, *k;
    if (dim)
        k = &j;
    else
        k = &i;

    for (i = 0; i < number; i++) {
        row = seq + i * length;
        for (j = 0; j < length; j++) {
            ch = row[j];
            if ((64 < ch && ch < 91) || (96 < ch && ch < 123))
                cnt[*k]++;
        }
    }

    if (!count) {
        double divisor;
        if (dim)
            divisor = 1. * number;
        else
            divisor = 1. * length;
        for (i = 0; i < PyArray_DIMS(msa)[dim]; i++)
            cnt[i] /= divisor;
    }
    return Py_BuildValue("O", occ);
}


static double calcOMES(double **joint, double **probs, long i, long j, int n) {

    /* Calculate OMES for a pair of columns in MSA. */

    int k, l;
    double *jrow, *iprb = probs[i], *jprb = probs[j], jp, omes = 0, inside;
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        for (l = 0; l < NUMCHARS; l++) {
            jp = jrow[l];
            inside = iprb[k] * jprb[l];
            if (inside != 0)
                omes += n * (jp - inside) * (jp - inside) / inside;
        }
    }
    return omes;
}


static PyObject *msaomes(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *omes;
    int ambiguity = 1, turbo = 1, debug = 0;

    static char *kwlist[] = {"msa", "omes",
                             "ambiguity", "turbo", "debug", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iii", kwlist,
                                     &msa, &omes,
                                     &ambiguity, &turbo, &debug))
        return NULL;

    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);

    /* check dimensions */
    long number = PyArray_DIMS(msa)[0], length = PyArray_DIMS(msa)[1];

    /* get pointers to data */
    char *seq = (char *) PyArray_DATA(msa); /*size: number x length */
    double *data = (double *) PyArray_DATA(omes);

    long i, j;
    /* allocate memory */
    unsigned char *iseq = malloc(number * sizeof(unsigned char));
    if (!iseq)
        return PyErr_NoMemory();

    /* hold transpose of the sorted character array */
    unsigned char **trans = malloc(length * sizeof(unsigned char *));
    if (!trans) {
        turbo = 0;
    }

    if (turbo) {
        /* allocate rows that will store columns of MSA */
        trans[0] = iseq;
        for (i = 1; i < length; i++) {
            trans[i] = malloc(number * sizeof(unsigned char));
            if (!trans[i]) {
                for (j = 1; j < i; j++)
                    free(trans[j]);
                free(trans);
                turbo = 0;
            }
        }
    }
    unsigned char *jseq = iseq; /* so that we don't get uninitialized warning*/

    /* length*27, a row for each column in the MSA */
    double **probs = malloc(length * sizeof(double *)), *prow;
    if (!probs) {
        if (turbo)
            for (j = 1; j < length; j++)
                free(trans[j]);
        free(trans);
        free(iseq);
        return PyErr_NoMemory();
    }

    /* 27x27, alphabet characters and a gap*/
    double **joint = malloc(NUMCHARS * sizeof(double *)), *jrow;
    if (!joint) {
        if (turbo)
            for (j = 1; j < length; j++)
                free(trans[j]);
        free(trans);
        free(iseq);
        free(probs);
        return PyErr_NoMemory();
    }

    for (i = 0; i < length; i++) {
        prow = malloc(NUMCHARS * sizeof(double));
        if (!prow) {
            for (j = 0; j < i; j++)
                free(probs[j]);
            free(probs);
            free(joint);
            if (turbo)
                for (j = 1; j < length; j++)
                    free(trans[j]);
            free(trans);
            free(iseq);
            return PyErr_NoMemory();
        }
        probs[i] = prow;
        for (j = 0; j < NUMCHARS; j++)
            prow[j] = 0;
    }

    for (i = 0; i < NUMCHARS; i++)  {
        joint[i] = malloc(NUMCHARS * sizeof(double));
        if (!joint[i]) {
            for (j = 0; j < i; j++)
                free(joint[j]);
            free(joint);
            for (j = 0; j < length; j++)
                free(probs[j]);
            free(probs);
            if (turbo)
                for (j = 1; j < length; j++)
                    free(trans[j]);
            free(trans);
            free(iseq);
            return PyErr_NoMemory();
        }
    }

    if (debug)
        printProbs(probs, length);

    unsigned char a, b;
    long k, l, diff, offset;
    double p_incr = 1. / number;
    double prb = 0;
    prow = probs[0];

    /* START OMES calculation */
    /* calculate first row of OMES matrix and all column probabilities */
    i = 0;
    data[0] = 0;
    for (j = 1; j < length; j++) {
        data[j * length + j] = 0; /* using empty, so needed for diagonal */
        jrow = probs[j];
        zeroJoint(joint);
        diff = j - 1;
        if (turbo) /* in turbo mode, there is a row for refined sequences */
            jseq = trans[j];
        for (k = 0; k < number; k++) {
            offset = k * length;
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
                prow[a] += p_incr;
            }

            b = (unsigned char) seq[offset + j];
            if (b > 90)
                b -= 96;
            else
                b -= 64;
            if (b < 1 || b > 26)
                b = 0; /* gap character */
            if (turbo)  /* we keep the refined chars for all sequences*/
                jseq[k] = b;
            joint[a][b] += p_incr;
            jrow[b] += p_incr;
        }

        if (ambiguity) {

            if (debug)
                printProbs(probs, length);
            if (diff)
                k = j;
            else
                k = 0;
            for (; k <= j; k++) {
                prow = probs[k];
                prb = prow[2];
                if (prb > 0) { /* B -> D, N  */
                    prb = prb / 2.;
                    prow[4] += prb;
                    prow[14] += prb;
                    prow[2] = 0;
                }
                prb = prow[10];
                if (prb > 0) { /* J -> I, L  */
                    prb = prb / 2.;
                    prow[9] += prb;
                    prow[12] += prb;
                    prow[10] = 0;
                }
                prb = prow[26];
                if (prb > 0) { /* Z -> E, Q  */
                    prb = prb / 2.;
                    prow[5] += prb;
                    prow[17] += prb;
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
                printProbs(probs, length);
            if (debug)
                printJoint(joint, i, j);
            sortJoint(joint);
            if (debug)
                printJoint(joint, i, j);
        }
        data[j] = data[length * j] = calcOMES(joint, probs, i, j, number);
    }
    if (debug)
        printProbs(probs, length);
    if (turbo)
        free(iseq);

    /* calculate rest of OMES matrix */
    long ioffset;
    for (i = 1; i < length; i++) {
        ioffset = i * length;
        if (turbo)
            iseq = trans[i];

        for (j = i + 1; j < length; j++) {
            zeroJoint(joint);

            if (turbo) {
                jseq = trans[j];
                for (k = 0; k < number; k++)
                    joint[iseq[k]][jseq[k]] += p_incr;

            } else {
                diff = j - i - 1;
                for (k = 0; k < number; k++) {
                    offset = k * length;
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
            if (ambiguity)
                sortJoint(joint);
        data[ioffset + j] = data[i + length * j] =
        calcOMES(joint, probs, i, j, number);
        }
    }

    /* free memory */
    for (i = 0; i < length; i++){
        free(probs[i]);
    }
    free(probs);
    for (i = 0; i < NUMCHARS; i++){
        free(joint[i]);
    }
    free(joint);
    if (turbo)
        for (j = 1; j < length; j++)
            free(trans[j]);
    free(trans);

    return Py_BuildValue("O", omes);
}


static PyObject *msasca(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *scainfo;
    int turbo = 1;
    static char *kwlist[] = {"msa", "sca", "turbo", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|i", kwlist,
                                     &msa, &scainfo, &turbo))
        return NULL;
    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);
    /* check dimensions */
    long number = PyArray_DIMS(msa)[0], length = PyArray_DIMS(msa)[1];
    /* get pointers to data */
    char *seq = (char *) PyArray_DATA(msa); /*size: number x length */
    double *sca = (double *) PyArray_DATA(scainfo);

    long i, j, k;
    double q[NUMCHARS] = {0., 0.073, 0., 0.025, 0.05, 0.061, 0.042, 0.072,
        0.023, 0.053, 0., 0.064, 0.089, 0.023, 0.043, 0., 0.052, 0.04, 0.052,
        0.073, 0.056, 0., 0.063, 0.013, 0., 0.033, 0.};
    int qlist[21] = {0, 1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13,
                        14, 16, 17, 18, 19, 20, 22, 23, 25};

    /* weighted probability matrix length*27 */
    double **wprob = malloc(length * sizeof(double *));
    if (!wprob)
        return PyErr_NoMemory();

    /* each row of weighted probability */
    for (i = 0; i < length; i++) {
        wprob[i] = malloc(NUMCHARS * sizeof(double));
        if (!wprob[i]) {
            for (j = 0; j < i; j++)
                free(wprob[j]);
            free(wprob);
            return PyErr_NoMemory();
        }
        for (j = 0; j < NUMCHARS; j++)
            wprob[i][j] = 0;
    }

    /* single column probability */
    double *prob;

    /* weighted x~ matrix array */
    double **wx = malloc(length * sizeof(double *));
    if (!turbo)
        free(wx);
    if (!wx)
        turbo = 0;
    if (turbo) {
        for (i = 0; i < length; i++) {
            wx[i] = malloc(number * sizeof(double));
            if (!wx[i]) {
                for (j = 0; j < i; j++)
                    free(wx[j]);
                free(wx);
                turbo = 0;
            }
        }
    }

    /* build weighted probability prob */
    for (i = 0; i < length; i++){
        prob = wprob[i];
        double phi[NUMCHARS];
        for (j = 0; j < NUMCHARS; j++){
            prob[j] = 0.0;
            phi[i] = 0.0;
        }
        for (j=0; j<number; j++){
            int temp = seq[j * length + i];
            temp = (temp > 96) ?  temp - 97 : temp - 65;
            if ((temp >= 0) && (temp <= 25))
                prob[temp + 1] += 1.0 ;
        }
        for (j=0; j<NUMCHARS; j++){
            prob[j] = prob[j] / number;
        }
        if (prob[2] > 0){ /* B -> D, N  */
            prob[4] += prob[2] / 2.;
            prob[14] += prob[2] / 2.;
            prob[2] = 0.;
        }
        if (prob[10] > 0){ /* J -> I, L  */
            prob[9] += prob[10] / 2.;
            prob[12] += prob[10] / 2.;
            prob[10] = 0.;
        }
        if (prob[26] > 0){ /* Z -> E, Q  */
            prob[4] += prob[26] / 2.;
            prob[17] += prob[26] / 2.;
            prob[26] = 0.;
        }
        if (prob[24] > 0) { /* X -> 20 AA */
            for (k = 0; k < 20; k++)
                prob[twenty[k]] += prob[24] / 20.;
            prob[24] = 0.;
        }
        double sum=0.0;
        for (j = 0; j < 21; j++){
            phi[qlist[j]] = (prob[qlist[j]] == 0.0 || q[qlist[j]] == 0.0
                            || prob[qlist[j]] == 1.0 || q[qlist[j]] == 1.0)
                    ? 0.0
                    : log(prob[qlist[j]] * (1 - q[qlist[j]]) /
                        (1 - prob[qlist[j]]) / q[qlist[j]]);
            phi[qlist[j]] = (phi[qlist[j]] >= 0.) ?
                phi[qlist[j]] : -phi[qlist[j]];
            prob[qlist[j]] = prob[qlist[j]] * phi[qlist[j]];
            sum += prob[qlist[j]] * prob[qlist[j]];
            prob[qlist[j]] = prob[qlist[j]] * phi[qlist[j]];
        }
        sum = sqrt(sum);
        if (sum == 0.)
            for (j = 0; j < 21; j++){
                prob[qlist[j]] = 0.;
            }
        else
            for (j = 0; j < 21; j++){
                prob[qlist[j]] = prob[qlist[j]] / sum;
            }
        prob[2] = (prob[4] + prob[14]) /2.0;
        prob[10] = (prob[9] + prob[12]) /2.0;
        prob[26] = (prob[4] + prob[17]) /2.0;
        sum =0.0;
        for (k = 0; k < 20; k++)
            sum += prob[twenty[k]];
        sum = sum / 20.0;
        prob[24] = sum;
        if (turbo){
            for (j = 0; j < number; j++){
                int temp = seq[j * length + i];
                temp = (temp > 96) ? temp - 97 : temp - 65;
                if (temp >= 0 && temp <= 25)
                    wx[i][j] = prob[temp + 1];
                else
                    wx[i][j] = 0.0;
            }
        }
    }

    /* Calculate SCA Matrix*/
    for (i=0;i<length;i++){
        for (j = i;j<length;j++){
            double *icol, *jcol, sumi=0.0, sumj=0.0, sum=0.0;
            if (turbo){
                icol=wx[i];
                jcol=wx[j];
                for (k=0; k< number; k++){
                    sumi += icol[k];
                    sumj += jcol[k];
                    sum += icol[k]*jcol[k];
                }
            }
            else{
                for (k = 0; k < number; k++){
                    int tempi = (seq[k*length + i] > 96) ?
                    seq[k * length + i] - 97 : seq[k * length + i] - 65;
                    double xi = (tempi >= 0 && tempi <= 25) ?
                        wprob[i][tempi + 1] : wprob[i][0];
                    int tempj = (seq[k * length + j] > 96) ?
                        seq[k * length + j] - 97 : seq[k * length + j] - 65;
                    double xj = (tempj >= 0 && tempj <= 25) ?
                        wprob[j][tempj + 1] : wprob[j][0];
                    sumi += xi;
                    sumj += xj;
                    sum += xi * xj;
                }
            }
            sum /= number;
            sumj /= number;
            sumi /= number;
            sum = sum - sumi * sumj;
            sum = sum >= 0 ? sum : -sum ;
            sca[i * length + j] = sca[j * length + i] = sum;
        }
    }

    /* free memory */
    for (j = 1; j < length; j++)
        free(wprob[j]);
    free(wprob);
    if (turbo){
        for (j = 1; j < length; j++)
            free(wx[j]);
        free(wx);
    }

    return Py_BuildValue("O", scainfo);
}


static PyObject *msameff(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa,*pythonw;
    double theta = 0.0;
    int meff_only = 1, refine = 0;
    int alignlist[26] = {1, 0, 2, 3, 4, 5, 6, 7, 8, 0, 9, 10, 11, 12,
             0, 13, 14, 15, 16, 17, 0, 18, 19, 0, 20, 0};
    static char *kwlist[] = {"msa", "theta", "meff_only", "refine", "w", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Odii|O", kwlist,
                                     &msa, &theta, &meff_only, &refine,
                                     &pythonw))
        return NULL;
    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);
    /* check dimensions */
    long number = PyArray_DIMS(msa)[0], length = PyArray_DIMS(msa)[1];
    long i, j, k, l = 0;
    /* get pointers to data */
    char *seq = (char *) PyArray_DATA(msa); /*size: number x length */

    /*Set ind and get l first.*/
    int *ind = malloc(length * sizeof(int));
    if (!ind) {
        return PyErr_NoMemory();
    }

    if (!refine){
        for (i = 0; i < length; i++){
            l += 1;
            ind[i] = l;
        }
    }
    else{
        for (i = 0; i < length; i++){
            if (seq[i] <= 90 && seq[i] >= 65){
                l += 1;
                ind[i] = l;
            }
            else
                ind[i] = 0;
        }
    }

    /*Use l to set align and w size.*/
    int *align = malloc(number * l * sizeof(int));
    if (!align) {
        free(ind);
        return PyErr_NoMemory();
    }
    for (i = 0; i < number * l; i++){
        align[i] = 0;
    }
    double *w = malloc(number * sizeof(double));
    if (!w) {
        free(ind);
        free(align);
        return PyErr_NoMemory();
    }

    #define align(x,y) align[(x)*l+(y)]

    /*Set align matrix*/
    for (i = 0; i < number; i++){
        for (j = 0; j < length; j++){
            if (ind[j] != 0){
                if (seq[i*length+j] >= 65 && seq[i*length+j] <= 90)
                    align(i,ind[j]-1) = alignlist[seq[i*length+j] - 65];
                else
                    align(i,ind[j]-1) = 0;
            }
        }
    }

    /*Calculate weight(w) for each sequence, sum of w is Meff*/
    for (i = 0; i < number; i++)
        w[i] = 1.;
    for (i = 0; i < number; i++)
        for (j = i+1; j < number; j++){
            double temp = 0.;
            for (k = 0; k < l; k++){
                if (align(i,k) != align(j,k))
                    temp += 1.;
            }
            temp /= l;
            if (temp < theta){
                w[i] += 1.;
                w[j] += 1.;
            }
        }
    double meff = 0.0;
    for (i = 0; i < number; i++){
        w[i] = 1./ w[i];
        meff += w[i];
    }

    #undef align

    /*Clean up memory.*/
    free(ind);
    if (meff_only == 1){
        free(align);
        free(w);
        return Py_BuildValue("d", meff);
    }
    else if (meff_only == 2){
        for (i = 0; i < number; i++)
            w[i] /= meff;
        return Py_BuildValue("dllll", meff, number, l , w, align);
    }
    else {
        free(align);
        pythonw = PyArray_GETCONTIGUOUS(pythonw);
        double *pw = (double *) PyArray_DATA(pythonw);
        for (i = 0; i < number; i++){
            pw[i]=w[i];
        }
        free(w);
        return Py_BuildValue("dO",meff,pythonw);
    }
}


static PyObject *msadipretest(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyArrayObject *msa;
    int refine = 0;
    int alignlist[26] = {1, 0, 2, 3, 4, 5, 6, 7, 8, 0, 9, 10, 11, 12,
             0, 13, 14, 15, 16, 17, 0, 18, 19, 0, 20, 0};
    static char *kwlist[] = {"msa", "refine", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Oi", kwlist,
                                     &msa, &refine))
        return NULL;
    msa = PyArray_GETCONTIGUOUS(msa);
    long number = PyArray_DIMS(msa)[0], length = PyArray_DIMS(msa)[1];
    char *seq = (char *) PyArray_DATA(msa);
    long i, j, k = 0, l = 0;
    int *ind = malloc(length * sizeof(int));
    if (!ind)
        return PyErr_NoMemory();
    if (!refine){
        for (i = 0; i < length; i++)
            ind[i] = i + 1;
        l = length;
    }
    else
        for (i = 0; i < length; i++)
            if (seq[i] <= 90 && seq[i] >= 65){
                l += 1;
                ind[i] = l;
            }
            else
                ind[i] = 0;
    for (i = 0; i < number; i++)
        for (j = 0; j < length; j++)
            if (ind[j])
                if (seq[i*length+j] >= 65 && seq[i*length+j] <= 90)
                    k = alignlist[seq[i*length+j]-65]>k?
                        alignlist[seq[i*length+j]-65]:k;
    free(ind);
    return Py_BuildValue("ii",l,k);
}


static PyObject *msadirectinfo1(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *cinfo, *pinfo;
    double theta = 0.2, pseudocount_weight = 0.5;
    int refine = 0, q = 0;
    static char *kwlist[] = {"msa", "c", "prob", "theta", "pseudocount_weight",
                             "refine", "q", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOddi|i", kwlist,
                                     &msa, &cinfo, &pinfo, &theta,
                                     &pseudocount_weight, &refine, &q))
        return NULL;
    long i, j, k, k1, k2;
    cinfo = PyArray_GETCONTIGUOUS(cinfo);
    pinfo = PyArray_GETCONTIGUOUS(pinfo);
    double *c = (double *) PyArray_DATA(cinfo);
    double *prob = (double *) PyArray_DATA(pinfo);

    /*Calculate meff, w and align.*/
    double meff = -1.;
    long number = 0, l = 0;
    int *align = NULL;
    double *w = NULL;
    PyObject *meffinfo;
    meffinfo = msameff(NULL, Py_BuildValue("(O)", msa),
             Py_BuildValue("{s:d,s:i,s:i}", "theta", theta, "meff_only", 2,
                 "refine", refine));
    if (!PyArg_ParseTuple(meffinfo, "dllll", &meff, &number, &l, &w, &align))
        return NULL;

    /*Build single probablity. use pseudocount_weight to weight it.*/
    double pse_weight_val = pseudocount_weight / q;
    double pro_weight = 1. - pseudocount_weight;
    for (i = 0; i < q*l; i++)
        prob[i] = pse_weight_val;
    #define prob(x,y) prob[(x)*q + (y)]
    #define align(x,y) align[(x)*l + (y)]
    for (i = 0; i < number; i++)
        for (j = 0; j < l; j++)
            prob(j, align(i,j)) += pro_weight * w[i];

    /*Calculate C matrix.*/
    double *joint = malloc(q*q*sizeof(double));
    if (!joint){
        free(w);
        free(align);
        return PyErr_NoMemory();
    }
    #define joint(x,y) joint[(x)*q + (y)]
    #define c(x,y) c[(x)*l*(q-1) + (y)]
    for (i = 0; i < l; i++){
        for (j = i; j < l; j++){

            if (i==j){
                for (k = 0; k < q*q; k++)
                    joint[k] = 0.;
                pse_weight_val = pseudocount_weight / q;
                for (k = 0; k < q; k++)
                    joint(k,k) = pse_weight_val;
            }
            else{
                pse_weight_val = pseudocount_weight / q / q;
                for (k = 0; k < q*q; k++)
                    joint[k] = pse_weight_val;
            }

            for (k = 0; k < number; k++){
                joint(align(k,i), align(k,j)) += pro_weight * w[k];
            }

            for (k1 = 0; k1 < q-1; k1++){
                for(k2 = 0; k2 < q-1; k2++){
                    c((q-1)*j+k2, (q-1)*i+k1) = c((q-1)*i+k1, (q-1)*j+k2) = joint(k1,k2) - prob(i,k1) * prob(j,k2);
                    // c((q-1)*j+k2, (q-1)*i+k1) = c((q-1)*i+k1, (q-1)*j+k2);
                }
            }
        }
    }

    free(w);
    free(align);
    free(joint);
    #undef prob
    #undef align
    #undef joint
    #undef c
    return Py_BuildValue("dllOO", meff, number, l, cinfo, pinfo);
}


static PyObject *msadirectinfo2(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *cinfo, *pinfo, *diinfo;
    long number = 0, l = 0, q = 0;
    static char *kwlist[] = {"n", "l", "c", "p", "di", "q", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "llOOOl", kwlist,
                                     &number, &l, &cinfo, &pinfo, &diinfo, &q))
        return NULL;
    cinfo = PyArray_GETCONTIGUOUS(cinfo);
    pinfo = PyArray_GETCONTIGUOUS(pinfo);
    diinfo = PyArray_GETCONTIGUOUS(diinfo);
    double *c = (double *) PyArray_DATA(cinfo);
    double *prob = (double *) PyArray_DATA(pinfo);
    double *di = (double *) PyArray_DATA(diinfo);

    long i, j, k1, k2;
    double *w = malloc(q * q * sizeof(double));
    if (!w)
        return NULL;
    for (i = 0; i < q*q; i++){
        w[i] = 0.0;
    }

    #define w(x,y) w[(x)*q+(y)]
    #define c(x,y) c[(x)*l*(q-1) + (y)]
    #define prob(x,y) prob[(x)*q + (y)]
    #define di(x,y) di[(x)*l + (y)]

    double epsilon = 1e-4, tiny = 1.0e-100;
    double diff = 1.0, sum1 = 0.0, sum2 = 0.0, sumpdir = 0.0, sumdi = 0.0;
    double *mu1 = malloc(q*sizeof(double)), *mu2 = malloc(q*sizeof(double));
    double *scra1 = malloc(q*sizeof(double)), *scra2 = malloc(q*sizeof(double));
    for (i = 0; i < l; i++){
        di(i,i) = 0.0;
        for (j = i+1; j < l; j++){
            for (k1 = 0; k1 < q-1; k1++){
                for (k2 = 0; k2 < q-1; k2++){
                    w(k1,k2) = exp(- c((q-1)*i + k1, (q-1)*j + k2));
                }
            }
            for (k1 = 0; k1 < q; k1++){
                w(q-1, k1) = w(k1, q-1) = 1.;
            }
            for (k1 = 0; k1 < q; k1++){
                mu1[k1] = 1./q;
                mu2[k1] = 1./q;
            }
            diff = 1.0;
            while (diff > epsilon){
                for (k1 = 0; k1 < q; k1++){
                    scra1[k1] = 0.0;
                    scra2[k1] = 0.0;
                }
                for (k1 = 0; k1 < q; k1++){
                    for (k2 = 0; k2 < q; k2++){
                        scra1[k1] += mu2[k2] * w(k1, k2);
                        scra2[k1] += mu1[k2] * w(k2, k1);
                    }
                }
                sum1 = 0.0;
                sum2 = 0.0;
                for (k1 = 0; k1 < q; k1++){
                    scra1[k1] = prob(i, k1) / scra1[k1];
                    sum1 += scra1[k1];
                    scra2[k1] = prob(j, k1) / scra2[k1];
                    sum2 += scra2[k1];
                }
                for (k1 = 0; k1 < q; k1++){
                    scra1[k1] /= sum1;
                    scra2[k1] /= sum2;
                }
                diff = -1.0;
                for (k1 = 0; k1 < q; k1++){
                    if (fabs(mu1[k1] - scra1[k1]) > diff)
                        diff = fabs(mu1[k1] - scra1[k1]);
                    if (fabs(mu2[k1] - scra2[k1]) > diff)
                        diff = fabs(mu2[k1] - scra2[k1]);
                    mu1[k1] = scra1[k1];
                    mu2[k1] = scra2[k1];
                }
            }

            sumpdir = 0.0;
            for (k1 = 0; k1 < q; k1++){
                for (k2 = 0; k2 < q; k2++){
                    w(k1,k2) = w(k1, k2) * mu1[k1] * mu2[k2];
                    sumpdir += w(k1,k2);
                }
            }

            sumdi = 0.0;
            for (k1 = 0; k1 < q; k1++){
                for (k2 = 0; k2 < q; k2++){
                    w(k1,k2) /= sumpdir;
                    sumdi += w(k1,k2) * log((w(k1,k2) + tiny) / (prob(i,k1) * prob(j,k2) +tiny));
                }
            }

            di(i,j) = di(j,i) = sumdi;
        }
    }

    #undef w
    #undef c
    #undef prob
    #undef di
    free(w);
    free(c);
    free(prob);
    free(mu1);
    free(mu2);
    free(scra1);
    free(scra2);
    return Py_BuildValue("O", diinfo);
}


static PyMethodDef msatools_methods[] = {

    {"msaentropy",  (PyCFunction)msaentropy,
     METH_VARARGS | METH_KEYWORDS,
     "Return an Shannon entropy array calculated for given character \n"
     "array that contains an MSA."},

    {"msamutinfo",  (PyCFunction)msamutinfo, METH_VARARGS | METH_KEYWORDS,
     "Return mutual information matrix calculated for given character \n"
     "array that contains an MSA."},

    {"msaocc",  (PyCFunction)msaocc, METH_VARARGS | METH_KEYWORDS,
     "Return occupancy (or count) array calculated for MSA rows or columns."},

    {"msaomes",  (PyCFunction)msaomes, METH_VARARGS | METH_KEYWORDS,
     "Return OMES matrix calculated for given character array that contains\n"
     "an MSA."},

    {"msasca",  (PyCFunction)msasca, METH_VARARGS | METH_KEYWORDS,
     "Return SCA matrix calculated for given character array that contains\n"
     "an MSA."},

    {"msameff",  (PyCFunction)msameff, METH_VARARGS | METH_KEYWORDS,
     "Return Meff calculated for given character array that contains\n"
     "an MSA."},

    {"msadirectinfo1",  (PyCFunction)msadirectinfo1, METH_VARARGS | METH_KEYWORDS,
     "Return DI correlation matrix to python for matrix inverse calculated\n"
     "for given character array that contains an MSA."},

    {"msadirectinfo2",  (PyCFunction)msadirectinfo2, METH_VARARGS | METH_KEYWORDS,
     "Return DI correlation matrix from inversed matrix."},

    {"msadipretest",  (PyCFunction)msadipretest, METH_VARARGS | METH_KEYWORDS,
     "Return some DI parameter to set array size."},

    {NULL, NULL, 0, NULL}
};



#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef msatools = {
        PyModuleDef_HEAD_INIT,
        "msatools",
        "Multiple sequence alignment analysis tools.",
        -1,
        msatools_methods,
};
PyMODINIT_FUNC PyInit_msatools(void) {
    import_array();
    return PyModule_Create(&msatools);
}
#else
PyMODINIT_FUNC initmsatools(void) {

    Py_InitModule3("msatools", msatools_methods,
        "Multiple sequence alignment analysis tools.");

    import_array();
}
#endif
