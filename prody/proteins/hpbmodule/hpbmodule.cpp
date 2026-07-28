/*
  hpb -- geometric methods for hydrophobic interactions and solvent accessible surface area
 
  Copyright (c) 2023, Xin Cao, Michelle H. Hummel, Bihua Yu, Evangelos A. Coutsias

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the “Software”), to deal
  in the Software without restriction, including without limitation the rights to
  use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
  the Software, and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

/////////////////////////////////////////////////////////////////////////
#define HPB_MODULE

#include "DT.h"
#include "find_vol_mod.h"
#include "hpbmodule.h"

/////////////////////////////////////////////////////////////////////////
//
// drivers in hpb module
//
/////////////////////////////////////////////////////////////////////////

static PyObject*
hpb_hpb(PyObject* self, PyObject* args)
{
    PyObject *listA; // the coordinates
    PyObject *listB; // the atom types and residues
	
    if ( ! PyArg_ParseTuple(args, "(OO)", &listA, &listB) ) {
	printf("Could not unparse objects\n");
        return NULL;
    }

    // read the two lists from python
    Py_INCREF(listA);
    Py_INCREF(listB);

    // test if selections are empty
    const int lenA = PyList_Size(listA);
    if ( lenA < 1 ) {
        printf("HPB ERROR: First selection didn't contain XYZ coordinates of any atoms. Please check!\n");
        // remove the lists in python
        Py_DECREF(listA);
        Py_DECREF(listB);
        return NULL;
    }
    const int lenB = PyList_Size(listB);
    if ( lenB < 1 ) {
        printf("HPB ERROR: Second selection didn't contain atom types. Please check!\n");
        // remove the lists in python
        Py_DECREF(listA);
        Py_DECREF(listB);
        return NULL;
    }
    
    // get the coodinates from the Python objects
    double* coordsA = getCoords( listA, lenA );
    
    // get the residue index and atom radii
    std::string* resIdx = new std::string[lenB];
    std::string* atmIdx = new std::string[lenB];
    std::string* chains = new std::string[lenB];
    double* radii = new double[lenB];
    getAtomtypes( listB, lenB, resIdx, atmIdx, radii, chains );

    // calculate the hydrophobic interactions, calcHPh
    PyObject* hydro;
    hydro = calcHPh( coordsA, lenA, resIdx, atmIdx, chains );

    // remove the lists in python
    Py_DECREF(listA);
    Py_DECREF(listB);
	
    // release memory
    delete[] coordsA;
    delete[] resIdx;
    delete[] atmIdx;
    delete[] radii;
    delete[] chains;

    return hydro;
}

static PyObject*
hpb_sasa(PyObject* self, PyObject* args)
{
    PyObject *listA; // the coordinates
    PyObject *listB; // the atom types and residues
    
    if ( ! PyArg_ParseTuple(args, "(OO)", &listA, &listB) ) {
        printf("Could not unparse objects\n");
        return NULL;
    }

    // read the two lists from python
    Py_INCREF(listA);
    Py_INCREF(listB);

    // test if selections are empty
    const int lenA = PyList_Size(listA);
    if ( lenA < 1 ) {
        printf("HPB ERROR: First selection didn't contain XYZ coordinates of any atoms. Please check!\n");
        // remove the lists in python
        Py_DECREF(listA);
        Py_DECREF(listB);
        return NULL;
    }
    const int lenB = PyList_Size(listB);
    if ( lenB < 1 ) {
        printf("HPB ERROR: Second selection didn't contain atom types. Please check!\n");
        // remove the lists in python
        Py_DECREF(listA);
        Py_DECREF(listB);
        return NULL;
    }
    
    // get the coodinates from the first list
    double* coordsA = getCoords( listA, lenA );

    // get the residue index and atom radii
    std::string* resIdx = new std::string[lenB];
    std::string* atmIdx = new std::string[lenB];
    std::string* chains = new std::string[lenB];
    double* radii = new double[lenB];
    getAtomtypes( listB, lenB, resIdx, atmIdx, radii, chains );
    
    // calculate the solvent accessible surface area, calcSASA
    PyObject* sasaV;
    sasaV = calcSASA( coordsA, lenA, resIdx, atmIdx, radii, chains );

    // remove the lists in python
    Py_DECREF(listA);
    Py_DECREF(listB);
    
    // release memory
    delete[] coordsA;
    delete[] resIdx;
    delete[] atmIdx;
    delete[] radii;
    delete[] chains;
    
    return sasaV;
}

static PyObject*
hpb_volume(PyObject* self, PyObject* args)
{
    PyObject *listA; // the coordinates
    PyObject *listB; // the atom types and residues

    if ( ! PyArg_ParseTuple(args, "(OO)", &listA, &listB) ) {
        printf("Could not unparse objects\n");
        return NULL;
    }

    // read the two lists from python
    Py_INCREF(listA);
    Py_INCREF(listB);

    // test if selections are empty
    const int lenA = PyList_Size(listA);
    if ( lenA < 1 ) {
        printf("HPB ERROR: First selection didn't contain XYZ coordinates of any atoms. Please check!\n");
        // remove the lists in python
        Py_DECREF(listA);
        Py_DECREF(listB);
        return NULL;
    }
    const int lenB = PyList_Size(listB);
    if ( lenB < 1 ) {
        printf("HPB ERROR: Second selection didn't contain atom types. Please check!\n");

        Py_DECREF(listA);
        Py_DECREF(listB);
        return NULL;
    }

    // get the coodinates from the first list
    double* coordsA = getCoords( listA, lenA );

    // get the residue index and atom radii
    std::string* resIdx = new std::string[lenB];
    std::string* atmIdx = new std::string[lenB];
    std::string* chains = new std::string[lenB];
    double* radii = new double[lenB];
    getAtomtypes( listB, lenB, resIdx, atmIdx, radii, chains );

    // calculate the volume, calcVolume
    PyObject* VolumeV;
    VolumeV = calcVolume( coordsA, lenA, resIdx, atmIdx, radii, chains );

    // remove the lists in python
    Py_DECREF(listA);
    Py_DECREF(listB);

    // release memory
    delete[] coordsA;
    delete[] resIdx;
    delete[] atmIdx;
    delete[] radii;
    delete[] chains;

    return VolumeV;
}

/////////////////////////////////////////////////////////////////////////
//
// converter: python object -> data arrays for C++
//
/////////////////////////////////////////////////////////////////////////

void getAtomtypes( PyObject* L, int length, std::string* resIdx, std::string* atmIdx, 
                   double* radii, std::string* chains )
{
    // loop through the arguments, pulling out the atom types and assign radii
    for ( int i = 0; i < length; i++ ) {
        PyObject* curAtom = PyList_GetItem(L,i);
        Py_INCREF(curAtom);

        PyObject* curVal = PyList_GetItem(curAtom,0);
        Py_INCREF(curVal);
#if PY_MAJOR_VERSION >= 3
        resIdx[i] = PyUnicode_AsUTF8(curVal);
#else
        resIdx[i] = PyString_AsString(curVal);
#endif
        Py_DECREF(curVal);
        
        curVal = PyList_GetItem(curAtom,1);
        Py_INCREF(curVal);
#if PY_MAJOR_VERSION >= 3
        atmIdx[i] = PyUnicode_AsUTF8(curVal);
#else
        atmIdx[i] = PyString_AsString(curVal);
#endif
        Py_DECREF(curVal);

        curVal = PyList_GetItem(curAtom,2);
        Py_INCREF(curVal);
#if PY_MAJOR_VERSION >= 3
        chains[i] = PyUnicode_AsUTF8(curVal);
#else
        chains[i] = PyString_AsString(curVal);
#endif
        Py_DECREF(curVal);
        
        if (atmIdx[i].substr(0,1) == "C") {
            radii[i] = 1.70+1.40;
        } else if (atmIdx[i].substr(0,1) == "N") {
            radii[i] = 1.55+1.40;
        } else if (atmIdx[i].substr(0,1) == "O") {
            radii[i] = 1.50+1.40;
        } else if (atmIdx[i].substr(0,1) == "S") {
            radii[i] = 1.80+1.40;
        } else if (atmIdx[i].substr(0,1) == "P") {
            radii[i] = 1.80+1.40;
        } else radii[i] = 0+1.40;

        Py_DECREF(curAtom);
    }
}

double* getCoords( PyObject* L, int length )
{
    // make space for the current coords
    double* coords = new double[3*length];

    // loop through the arguments, pulling out the XYZ coordinates.
    for ( int i = 0; i < length; i++ ) {
        PyObject* curCoord = PyList_GetItem(L,i);
        Py_INCREF(curCoord);

        PyObject* curVal = PyList_GetItem(curCoord,0);
        Py_INCREF(curVal);
        coords[3*i] = PyFloat_AsDouble(curVal);
        Py_DECREF(curVal);

        curVal = PyList_GetItem(curCoord,1);
        Py_INCREF(curVal);
        coords[3*i+1] = PyFloat_AsDouble(curVal);
        Py_DECREF(curVal);

        curVal = PyList_GetItem(curCoord,2);
        Py_INCREF(curVal);
        coords[3*i+2] = PyFloat_AsDouble(curVal);
        Py_DECREF(curVal);

        Py_DECREF(curCoord);
    }

    return coords;
}

/////////////////////////////////////////////////////////////////////////
//
// calcHPh
//
/////////////////////////////////////////////////////////////////////////

PyObject* calcHPh( double* coordsA, int length, std::string* resIdx, 
                   std::string* atmIdx, std::string* chains )
{
    PyObject* rVal = PyList_New(0);
    Py_INCREF(rVal);

    double* weights = (double*)malloc(sizeof(double)*length);
    double** coordc = new double*[length];
    
    for(int i = 0; i < length; ++ i) {
        coordc[i] = new double[3];
        weights[i] = 2.5*2.5; // initial radius for hydrophobic interactions
        for(int j = 0; j < 3; ++ j) {
            coordc[i][j] = coordsA[3*i+j];
        }
    }
    
    double* LIV_atom, *S_atom, *LIS_atom;
    LIV_atom = new double[length];
    S_atom = new double[length];
    LIS_atom = new double[length];

    for (int i = 0; i < length; ++ i) {
        LIV_atom[i] = 0;
        LIS_atom[i] = 0;
        S_atom[i] = 0;
    }
    
    double* allvalues = new double[9*length];
    int PR = 1;
    find_vol_intersection( length, coordc, weights, resIdx, allvalues, LIV_atom, S_atom, LIS_atom, PR );

    // save values of HPh interactions
    for (int i = 0; i < PR; i++) {
        int Idx1 = int(allvalues[3*i])-1;
        int Idx2 = int(allvalues[3*i+1])-1;
        std::string res1 = resIdx[ Idx1 ];
        std::string res2 = resIdx[ Idx2 ];
        std::string atom1 = atmIdx[ Idx1 ];
        std::string atom2 = atmIdx[ Idx2 ];
        std::string cha1 = chains[ Idx1 ];
        std::string cha2 = chains[ Idx2 ];
        
        PyObject* curPair = Py_BuildValue( "[s,s,s,s,s,s,d]", res1.c_str(), atom1.c_str(), cha1.c_str(), 
                                           res2.c_str(), atom2.c_str(), cha2.c_str(), allvalues[3*i+2] );
        Py_INCREF(curPair);
        
        PyList_Append(rVal, curPair);
    }

    for(int i = 0; i < length ; ++ i) delete[] coordc[i];
    coordc = NULL;
    delete[] allvalues;
    delete[] weights;
    delete[] LIV_atom;
    delete[] LIS_atom;
    delete[] S_atom;

    return rVal;
}

/////////////////////////////////////////////////////////////////////////
//
// calcSASA
//
/////////////////////////////////////////////////////////////////////////

PyObject* calcSASA( double* coordsA, int length, std::string* resIdx, 
                    std::string* atmIdx, double* radii, std::string* chains )
{
    PyObject* rVal = PyList_New(0);
    Py_INCREF(rVal);

    double* weights = (double*)malloc(sizeof(double)*length);
    double** coordc = new double*[length];
    
    for(int i = 0; i < length; ++ i) {
        coordc[i] = new double[3];
        weights[i] = radii[i]*radii[i];
        for(int j = 0; j < 3; ++ j) {
            coordc[i][j] = coordsA[3*i+j];
        }
    }
    
    double* LIV_atom, *S_atom, *LIS_atom;
    LIV_atom = new double[length];
    S_atom = new double[length];
    LIS_atom = new double[length];

    for (int i = 0; i < length; ++ i) {
        LIV_atom[i] = 0;
        LIS_atom[i] = 0;
        S_atom[i] = 0;
    }
    
    double* allvalues;
    int PR = 2;
    find_vol_intersection( length, coordc, weights, resIdx, allvalues, LIV_atom, S_atom, LIS_atom, PR );

    // save values of SASA
    double s_cur = 0;
    for (int i = 0; i < length-1; i++) {
        std::string res1 = resIdx[ i ];
        std::string res2 = resIdx[ i+1 ];
        
        s_cur = s_cur + S_atom[i];
        if (res1 != res2) {
            PyObject* curPair = Py_BuildValue( "[s,s,d]", res1.c_str(), chains[i].c_str(), s_cur );
            Py_INCREF(curPair);
            PyList_Append(rVal, curPair);
            s_cur = 0;
        }
        
        if (i == length-2) {
            s_cur = s_cur + S_atom[i+1];
            PyObject* curPair = Py_BuildValue( "[s,s,d]", res1.c_str(), chains[i].c_str(), s_cur );
            Py_INCREF(curPair);
            PyList_Append(rVal, curPair);
        }
    }

    for(int i = 0; i < length ; ++ i) delete[] coordc[i];
    coordc = NULL;
    delete[] allvalues;
    delete[] weights;
    delete[] LIV_atom;
    delete[] LIS_atom;
    delete[] S_atom;

    return rVal;
}

/////////////////////////////////////////////////////////////////////////
//
// calcVolume: enclosed by SASA
//
/////////////////////////////////////////////////////////////////////////

PyObject* calcVolume( double* coordsA, int length, std::string* resIdx,
                    std::string* atmIdx, double* radii, std::string* chains )
{
    PyObject* rVal = PyList_New(0);
    Py_INCREF(rVal);

    double* weights = (double*)malloc(sizeof(double)*length);
    double** coordc = new double*[length];

    for(int i = 0; i < length; ++ i) {
        coordc[i] = new double[3];
        weights[i] = radii[i]*radii[i];
        for(int j = 0; j < 3; ++ j) {
            coordc[i][j] = coordsA[3*i+j];
        }
    }

    double* LIV_atom, *S_atom, *LIS_atom;
    LIV_atom = new double[length];
    S_atom = new double[length];
    LIS_atom = new double[length];

    for (int i = 0; i < length; ++ i) {
        LIV_atom[i] = 0;
        LIS_atom[i] = 0;
        S_atom[i] = 0;
    }

    double* allvalues;
    int PR = 2;
    find_vol_intersection( length, coordc, weights, resIdx, allvalues, LIV_atom, S_atom, LIS_atom, PR );

    // save volume values
    double v_cur = 0;
    for (int i = 0; i < length-1; i++) {
        std::string res1 = resIdx[ i ];
        std::string res2 = resIdx[ i+1 ];

        v_cur = v_cur + LIV_atom[i];
        if (res1 != res2) {
            PyObject* curPair = Py_BuildValue( "[s,s,d]", res1.c_str(), chains[i].c_str(), v_cur );
            Py_INCREF(curPair);
            PyList_Append(rVal, curPair);
            v_cur = 0;
        }

        if (i == length-2) {
            v_cur = v_cur + LIV_atom[i+1];
            PyObject* curPair = Py_BuildValue( "[s,s,d]", res1.c_str(), chains[i].c_str(), v_cur );
            Py_INCREF(curPair);
            PyList_Append(rVal, curPair);
        }
    }

    for(int i = 0; i < length ; ++ i) delete[] coordc[i];
    coordc = NULL;
    delete[] allvalues;
    delete[] weights;
    delete[] LIV_atom;
    delete[] LIS_atom;
    delete[] S_atom;

    return rVal;
}

/////////////////////////////////////////////////////////////////////////
//
// Python setup for the module
//
/////////////////////////////////////////////////////////////////////////

static PyMethodDef hpbMethods[] = {
	{"hpb", hpb_hpb, METH_VARARGS, "Hydrophobic interactions."},
        {"sasa", hpb_sasa, METH_VARARGS, "Solvent accessible surface area."},
        {"volume", hpb_volume, METH_VARARGS, "Volume for residue."},
	{NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef hpb = {
        PyModuleDef_HEAD_INIT,
        "hpb",
        "Hydrophobic, SASA and Volume.",
        -1,
        hpbMethods,
};

PyMODINIT_FUNC PyInit_hpb(void) {
    return PyModule_Create(&hpb);
}
#else
PyMODINIT_FUNC inithpb(void)
{
	(void) Py_InitModule3("hpb", hpbMethods, "Hydrophobic, SASA and Volume.");
}
#endif
