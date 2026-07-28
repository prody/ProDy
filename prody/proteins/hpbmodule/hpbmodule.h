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

#ifndef Py_HPBMODULE_H
#define Py_HPBMODULE_H

/////////////////////////////////////////////////////////////////////////
#include <Python.h>

////////////////////////////////////////////////////////////////////////
//
// Function Declarations
//
////////////////////////////////////////////////////////////////////////
//
// drivers in hpb module
//
static PyObject* hpb_hpb(PyObject* self, PyObject* args);

static PyObject* sasa_sasa(PyObject* self, PyObject* args);

//
// converter: python object -> data arrays for C++
//
void getAtomtypes( PyObject* L, int length, std::string* resIdx, std::string* atmIdx, 
                   double* radii, std::string* chains );
double* getCoords( PyObject* L, int len );

//
// functions to calculate hydrophobic, SASA and volume
//
PyObject* calcHPh( double* coordsA, int length, std::string* resIdx,
                   std::string* atmIdx, std::string* chains );
PyObject* calcSASA( double* coordsA, int length, std::string* resIdx,
                    std::string* atmIdx, double* radii, std::string* chains );
PyObject* calcVolume( double* coordsA, int length, std::string* resIdx,
                    std::string* atmIdx, double* radii, std::string* chains );

#endif
