/*
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

#ifndef FUN_H
#define FUN_H
#include <iostream>
#ifdef __cplusplus
extern"C" {
#endif
  extern void regtet_( double*, double*, double*, double*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
              int&, int&, int&, int&, const int&, const int&, const int&, double&, double&, double&, double&,
              int&, int&, int&, int&, int&, int&, int&, double&, int&, int&, int&, int&, int&, int&, int& );
#ifdef __cplusplus
}
#endif

void Delaunay_Tetrahedrize(int& n_v, double**& xyz, double*&weights, int& n_t, int**& tetra_node, int**&tetra_nbrs)
{
    const int nmax=150000, nvmax= 7*nmax, nhmax=1500;

    double* x = new double[nmax];
    double* y = new double[nmax];
    double* z = new double[nmax];
    double* w = new double[nmax];
    float* v = new float[nmax];
    int* ix = new int[nmax];
    int* iy = new int[nmax];
    int* iz = new int[nmax];
    int* iw = new int[nmax];
    int* ix2 = new int[nmax];
    int* iy2 = new int[nmax];
    int* iz2 = new int[nmax];
    int* iw2 = new int[nmax];
    int* is = new int[nmax];
    int* ifl = new int[nvmax];
    int* io = new int[nmax];
    int* id = new int[nvmax];
    int ih[nhmax];
    int nv, nw, nt, nd, naddl, isu, jsu, ksu, nsu, icfig, iwfig;
    double wlenx, wleny, wlenz, wlenw, epz;
    /*
     * ATTENTION!!!!!!
     *
     * BE CAREFUL ABOU THE BOOL AND LOGICAL
     * in old fortran logical(4) true is 1, false is 0;
     * how to convert a c++ bool or int to fortran logical?
     *
     */
    //bool delaun, pntoff, flphis, artfcl, random, reccor, redchk;
    int delaun, pntoff, flphis, artfcl, random, reccor, redchk;

    bool bigbox;
    int* icon = new int[nvmax*8];
    for(int i = 0; i < nvmax; ++ i)
        for(int j = 0; j < 8; ++ j)
          icon[i*8+j] = 0;

    nv = n_v;
    if(nv > nmax) {
        std::cerr << "Too many vertices for Delaunay tetrahedrization:" << nv << std::endl;
        return;
    }

    for(int i = 0; i < nv; ++ i) w[i] = weights[i];
    //if(w[0] == 0) delaun = true;//A Delaunay tetrahedralization (no weights) will be computed.
    if(w[0] == 0) delaun = 1;//A Delaunay tetrahedralization (no weights) will be computed.
    else delaun = 0;//A Regular tetrahedralization will be computed.
    //else delaun = false;//A Regular tetrahedralization will be computed.
    /*
     * ----------------------------------------------------------------------
     *
     *     pntoff =.true.
     *  Some input points will be inactive during tetrahedralization computation.
     */
    pntoff = 0;
    //pntoff = false;
     /*  All input points will be active during tetrahedralization computation.
     *
     * ----------------------------------------------------------------------
     *
     *  Note that in order to use flipping history, parameter nvmax in the program
     *  should be set to about 55*nmax, otherwise to about 7*nmax.
     *
     *     flphis=.true.
     *  Tetrahedron list with flipping history will be used.
     */
    //flphis = false;
    flphis = 0;
     /*  Shishkebab method will be used.
     *
     * ----------------------------------------------------------------------
     *
     *     artfcl=.true.
     *  Output tetrahedron list will include both real and artificial tetrahedra in
     *  the final tetrahedralization together with flipping history tetrahedra if
     *  flipping history was used for locating points.'
     */
    //artfcl = false;
    artfcl = 0;
     /*  Output tetrahedron list will only include real tetrahedra in the final
     *  tetrahedralization.
     *
     * ----------------------------------------------------------------------
     *
     *       random = .true.
     *  Input points will be inserted in a random order.
     *
     *       if(prompt)then
     *          read(9,*) isu, jsu, ksu, nsu
     *          write(*,*)'The four seeds for randomizing are:'
     *          write(*,*)isu, jsu, ksu, nsu
     *       else
     *          write(*,*)'Enter seeds for randomizing (4 integers):'
     *          read(5,*) isu, jsu, ksu, nsu
     */
    isu = 521288629;
    jsu = 362436069;
    ksu = 16163801;
    nsu = 131199299;
    random = 0;
    //random = false;
    /*
     *  Input points will be inserted in their current order.
     *
     * ----------------------------------------------------------------------
     *
     *       reccor=.true.
     *  Points that define a rectangular regular grid on the surface of a rectangular
     *  polyhedron that contains the set of input points will be part of the set for
     *  which tetrahedralization is to be computed; dimensions of polyhedron and
     *  choice of points that define grid are according to specifications provided by
     *  user.
     *
     *  If xmax, ymax, zmax are the maximum values of the x, y, and z coordinates of
     *  the input points, respectively, and xmin, ymin, zmin are the minimum values,
     *  then for positive numbers wlenx, wleny, wlenz (provided by user), the eight
     *  vertices of the the polyhedron will be: '
     *
     *  (xmin-wlenx, ymin-wleny, zmin-wlenz), (xmax+wlenx, ymin-wleny, zmin-wlenz),
     *  (xmax+wlenx, ymax+wleny  zmin-wlenz), (xmin-wlenx, ymax+wleny, zmin-wlenz),
     *  (xmin-wlenx, ymin-wleny, zmax+wlenz), (xmax+wlenx, ymin-wleny, zmax+wlenz),
     *  (xmax+wlenx, ymax+wleny  zmax+wlenz), (xmin-wlenx, ymax+wleny, zmax+wlenz).
     *
     *  For positive integer naddl (provided by user) for each facet of the
     *  polyhedron a set of naddl x naddl points is generated by program; this set
     *  defines a rectangular regular grid on the facet and contains the four
     *  vertices of the facet; the points in the union of the six sets thus generated
     *  define the rectangular grid on the surface of the polyhedron; naddl can not
     *  be less than 2; if it equals 2 then the grid is defined exactly by the 8
     *  vertices of the polyhedron.
     *
     *       if(.not.delaun) then
     *  If wmin is the minimum value of the weights of the input points then for a
     *  real number wlenw (provided by user) a weight equal to wmin - wlenw is
     *  assigned by the program to each point in the rectangular grid on the surface
     *  of the polyhedron.
     *       endif
     *
     *       if(prompt)then
     *          read(9,*) wlenx, wleny, wlenz
     *          if(.not.delaun) read(9,*) wlenw
     *          read(9,*) naddl
     *          write(*,*)'The values of wlenx, wleny, wlenz are:'
     *          write(*,*) wlenx, wleny, wlenz
     *          if(.not.delaun) then
     *             write(*,*)'The value of wlenw is:'
     *             write(*,*) wlenw
     *          endif
     *          write(*,*)'The value of naddl is:'
     *          write(*,*) naddl
     *       else
     *          write(*,*)'Enter wlenx, wleny, wlenz (3 positive real numbers):'
     *          read(5,*) wlenx, wleny, wlenz
     *          if(.not.delaun) then
     *             write(*,*)'Enter wlenw (a real number):'
     *             read(5,*) wlenw
     *          endif
     *          write(*,*)'Enter naddl (an integer greater than 1):'
     *          read(5,*) naddl
     *          write(10,*) wlenx, wleny, wlenz
     *          if(.not.delaun) write(10,*) wlenw
     *          write(10,*) naddl
     *       endif
     */
    reccor= 0;
    //reccor= false;
    /*
     * Points that define a rectangular regular grid on the surface of a rectangular
     * polyhedron that contains a set of input points will not be part of the set
     * for which a tetrahedralization is to be computed.
     *
     * ----------------------------------------------------------------------
     */
    //redchk = false;//interface_nbrs.f90 continue makes sense?
    redchk = 0;
    bigbox = 0;
    //bigbox = false;
    icfig = 9;
    //if(!delaun) iwfig = 9;
    if(delaun == 0) iwfig = 9;
    epz = 0.001;
    for(int i = 0; i < nv; ++ i) {
        x[i] = xyz[i][0];
        y[i] = xyz[i][1];
        z[i] = xyz[i][2];
        is[i] = 1;
    }
    
    regtet_(x, y, z, w, v, ix, iy, iz, iw, ix2, iy2, iz2, iw2,
                      icon, is, ifl, io, id, ih, nv, nw, nt, nd, nmax,
                      nvmax, nhmax, wlenx, wleny, wlenz, wlenw, naddl,
                      isu, jsu, ksu, nsu, icfig, iwfig, epz, delaun,
                      pntoff, flphis, artfcl, random, reccor, redchk);

    n_t = nt;

    tetra_node = new int*[nt];
    for(int i = 0; i < nt; ++ i) tetra_node[i] = new int[4];

    tetra_nbrs = new int*[nt];
    for(int i = 0; i < nt; ++ i) tetra_nbrs[i] = new int[4];

    for(int i = 0; i < nt; ++ i){
      tetra_node[i][0] = icon[i*8+4];
      tetra_node[i][1] = icon[i*8+5];
      tetra_node[i][2] = icon[i*8+6];
      tetra_node[i][3] = icon[i*8+7];
      tetra_nbrs[i][0] = icon[i*8+0];
      tetra_nbrs[i][1] = icon[i*8+1];
      tetra_nbrs[i][2] = icon[i*8+2];
      tetra_nbrs[i][3] = icon[i*8+3];
    }

  delete[] icon; delete[] x; delete[] y; delete[] z; delete[] v; delete[] ix; delete[] iy; delete[] iz;
  delete[] iw; delete[] ix2; delete[] iy2; delete[] iz2; delete[] iw2; delete[] is; delete[] ifl; delete[] id;
  icon = NULL; x = NULL; y = NULL; z = NULL; v = NULL; ix = NULL; iy = NULL; iz = NULL;
  iw = NULL; ix2 = NULL; iy2 = NULL; iz2 = NULL; iw2 = NULL; is = NULL; ifl = NULL; id = NULL;
}

#endif // FUN_H
