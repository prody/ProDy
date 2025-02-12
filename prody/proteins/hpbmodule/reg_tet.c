/* reg_tet.f -- translated by f2c (version 20230428).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__9 = 9;
static integer c__30 = 30;

/* c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12 */
/* c */
/* c     Program regtet */
/* c */
/* c     This program will compute a Regular tetrahedralization for */
/* c     a set of points in 3-dimensional space. */
/* c */
/* c     A Regular tetrahedralization is the dual of a Power diagram. */
/* c     It is essentially a Delaunay tetrahedralization with weights. */
/* c */
/* c     If no weights are present the program will compute a Delaunay */
/* c     tetrahedralization. It is noted that if weights are present */
/* c     and degeneracies exist there may be points whose Power cells */
/* c     are not 3-dimensional. */
/* c */
/* c     The program is based on an algorithm for constructing Regular */
/* c     tetrahedralizations with incremental topological flipping. */
/* c */
/* c     Computations in this program for the purpose of computing */
/* c     the tetrahedralization are done in exact arithmetic whenever */
/* c     floating point arithmetic (done in double precision) does not */
/* c     seem appropriate. */
/* c */
/* c     Subroutine regtet is the driver routine of the program and is */
/* c     called from the main routine. */
/* c */
/* c     Documentation appears below in driver routine regtet. */
/* c */
/* c     Author: Javier Bernal */
/* c             National Institute of Standards and Technology (NIST) */
/* c */
/* c     Disclaimer: */
/* c */
/* c     This software was developed at the National Institute of Standards */
/* c     and Technology by employees of the Federal Government in the */
/* c     course of their official duties. Pursuant to title 17 Section 105 */
/* c     of the United States Code this software is not subject to */
/* c     copyright protection and is in the public domain. This software is */
/* c     experimental. NIST assumes no responsibility whatsoever for its */
/* c     use by other parties, and makes no guarantees, expressed or */
/* c     implied, about its quality, reliability, or any other */
/* c     characteristic. We would appreciate acknowledgement if the */
/* c     software is used. */
/* c */
/* *MAIN */
/* c */
/*      program main */
/* c */
/* c     Setting of parameters: */
/* c */
/* c     1. Flipping history used for locating points: */
/* c */
/* c     integer nmax, nvmax, nhmax */
/* c     parameter (nmax=150000, nvmax=55*nmax, nhmax=1500) */
/* c */
/* c     2. Flipping history not used for locating points: */
/* c */
/*      integer nmax, nvmax, nhmax */
/*      parameter (nmax=150000, nvmax= 7*nmax, nhmax=1500) */
/* c */
/*      double precision x(nmax), y(nmax), z(nmax), w(nmax) */
/*      real v(nmax) */
/*      integer ix(nmax), iy(nmax), iz(nmax), iw(nmax) */
/*      integer ix2(nmax), iy2(nmax), iz2(nmax), iw2(nmax) */
/*      integer icon(8,nvmax), is(nmax), ifl(nvmax), io(nmax) */
/*      integer id(nvmax), ih(nhmax) */
/*      integer nv, nw, nt, nd, naddl, isu, jsu, ksu, nsu, icfig, iwfig */
/*      double precision wlenx, wleny, wlenz, wlenw, epz */
/*      logical delaun, pntoff, flphis, artfcl, random, reccor, redchk */
/* c */
/*      logical prompt, bigbox */
/*      character*1 answ */
/*      double precision xcor, ycor, zcor, wght */
/*      integer np, i, j, iric, irec */
/*      integer ideli, ipnti, iflpi, iarti, irani, ireci, iredi */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      write(*,*)'This program is for computing a Regular ', */
/*     *          'tetrahedralization for a set' */
/*      write(*,*)'of points in 3-d space, ', */
/*     *          'i. e. a Delaunay tetrahedralization with weights.' */
/*      write(*,*)' ' */
/*      write(*,*)'If no weights are present a Delaunay ', */
/*     *          'tetrahedralization is computed.' */
/*      write(*,*)'Note that if weights are present and degeneracies ', */
/*     *          'exist there may' */
/*      write(*,*)'be points whose Power cells are ', */
/*     *          'not 3-dimensional.' */
/*      write(*,*)' ' */
/*      write(*,*)'Computations in this program are done in exact ', */
/*     *          'arithmetic whenever' */
/*      write(*,*)'floating point arithmetic (done in double precision) ', */
/*     *          'does not seem' */
/*      write(*,*)'appropriate.' */
/*      write(*,*)' ' */
/*      write(*,*)'Documentation about input/output variables, etc. ', */
/*     *          'appears in driver' */
/*      write(*,*)'routine regtet' */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/* c     inquire about existence of prompts file */
/* c */
/*      inquire(file = 'prompts.data', exist = prompt) */
/*      write(*,*)' ' */
/*      if(prompt)then */
/*         write(*,*)'User intervention not required: responses to ', */
/*     *             'would-be program generated' */
/*         write(*,*)'prompts will be read from file prompts.data ', */
/*     *             '(prompts will be suppressed).' */
/*         open (9, file = 'prompts.data') */
/*      else */
/*         write(*,*)'User intervention required: file named ', */
/*     *             'prompts.data does not exist, thus' */
/*         write(*,*)'program will generate prompts for which ', */
/*     *             'user must provide responses (in' */
/*         write(*,*)'order to avoid user intervention a file ', */
/*     *             'named prompts.data must exist that' */
/*         write(*,*)'contains the responses to the would-be ', */
/*     *             'prompts; prompts are then suppressed;' */
/*         write(*,*)'for this purpose the program will generate ', */
/*     *             'a file named prompts.tentative' */
/*         write(*,*)'that will contain the responses provided ', */
/*     *             'by the user during the current run,' */
/*         write(*,*)'and which can then be renamed prompts.data ', */
/*     *             ' by the user for future runs).' */
/*         open (10, file = 'prompts.tentative') */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*   10 format (a1) */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'A Delaunay tetrahedralization (no weights) to be ', */
/*     *             'computed?(y/n)' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         delaun=.true. */
/*         write(*,*)'A Delaunay tetrahedralization will be computed.' */
/*      else */
/*         delaun=.false. */
/*         write(*,*)'A Regular tetrahedralization will be computed.' */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'Any input points to be inactive during ', */
/*     *             'tetrahedralization computation?(y/n)' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         pntoff =.true. */
/*         write(*,*)'Some input points will be inactive during ', */
/*     *             'tetrahedralization computation.' */
/*      else */
/*         pntoff =.false. */
/*         write(*,*)'All input points will be active during ', */
/*     *             'tetrahedralization computation.' */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'Tetrahedron list with flipping history to be ', */
/*     *             'used for locating points?(y/n)' */
/*         write(*,*)'(Otherwise shishkebab method will be used. Note ', */
/*     *             'that in order' */
/*         write(*,*)'to use flipping history, parameter nvmax in the ', */
/*     *             'program should' */
/*         write(*,*)'be set to about 55*nmax, otherwise to about ', */
/*     *             '7*nmax).' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         flphis=.true. */
/*         write(*,*)'Tetrahedron list with flipping history will be ', */
/*     *             'used.' */
/*      else */
/*         flphis=.false. */
/*         write(*,*)'Shishkebab method will be used.' */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'Do you want output tetrahedron list to include ', */
/*     *             'both real and artificial' */
/*         write(*,*)'tetrahedra in the final tetrahedralization?(y/n)' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         artfcl=.true. */
/*         write(*,*)'Output tetrahedron list will include both real ', */
/*     *             'and artificial tetrahedra' */
/*         write(*,*)'in the final tetrahedralization together with ', */
/*     *             'flipping history tetrahedra' */
/*         write(*,*)'if flipping history was used for locating points.' */
/*      else */
/*         artfcl=.false. */
/*         write(*,*)'Output tetrahedron list will only include ', */
/*     *             'real tetrahedra in the final' */
/*         write(*,*)'tetrahedralization.' */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'Do you want input points to be inserted in ', */
/*     *             'a random order?(y/n)' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         random = .true. */
/*         write(*,*)'Input points will be inserted in a random order.' */
/*         write(*,*)' ' */
/*         if(prompt)then */
/*            read(9,*) isu, jsu, ksu, nsu */
/*            write(*,*)'The four seeds for randomizing are:' */
/*            write(*,*)isu, jsu, ksu, nsu */
/*         else */
/*            write(*,*)'Enter seeds for randomizing (4 integers):' */
/*            read(5,*) isu, jsu, ksu, nsu */
/* c */
/* c           isu = 521288629 */
/* c           jsu = 362436069 */
/* c           ksu = 16163801 */
/* c           nsu = 131199299 */
/* c */
/*            write(10,*) isu, jsu, ksu, nsu */
/*         endif */
/*      else */
/*         random = .false. */
/*         write(*,*)'Input points will be inserted in their ', */
/*     *             'current order.' */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'Do you want points that define a rectangular ', */
/*     *             'regular grid on the surface' */
/*         write(*,*)'of a rectangular polyhedron that contains ', */
/*     *             'set of input points to become,' */
/*         write(*,*)'together with set of input points and ', */
/*     *             ' artificial points, the set for' */
/*         write(*,*)'which a tetrahedralization is to be computed?(y/n)' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         reccor=.true. */
/*         write(*,*)'Points that define a rectangular regular grid on ', */
/*     *             'the surface of a' */
/*         write(*,*)'rectangular polyhedron that contains set of input ', */
/*     *             'points will be part' */
/*         write(*,*)'of the set for which tetrahedralization is to be ', */
/*     *             'computed; dimensions' */
/*         write(*,*)'of polyhedron and choice of points that define ', */
/*     *             'grid are according to' */
/*         write(*,*)'specifications provided by user.' */
/*         write(*,*)' ' */
/*         write(*,*)'If xmax, ymax, zmax are the maximum values of the ', */
/*     *             'x, y, and z coordinates' */
/*         write(*,*)'of the input points, respectively, and xmin, ', */
/*     *             'ymin, zmin are the minimum' */
/*         write(*,*)'values, then for positive numbers ', */
/*     *             'wlenx, wleny, wlenz (provided by user),' */
/*         write(*,*)'the eight vertices of the the polyhedron will be: ' */
/*         write(*,*)'(xmin-wlenx, ymin-wleny, zmin-wlenz), ', */
/*     *             '(xmax+wlenx, ymin-wleny, zmin-wlenz),' */
/*         write(*,*)'(xmax+wlenx, ymax+wleny  zmin-wlenz), ', */
/*     *             '(xmin-wlenx, ymax+wleny, zmin-wlenz),' */
/*         write(*,*)'(xmin-wlenx, ymin-wleny, zmax+wlenz), ', */
/*     *             '(xmax+wlenx, ymin-wleny, zmax+wlenz),' */
/*         write(*,*)'(xmax+wlenx, ymax+wleny  zmax+wlenz), ', */
/*     *             '(xmin-wlenx, ymax+wleny, zmax+wlenz).' */
/*         write(*,*)' ' */
/*         write(*,*)'For positive integer naddl (provided by user) ', */
/*     *             'for each facet of the' */
/*         write(*,*)'polyhedron a set of naddl x naddl points is ', */
/*     *             'generated by program;' */
/*         write(*,*)'this set defines a rectangular regular grid ', */
/*     *             'on the facet and contains' */
/*         write(*,*)'the four vertices of the facet; the points in ', */
/*     *             'the union of the six' */
/*         write(*,*)'sets thus generated define the rectangular ', */
/*     *             'grid on the surface of' */
/*         write(*,*)'the polyhedron; naddl can not be less than 2; ', */
/*     *             'if it equals 2 then' */
/*         write(*,*)'the grid is defined exactly by the 8 vertices ', */
/*     *             'of the polyhedron.' */
/*         if(.not.delaun) then */
/*            write(*,*)' ' */
/*            write(*,*)'If wmin is the minimum value of the weights of ', */
/*     *                'the input points then' */
/*            write(*,*)'for a real number wlenw (provided ', */
/*     *                'by user) a weight equal to wmin - wlenw' */
/*            write(*,*)'is assigned by the program to ', */
/*     *                'each point in the rectangular grid on the' */
/*            write(*,*)'surface of the polyhedron.' */
/*         endif */
/*         write(*,*)' ' */
/*         if(prompt)then */
/*            read(9,*) wlenx, wleny, wlenz */
/*            if(.not.delaun) read(9,*) wlenw */
/*            read(9,*) naddl */
/*            write(*,*)'The values of wlenx, wleny, wlenz are:' */
/*            write(*,*) wlenx, wleny, wlenz */
/*            if(.not.delaun) then */
/*               write(*,*)'The value of wlenw is:' */
/*               write(*,*) wlenw */
/*            endif */
/*            write(*,*)'The value of naddl is:' */
/*            write(*,*) naddl */
/*         else */
/*            write(*,*)'Enter wlenx, wleny, wlenz ', */
/*     *      '(3 positive real numbers):' */
/*            read(5,*) wlenx, wleny, wlenz */
/*            if(.not.delaun) then */
/*               write(*,*)'Enter wlenw (a real number):' */
/*               read(5,*) wlenw */
/*            endif */
/*            write(*,*)'Enter naddl (an integer greater than 1):' */
/*            read(5,*) naddl */
/*            write(10,*) wlenx, wleny, wlenz */
/*            if(.not.delaun) write(10,*) wlenw */
/*            write(10,*) naddl */
/*         endif */
/*      else */
/*         reccor=.false. */
/*         write(*,*)'Points that define a rectangular regular grid on ', */
/*     *             'the surface of a' */
/*         write(*,*)'rectangular polyhedron that ', */
/*     *             'contains set of input points will not' */
/*         write(*,*)'be part of the set for which a tetrahedralization ', */
/*     *             'is to be computed.' */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      if(delaun) go to 50 */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'Do you want redundant points to be checked for ', */
/*     *             'redundancy after ' */
/*         write(*,*)'final tetrahedralization has been computed?(y/n)' */
/*         write(*,*)'(doing so may require some additional time).' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         redchk=.true. */
/*         write(*,*)'Redundant points will be checked for redundancy.' */
/*      else */
/*         redchk=.false. */
/*         write(*,*)'Redundant points will not be checked for ', */
/*     *             'redundancy.' */
/*      endif */
/*   50 continue */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      bigbox=.false. */
/*      if(.not.reccor) go to 55 */
/*      write(*,*)' ' */
/*      answ = ' ' */
/*      if(prompt)then */
/*         read(9,10) answ */
/*      else */
/*         write(*,*)'Do you want points that define rectangular grid ', */
/*     *             'as described above' */
/*         write(*,*)'to be saved in a file?(y/n)' */
/*         read(5,10) answ */
/*         write(10,10) answ */
/*      endif */
/*      if(answ.eq.'y'.or.answ.eq.'Y') then */
/*         bigbox=.true. */
/*         write(*,*)'Rectangular grid will be saved in a file.' */
/*      else */
/*         write(*,*)'Rectangular grid will not be saved in a file.' */
/*      endif */
/*   55 continue */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      if(prompt)then */
/*         read(9,*) icfig */
/*         write(*,*)'The number of significant figures of decimal part ', */
/*     *             'of point coordinates is ', icfig */
/*      else */
/*         write(*,*)'Enter icfig, i. e. number of significant figures ', */
/*     *             'of decimal part' */
/*         write(*,*)'of point coordinates, -1 < icfig < 10 (total ', */
/*     *             'number of significant' */
/*         write(*,*)'figures should be at most 14 with at most 9 to ', */
/*     *             'either the left or' */
/*         write(*,*)'the right of the decimal point):' */
/*         read(5,*) icfig */
/*         write(10,*) icfig */
/*      endif */
/*      if(delaun) go to 60 */
/*      write(*,*)' ' */
/*      if(prompt)then */
/*         read(9,*) iwfig */
/*         write(*,*)'The number of significant figures of decimal part ', */
/*     *             'of point weights is ', iwfig */
/*      else */
/*         write(*,*)'Enter iwfig, i. e. number of significant figures ', */
/*     *             'of decimal part' */
/*         write(*,*)'of point weights, -1 < iwfig < 10, ', */
/*     *             '-1 < 2*icfig - iwfig < 10 (total' */
/*         write(*,*)'number of significant figures should ', */
/*     *             'be at most 14 with at most 9' */
/*         write(*,*)'to either the left or the right of the ', */
/*     *             'decimal point):' */
/*         read(5,*) iwfig */
/*         write(10,*) iwfig */
/*      endif */
/*   60 continue */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/* c     open files */
/* c */
/*      open (11, file = 'pnts-wghts') */
/*      open (12, file = 'tetrahedra') */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/* c     set tolerance */
/* c */
/*      epz = 0.001d0 */
/* c */
/* c     read vertex data */
/* c */
/*      write(*,*)' ' */
/*      write(*,*)'Reading vertex data ...' */
/*      nv = 0 */
/*      if(delaun) go to 130 */
/* c */
/* c     read vertex data with weights */
/* c */
/*  100 continue */
/*      read (11, *, end = 120) xcor, ycor, zcor, wght */
/*      nv = nv + 1 */
/*      if (nv .gt. nmax) stop 10 */
/*      x(nv) = xcor */
/*      y(nv) = ycor */
/*      z(nv) = zcor */
/*      w(nv) = wght */
/*      go to 100 */
/*  120 continue */
/*      go to 140 */
/* c */
/* c     read vertex data without weights */
/* c */
/*  130 continue */
/*      read (11, *, end = 140) xcor, ycor, zcor */
/*      nv = nv + 1 */
/*      if (nv .gt. nmax) stop 20 */
/*      x(nv) = xcor */
/*      y(nv) = ycor */
/*      z(nv) = zcor */
/*      w(nv) = 0.0d0 */
/*      go to 130 */
/*  140 continue */
/* c */
/* c     read off-on information about input points if there are any that */
/* c     are to be inactive during the tetrahedralization computation */
/* c */
/*      if(pntoff) then */
/*         open (13, file = 'points-off') */
/*         read (13,*) np */
/*         if(np .ne. nv) stop 30 */
/*         read (13,150) (is(i), i = 1, np) */
/*      endif */
/*  150 format (40(1x,i1)) */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      write(*,*)'Computation of tetrahedralization to begin: driver ', */
/*     *          'subroutine regtet' */
/*      write(*,*)'will be called from the ', */
/*     *          'main routine and points will be processed.' */
/*      write(*,*)' ' */
/*      write(*,*)'Please wait ...' */
/*      write(*,*)' ' */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/* c     call regtet to compute tetrahedralization */
/* c */
/*      call regtet(x, y, z, w, v, ix, iy, iz, iw, ix2, iy2, iz2, iw2, */
/*     *            icon, is, ifl, io, id, ih, nv, nw, nt, nd, nmax, */
/*     *            nvmax, nhmax, wlenx, wleny, wlenz, wlenw, naddl, */
/*     *            isu, jsu, ksu, nsu, icfig, iwfig, epz, delaun, */
/*     *            pntoff, flphis, artfcl, random, reccor, redchk) */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      write(*,*)'(Back to the main routine).' */
/*      write(*,*)' ' */
/*      write(*,*)'Computation of tetrahedralization has been completed.' */
/*      write(*,*)' ' */
/*      write(*,*)'Number of input vertices = ',nv */
/*      write(*,*)' ' */
/*      write(*,*)'Number of output vertices = ',nw */
/*      write(*,*)' ' */
/*      write(*,*)'Length of final tetrahedron list = ',nt */
/*      write(*,*)' ' */
/*      write(*,*)'Number of real tetrahedra = ',nd */
/*      write(*,*)' ' */
/*      write(*,*)'(The output vertices are the vertices of tetrahedra ', */
/*     *          'in the final' */
/*      write(*,*)'tetrahedron list).' */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/* c     write tetrahedralization information in files */
/* c */
/*      ideli = 0 */
/*      ipnti = 0 */
/*      iflpi = 0 */
/*      iarti = 0 */
/*      irani = 0 */
/*      ireci = 0 */
/*      iredi = 0 */
/*      if(delaun) ideli = 1 */
/*      if(pntoff) ipnti = 1 */
/*      if(flphis) iflpi = 1 */
/*      if(artfcl) iarti = 1 */
/*      if(random) irani = 1 */
/*      if(reccor) ireci = 1 */
/*      if(redchk) iredi = 1 */
/*      write(*,*)' ' */
/*      write(*,*)'Saving tetrahedralization data in output file ...' */
/*      if(nt.ne.0)then */
/*         write (12,480) ideli, ipnti, iflpi, iarti, irani, ireci, iredi */
/*         if(.not.delaun) then */
/*            write (12,*) nw, nt, icfig, iwfig */
/*         else */
/*            write (12,*) nw, nt, icfig */
/*         endif */
/*         write (12,500) ((icon(i,j), i = 1, 8), j = 1, nt) */
/*         write (12,550) (is(i), i = 1, nw) */
/*         if(reccor .and. .not.delaun) then */
/*            write (12,*) wlenx, wleny, wlenz, wlenw */
/*            write (12,*) naddl */
/*         elseif(reccor) then */
/*            write (12,*) wlenx, wleny, wlenz */
/*            write (12,*) naddl */
/*         endif */
/*      else */
/*         write (*,*)'warning: no real tetrahedra were created.' */
/*         write (*,*)' ' */
/*         write (12,*)'warning: no real tetrahedra were created.' */
/*         write (12,*)' ' */
/*      endif */
/* c */
/*  480 format (7(1x,i1)) */
/*  500 format (8i10) */
/*  550 format (7i10) */
/* c */
/*      if(bigbox) then */
/*         open (14, file = 'bgbox-pnts') */
/*         iric = 1 */
/*         if(artfcl) iric = 9 */
/*         irec = iric + 6*(naddl**2) - 12*naddl + 7 */
/*         if(delaun) then */
/*            do 580 i = iric, irec */
/*               write (14,*) x(i), y(i), z(i) */
/*  580       continue */
/*         else */
/*            do 590 i = iric, irec */
/*               write (14,*) x(i), y(i), z(i), w(i) */
/*  590       continue */
/*         endif */
/*      endif */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/* c     write information in file in a readable form */
/* c */
/* c     open (23, file = 'exeinfor') */
/* c     write (23,*)' ' */
/* c     write (23,*)'Number of vertices   =',nw */
/* c     write (23,*)'Number of tetrahedra =',nt */
/* c     write (23,*)' ' */
/* c     if(nt.ne.0)then */
/* c        do 600 i = 1, nt */
/* c           write (23,800) i, (icon(j,i), j=1,8) */
/* c 600    continue */
/* c        do 650 i = 1, nw */
/* c           write (23,850) i, is(i) */
/* c 650    continue */
/* c     else */
/* c        write (23,*)'warning: no real tetrahedra were created.' */
/* c        write (23,*)' ' */
/* c     endif */
/* c 800 format ('tetr: ', i4, ' info: ', 8(1x,i4)) */
/* c 850 format ('vert: ', i4, ' is: ',i4) */
/* c */
/* c---------------------------------------------------------------------- */
/* c */
/*      write(*,*)' ' */
/*      stop */
/*      end */
/* *REGTET */
/* ********************************************************************** */

/*     Driver subroutine of Fortran 77 program REGTET [3] for computing */
/*     a Regular tetrahedralization for a set of points in 3-dimensional */
/*     space. */

/*     A Regular tetrahedralization is the dual of a Power diagram. */
/*     It is essentially a Delaunay tetrahedralization with weights. */
/*     In the absence of weights the program simply computes a Delaunay */
/*     tetrahedralization. It is noted that if weights are present */
/*     and degeneracies exist there may be points whose Power cells */
/*     are not 3-dimensional. */

/*     Computations in this program for the purpose of computing */
/*     the tetrahedralization are done in exact arithmetic whenever */
/*     floating point arithmetic (done in double precision) does not */
/*     seem appropriate. */

/*     The program is based on an algorithm by Edelsbrunner and Shah [1] */
/*     for constructing Regular tetrahedralizations with incremental */
/*     topological flipping. At the start of the execution of the */
/*     program a Regular tetrahedralization for the vertices of an */
/*     artificial cube that contains the input points is constructed. */
/*     Throughout the execution the vertices of this cube (artificial */
/*     points) are treated in the proper lexicographical manner [2] so */
/*     that the final tetrahedralization is correct. The program has the */
/*     capability of maintaining at all times during the execution of */
/*     the program a list of all tetrahedra in the current and previous */
/*     tetrahedralizations. This list is in the form of a directed */
/*     acyclic graph that represents the history of the flips the */
/*     program has performed, and it is used by the program to identify */
/*     a tetrahedron in the current tetrahedralization that contains a */
/*     new input point. Finally, the program has the capability of */
/*     adding the input points in a random sequence. */

/*     [1] H. Edelsbrunner and N. R. Shah, Incremental topological */
/*         flipping works for regular triangulations, Algorithmica 15(3), */
/*         223-241 (1996). */

/*     [2] J. Bernal, Lexicographical manipulations for correctly */
/*         computing regular tetrahedralizations with incremental */
/*         topological flipping, NISTIR 6335 (1999). */

/*     [3] J. Bernal, REGTET: A program for computing regular */
/*         tetrahedralizations (long version), NISTIR 6786 (2001). */

/*         Author: Javier Bernal */

/* ********************************************************************** */

/*     The following are examples of how parameters that are dimensions */
/*     of arrays can be defined in main routine: */

/*     1. Flipping history used for locating points: */

/*     integer nmax, nvmax, nhmax */
/*     parameter (nmax=150000, nvmax=55*nmax, nhmax=1500) */

/*     2. Flipping history not used for locating points: */

/*     integer nmax, nvmax, nhmax */
/*     parameter (nmax=150000, nvmax= 7*nmax, nhmax=1500) */

/*     Arrays and logical variables should be defined in main routine */
/*     as follows: */

/*     double precision x(nmax), y(nmax), z(nmax), w(nmax) */
/*     real v(nmax) */
/*     integer ix(nmax), iy(nmax), iz(nmax), iw(nmax) */
/*     integer ix2(nmax), iy2(nmax), iz2(nmax), iw2(nmax) */
/*     integer icon(8,nvmax), is(nmax), ifl(nvmax), io(nmax) */
/*     integer id(nvmax), ih(nhmax) */
/*     integer nv, nw, nt, nd, naddl, isu, jsu, ksu, nsu, icfig, iwfig */
/*     double precision wlenx, wleny, wlenz, wlenw, epz */
/*     logical delaun, pntoff, flphis, artfcl, random, reccor, redchk */

/*     Subroutine regtet should be called in main routine as follows: */

/*     call regtet(x, y, z, w, v, ix, iy, iz, iw, ix2, iy2, iz2, iw2, */
/*    *            icon, is, ifl, io, id, ih, nv, nw, nt, nd, nmax, */
/*    *            nvmax, nhmax, wlenx, wleny, wlenz, wlenw, naddl, */
/*    *            isu, jsu, ksu, nsu, icfig, iwfig, epz, delaun, */
/*    *            pntoff, flphis, artfcl, random, reccor, redchk) */

/* ********************************************************************** */

/*     Description of variables in calling sequence */
/*     -------------------------------------------- */

/*     For the purpose of describing these variables six sets of points */
/*     are first described. It should be noted that the descriptions of */
/*     some of these sets may refer to variables that have not yet been */
/*     described. */

/*     SI     Set of points on input, i. e. initial set of points */
/*            provided by the user. Input integer variable nv (described */
/*            below) equals the number of points in this set. */

/*     SP     Set of points that define a rectangular regular grid on */
/*            the surface of a rectangular polyhedron that contains SI */
/*            in its interior and which are internally generated by the */
/*            program if logical variable reccor (described below) is */
/*            .true.. This polyhedron is big enough to contain SI and */
/*            its size and points in SP are determined by the values */
/*            of variables wlenx, wleny, wlenz and naddl (described */
/*            below) provided by the user. If logical variable delaun */
/*            (described below) is .false. and if wmin is the minimum */
/*            weight of points in SI then a weight equal to wmin - wlenw */
/*            is assigned by the program to each point in SP where the */
/*            value of variable wlenw (described below) is provided by */
/*            the user. */

/*     SA     Set of eight points which are the vertices of an artificial */
/*            cube. These points are essentially at infinity and must be */
/*            treated by the program in the proper lexicographical */
/*            manner. These points are the artificial points. */

/*     SD     Set of points for which a Regular tetrahedralization is */
/*            desired. Set is defined according to the value of logical */
/*            variable reccor (described below) as follows: */
/*            if reccor is .false. then SD equals SI; */
/*            if reccor is .true. then SD equals the union of SP and SI. */
/*            If a point in SD has a weight as a member of either SP or */
/*            SI, it gets the same weight in SD. */
/*            The points in SD are the real points. */

/*     SU     Set of points for which the program actually computes a */
/*            Regular tetrahedralization. */
/*            Set equals the union of SA and SD. */
/*            If a point in SU is in either SP or SI and has a weight as */
/*            a member of either set, it gets the same weight in SU. */
/*            The final tetrahedralization computed by the program is */
/*            Regular for SU and contains a Regular tetrahedralization */
/*            for SD. */

/*     SO     Set of points on output. Points in this set are the */
/*            vertices of tetrahedra in output tetrahedron list. Output */
/*            integer variable nw (described below) equals the number of */
/*            points in this set. Set is defined according to the value */
/*            of logical variable artfcl (described below) as follows: */
/*            if artfcl is .false. then SO equals SD; */
/*            if artfcl is .true. then SO equals SU. */
/*            (As described above SD and SU are defined according to the */
/*            value of logical variable reccor (described below)). */
/*            If a point in SO is in either SP or SI and has a weight as */
/*            a member of either set, it gets the same weight in SO. */

/*            Depending on how SO is defined points in SO are ordered as */
/*            follows (here ng is the number of points in SP): */

/*            case 1: SO equals SI: for each i, i = 1, ..., nv, the */
/*            ith point in SI is the ith point in SO; */
/*            case 2: SO equals the union of SP and SI: for each i, */
/*            i = 1, ..., ng, the ith point in SP is the ith point in SO, */
/*            and for each i, i = 1, ..., nv, the ith point in SI is the */
/*            (i+ng)th point in SO; */
/*            case 3: SO equals the union of SA and SI: for each i, */
/*            i = 1, ..., 8, the ith point in SA is the ith point in SO, */
/*            and for each i, i = 1, ..., nv, the ith point in SI is the */
/*            (i+8)th point in SO; */
/*            case 4: SO equals the union of SA, SP, and SI: for each i, */
/*            i = 1, ..., 8, the ith point in SA is the ith point in SO, */
/*            for each i, i = 1, ..., ng, the ith point in SP is the */
/*            (i+8)th point in SO, and for each i, i = 1, ..., nv, the */
/*            ith point in SI is the (i+ng+8)th point in SO. */

/*            In what follows for each i, i = 1, ..., nw, the ith point */
/*            in SO is also refered to as point i. */

/*     The description of the variables follows: */

/*     delaun Input logical flag; */
/*            if .true. (no weights) then a Delaunay tetrahedralization */
/*            is to be computed; */
/*            if .false. then a Regular tetrahedralization is to be */
/*            computed. */

/*     pntoff Input logical flag; */
/*            if .true. then some input points will be inactive during */
/*            tetrahedralization computation; */
/*            if .false. then all input points will be active during */
/*            tetrahedralization computation. */

/*     flphis Input logical flag; */
/*            if .true. then tetrahedron list with flipping history will */
/*            be used for locating points; */
/*            if .false. then tetrahedron list with flipping history will */
/*            not be used for locating points; a shishkebab method that */
/*            locates points by checking tetrahedra in the current */
/*            tetrahedralization will be used. */

/*            The use of the tetrahedron list with flipping history */
/*            requires a lot more space than the alternative; parameter */
/*            nvmax (described below) must be set to about 55 times */
/*            the value of parameter nmax (described below) if flphis */
/*            equals .true.; otherwise it must be set to about 7 times */
/*            the value of parameter nmax. */

/*     artfcl Input logical flag; */
/*            if .true. and if logical variable flphis (described above) */
/*            is set to .true. then output tetrahedron list will include */
/*            the final tetrahedralization for SO (equal to SU) together */
/*            with the flipping history tetrahedra, i.e. tetrahedra that */
/*            at some point during the execution of the program were part */
/*            of a tetrahedralization for SO but that are not in the */
/*            final tetrahedralization; */
/*            if .true. and if flphis equals .false. then output */
/*            tetrahedron list will just be the final tetrahedralization */
/*            for SO (equal to SU). */
/*            if .false. then output tetrahedron list will only include */
/*            the final tetrahedralization for SO (equal to SD)). */

/*     random Input logical flag; */
/*            if .true. then the points in SI are to be added by the */
/*            program in a random fashion; */
/*            if .false. then the points in SI are to be added by the */
/*            program in their original order. */

/*            If the points in SI are already randomized on input then */
/*            there is no need for the program to randomize them again */
/*            so that logical variable random should be set to .false.. */
/*            If the points are in a nice order and the shishkebab */
/*            method is to be used for locating points then random */
/*            should be set to .false.. */

/*     reccor Input logical flag; */
/*            if .true. then SD includes SP; */
/*            if .false. then SD does not include SP. */

/*            Including SP in SD is recommended if it is not desirable */
/*            to have points from SI in the boundary of the convex hull */
/*            of SD. However the final tetrahedralization for SO, even */
/*            though it is Regular and equals or contains a Regular */
/*            tetrahedralization for SD, will not necessarily contain a */
/*            Regular tetrahedralization for SI. */

/*     redchk Input logical flag; */
/*            used only if logical variable delaun (described above) is */
/*            set to .false. (weights are being used); */
/*            if .true. then redundant points are to be tested for */
/*            redundancy after final tetrahedralization has been */
/*            computed; */
/*            if .false. then they are not to be tested for redundancy. */

/*     nmax   Input integer variable; */
/*            must be defined in a parameter statement in main routine; */
/*            dimension of single arrays x, y, z, w, v, ix, iy, iz, iw, */
/*            ix2, iy2, iz2, iw2, is, io (all are described below). */

/*     nvmax  Input integer variable; */
/*            must be defined in a parameter statement in main routine; */
/*            second dimension of double array icon and dimension of */
/*            single arrays ifl, id (these arrays are described below); */
/*            nvmax should be set to 55*nmax if logical variable flphis */
/*            (described above) is set to .true.; to 7*nmax otherwise. */

/*     nhmax  Input integer variable; */
/*            must be defined in a parameter statement in main routine; */
/*            dimension of single array ih (described below). */

/*     nv     Input integer variable that can not exceed nmax; */
/*            number of points or vertices in SI; */
/*            same value on output. */

/*     nw     Output integer variable that can not exceed nmax; */
/*            number of points or vertices in SO. */

/*     nt     Output integer variable that can not exceed nvmax; */
/*            number of tetrahedra in final tetrahedron list; */
/*            if logical variable artfcl is .false. then nt equals */
/*            the number of tetrahedra in the final tetrahedralization */
/*            for SO (equal to SD); */
/*            if artfcl is .true. and logical variable flphis is .true. */
/*            then nt equals the number of tetrahedra in the  final */
/*            tetrahedralization for SO (equal to SU) plus the */
/*            number of flipping history tetrahedra, i. e. tetrahedra */
/*            that at some point during the execution of the program were */
/*            part of a tetrahedralization for SO but that are not part */
/*            of the final tetrahedralization (a tetrahedron in a */
/*            previous tetrahedralization will not be in the final */
/*            tetrahedralization if at some time during the execution of */
/*            the program it has been eliminated and replaced by other */
/*            tetrahedra through a flip); */
/*            if artfcl is .true. and flphis is .false. then nt equals */
/*            the number of tetrahedra in the final tetrahedralization */
/*            for SO (equal to SU). */

/*     nd     Output integer variable that can not exceed nvmax; */
/*            number of real tetrahedra in final tetrahedron list; */
/*            i. e. number of tetrahedra in Regular tetrahedralization */
/*            for SD (the real points) which is contained in final */
/*            Regular tetrahedralization computed by the program for SU. */

/*     x      Input/output real*8 single array of dimension nmax; */
/*            on input for each i, i = 1, ..., nv, x(i) is the */
/*            x-coordinate of the ith point in SI; */
/*            with icfig as decribed below, on output for each i, */
/*            i = 1, ..., nw, if the ith point in SO is not in SA */
/*            then x(i) is the x-coordinate of the ith point in SO */
/*            rounded off so that its decimal part has icfig */
/*            significant figures; if it is in SA then x(i) is a */
/*            program generated value associated internally by the */
/*            program with the ith point in SO. */

/*     y      Input/output real*8 single array of dimension nmax; */
/*            on input for each i, i = 1, ..., nv, y(i) is the */
/*            y-coordinate of the ith point in SI; */
/*            with icfig as decribed below, on output for each i, */
/*            i = 1, ..., nw, if the ith point in SO is not in SA */
/*            then y(i) is the y-coordinate of the ith point in SO */
/*            rounded off so that its decimal part has icfig */
/*            significant figures; if it is in SA then y(i) is a */
/*            program generated value associated internally by the */
/*            program with the ith point in SO. */

/*     z      Input/output real*8 single array of dimension nmax; */
/*            on input for each i, i = 1, ..., nv, z(i) is the */
/*            z-coordinate of the ith point in SI; */
/*            with icfig as decribed below, on output for each i, */
/*            i = 1, ..., nw, if the ith point in SO is not in SA */
/*            then z(i) is the z-coordinate of the ith point in SO */
/*            rounded off so that its decimal part has icfig */
/*            significant figures; if it is in SA then z(i) is a */
/*            program generated value associated internally by the */
/*            program with the ith point in SO. */

/*     w      Input/output real*8 single array of dimension nmax; */
/*            on input for each i, i = 1, ..., nv, w(i) is the */
/*            weight of the ith point in SI; */
/*            with iwfig as decribed below, on output for each i, */
/*            i = 1, ..., nw, if the ith point in SO is not in SA */
/*            then w(i) is the weight of the ith point in SO */
/*            rounded off so that its decimal part has iwfig */
/*            significant figures; if it is in SA then w(i) is a */
/*            program generated value associated internally by the */
/*            program with the ith point in SO. */

/*     v      Real single array of dimension nmax; */
/*            internally used by program. */

/*     ix     Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     iy     Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     iz     Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     iw     Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     ix2    Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     iy2    Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     iz2    Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     iw2    Integer single array of dimension nmax; */
/*            internally used by the program. */

/*     icon   Output integer double array of dimensions 8 and nvmax; */
/*            this is the tetrahedron list; */
/*            actually this is a list of 8 x nt integers, and it is a */
/*            list of tetrahedra in the sense that for each j, */
/*            j = 1, ..., nt, the 8 integers icon(i,j), i = 1, .., 8, are */
/*            associated with a tetrahedron, the jth tetrahedron or */
/*            tetrahedron j, as will be described below; */
/*            if logical variables artfcl and flphis are both .true. */
/*            then tetrahedra in this list are those in the final */
/*            tetrahedralization for SO (equal to SU) together with the */
/*            flipping history tetrahedra, i. e. tetrahedra that at some */
/*            point during the execution of the program were part of a */
/*            tetrahedralization for SO but that are not part  of the */
/*            final tetrahedralization; */
/*            if artfcl is .true. and flphis is .false. then tetrahedra */
/*            in this list are just those in the final tetrahedralization */
/*            for SO (equal to SU); */
/*            if artfcl is .false. then tetrahedra in this list are those */
/*            in the final tetrahedralization for SO (equal to SD); */
/*            for each j, j = 1, ..., nt, if icon(5,j) is negative (can */
/*            only happen if artfcl equals .true.) then tetrahedron j */
/*            is not in the final tetrahedralization for SO, it was in a */
/*            previous tetrahedralization (it was eliminated) and its */
/*            vertices are the points -icon(5,j), icon(6,j), icon(7,j), */
/*            and icon(8,j) in SO; in addition if flphis is .true. the */
/*            tetrahedra by which tetrahedron j was replaced through a */
/*            flip can be identified  as follows: for each i, */
/*            i = 1, ..., 4, if icon(i,j) is positive then tetrahedron */
/*            icon(i,j) is one of those tetrahedra; */
/*            for each j, j = 1, ..., nt, if icon(5,j) is positive then */
/*            tetrahedron j is in the final tetrahedralization for SO, */
/*            and its vertices are the points icon(5,j), icon(6,j), */
/*            icon(7,j), and icon(8,j) in SO; in addition the tetrahedra */
/*            in the final tetrahedralization that share a facet with */
/*            tetrahedron j can be identified as follows: for each i, */
/*            i = 1, ..., 4, if icon(i,j) is positive then tetrahedron */
/*            icon(i,j) is one of those tetrahedra; */
/*            for each j, j = 1, ..., nt, the vertices of tetrahedron j */
/*            are ordered as follows: when viewed from vertex icon(5,j) */
/*            (-icon(5,j) if icon(5,j) is negative) the other three */
/*            vertices icon(6,j), icon(7,j), icon(8,j) appear in this */
/*            order in a clockwise direction around the circle that */
/*            contains them; */
/*            for each j, j = 1, ..., nt, if tetrahedron j is in the */
/*            final tetrahedralization, i. e. icon(5,j) is positive, */
/*            then the tetrahedra in the final tetrahedralization that */
/*            share a facet with tetrahedron j are ordered as follows: */
/*            for each i, i = 1, ..., 4, if icon(i,j) is positive, */
/*            tetrahedron j shares with tetrahedron icon(i,j) the */
/*            facet of tetrahedron j that does not contain vertex */
/*            icon(i+4,j). */

/*     is     Input/output integer single array of dimension nmax; */
/*            on input if logical variable pntoff is .true. then for */
/*            each i, i = 1, ..., nv, if the value of is(i) equals 1 */
/*            then the ith point in SI is to be active during the */
/*            tetrahedralization computation; if it equals 0 then the */
/*            ith point in SI is to be inactive during the computation; */
/*            on input if logical variable pntoff is .false. then */
/*            all points in SI are to be active during the */
/*            tetrahedralization computation and array is does not */
/*            have to be set to any values; */
/*            on output for each i, i = 1, ..., nw, the value of is(i) */
/*            indicates how the ith point in SO was processed by the */
/*            program as follows: */
/*            if is(i) is zero then point i was not considered as a */
/*            vertex for tetrahedralization; */
/*            if is(i) is positive then point i is part of the final */
/*            tetrahedralization for SO, i. e. there is at least one */
/*            tetrahedron in the final tetrahedralization with point i */
/*            as a vertex, and tetrahedron is(i) is one such tetrahedron */
/*            (actually if point i is in SA then is(i) is always */
/*            positive); */
/*            if is(i) is less than -8 then point i was found to be */
/*            redundant as the program was trying to insert it into the */
/*            current tetrahedralization because a point previously */
/*            inserted (point -is(i) in SO if artfcl is .true. (SO equals */
/*            SU), point -is(i)-8 in SO if artfcl is .false. (SO equals */
/*            SD)) was identical to it and either the weight of the */
/*            previously inserted point was larger or equal to the weight */
/*            of point i or there were no weights; */
/*            if is(i) equals -2 then point i had been inserted by the */
/*            program into the tetrahedralization but was found to be */
/*            redundant because another point was later inserted by the */
/*            program that was identical to point i and whose weight was */
/*            larger than that of point i (this case is not possible if */
/*            there are no weights); */
/*            if is(i) equals -3 then point i was found to be redundant */
/*            in the sense of a Regular tetrahedralization as the program */
/*            was trying to insert it into the current tetrahedralization */
/*            because of its weight as compared to the weights of the */
/*            vertices of the tetrahedron in the current */
/*            tetrahedralization that contains it even though it was not */
/*            identical to a previously inserted point (this case is not */
/*            possible if there are no weights); */
/*            if is(i) equals -4 then point i had been inserted by the */
/*            program into the tetrahedralization but was found to be */
/*            redundant in the sense of a Regular tetrahedralization */
/*            because of the weight of another point, not identical to */
/*            point i, that was later inserted by the program together */
/*            with the weights of three other previously inserted points */
/*            as compared to the weight of point i (this case is not */
/*            possible if there are no weights). */

/*     ifl    Integer single array of dimension nvmax; */
/*            internally used by program. */

/*     io     Integer single array of dimension nmax; */
/*            internally used by program. */

/*     id     Integer single array of dimension nvmax; */
/*            internally used by program. */

/*     ih     Integer single array of dimension nhmax; */
/*            internally used by program. */

/*     wlenx */
/*     wleny */
/*     wlenz  Input real*8 variables; */
/*            If reccor is .true. then these are three positive real */
/*            numbers provided by the user to be used by the program */
/*            to identify a rectangular polyhedron that contains SI in */
/*            its interior. This is the polyhedron whose surface will */
/*            contain the set SP. If xmax, ymax, zmax are the maximum */
/*            values of the x, y, and z coordinates of the points in SI, */
/*            respectively, and xmin, ymin, zmin are the minimum values, */
/*            then the eight vertices of the polyhedron will be: */
/*            (xmin-wlenx, ymin-wleny, zmin-wlenz), */
/*            (xmax+wlenx, ymin-wleny, zmin-wlenz), */
/*            (xmax+wlenx, ymax+wleny  zmin-wlenz), */
/*            (xmin-wlenx, ymax+wleny, zmin-wlenz), */
/*            (xmin-wlenx, ymin-wleny, zmax+wlenz), */
/*            (xmax+wlenx, ymin-wleny, zmax+wlenz), */
/*            (xmax+wlenx, ymax+wleny  zmax+wlenz), */
/*            (xmin-wlenx, ymax+wleny, zmax+wlenz). */

/*     wlenw  Input real*8 variable; */
/*            If reccor is .true. and delaun is .false. then this is a */
/*            real number provided  by the user to be used by the program */
/*            to determine a weight to be assigned to each point in SP. */
/*            If wmin is the minimum value of the weights of the points */
/*            in SI then this weight will be wmin - wlenw. */

/*     naddl  Input integer variable; */
/*            If reccor is .true. then this is a positive integer */
/*            greater than 1 provided by the user to be used by the */
/*            program to determine the set SP. The points in SP define */
/*            a rectangular regular grid on the surface of a rectangular */
/*            polyhedron (described above) that contains SI in its */
/*            interior. For each facet of the polyhedron a set of naddl */
/*            x naddl points is generated by the program that defines a */
/*            rectangular regular grid on the facet and that contains */
/*            the four vertices of the facet. SP is then the union of */
/*            the six sets thus generated (one per facet). It then */
/*            follows that the number of points in SP must be */
/*            6(naddl-2)(naddl-2)+12(naddl-2)+8 which reduces to */
/*            6(naddl)(naddl)-12(naddl)+8. */
/*            It also follows that if naddl equals 2 then the points in */
/*            SP are exactly the 8 vertices of the polyhedron. */

/*     isu */
/*     jsu */
/*     ksu */
/*     nsu    Input integer variables; */
/*            If random is .true. then these are four integers provided */
/*            by the user to be used by program as seeds for identifying */
/*            random order in which points in SI are to be added; */
/*            they can be any four integers. */

/*     icfig  Input integer variable; */
/*            value is the number of significant figures of decimal part */
/*            of coordinates of input points; value should be nonnegative */
/*            and less than 10. */

/*     iwfig  Input integer variable; */
/*            value is the number of significant figures of decimal part */
/*            of weights (if any) of input points; value should be */
/*            nonnegative, less than 10, and not greater than twice the */
/*            value of variable icfig (described above). */

/*     epz    Input real*8 variable; */
/*            tolerance used by the program to switch from floating */
/*            point arithmetic to exact arithmetic by testing against */
/*            this tolerance whether certain quantities are too close */
/*            to zero; setting it equal to numbers such as 0.1, 0.01 */
/*            has worked well so far. */

/* ********************************************************************** */

/*     Examples of settings for logical variables delaun, pntoff, */
/*     flphis, artfcl, random, reccor, and redchk: */

/*     Delaunay tetrahedralization for set of randomized input points */
/*     is desired and nothing else (number of input points equals number */
/*     of output points, all input points are to be active during the */
/*     tetrahedralization computation, and tetrahedron list is exactly */
/*     a list of the tetrahedra in the final tetrahedralization for the */
/*     set of input points): */

/*     delaun = .true. */
/*     pntoff = .false. */
/*     flphis = .true. */
/*     artfcl = .false. */
/*     random = .false. */
/*     reccor = .false. */

/*     The same as above but a Regular tetrahedralization is desired, */
/*     the input points are not randomized, the flipping history */
/*     is to be used for locating points, and redundant points are not */
/*     to be tested for redundancy after the final tetrahedralization */
/*     has been computed: */

/*     delaun = .false. */
/*     pntoff = .false. */
/*     flphis = .true. */
/*     artfcl = .false. */
/*     random = .true. */
/*     reccor = .false. */
/*     redchk = .false. */

/*     The same as above but a Regular tetrahedralization is desired, */
/*     the input points are not randomized, the output tetrahedron */
/*     list is to include artificial tetrahedra information, and the */
/*     shishkebab method is to be used for locating points: */

/*     delaun = .false. */
/*     pntoff = .false. */
/*     flphis = .false. */
/*     artfcl = .true. */
/*     random = .true. */
/*     reccor = .false. */
/*     redchk = .false. */

/* ********************************************************************** */

/* Subroutine */ int regtet_(doublereal *x, doublereal *y, doublereal *z__, 
	doublereal *w, real *v, integer *ix, integer *iy, integer *iz, 
	integer *iw, integer *ix2, integer *iy2, integer *iz2, integer *iw2, 
	integer *icon, integer *is, integer *ifl, integer *io, integer *id, 
	integer *ih, integer *nv, integer *nw, integer *nt, integer *nd, 
	integer *nmax, integer *nvmax, integer *nhmax, doublereal *wlenx, 
	doublereal *wleny, doublereal *wlenz, doublereal *wlenw, integer *
	naddl, integer *isu, integer *jsu, integer *ksu, integer *nsu, 
	integer *icfig, integer *iwfig, doublereal *epz, logical *delaun, 
	logical *pntoff, logical *flphis, logical *artfcl, logical *random, 
	logical *reccor, logical *redchk)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, ni, no, nu, nv1, ihn, irec, ipre, jpre, npre;
    doublereal wmin, xmin, ymin, xmax, ymax, zmin, zmax, wmax, xpre, ypre, 
	    zpre, wpre;
    integer inow, jnow;
    doublereal wnow, xnow, ynow, znow;
    integer irec1, ired0, ired1, ired2, ired3, ired4, mhalf, iftal, iredp, 
	    isclp[2], isclr[2], iredx, tetra, mfull, isclw[2], iconx, iorix, 
	    isphx, tetru;
    extern /* Subroutine */ int delchk_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, logical *, integer *, integer *, integer *,
	     integer *, integer *, doublereal *), rdmord_(real *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), consis_(
	    integer *, integer *, integer *, integer *, integer *), orient_(
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *), convex_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *);
    integer mxlook;
    extern /* Subroutine */ int poltri_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, logical *, logical *, logical *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *), revtet_(integer *, integer *, integer *, integer *,
	     integer *, integer *, logical *), ruvtet_(integer *, integer *, 
	    integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___48 = { 0, 6, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };




/*     integer nva */

/*     initialize Fortran 77 word lengths */

    /* Parameter adjustments */
    --ih;
    --id;
    --io;
    --ifl;
    --is;
    icon -= 9;
    --iw2;
    --iz2;
    --iy2;
    --ix2;
    --iw;
    --iz;
    --iy;
    --ix;
    --v;
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    mhalf = 32768;
    mfull = 1073741824;

/*     testing parameters and number of input points or vertices */

    if (*nv < 1 || *nv > *nmax || *nvmax < 12) {
	s_wsle(&io___3);
	do_lio(&c__3, &c__1, (char *)&(*nv), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&(*nvmax), (ftnlen)sizeof(integer));
	e_wsle();
	s_stop("110", (ftnlen)3);
    }

/*     initialize arrays ih, w, is and id */

    i__1 = *nhmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ih[i__] = 0;
/* L50: */
    }
    if (*delaun) {
	i__1 = *nmax;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    w[i__] = 0.;
/* L60: */
	}
    }
    if (! (*pntoff)) {
	i__1 = *nmax;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    is[i__] = 1;
/* L80: */
	}
    }
    if (! (*flphis)) {
	ihn = 0;
	iftal = 0;
	i__1 = *nvmax;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    id[i__] = 0;
/* L100: */
	}
    }

/*     test variables associated with a possible rectangular polyhedron */

    if (*reccor) {
	if (*wlenx <= 0. || *wleny <= 0. || *wlenz <= 0.) {
	    s_stop("120", (ftnlen)3);
	}
	if (*naddl < 2) {
	    s_stop("130", (ftnlen)3);
	}
    } else {
	*wlenx = 0.;
	*wleny = 0.;
	*wlenz = 0.;
	*wlenw = 0.;
	*naddl = 0;
    }

/*     calculating min and max */

    xmax = x[1];
    xmin = x[1];
    ymax = y[1];
    ymin = y[1];
    zmax = z__[1];
    zmin = z__[1];
    wmax = w[1];
    wmin = w[1];
    i__1 = *nv;
    for (no = 1; no <= i__1; ++no) {
	if (x[no] > xmax) {
	    xmax = x[no];
	}
	if (x[no] < xmin) {
	    xmin = x[no];
	}
	if (y[no] > ymax) {
	    ymax = y[no];
	}
	if (y[no] < ymin) {
	    ymin = y[no];
	}
	if (z__[no] > zmax) {
	    zmax = z__[no];
	}
	if (z__[no] < zmin) {
	    zmin = z__[no];
	}
	if (w[no] > wmax) {
	    wmax = w[no];
	}
	if (w[no] < wmin) {
	    wmin = w[no];
	}
/* L150: */
    }

/*     if randomizing input data then obtain random order of integers */
/*     from 1 to nv and randomize data */

    if (! (*random)) {
	goto L185;
    }
    rdmord_(&v[1], &io[1], nv, isu, jsu, ksu, nsu);
    i__1 = *nv;
    for (no = 1; no <= i__1; ++no) {
	if (io[no] < 0) {
	    goto L180;
	}
	ipre = no;
	xpre = x[ipre];
	ypre = y[ipre];
	zpre = z__[ipre];
	wpre = w[ipre];
	jpre = is[ipre];
L170:
	inow = io[ipre];
	io[ipre] = -inow;
	xnow = x[inow];
	ynow = y[inow];
	znow = z__[inow];
	wnow = w[inow];
	jnow = is[inow];
	x[inow] = xpre;
	y[inow] = ypre;
	z__[inow] = zpre;
	w[inow] = wpre;
	is[inow] = jpre;
	if (inow == no) {
	    goto L180;
	}
	ipre = inow;
	xpre = xnow;
	ypre = ynow;
	zpre = znow;
	wpre = wnow;
	jpre = jnow;
	if (io[ipre] < 0) {
	    s_stop("140", (ftnlen)3);
	}
	goto L170;
L180:
	;
    }

/*     OPEN(22,FILE='ran.dat') */
/*     DO 182 I=1,NV */
/*        WRITE(22,*)X(I),Y(I),Z(I),W(I) */
/* 182 CONTINUE */

/*     shift data */

L185:
    irec = 8;
    if (*reccor) {
/* Computing 2nd power */
	i__1 = *naddl;
	irec = irec + i__1 * i__1 * 6 - *naddl * 12 + 8;
    }
    irec1 = irec + 1;
    nv1 = *nv;
    *nv += irec;
    if (*nv > *nmax) {
	s_stop("150", (ftnlen)3);
    }
    i__1 = irec1;
    for (no = *nv; no >= i__1; --no) {
	nu = no - irec;
	x[no] = x[nu];
	y[no] = y[nu];
	z__[no] = z__[nu];
	w[no] = w[nu];
	is[no] = is[nu];
	if (*random) {
	    if (io[nu] >= 0) {
		s_stop("160", (ftnlen)3);
	    }
	    io[nu] = -io[nu] + irec;
	}
/* L190: */
    }

/*     initialize is for additional data */

    i__1 = irec;
    for (i__ = 1; i__ <= i__1; ++i__) {
	is[i__] = 1;
/* L200: */
    }

/*     write(*,*)' ' */
/*     write(*,*)'Entering poltri ...' */
    poltri_(&x[1], &y[1], &z__[1], &w[1], &ix[1], &iy[1], &iz[1], &iw[1], &
	    ix2[1], &iy2[1], &iz2[1], &iw2[1], &icon[9], &is[1], &ifl[1], &id[
	    1], &ih[1], &ihn, &xmin, &ymin, &zmin, &wmin, &xmax, &ymax, &zmax,
	     &wmax, &iftal, nv, nmax, nvmax, nhmax, wlenx, wleny, wlenz, 
	    wlenw, &tetra, &mxlook, &irec, naddl, &iredx, delaun, flphis, 
	    redchk, icfig, iwfig, &mhalf, &mfull, isclp, isclw, isclr, epz);
/*     write(*,*)' ' */
/*     write(*,*)'Checking tetrahedralization ...' */
/*     WRITE(*,*)' ' */
/*     WRITE(*,*)'MAXLOOK=',MAXLOOK,' IHN=',IHN */
/*     write (*,*)' ' */
/*     write (*,*)'Leaving poltri ...' */

/*     write(*,*)' ' */
/*     write(*,*)'Entering consis ...' */
    consis_(&icon[9], &is[1], &ifl[1], nv, &tetra);

/*     write(*,*)' ' */
/*     write(*,*)'Entering convex ...' */
    convex_(&icon[9], &tetra, &ifl[1], &x[1], &y[1], &z__[1], &ix[1], &iy[1], 
	    &iz[1], &ix2[1], &iy2[1], &iz2[1], &iconx, &mhalf, &mfull, isclp, 
	    epz);

/*     write(*,*)' ' */
/*     write(*,*)'Entering orient ...' */
    orient_(&tetra, &icon[9], &ifl[1], &x[1], &y[1], &z__[1], &ix[1], &iy[1], 
	    &iz[1], &ix2[1], &iy2[1], &iz2[1], &iorix, &mhalf, &mfull, isclp, 
	    epz);

/*     write(*,*)' ' */
/*     write(*,*)'Entering delchk ...' */
    delchk_(&tetra, &icon[9], &ifl[1], &x[1], &y[1], &z__[1], &w[1], &ix[1], &
	    iy[1], &iz[1], &iw[1], &ix2[1], &iy2[1], &iz2[1], &iw2[1], &isphx,
	     delaun, &mhalf, &mfull, isclp, isclw, isclr, epz);

/*     checking for possible warnings */

    if (iredx != 0) {
	s_wsle(&io___41);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___42);
	do_lio(&c__9, &c__1, "Warning: redundancy violations detected", (
		ftnlen)39);
	e_wsle();
	s_wsle(&io___43);
	do_lio(&c__9, &c__1, "Number of violations = ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iredx, (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (iconx != 0) {
	s_wsle(&io___44);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___45);
	do_lio(&c__9, &c__1, "Warning: convexity violations detected", (
		ftnlen)38);
	e_wsle();
	s_wsle(&io___46);
	do_lio(&c__9, &c__1, "Number of violations = ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iconx, (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (iorix != 0) {
	s_wsle(&io___47);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___48);
	do_lio(&c__9, &c__1, "Warning: tetrahedra orientation violations det\
ected", (ftnlen)51);
	e_wsle();
	s_wsle(&io___49);
	do_lio(&c__9, &c__1, "Number of violations = ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&iorix, (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (isphx != 0) {
	s_wsle(&io___50);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___51);
	do_lio(&c__9, &c__1, "Warning: sphere criterion violations detected", 
		(ftnlen)45);
	e_wsle();
	s_wsle(&io___52);
	do_lio(&c__9, &c__1, "Number of violations = ", (ftnlen)23);
	do_lio(&c__3, &c__1, (char *)&isphx, (ftnlen)sizeof(integer));
	e_wsle();
    }
    if (iredx != 0 || iconx != 0 || iorix != 0 || isphx != 0) {
	s_wsle(&io___53);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___54);
	do_lio(&c__9, &c__1, "Increasing tolerance epz could improve situati\
on", (ftnlen)48);
	e_wsle();
    }

/*     readjust data structure for randomizing */

    if (! (*random)) {
	goto L290;
    }
/*     write(*,*)' ' */
/*     write(*,*)'Readjusting data structure for randomizing ...' */
    if (*nv > *nvmax) {
	s_stop("170", (ftnlen)3);
    }
    nu = *nv - irec;
    i__1 = nu;
    for (no = 1; no <= i__1; ++no) {
	ifl[io[no]] = no + irec;
/* L230: */
    }
    i__1 = tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 5; j <= 8; ++j) {
	    if (icon[j + (i__ << 3)] > irec) {
		icon[j + (i__ << 3)] = ifl[icon[j + (i__ << 3)]];
	    } else if (icon[j + (i__ << 3)] < -irec) {
		icon[j + (i__ << 3)] = -ifl[-icon[j + (i__ << 3)]];
	    }
/* L240: */
	}
/* L250: */
    }
    i__1 = *nv;
    for (i__ = irec1; i__ <= i__1; ++i__) {
	if (is[i__] < -8) {
	    is[i__] = -ifl[-is[i__]];
	}
/* L255: */
    }
    i__1 = *nv;
    for (i__ = irec1; i__ <= i__1; ++i__) {
	ifl[i__] = is[i__];
/* L260: */
    }
    i__1 = *nv;
    for (i__ = irec1; i__ <= i__1; ++i__) {
	is[i__] = ifl[io[i__ - irec]];
/* L270: */
    }
    i__1 = nv1;
    for (no = 1; no <= i__1; ++no) {
	if (io[no] < 0) {
	    goto L280;
	}
	nu = no + irec;
	ipre = nu;
	npre = no;
	xpre = x[ipre];
	ypre = y[ipre];
	zpre = z__[ipre];
	wpre = w[ipre];
L275:
	inow = io[npre];
	io[npre] = -inow;
	if (inow == nu) {
	    x[ipre] = xpre;
	    y[ipre] = ypre;
	    z__[ipre] = zpre;
	    w[ipre] = wpre;
	    goto L280;
	}
	x[ipre] = x[inow];
	y[ipre] = y[inow];
	z__[ipre] = z__[inow];
	w[ipre] = w[inow];
	ipre = inow;
	npre = ipre - irec;
	if (io[npre] < 0) {
	    s_stop("180", (ftnlen)3);
	}
	goto L275;
L280:
	;
    }

/*     write(*,*)' ' */
/*     write(*,*)'Entering consis ...' */
    consis_(&icon[9], &is[1], &ifl[1], nv, &tetra);

L290:
    nu = *nv - 8;
    if (! (*artfcl)) {
/*        write(*,*)' ' */
/*        write(*,*)'Entering revtet ...' */
	revtet_(&tetra, &tetru, &icon[9], nv, &is[1], &ifl[1], flphis);
	i__1 = nu;
	for (no = 1; no <= i__1; ++no) {
	    x[no] = x[no + 8];
	    y[no] = y[no + 8];
	    z__[no] = z__[no + 8];
	    w[no] = w[no + 8];
/* L293: */
	}
	if (tetru == 0) {
	    goto L300;
	}

/*        write(*,*)' ' */
/*        write(*,*)'Entering consis ...' */
	consis_(&icon[9], &is[1], &ifl[1], &nu, &tetru);
    } else if (! (*flphis)) {
	ruvtet_(&tetra, &tetru, &icon[9], &is[1], &ifl[1]);
	consis_(&icon[9], &is[1], &ifl[1], nv, &tetra);
    } else {

/*        count true tetrahedra */

	tetru = 0;
	i__1 = tetra;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (icon[(i__ << 3) + 6] <= 8 || icon[(i__ << 3) + 7] <= 8 || 
		    icon[(i__ << 3) + 8] <= 8 || icon[(i__ << 3) + 5] <= 8) {
		goto L295;
	    }
	    ++tetru;
L295:
	    ;
	}
    }

/*     count redundant vertices */

L300:
    *nd = tetru;
    iredp = 0;
    ired0 = 0;
    ired1 = 0;
    ired2 = 0;
    ired3 = 0;
    ired4 = 0;
    if (*artfcl) {
	ni = 9;
	*nw = *nv;
	*nt = tetra;
    } else {
	ni = 1;
	*nw = nu;
	*nt = tetru;
    }
    i__1 = *nw;
    for (i__ = ni; i__ <= i__1; ++i__) {
	if (is[i__] > 0) {
	    ++iredp;
	} else if (is[i__] == 0) {
	    ++ired0;
	} else if (is[i__] < -8) {
	    ++ired1;
	} else if (is[i__] == -2) {
	    ++ired2;
	} else if (is[i__] == -3) {
	    ++ired3;
	} else if (is[i__] == -4) {
	    ++ired4;
	} else {
	    s_stop("190", (ftnlen)3);
	}
/* L400: */
    }

/*     OPEN(23,FILE='unr.dat') */
/*     DO 500 I=1,NW */
/*        WRITE(23,*)X(I),Y(I),Z(I),W(I) */
/* 500 CONTINUE */

/*     nva = nv */
    *nv = nv1;

/*     write info to screen */

/*     wtenv=float(tetra)/float(nva) */
/*     wtena=float(tetra)/float(iredp) */
/*     wtuna=float(tetru)/float(iredp) */
/*     write (*,*) ' ' */
/*     write (*,*) 'Tetrahedralization data ...' */
/*     write (*,*) ' ' */
/*     write (*,*) 'minimum weight = ',wmin */
/*     write (*,*) 'maximum weight = ',wmax */
/*     write (*,*) 'number of true vertices: ', nu */
/*     write (*,*) 'number of active vertices: ',iredp */
/*     write (*,*) 'maximum number of vertices parameter = ', nmax */
/*     write (*,*) 'maximum number of tetrahed parameter = ', nvmax */
/*     write (*,*) 'number of tetrahedra of all kinds: ', tetra */
/*     write (*,*) 'all tetrahedra-all vertices ratio: ',wtenv */
/*     write (*,*) 'number of true tetrahedra: ', tetru */
/*     write (*,*) ' all tetrahedra-active vertices ratio: ',wtena */
/*     write (*,*) 'true tetrahedra-active vertices ratio: ',wtuna */
/*     write (*,*) 'max levels gone down in hierarchy = ', mxlook */
/*     write (*,*) 'points active at the end of current run   = ',iredp */
/*     write (*,*) 'points inactive at the end of current run = ',ired0 */
/*     write (*,*) 'points redundant by initial substitution  = ',ired1 */
/*     write (*,*) 'points redundant by eventual substitution = ',ired2 */
/*     write (*,*) 'points redundant by initial elimination   = ',ired3 */
/*     write (*,*) 'points redundant by eventual elimination  = ',ired4 */

    return 0;
} /* regtet_ */

/* POLTRI */

/*     This subroutine will obtain initial cube and will divide it */
/*     into 12 tetrahedra; insert points into tetrahedralization */

/* Subroutine */ int poltri_(doublereal *x, doublereal *y, doublereal *z__, 
	doublereal *w, integer *ix, integer *iy, integer *iz, integer *iw, 
	integer *ix2, integer *iy2, integer *iz2, integer *iw2, integer *icon,
	 integer *is, integer *ifl, integer *id, integer *ih, integer *ihn, 
	doublereal *xmin, doublereal *ymin, doublereal *zmin, doublereal *
	wmin, doublereal *xmax, doublereal *ymax, doublereal *zmax, 
	doublereal *wmax, integer *iftal, integer *nv, integer *nmax, integer 
	*nvmax, integer *nhmax, doublereal *wlenx, doublereal *wleny, 
	doublereal *wlenz, doublereal *wlenw, integer *tetra, integer *mxlook,
	 integer *irec, integer *naddl, integer *iredx, logical *delaun, 
	logical *flphis, logical *redchk, integer *icsfig, integer *iwsfig, 
	integer *mhalf, integer *mfull, integer *isclp, integer *isclw, 
	integer *isclr, doublereal *epz)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer i_dnnt(doublereal *), s_wsle(cilist *), do_lio(integer *, integer 
	    *, char *, ftnlen), e_wsle();
    double d_int(doublereal *);

    /* Local variables */
    integer i__, j, k, i9, i10, ng;
    doublereal xc[8], yc[8], zc[8];
    integer ixc[8], iyc[8], izc[8];
    doublereal xmn, ymn, zmn, xmx, ymx, zmx, wbig, wlan, wman, xlan, wlen, 
	    wmen, xlen, ylen, zlen, xctr, yctr, zctr, ylan, zlan, xint, yint, 
	    zint, xcor, ycor, zcor;
    integer itcn1, itcn2, itcn3, itcn4, naddm;
    doublereal decml, dscle, dfill, dscli;
    integer isgcl, itchk;
    doublereal dfull;
    integer iredu, isclu, issin, icsfi2, ibsfig;
    extern /* Subroutine */ int decomp_(integer *, integer *, integer *, 
	    integer *);
    integer irsfig, itsfig;
    extern /* Subroutine */ int pntred_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, integer *, integer *, integer *), 
	    pntins_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *);

    /* Fortran I/O blocks */
    static cilist io___101 = { 0, 6, 0, 0, 0 };
    static cilist io___102 = { 0, 6, 0, 0, 0 };
    static cilist io___103 = { 0, 6, 0, 0, 0 };
    static cilist io___105 = { 0, 6, 0, 0, 0 };
    static cilist io___106 = { 0, 6, 0, 0, 0 };
    static cilist io___107 = { 0, 6, 0, 0, 0 };
    static cilist io___108 = { 0, 6, 0, 0, 0 };
    static cilist io___109 = { 0, 6, 0, 0, 0 };
    static cilist io___110 = { 0, 6, 0, 0, 0 };
    static cilist io___111 = { 0, 6, 0, 0, 0 };
    static cilist io___112 = { 0, 6, 0, 0, 0 };
    static cilist io___113 = { 0, 6, 0, 0, 0 };
    static cilist io___114 = { 0, 6, 0, 0, 0 };
    static cilist io___115 = { 0, 6, 0, 0, 0 };
    static cilist io___116 = { 0, 6, 0, 0, 0 };
    static cilist io___117 = { 0, 6, 0, 0, 0 };
    static cilist io___118 = { 0, 6, 0, 0, 0 };
    static cilist io___119 = { 0, 6, 0, 0, 0 };
    static cilist io___120 = { 0, 6, 0, 0, 0 };
    static cilist io___121 = { 0, 6, 0, 0, 0 };
    static cilist io___122 = { 0, 6, 0, 0, 0 };
    static cilist io___125 = { 0, 6, 0, 0, 0 };
    static cilist io___126 = { 0, 6, 0, 0, 0 };
    static cilist io___127 = { 0, 6, 0, 0, 0 };
    static cilist io___128 = { 0, 6, 0, 0, 0 };
    static cilist io___129 = { 0, 6, 0, 0, 0 };
    static cilist io___130 = { 0, 6, 0, 0, 0 };
    static cilist io___131 = { 0, 6, 0, 0, 0 };
    static cilist io___132 = { 0, 6, 0, 0, 0 };
    static cilist io___138 = { 0, 6, 0, 0, 0 };
    static cilist io___139 = { 0, 6, 0, 0, 0 };
    static cilist io___140 = { 0, 6, 0, 0, 0 };
    static cilist io___141 = { 0, 6, 0, 0, 0 };
    static cilist io___154 = { 0, 6, 0, 0, 0 };
    static cilist io___155 = { 0, 6, 0, 0, 0 };




/*     initialize */

    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --ih;
    --id;
    --ifl;
    --is;
    icon -= 9;
    --iw2;
    --iz2;
    --iy2;
    --ix2;
    --iw;
    --iz;
    --iy;
    --ix;
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    *mxlook = 0;

/*     test parameters */

    if (*nv > *nmax) {
	s_stop("210", (ftnlen)3);
    }
    if (*nv < 9) {
	s_stop("220", (ftnlen)3);
    }
    if (*nvmax < 12) {
	s_stop("230", (ftnlen)3);
    }

/*     construct cube */

    xlen = *xmax - *xmin;
    ylen = *ymax - *ymin;
    zlen = *zmax - *zmin;
/* Computing MAX */
    d__1 = max(xlen,ylen);
/* Computing MAX */
    d__2 = max(*wlenx,*wleny);
    wlan = max(d__1,zlen) / 2. + max(d__2,*wlenz);
    wlen = wlan + 4.;
    if (wlen <= wlan) {
	s_stop("235", (ftnlen)3);
    }

    xctr = (*xmax + *xmin) / 2.;
    yctr = (*ymax + *ymin) / 2.;
    zctr = (*zmax + *zmin) / 2.;
/*     WRITE(*,*)'XYCTR=',XCTR,YCTR,ZCTR,' WLEN=',WLEN */

/*     identify cube corner direction vectors */

    xc[0] = -1.;
    yc[0] = -1.;
    zc[0] = 1.;

    xc[1] = -1.;
    yc[1] = 1.;
    zc[1] = 1.;

    xc[2] = 1.;
    yc[2] = 1.;
    zc[2] = 1.;

    xc[3] = 1.;
    yc[3] = -1.;
    zc[3] = 1.;

    xc[4] = -1.;
    yc[4] = -1.;
    zc[4] = -1.;

    xc[5] = -1.;
    yc[5] = 1.;
    zc[5] = -1.;

    xc[6] = 1.;
    yc[6] = 1.;
    zc[6] = -1.;

    xc[7] = 1.;
    yc[7] = -1.;
    zc[7] = -1.;

/*     compute corners of cube */

    for (i__ = 1; i__ <= 8; ++i__) {
	x[i__] = xctr + wlen * xc[i__ - 1];
	y[i__] = yctr + wlen * yc[i__ - 1];
	z__[i__] = zctr + wlen * zc[i__ - 1];
	if (x[i__] >= *xmin && x[i__] <= *xmax || y[i__] >= *ymin && y[i__] <=
		 *ymax || z__[i__] >= *zmin && z__[i__] <= *zmax) {
	    s_stop("240", (ftnlen)3);
	}
/* L50: */
    }

/*     make coordinates of corners of cube into whole numbers */

    dfull = (doublereal) (*mfull);
    if ((d__1 = x[3] + 1., abs(d__1)) >= dfull) {
	s_stop("242", (ftnlen)3);
    }
    if ((d__1 = y[3] + 1., abs(d__1)) >= dfull) {
	s_stop("243", (ftnlen)3);
    }
    if ((d__1 = z__[3] + 1., abs(d__1)) >= dfull) {
	s_stop("244", (ftnlen)3);
    }
    if ((d__1 = x[5] - 1., abs(d__1)) >= dfull) {
	s_stop("245", (ftnlen)3);
    }
    if ((d__1 = y[5] - 1., abs(d__1)) >= dfull) {
	s_stop("246", (ftnlen)3);
    }
    if ((d__1 = z__[5] - 1., abs(d__1)) >= dfull) {
	s_stop("247", (ftnlen)3);
    }

    d__1 = x[3] + 1.;
    xmx = (doublereal) i_dnnt(&d__1);
    d__1 = y[3] + 1.;
    ymx = (doublereal) i_dnnt(&d__1);
    d__1 = z__[3] + 1.;
    zmx = (doublereal) i_dnnt(&d__1);
    d__1 = x[5] - 1.;
    xmn = (doublereal) i_dnnt(&d__1);
    d__1 = y[5] - 1.;
    ymn = (doublereal) i_dnnt(&d__1);
    d__1 = z__[5] - 1.;
    zmn = (doublereal) i_dnnt(&d__1);

    xlan = xmx - xmn;
    ylan = ymx - ymn;
    zlan = zmx - zmn;
/* Computing MAX */
    d__1 = max(xlan,ylan);
    wlan = max(d__1,zlan);

    x[1] = xmn;
    y[1] = ymn;
    z__[1] = zmn + wlan;

    x[2] = xmn;
    y[2] = ymn + wlan;
    z__[2] = zmn + wlan;

    x[3] = xmn + wlan;
    y[3] = ymn + wlan;
    z__[3] = zmn + wlan;

    x[4] = xmn + wlan;
    y[4] = ymn;
    z__[4] = zmn + wlan;

    x[5] = xmn;
    y[5] = ymn;
    z__[5] = zmn;

    x[6] = xmn;
    y[6] = ymn + wlan;
    z__[6] = zmn;

    x[7] = xmn + wlan;
    y[7] = ymn + wlan;
    z__[7] = zmn;

    x[8] = xmn + wlan;
    y[8] = ymn;
    z__[8] = zmn;

    for (i__ = 1; i__ <= 8; ++i__) {
	if (x[i__] >= *xmin && x[i__] <= *xmax || y[i__] >= *ymin && y[i__] <=
		 *ymax || z__[i__] >= *zmin && z__[i__] <= *zmax) {
	    s_stop("250", (ftnlen)3);
	}
/* L55: */
    }

    if (*irec == 8) {
	goto L155;
    }

/*     compute corners of rectangular polyhedron */

    x[9] = *xmin - *wlenx;
    y[9] = *ymin - *wleny;
    z__[9] = *zmax + *wlenz;

    x[10] = *xmin - *wlenx;
    y[10] = *ymax + *wleny;
    z__[10] = *zmax + *wlenz;

    x[11] = *xmax + *wlenx;
    y[11] = *ymax + *wleny;
    z__[11] = *zmax + *wlenz;

    x[12] = *xmax + *wlenx;
    y[12] = *ymin - *wleny;
    z__[12] = *zmax + *wlenz;

    x[13] = *xmin - *wlenx;
    y[13] = *ymin - *wleny;
    z__[13] = *zmin - *wlenz;

    x[14] = *xmin - *wlenx;
    y[14] = *ymax + *wleny;
    z__[14] = *zmin - *wlenz;

    x[15] = *xmax + *wlenx;
    y[15] = *ymax + *wleny;
    z__[15] = *zmin - *wlenz;

    x[16] = *xmax + *wlenx;
    y[16] = *ymin - *wleny;
    z__[16] = *zmin - *wlenz;

    for (i__ = 9; i__ <= 16; ++i__) {
	if (x[i__] >= *xmin && x[i__] <= *xmax || y[i__] >= *ymin && y[i__] <=
		 *ymax || z__[i__] >= *zmin && z__[i__] <= *zmax) {
	    s_stop("260", (ftnlen)3);
	}
/* L60: */
    }
    if (x[1] >= x[9] || y[1] >= y[9] || z__[1] <= z__[9] || x[2] >= x[10] || 
	    y[2] <= y[10] || z__[2] <= z__[10] || x[3] <= x[11] || y[3] <= y[
	    11] || z__[3] <= z__[11] || x[4] <= x[12] || y[4] >= y[12] || z__[
	    4] <= z__[12] || x[5] >= x[13] || y[5] >= y[13] || z__[5] >= z__[
	    13] || x[6] >= x[14] || y[6] <= y[14] || z__[6] >= z__[14] || x[7]
	     <= x[15] || y[7] <= y[15] || z__[7] >= z__[15] || x[8] <= x[16] 
	    || y[8] >= y[16] || z__[8] >= z__[16]) {
	s_stop("270", (ftnlen)3);
    }

    *xmin -= *wlenx;
    *ymin -= *wleny;
    *zmin -= *wlenz;
    *xmax += *wlenx;
    *ymax += *wleny;
    *zmax += *wlenz;

    if (*naddl == 2) {
	goto L155;
    }

/*     compute other points in grid on surface of polyhedron */

    naddm = *naddl - 2;
    xint = (*xmax - *xmin) / (doublereal) (*naddl - 1);
    yint = (*ymax - *ymin) / (doublereal) (*naddl - 1);
    zint = (*zmax - *zmin) / (doublereal) (*naddl - 1);
    ng = 16;

    i__1 = naddm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xcor = *xmin + (doublereal) i__ * xint;
	ng += 4;
	x[ng - 3] = xcor;
	y[ng - 3] = *ymin;
	z__[ng - 3] = *zmin;
	x[ng - 2] = xcor;
	y[ng - 2] = *ymin;
	z__[ng - 2] = *zmax;
	x[ng - 1] = xcor;
	y[ng - 1] = *ymax;
	z__[ng - 1] = *zmin;
	x[ng] = xcor;
	y[ng] = *ymax;
	z__[ng] = *zmax;
/* L119: */
    }

    i__1 = naddm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ycor = *ymin + (doublereal) i__ * yint;
	ng += 4;
	y[ng - 3] = ycor;
	z__[ng - 3] = *zmin;
	x[ng - 3] = *xmin;
	y[ng - 2] = ycor;
	z__[ng - 2] = *zmin;
	x[ng - 2] = *xmax;
	y[ng - 1] = ycor;
	z__[ng - 1] = *zmax;
	x[ng - 1] = *xmin;
	y[ng] = ycor;
	z__[ng] = *zmax;
	x[ng] = *xmax;
/* L120: */
    }

    i__1 = naddm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zcor = *zmin + (doublereal) i__ * zint;
	ng += 4;
	z__[ng - 3] = zcor;
	x[ng - 3] = *xmin;
	y[ng - 3] = *ymin;
	z__[ng - 2] = zcor;
	x[ng - 2] = *xmin;
	y[ng - 2] = *ymax;
	z__[ng - 1] = zcor;
	x[ng - 1] = *xmax;
	y[ng - 1] = *ymin;
	z__[ng] = zcor;
	x[ng] = *xmax;
	y[ng] = *ymax;
/* L121: */
    }

    i__1 = naddm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xcor = *xmin + (doublereal) i__ * xint;
	i__2 = naddm;
	for (j = 1; j <= i__2; ++j) {
	    ycor = *ymin + (doublereal) j * yint;
	    ng += 2;
	    x[ng - 1] = xcor;
	    y[ng - 1] = ycor;
	    z__[ng - 1] = *zmin;
	    x[ng] = xcor;
	    y[ng] = ycor;
	    z__[ng] = *zmax;
/* L125: */
	}
/* L130: */
    }

    i__1 = naddm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ycor = *ymin + (doublereal) i__ * yint;
	i__2 = naddm;
	for (j = 1; j <= i__2; ++j) {
	    zcor = *zmin + (doublereal) j * zint;
	    ng += 2;
	    y[ng - 1] = ycor;
	    z__[ng - 1] = zcor;
	    x[ng - 1] = *xmin;
	    y[ng] = ycor;
	    z__[ng] = zcor;
	    x[ng] = *xmax;
/* L135: */
	}
/* L140: */
    }

    i__1 = naddm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zcor = *zmin + (doublereal) i__ * zint;
	i__2 = naddm;
	for (j = 1; j <= i__2; ++j) {
	    xcor = *xmin + (doublereal) j * xint;
	    ng += 2;
	    z__[ng - 1] = zcor;
	    x[ng - 1] = xcor;
	    y[ng - 1] = *ymin;
	    z__[ng] = zcor;
	    x[ng] = xcor;
	    y[ng] = *ymax;
/* L145: */
	}
/* L150: */
    }

    if (ng != *irec) {
	s_stop("280", (ftnlen)3);
    }

L155:

/*     find i9 */

    i__1 = *nv;
    for (i__ = 9; i__ <= i__1; ++i__) {
	if (is[i__] != 0) {
	    goto L158;
	}
/* L157: */
    }
    s_stop("282", (ftnlen)3);
L158:
    i9 = i__;

/*     divide cube into 12 tetrahedra */

    icon[9] = 0;
    icon[10] = 5;
    icon[11] = 9;
    icon[12] = 2;
    icon[13] = i9;
    icon[14] = 7;
    icon[15] = 3;
    icon[16] = 4;

    icon[17] = 0;
    icon[18] = 1;
    icon[19] = 12;
    icon[20] = 3;
    icon[21] = i9;
    icon[22] = 2;
    icon[23] = 3;
    icon[24] = 7;

    icon[25] = 0;
    icon[26] = 5;
    icon[27] = 2;
    icon[28] = 4;
    icon[29] = i9;
    icon[30] = 2;
    icon[31] = 1;
    icon[32] = 3;

    icon[33] = 0;
    icon[34] = 3;
    icon[35] = 11;
    icon[36] = 6;
    icon[37] = i9;
    icon[38] = 5;
    icon[39] = 1;
    icon[40] = 2;

    icon[41] = 0;
    icon[42] = 6;
    icon[43] = 1;
    icon[44] = 3;
    icon[45] = i9;
    icon[46] = 3;
    icon[47] = 1;
    icon[48] = 4;

    icon[49] = 0;
    icon[50] = 4;
    icon[51] = 7;
    icon[52] = 5;
    icon[53] = i9;
    icon[54] = 4;
    icon[55] = 1;
    icon[56] = 5;

    icon[57] = 0;
    icon[58] = 6;
    icon[59] = 8;
    icon[60] = 9;
    icon[61] = i9;
    icon[62] = 8;
    icon[63] = 4;
    icon[64] = 5;

    icon[65] = 0;
    icon[66] = 10;
    icon[67] = 9;
    icon[68] = 7;
    icon[69] = i9;
    icon[70] = 8;
    icon[71] = 5;
    icon[72] = 7;

    icon[73] = 0;
    icon[74] = 7;
    icon[75] = 8;
    icon[76] = 1;
    icon[77] = i9;
    icon[78] = 7;
    icon[79] = 4;
    icon[80] = 8;

    icon[81] = 0;
    icon[82] = 11;
    icon[83] = 12;
    icon[84] = 8;
    icon[85] = i9;
    icon[86] = 7;
    icon[87] = 5;
    icon[88] = 6;

    icon[89] = 0;
    icon[90] = 4;
    icon[91] = 12;
    icon[92] = 10;
    icon[93] = i9;
    icon[94] = 6;
    icon[95] = 5;
    icon[96] = 2;

    icon[97] = 0;
    icon[98] = 2;
    icon[99] = 10;
    icon[100] = 11;
    icon[101] = i9;
    icon[102] = 6;
    icon[103] = 2;
    icon[104] = 7;

    *tetra = 12;

    wmen = *wmin;
    if (! (*delaun) && *irec > 8) {
	wmen -= *wlenw;
	i__1 = *irec;
	for (i__ = 9; i__ <= i__1; ++i__) {
	    w[i__] = wmen;
/* L160: */
	}
    }
    if (wmen < *wmin) {
	*wmin = wmen;
    }
    if (wmen > *wmax) {
	*wmax = wmen;
    }
    wman = *wmin - 10.;
    for (i__ = 1; i__ <= 8; ++i__) {
	w[i__] = wman;
/* L170: */
    }

    is[1] = 5;
    is[2] = 12;
    is[3] = 1;
    is[4] = 9;
    is[5] = 10;
    is[6] = 12;
    is[7] = 12;
    is[8] = 9;
    is[i9] = 12;

/*     test # of significant figures of nondecimal part of coordinates */

    wbig = 0.;
    if (wbig < abs(*xmax)) {
	wbig = abs(*xmax);
    }
    if (wbig < abs(*xmin)) {
	wbig = abs(*xmin);
    }
    if (wbig < abs(*ymax)) {
	wbig = abs(*ymax);
    }
    if (wbig < abs(*ymin)) {
	wbig = abs(*ymin);
    }
    if (wbig < abs(*zmax)) {
	wbig = abs(*zmax);
    }
    if (wbig < abs(*zmin)) {
	wbig = abs(*zmin);
    }
    wbig += *epz;
/*     WRITE(*,*)'COORDINATES WBIG=',WBIG */
    ibsfig = 0;
L180:
    ++ibsfig;
    wbig /= 10.;
    if (wbig >= 1.) {
	goto L180;
    }
    if (ibsfig > 9) {
	s_wsle(&io___101);
	do_lio(&c__9, &c__1, "Number of significant figures of largest ", (
		ftnlen)41);
	do_lio(&c__9, &c__1, "nondecimal part of", (ftnlen)18);
	e_wsle();
	s_wsle(&io___102);
	do_lio(&c__9, &c__1, "a point coordinate appears to be greater than \
9.", (ftnlen)48);
	e_wsle();
	s_wsle(&io___103);
	do_lio(&c__9, &c__1, "Program is terminated.", (ftnlen)22);
	e_wsle();
	s_stop("286", (ftnlen)3);
    }
    itsfig = ibsfig + *icsfig;
/*     WRITE(*,*)'ITSFIG=',ITSFIG,' IBSFIG=',IBSFIG,' ICSFIG=',ICSFIG */
    if (itsfig > 14) {
	s_wsle(&io___105);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___106);
	do_lio(&c__9, &c__1, "For this execution of the program the largest ",
		 (ftnlen)46);
	do_lio(&c__9, &c__1, "total number of", (ftnlen)15);
	e_wsle();
	s_wsle(&io___107);
	do_lio(&c__9, &c__1, "significant figures ", (ftnlen)20);
	do_lio(&c__9, &c__1, "that a coordinate requires appears to be ", (
		ftnlen)41);
	do_lio(&c__3, &c__1, (char *)&itsfig, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___108);
	do_lio(&c__9, &c__1, "Since the maximum allowed is 14, the number of "
		, (ftnlen)47);
	do_lio(&c__9, &c__1, "significant", (ftnlen)11);
	e_wsle();
	s_wsle(&io___109);
	do_lio(&c__9, &c__1, "figures of the decimal part of the coordinates "
		, (ftnlen)47);
	do_lio(&c__9, &c__1, "for this run is ", (ftnlen)16);
	e_wsle();
	s_wsle(&io___110);
	do_lio(&c__9, &c__1, "decreased accordingly.", (ftnlen)22);
	e_wsle();
	*icsfig = 14 - ibsfig;
	s_wsle(&io___111);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___112);
	do_lio(&c__9, &c__1, "Now icfig = ", (ftnlen)12);
	do_lio(&c__3, &c__1, (char *)&(*icsfig), (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___113);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
    }

/*     if a Regular tetrahedralization test # of significant figures */
/*     of nondecimal part of weights */

    if (*delaun) {
	goto L200;
    }
    wbig = 0.;
    if (wbig < abs(*wmax)) {
	wbig = abs(*wmax);
    }
    if (wbig < abs(*wmin)) {
	wbig = abs(*wmin);
    }
    wbig += *epz;
/*     WRITE(*,*)'COORDINATES WBIG=',WBIG */
    ibsfig = 0;
L190:
    ++ibsfig;
    wbig /= 10.;
    if (wbig >= 1.) {
	goto L190;
    }
    if (ibsfig > 9) {
	s_wsle(&io___114);
	do_lio(&c__9, &c__1, "Number of significant figures of largest ", (
		ftnlen)41);
	do_lio(&c__9, &c__1, "nondecimal part of", (ftnlen)18);
	e_wsle();
	s_wsle(&io___115);
	do_lio(&c__9, &c__1, "a weight appears to be greater than 9.", (
		ftnlen)38);
	e_wsle();
	s_wsle(&io___116);
	do_lio(&c__9, &c__1, "Program is terminated.", (ftnlen)22);
	e_wsle();
	s_stop("288", (ftnlen)3);
    }
    itsfig = ibsfig + *iwsfig;
/*     WRITE(*,*)'ITSFIG=',ITSFIG,' IBSFIG=',IBSFIG,' IWSFIG=',IWSFIG */
    if (itsfig > 14) {
	s_wsle(&io___117);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___118);
	do_lio(&c__9, &c__1, "For this execution of the program the largest ",
		 (ftnlen)46);
	do_lio(&c__9, &c__1, "total number of", (ftnlen)15);
	e_wsle();
	s_wsle(&io___119);
	do_lio(&c__9, &c__1, "significant figures ", (ftnlen)20);
	do_lio(&c__9, &c__1, "that a weight requires appears to be ", (ftnlen)
		37);
	do_lio(&c__3, &c__1, (char *)&itsfig, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___120);
	do_lio(&c__9, &c__1, "Since the maximum allowed is 14, the number of "
		, (ftnlen)47);
	do_lio(&c__9, &c__1, "significant", (ftnlen)11);
	e_wsle();
	s_wsle(&io___121);
	do_lio(&c__9, &c__1, "figures of the decimal part of the weights for "
		, (ftnlen)47);
	do_lio(&c__9, &c__1, "this run is ", (ftnlen)12);
	e_wsle();
	s_wsle(&io___122);
	do_lio(&c__9, &c__1, "decreased accordingly.", (ftnlen)22);
	e_wsle();
	*iwsfig = 14 - ibsfig;
	icsfi2 = *icsfig << 1;
	irsfig = icsfi2 - *iwsfig;
	if (irsfig > 9) {
	    s_wsle(&io___125);
	    do_lio(&c__9, &c__1, "In order to make this number compatible wi\
th ", (ftnlen)45);
	    do_lio(&c__9, &c__1, "that of the", (ftnlen)11);
	    e_wsle();
	    s_wsle(&io___126);
	    do_lio(&c__9, &c__1, "coordinates, the latter is also decreased ",
		     (ftnlen)42);
	    do_lio(&c__9, &c__1, "accordingly.", (ftnlen)12);
	    e_wsle();
	    *icsfig = (*iwsfig + 9) / 2;
	}
	s_wsle(&io___127);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___128);
	do_lio(&c__9, &c__1, "Now icfig = ", (ftnlen)12);
	do_lio(&c__3, &c__1, (char *)&(*icsfig), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " iwfig = ", (ftnlen)9);
	do_lio(&c__3, &c__1, (char *)&(*iwsfig), (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___129);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
    }
L200:

/*     test number of significant figures of decimal part of coordinates */

    if (*icsfig < 0 || *icsfig > 9) {
	s_wsle(&io___130);
	do_lio(&c__9, &c__1, "Number of significant figures of decimal", (
		ftnlen)40);
	e_wsle();
	s_wsle(&io___131);
	do_lio(&c__9, &c__1, "part of coordinates is out of range.", (ftnlen)
		36);
	e_wsle();
	s_wsle(&io___132);
	do_lio(&c__9, &c__1, "Program is terminated.", (ftnlen)22);
	e_wsle();
	s_stop("290", (ftnlen)3);
    }
    isclu = 1;
    dscle = 1.;
    if (*icsfig == 0) {
	goto L220;
    }
    i__1 = *icsfig;
    for (i__ = 1; i__ <= i__1; ++i__) {
	isclu *= 10;
	dscle *= 10.;
/* L210: */
    }
L220:
    if (abs(isclu) >= *mfull) {
	s_stop("295", (ftnlen)3);
    }
    decomp_(&isclp[1], &isgcl, &isclu, mhalf);
    if (isgcl != 1) {
	s_stop("297", (ftnlen)3);
    }

/*     test lengths of x, y, z-coordinates, shift and make them integers */

    dfull = (doublereal) (*mfull);
    dfill = dfull / dscle;
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ix2[i__] = 0;
	iy2[i__] = 0;
	iz2[i__] = 0;
	if ((d__1 = x[i__], abs(d__1)) < dfill) {
	    d__1 = dscle * x[i__];
	    ix[i__] = i_dnnt(&d__1);
	    if ((i__2 = ix[i__], abs(i__2)) < *mfull) {
		x[i__] = (doublereal) ix[i__] / dscle;
		goto L225;
	    }
	}
	if ((d__1 = x[i__], abs(d__1)) >= dfull) {
	    s_stop("305", (ftnlen)3);
	}
	ix[i__] = (integer) x[i__];
	if ((i__2 = ix[i__], abs(i__2)) >= *mfull) {
	    s_stop("310", (ftnlen)3);
	}
	decml = (x[i__] - d_int(&x[i__])) * dscle;
	if (abs(decml) >= dfull) {
	    s_stop("312", (ftnlen)3);
	}
	ix2[i__] = i_dnnt(&decml);
	if ((i__2 = ix2[i__], abs(i__2)) >= *mfull) {
	    s_stop("315", (ftnlen)3);
	}
	if ((i__2 = ix2[i__], abs(i__2)) == 0) {
	    x[i__] = (doublereal) ix[i__];
	    ix2[i__] = *mfull;
	} else {
	    x[i__] = (doublereal) ix[i__] + (doublereal) ix2[i__] / dscle;
	}
L225:
	if ((d__1 = y[i__], abs(d__1)) < dfill) {
	    d__1 = dscle * y[i__];
	    iy[i__] = i_dnnt(&d__1);
	    if ((i__2 = iy[i__], abs(i__2)) < *mfull) {
		y[i__] = (doublereal) iy[i__] / dscle;
		goto L230;
	    }
	}
	if ((d__1 = y[i__], abs(d__1)) >= dfull) {
	    s_stop("320", (ftnlen)3);
	}
	iy[i__] = (integer) y[i__];
	if ((i__2 = iy[i__], abs(i__2)) >= *mfull) {
	    s_stop("325", (ftnlen)3);
	}
	decml = (y[i__] - d_int(&y[i__])) * dscle;
	if (abs(decml) >= dfull) {
	    s_stop("327", (ftnlen)3);
	}
	iy2[i__] = i_dnnt(&decml);
	if ((i__2 = iy2[i__], abs(i__2)) >= *mfull) {
	    s_stop("330", (ftnlen)3);
	}
	if ((i__2 = iy2[i__], abs(i__2)) == 0) {
	    y[i__] = (doublereal) iy[i__];
	    iy2[i__] = *mfull;
	} else {
	    y[i__] = (doublereal) iy[i__] + (doublereal) iy2[i__] / dscle;
	}
L230:
	if ((d__1 = z__[i__], abs(d__1)) < dfill) {
	    d__1 = dscle * z__[i__];
	    iz[i__] = i_dnnt(&d__1);
	    if ((i__2 = iz[i__], abs(i__2)) < *mfull) {
		z__[i__] = (doublereal) iz[i__] / dscle;
		goto L235;
	    }
	}
	if ((d__1 = z__[i__], abs(d__1)) >= dfull) {
	    s_stop("335", (ftnlen)3);
	}
	iz[i__] = (integer) z__[i__];
	if ((i__2 = iz[i__], abs(i__2)) >= *mfull) {
	    s_stop("340", (ftnlen)3);
	}
	decml = (z__[i__] - d_int(&z__[i__])) * dscle;
	if (abs(decml) >= dfull) {
	    s_stop("342", (ftnlen)3);
	}
	iz2[i__] = i_dnnt(&decml);
	if ((i__2 = iz2[i__], abs(i__2)) >= *mfull) {
	    s_stop("345", (ftnlen)3);
	}
	if ((i__2 = iz2[i__], abs(i__2)) == 0) {
	    z__[i__] = (doublereal) iz[i__];
	    iz2[i__] = *mfull;
	} else {
	    z__[i__] = (doublereal) iz[i__] + (doublereal) iz2[i__] / dscle;
	}
L235:
	;
    }

/*     if a Regular tetrahedralization test number of significant figures */
/*     of decimal part of weights, test lengths of weights, shift and */
/*     make them integers */

    if (*delaun) {
	goto L300;
    }
    icsfi2 = *icsfig << 1;
    irsfig = icsfi2 - *iwsfig;
    if (*iwsfig < 0 || *iwsfig > 9 || irsfig < 0 || irsfig > 9) {
	s_wsle(&io___138);
	do_lio(&c__9, &c__1, "Either number of significant figures of decimal"
		, (ftnlen)47);
	e_wsle();
	s_wsle(&io___139);
	do_lio(&c__9, &c__1, "part of weights is out of range or it is not", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___140);
	do_lio(&c__9, &c__1, "compatible with that of point coordinates.", (
		ftnlen)42);
	e_wsle();
	s_wsle(&io___141);
	do_lio(&c__9, &c__1, "Program is terminated.", (ftnlen)22);
	e_wsle();
	s_stop("350", (ftnlen)3);
    }
    isclu = 1;
    dscli = 1.;
    if (*iwsfig == 0) {
	goto L250;
    }
    i__1 = *iwsfig;
    for (i__ = 1; i__ <= i__1; ++i__) {
	isclu *= 10;
	dscli *= 10.;
/* L240: */
    }
L250:
    if (abs(isclu) >= *mfull) {
	s_stop("360", (ftnlen)3);
    }
    decomp_(&isclw[1], &isgcl, &isclu, mhalf);
    if (isgcl != 1) {
	s_stop("363", (ftnlen)3);
    }
    dfill = dfull / dscli;
    i__1 = *nv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iw2[i__] = 0;
	if ((d__1 = w[i__], abs(d__1)) < dfill) {
	    d__1 = dscli * w[i__];
	    iw[i__] = i_dnnt(&d__1);
	    if ((i__2 = iw[i__], abs(i__2)) < *mfull) {
		w[i__] = (doublereal) iw[i__] / dscli;
		goto L260;
	    }
	}
	if ((d__1 = w[i__], abs(d__1)) >= dfull) {
	    s_stop("365", (ftnlen)3);
	}
	iw[i__] = (integer) w[i__];
	if ((i__2 = iw[i__], abs(i__2)) >= *mfull) {
	    s_stop("370", (ftnlen)3);
	}
	decml = (w[i__] - d_int(&w[i__])) * dscli;
	if (abs(decml) >= dfull) {
	    s_stop("372", (ftnlen)3);
	}
	iw2[i__] = i_dnnt(&decml);
	if ((i__2 = iw2[i__], abs(i__2)) >= *mfull) {
	    s_stop("375", (ftnlen)3);
	}
	if ((i__2 = iw2[i__], abs(i__2)) == 0) {
	    w[i__] = (doublereal) iw[i__];
	    iw2[i__] = *mfull;
	} else {
	    w[i__] = (doublereal) iw[i__] + (doublereal) iw2[i__] / dscli;
	}
L260:
	;
    }
    isclu = 1;
    if (irsfig == 0) {
	goto L290;
    }
    i__1 = irsfig;
    for (i__ = 1; i__ <= i__1; ++i__) {
	isclu *= 10;
/* L270: */
    }
L290:
    if (abs(isclu) >= *mfull) {
	s_stop("385", (ftnlen)3);
    }
    decomp_(&isclr[1], &isgcl, &isclu, mhalf);
    if (isgcl != 1) {
	s_stop("390", (ftnlen)3);
    }
L300:

/*     get cube corner directions in their integer form */

    for (i__ = 1; i__ <= 8; ++i__) {
	ixc[i__ - 1] = i_dnnt(&xc[i__ - 1]);
	iyc[i__ - 1] = i_dnnt(&yc[i__ - 1]);
	izc[i__ - 1] = i_dnnt(&zc[i__ - 1]);
/* L320: */
    }

/*     add all points to tetrahedralization */

    i10 = i9 + 1;
    if (*nv < i10) {
	goto L400;
    }
/*     write(*,*)' ' */
/*     write(*,*)'Adding all points to tetrahedralization ...' */
    itchk = 0;
    itcn1 = 0;
    itcn2 = 0;
    itcn3 = 0;
    itcn4 = 0;
    issin = i9;
    i__1 = *nv;
    for (k = i10; k <= i__1; ++k) {
/*         if(k.le.(k/1000)*1000) */
/*     *   write(*,*)'Number of points processed = ',k */
	if (is[k] == 0) {
	    goto L380;
	}
	pntins_(&x[1], &y[1], &z__[1], &w[1], &ix[1], &iy[1], &iz[1], &iw[1], 
		&ix2[1], &iy2[1], &iz2[1], &iw2[1], &icon[9], &is[1], &ifl[1],
		 &id[1], &ih[1], ihn, &k, nvmax, nhmax, tetra, mxlook, ixc, 
		iyc, izc, &issin, iftal, delaun, flphis, mhalf, mfull, &isclp[
		1], &isclw[1], &isclr[1], epz, &itchk, &itcn1, &itcn2, &itcn3,
		 &itcn4);
L380:
	;
    }
/*     WRITE(*,*)' ' */
/*     WRITE(*,*)'POINTINS ITCHK=',ITCHK */
/*     WCPNT = REAL(ITCHK)/REAL(NV-9) */
/*     WRITE(*,*)' ' */
/*     WRITE(*,*)'DURING INSERTION AVERAGE NUMBER OF ', */
/*    *          'TETRAHEDRA CHECKED PER POINT WAS ' */
/*     WRITE(*,*) WCPNT */
/*     WRITE(*,*)' ' */
/*     WRITE(*,*)'ITCN1234=',ITCN1,ITCN2,ITCN3,ITCN4 */
L400:

/*     test redundant points */

    *iredx = 0;
    if (*delaun) {
	goto L1000;
    }
    if (*nv < 10) {
	goto L1000;
    }
    i__1 = *nv;
    for (k = 9; k <= i__1; ++k) {
	if (is[k] <= 0) {
	    goto L420;
	}
L420:
	;
    }
    if (*redchk) {
	s_wsle(&io___154);
	do_lio(&c__9, &c__1, " ", (ftnlen)1);
	e_wsle();
	s_wsle(&io___155);
	do_lio(&c__9, &c__1, "Testing redundant points ...", (ftnlen)28);
	e_wsle();
    }
    itchk = 0;
    itcn1 = 0;
    itcn2 = 0;
    itcn3 = 0;
    itcn4 = 0;
    iredu = 0;
    i__1 = *nv;
    for (k = 9; k <= i__1; ++k) {
	if (is[k] >= 0) {
	    goto L450;
	}
	++iredu;
	if (! (*redchk)) {
	    goto L500;
	}
	pntred_(&x[1], &y[1], &z__[1], &w[1], &ix[1], &iy[1], &iz[1], &iw[1], 
		&ix2[1], &iy2[1], &iz2[1], &iw2[1], &icon[9], &is[1], &id[1], 
		&k, tetra, ixc, iyc, izc, iredx, &issin, iftal, delaun, 
		flphis, mhalf, mfull, &isclp[1], &isclw[1], &isclr[1], epz, &
		itchk, &itcn1, &itcn2, &itcn3, &itcn4);
L450:
	if (is[k] <= 0) {
	    goto L500;
	}
	issin = k;
L500:
	;
    }
/*     WRITE(*,*)' ' */
/*     WRITE(*,*)'POINTRED ITCHK=',ITCHK */
/*     WRITE(*,*)' ' */
/*     WRITE(*,*)'ITCN1234=',ITCN1,ITCN2,ITCN3,ITCN4 */
/*     write(*,*)' ' */
/*     write(*,*)'Number of redundant points = ',iredu */

L1000:
    return 0;
} /* poltri_ */

/* PNTINS */

/*     This subroutine will find location of new point */

/* This routine also calls routine 'sphere' for the purpose of optimizing */
/* for the locally Regular property */

/* Subroutine */ int pntins_(doublereal *xi, doublereal *yi, doublereal *zi, 
	doublereal *wi, integer *x, integer *y, integer *z__, integer *w, 
	integer *x2, integer *y2, integer *z2, integer *w2, integer *icon, 
	integer *is, integer *ifl, integer *id, integer *ih, integer *ihn, 
	integer *k, integer *nvmax, integer *nhmax, integer *tetra, integer *
	mxlook, integer *xc, integer *yc, integer *zc, integer *issin, 
	integer *iftal, logical *delaun, logical *flphis, integer *mhalf, 
	integer *mfull, integer *isclp, integer *isclw, integer *isclr, 
	doublereal *epz, integer *itchk, integer *itcn1, integer *itcn2, 
	integer *itcn3, integer *itcn4)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer ired, look, curr, side1, side2, itype;
    extern /* Subroutine */ int edgins_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, logical *, integer *, integer *, integer *, 
	    integer *, integer *), sphere_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    logical *, logical *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *), shishk_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *), sidins_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, logical *, integer *, integer *, integer *, 
	    integer *, integer *), gettet_(integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *), lkdown_(integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *), intins_(integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , logical *, logical *, integer *, integer *, integer *, integer *
	    , integer *, doublereal *);
    integer newtts;
    extern /* Subroutine */ int vrtins_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, logical *, 
	    logical *, integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --zc;
    --yc;
    --xc;
    --ih;
    --id;
    --ifl;
    --is;
    icon -= 9;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;
    --wi;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    if (*flphis) {
	itype = -1;
	look = 1;
	gettet_(&itype, k, &xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[
		1], &y2[1], &z2[1], &icon[9], &curr, &side1, &side2, &xc[1], &
		yc[1], &zc[1], mhalf, mfull, &isclp[1], epz, itchk);

L50:
	if (icon[(curr << 3) + 5] < 0) {
	    lkdown_(&icon[9], &curr, &xi[1], &yi[1], &zi[1], &x[1], &y[1], &
		    z__[1], &x2[1], &y2[1], &z2[1], &itype, k, &side1, &side2,
		     &xc[1], &yc[1], &zc[1], mhalf, mfull, &isclp[1], epz, 
		    itchk);
	    ++look;
	    if (itype == -1) {
		s_stop("410", (ftnlen)3);
	    }
	    goto L50;
	}
    } else {
	shishk_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1],
		 &z2[1], &is[1], &icon[9], &id[1], issin, k, &side1, &side2, &
		curr, iftal, &itype, tetra, &xc[1], &yc[1], &zc[1], mhalf, 
		mfull, &isclp[1], epz, itchk);
	if (itype == -1) {
	    s_stop("420", (ftnlen)3);
	}
    }

    if (itype == 1) {
	vrtins_(k, &w[1], &w2[1], &icon[9], nvmax, tetra, &curr, &is[1], &id[
		1], iftal, &side1, &ifl[1], &newtts, &ired, delaun, flphis, 
		mhalf, mfull, &isclw[1]);
	++(*itcn1);
    } else if (itype == 2) {
	intins_(k, &xi[1], &yi[1], &zi[1], &wi[1], &x[1], &y[1], &z__[1], &w[
		1], &x2[1], &y2[1], &z2[1], &w2[1], &icon[9], &ih[1], ihn, 
		nhmax, nvmax, tetra, &curr, &is[1], &ifl[1], &newtts, &ired, 
		delaun, flphis, mhalf, mfull, &isclp[1], &isclw[1], &isclr[1],
		 epz);
	++(*itcn2);
    } else if (itype == 3) {
	edgins_(k, &x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		1], &icon[9], &ih[1], ihn, nhmax, nvmax, tetra, &curr, &is[1],
		 &side1, &side2, &ifl[1], &newtts, &ired, delaun, flphis, 
		mhalf, mfull, &isclp[1], &isclw[1], &isclr[1]);
	++(*itcn3);
    } else if (itype == 4) {
	sidins_(k, &x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		1], &icon[9], &ih[1], ihn, nhmax, nvmax, tetra, &curr, &is[1],
		 &side1, &ifl[1], &newtts, &ired, delaun, flphis, mhalf, 
		mfull, &isclp[1], &isclw[1], &isclr[1]);
	++(*itcn4);
    } else {
	s_stop("430", (ftnlen)3);
    }
    if (ired == 1) {
	goto L1000;
    }
    *issin = *k;

/*     optimize for Regular/Delaunay property */

    sphere_(k, &icon[9], &ifl[1], &ih[1], ihn, nhmax, &newtts, &xi[1], &yi[1],
	     &zi[1], &wi[1], &x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &
	    z2[1], &w2[1], tetra, &is[1], nvmax, &xc[1], &yc[1], &zc[1], 
	    delaun, flphis, mhalf, mfull, &isclp[1], &isclw[1], &isclr[1], 
	    epz);
L1000:

    if (look > *mxlook) {
	*mxlook = look;
    }

    return 0;
} /* pntins_ */

/* PNTRED */

/* This subroutine will test a redundant point for redundancy */


/* Subroutine */ int pntred_(doublereal *xi, doublereal *yi, doublereal *zi, 
	doublereal *wi, integer *x, integer *y, integer *z__, integer *w, 
	integer *x2, integer *y2, integer *z2, integer *w2, integer *icon, 
	integer *is, integer *id, integer *k, integer *tetra, integer *xc, 
	integer *yc, integer *zc, integer *idmax, integer *issin, integer *
	iftal, logical *delaun, logical *flphis, integer *mhalf, integer *
	mfull, integer *isclp, integer *isclw, integer *isclr, doublereal *
	epz, integer *itchk, integer *itcn1, integer *itcn2, integer *itcn3, 
	integer *itcn4)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer a, b, c__, d__, i__, curr;
    doublereal xctr, yctr, zctr;
    integer side1, side2, site1, site2, site3, site4;
    extern /* Subroutine */ int ctrad_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, logical 
	    *, integer *);
    integer itide;
    doublereal tdist;
    integer itype;
    extern /* Subroutine */ int iqsig1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, integer *), iqsig2_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, logical *, integer *);
    integer fndsit;
    extern /* Subroutine */ int bisphr_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *), shishk_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *), iqsign_(integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, logical *, integer *), gettet_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *), 
	    lkdown_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *), iwsign_(integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *), sitord_(
	    integer *, integer *, integer *);
    integer ipossi;



    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --zc;
    --yc;
    --xc;
    --id;
    --is;
    icon -= 9;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;
    --wi;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    if (*flphis) {
	itype = -1;
	gettet_(&itype, k, &xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[
		1], &y2[1], &z2[1], &icon[9], &curr, &side1, &side2, &xc[1], &
		yc[1], &zc[1], mhalf, mfull, &isclp[1], epz, itchk);

L50:
	if (icon[(curr << 3) + 5] < 0) {
	    lkdown_(&icon[9], &curr, &xi[1], &yi[1], &zi[1], &x[1], &y[1], &
		    z__[1], &x2[1], &y2[1], &z2[1], &itype, k, &side1, &side2,
		     &xc[1], &yc[1], &zc[1], mhalf, mfull, &isclp[1], epz, 
		    itchk);
	    if (itype == -1) {
		s_stop("510", (ftnlen)3);
	    }
	    goto L50;
	}
    } else {
	shishk_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1],
		 &z2[1], &is[1], &icon[9], &id[1], issin, k, &side1, &side2, &
		curr, iftal, &itype, tetra, &xc[1], &yc[1], &zc[1], mhalf, 
		mfull, &isclp[1], epz, itchk);
	if (itype == -1) {
	    s_stop("515", (ftnlen)3);
	}
    }

    if (itype == 1) {
	++(*itcn1);
	site1 = icon[side1 + 4 + (curr << 3)];
/*        itide = 1 */
/*        if(w(k) .gt. w(site1)) itide =-1 */
	iwsign_(&w[1], &w2[1], &site1, k, mhalf, mfull, &isclw[1], &itide);
    } else if (itype == 2) {
	++(*itcn2);
	a = icon[(curr << 3) + 5];
	b = icon[(curr << 3) + 6];
	c__ = icon[(curr << 3) + 7];
	d__ = icon[(curr << 3) + 8];
	if (a <= 8 || b <= 8 || c__ <= 8 || d__ <= 8) {
	    s_stop("520", (ftnlen)3);
	}
	ctrad_(&xi[1], &yi[1], &zi[1], &wi[1], &xctr, &yctr, &zctr, &a, &b, &
		c__, &d__, epz, delaun, &ipossi);
	if (ipossi == 1) {
	    goto L60;
	}
	bisphr_(&xi[1], &yi[1], &zi[1], &wi[1], k, &a, epz, &xctr, &yctr, &
		zctr, &tdist, delaun, &ipossi);
	if (ipossi == 1) {
	    goto L60;
	}
	itide = 1;
	if (tdist > 0.) {
	    itide = -1;
	}
	goto L1000;
L60:
	iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&a, &b, &c__, &d__, k, mhalf, mfull, &isclp[1], &isclw[1], &
		isclr[1], delaun, &itide);
    } else if (itype == 3) {
	++(*itcn3);
	fndsit = 0;
	for (i__ = 5; i__ <= 8; ++i__) {
	    if (fndsit == 0) {
		if (i__ == side1 + 4 || i__ == side2 + 4) {
		    goto L100;
		} else {
		    site1 = icon[i__ + (curr << 3)];
		    fndsit = 1;
		}
	    } else {
		if (i__ == side1 + 4 || i__ == side2 + 4) {
		    goto L100;
		} else {
		    site2 = icon[i__ + (curr << 3)];
		    goto L150;
		}
	    }
L100:
	    ;
	}
	s_stop("530", (ftnlen)3);
L150:

	if (site1 <= 8 || site2 <= 8) {
	    s_stop("540", (ftnlen)3);
	}
	iqsig1_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&site1, &site2, k, mhalf, mfull, &isclp[1], &isclw[1], &isclr[
		1], delaun, &itide);
    } else if (itype == 4) {
	++(*itcn4);
	site1 = icon[side1 + 4 + (curr << 3)];
	sitord_(&icon[9], &site1, &curr);

	site2 = icon[(curr << 3) + 6];
	site3 = icon[(curr << 3) + 7];
	site4 = icon[(curr << 3) + 8];

	if (site2 <= 8 || site3 <= 8 || site4 <= 8) {
	    s_stop("550", (ftnlen)3);
	}
	iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&site2, &site3, &site4, k, mhalf, mfull, &isclp[1], &isclw[1],
		 &isclr[1], delaun, &itide);
    } else {
	s_stop("560", (ftnlen)3);
    }

L1000:
    if (itide < 0) {
	++(*idmax);
    }

    return 0;
} /* pntred_ */

/* GETTET */

/* This subroutine will test each of the 1st 12 tetrahedra to find */
/* where new point is located.  It'll do so by calling 'gette2'. */

/* Subroutine */ int gettet_(integer *itype, integer *k, doublereal *xi, 
	doublereal *yi, doublereal *zi, integer *x, integer *y, integer *z__, 
	integer *x2, integer *y2, integer *z2, integer *icon, integer *curr, 
	integer *side1, integer *side2, integer *xc, integer *yc, integer *zc,
	 integer *mhalf, integer *mfull, integer *isclp, doublereal *epz, 
	integer *itchk)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer a, b, c__, d__, flag__;
    extern /* Subroutine */ int gette2_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *), 
	    vrtord_(integer *, integer *, integer *, integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    icon -= 9;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    for (*curr = 1; *curr <= 12; ++(*curr)) {
	++(*itchk);
	a = (i__1 = icon[(*curr << 3) + 5], abs(i__1));
	b = icon[(*curr << 3) + 6];
	c__ = icon[(*curr << 3) + 7];
	d__ = icon[(*curr << 3) + 8];

	flag__ = icon[(*curr << 3) + 5];
	icon[(*curr << 3) + 5] = a;

	vrtord_(&icon[9], curr, &a, &b, &c__, &d__);
	if (flag__ < 0) {
	    icon[(*curr << 3) + 5] = -a;
	}
	if (a <= 8 || b > 8 || c__ > 8 || d__ > 8) {
	    s_stop("610", (ftnlen)3);
	}

	gette2_(&a, &b, &c__, &d__, &xi[1], &yi[1], &zi[1], &x[1], &y[1], &
		z__[1], &x2[1], &y2[1], &z2[1], itype, k, side1, side2, &
		flag__, &xc[1], &yc[1], &zc[1], mhalf, mfull, &isclp[1], epz);
	if (*itype != -1) {
	    goto L200;
	}

/* L100: */
    }
    s_stop("620", (ftnlen)3);
L200:
    return 0;
} /* gettet_ */

/* LKDOWN */

/* This subroutine will traverse thru children of curr until point */
/* is found to be in one of them */

/* Subroutine */ int lkdown_(integer *icon, integer *curr, doublereal *xi, 
	doublereal *yi, doublereal *zi, integer *x, integer *y, integer *z__, 
	integer *x2, integer *y2, integer *z2, integer *itype, integer *k, 
	integer *side1, integer *side2, integer *xc, integer *yc, integer *zc,
	 integer *mhalf, integer *mfull, integer *isclp, doublereal *epz, 
	integer *itchk)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer a, b, c__, d__, i__, flag__;
    extern /* Subroutine */ int gette2_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *);
    integer newcur;
    extern /* Subroutine */ int vrtord_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);



/*     test children of current tetrahedron */

    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;
    icon -= 9;

    /* Function Body */
    *itype = -1;
    for (i__ = 1; i__ <= 4; ++i__) {
	newcur = icon[i__ + (*curr << 3)];
	if (newcur <= 0) {
	    goto L100;
	}
	++(*itchk);

	a = (i__1 = icon[(newcur << 3) + 5], abs(i__1));
	b = icon[(newcur << 3) + 6];
	c__ = icon[(newcur << 3) + 7];
	d__ = icon[(newcur << 3) + 8];

	flag__ = icon[(newcur << 3) + 5];
	icon[(newcur << 3) + 5] = a;

	vrtord_(&icon[9], &newcur, &a, &b, &c__, &d__);
	if (flag__ < 0) {
	    icon[(newcur << 3) + 5] = -a;
	}
	if (a <= 8) {
	    s_stop("710", (ftnlen)3);
	}

	gette2_(&a, &b, &c__, &d__, &xi[1], &yi[1], &zi[1], &x[1], &y[1], &
		z__[1], &x2[1], &y2[1], &z2[1], itype, k, side1, side2, &
		flag__, &xc[1], &yc[1], &zc[1], mhalf, mfull, &isclp[1], epz);
	if (*itype == -1) {
	    goto L100;
	}
	*curr = newcur;
	goto L1000;
L100:
	;
    }

L1000:
    return 0;
} /* lkdown_ */

/* GETTE2 */

/* This subroutine will check for each tetra, if the point is equal to an */
/* existing vertex, inside (interior, edge, side), or outside curr tetra. */

/* Subroutine */ int gette2_(integer *a, integer *b, integer *c__, integer *
	d__, doublereal *xi, doublereal *yi, doublereal *zi, integer *x, 
	integer *y, integer *z__, integer *x2, integer *y2, integer *z2, 
	integer *itype, integer *k, integer *side1, integer *side2, integer *
	flag__, integer *xc, integer *yc, integer *zc, integer *mhalf, 
	integer *mfull, integer *isclp, doublereal *epz)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer ifn, iside[4], ipout;
    extern /* Subroutine */ int ipsig3_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), ipsig4_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), ipsig6_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), ipsign_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     irsign_(doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *);



/*     determine position of point k relative to facets of tetrahedron */

    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    if (*b <= 8 || *c__ <= 8 || *d__ <= 8) {
	goto L100;
    }

    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, b, d__, c__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[0] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, a, c__, d__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[1] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, a, d__, b, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[2] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, a, b, c__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[3] = ipout;

/*     point k is not in tetrahedron */

L50:
    if (iside[0] < 0 || iside[1] < 0 || iside[2] < 0 || iside[3] < 0) {
	goto L1000;
    }

    if (*flag__ < 0) {
	*itype = 0;
	goto L1000;
    }

/*     point k is in the interior of tetrahedron */

    if (iside[0] > 0 && iside[1] > 0 && iside[2] > 0 && iside[3] > 0) {
	*itype = 2;
	goto L1000;
    }

/*     unacceptable situation */

    if (iside[0] == 0 && iside[1] == 0 && iside[2] == 0 && iside[3] == 0) {
	s_stop("805", (ftnlen)3);
    }

/*     point k is a vertex of tetrahedron */

    if (iside[0] == 0 && iside[1] == 0 && iside[2] == 0) {
	*itype = 1;
	*side1 = 4;
	goto L1000;
    } else if (iside[0] == 0 && iside[1] == 0 && iside[3] == 0) {
	*itype = 1;
	*side1 = 3;
	goto L1000;
    } else if (iside[0] == 0 && iside[2] == 0 && iside[3] == 0) {
	*itype = 1;
	*side1 = 2;
	goto L1000;
    } else if (iside[1] == 0 && iside[2] == 0 && iside[3] == 0) {
	*itype = 1;
	*side1 = 1;
	goto L1000;
    }

/*     point k is in the interior of an edge of tetrahedron */

    if (iside[0] == 0 && iside[1] == 0) {
	*itype = 3;
	*side1 = 1;
	*side2 = 2;
	goto L1000;
    } else if (iside[0] == 0 && iside[2] == 0) {
	*itype = 3;
	*side1 = 1;
	*side2 = 3;
	goto L1000;
    } else if (iside[0] == 0 && iside[3] == 0) {
	*itype = 3;
	*side1 = 1;
	*side2 = 4;
	goto L1000;
    } else if (iside[1] == 0 && iside[2] == 0) {
	*itype = 3;
	*side1 = 2;
	*side2 = 3;
	goto L1000;
    } else if (iside[1] == 0 && iside[3] == 0) {
	*itype = 3;
	*side1 = 2;
	*side2 = 4;
	goto L1000;
    } else if (iside[2] == 0 && iside[3] == 0) {
	*itype = 3;
	*side1 = 3;
	*side2 = 4;
	goto L1000;
    }

/*     point k is in the interior of a facet of tetrahedron */

    *itype = 4;
    if (iside[0] == 0) {
	*side1 = 1;
    } else if (iside[1] == 0) {
	*side1 = 2;
    } else if (iside[2] == 0) {
	*side1 = 3;
    } else if (iside[3] == 0) {
	*side1 = 4;
    } else {
	s_stop("807", (ftnlen)3);
    }
    goto L1000;

/*     there is at least one artificial vertex */

L100:
    if (*b <= 8) {
	iside[0] = 1;
	goto L120;
    } else if (*d__ <= 8 && *c__ <= 8) {
	ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], d__, c__, k, b, mhalf, mfull, &isclp[1], &ipout);
	iside[0] = ipout;
	if (iside[0] != 0) {
	    goto L120;
	}
	ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], d__, c__, k, b, mhalf, mfull, &isclp[1], &ipout);
	iside[0] = ipout;
	goto L120;
    } else if (*d__ <= 8) {
	ifn = 0;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], b, c__, k, d__, mhalf, mfull, &isclp[1], &ifn, &
		ipout);
	iside[0] = ipout;
    } else {
	ifn = 1;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], b, d__, k, c__, mhalf, mfull, &isclp[1], &ifn, &
		ipout);
	iside[0] = ipout;
    }
    if (iside[0] != 0) {
	goto L120;
    }
    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], b, d__, c__, k, 
	    mhalf, mfull, &isclp[1], &ipout);
    iside[0] = ipout;
L120:

    if (*c__ <= 8 && *d__ <= 8) {
	ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], c__, d__, k, a, mhalf, mfull, &isclp[1], &ipout);
	iside[1] = ipout;
	if (iside[1] != 0) {
	    goto L140;
	}
	ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], c__, d__, k, a, mhalf, mfull, &isclp[1], &ipout);
	iside[1] = ipout;
	goto L140;
    } else if (*c__ <= 8) {
	ifn = 0;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], a, d__, k, c__, mhalf, mfull, &isclp[1], &ifn, &
		ipout);
	iside[1] = ipout;
    } else {
	ifn = 1;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], a, c__, k, d__, mhalf, mfull, &isclp[1], &ifn, &
		ipout);
	iside[1] = ipout;
    }
    if (iside[1] != 0) {
	goto L140;
    }
    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], a, c__, d__, k, 
	    mhalf, mfull, &isclp[1], &ipout);
    iside[1] = ipout;
L140:

    if (*d__ > 8) {
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], a, d__, b, k, 
		mhalf, mfull, &isclp[1], &ipout);
	iside[2] = ipout;
	goto L160;
    } else if (*b <= 8) {
	ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], d__, b, k, a, mhalf, mfull, &isclp[1], &ipout);
	iside[2] = ipout;
	if (iside[2] != 0) {
	    goto L160;
	}
	ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], d__, b, k, a, mhalf, mfull, &isclp[1], &ipout);
	iside[2] = ipout;
	goto L160;
    } else {
	ifn = 0;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], a, b, k, d__, mhalf, mfull, &isclp[1], &ifn, &ipout);
	iside[2] = ipout;
    }
    if (iside[2] != 0) {
	goto L160;
    }
    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], a, d__, b, k, 
	    mhalf, mfull, &isclp[1], &ipout);
    iside[2] = ipout;
L160:

    if (*c__ > 8) {
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], a, b, c__, k, 
		mhalf, mfull, &isclp[1], &ipout);
	iside[3] = ipout;
	goto L180;
    } else if (*b <= 8) {
	ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], b, c__, k, a, mhalf, mfull, &isclp[1], &ipout);
	iside[3] = ipout;
	if (iside[3] != 0) {
	    goto L180;
	}
	ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], b, c__, k, a, mhalf, mfull, &isclp[1], &ipout);
	iside[3] = ipout;
	goto L180;
    } else {
	ifn = 1;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], a, b, k, c__, mhalf, mfull, &isclp[1], &ifn, &ipout);
	iside[3] = ipout;
    }
    if (iside[3] != 0) {
	goto L180;
    }
    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], a, b, c__, k, 
	    mhalf, mfull, &isclp[1], &ipout);
    iside[3] = ipout;
L180:
    goto L50;

L1000:
    return 0;
} /* gette2_ */

/* SHISHK */

/*     shishkebab routines */

/*     subroutine shishk to - */

/*     move from a vertex in the tetrahedralization to a tetrahedron */
/*     that contains a point, and identify the type of location of */
/*     the point with respect to the tetrahedron */

/* Subroutine */ int shishk_(doublereal *xi, doublereal *yi, doublereal *zi, 
	integer *x, integer *y, integer *z__, integer *x2, integer *y2, 
	integer *z2, integer *is, integer *icon, integer *id, integer *ileft, 
	integer *k, integer *side1, integer *side2, integer *iscur, integer *
	iftal, integer *itype, integer *ivnxt, integer *xc, integer *yc, 
	integer *zc, integer *mhalf, integer *mfull, integer *isclp, 
	doublereal *epz, integer *itchk)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer a, b, c__, d__, i__, site0, site1, isadj, ilift, isini, imist, 
	    islst;
    extern /* Subroutine */ int fcedge_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), fcfind_(doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *), fctest_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *), 
	    reordr_(integer *, integer *, integer *, integer *), sitord_(
	    integer *, integer *, integer *);



/*     reinitialize array id if necessary */

    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --id;
    icon -= 9;
    --is;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    if (*iftal > 10000000) {
	*iftal = 0;
	i__1 = *ivnxt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    id[i__] = 0;
/* L50: */
	}
    }

    if (*ileft <= 8) {
	s_stop("911", (ftnlen)3);
    }
    a = *ileft;
L100:

/*     find tetrahedron with point a as a vertex for which the ray with */
/*     origin point a and through point k intersects the facet of the */
/*     tetrahedron opposite to point a */

    *itype = 0;
    ++(*iftal);
    *iscur = is[a];
    if (*iscur <= 0 || *iscur > *ivnxt) {
	s_stop("912", (ftnlen)3);
    }
    isini = *iscur;

/*     reorder isini so that vertex a equals icon(5,isini) */

    sitord_(&icon[9], &a, &isini);

/*     test current facet */

L400:
    b = icon[(*iscur << 3) + 6];
    c__ = icon[(*iscur << 3) + 7];
    d__ = icon[(*iscur << 3) + 8];
    id[*iscur] = *iftal;

    ++(*itchk);
    fctest_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], &xc[1], &yc[1], &zc[1], itype, k, &imist, &a, &b, &c__, &
	    d__, side1, side2, mhalf, mfull, &isclp[1], epz);
    if (*itype > 0) {
	goto L9000;
    }
    if (*itype == 0) {
	goto L500;
    }
    if (*itype == -2) {
	goto L1100;
    }
    if (*itype == -3) {
	a = imist;
	goto L100;
    } else if (*itype == -4) {
	site0 = a;
	site1 = imist;
	goto L2000;
    } else {
	s_stop("913", (ftnlen)3);
    }

/*     obtain next tetrahedron with point a as a vertex */

L500:
    isadj = (i__1 = icon[(*iscur << 3) + 2], abs(i__1));
    if (isadj <= 0 || isadj > *ivnxt) {
	s_stop("914", (ftnlen)3);
    }
    if (id[isadj] == *iftal) {
	goto L600;
    }
    ilift = icon[(*iscur << 3) + 8];
    goto L900;
L600:
    isadj = (i__1 = icon[(*iscur << 3) + 3], abs(i__1));
    if (isadj <= 0 || isadj > *ivnxt) {
	s_stop("915", (ftnlen)3);
    }
    if (id[isadj] == *iftal) {
	goto L700;
    }
    ilift = icon[(*iscur << 3) + 6];
    goto L900;
L700:
    isadj = (i__1 = icon[(*iscur << 3) + 4], abs(i__1));
    if (isadj <= 0 || isadj > *ivnxt) {
	s_stop("916", (ftnlen)3);
    }
    if (*iscur == isini) {
	goto L800;
    }
    if ((i__1 = icon[(isadj << 3) + 3], abs(i__1)) == *iscur) {
	*iscur = isadj;
	goto L700;
    } else if ((i__1 = icon[(isadj << 3) + 2], abs(i__1)) == *iscur) {
	*iscur = isadj;
	goto L600;
    } else if ((i__1 = icon[(isadj << 3) + 4], abs(i__1)) == *iscur) {
	if (isadj != isini) {
	    s_stop("917", (ftnlen)3);
	}
	goto L1000;
    } else {
	s_stop("918", (ftnlen)3);
    }
L800:
    if (id[isadj] == *iftal) {
	goto L1000;
    }
    ilift = icon[(*iscur << 3) + 7];

/*     reorder isadj so that a equals icon(5,isadj) and ilift */
/*     equals icon(6,isadj) */

L900:
    reordr_(&icon[9], &a, &ilift, &isadj);
    *iscur = isadj;
    goto L400;

/*     can not find intersected tetrahedron */

L1000:
    s_stop("919", (ftnlen)3);

/*     obtain next tetrahedron along line segment as it crosses a facet */

L1100:
    islst = *iscur;
    isadj = (i__1 = icon[(*iscur << 3) + 1], abs(i__1));
    if (isadj <= 0 || isadj > *ivnxt) {
	s_stop("921", (ftnlen)3);
    }
    *iscur = isadj;
    if ((i__1 = icon[(*iscur << 3) + 1], abs(i__1)) == islst) {
	ilift = icon[(*iscur << 3) + 5];
    } else if ((i__1 = icon[(*iscur << 3) + 2], abs(i__1)) == islst) {
	ilift = icon[(*iscur << 3) + 6];
    } else if ((i__1 = icon[(*iscur << 3) + 3], abs(i__1)) == islst) {
	ilift = icon[(*iscur << 3) + 7];
    } else if ((i__1 = icon[(*iscur << 3) + 4], abs(i__1)) == islst) {
	ilift = icon[(*iscur << 3) + 8];
    } else {
	s_stop("922", (ftnlen)3);
    }

/*     obtain opposite facet of tetrahedron intersected by line */
/*     segment */

    ++(*itchk);
    fcfind_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], &xc[1], &yc[1], &zc[1], itype, ileft, k, &ilift, &imist, &
	    b, &c__, &d__, side1, side2, mhalf, mfull, &isclp[1], epz);
    if (*itype > 0) {
	reordr_(&icon[9], &ilift, &b, iscur);
	goto L9000;
    } else if (*itype == -2) {
	sitord_(&icon[9], &imist, iscur);
	b = icon[(*iscur << 3) + 6];
	c__ = icon[(*iscur << 3) + 7];
	d__ = icon[(*iscur << 3) + 8];
	goto L1100;
    } else if (*itype == -3) {
	a = ilift;
	goto L100;
    } else if (*itype == -4) {
	if (imist == b) {
	    site0 = c__;
	} else if (imist == c__) {
	    site0 = d__;
	} else {
	    site0 = b;
	}
	site1 = imist;
	goto L2000;
    } else {
	s_stop("923", (ftnlen)3);
    }

/*     obtain next tetrahedron along line segment as it crosses an edge */

L2000:
    fcedge_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], itype, ileft, k, &icon[9], iscur, &imist, ivnxt, &site0, &
	    site1, side1, side2, mhalf, mfull, &isclp[1], itchk);
    if (*itype > 0) {
	goto L9000;
    }
    if (*itype == -2) {
	sitord_(&icon[9], &imist, iscur);
	b = icon[(*iscur << 3) + 6];
	c__ = icon[(*iscur << 3) + 7];
	d__ = icon[(*iscur << 3) + 8];
	goto L1100;
    } else if (*itype == -3) {
	a = imist;
	goto L100;
    } else if (*itype == -4) {
	if (imist == site1) {
	    site0 = icon[(*iscur << 3) + 7];
	} else if (imist == site0) {
	    site0 = site1;
	}
	site1 = imist;
	goto L2000;
    } else {
	s_stop("924", (ftnlen)3);
    }

L9000:
    return 0;
} /* shishk_ */

/* IRSIGN */

/*     subroutine to determine position of point site0 with respect */
/*     to the plane spanned by points site1, site2, site3 */

/* Subroutine */ int irsign_(doublereal *xi, doublereal *yi, doublereal *zi, 
	integer *x, integer *y, integer *z__, integer *x2, integer *y2, 
	integer *z2, integer *site0, integer *site1, integer *site2, integer *
	site3, integer *mhalf, integer *mfull, integer *isclp, doublereal *
	epz, integer *ipout)
{
    doublereal dist;
    extern /* Subroutine */ int dstnce_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *), ipsign_(integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, integer *, integer *)
	    ;
    integer ipossi;



    /* Parameter adjustments */
    --isclp;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    dstnce_(&xi[1], &yi[1], &zi[1], site1, site2, site3, epz, site0, &dist, &
	    ipossi);
    if (ipossi == 0) {
	*ipout = 1;
	if (dist < 0.) {
	    *ipout = -1;
	}
    } else {
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], site1, site2, 
		site3, site0, mhalf, mfull, &isclp[1], ipout);
    }

    return 0;
} /* irsign_ */

/* FCTEST */

/*     This subroutine will test whether a ray with origin a vertex of */
/*     a tetrahedron intersects the facet opposite the vertex of the */
/*     tetrahedron and whether a point in the interior of the ray is */
/*     contained in the tetrahedron */

/* Subroutine */ int fctest_(doublereal *xi, doublereal *yi, doublereal *zi, 
	integer *x, integer *y, integer *z__, integer *x2, integer *y2, 
	integer *z2, integer *xc, integer *yc, integer *zc, integer *itype, 
	integer *k, integer *imist, integer *a, integer *b, integer *c__, 
	integer *d__, integer *side1, integer *side2, integer *mhalf, integer 
	*mfull, integer *isclp, doublereal *epz)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer iside[4], ipout, iasign;
    extern /* Subroutine */ int artsig_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), irsign_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *), pntype_(
	    integer *, integer *, integer *, integer *);



/*     determine whether ray with origin point a and through point k */
/*     intersects facet of current tetrahedron opposite to point a */

    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    if (*b <= 8 || *c__ <= 8 || *d__ <= 8) {
	goto L100;
    }

    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, a, c__, d__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[1] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, a, d__, b, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[2] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, a, b, c__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[3] = ipout;

    if (iside[1] < 0 || iside[2] < 0 || iside[3] < 0) {
	goto L1000;
    }

/*     determine whether point k is in tetrahedron */

    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, b, d__, c__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[0] = ipout;
    if (iside[0] < 0) {
	goto L500;
    }

L50:
    pntype_(iside, itype, side1, side2);
    goto L1000;

/*     there is at least one artificial vertex */

L100:
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], a, c__, d__, k, mhalf, mfull, &isclp[1], &iasign);
    iside[1] = iasign;
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], a, d__, b, k, mhalf, mfull, &isclp[1], &iasign);
    iside[2] = iasign;
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], a, b, c__, k, mhalf, mfull, &isclp[1], &iasign);
    iside[3] = iasign;

    if (iside[1] < 0 || iside[2] < 0 || iside[3] < 0) {
	goto L1000;
    }

/*     determine whether point k is in tetrahedron */

    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], b, d__, c__, k, mhalf, mfull, &isclp[1], &iasign);
    iside[0] = iasign;
    if (iside[0] < 0) {
	goto L500;
    }
    goto L50;

/*     ray intersects facet but point k is not in tetrahedron */

L500:

/*     ray intersects interior of facet */

    if (iside[1] > 0 && iside[2] > 0 && iside[3] > 0) {
	*itype = -2;
	goto L1000;
    }

    if (iside[1] == 0 && iside[2] == 0 && iside[3] == 0) {
	s_stop("931", (ftnlen)3);
    }

/*     ray intersects a vertex of facet */

    if (iside[1] == 0 && iside[2] == 0) {
	*itype = -3;
	*imist = *d__;
	goto L1000;
    } else if (iside[1] == 0 && iside[3] == 0) {
	*itype = -3;
	*imist = *c__;
	goto L1000;
    } else if (iside[2] == 0 && iside[3] == 0) {
	*itype = -3;
	*imist = *b;
	goto L1000;
    }

/*     ray intersects the interior of an edge of facet */

    *itype = -4;
    if (iside[1] == 0) {
	*imist = *c__;
    } else if (iside[2] == 0) {
	*imist = *d__;
    } else if (iside[3] == 0) {
	*imist = *b;
    } else {
	s_stop("932", (ftnlen)3);
    }

L1000:
    return 0;
} /* fctest_ */

/* FCFIND */

/*     This subroutine tests whether a point on a ray that intersects */
/*     the interior of a facet of a tetrahedron is in the tetrahedron */
/*     and if not finds other facet of the tetrahedron intersected by */
/*     the ray */

/* Subroutine */ int fcfind_(doublereal *xi, doublereal *yi, doublereal *zi, 
	integer *x, integer *y, integer *z__, integer *x2, integer *y2, 
	integer *z2, integer *xc, integer *yc, integer *zc, integer *itype, 
	integer *ileft, integer *k, integer *ilift, integer *imist, integer *
	b, integer *c__, integer *d__, integer *side1, integer *side2, 
	integer *mhalf, integer *mfull, integer *isclp, doublereal *epz)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer idot1, idot2, idot3, idut1, idut2, idut3, iside[4], ipout, iasign;
    extern /* Subroutine */ int artsig_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), irsign_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *), pntype_(
	    integer *, integer *, integer *, integer *);



/*     determine whether point k is in tetrahedron */

    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    if (*b <= 8 || *c__ <= 8 || *d__ <= 8 || *ilift <= 8) {
	goto L100;
    }

    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, ilift, d__, c__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[1] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, ilift, c__, b, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[2] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, ilift, b, d__, mhalf, mfull, &isclp[1], epz, &ipout);
    iside[3] = ipout;

    if (iside[1] < 0 || iside[2] < 0 || iside[3] < 0) {
	goto L200;
    }

L50:
    iside[0] = 1;
    pntype_(iside, itype, side1, side2);
    goto L1000;

/*     there is at least one artificial vertex */

L100:
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ilift, d__, c__, k, mhalf, mfull, &isclp[1], &iasign);
    iside[1] = iasign;
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ilift, b, d__, k, mhalf, mfull, &isclp[1], &iasign);
    iside[3] = iasign;
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ilift, c__, b, k, mhalf, mfull, &isclp[1], &iasign);
    iside[2] = iasign;

    if (iside[1] < 0 || iside[2] < 0 || iside[3] < 0) {
	goto L200;
    }
    goto L50;

/*     k is not in tetrahedron */

/*     determine position of ilift with repect to current situation */

L200:
    if (*b <= 8 || *c__ <= 8 || *d__ <= 8 || *ilift <= 8) {
	goto L300;
    }

    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], ileft, ilift, c__, b, mhalf, mfull, &isclp[1], epz, &idut1)
	    ;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], ileft, ilift, d__, c__, mhalf, mfull, &isclp[1], epz, &
	    idut2);
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], ileft, ilift, b, d__, mhalf, mfull, &isclp[1], epz, &idut3)
	    ;

    if (idut1 <= 0 || idut2 <= 0 || idut3 <= 0) {
	goto L700;
    }
    goto L400;

/*     there is at least one artificial vertex */

L300:
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ilift, c__, b, ileft, mhalf, mfull, &isclp[1], &idut1);
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ilift, d__, c__, ileft, mhalf, mfull, &isclp[1], &idut2);
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ilift, b, d__, ileft, mhalf, mfull, &isclp[1], &idut3);

    if (idut1 <= 0 || idut2 <= 0 || idut3 <= 0) {
	goto L700;
    }

/*     ilift, ileft, b, d, c, form a strictly convex set */

L400:
    if (*b <= 8 || *c__ <= 8 || *d__ <= 8 || *ilift <= 8) {
	goto L500;
    }

    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, ileft, b, ilift, mhalf, mfull, &isclp[1], epz, &idot1);
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, ileft, c__, ilift, mhalf, mfull, &isclp[1], epz, &idot2)
	    ;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], k, ileft, d__, ilift, mhalf, mfull, &isclp[1], epz, &idot3)
	    ;
    goto L600;

L500:
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ileft, b, ilift, k, mhalf, mfull, &isclp[1], &idot1);
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ileft, c__, ilift, k, mhalf, mfull, &isclp[1], &idot2);
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], ileft, d__, ilift, k, mhalf, mfull, &isclp[1], &idot3);

L600:
    *itype = -2;
    if (idot1 < 0 && idot2 > 0) {
	*imist = *d__;
    } else if (idot2 < 0 && idot3 > 0) {
	*imist = *b;
    } else if (idot3 < 0 && idot1 > 0) {
	*imist = *c__;
    } else if (idot1 == 0 && idot2 == 0 && idot3 == 0) {
	*itype = -3;
    } else if (idot1 == 0) {
	*itype = -4;
	*imist = *b;
    } else if (idot2 == 0) {
	*itype = -4;
	*imist = *c__;
    } else if (idot3 == 0) {
	*itype = -4;
	*imist = *d__;
    } else {
	s_stop("951", (ftnlen)3);
    }
    goto L1000;

L700:
    if (idut1 <= 0 && idut2 <= 0 && idut3 <= 0) {
	s_stop("952", (ftnlen)3);
    }
    *itype = -2;
    if (idut1 <= 0 && idut2 <= 0) {
	*imist = *c__;
    } else if (idut2 <= 0 && idut3 <= 0) {
	*imist = *d__;
    } else if (idut3 <= 0 && idut1 <= 0) {
	*imist = *b;
    } else if (idut1 <= 0) {
	if (*d__ > 8 && *ilift > 8) {
	    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &
		    y2[1], &z2[1], k, ileft, d__, ilift, mhalf, mfull, &isclp[
		    1], epz, &idot3);
	} else {
	    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], ileft, d__, ilift, k, mhalf, mfull, &isclp[
		    1], &idot3);
	}
	if (idot3 > 0) {
	    *imist = *b;
	} else if (idot3 < 0) {
	    *imist = *c__;
	} else {
	    *itype = -4;
	    *imist = *d__;
	}
    } else if (idut2 <= 0) {
	if (*b > 8 && *ilift > 8) {
	    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &
		    y2[1], &z2[1], k, ileft, b, ilift, mhalf, mfull, &isclp[1]
		    , epz, &idot1);
	} else {
	    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], ileft, b, ilift, k, mhalf, mfull, &isclp[1]
		    , &idot1);
	}
	if (idot1 > 0) {
	    *imist = *c__;
	} else if (idot1 < 0) {
	    *imist = *d__;
	} else {
	    *itype = -4;
	    *imist = *b;
	}
    } else {
	if (*c__ > 8 && *ilift > 8) {
	    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &
		    y2[1], &z2[1], k, ileft, c__, ilift, mhalf, mfull, &isclp[
		    1], epz, &idot2);
	} else {
	    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], ileft, c__, ilift, k, mhalf, mfull, &isclp[
		    1], &idot2);
	}
	if (idot2 > 0) {
	    *imist = *d__;
	} else if (idot2 < 0) {
	    *imist = *b;
	} else {
	    *itype = -4;
	    *imist = *c__;
	}
    }

L1000:
    return 0;
} /* fcfind_ */

/* FCEDGE */

/*     This subroutine will test whether a ray through an edge of a */
/*     tetrahedron intersects either of the facets of the tetrahedron */
/*     opposite the edge and whether a point in the interior of the */
/*     ray is contained in the tetrahedron */

/* Subroutine */ int fcedge_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *xc, integer *yc, integer *zc, 
	integer *itype, integer *ileft, integer *k, integer *icon, integer *
	iscur, integer *imist, integer *ivnxt, integer *site0, integer *site1,
	 integer *side1, integer *side2, integer *mhalf, integer *mfull, 
	integer *isclp, integer *itchk)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer idut, idot0, site2, site3, iside[4], isnow, ipout, iasign;
    extern /* Subroutine */ int ipsign_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), artsig_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), reordr_(
	    integer *, integer *, integer *, integer *), pntype_(integer *, 
	    integer *, integer *, integer *);



/*     find intersecting facet */

    /* Parameter adjustments */
    --isclp;
    icon -= 9;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    reordr_(&icon[9], site0, site1, iscur);
    site2 = icon[(*iscur << 3) + 7];
    *site0 = icon[(*iscur << 3) + 8];
    isnow = (i__1 = icon[(*iscur << 3) + 1], abs(i__1));

L300:
    ++(*itchk);
    if (isnow <= 0 || isnow > *ivnxt) {
	s_stop("961", (ftnlen)3);
    }
    if (isnow == *iscur) {
	s_stop("963", (ftnlen)3);
    }
    reordr_(&icon[9], site0, site1, &isnow);
    site3 = icon[(isnow << 3) + 8];
    if (*site1 > 8 && site2 > 8 && site3 > 8) {
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], site1, &site3, 
		&site2, k, mhalf, mfull, &isclp[1], &idut);
    } else {
	artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], site1, &site3, &site2, k, mhalf, mfull, &isclp[1], &
		idut);
    }
    if (idut >= 0) {
	goto L400;
    }
    *site0 = site3;
    isnow = (i__1 = icon[(isnow << 3) + 1], abs(i__1));
    goto L300;

L400:
    *iscur = isnow;

/*     determine whether point k is in tetrahedron */

    if (*site0 <= 8 || *site1 <= 8 || site2 <= 8 || site3 <= 8) {
	goto L500;
    }

    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], site1, site0, &
	    site3, k, mhalf, mfull, &isclp[1], &ipout);
    iside[2] = ipout;
    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &site2, &site3, 
	    site0, k, mhalf, mfull, &isclp[1], &ipout);
    iside[1] = ipout;

    if (iside[1] < 0 || iside[2] < 0) {
	goto L600;
    }

L450:
    iside[0] = idut;
    iside[3] = 1;
    pntype_(iside, itype, side1, side2);
    goto L1000;

/*     there is at least one artificial vertex */

L500:
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], site1, site0, &site3, k, mhalf, mfull, &isclp[1], &iasign);
    iside[2] = iasign;
    artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1], &
	    zc[1], &site2, &site3, site0, k, mhalf, mfull, &isclp[1], &iasign)
	    ;
    iside[1] = iasign;

    if (iside[1] < 0 || iside[2] < 0) {
	goto L600;
    }
    goto L450;

/*     k is not in tetrahedron but ray intersects one of the facets */
/*     of the tetrahedron opposite the edge */

L600:
    if (*site0 > 8 && site3 > 8) {
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], ileft, site0, &
		site3, k, mhalf, mfull, &isclp[1], &idot0);
    } else {
	artsig_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], ileft, site0, &site3, k, mhalf, mfull, &isclp[1], &
		idot0);
    }
    if (idot0 > 0) {
	if (idut > 0) {
	    *itype = -2;
	    *imist = *site1;
	} else {
	    *itype = -4;
	    *imist = site2;
	}
    } else if (idot0 < 0) {
	if (idut > 0) {
	    *itype = -2;
	    *imist = site2;
	} else {
	    *itype = -4;
	    *imist = *site1;
	}
    } else {
	if (idut > 0) {
	    *itype = -4;
	    *imist = *site0;
	} else {
	    *itype = -3;
	    *imist = site3;
	}
    }

L1000:
    return 0;
} /* fcedge_ */

/* PNTYPE */

/*     This subroutine determines point type with respect to a */
/*     tetrahedron that contains a point */

/* Subroutine */ int pntype_(integer *iside, integer *itype, integer *side1, 
	integer *side2)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);



/*     point is in the interior of tetrahedron */

    /* Parameter adjustments */
    --iside;

    /* Function Body */
    if (iside[1] > 0 && iside[2] > 0 && iside[3] > 0 && iside[4] > 0) {
	*itype = 2;
	goto L1000;
    }

/*     unacceptable situation */

    if (iside[1] == 0 && iside[2] == 0 && iside[3] == 0 && iside[4] == 0) {
	s_stop("971", (ftnlen)3);
    }

/*     point is a vertex of tetrahedron */

    if (iside[1] == 0 && iside[2] == 0 && iside[3] == 0) {
	*itype = 1;
	*side1 = 4;
	goto L1000;
    } else if (iside[1] == 0 && iside[2] == 0 && iside[4] == 0) {
	*itype = 1;
	*side1 = 3;
	goto L1000;
    } else if (iside[1] == 0 && iside[3] == 0 && iside[4] == 0) {
	*itype = 1;
	*side1 = 2;
	goto L1000;
    } else if (iside[2] == 0 && iside[3] == 0 && iside[4] == 0) {
	*itype = 1;
	*side1 = 1;
	goto L1000;
    }

/*     point is in the interior of an edge of tetrahedron */

    if (iside[1] == 0 && iside[2] == 0) {
	*itype = 3;
	*side1 = 1;
	*side2 = 2;
	goto L1000;
    } else if (iside[1] == 0 && iside[3] == 0) {
	*itype = 3;
	*side1 = 1;
	*side2 = 3;
	goto L1000;
    } else if (iside[1] == 0 && iside[4] == 0) {
	*itype = 3;
	*side1 = 1;
	*side2 = 4;
	goto L1000;
    } else if (iside[2] == 0 && iside[3] == 0) {
	*itype = 3;
	*side1 = 2;
	*side2 = 3;
	goto L1000;
    } else if (iside[2] == 0 && iside[4] == 0) {
	*itype = 3;
	*side1 = 2;
	*side2 = 4;
	goto L1000;
    } else if (iside[3] == 0 && iside[4] == 0) {
	*itype = 3;
	*side1 = 3;
	*side2 = 4;
	goto L1000;
    }

/*     point is in the interior of a facet of tetrahedron */

    *itype = 4;
    if (iside[1] == 0) {
	*side1 = 1;
    } else if (iside[2] == 0) {
	*side1 = 2;
    } else if (iside[3] == 0) {
	*side1 = 3;
    } else if (iside[4] == 0) {
	*side1 = 4;
    } else {
	s_stop("972", (ftnlen)3);
    }

L1000:
    return 0;
} /* pntype_ */

/* ARTSIG */

/*     This subroutine determines the position of a point with respect */
/*     to a plane when artificial points may be involved */

/* Subroutine */ int artsig_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *xc, integer *yc, integer *zc, 
	integer *ib, integer *id, integer *ic, integer *k, integer *mhalf, 
	integer *mfull, integer *isclp, integer *iasign)
{
    integer b, c__, d__, ifn;
    extern /* Subroutine */ int ipsig3_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), ipsig4_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), ipsig6_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), ipsign_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     vrtarr_(integer *, integer *, integer *, integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    vrtarr_(ib, id, ic, &b, &d__, &c__);

    if (d__ > 8 && c__ > 8) {
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &b, &d__, &c__,
		 k, mhalf, mfull, &isclp[1], iasign);
	goto L100;
    } else if (b <= 8) {
	*iasign = 1;
	goto L100;
    } else if (d__ <= 8 && c__ <= 8) {
	ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], &d__, &c__, k, &b, mhalf, mfull, &isclp[1], iasign);
	if (*iasign != 0) {
	    goto L100;
	}
	ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], &d__, &c__, k, &b, mhalf, mfull, &isclp[1], iasign);
	goto L100;
    } else if (d__ <= 8) {
	ifn = 0;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], &b, &c__, k, &d__, mhalf, mfull, &isclp[1], &ifn, 
		iasign);
    } else {
	ifn = 1;
	ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &yc[1],
		 &zc[1], &b, &d__, k, &c__, mhalf, mfull, &isclp[1], &ifn, 
		iasign);
    }
    if (*iasign != 0) {
	goto L100;
    }
    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &b, &d__, &c__, k, 
	    mhalf, mfull, &isclp[1], iasign);

L100:
    return 0;
} /* artsig_ */

/* VRTINS */

/* This subroutine will insert a point located at a vertex of the current */
/* tetrahedron */

/* Subroutine */ int vrtins_(integer *k, integer *w, integer *w2, integer *
	icon, integer *nvmax, integer *tetra, integer *curr, integer *is, 
	integer *id, integer *iftal, integer *side, integer *ifl, integer *
	newtts, integer *ired, logical *delaun, logical *flphis, integer *
	mhalf, integer *mfull, integer *isclw)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, indx, site1, site2, isadj, itide, ilift, isini, iscur, 
	    isfrt, islst;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *), iwsign_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), sitord_(integer *, 
	    integer *, integer *);



    /* Parameter adjustments */
    --isclw;
    --ifl;
    --id;
    --is;
    icon -= 9;
    --w2;
    --w;

    /* Function Body */
    site1 = icon[*side + 4 + (*curr << 3)];
    if (site1 <= 8) {
	s_stop("1010", (ftnlen)4);
    }

    if (*delaun) {
	goto L50;
    }
    iwsign_(&w[1], &w2[1], k, &site1, mhalf, mfull, &isclw[1], &itide);
    if (itide > 0) {
	goto L100;
    }
L50:
    is[*k] = -site1;
    *ired = 1;
    goto L2000;

L100:
    if (! (*flphis)) {
	goto L1000;
    }
    *ired = 0;
    is[site1] = -2;
    sitord_(&icon[9], &site1, curr);
    isini = *curr;
    ++(*tetra);
    if (*tetra > *nvmax) {
	s_stop("1020", (ftnlen)4);
    }
    isfrt = *tetra;
    icon[(isini << 3) + 5] = -(*tetra);
    icon[(*tetra << 3) + 5] = isini;
    islst = isini;
    iscur = icon[(islst << 3) + 2];
    if (iscur <= 0) {
	s_stop("1030", (ftnlen)4);
    }
    site2 = icon[(islst << 3) + 8];
    is[*k] = *tetra;
    is[site2] = *tetra;
    is[icon[(islst << 3) + 6]] = *tetra;
    is[icon[(islst << 3) + 7]] = *tetra;
    *newtts = 1;
    ifl[1] = *tetra;

L200:
    reordr_(&icon[9], &site1, &site2, &iscur);
    if (icon[(iscur << 3) + 4] != islst) {
	s_stop("1040", (ftnlen)4);
    }
    ++(*tetra);
    if (*tetra > *nvmax) {
	s_stop("1050", (ftnlen)4);
    }
    icon[(iscur << 3) + 5] = -(*tetra);
    icon[(*tetra << 3) + 5] = iscur;
    islst = iscur;
    indx = 2;
    iscur = icon[(islst << 3) + 2];
    if (iscur <= 0) {
	s_stop("1060", (ftnlen)4);
    }
    site2 = icon[(islst << 3) + 8];
    is[site2] = *tetra;
    ++(*newtts);
    ifl[*newtts] = *tetra;
    if (icon[(iscur << 3) + 5] > 0) {
	goto L200;
    }

L500:
    if (indx == 2) {
	indx = 3;
	iscur = icon[(islst << 3) + 3];
	if (iscur <= 0) {
	    s_stop("1070", (ftnlen)4);
	}
	site2 = icon[(islst << 3) + 6];
	if (icon[(iscur << 3) + 5] > 0) {
	    goto L200;
	}
	goto L500;
    } else if (indx == 3) {
	if (islst != isini) {
	    iscur = islst;
	    islst = icon[(iscur << 3) + 4];
	    if (icon[(islst << 3) + 2] == iscur) {
		indx = 2;
	    } else if (icon[(islst << 3) + 3] == iscur) {
		indx = 3;
	    } else if (icon[(islst << 3) + 4] == iscur) {
		indx = 4;
	    } else {
		s_stop("1080", (ftnlen)4);
	    }
	    goto L500;
	} else {
	    indx = 4;
	    iscur = icon[(islst << 3) + 4];
	    if (iscur <= 0) {
		s_stop("1090", (ftnlen)4);
	    }
	    site2 = icon[(islst << 3) + 7];
	    if (icon[(iscur << 3) + 5] > 0) {
		goto L200;
	    }
	    goto L500;
	}
    }
    if (islst != isini) {
	s_stop("1110", (ftnlen)4);
    }

    i__1 = *tetra;
    for (i__ = isfrt; i__ <= i__1; ++i__) {
	iscur = icon[(i__ << 3) + 5];
	for (j = 2; j <= 4; ++j) {
	    icon[j + (i__ << 3)] = -icon[(icon[j + (iscur << 3)] << 3) + 5];
/* L600: */
	}
	for (j = 6; j <= 8; ++j) {
	    icon[j + (i__ << 3)] = icon[j + (iscur << 3)];
/* L700: */
	}
/* L800: */
    }

    i__1 = *tetra;
    for (i__ = isfrt; i__ <= i__1; ++i__) {
	iscur = icon[(i__ << 3) + 5];
	icon[(i__ << 3) + 5] = *k;
	icon[(iscur << 3) + 5] = -site1;
	isadj = icon[(iscur << 3) + 1];
	icon[(iscur << 3) + 1] = i__;
	icon[(iscur << 3) + 2] = 0;
	icon[(iscur << 3) + 3] = 0;
	icon[(iscur << 3) + 4] = 0;
	icon[(i__ << 3) + 1] = isadj;
	if (isadj == 0) {
	    goto L900;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (isadj << 3)] == iscur) {
		goto L860;
	    }
/* L840: */
	}
	s_stop("1120", (ftnlen)4);
L860:
	icon[j + (isadj << 3)] = i__;
L900:
	;
    }
    goto L2000;

L1000:
    *ired = 0;
    ++(*iftal);
    iscur = is[site1];
    if (iscur <= 0 || iscur > *tetra) {
	s_stop("1130", (ftnlen)4);
    }
    isini = iscur;
    sitord_(&icon[9], &site1, &isini);
    is[site1] = -2;
    is[*k] = isini;
    *newtts = 0;

L1400:
    ++(*newtts);
    ifl[*newtts] = iscur;
    id[iscur] = *iftal;
    icon[(iscur << 3) + 5] = *k;

    isadj = icon[(iscur << 3) + 2];
    if (isadj <= 0 || isadj > *tetra) {
	s_stop("1140", (ftnlen)4);
    }
    if (id[isadj] == *iftal) {
	goto L1600;
    }
    ilift = icon[(iscur << 3) + 8];
    goto L1900;
L1600:
    isadj = icon[(iscur << 3) + 3];
    if (isadj <= 0 || isadj > *tetra) {
	s_stop("1150", (ftnlen)4);
    }
    if (id[isadj] == *iftal) {
	goto L1700;
    }
    ilift = icon[(iscur << 3) + 6];
    goto L1900;
L1700:
    isadj = icon[(iscur << 3) + 4];
    if (isadj <= 0 || isadj > *tetra) {
	s_stop("1160", (ftnlen)4);
    }
    if (iscur == isini) {
	goto L1800;
    }
    if (icon[(isadj << 3) + 3] == iscur) {
	iscur = isadj;
	goto L1700;
    } else if (icon[(isadj << 3) + 2] == iscur) {
	iscur = isadj;
	goto L1600;
    } else if (icon[(isadj << 3) + 4] == iscur) {
	if (isadj != isini) {
	    s_stop("1170", (ftnlen)4);
	}
	goto L2000;
    } else {
	s_stop("1180", (ftnlen)4);
    }
L1800:
    if (id[isadj] == *iftal) {
	goto L2000;
    }
    ilift = icon[(iscur << 3) + 7];

L1900:
    reordr_(&icon[9], &site1, &ilift, &isadj);
    iscur = isadj;
    goto L1400;

L2000:
    return 0;
} /* vrtins_ */

/* INTINS */

/* This subroutine will insert a point located in the interior of the */
/* current tetrahedron.  Four new tetra will be created. */

/* Subroutine */ int intins_(integer *k, doublereal *xi, doublereal *yi, 
	doublereal *zi, doublereal *wi, integer *x, integer *y, integer *z__, 
	integer *w, integer *x2, integer *y2, integer *z2, integer *w2, 
	integer *icon, integer *ih, integer *ihn, integer *nhmax, integer *
	nvmax, integer *tetra, integer *curr, integer *is, integer *ifl, 
	integer *newtts, integer *ired, logical *delaun, logical *flphis, 
	integer *mhalf, integer *mfull, integer *isclp, integer *isclw, 
	integer *isclr, doublereal *epz)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer a, b, c__, d__, i__, j, adj, new1, new2, new3, new4, newi[4];
    doublereal xctr, yctr, zctr;
    extern /* Subroutine */ int ctrad_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, logical 
	    *, integer *);
    integer itide;
    doublereal tdist;
    extern /* Subroutine */ int bisphr_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *), iqsign_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *), sitord_(integer *, 
	    integer *, integer *);
    integer ipossi;



    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --ifl;
    --is;
    --ih;
    icon -= 9;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;
    --wi;
    --zi;
    --yi;
    --xi;

    /* Function Body */
    a = icon[(*curr << 3) + 5];
    b = icon[(*curr << 3) + 6];
    c__ = icon[(*curr << 3) + 7];
    d__ = icon[(*curr << 3) + 8];

    if (*delaun) {
	goto L30;
    }
    if (a <= 8 || b <= 8 || c__ <= 8 || d__ <= 8) {
	goto L30;
    }
    ctrad_(&xi[1], &yi[1], &zi[1], &wi[1], &xctr, &yctr, &zctr, &a, &b, &c__, 
	    &d__, epz, delaun, &ipossi);
    if (ipossi == 1) {
	goto L20;
    }
    bisphr_(&xi[1], &yi[1], &zi[1], &wi[1], k, &a, epz, &xctr, &yctr, &zctr, &
	    tdist, delaun, &ipossi);
    if (ipossi == 1) {
	goto L20;
    }
    itide = 1;
    if (tdist > 0.) {
	itide = -1;
    }
    goto L25;
L20:
    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], &a, 
	    &b, &c__, &d__, k, mhalf, mfull, &isclp[1], &isclw[1], &isclr[1], 
	    delaun, &itide);
L25:
    if (itide <= 0) {
	goto L30;
    }
    is[*k] = -3;
    *ired = 1;
    goto L1000;

L30:
    *ired = 0;

    if (*flphis) {
	new1 = *tetra + 1;
	new2 = *tetra + 2;
	new3 = *tetra + 3;
	new4 = *tetra + 4;
	*tetra = new4;
	if (*tetra > *nvmax) {
	    s_stop("1210", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    *tetra = new1;
	    if (*tetra > *nvmax) {
		s_stop("1220", (ftnlen)4);
	    }
	} else {
	    new1 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new2 = *tetra + 1;
	    *tetra = new2;
	    if (*tetra > *nvmax) {
		s_stop("1230", (ftnlen)4);
	    }
	} else {
	    new2 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new3 = *tetra + 1;
	    *tetra = new3;
	    if (*tetra > *nvmax) {
		s_stop("1240", (ftnlen)4);
	    }
	} else {
	    new3 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new4 = *tetra + 1;
	    *tetra = new4;
	    if (*tetra > *nvmax) {
		s_stop("1250", (ftnlen)4);
	    }
	} else {
	    new4 = ih[*ihn];
	    --(*ihn);
	}
    }

    ifl[1] = new1;
    ifl[2] = new2;
    ifl[3] = new3;
    ifl[4] = new4;

    *newtts = 4;

    for (i__ = 1; i__ <= 8; ++i__) {
	icon[i__ + (new1 << 3)] = icon[i__ + (*curr << 3)];
	icon[i__ + (new2 << 3)] = icon[i__ + (*curr << 3)];
	icon[i__ + (new3 << 3)] = icon[i__ + (*curr << 3)];
	icon[i__ + (new4 << 3)] = icon[i__ + (*curr << 3)];
/* L50: */
    }

/*     update new1 */

    icon[(new1 << 3) + 2] = new2;
    icon[(new1 << 3) + 3] = new3;
    icon[(new1 << 3) + 4] = new4;
    icon[(new1 << 3) + 5] = *k;

/*     update new2 */

    icon[(new2 << 3) + 1] = new1;
    icon[(new2 << 3) + 3] = new3;
    icon[(new2 << 3) + 4] = new4;
    icon[(new2 << 3) + 6] = *k;
    sitord_(&icon[9], k, &new2);

/*     update new3 */

    icon[(new3 << 3) + 1] = new1;
    icon[(new3 << 3) + 2] = new2;
    icon[(new3 << 3) + 4] = new4;
    icon[(new3 << 3) + 7] = *k;
    sitord_(&icon[9], k, &new3);

/*     update new4 */

    icon[(new4 << 3) + 1] = new1;
    icon[(new4 << 3) + 2] = new2;
    icon[(new4 << 3) + 3] = new3;
    icon[(new4 << 3) + 8] = *k;
    sitord_(&icon[9], k, &new4);

/*     update is(*) */

    is[*k] = new1;

    is[a] = new2;
    is[b] = new1;
    is[c__] = new1;
    is[d__] = new1;

/*     update neighboring tetra */

    newi[0] = new1;
    newi[1] = new2;
    newi[2] = new3;
    newi[3] = new4;
    for (i__ = 1; i__ <= 4; ++i__) {
	adj = icon[i__ + (*curr << 3)];
	if (adj != 0) {
	    for (j = 1; j <= 4; ++j) {
		if (icon[j + (adj << 3)] == *curr) {
		    icon[j + (adj << 3)] = newi[i__ - 1];
		    goto L200;
		}
/* L175: */
	    }
	    s_stop("1260", (ftnlen)4);
	}
L200:
	;
    }

/*     flag current tetra to denote that it has children */

    if (*flphis) {
	icon[(*curr << 3) + 5] = -icon[(*curr << 3) + 5];
	icon[(*curr << 3) + 1] = new1;
	icon[(*curr << 3) + 2] = new2;
	icon[(*curr << 3) + 3] = new3;
	icon[(*curr << 3) + 4] = new4;
    } else {
	icon[(*curr << 3) + 5] = -icon[(*curr << 3) + 5];
	++(*ihn);
	if (*ihn > *nhmax) {
	    s_stop("1270", (ftnlen)4);
	}
	ih[*ihn] = *curr;
    }

L1000:
    return 0;
} /* intins_ */

/* SIDINS */

/* This subroutine will insert a point which is on side */
/* of curr (and adj). Six new tetra will be created. */

/* Subroutine */ int sidins_(integer *k, integer *x, integer *y, integer *z__,
	 integer *w, integer *x2, integer *y2, integer *z2, integer *w2, 
	integer *icon, integer *ih, integer *ihn, integer *nhmax, integer *
	nvmax, integer *tetra, integer *curr, integer *is, integer *side, 
	integer *ifl, integer *newtts, integer *ired, logical *delaun, 
	logical *flphis, integer *mhalf, integer *mfull, integer *isclp, 
	integer *isclw, integer *isclr)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, adj, new1, new2, new3, new4, new5, new6, newi[6], temp, 
	    site1, site2, site3, site4, itide;
    extern /* Subroutine */ int iqsig2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *), reordr_(integer *, 
	    integer *, integer *, integer *), sitord_(integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --ifl;
    --is;
    --ih;
    icon -= 9;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    adj = icon[*side + (*curr << 3)];
    if (adj == 0) {
	s_stop("1310", (ftnlen)4);
    }

/*     rearrange curr */

    site1 = icon[*side + 4 + (*curr << 3)];
    sitord_(&icon[9], &site1, curr);

    site2 = icon[(*curr << 3) + 6];
    site3 = icon[(*curr << 3) + 7];
    site4 = icon[(*curr << 3) + 8];

    if (*delaun) {
	goto L30;
    }
    if (site2 <= 8 || site3 <= 8 || site4 <= 8) {
	goto L30;
    }
    iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], &
	    site2, &site3, &site4, k, mhalf, mfull, &isclp[1], &isclw[1], &
	    isclr[1], delaun, &itide);
    if (itide <= 0) {
	goto L30;
    }
    is[*k] = -3;
    *ired = 1;
    goto L1000;

L30:
    *ired = 0;

    if (*flphis) {
	new1 = *tetra + 1;
	new2 = *tetra + 2;
	new3 = *tetra + 3;
	new4 = *tetra + 4;
	new5 = *tetra + 5;
	new6 = *tetra + 6;
	*tetra = new6;
	if (*tetra > *nvmax) {
	    s_stop("1320", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    *tetra = new1;
	    if (*tetra > *nvmax) {
		s_stop("1330", (ftnlen)4);
	    }
	} else {
	    new1 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new2 = *tetra + 1;
	    *tetra = new2;
	    if (*tetra > *nvmax) {
		s_stop("1340", (ftnlen)4);
	    }
	} else {
	    new2 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new3 = *tetra + 1;
	    *tetra = new3;
	    if (*tetra > *nvmax) {
		s_stop("1350", (ftnlen)4);
	    }
	} else {
	    new3 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new4 = *tetra + 1;
	    *tetra = new4;
	    if (*tetra > *nvmax) {
		s_stop("1360", (ftnlen)4);
	    }
	} else {
	    new4 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new5 = *tetra + 1;
	    *tetra = new5;
	    if (*tetra > *nvmax) {
		s_stop("1370", (ftnlen)4);
	    }
	} else {
	    new5 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new6 = *tetra + 1;
	    *tetra = new6;
	    if (*tetra > *nvmax) {
		s_stop("1380", (ftnlen)4);
	    }
	} else {
	    new6 = ih[*ihn];
	    --(*ihn);
	}
    }

    ifl[1] = new1;
    ifl[2] = new2;
    ifl[3] = new3;
    ifl[4] = new4;
    ifl[5] = new5;
    ifl[6] = new6;

    *newtts = 6;

/*     create new1, new2, new3 */

    for (i__ = 1; i__ <= 8; ++i__) {
	icon[i__ + (new1 << 3)] = icon[i__ + (*curr << 3)];
	icon[i__ + (new2 << 3)] = icon[i__ + (*curr << 3)];
	icon[i__ + (new3 << 3)] = icon[i__ + (*curr << 3)];
/* L50: */
    }

    icon[(new1 << 3) + 1] = new4;
    icon[(new1 << 3) + 3] = new2;
    icon[(new1 << 3) + 4] = new3;
    icon[(new1 << 3) + 6] = *k;
    sitord_(&icon[9], k, &new1);

    icon[(new2 << 3) + 1] = new6;
    icon[(new2 << 3) + 2] = new1;
    icon[(new2 << 3) + 4] = new3;
    icon[(new2 << 3) + 7] = *k;
    sitord_(&icon[9], k, &new2);

    icon[(new3 << 3) + 1] = new5;
    icon[(new3 << 3) + 2] = new1;
    icon[(new3 << 3) + 3] = new2;
    icon[(new3 << 3) + 8] = *k;
    sitord_(&icon[9], k, &new3);

/*     update is(*) */

    is[*k] = new1;

    is[site1] = new1;
    is[site2] = new2;
    is[site3] = new1;
    is[site4] = new1;

/*     update curr's neighbors */

    newi[0] = new1;
    newi[1] = new2;
    newi[2] = new3;
    for (i__ = 2; i__ <= 4; ++i__) {
	temp = icon[i__ + (*curr << 3)];
	if (temp != 0) {
	    for (j = 1; j <= 4; ++j) {
		if (icon[j + (temp << 3)] == *curr) {
		    icon[j + (temp << 3)] = newi[i__ - 2];
		    goto L200;
		}
/* L150: */
	    }
	    s_stop("1390", (ftnlen)4);
	}
L200:
	;
    }

/*     flag curr to show its children */

    if (*flphis) {
	icon[(*curr << 3) + 5] = -icon[(*curr << 3) + 5];
	icon[(*curr << 3) + 1] = new1;
	icon[(*curr << 3) + 2] = new2;
	icon[(*curr << 3) + 3] = new3;
	icon[(*curr << 3) + 4] = -adj;
    } else {
	icon[(*curr << 3) + 5] = -icon[(*curr << 3) + 5];
	++(*ihn);
	if (*ihn > *nhmax) {
	    s_stop("1410", (ftnlen)4);
	}
	ih[*ihn] = *curr;
    }

/*     update 2nd tet (adj) */

    for (i__ = 1; i__ <= 4; ++i__) {
	if (icon[i__ + (adj << 3)] == *curr) {
	    site1 = icon[i__ + 4 + (adj << 3)];
	    goto L325;
	}
/* L300: */
    }
    s_stop("1420", (ftnlen)4);
L325:

    reordr_(&icon[9], &site1, &site2, &adj);

/*     create new4, new5, new6 */

    for (i__ = 1; i__ <= 8; ++i__) {
	icon[i__ + (new4 << 3)] = icon[i__ + (adj << 3)];
	icon[i__ + (new5 << 3)] = icon[i__ + (adj << 3)];
	icon[i__ + (new6 << 3)] = icon[i__ + (adj << 3)];
/* L350: */
    }

    icon[(new4 << 3) + 1] = new1;
    icon[(new4 << 3) + 3] = new5;
    icon[(new4 << 3) + 4] = new6;
    icon[(new4 << 3) + 6] = *k;
    sitord_(&icon[9], k, &new4);

    icon[(new5 << 3) + 1] = new3;
    icon[(new5 << 3) + 2] = new4;
    icon[(new5 << 3) + 4] = new6;
    icon[(new5 << 3) + 7] = *k;
    sitord_(&icon[9], k, &new5);

    icon[(new6 << 3) + 1] = new2;
    icon[(new6 << 3) + 2] = new4;
    icon[(new6 << 3) + 3] = new5;
    icon[(new6 << 3) + 8] = *k;
    sitord_(&icon[9], k, &new6);

/*     update is(*) */

    is[site1] = new4;

/*     update adj's neighbors */

    newi[3] = new4;
    newi[4] = new5;
    newi[5] = new6;
    for (i__ = 2; i__ <= 4; ++i__) {
	temp = icon[i__ + (adj << 3)];
	if (temp != 0) {
	    for (j = 1; j <= 4; ++j) {
		if (icon[j + (temp << 3)] == adj) {
		    icon[j + (temp << 3)] = newi[i__ + 1];
		    goto L500;
		}
/* L450: */
	    }
	    s_stop("1430", (ftnlen)4);
	}
L500:
	;
    }

/*     flag adj to show its children */

    if (*flphis) {
	icon[(adj << 3) + 5] = -icon[(adj << 3) + 5];
	icon[(adj << 3) + 1] = new4;
	icon[(adj << 3) + 2] = new5;
	icon[(adj << 3) + 3] = new6;
	icon[(adj << 3) + 4] = -(*curr);
    } else {
	icon[(adj << 3) + 5] = -icon[(adj << 3) + 5];
	++(*ihn);
	if (*ihn > *nhmax) {
	    s_stop("1440", (ftnlen)4);
	}
	ih[*ihn] = adj;
    }

L1000:
    return 0;
} /* sidins_ */

/* EDGINS */

/* This subroutine will insert point on edge of curr tetra. */

/* Subroutine */ int edgins_(integer *k, integer *x, integer *y, integer *z__,
	 integer *w, integer *x2, integer *y2, integer *z2, integer *w2, 
	integer *icon, integer *ih, integer *ihn, integer *nhmax, integer *
	nvmax, integer *tetra, integer *curr, integer *is, integer *side1, 
	integer *side2, integer *ifl, integer *newtts, integer *ired, logical 
	*delaun, logical *flphis, integer *mhalf, integer *mfull, integer *
	isclp, integer *isclw, integer *isclr)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, adj, nel1, nel2, new1, new2, now1, now2, site1, site2, 
	    itide;
    extern /* Subroutine */ int iqsig1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, integer *);
    integer fndsit, firtet;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *), sitord_(integer *, integer *, integer *);
    integer prvtet;



/*     find endpoints of edge */

    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --ifl;
    --is;
    --ih;
    icon -= 9;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    fndsit = 0;
    for (i__ = 5; i__ <= 8; ++i__) {
	if (fndsit == 0) {
	    if (i__ == *side1 + 4 || i__ == *side2 + 4) {
		goto L100;
	    } else {
		site1 = icon[i__ + (*curr << 3)];
		fndsit = 1;
	    }
	} else {
	    if (i__ == *side1 + 4 || i__ == *side2 + 4) {
		goto L100;
	    } else {
		site2 = icon[i__ + (*curr << 3)];
		goto L150;
	    }
	}
L100:
	;
    }
    s_stop("1510", (ftnlen)4);
L150:

    if (*delaun) {
	goto L160;
    }
    if (site1 <= 8 || site2 <= 8) {
	goto L160;
    }
    iqsig1_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], &
	    site1, &site2, k, mhalf, mfull, &isclp[1], &isclw[1], &isclr[1], 
	    delaun, &itide);
    if (itide <= 0) {
	goto L160;
    }
    is[*k] = -3;
    *ired = 1;
    goto L1000;

/*     order vertices of tetrahedra around edge */

L160:
    firtet = *curr;

L163:
    reordr_(&icon[9], &site1, &site2, curr);
    *curr = icon[(*curr << 3) + 3];
    if (*curr != firtet) {
	goto L163;
    }

    *ired = 0;
    if (*flphis) {
	new1 = *tetra + 1;
	new2 = *tetra + 2;
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    new2 = *tetra + 2;
	} else if (*ihn == 1) {
	    new1 = ih[1];
	    new2 = *tetra + 1;
	} else {
	    new1 = ih[*ihn];
	    new2 = ih[*ihn - 1];
	}
    }
    nel1 = new1;
    nel2 = new2;
    *newtts = 0;

    is[*k] = new1;
    is[site1] = new2;
    is[site2] = new1;

/*     create 2 new tetra */

L175:
    if (*flphis) {
	now1 = *tetra + 1;
	now2 = *tetra + 2;
	*tetra = now2;
	if (*tetra > *nvmax) {
	    s_stop("1520", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    now1 = *tetra + 1;
	    *tetra = now1;
	    if (*tetra > *nvmax) {
		s_stop("1530", (ftnlen)4);
	    }
	} else {
	    now1 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    now2 = *tetra + 1;
	    *tetra = now2;
	    if (*tetra > *nvmax) {
		s_stop("1540", (ftnlen)4);
	    }
	} else {
	    now2 = ih[*ihn];
	    --(*ihn);
	}
    }

    ++(*newtts);
    ifl[*newtts] = now1;
    ++(*newtts);
    ifl[*newtts] = now2;

    for (i__ = 1; i__ <= 8; ++i__) {
	icon[i__ + (now1 << 3)] = icon[i__ + (*curr << 3)];
	icon[i__ + (now2 << 3)] = icon[i__ + (*curr << 3)];
/* L180: */
    }

    icon[(now1 << 3) + 2] = now2;
    icon[(now1 << 3) + 5] = *k;
    icon[(now2 << 3) + 1] = now1;
    icon[(now2 << 3) + 6] = *k;

    icon[(nel1 << 3) + 3] = now1;
    icon[(nel2 << 3) + 3] = now2;
    icon[(now1 << 3) + 4] = nel1;
    icon[(now2 << 3) + 4] = nel2;
    sitord_(&icon[9], k, &now1);
    sitord_(&icon[9], k, &now2);

/*     update is(*) */

    is[icon[(*curr << 3) + 7]] = now1;

/*     update neighbors of curr */

    for (i__ = 1; i__ <= 2; ++i__) {
	adj = icon[i__ + (*curr << 3)];
	if (adj != 0) {
	    for (j = 1; j <= 4; ++j) {
		if (icon[j + (adj << 3)] == *curr) {
		    if (i__ == 1) {
			icon[j + (adj << 3)] = now1;
		    } else {
			icon[j + (adj << 3)] = now2;
		    }
		    goto L300;
		}
/* L250: */
	    }
	    s_stop("1550", (ftnlen)4);
	}
L300:
	;
    }

    prvtet = *curr;
    *curr = icon[(*curr << 3) + 3];

/*     show children of old tetra */

    if (*flphis) {
	icon[(prvtet << 3) + 5] = -icon[(prvtet << 3) + 5];
	icon[(prvtet << 3) + 1] = now1;
	icon[(prvtet << 3) + 2] = now2;
	icon[(prvtet << 3) + 3] = -(*curr);
	icon[(prvtet << 3) + 4] = -icon[(prvtet << 3) + 4];
    } else {
	icon[(prvtet << 3) + 5] = -icon[(prvtet << 3) + 5];
	++(*ihn);
	if (*ihn > *nhmax) {
	    s_stop("1560", (ftnlen)4);
	}
	ih[*ihn] = prvtet;
    }

/*     go to next tetrahedron until we're back at firtet */

    if (*curr != firtet) {
	nel1 = now1;
	nel2 = now2;
	goto L175;
    } else {
	icon[(new1 << 3) + 4] = now1;
	icon[(new2 << 3) + 2] = now2;
	icon[(now1 << 3) + 3] = new1;
	icon[(now2 << 3) + 3] = new2;
    }

L1000:
    return 0;
} /* edgins_ */

/* SPHERE */

/* This subroutine will optimize locally at point k for the */
/* Regular/Delaunay property */

/* Subroutine */ int sphere_(integer *k, integer *icon, integer *ifl, integer 
	*ih, integer *ihn, integer *nhmax, integer *newtts, doublereal *xi, 
	doublereal *yi, doublereal *zi, doublereal *wi, integer *x, integer *
	y, integer *z__, integer *w, integer *x2, integer *y2, integer *z2, 
	integer *w2, integer *tetra, integer *is, integer *nvmax, integer *xc,
	 integer *yc, integer *zc, logical *delaun, logical *flphis, integer *
	mhalf, integer *mfull, integer *isclp, integer *isclw, integer *isclr,
	 doublereal *epz)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();

    /* Local variables */
    integer a, b, c__, d__, i__, j, ib, ic, id, adj, ifn, now;
    doublereal xctr, yctr, zctr;
    integer istt, iside[4], itide;
    extern /* Subroutine */ int ctrad_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, logical 
	    *, integer *), flip23_(integer *, integer *, integer *, integer *,
	     integer *, integer *, integer *, integer *, logical *, integer *,
	     integer *, integer *, integer *), flip32_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    logical *, integer *, integer *, integer *, integer *, integer *);
    integer isodd;
    extern /* Subroutine */ int flip22_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, logical *, 
	    integer *, integer *, integer *, integer *, integer *), flip41_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *, integer *, integer *, 
	    integer *, integer *), flip21_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, logical *, 
	    integer *, integer *, integer *, integer *, integer *), flip31_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *, integer *, integer *, 
	    integer *, integer *);
    integer isite;
    doublereal tdist;
    integer ipout;
    extern /* Subroutine */ int ipsig1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), ipsig2_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     ipsig3_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), ipsig4_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), iqsig2_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, logical *, integer *), ipsig6_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), iqsig1_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *);
    integer oddsid;
    extern /* Subroutine */ int bisphr_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *), ipsign_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), iqsign_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, logical *, 
	    integer *), irsign_(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *), reordr_(integer *, integer *,
	     integer *, integer *);
    integer sidist[4], ipossi, opvert;
    extern /* Subroutine */ int vrtarr_(integer *, integer *, integer *, 
	    integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___317 = { 0, 6, 0, 0, 0 };
    static cilist io___318 = { 0, 6, 0, 0, 0 };




    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --zc;
    --yc;
    --xc;
    --is;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;
    --wi;
    --zi;
    --yi;
    --xi;
    --ih;
    --ifl;
    icon -= 9;

    /* Function Body */
    i__ = 1;
L100:
    if (i__ > *newtts) {
	goto L1500;
    }

    now = ifl[i__];
    if (icon[(now << 3) + 5] < 0) {
	goto L1000;
    }
    a = icon[(now << 3) + 5];
    if (a != *k) {
	s_stop("1610", (ftnlen)4);
    }

/*     look at adj tet */

    adj = icon[(now << 3) + 1];
    if (adj == 0) {
	goto L1000;
    }
    if (icon[(adj << 3) + 5] == *k) {
	s_stop("1620", (ftnlen)4);
    }

    b = icon[(now << 3) + 6];
    c__ = icon[(now << 3) + 7];
    d__ = icon[(now << 3) + 8];
    if (b < c__) {
	b = c__;
    }
    if (b < d__) {
	b = d__;
    }
    if (b <= 8) {
	s_stop("1630", (ftnlen)4);
    }
    reordr_(&icon[9], &a, &b, &now);
    c__ = icon[(now << 3) + 7];
    d__ = icon[(now << 3) + 8];

/*     reorder adj */

    reordr_(&icon[9], &b, &c__, &adj);

    if (icon[(adj << 3) + 7] != d__) {
	s_stop("1640", (ftnlen)4);
    }
    if (icon[(adj << 3) + 4] != now) {
	s_stop("1650", (ftnlen)4);
    }

    opvert = icon[(adj << 3) + 8];
    if (opvert <= 8 || c__ <= 8 || d__ <= 8) {
	goto L300;
    }

/*     test whether now and adj form a Regular configuration */

    ctrad_(&xi[1], &yi[1], &zi[1], &wi[1], &xctr, &yctr, &zctr, &b, &c__, &
	    d__, &opvert, epz, delaun, &ipossi);
    if (ipossi == 1) {
	goto L150;
    }
    bisphr_(&xi[1], &yi[1], &zi[1], &wi[1], &a, &b, epz, &xctr, &yctr, &zctr, 
	    &tdist, delaun, &ipossi);
    if (ipossi == 1) {
	goto L150;
    }
    if (tdist <= 0.) {
	goto L1000;
    }
    goto L170;
L150:
    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], &a, 
	    &b, &c__, &d__, &opvert, mhalf, mfull, &isclp[1], &isclw[1], &
	    isclr[1], delaun, &itide);
    if (itide >= 0) {
	goto L1000;
    }

/*     compute distances from opvert to facets of now */

L170:
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], &a, &c__, &opvert, &d__, mhalf, mfull, &isclp[1], epz, &
	    ipout);
    iside[1] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], &a, &b, &d__, &opvert, mhalf, mfull, &isclp[1], epz, &
	    ipout);
    iside[2] = ipout;
    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1], &
	    z2[1], &a, &b, &opvert, &c__, mhalf, mfull, &isclp[1], epz, &
	    ipout);
    iside[3] = ipout;

/*     set sidist array */

    for (j = 2; j <= 4; ++j) {
	if (iside[j - 1] > 0) {
	    sidist[j - 1] = 0;
	} else if (iside[j - 1] < 0) {
	    sidist[j - 1] = -1;
	} else {
	    sidist[j - 1] = 1;
	}
/* L200: */
    }

/*     flip according to type of flip */

    if (sidist[1] == 0 && sidist[2] == 0 && sidist[3] == 0) {
	flip23_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], nvmax);
	goto L1000;
    } else if (sidist[1] == -1 && sidist[2] == 0 && sidist[3] == 0) {
	oddsid = 2;
	flip32_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == -1 && sidist[3] == 0) {
	oddsid = 3;
	flip32_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == 0 && sidist[3] == -1) {
	oddsid = 4;
	flip32_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 1 && sidist[2] == 0 && sidist[3] == 0) {
	oddsid = 2;
	flip22_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == 1 && sidist[3] == 0) {
	oddsid = 3;
	flip22_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == 0 && sidist[3] == 1) {
	oddsid = 4;
	flip22_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (*delaun) {
	s_wsle(&io___317);
	do_lio(&c__9, &c__1, "Warning: Delaunay sphere violation", (ftnlen)34)
		;
	e_wsle();
	goto L1000;
    }

    if (sidist[1] == 0 && sidist[2] == -1 && sidist[3] == -1) {
	oddsid = 2;
	flip41_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == -1 && sidist[2] == 0 && sidist[3] == -1) {
	oddsid = 3;
	flip41_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == -1 && sidist[2] == -1 && sidist[3] == 0) {
	oddsid = 4;
	flip41_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 0 && sidist[2] == 1 && sidist[3] == 1) {
	oddsid = 2;
	flip21_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == 0 && sidist[3] == 1) {
	oddsid = 3;
	flip21_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == 1 && sidist[3] == 0) {
	oddsid = 4;
	flip21_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 0 && sidist[2] == -1 && sidist[3] == 1) {
	oddsid = 2;
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == 0 && sidist[3] == -1) {
	oddsid = 3;
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == -1 && sidist[2] == 1 && sidist[3] == 0) {
	oddsid = 4;
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 0 && sidist[2] == 1 && sidist[3] == -1) {
	oddsid = -2;
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == -1 && sidist[2] == 0 && sidist[3] == 1) {
	oddsid = -3;
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == -1 && sidist[3] == 0) {
	oddsid = -4;
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else {
	s_wsle(&io___318);
	do_lio(&c__9, &c__1, "Warning: Regular sphere violation", (ftnlen)33);
	e_wsle();
    }
    goto L1000;

L300:
    if (opvert <= 8 && c__ > 8 && d__ > 8) {
	goto L1000;
    }

/*     determine signs of distances from opvert to facets of */
/*     tetrahedron now */

    if (opvert > 8) {
	if (c__ <= 8 && d__ <= 8) {
	    ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &c__, &a, &opvert, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[1] = ipout;
	    if (iside[1] != 0) {
		goto L310;
	    }
	    ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &c__, &a, &opvert, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[1] = ipout;
	    goto L310;
	}
	vrtarr_(&c__, &d__, &opvert, &ib, &ic, &id);
	if (ic > 8) {
	    ifn = 0;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &ib, &ic, &a, &id, mhalf, mfull, &isclp[1],
		     &ifn, &ipout);
	    iside[1] = ipout;
	} else {
	    ifn = 1;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &ib, &id, &a, &ic, mhalf, mfull, &isclp[1],
		     &ifn, &ipout);
	    iside[1] = ipout;
	}
	if (iside[1] != 0) {
	    goto L310;
	}
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &ib, &id, &ic, 
		&a, mhalf, mfull, &isclp[1], &ipout);
	iside[1] = ipout;
L310:

	vrtarr_(&b, &opvert, &d__, &ib, &ic, &id);
	if (d__ > 8) {
	    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &ib, &id, &
		    ic, &a, mhalf, mfull, &isclp[1], &ipout);
	    iside[2] = ipout;
	    goto L320;
	}
	if (ic > 8) {
	    ifn = 0;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &ib, &ic, &a, &id, mhalf, mfull, &isclp[1],
		     &ifn, &ipout);
	    iside[2] = ipout;
	} else {
	    ifn = 1;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &ib, &id, &a, &ic, mhalf, mfull, &isclp[1],
		     &ifn, &ipout);
	    iside[2] = ipout;
	}
	if (iside[2] != 0) {
	    goto L320;
	}
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &ib, &id, &ic, 
		&a, mhalf, mfull, &isclp[1], &ipout);
	iside[2] = ipout;
L320:

	vrtarr_(&b, &c__, &opvert, &ib, &ic, &id);
	if (c__ > 8) {
	    ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &ib, &id, &
		    ic, &a, mhalf, mfull, &isclp[1], &ipout);
	    iside[3] = ipout;
	    goto L330;
	}
	if (ic > 8) {
	    ifn = 0;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &ib, &ic, &a, &id, mhalf, mfull, &isclp[1],
		     &ifn, &ipout);
	    iside[3] = ipout;
	} else {
	    ifn = 1;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &ib, &id, &a, &ic, mhalf, mfull, &isclp[1],
		     &ifn, &ipout);
	    iside[3] = ipout;
	}
	if (iside[3] != 0) {
	    goto L330;
	}
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &ib, &id, &ic, 
		&a, mhalf, mfull, &isclp[1], &ipout);
	iside[3] = ipout;
L330:

	;
    } else {
	if (c__ <= 8 && d__ <= 8) {
	    iside[1] = 1;
	} else if (c__ > 8) {
	    ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &d__, &a, &c__, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[1] = ipout;
	    if (iside[1] != 0) {
		goto L340;
	    }
	    ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &d__, &a, &c__, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[1] = ipout;
	} else {
	    ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &c__, &opvert, &a, &d__, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[1] = ipout;
	    if (iside[1] != 0) {
		goto L340;
	    }
	    ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &c__, &opvert, &a, &d__, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[1] = ipout;
	}
L340:

	if (d__ > 8) {
	    ifn = 1;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &b, &d__, &a, &opvert, mhalf, mfull, &
		    isclp[1], &ifn, &ipout);
	    iside[2] = ipout;
	} else {
	    ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &opvert, &a, &b, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[2] = ipout;
	    if (iside[2] != 0) {
		goto L350;
	    }
	    ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &opvert, &a, &b, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[2] = ipout;
	    goto L350;
	}
	if (iside[2] != 0) {
	    goto L350;
	}
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &b, &d__, &
		opvert, &a, mhalf, mfull, &isclp[1], &ipout);
	iside[2] = ipout;
L350:

	if (c__ > 8) {
	    ifn = 0;
	    ipsig3_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &b, &c__, &a, &opvert, mhalf, mfull, &
		    isclp[1], &ifn, &ipout);
	    iside[3] = ipout;
	} else {
	    ipsig4_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &c__, &a, &b, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[3] = ipout;
	    if (iside[3] != 0) {
		goto L360;
	    }
	    ipsig6_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &c__, &a, &b, mhalf, mfull, &
		    isclp[1], &ipout);
	    iside[3] = ipout;
	    goto L360;
	}
	if (iside[3] != 0) {
	    goto L360;
	}
	ipsign_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &b, &opvert, &
		c__, &a, mhalf, mfull, &isclp[1], &ipout);
	iside[3] = ipout;
L360:
	;
    }

/*     set sidist array */

    for (j = 2; j <= 4; ++j) {
	if (iside[j - 1] > 0) {
	    sidist[j - 1] = 0;
	} else if (iside[j - 1] < 0) {
	    sidist[j - 1] = -1;
	} else {
	    sidist[j - 1] = 1;
	}
/* L400: */
    }

/*     flip according to type of flip if possible */

    if (sidist[1] == 0 && sidist[2] == 0 && sidist[3] == 0) {
	if (opvert > 8) {
	    goto L420;
	}
	if (c__ <= 8 && d__ <= 8) {
	    ipsig2_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &d__, &opvert, &c__, &a, &b, 
		    mhalf, mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L420;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else if (c__ <= 8) {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &c__, &opvert, &b, &d__, &a, &b, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L420;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &opvert, &c__, &b, &a, &c__, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L420;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	}
	if (itide >= 0) {
	    goto L1000;
	}
L420:
	flip23_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], nvmax);
	goto L1000;
    } else if (sidist[1] == -1 && sidist[2] == 0 && sidist[3] == 0) {
	if (c__ <= 8 && d__ <= 8) {
	    s_stop("1660", (ftnlen)4);
	}
	if (opvert > 8) {
	    goto L440;
	}
	if (c__ <= 8) {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &c__, &opvert, &b, &d__, &a, &b, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L440;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &opvert, &c__, &b, &a, &c__, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L440;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	}
	if (itide >= 0) {
	    goto L1000;
	}
L440:
	oddsid = 2;
	flip32_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == -1 && sidist[3] == 0) {
	if (d__ > 8) {
	    goto L1000;
	}
	if (opvert > 8 && c__ > 8) {
	    goto L460;
	}
	if (opvert > 8) {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &c__, &b, &opvert, &a, &b, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L460;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else if (c__ <= 8) {
	    ipsig2_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &d__, &opvert, &c__, &a, &b, 
		    mhalf, mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L460;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &opvert, &c__, &b, &a, &c__, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L460;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	}
	if (itide >= 0) {
	    goto L1000;
	}
L460:
	oddsid = 3;
	flip32_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == 0 && sidist[3] == -1) {
	if (c__ > 8) {
	    goto L1000;
	}
	if (opvert > 8 && d__ > 8) {
	    goto L480;
	}
	if (opvert > 8) {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &c__, &b, &opvert, &a, &b, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L480;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else if (d__ <= 8) {
	    ipsig2_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &d__, &opvert, &c__, &a, &b, 
		    mhalf, mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L480;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &c__, &opvert, &b, &d__, &a, &b, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L480;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	}
	if (itide >= 0) {
	    goto L1000;
	}
L480:
	oddsid = 4;
	flip32_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 1 && sidist[2] == 0 && sidist[3] == 0) {
	if (c__ <= 8 && d__ <= 8) {
	    s_stop("1670", (ftnlen)4);
	}
	if (opvert > 8) {
	    goto L500;
	}
	if (c__ <= 8) {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &c__, &opvert, &b, &d__, &a, &b, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L500;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &opvert, &c__, &b, &a, &c__, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L500;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	}
	if (itide >= 0) {
	    goto L1000;
	}
L500:
	oddsid = 2;
	flip22_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == 1 && sidist[3] == 0) {
	if (opvert > 8 && d__ > 8) {
	    iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &b, &d__, &opvert, &a, mhalf, mfull, &isclp[1], &
		    isclw[1], &isclr[1], delaun, &itide);
	    if (itide >= 0) {
		goto L1000;
	    }
	    goto L520;
	}
	if (opvert > 8) {
	    goto L520;
	}
	if (d__ > 8) {
	    goto L1000;
	}
	if (c__ > 8) {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &d__, &opvert, &c__, &b, &a, &c__, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L520;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else {
	    ipsig2_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &d__, &opvert, &c__, &a, &b, 
		    mhalf, mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L520;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	}
	if (itide >= 0) {
	    goto L1000;
	}
L520:
	oddsid = 3;
	flip22_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (sidist[1] == 0 && sidist[2] == 0 && sidist[3] == 1) {
	if (opvert > 8 && c__ > 8) {
	    iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &b, &opvert, &c__, &a, mhalf, mfull, &isclp[1], &
		    isclw[1], &isclr[1], delaun, &itide);
	    if (itide >= 0) {
		goto L1000;
	    }
	    goto L540;
	}
	if (opvert > 8) {
	    goto L540;
	}
	if (c__ > 8) {
	    goto L1000;
	}
	if (d__ > 8) {
	    ipsig1_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &c__, &opvert, &b, &d__, &a, &b, mhalf, 
		    mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L540;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	} else {
	    ipsig2_(&x[1], &y[1], &z__[1], &x2[1], &y2[1], &z2[1], &xc[1], &
		    yc[1], &zc[1], &opvert, &d__, &opvert, &c__, &a, &b, 
		    mhalf, mfull, &isclp[1], &istt);
	    if (istt > 0) {
		goto L540;
	    }
	    if (istt < 0) {
		goto L1000;
	    }
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &opvert, &d__, &c__, &b, &a, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	}
	if (itide >= 0) {
	    goto L1000;
	}
L540:
	oddsid = 4;
	flip22_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
	goto L1000;
    } else if (*delaun) {
	goto L1000;
    }

    if (sidist[1] == 0 && sidist[2] == -1 && sidist[3] == -1) {
	oddsid = 2;
	goto L900;
    } else if (sidist[1] == -1 && sidist[2] == 0 && sidist[3] == -1) {
	oddsid = 3;
	goto L900;
    } else if (sidist[1] == -1 && sidist[2] == -1 && sidist[3] == 0) {
	oddsid = 4;
	goto L900;
    } else if (sidist[1] == 0 && sidist[2] == 1 && sidist[3] == 1) {
	oddsid = 2;
	if (opvert <= 8) {
	    goto L900;
	}
	isodd = icon[(now << 3) + 6];
	if (isodd <= 8) {
	    s_stop("1680", (ftnlen)4);
	}
	iqsig1_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&isodd, &opvert, k, mhalf, mfull, &isclp[1], &isclw[1], &
		isclr[1], delaun, &itide);
	if (itide >= 0) {
	    goto L1000;
	}
	flip21_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == 0 && sidist[3] == 1) {
	oddsid = 3;
	if (opvert <= 8) {
	    goto L900;
	}
	isodd = icon[(now << 3) + 7];
	if (isodd <= 8) {
	    s_stop("1690", (ftnlen)4);
	}
	iqsig1_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&isodd, &opvert, k, mhalf, mfull, &isclp[1], &isclw[1], &
		isclr[1], delaun, &itide);
	if (itide >= 0) {
	    goto L1000;
	}
	flip21_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == 1 && sidist[3] == 0) {
	oddsid = 4;
	if (opvert <= 8) {
	    goto L900;
	}
	isodd = icon[(now << 3) + 8];
	if (isodd <= 8) {
	    s_stop("1710", (ftnlen)4);
	}
	iqsig1_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&isodd, &opvert, k, mhalf, mfull, &isclp[1], &isclw[1], &
		isclr[1], delaun, &itide);
	if (itide >= 0) {
	    goto L1000;
	}
	flip21_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 0 && sidist[2] == -1 && sidist[3] == 1) {
	oddsid = 2;
	isite = icon[(now << 3) + 7];
	if (opvert <= 8 || isite <= 8) {
	    goto L900;
	}
	isodd = icon[(now << 3) + 6];
	if (isodd <= 8) {
	    s_stop("1720", (ftnlen)4);
	}
	iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&isodd, &isite, &opvert, k, mhalf, mfull, &isclp[1], &isclw[1]
		, &isclr[1], delaun, &itide);
	if (itide >= 0) {
	    goto L1000;
	}
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == 0 && sidist[3] == -1) {
	oddsid = 3;
	isite = icon[(now << 3) + 8];
	if (opvert <= 8 || isite <= 8) {
	    goto L900;
	}
	s_stop("1730", (ftnlen)4);
    } else if (sidist[1] == -1 && sidist[2] == 1 && sidist[3] == 0) {
	oddsid = 4;
	isite = icon[(now << 3) + 6];
	if (opvert <= 8 || isite <= 8) {
	    goto L900;
	}
	isodd = icon[(now << 3) + 8];
	if (isodd <= 8) {
	    s_stop("1740", (ftnlen)4);
	}
	iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&isodd, &isite, &opvert, k, mhalf, mfull, &isclp[1], &isclw[1]
		, &isclr[1], delaun, &itide);
	if (itide >= 0) {
	    goto L1000;
	}
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 0 && sidist[2] == 1 && sidist[3] == -1) {
	oddsid = -2;
	isite = icon[(now << 3) + 8];
	if (opvert <= 8 || isite <= 8) {
	    goto L800;
	}
	isodd = icon[(now << 3) + 6];
	if (isodd <= 8) {
	    s_stop("1750", (ftnlen)4);
	}
	iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&isodd, &isite, &opvert, k, mhalf, mfull, &isclp[1], &isclw[1]
		, &isclr[1], delaun, &itide);
	if (itide >= 0) {
	    goto L1000;
	}
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == -1 && sidist[2] == 0 && sidist[3] == 1) {
	oddsid = -3;
	isite = icon[(now << 3) + 6];
	if (opvert <= 8 || isite <= 8) {
	    goto L800;
	}
	isodd = icon[(now << 3) + 7];
	if (isodd <= 8) {
	    s_stop("1760", (ftnlen)4);
	}
	iqsig2_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[1], 
		&isodd, &isite, &opvert, k, mhalf, mfull, &isclp[1], &isclw[1]
		, &isclr[1], delaun, &itide);
	if (itide >= 0) {
	    goto L1000;
	}
	flip31_(&icon[9], k, &now, &adj, &ifl[1], &ih[1], ihn, nhmax, flphis, 
		tetra, newtts, &is[1], &oddsid, nvmax);
    } else if (sidist[1] == 1 && sidist[2] == -1 && sidist[3] == 0) {
	oddsid = -4;
	isite = icon[(now << 3) + 7];
	if (opvert <= 8 || isite <= 8) {
	    goto L800;
	}
	s_stop("1770", (ftnlen)4);
    }
    goto L1000;

L800:
    oddsid = -oddsid;
L900:
    isite = icon[oddsid + 4 + (now << 3)];
    if (isite <= 8) {
	s_stop("1780", (ftnlen)4);
    }

L1000:
    ++i__;
    goto L100;

L1500:

    return 0;
} /* sphere_ */

/* FLIP23 */

/*     This subroutine will perform a 2->3 flip. */

/* Subroutine */ int flip23_(integer *icon, integer *k, integer *now, integer 
	*adj, integer *ifl, integer *ih, integer *ihn, integer *nhmax, 
	logical *flphis, integer *tetra, integer *newtts, integer *is, 
	integer *nvmax)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, new1, new2, new3, next;



    /* Parameter adjustments */
    --is;
    --ih;
    --ifl;
    icon -= 9;

    /* Function Body */
    if (*flphis) {
	new1 = *tetra + 1;
	new2 = *tetra + 2;
	new3 = *tetra + 3;
	*tetra = new3;
	if (*tetra > *nvmax) {
	    s_stop("1810", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    *tetra = new1;
	    if (*tetra > *nvmax) {
		s_stop("1820", (ftnlen)4);
	    }
	} else {
	    new1 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new2 = *tetra + 1;
	    *tetra = new2;
	    if (*tetra > *nvmax) {
		s_stop("1830", (ftnlen)4);
	    }
	} else {
	    new2 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new3 = *tetra + 1;
	    *tetra = new3;
	    if (*tetra > *nvmax) {
		s_stop("1840", (ftnlen)4);
	    }
	} else {
	    new3 = ih[*ihn];
	    --(*ihn);
	}
    }

    *newtts += 3;
    if (*newtts > *nvmax) {
	s_stop("1845", (ftnlen)4);
    }
    ifl[*newtts - 2] = new1;
    ifl[*newtts - 1] = new2;
    ifl[*newtts] = new3;

/*     create new1 */

    icon[(new1 << 3) + 1] = icon[(*adj << 3) + 1];
    icon[(new1 << 3) + 2] = new2;
    icon[(new1 << 3) + 3] = new3;
    icon[(new1 << 3) + 4] = icon[(*now << 3) + 2];
    icon[(new1 << 3) + 5] = *k;
    icon[(new1 << 3) + 6] = icon[(*adj << 3) + 6];
    icon[(new1 << 3) + 7] = icon[(*adj << 3) + 7];
    icon[(new1 << 3) + 8] = icon[(*adj << 3) + 8];

/*     create new2 */

    icon[(new2 << 3) + 1] = icon[(*adj << 3) + 2];
    icon[(new2 << 3) + 2] = new3;
    icon[(new2 << 3) + 3] = new1;
    icon[(new2 << 3) + 4] = icon[(*now << 3) + 3];
    icon[(new2 << 3) + 5] = *k;
    icon[(new2 << 3) + 6] = icon[(*adj << 3) + 7];
    icon[(new2 << 3) + 7] = icon[(*adj << 3) + 5];
    icon[(new2 << 3) + 8] = icon[(*adj << 3) + 8];

/*     create new3 */

    icon[(new3 << 3) + 1] = icon[(*adj << 3) + 3];
    icon[(new3 << 3) + 2] = new2;
    icon[(new3 << 3) + 3] = icon[(*now << 3) + 4];
    icon[(new3 << 3) + 4] = new1;
    icon[(new3 << 3) + 5] = *k;
    icon[(new3 << 3) + 6] = icon[(*adj << 3) + 6];
    icon[(new3 << 3) + 7] = icon[(*adj << 3) + 8];
    icon[(new3 << 3) + 8] = icon[(*adj << 3) + 5];

/*     update neighboring tetrahedra */

    for (i__ = 2; i__ <= 4; ++i__) {
	next = icon[i__ + (*now << 3)];
	if (next == 0) {
	    goto L100;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (next << 3)] == *now) {
		goto L60;
	    }
/* L50: */
	}
	s_stop("1850", (ftnlen)4);
L60:
	if (i__ == 2) {
	    icon[j + (next << 3)] = new1;
	} else if (i__ == 3) {
	    icon[j + (next << 3)] = new2;
	} else {
	    icon[j + (next << 3)] = new3;
	}
L100:
	;
    }

    for (i__ = 1; i__ <= 3; ++i__) {
	next = icon[i__ + (*adj << 3)];
	if (next == 0) {
	    goto L200;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (next << 3)] == *adj) {
		goto L160;
	    }
/* L150: */
	}
	s_stop("1860", (ftnlen)4);
L160:
	if (i__ == 1) {
	    icon[j + (next << 3)] = new1;
	} else if (i__ == 2) {
	    icon[j + (next << 3)] = new2;
	} else {
	    icon[j + (next << 3)] = new3;
	}
L200:
	;
    }

/*     update is(*) */

    is[icon[(*now << 3) + 5]] = new3;
    is[icon[(*now << 3) + 6]] = new3;
    is[icon[(*now << 3) + 7]] = new3;
    is[icon[(*now << 3) + 8]] = new1;
    is[icon[(*adj << 3) + 8]] = new3;

/*     mark 2 old tetra to show children */

    if (*flphis) {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*now << 3) + 1] = new1;
	icon[(*now << 3) + 2] = new2;
	icon[(*now << 3) + 3] = new3;
	icon[(*now << 3) + 4] = 0;

	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(*adj << 3) + 1] = new1;
	icon[(*adj << 3) + 2] = new2;
	icon[(*adj << 3) + 3] = new3;
	icon[(*adj << 3) + 4] = 0;
    } else {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	*ihn += 2;
	if (*ihn > *nhmax) {
	    s_stop("1870", (ftnlen)4);
	}
	ih[*ihn] = *now;
	ih[*ihn - 1] = *adj;
    }

    return 0;
} /* flip23_ */

/* FLIP32 */

/* This subroutine will perform a 3->2 flip.  As a result, 3 new */
/* tetra will be created. */

/* Subroutine */ int flip32_(integer *icon, integer *k, integer *now, integer 
	*adj, integer *ifl, integer *ih, integer *ihn, integer *nhmax, 
	logical *flphis, integer *tetra, integer *newtts, integer *is, 
	integer *negsid, integer *nvmax)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, new1, new2, site1, site2, neigh, nxtadj;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *), sitord_(integer *, integer *, integer *);
    integer nxtnow;



/*     reorder now */

    /* Parameter adjustments */
    --is;
    --ih;
    --ifl;
    icon -= 9;

    /* Function Body */
    site1 = icon[*negsid + 4 + (*now << 3)];
    reordr_(&icon[9], k, &site1, now);

/*     check if now & adj have same neighbor, reorder adj */

    nxtnow = icon[(*now << 3) + 2];
    if (icon[(nxtnow << 3) + 5] != *k) {
	s_stop("1910", (ftnlen)4);
    }

    sitord_(&icon[9], &site1, adj);
    nxtadj = icon[(*adj << 3) + 1];
    if (nxtnow != nxtadj) {
	goto L1000;
    }

    for (i__ = 1; i__ <= 4; ++i__) {
	if (icon[i__ + (*adj << 3)] == *now) {
	    site2 = icon[i__ + 4 + (*adj << 3)];
	    goto L215;
	}
/* L210: */
    }
    s_stop("1920", (ftnlen)4);
L215:
    reordr_(&icon[9], &site1, &site2, adj);

/*     reorder nxtnow */

    reordr_(&icon[9], k, &site2, &nxtnow);

    if (*flphis) {
	new1 = *tetra + 1;
	new2 = *tetra + 2;
	*tetra = new2;
	if (*tetra > *nvmax) {
	    s_stop("1930", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    *tetra = new1;
	    if (*tetra > *nvmax) {
		s_stop("1940", (ftnlen)4);
	    }
	} else {
	    new1 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new2 = *tetra + 1;
	    *tetra = new2;
	    if (*tetra > *nvmax) {
		s_stop("1950", (ftnlen)4);
	    }
	} else {
	    new2 = ih[*ihn];
	    --(*ihn);
	}
    }

    *newtts += 2;
    if (*newtts > *nvmax) {
	s_stop("1955", (ftnlen)4);
    }
    ifl[*newtts - 1] = new1;
    ifl[*newtts] = new2;

/*     create new1 */

    icon[(new1 << 3) + 1] = icon[(*adj << 3) + 4];
    icon[(new1 << 3) + 2] = icon[(nxtnow << 3) + 3];
    icon[(new1 << 3) + 3] = new2;
    icon[(new1 << 3) + 4] = icon[(*now << 3) + 4];
    icon[(new1 << 3) + 5] = *k;
    icon[(new1 << 3) + 6] = icon[(*now << 3) + 6];
    icon[(new1 << 3) + 7] = icon[(*now << 3) + 7];
    icon[(new1 << 3) + 8] = icon[(*adj << 3) + 6];

/*     create new2 */

    icon[(new2 << 3) + 1] = icon[(*adj << 3) + 3];
    icon[(new2 << 3) + 2] = icon[(nxtnow << 3) + 4];
    icon[(new2 << 3) + 3] = icon[(*now << 3) + 3];
    icon[(new2 << 3) + 4] = new1;
    icon[(new2 << 3) + 5] = *k;
    icon[(new2 << 3) + 6] = icon[(*now << 3) + 6];
    icon[(new2 << 3) + 7] = icon[(*adj << 3) + 6];
    icon[(new2 << 3) + 8] = icon[(*now << 3) + 8];

/*     update neighboring tetrahedra */

    for (i__ = 3; i__ <= 4; ++i__) {
	neigh = icon[i__ + (*now << 3)];
	if (neigh == 0) {
	    goto L400;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (neigh << 3)] == *now) {
		goto L375;
	    }
/* L350: */
	}
	s_stop("1960", (ftnlen)4);
L375:
	if (i__ == 3) {
	    icon[j + (neigh << 3)] = new2;
	} else {
	    icon[j + (neigh << 3)] = new1;
	}
L400:
	;
    }

    for (i__ = 3; i__ <= 4; ++i__) {
	neigh = icon[i__ + (*adj << 3)];
	if (neigh == 0) {
	    goto L500;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (neigh << 3)] == *adj) {
		goto L475;
	    }
/* L450: */
	}
	s_stop("1970", (ftnlen)4);
L475:
	if (i__ == 3) {
	    icon[j + (neigh << 3)] = new2;
	} else {
	    icon[j + (neigh << 3)] = new1;
	}
L500:
	;
    }

    for (i__ = 3; i__ <= 4; ++i__) {
	neigh = icon[i__ + (nxtnow << 3)];
	if (neigh == 0) {
	    goto L600;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (neigh << 3)] == nxtnow) {
		goto L575;
	    }
/* L550: */
	}
	s_stop("1980", (ftnlen)4);
L575:
	if (i__ == 3) {
	    icon[j + (neigh << 3)] = new1;
	} else {
	    icon[j + (neigh << 3)] = new2;
	}
L600:
	;
    }

/*     update is(*) */

    is[icon[(*now << 3) + 5]] = new1;
    is[icon[(*now << 3) + 6]] = new1;
    is[icon[(*now << 3) + 7]] = new1;
    is[icon[(*now << 3) + 8]] = new2;
    is[icon[(*adj << 3) + 6]] = new1;

/*     show children of adj, now, nxtnow */

    if (*flphis) {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*now << 3) + 1] = new1;
	icon[(*now << 3) + 2] = new2;
	icon[(*now << 3) + 3] = 0;
	icon[(*now << 3) + 4] = 0;

	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(*adj << 3) + 1] = new1;
	icon[(*adj << 3) + 2] = new2;
	icon[(*adj << 3) + 3] = 0;
	icon[(*adj << 3) + 4] = 0;

	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	icon[(nxtnow << 3) + 1] = new1;
	icon[(nxtnow << 3) + 2] = new2;
	icon[(nxtnow << 3) + 3] = 0;
	icon[(nxtnow << 3) + 4] = 0;
    } else {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	*ihn += 3;
	if (*ihn > *nhmax) {
	    s_stop("1990", (ftnlen)4);
	}
	ih[*ihn] = *now;
	ih[*ihn - 1] = *adj;
	ih[*ihn - 2] = nxtnow;
    }

L1000:
    return 0;
} /* flip32_ */

/* FLIP22 */

/* This subroutine will perform a 2->2 flip.  Four new tetra will be */
/* created. */

/* Subroutine */ int flip22_(integer *icon, integer *k, integer *now, integer 
	*adj, integer *ifl, integer *ih, integer *ihn, integer *nhmax, 
	logical *flphis, integer *tetra, integer *newtts, integer *is, 
	integer *zersid, integer *nvmax)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, new1, new2, new3, new4, site1, site2, neigh, nxtadj;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *), sitord_(integer *, integer *, integer *);
    integer nxtnow;



/*     reorder now */

    /* Parameter adjustments */
    --is;
    --ih;
    --ifl;
    icon -= 9;

    /* Function Body */
    site1 = icon[*zersid + 4 + (*now << 3)];
    reordr_(&icon[9], k, &site1, now);

/*     define nxtnow */

    nxtnow = icon[(*now << 3) + 2];
    if (icon[(nxtnow << 3) + 5] != *k) {
	s_stop("2010", (ftnlen)4);
    }

/*     reorder adj */

    sitord_(&icon[9], &site1, adj);

/*     define nxtadj */

    nxtadj = icon[(*adj << 3) + 1];

/*     are nxtnow and nxtadj neighbors? */

    for (i__ = 1; i__ <= 4; ++i__) {
	if (icon[i__ + (nxtnow << 3)] == nxtadj) {
	    goto L6;
	}
/* L5: */
    }
    goto L2000;
L6:

    if (*flphis) {
	new1 = *tetra + 1;
	new2 = *tetra + 2;
	new3 = *tetra + 3;
	new4 = *tetra + 4;
	*tetra = new4;
	if (*tetra > *nvmax) {
	    s_stop("2020", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    *tetra = new1;
	    if (*tetra > *nvmax) {
		s_stop("2030", (ftnlen)4);
	    }
	} else {
	    new1 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new2 = *tetra + 1;
	    *tetra = new2;
	    if (*tetra > *nvmax) {
		s_stop("2040", (ftnlen)4);
	    }
	} else {
	    new2 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new3 = *tetra + 1;
	    *tetra = new3;
	    if (*tetra > *nvmax) {
		s_stop("2050", (ftnlen)4);
	    }
	} else {
	    new3 = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    new4 = *tetra + 1;
	    *tetra = new4;
	    if (*tetra > *nvmax) {
		s_stop("2060", (ftnlen)4);
	    }
	} else {
	    new4 = ih[*ihn];
	    --(*ihn);
	}
    }

    *newtts += 4;
    if (*newtts > *nvmax) {
	s_stop("2065", (ftnlen)4);
    }
    ifl[*newtts - 3] = new1;
    ifl[*newtts - 2] = new2;
    ifl[*newtts - 1] = new3;
    ifl[*newtts] = new4;

/*     reorder adj, nxtnow, nxtadj */

    for (i__ = 2; i__ <= 4; ++i__) {
	if (icon[i__ + (*adj << 3)] == *now) {
	    site2 = icon[i__ + 4 + (*adj << 3)];
	    goto L15;
	}
/* L10: */
    }
    s_stop("2070", (ftnlen)4);
L15:
    reordr_(&icon[9], &site1, &site2, adj);

    for (i__ = 2; i__ <= 4; ++i__) {
	if (icon[i__ + (nxtnow << 3)] == *now) {
	    site1 = icon[i__ + 4 + (nxtnow << 3)];
	    goto L25;
	}
/* L20: */
    }
    s_stop("2080", (ftnlen)4);
L25:
    reordr_(&icon[9], k, &site1, &nxtnow);

    reordr_(&icon[9], &site1, &site2, &nxtadj);

/*     create new1 */

    icon[(new1 << 3) + 1] = icon[(*adj << 3) + 4];
    icon[(new1 << 3) + 2] = new3;
    icon[(new1 << 3) + 3] = new2;
    icon[(new1 << 3) + 4] = icon[(*now << 3) + 4];
    icon[(new1 << 3) + 5] = *k;
    icon[(new1 << 3) + 6] = icon[(*now << 3) + 6];
    icon[(new1 << 3) + 7] = icon[(*now << 3) + 7];
    icon[(new1 << 3) + 8] = icon[(*adj << 3) + 6];

/*     create new2 */

    icon[(new2 << 3) + 1] = icon[(*adj << 3) + 3];
    icon[(new2 << 3) + 2] = new4;
    icon[(new2 << 3) + 3] = icon[(*now << 3) + 3];
    icon[(new2 << 3) + 4] = new1;
    icon[(new2 << 3) + 5] = *k;
    icon[(new2 << 3) + 6] = icon[(*now << 3) + 6];
    icon[(new2 << 3) + 7] = icon[(*adj << 3) + 6];
    icon[(new2 << 3) + 8] = icon[(*now << 3) + 8];

/*     create new3 */

    icon[(new3 << 3) + 1] = icon[(nxtadj << 3) + 3];
    icon[(new3 << 3) + 2] = new4;
    icon[(new3 << 3) + 3] = new1;
    icon[(new3 << 3) + 4] = icon[(nxtnow << 3) + 3];
    icon[(new3 << 3) + 5] = *k;
    icon[(new3 << 3) + 6] = icon[(*now << 3) + 7];
    icon[(new3 << 3) + 7] = icon[(nxtnow << 3) + 6];
    icon[(new3 << 3) + 8] = icon[(*adj << 3) + 6];

/*     create new4 */

    icon[(new4 << 3) + 1] = icon[(nxtadj << 3) + 4];
    icon[(new4 << 3) + 2] = new2;
    icon[(new4 << 3) + 3] = new3;
    icon[(new4 << 3) + 4] = icon[(nxtnow << 3) + 4];
    icon[(new4 << 3) + 5] = *k;
    icon[(new4 << 3) + 6] = icon[(nxtnow << 3) + 6];
    icon[(new4 << 3) + 7] = icon[(*now << 3) + 8];
    icon[(new4 << 3) + 8] = icon[(*adj << 3) + 6];

/*     update is(*) */

    is[icon[(*now << 3) + 5]] = new1;
    is[icon[(*now << 3) + 6]] = new1;
    is[icon[(*now << 3) + 7]] = new1;
    is[icon[(*now << 3) + 8]] = new2;
    is[icon[(*adj << 3) + 6]] = new1;
    is[icon[(nxtnow << 3) + 6]] = new3;

/*     update neighbors of now */

    for (i__ = 3; i__ <= 4; ++i__) {
	neigh = icon[i__ + (*now << 3)];
	if (neigh == 0) {
	    goto L100;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (neigh << 3)] == *now) {
		goto L75;
	    }
/* L50: */
	}
	s_stop("2090", (ftnlen)4);
L75:
	if (i__ == 3) {
	    icon[j + (neigh << 3)] = new2;
	} else {
	    icon[j + (neigh << 3)] = new1;
	}
L100:
	;
    }

/*     update neighbors of adj */

    for (i__ = 3; i__ <= 4; ++i__) {
	neigh = icon[i__ + (*adj << 3)];
	if (neigh == 0) {
	    goto L400;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (neigh << 3)] == *adj) {
		goto L375;
	    }
/* L350: */
	}
	s_stop("2110", (ftnlen)4);
L375:
	if (i__ == 3) {
	    icon[j + (neigh << 3)] = new2;
	} else {
	    icon[j + (neigh << 3)] = new1;
	}
L400:
	;
    }

/*     update neighbors of nxtnow */

    for (i__ = 3; i__ <= 4; ++i__) {
	neigh = icon[i__ + (nxtnow << 3)];
	if (neigh == 0) {
	    goto L600;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (neigh << 3)] == nxtnow) {
		goto L590;
	    }
/* L575: */
	}
	s_stop("2120", (ftnlen)4);
L590:
	if (i__ == 3) {
	    icon[j + (neigh << 3)] = new3;
	} else {
	    icon[j + (neigh << 3)] = new4;
	}
L600:
	;
    }

/*     update neighbors of nxtadj */

    for (i__ = 3; i__ <= 4; ++i__) {
	neigh = icon[i__ + (nxtadj << 3)];
	if (neigh == 0) {
	    goto L900;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (neigh << 3)] == nxtadj) {
		goto L890;
	    }
/* L875: */
	}
	s_stop("2130", (ftnlen)4);
L890:
	if (i__ == 3) {
	    icon[j + (neigh << 3)] = new3;
	} else {
	    icon[j + (neigh << 3)] = new4;
	}
L900:
	;
    }

/*     show children of old tetra */

    if (*flphis) {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*now << 3) + 1] = new1;
	icon[(*now << 3) + 2] = new2;
	icon[(*now << 3) + 3] = -nxtnow;
	icon[(*now << 3) + 4] = 0;

	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(*adj << 3) + 1] = new1;
	icon[(*adj << 3) + 2] = new2;
	icon[(*adj << 3) + 3] = -nxtadj;
	icon[(*adj << 3) + 4] = 0;

	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	icon[(nxtnow << 3) + 1] = new3;
	icon[(nxtnow << 3) + 2] = new4;
	icon[(nxtnow << 3) + 3] = -(*now);
	icon[(nxtnow << 3) + 4] = 0;

	icon[(nxtadj << 3) + 5] = -icon[(nxtadj << 3) + 5];
	icon[(nxtadj << 3) + 1] = new3;
	icon[(nxtadj << 3) + 2] = new4;
	icon[(nxtadj << 3) + 3] = -(*adj);
	icon[(nxtadj << 3) + 4] = 0;
    } else {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	icon[(nxtadj << 3) + 5] = -icon[(nxtadj << 3) + 5];
	*ihn += 4;
	if (*ihn > *nhmax) {
	    s_stop("2140", (ftnlen)4);
	}
	ih[*ihn] = *now;
	ih[*ihn - 1] = *adj;
	ih[*ihn - 2] = nxtnow;
	ih[*ihn - 3] = nxtadj;
    }

L2000:
    return 0;
} /* flip22_ */

/* FLIP41 */

/* This subroutine will perform a 4->1 flip. 1 tetrahedron will be */
/* created from 4 tetrahedra */

/* Subroutine */ int flip41_(integer *icon, integer *k, integer *now, integer 
	*adj, integer *ifl, integer *ih, integer *ihn, integer *nhmax, 
	logical *flphis, integer *tetra, integer *newtts, integer *is, 
	integer *zersid, integer *nvmax)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer j, new1, site1, site2, neigh, nxtadj;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *);
    integer nxtnow;



/*     reorder now */

    /* Parameter adjustments */
    --is;
    --ih;
    --ifl;
    icon -= 9;

    /* Function Body */
    site1 = icon[*zersid + 4 + (*now << 3)];
    reordr_(&icon[9], k, &site1, now);
    site2 = icon[(*now << 3) + 7];

/*     reorder adj */

    reordr_(&icon[9], &site1, &site2, adj);

/*     define nxtnow and nxtadj */

    nxtnow = icon[(*now << 3) + 4];
    nxtadj = icon[(*now << 3) + 3];

/*     do now, adj, nxtnow, nxtadj form a tetrahedron? */

    if (icon[(*adj << 3) + 3] != nxtnow || icon[(*adj << 3) + 2] != nxtadj) {
	goto L2000;
    }

/*     flip */

    if (*flphis) {
	new1 = *tetra + 1;
	*tetra = new1;
	if (*tetra > *nvmax) {
	    s_stop("2210", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    *tetra = new1;
	    if (*tetra > *nvmax) {
		s_stop("2220", (ftnlen)4);
	    }
	} else {
	    new1 = ih[*ihn];
	    --(*ihn);
	}
    }

    ++(*newtts);
    if (*newtts > *nvmax) {
	s_stop("2225", (ftnlen)4);
    }
    ifl[*newtts] = new1;

/*     reorder nxtnow and nxtadj */

    reordr_(&icon[9], k, &site2, &nxtnow);
    reordr_(&icon[9], k, &site1, &nxtadj);
    if (icon[(nxtnow << 3) + 2] != nxtadj) {
	s_stop("2230", (ftnlen)4);
    }

/*     create tetra */

    icon[(new1 << 3) + 1] = icon[(*adj << 3) + 1];
    icon[(new1 << 3) + 2] = icon[(nxtadj << 3) + 2];
    icon[(new1 << 3) + 3] = icon[(nxtnow << 3) + 3];
    icon[(new1 << 3) + 4] = icon[(*now << 3) + 2];
    icon[(new1 << 3) + 5] = *k;
    icon[(new1 << 3) + 6] = site2;
    icon[(new1 << 3) + 7] = icon[(*now << 3) + 8];
    icon[(new1 << 3) + 8] = icon[(*adj << 3) + 8];

/*     update is(*) */

    is[site1] = -4;
    is[icon[(new1 << 3) + 5]] = new1;
    is[icon[(new1 << 3) + 6]] = new1;
    is[icon[(new1 << 3) + 7]] = new1;
    is[icon[(new1 << 3) + 8]] = new1;

/*     update neighbor of now */

    neigh = icon[(*now << 3) + 2];
    if (neigh == 0) {
	goto L200;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == *now) {
	    goto L160;
	}
/* L140: */
    }
    s_stop("2240", (ftnlen)4);
L160:
    icon[j + (neigh << 3)] = new1;
L200:

/*     update neighbor of adj */

    neigh = icon[(*adj << 3) + 1];
    if (neigh == 0) {
	goto L300;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == *adj) {
	    goto L260;
	}
/* L240: */
    }
    s_stop("2250", (ftnlen)4);
L260:
    icon[j + (neigh << 3)] = new1;
L300:

/*     update neighbor of nxtnow */

    neigh = icon[(nxtnow << 3) + 3];
    if (neigh == 0) {
	goto L400;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nxtnow) {
	    goto L360;
	}
/* L340: */
    }
    s_stop("2260", (ftnlen)4);
L360:
    icon[j + (neigh << 3)] = new1;
L400:

/*     update neighbor of nxtadj */

    neigh = icon[(nxtadj << 3) + 2];
    if (neigh == 0) {
	goto L500;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nxtadj) {
	    goto L460;
	}
/* L440: */
    }
    s_stop("2270", (ftnlen)4);
L460:
    icon[j + (neigh << 3)] = new1;
L500:

/*     show children of old tetrahedra */

    if (*flphis) {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*now << 3) + 1] = new1;
	icon[(*now << 3) + 2] = 0;
	icon[(*now << 3) + 3] = 0;
	icon[(*now << 3) + 4] = 0;

	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(*adj << 3) + 1] = new1;
	icon[(*adj << 3) + 2] = 0;
	icon[(*adj << 3) + 3] = 0;
	icon[(*adj << 3) + 4] = 0;

	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	icon[(nxtnow << 3) + 1] = new1;
	icon[(nxtnow << 3) + 2] = 0;
	icon[(nxtnow << 3) + 3] = 0;
	icon[(nxtnow << 3) + 4] = 0;

	icon[(nxtadj << 3) + 5] = -icon[(nxtadj << 3) + 5];
	icon[(nxtadj << 3) + 1] = new1;
	icon[(nxtadj << 3) + 2] = 0;
	icon[(nxtadj << 3) + 3] = 0;
	icon[(nxtadj << 3) + 4] = 0;
    } else {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	icon[(nxtadj << 3) + 5] = -icon[(nxtadj << 3) + 5];
	*ihn += 4;
	if (*ihn > *nhmax) {
	    s_stop("2280", (ftnlen)4);
	}
	ih[*ihn] = *now;
	ih[*ihn - 1] = *adj;
	ih[*ihn - 2] = nxtnow;
	ih[*ihn - 3] = nxtadj;
    }

L2000:
    return 0;
} /* flip41_ */

/* FLIP21 */

/* This subroutine will perform 2->1 flips. 1 tetrahedron will */
/* be created from 2 tetrahedra for each flip */

/* Subroutine */ int flip21_(integer *icon, integer *k, integer *now, integer 
	*adj, integer *ifl, integer *ih, integer *ihn, integer *nhmax, 
	logical *flphis, integer *tetra, integer *newtts, integer *is, 
	integer *zersid, integer *nvmax)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer j, nel1, new1, site1, site2, site3, site4, neigh, lstadj, nxtadj, 
	    initet;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *);
    integer lstnow, nxtnow;



/*     reorder now */

    /* Parameter adjustments */
    --is;
    --ih;
    --ifl;
    icon -= 9;

    /* Function Body */
    site1 = icon[*zersid + 4 + (*now << 3)];
    reordr_(&icon[9], k, &site1, now);
    site2 = icon[(*now << 3) + 7];
    site3 = icon[(*now << 3) + 8];

/*     reordr adj */

    reordr_(&icon[9], &site1, &site2, adj);

    nxtnow = icon[(*now << 3) + 3];
    nxtadj = icon[(*adj << 3) + 2];
    if (nxtnow == nxtadj) {
	s_stop("2310", (ftnlen)4);
    }

/*     go around edge to test for flipping */

L100:
    reordr_(&icon[9], k, &site1, &nxtnow);
    reordr_(&icon[9], &site1, &site3, &nxtadj);
    if (icon[(nxtnow << 3) + 1] != nxtadj) {
	goto L2000;
    }
    site3 = icon[(nxtnow << 3) + 8];
    nxtnow = icon[(nxtnow << 3) + 3];
    nxtadj = icon[(nxtadj << 3) + 2];
    if (nxtnow == *now) {
	goto L200;
    }
    if (nxtadj == *adj) {
	s_stop("2320", (ftnlen)4);
    }
    goto L100;

/*     flip */

L200:

    if (nxtadj != *adj) {
	s_stop("2330", (ftnlen)4);
    }
    if (*flphis) {
	initet = *tetra + 1;
    } else {
	if (*ihn == 0) {
	    initet = *tetra + 1;
	} else {
	    initet = ih[*ihn];
	}
    }
    new1 = initet;
    site4 = icon[(*adj << 3) + 8];

/*     go around edge for creating new tetrahedra */

L300:
    nel1 = new1;
    if (*flphis) {
	new1 = *tetra + 1;
	*tetra = new1;
	if (*tetra > *nvmax) {
	    s_stop("2340", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    new1 = *tetra + 1;
	    *tetra = new1;
	    if (*tetra > *nvmax) {
		s_stop("2350", (ftnlen)4);
	    }
	} else {
	    new1 = ih[*ihn];
	    --(*ihn);
	}
    }

    ++(*newtts);
    if (*newtts > *nvmax) {
	s_stop("2355", (ftnlen)4);
    }
    ifl[*newtts] = new1;

/*     create tetra */

    icon[(new1 << 3) + 1] = icon[(nxtadj << 3) + 1];
    icon[(new1 << 3) + 2] = icon[(nxtnow << 3) + 2];
    icon[(nel1 << 3) + 3] = new1;
    icon[(new1 << 3) + 4] = nel1;
    icon[(new1 << 3) + 5] = *k;
    icon[(new1 << 3) + 6] = site4;
    icon[(new1 << 3) + 7] = icon[(nxtnow << 3) + 7];
    icon[(new1 << 3) + 8] = icon[(nxtnow << 3) + 8];

/*     update is(*) */

    is[icon[(new1 << 3) + 8]] = new1;

/*     update neighbor of nxtnow */

    neigh = icon[(nxtnow << 3) + 2];
    if (neigh == 0) {
	goto L400;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nxtnow) {
	    goto L360;
	}
/* L340: */
    }
    s_stop("2360", (ftnlen)4);
L360:
    icon[j + (neigh << 3)] = new1;
L400:

/*     update neighbor of nxtadj */

    neigh = icon[(nxtadj << 3) + 1];
    if (neigh == 0) {
	goto L500;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nxtadj) {
	    goto L460;
	}
/* L440: */
    }
    s_stop("2370", (ftnlen)4);
L460:
    icon[j + (neigh << 3)] = new1;
L500:

/*     show children of nxtnow, nxtadj */

    lstnow = nxtnow;
    lstadj = nxtadj;
    nxtnow = icon[(lstnow << 3) + 3];
    nxtadj = icon[(lstadj << 3) + 2];

    if (*flphis) {
	icon[(lstnow << 3) + 5] = -icon[(lstnow << 3) + 5];
	icon[(lstnow << 3) + 1] = new1;
	icon[(lstnow << 3) + 2] = 0;
	icon[(lstnow << 3) + 3] = -nxtnow;
	icon[(lstnow << 3) + 4] = -icon[(lstnow << 3) + 4];

	icon[(lstadj << 3) + 5] = -icon[(lstadj << 3) + 5];
	icon[(lstadj << 3) + 1] = new1;
	icon[(lstadj << 3) + 2] = -nxtadj;
	icon[(lstadj << 3) + 3] = -icon[(lstadj << 3) + 3];
	icon[(lstadj << 3) + 4] = 0;
    } else {
	icon[(lstnow << 3) + 5] = -icon[(lstnow << 3) + 5];
	icon[(lstadj << 3) + 5] = -icon[(lstadj << 3) + 5];
	*ihn += 2;
	if (*ihn > *nhmax) {
	    s_stop("2380", (ftnlen)4);
	}
	ih[*ihn] = lstnow;
	ih[*ihn - 1] = lstadj;
    }

    if (nxtnow != *now) {
	goto L300;
    }
    icon[(new1 << 3) + 3] = initet;
    icon[(initet << 3) + 4] = new1;
    is[*k] = new1;
    is[site4] = new1;
    is[site1] = -4;

L2000:
    return 0;
} /* flip21_ */

/* FLIP31 */

/* This subroutine will perform 3->1 flips. 1 tetrahedron will */
/* be created from 3 tetrahedra for each flip */

/* Subroutine */ int flip31_(integer *icon, integer *k, integer *now, integer 
	*adj, integer *ifl, integer *ih, integer *ihn, integer *nhmax, 
	logical *flphis, integer *tetra, integer *newtts, integer *is, 
	integer *zersid, integer *nvmax)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer j, nrd, nrt, site1, site2, site3, site4, site5, neigh, nside, 
	    nitat, nitet, nxtadj;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *);
    integer nxtnrd, nxtnow;



    /* Parameter adjustments */
    --is;
    --ih;
    --ifl;
    icon -= 9;

    /* Function Body */
    nside = abs(*zersid);

/*     reorder now */

    site1 = icon[nside + 4 + (*now << 3)];
    reordr_(&icon[9], k, &site1, now);

    if (*zersid > 0) {
	goto L10;
    }

/*     define nrt */

    nrt = icon[(*now << 3) + 4];

/*     reordr nrt */

    reordr_(&icon[9], k, &site1, &nrt);
    if (icon[(nrt << 3) + 1] != *adj) {
	goto L2000;
    }

/*     define nrd, redefine now */

    nrd = *now;
    *now = nrt;
    goto L20;

/*     define nrd */

L10:
    nrd = icon[(*now << 3) + 3];

/*     reorder nrd */

    reordr_(&icon[9], k, &site1, &nrd);
    if (icon[(nrd << 3) + 1] != *adj) {
	goto L2000;
    }

L20:

/*     reorder adj */

    site2 = icon[(*now << 3) + 7];
    reordr_(&icon[9], &site1, &site2, adj);

/*     define nxtnow, nxtadj, nxtnrd, and reorder */

    nxtnow = icon[(*now << 3) + 4];
    nxtadj = icon[(*adj << 3) + 3];
    nxtnrd = icon[(nrd << 3) + 3];
    reordr_(&icon[9], k, &site1, &nxtnow);
    reordr_(&icon[9], &site1, &site2, &nxtadj);
    reordr_(&icon[9], k, &site1, &nxtnrd);
    site5 = icon[(nxtnow << 3) + 7];
    if (icon[(nxtadj << 3) + 8] != site5 || icon[(nxtnrd << 3) + 8] != site5) 
	    {
	goto L2000;
    }

/*     flip */

    site3 = icon[(*now << 3) + 8];
    site4 = icon[(*adj << 3) + 8];
    if (*flphis) {
	nitet = *tetra + 1;
	nitat = *tetra + 2;
	*tetra = nitat;
	if (*tetra > *nvmax) {
	    s_stop("2410", (ftnlen)4);
	}
    } else {
	if (*ihn == 0) {
	    nitet = *tetra + 1;
	    *tetra = nitet;
	    if (*tetra > *nvmax) {
		s_stop("2420", (ftnlen)4);
	    }
	} else {
	    nitet = ih[*ihn];
	    --(*ihn);
	}
	if (*ihn == 0) {
	    nitat = *tetra + 1;
	    *tetra = nitat;
	    if (*tetra > *nvmax) {
		s_stop("2430", (ftnlen)4);
	    }
	} else {
	    nitat = ih[*ihn];
	    --(*ihn);
	}
    }

    *newtts += 2;
    if (*newtts > *nvmax) {
	s_stop("2435", (ftnlen)4);
    }
    ifl[*newtts - 1] = nitet;
    ifl[*newtts] = nitat;

/*     create new tetrahedra */

    icon[(nitet << 3) + 1] = icon[(*adj << 3) + 1];
    icon[(nitet << 3) + 2] = icon[(nrd << 3) + 2];
    icon[(nitet << 3) + 3] = nitat;
    icon[(nitet << 3) + 4] = icon[(*now << 3) + 2];
    icon[(nitet << 3) + 5] = *k;
    icon[(nitet << 3) + 6] = site2;
    icon[(nitet << 3) + 7] = site3;
    icon[(nitet << 3) + 8] = site4;

    icon[(nitat << 3) + 1] = icon[(nxtadj << 3) + 1];
    icon[(nitat << 3) + 2] = icon[(nxtnrd << 3) + 2];
    icon[(nitat << 3) + 3] = icon[(nxtnow << 3) + 2];
    icon[(nitat << 3) + 4] = nitet;
    icon[(nitat << 3) + 5] = *k;
    icon[(nitat << 3) + 6] = site2;
    icon[(nitat << 3) + 7] = site4;
    icon[(nitat << 3) + 8] = site5;

/*     update is(*) */

    is[*k] = nitat;
    is[site1] = -4;
    is[site2] = nitat;
    is[site3] = nitet;
    is[site4] = nitat;
    is[site5] = nitat;

/*     update neighbors of adj, nrd, now */

    neigh = icon[(*now << 3) + 2];
    if (neigh == 0) {
	goto L200;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == *now) {
	    goto L160;
	}
/* L140: */
    }
    s_stop("2440", (ftnlen)4);
L160:
    icon[j + (neigh << 3)] = nitet;
L200:

    neigh = icon[(*adj << 3) + 1];
    if (neigh == 0) {
	goto L300;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == *adj) {
	    goto L260;
	}
/* L240: */
    }
    s_stop("2450", (ftnlen)4);
L260:
    icon[j + (neigh << 3)] = nitet;
L300:

    neigh = icon[(nrd << 3) + 2];
    if (neigh == 0) {
	goto L330;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nrd) {
	    goto L320;
	}
/* L310: */
    }
    s_stop("2460", (ftnlen)4);
L320:
    icon[j + (neigh << 3)] = nitet;
L330:

/*     update neighbors of nxtnow, nxtadj, nxtnrd */

    neigh = icon[(nxtnow << 3) + 2];
    if (neigh == 0) {
	goto L400;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nxtnow) {
	    goto L360;
	}
/* L340: */
    }
    s_stop("2470", (ftnlen)4);
L360:
    icon[j + (neigh << 3)] = nitat;
L400:

    neigh = icon[(nxtadj << 3) + 1];
    if (neigh == 0) {
	goto L500;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nxtadj) {
	    goto L460;
	}
/* L440: */
    }
    s_stop("2480", (ftnlen)4);
L460:
    icon[j + (neigh << 3)] = nitat;
L500:

    neigh = icon[(nxtnrd << 3) + 2];
    if (neigh == 0) {
	goto L600;
    }
    for (j = 1; j <= 4; ++j) {
	if (icon[j + (neigh << 3)] == nxtnrd) {
	    goto L560;
	}
/* L540: */
    }
    s_stop("2490", (ftnlen)4);
L560:
    icon[j + (neigh << 3)] = nitat;
L600:

/*     show children of old tetrahedra */

    if (*flphis) {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*now << 3) + 1] = nitet;
	icon[(*now << 3) + 2] = -nxtnow;
	icon[(*now << 3) + 3] = 0;
	icon[(*now << 3) + 4] = 0;

	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(*adj << 3) + 1] = nitet;
	icon[(*adj << 3) + 2] = -nxtadj;
	icon[(*adj << 3) + 3] = 0;
	icon[(*adj << 3) + 4] = 0;

	icon[(nrd << 3) + 5] = -icon[(nrd << 3) + 5];
	icon[(nrd << 3) + 1] = nitet;
	icon[(nrd << 3) + 2] = -nxtnrd;
	icon[(nrd << 3) + 3] = 0;
	icon[(nrd << 3) + 4] = 0;

	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	icon[(nxtnow << 3) + 1] = nitat;
	icon[(nxtnow << 3) + 2] = -(*now);
	icon[(nxtnow << 3) + 3] = 0;
	icon[(nxtnow << 3) + 4] = 0;

	icon[(nxtadj << 3) + 5] = -icon[(nxtadj << 3) + 5];
	icon[(nxtadj << 3) + 1] = nitat;
	icon[(nxtadj << 3) + 2] = -(*adj);
	icon[(nxtadj << 3) + 3] = 0;
	icon[(nxtadj << 3) + 4] = 0;

	icon[(nxtnrd << 3) + 5] = -icon[(nxtnrd << 3) + 5];
	icon[(nxtnrd << 3) + 1] = nitat;
	icon[(nxtnrd << 3) + 2] = -nrd;
	icon[(nxtnrd << 3) + 3] = 0;
	icon[(nxtnrd << 3) + 4] = 0;
    } else {
	icon[(*now << 3) + 5] = -icon[(*now << 3) + 5];
	icon[(*adj << 3) + 5] = -icon[(*adj << 3) + 5];
	icon[(nrd << 3) + 5] = -icon[(nrd << 3) + 5];
	icon[(nxtnow << 3) + 5] = -icon[(nxtnow << 3) + 5];
	icon[(nxtadj << 3) + 5] = -icon[(nxtadj << 3) + 5];
	icon[(nxtnrd << 3) + 5] = -icon[(nxtnrd << 3) + 5];
	*ihn += 6;
	if (*ihn > *nhmax) {
	    s_stop("2495", (ftnlen)4);
	}
	ih[*ihn] = *now;
	ih[*ihn - 1] = *adj;
	ih[*ihn - 2] = nrd;
	ih[*ihn - 3] = nxtnow;
	ih[*ihn - 4] = nxtadj;
	ih[*ihn - 5] = nxtnrd;
    }

L2000:
    return 0;
} /* flip31_ */

/* CONVEX */

/* This subroutine will classify all tetra in array ifl, where: */
/* 	ifl(curr) = -1 -> tetra has children */
/* 	ifl(curr) =  0 -> tetra is outside convex hull */
/* 	ifl(curr) =  1 -> tetra is inside convex hull */
/* Then, a verification of the surface's convexity will be run. */

/* Subroutine */ int convex_(integer *icon, integer *tetra, integer *ifl, 
	doublereal *xi, doublereal *yi, doublereal *zi, integer *x, integer *
	y, integer *z__, integer *x2, integer *y2, integer *z2, integer *
	idmin, integer *mhalf, integer *mfull, integer *isclp, doublereal *
	epz)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer a, b, c__, i__, j, adj, nop, now, ikon[8]	/* was [8][1] */, 
	    curr, vert, site1, site2, site3, iside;
    extern /* Subroutine */ int irsign_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *), reordr_(
	    integer *, integer *, integer *, integer *), sitord_(integer *, 
	    integer *, integer *);



    /* Parameter adjustments */
    --isclp;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;
    --ifl;
    icon -= 9;

    /* Function Body */
    *idmin = 0;

/*     classify all tetra */

    i__1 = *tetra;
    for (curr = 1; curr <= i__1; ++curr) {
	if (icon[(curr << 3) + 5] < 0) {
	    ifl[curr] = -1;
	} else if (icon[(curr << 3) + 5] <= 8 || icon[(curr << 3) + 6] <= 8 ||
		 icon[(curr << 3) + 7] <= 8 || icon[(curr << 3) + 8] <= 8) {
	    ifl[curr] = 0;
	} else {
	    ifl[curr] = 1;
	}
/* L100: */
    }

/*     take all tetra s.t. ifl(tetra) = 1 and check convexity */

    i__1 = *tetra;
    for (curr = 1; curr <= i__1; ++curr) {
	if (ifl[curr] != 1) {
	    goto L300;
	}
	for (i__ = 1; i__ <= 4; ++i__) {
	    adj = icon[i__ + (curr << 3)];
	    if (adj == 0) {
		s_stop("2510", (ftnlen)4);
	    }
	    if (ifl[adj] != 0) {
		goto L200;
	    }
	    for (j = 1; j <= 8; ++j) {
		ikon[j - 1] = icon[j + (curr << 3)];
/* L150: */
	    }
	    site1 = icon[i__ + 4 + (curr << 3)];
	    sitord_(ikon, &site1, &c__1);
	    a = ikon[5];
	    b = ikon[6];
	    c__ = ikon[7];
	    site1 = a;
	    site2 = b;
	    site3 = c__;
	    reordr_(ikon, &site1, &site2, &c__1);
	    for (j = 1; j <= 3; ++j) {
		nop = curr;
		now = ikon[3];
L160:
		if (now == curr || now == 0) {
		    s_stop("2520", (ftnlen)4);
		}
		reordr_(&icon[9], &site1, &site2, &now);
		if (ifl[now] != 0) {
		    nop = now;
		    now = icon[(nop << 3) + 4];
		    goto L160;
		} else {
		    if (nop == curr) {
			goto L170;
		    }
		    vert = icon[(now << 3) + 8];
		    irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &
			    x2[1], &y2[1], &z2[1], &vert, &site1, &site2, &
			    site3, mhalf, mfull, &isclp[1], epz, &iside);
		    if (iside > 0) {
			++(*idmin);
		    }
		}

L170:
		if (j == 1) {
		    site1 = b;
		    site2 = c__;
		    site3 = a;
		    reordr_(ikon, &site1, &site2, &c__1);
		}
		if (j == 2) {
		    site1 = c__;
		    site2 = a;
		    site3 = b;
		    reordr_(ikon, &site1, &site2, &c__1);
		}
/* L175: */
	    }
L200:
	    ;
	}
L300:
	;
    }

    return 0;
} /* convex_ */

/* DELCHK */

/* This subroutine will test how well the Regular/Delaunay property is */
/* satisfied by the tetrahedra inside convex hull of point set */

/* Subroutine */ int delchk_(integer *tetra, integer *icon, integer *ifl, 
	doublereal *xi, doublereal *yi, doublereal *zi, doublereal *wi, 
	integer *x, integer *y, integer *z__, integer *w, integer *x2, 
	integer *y2, integer *z2, integer *w2, integer *idmax, logical *
	delaun, integer *mhalf, integer *mfull, integer *isclp, integer *
	isclw, integer *isclr, doublereal *epz)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer a, b, c__, d__, i__, j, k, adj, ikon[8]	/* was [8][1] */;
    doublereal xctr, yctr, zctr;
    integer site1, site2, site3;
    extern /* Subroutine */ int ctrad_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, logical 
	    *, integer *);
    integer itide, isite;
    doublereal tdist;
    extern /* Subroutine */ int bisphr_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *), iqsign_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, logical *, integer *), reordr_(integer *, 
	    integer *, integer *, integer *), sitord_(integer *, integer *, 
	    integer *);
    integer ipossi, opvert;



/*     initialize */

    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;
    --wi;
    --zi;
    --yi;
    --xi;
    --ifl;
    icon -= 9;

    /* Function Body */
    *idmax = 0;

/*     test all tetra with ifl(tetra)=1 for the Regular/Delaunay property */

    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ifl[i__] != 1) {
	    goto L200;
	}
	a = icon[(i__ << 3) + 5];
	b = icon[(i__ << 3) + 6];
	c__ = icon[(i__ << 3) + 7];
	d__ = icon[(i__ << 3) + 8];
	for (j = 1; j <= 4; ++j) {
	    adj = icon[j + (i__ << 3)];
	    if (adj == 0) {
		s_stop("2610", (ftnlen)4);
	    }
	    if (ifl[adj] != 1) {
		goto L100;
	    }
	    if (adj > i__) {
		goto L100;
	    }
	    for (k = 1; k <= 8; ++k) {
		ikon[k - 1] = icon[k + (i__ << 3)];
/* L50: */
	    }
	    isite = icon[j + 4 + (i__ << 3)];
	    sitord_(ikon, &isite, &c__1);

	    site1 = ikon[5];
	    site2 = ikon[6];
	    site3 = ikon[7];
	    reordr_(&icon[9], &site1, &site2, &adj);
	    if (icon[(adj << 3) + 7] != site3) {
		s_stop("2620", (ftnlen)4);
	    }
	    if (icon[(adj << 3) + 4] != i__) {
		s_stop("2630", (ftnlen)4);
	    }

	    opvert = icon[(adj << 3) + 8];
	    ctrad_(&xi[1], &yi[1], &zi[1], &wi[1], &xctr, &yctr, &zctr, &
		    site1, &site2, &site3, &opvert, epz, delaun, &ipossi);
	    if (ipossi == 1) {
		goto L70;
	    }
	    bisphr_(&xi[1], &yi[1], &zi[1], &wi[1], &isite, &site1, epz, &
		    xctr, &yctr, &zctr, &tdist, delaun, &ipossi);
	    if (ipossi == 1) {
		goto L70;
	    }
	    if (tdist <= 0.) {
		goto L100;
	    }
	    goto L80;
L70:
	    iqsign_(&x[1], &y[1], &z__[1], &w[1], &x2[1], &y2[1], &z2[1], &w2[
		    1], &a, &b, &c__, &d__, &opvert, mhalf, mfull, &isclp[1], 
		    &isclw[1], &isclr[1], delaun, &itide);
	    if (itide >= 0) {
		goto L100;
	    }
L80:
	    ++(*idmax);
L100:
	    ;
	}
L200:
	;
    }

    return 0;
} /* delchk_ */

/* REVTET */

/* This subroutine will compress data structure in order to save space */
/* by eliminating artificial and discarded tetrahedra */

/* Subroutine */ int revtet_(integer *tetra, integer *tetru, integer *icon, 
	integer *nv, integer *is, integer *ifl, logical *flphis)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, ii, ielm;



/*     identify true tetrahedra */

    /* Parameter adjustments */
    --ifl;
    --is;
    icon -= 9;

    /* Function Body */
    *tetru = 0;
    ielm = 0;
    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (icon[(i__ << 3) + 5] <= 0 || icon[(i__ << 3) + 6] <= 0 || icon[(
		i__ << 3) + 7] <= 0 || icon[(i__ << 3) + 8] <= 0) {
	    ++ielm;
	    ifl[i__] = 0;
	} else if (icon[(i__ << 3) + 5] <= 8 || icon[(i__ << 3) + 6] <= 8 || 
		icon[(i__ << 3) + 7] <= 8 || icon[(i__ << 3) + 8] <= 8) {
	    ifl[i__] = 0;
	} else {
	    ++(*tetru);
	    ifl[i__] = 1;
	}
/* L100: */
    }
    if (*tetru == 0) {
	goto L1000;
    }

/*     zero out nonexistent tetrahedra in icon */

    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ifl[i__] == 0) {
	    goto L300;
	}
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (i__ << 3)] <= 0 || icon[j + (i__ << 3)] > *tetra) {
		s_stop("2710", (ftnlen)4);
	    }
	    if (ifl[icon[j + (i__ << 3)]] == 0) {
		icon[j + (i__ << 3)] = 0;
	    }
/* L200: */
	}
L300:
	;
    }

/*     compress icon */

    ii = 0;
    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ifl[i__] == 0) {
	    goto L500;
	}
	++ii;
	ifl[i__] = ii;
	for (j = 1; j <= 8; ++j) {
	    icon[j + (ii << 3)] = icon[j + (i__ << 3)];
/* L400: */
	}
L500:
	;
    }

/*     update icon for tetrahedra and is for vertices */

    i__1 = *nv;
    for (i__ = 9; i__ <= i__1; ++i__) {
	if (is[i__] > 0) {
	    is[i__ - 8] = 1;
	} else {
	    is[i__ - 8] = is[i__];
	}
/* L550: */
    }
    i__1 = *tetru;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (i__ << 3)] == 0) {
		goto L600;
	    }
	    icon[j + (i__ << 3)] = ifl[icon[j + (i__ << 3)]];
L600:
	    ;
	}
	for (j = 5; j <= 8; ++j) {
	    icon[j + (i__ << 3)] += -8;
	    if (is[icon[j + (i__ << 3)]] <= 0) {
		s_stop("2720", (ftnlen)4);
	    }
	    is[icon[j + (i__ << 3)]] = i__;
/* L700: */
	}
/* L800: */
    }

L1000:
    if (! (*flphis)) {
	*tetra -= ielm;
    }

    return 0;
} /* revtet_ */

/* RUVTET */

/* This subroutine will compress data structure in order to save space */
/* by eliminating discarded tetrahedra while keeping artificial ones */

/* Subroutine */ int ruvtet_(integer *tetra, integer *tetru, integer *icon, 
	integer *is, integer *ifl)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, ii, ielm;



/*     identify true tetrahedra */

    /* Parameter adjustments */
    --ifl;
    --is;
    icon -= 9;

    /* Function Body */
    *tetru = 0;
    ielm = 0;
    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (icon[(i__ << 3) + 5] <= 0 || icon[(i__ << 3) + 6] <= 0 || icon[(
		i__ << 3) + 7] <= 0 || icon[(i__ << 3) + 8] <= 0) {
	    ++ielm;
	    ifl[i__] = 0;
	} else if (icon[(i__ << 3) + 5] <= 8 || icon[(i__ << 3) + 6] <= 8 || 
		icon[(i__ << 3) + 7] <= 8 || icon[(i__ << 3) + 8] <= 8) {
	    ifl[i__] = 1;
	} else {
	    ++(*tetru);
	    ifl[i__] = 1;
	}
/* L100: */
    }

/*     compress icon */

    ii = 0;
    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ifl[i__] == 0) {
	    goto L500;
	}
	++ii;
	ifl[i__] = ii;
	for (j = 1; j <= 8; ++j) {
	    icon[j + (ii << 3)] = icon[j + (i__ << 3)];
/* L400: */
	}
L500:
	;
    }

    *tetra -= ielm;

/*     update icon for tetrahedra and is for vertices */

    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 4; ++j) {
	    if (icon[j + (i__ << 3)] == 0) {
		goto L600;
	    }
	    icon[j + (i__ << 3)] = ifl[icon[j + (i__ << 3)]];
L600:
	    ;
	}
	for (j = 5; j <= 8; ++j) {
	    if (is[icon[j + (i__ << 3)]] <= 0) {
		s_stop("2730", (ftnlen)4);
	    }
	    is[icon[j + (i__ << 3)]] = i__;
/* L700: */
	}
/* L800: */
    }

    return 0;
} /* ruvtet_ */

/* CONSIS */

/*     subroutine consis to - */

/*     test consistency of diagram */

/*     May 1, 1989 */

/* Subroutine */ int consis_(integer *icon, integer *is, integer *ifl, 
	integer *n, integer *ivnxt)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, ikon[8]	/* was [8][1] */, indx, site0, site1, site2, 
	    site3, isini, isone, iscur, islst;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *);



/*     test initial tetrahedron for each site */

    /* Parameter adjustments */
    --ifl;
    --is;
    icon -= 9;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iscur = is[i__];
	if (iscur <= 0) {
	    goto L50;
	}
	if (icon[(iscur << 3) + 5] != i__ && icon[(iscur << 3) + 6] != i__ && 
		icon[(iscur << 3) + 7] != i__ && icon[(iscur << 3) + 8] != 
		i__) {
	    s_stop("2810", (ftnlen)4);
	}
L50:
	;
    }

/*     initialize */

    isone = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (is[i__] > 0) {
	    goto L80;
	}
/* L60: */
    }
    s_stop("2820", (ftnlen)4);
L80:
    islst = is[i__];
    isini = islst;

    i__1 = *ivnxt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ifl[i__] = 0;
/* L100: */
    }

    ifl[isini] = 1;
    indx = 1;
    iscur = icon[(isini << 3) + 1];
    if (iscur == 0) {
	goto L500;
    }
    site0 = icon[(isini << 3) + 5];
    site1 = icon[(isini << 3) + 6];
    site2 = icon[(isini << 3) + 7];
    site3 = icon[(isini << 3) + 8];

/*     reorder iscur relative to site1 and site2, and test */

L200:
    if (site0 == site1 || site0 == site2 || site0 == site3 || site1 == site2 
	    || site1 == site3 || site2 == site3) {
	s_stop("2830", (ftnlen)4);
    }
    reordr_(&icon[9], &site1, &site2, &iscur);
    if (icon[(iscur << 3) + 7] != site3) {
	s_stop("2840", (ftnlen)4);
    }
    if (icon[(iscur << 3) + 4] != islst) {
	s_stop("2850", (ftnlen)4);
    }
    if (icon[(iscur << 3) + 8] == site0) {
	s_stop("2855", (ftnlen)4);
    }
    ifl[iscur] = 1;

/*     obtain next tetrahedron */

    islst = iscur;
    indx = 1;
    iscur = icon[(islst << 3) + 1];
    if (iscur == 0) {
	goto L500;
    }
    site0 = icon[(islst << 3) + 5];
    site1 = icon[(islst << 3) + 6];
    site2 = icon[(islst << 3) + 7];
    site3 = icon[(islst << 3) + 8];
    if (ifl[iscur] != 1) {
	goto L200;
    }

/*     reorder iscur relative to site1 and site2, and test */

L300:
    if (site0 == site1 || site0 == site2 || site0 == site3 || site1 == site2 
	    || site1 == site3 || site2 == site3) {
	s_stop("2860", (ftnlen)4);
    }
    for (i__ = 1; i__ <= 8; ++i__) {
	ikon[i__ - 1] = icon[i__ + (iscur << 3)];
/* L400: */
    }
    reordr_(ikon, &site1, &site2, &isone);
    if (ikon[6] != site3) {
	s_stop("2865", (ftnlen)4);
    }
    if (ikon[3] != islst) {
	s_stop("2870", (ftnlen)4);
    }
    if (ikon[7] == site0) {
	s_stop("2875", (ftnlen)4);
    }

/*     obtain next tetrahedron */

L500:
    if (indx == 1) {
	indx = 2;
	iscur = icon[(islst << 3) + 2];
	if (iscur == 0) {
	    goto L500;
	}
	site0 = icon[(islst << 3) + 6];
	site1 = icon[(islst << 3) + 5];
	site2 = icon[(islst << 3) + 8];
	site3 = icon[(islst << 3) + 7];
	if (ifl[iscur] != 1) {
	    goto L200;
	}
	goto L300;
    } else if (indx == 2) {
	indx = 3;
	iscur = icon[(islst << 3) + 3];
	if (iscur == 0) {
	    goto L500;
	}
	site0 = icon[(islst << 3) + 7];
	site1 = icon[(islst << 3) + 5];
	site2 = icon[(islst << 3) + 6];
	site3 = icon[(islst << 3) + 8];
	if (ifl[iscur] != 1) {
	    goto L200;
	}
	goto L300;
    } else if (indx == 3) {
	if (islst != isini) {
	    iscur = islst;
	    islst = icon[(iscur << 3) + 4];
	    if (islst == 0) {
		s_stop("2880", (ftnlen)4);
	    }
	    if (icon[(islst << 3) + 1] == iscur) {
		indx = 1;
	    } else if (icon[(islst << 3) + 2] == iscur) {
		indx = 2;
	    } else if (icon[(islst << 3) + 3] == iscur) {
		indx = 3;
	    } else if (icon[(islst << 3) + 4] == iscur) {
		indx = 4;
	    } else {
		s_stop("2885", (ftnlen)4);
	    }
	    goto L500;
	} else {
	    indx = 4;
	    iscur = icon[(islst << 3) + 4];
	    if (iscur == 0) {
		goto L500;
	    }
	    site0 = icon[(islst << 3) + 8];
	    site1 = icon[(islst << 3) + 5];
	    site2 = icon[(islst << 3) + 7];
	    site3 = icon[(islst << 3) + 6];
	    if (ifl[iscur] != 1) {
		goto L200;
	    }
	    goto L300;
	}
    }
    if (islst != isini) {
	s_stop("2890", (ftnlen)4);
    }

/*     write (*,*) ' ' */
/*     write (*,*) '**************************************' */
/*     write (*,*) 'consistency check satisfied' */
/*     write (*,*) '**************************************' */
/*     write (*,*) ' ' */

    return 0;
} /* consis_ */

/* ORIENT */

/*     This subroutine will test the orientation of the tetrahedra */

/* Subroutine */ int orient_(integer *tetra, integer *icon, integer *ifl, 
	doublereal *xi, doublereal *yi, doublereal *zi, integer *x, integer *
	y, integer *z__, integer *x2, integer *y2, integer *z2, integer *
	idmin, integer *mhalf, integer *mfull, integer *isclp, doublereal *
	epz)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer a, b, c__, d__, i__, iside;
    extern /* Subroutine */ int irsign_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *);



/*     test all tetrahedra with ifl equal to 1 */

    /* Parameter adjustments */
    --isclp;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;
    --zi;
    --yi;
    --xi;
    --ifl;
    icon -= 9;

    /* Function Body */
    *idmin = 0;
    i__1 = *tetra;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ifl[i__] != 1) {
	    goto L200;
	}
	a = icon[(i__ << 3) + 5];
	b = icon[(i__ << 3) + 6];
	c__ = icon[(i__ << 3) + 7];
	d__ = icon[(i__ << 3) + 8];
	irsign_(&xi[1], &yi[1], &zi[1], &x[1], &y[1], &z__[1], &x2[1], &y2[1],
		 &z2[1], &d__, &a, &b, &c__, mhalf, mfull, &isclp[1], epz, &
		iside);
	if (iside <= 0) {
	    ++(*idmin);
	}
L200:
	;
    }

    return 0;
} /* orient_ */

/* REORDR */

/*     subroutine reordr to - */

/*     reorder icon(i,iscur), i = 1, ..., 8, so that site1 equals */
/*     icon(5,iscur) and site2 equals icon(6,iscur) */

/*     July 22, 1988 */

/* Subroutine */ int reordr_(integer *icon, integer *site1, integer *site2, 
	integer *iscur)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer itemp;



    /* Parameter adjustments */
    icon -= 9;

    /* Function Body */
    if (icon[(*iscur << 3) + 5] == *site1) {
	goto L200;
    }
    if (icon[(*iscur << 3) + 6] == *site1) {
	itemp = icon[(*iscur << 3) + 1];
	icon[(*iscur << 3) + 1] = icon[(*iscur << 3) + 2];
	icon[(*iscur << 3) + 2] = icon[(*iscur << 3) + 4];
	icon[(*iscur << 3) + 4] = itemp;
	itemp = icon[(*iscur << 3) + 5];
	icon[(*iscur << 3) + 5] = icon[(*iscur << 3) + 6];
	icon[(*iscur << 3) + 6] = icon[(*iscur << 3) + 8];
	icon[(*iscur << 3) + 8] = itemp;
    } else if (icon[(*iscur << 3) + 7] == *site1) {
	itemp = icon[(*iscur << 3) + 1];
	icon[(*iscur << 3) + 1] = icon[(*iscur << 3) + 3];
	icon[(*iscur << 3) + 3] = icon[(*iscur << 3) + 2];
	icon[(*iscur << 3) + 2] = itemp;
	itemp = icon[(*iscur << 3) + 5];
	icon[(*iscur << 3) + 5] = icon[(*iscur << 3) + 7];
	icon[(*iscur << 3) + 7] = icon[(*iscur << 3) + 6];
	icon[(*iscur << 3) + 6] = itemp;
    } else if (icon[(*iscur << 3) + 8] == *site1) {
	itemp = icon[(*iscur << 3) + 1];
	icon[(*iscur << 3) + 1] = icon[(*iscur << 3) + 4];
	icon[(*iscur << 3) + 4] = icon[(*iscur << 3) + 3];
	icon[(*iscur << 3) + 3] = itemp;
	itemp = icon[(*iscur << 3) + 5];
	icon[(*iscur << 3) + 5] = icon[(*iscur << 3) + 8];
	icon[(*iscur << 3) + 8] = icon[(*iscur << 3) + 7];
	icon[(*iscur << 3) + 7] = itemp;
    } else {
	s_stop("2910", (ftnlen)4);
    }
L200:

    if (icon[(*iscur << 3) + 6] == *site2) {
	goto L300;
    }
    if (icon[(*iscur << 3) + 7] == *site2) {
	itemp = icon[(*iscur << 3) + 2];
	icon[(*iscur << 3) + 2] = icon[(*iscur << 3) + 3];
	icon[(*iscur << 3) + 3] = icon[(*iscur << 3) + 4];
	icon[(*iscur << 3) + 4] = itemp;
	itemp = icon[(*iscur << 3) + 6];
	icon[(*iscur << 3) + 6] = icon[(*iscur << 3) + 7];
	icon[(*iscur << 3) + 7] = icon[(*iscur << 3) + 8];
	icon[(*iscur << 3) + 8] = itemp;
    } else if (icon[(*iscur << 3) + 8] == *site2) {
	itemp = icon[(*iscur << 3) + 2];
	icon[(*iscur << 3) + 2] = icon[(*iscur << 3) + 4];
	icon[(*iscur << 3) + 4] = icon[(*iscur << 3) + 3];
	icon[(*iscur << 3) + 3] = itemp;
	itemp = icon[(*iscur << 3) + 6];
	icon[(*iscur << 3) + 6] = icon[(*iscur << 3) + 8];
	icon[(*iscur << 3) + 8] = icon[(*iscur << 3) + 7];
	icon[(*iscur << 3) + 7] = itemp;
    } else {
	s_stop("2920", (ftnlen)4);
    }
L300:

    return 0;
} /* reordr_ */

/* SITORD */

/*     subroutine sitord to - */

/*     reorder icon(i,iscur), i = 1, ..., 8, so that site1 equals */
/*     icon(5,iscur) */

/*     July 22, 1988 */

/* Subroutine */ int sitord_(integer *icon, integer *site1, integer *iscur)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer itemp;



    /* Parameter adjustments */
    icon -= 9;

    /* Function Body */
    if (icon[(*iscur << 3) + 5] == *site1) {
	goto L200;
    }
    if (icon[(*iscur << 3) + 6] == *site1) {
	itemp = icon[(*iscur << 3) + 1];
	icon[(*iscur << 3) + 1] = icon[(*iscur << 3) + 2];
	icon[(*iscur << 3) + 2] = icon[(*iscur << 3) + 4];
	icon[(*iscur << 3) + 4] = itemp;
	itemp = icon[(*iscur << 3) + 5];
	icon[(*iscur << 3) + 5] = icon[(*iscur << 3) + 6];
	icon[(*iscur << 3) + 6] = icon[(*iscur << 3) + 8];
	icon[(*iscur << 3) + 8] = itemp;
    } else if (icon[(*iscur << 3) + 7] == *site1) {
	itemp = icon[(*iscur << 3) + 1];
	icon[(*iscur << 3) + 1] = icon[(*iscur << 3) + 3];
	icon[(*iscur << 3) + 3] = icon[(*iscur << 3) + 2];
	icon[(*iscur << 3) + 2] = itemp;
	itemp = icon[(*iscur << 3) + 5];
	icon[(*iscur << 3) + 5] = icon[(*iscur << 3) + 7];
	icon[(*iscur << 3) + 7] = icon[(*iscur << 3) + 6];
	icon[(*iscur << 3) + 6] = itemp;
    } else if (icon[(*iscur << 3) + 8] == *site1) {
	itemp = icon[(*iscur << 3) + 1];
	icon[(*iscur << 3) + 1] = icon[(*iscur << 3) + 4];
	icon[(*iscur << 3) + 4] = icon[(*iscur << 3) + 3];
	icon[(*iscur << 3) + 3] = itemp;
	itemp = icon[(*iscur << 3) + 5];
	icon[(*iscur << 3) + 5] = icon[(*iscur << 3) + 8];
	icon[(*iscur << 3) + 8] = icon[(*iscur << 3) + 7];
	icon[(*iscur << 3) + 7] = itemp;
    } else {
	s_stop("3010", (ftnlen)4);
    }
L200:
    return 0;
} /* sitord_ */

/* VRTORD */

/* This routine will order vertices a, b, c, d of a tetrahedron */
/* so that a>b, b>c, b>d. Data structure is updated */

/* Subroutine */ int vrtord_(integer *icon, integer *curr, integer *a, 
	integer *b, integer *c__, integer *d__)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer it;
    extern /* Subroutine */ int reordr_(integer *, integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    icon -= 9;

    /* Function Body */
    if (*a < *b) {
	it = *a;
	*a = *b;
	*b = it;
    }
    if (*a < *c__) {
	it = *a;
	*a = *c__;
	*c__ = it;
    }
    if (*a < *d__) {
	it = *a;
	*a = *d__;
	*d__ = it;
    }
    if (*b < *c__) {
	*b = *c__;
    }
    if (*b < *d__) {
	*b = *d__;
    }
    reordr_(&icon[9], a, b, curr);
    *c__ = icon[(*curr << 3) + 7];
    *d__ = icon[(*curr << 3) + 8];
    if (*b > *a || *c__ > *b || *d__ > *b) {
	s_stop("3110", (ftnlen)4);
    }

    return 0;
} /* vrtord_ */

/* VRTARR */

/* This routine will arrange vertices b, c, d of a tetrahedron */
/* so that b>c, b>d. Data structure is not updated */

/* Subroutine */ int vrtarr_(integer *i2, integer *i3, integer *i4, integer *
	b, integer *c__, integer *d__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer ix;



    *b = *i2;
    *c__ = *i3;
    *d__ = *i4;
/* Computing MAX */
    i__1 = max(*b,*c__);
    ix = max(i__1,*d__);
    if (*b == ix) {
	goto L100;
    }
    if (*c__ == ix) {
	ix = *b;
	*b = *c__;
	*c__ = *d__;
	*d__ = ix;
    } else {
	ix = *b;
	*b = *d__;
	*d__ = *c__;
	*c__ = ix;
    }
L100:
    if (*c__ > *b || *d__ > *b) {
	s_stop("3210", (ftnlen)4);
    }

    return 0;
} /* vrtarr_ */

/* RDMORD */

/*     subroutine to reorder randomly integers from 1 to n */


/* Subroutine */ int rdmord_(real *wr, integer *io, integer *n, integer *isu, 
	integer *jsu, integer *ksu, integer *nsu)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, iu, ju, ku, nu;
    extern /* Subroutine */ int umi_(real *, integer *, integer *, integer *, 
	    integer *);
    real rumi;
    extern /* Subroutine */ int mzrans_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), trsort_(
	    real *, integer *, integer *);



/*     initialize */

/*     isu = 521288629 */
/*     jsu = 362436069 */
/*     ksu = 16163801 */
/*     nsu = 131199299 */

    /* Parameter adjustments */
    --io;
    --wr;

    /* Function Body */
    mzrans_(isu, jsu, ksu, nsu, &iu, &ju, &ku, &nu);

/*     get numbers */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	umi_(&rumi, &iu, &ju, &ku, &nu);
	wr[i__] = rumi;
	io[i__] = i__;
/* L10: */
    }

/*     sort in increasing order thus obtaining random order */
/*     of integers from 1 to n */

/*     OPEN(21,FILE='umi.dat') */
/*     WRITE(21,*)(WR(I),I=1,N) */

    trsort_(&wr[1], &io[1], n);

    return 0;
} /* rdmord_ */

/* TRSORT */

/*     subroutine trsort to - */

/*     sort an array of real numbers in increasing order */
/*     in O(k log k) time */

/*     January 15, 1988 */

/* Subroutine */ int trsort_(real *var, integer *ji, integer *klen)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, k, jj, iis, iast;



/*     create initial tree in decreasing order */

    /* Parameter adjustments */
    --ji;
    --var;

    /* Function Body */
    iast = *klen;
    i__1 = *klen;
    for (k = 1; k <= i__1; ++k) {
	i__ = k;
L50:

/*     check if current node is as small as father */

	if (i__ == 1) {
	    goto L100;
	}
	if (var[ji[i__]] <= var[ji[i__ / 2]]) {
	    goto L100;
	}
	jj = ji[i__];
	ji[i__] = ji[i__ / 2];
	ji[i__ / 2] = jj;
	i__ /= 2;
	goto L50;
L100:
	;
    }
    if (iast == 1) {
	goto L160;
    }

/*     sort by shrinking tree: last element is moved to the */
/*     first position */

L102:
    i__ = 1;
    jj = ji[1];
    ji[1] = ji[iast];
    ji[iast] = jj;
    if (iast == 2) {
	goto L160;
    }
    --iast;

L105:
    iis = i__ << 1;

/*     check which sons exist */

    if ((i__1 = iis - iast) < 0) {
	goto L110;
    } else if (i__1 == 0) {
	goto L140;
    } else {
	goto L150;
    }

/*     both sons exist */

L110:
    if (var[ji[i__]] < var[ji[iis]]) {
	goto L120;
    }
    if (var[ji[i__]] >= var[ji[iis + 1]]) {
	goto L150;
    }
    goto L125;

/*     check which son to be switched */

L120:
    if (var[ji[iis]] >= var[ji[iis + 1]]) {
	goto L130;
    }

/*     adjust to switch with right son */

L125:
    ++iis;

/*     switch */

L130:
    jj = ji[i__];
    ji[i__] = ji[iis];
    ji[iis] = jj;
    i__ = iis;
    goto L105;

/*     only left son exists */

L140:
    if (var[ji[i__]] >= var[ji[iis]]) {
	goto L150;
    }
    jj = ji[i__];
    ji[i__] = ji[iis];
    ji[iis] = jj;

/*     no more switching needed */

L150:
    goto L102;

/*     sorting is finished */

L160:

    return 0;
} /* trsort_ */

/* UMI */

/*     subroutine umi to - */

/*     generate random numbers */

/* Subroutine */ int umi_(real *rumi, integer *i__, integer *j, integer *k, 
	integer *n)
{
    integer mzrn;
    extern /* Subroutine */ int mzran_(integer *, integer *, integer *, 
	    integer *, integer *);



    mzran_(i__, j, k, n, &mzrn);
    *rumi = mzrn * (float)2.328306e-10 + (float).5;

    return 0;
} /* umi_ */

/* MZRAN */

/*     subroutine mzran to - */

/*     do computations in order to generate random numbers */

/* Subroutine */ int mzran_(integer *i__, integer *j, integer *k, integer *n, 
	integer *mzrn)
{


    *mzrn = *i__ - *k;
    if (*mzrn < 0) {
	*mzrn += 2147483579;
    }
    *i__ = *j;
    *j = *k;
    *k = *mzrn;
    *n = *n * 69069 + 1013904243;
    *mzrn += *n;

    return 0;
} /* mzran_ */

/* MZRANS */

/*     subroutine mzrans to - */

/*     initialize in order to generate random numbers */

/* Subroutine */ int mzrans_(integer *is, integer *js, integer *ks, integer *
	ns, integer *i__, integer *j, integer *k, integer *n)
{

/*     save is,js,ks,ns */
/*     data is,js,ks,ns/521288629,362436069,16163801,1131199299/ */


    *i__ = abs(*is) + 1;
    *j = abs(*js) + 1;
    *k = abs(*ks) + 1;
    *n = *ns;

    return 0;
} /* mzrans_ */

/* DSTNCE */

/* This subroutine will compute the distance from a point to a facet of */
/* a tetrahedron. */

/* Subroutine */ int dstnce_(doublereal *x, doublereal *y, doublereal *z__, 
	integer *p, integer *q, integer *r__, doublereal *epz, integer *k, 
	doublereal *dist, integer *ipossi)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal dst1, dst2, dst3, dlen, dmax__, dlun, dstp, dotx, doty, dotz, 
	    xvec1, yvec1, zvec1, xvec2, yvec2, zvec2, xvec3, yvec3, zvec3, 
	    xvecp, yvecp, zvecp;



    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    *ipossi = 0;
    xvec1 = x[*q] - x[*p];
    yvec1 = y[*q] - y[*p];
    zvec1 = z__[*q] - z__[*p];
    xvec2 = x[*r__] - x[*p];
    yvec2 = y[*r__] - y[*p];
    zvec2 = z__[*r__] - z__[*p];
    xvec3 = x[*q] - x[*r__];
    yvec3 = y[*q] - y[*r__];
    zvec3 = z__[*q] - z__[*r__];
/* Computing 2nd power */
    d__1 = xvec1;
/* Computing 2nd power */
    d__2 = yvec1;
/* Computing 2nd power */
    d__3 = zvec1;
    dst1 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
    d__1 = xvec2;
/* Computing 2nd power */
    d__2 = yvec2;
/* Computing 2nd power */
    d__3 = zvec2;
    dst2 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
    d__1 = xvec3;
/* Computing 2nd power */
    d__2 = yvec3;
/* Computing 2nd power */
    d__3 = zvec3;
    dst3 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    if (dst1 < *epz || dst2 < *epz || dst3 < *epz) {
	*ipossi = 1;
	goto L1000;
    }
/* Computing MAX */
    d__1 = max(dst1,dst2);
    dmax__ = max(d__1,dst3);

    dotx = yvec1 * zvec2 - yvec2 * zvec1;
    doty = -xvec1 * zvec2 + xvec2 * zvec1;
    dotz = xvec1 * yvec2 - xvec2 * yvec1;
/* Computing 2nd power */
    d__1 = dotx;
/* Computing 2nd power */
    d__2 = doty;
/* Computing 2nd power */
    d__3 = dotz;
    dlen = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    if (dlen < *epz || dlen / dmax__ < *epz) {
	*ipossi = 1;
	goto L1000;
    }

    xvecp = x[*k] - x[*p];
    yvecp = y[*k] - y[*p];
    zvecp = z__[*k] - z__[*p];
/* Computing 2nd power */
    d__1 = xvecp;
/* Computing 2nd power */
    d__2 = yvecp;
/* Computing 2nd power */
    d__3 = zvecp;
    dstp = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    if (dstp < *epz) {
	*ipossi = 1;
	goto L1000;
    }

    dlun = dstp * dmax__;
    dlun = max(dlen,dlun);
    *dist = (xvecp * dotx + yvecp * doty + zvecp * dotz) / dlun;
    if (*dist > -(*epz) && *dist < *epz) {
	*ipossi = 1;
    }

L1000:
    return 0;
} /* dstnce_ */

/* CTRAD */

/* This subroutine will compute the orthogonal center of a tetrahedron */

/* Subroutine */ int ctrad_(doublereal *x, doublereal *y, doublereal *z__, 
	doublereal *w, doublereal *xctr, doublereal *yctr, doublereal *zctr, 
	integer *a, integer *b, integer *c__, integer *d__, doublereal *epz, 
	logical *delaun, integer *ipossi)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal xe, ye, ze, xl, xm, ym, zm, xn, yn, zn, xq, yq, xu, yu, zu, xv,
	     yv, zv, xw, yw, zw, zq, yl, zl, xt, yt, zt, dmax__, norm, denom, 
	    normu, normv, lambda;



/*     initialize */

    /* Parameter adjustments */
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    *ipossi = 0;

/*     find midpoints of edges ac and ab */

    xm = (x[*a] + x[*c__]) / 2.;
    ym = (y[*a] + y[*c__]) / 2.;
    zm = (z__[*a] + z__[*c__]) / 2.;

    xn = (x[*a] + x[*b]) / 2.;
    yn = (y[*a] + y[*b]) / 2.;
    zn = (z__[*a] + z__[*b]) / 2.;

/*     compute edge vectors u and v for edges ac and ab */

    xu = x[*c__] - x[*a];
    yu = y[*c__] - y[*a];
    zu = z__[*c__] - z__[*a];

    xv = x[*b] - x[*a];
    yv = y[*b] - y[*a];
    zv = z__[*b] - z__[*a];

/*     compute lengths of u and v */

/* Computing 2nd power */
    d__1 = xu;
/* Computing 2nd power */
    d__2 = yu;
/* Computing 2nd power */
    d__3 = zu;
    normu = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* Computing 2nd power */
    d__1 = xv;
/* Computing 2nd power */
    d__2 = yv;
/* Computing 2nd power */
    d__3 = zv;
    normv = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    if (normu < *epz || normv < *epz) {
	*ipossi = 1;
	goto L1000;
    }
    dmax__ = max(normu,normv);

/*     find perpendicular to facet abc of tetrahedron */

    xw = yu * zv - zu * yv;
    yw = -xu * zv + zu * xv;
    zw = xu * yv - yu * xv;

/*     test whether edges ac, ab are colinear */

/* Computing 2nd power */
    d__1 = xw;
/* Computing 2nd power */
    d__2 = yw;
/* Computing 2nd power */
    d__3 = zw;
    norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3) / dmax__;
    if (norm < *epz) {
	*ipossi = 1;
	goto L1000;
    }
    xw /= normu;
    yw /= normu;
    zw /= normu;

/*     normalize u and v */

    xu /= normu;
    yu /= normu;
    zu /= normu;
    xv /= normv;
    yv /= normv;
    zv /= normv;

/*     compute orthogonal center of edge ac */

    if (! (*delaun)) {
	lambda = (w[*a] - w[*c__]) / normu / 2.;
	xm += lambda * xu;
	ym += lambda * yu;
	zm += lambda * zu;
    }

/*     compute orthogonal center of edge ab */

    if (! (*delaun)) {
	lambda = (w[*a] - w[*b]) / normv / 2.;
	xn += lambda * xv;
	yn += lambda * yv;
	zn += lambda * zv;
    }

/*     find perpendicular to edge v in plane that contains facet abc */

    xq = yw * zv - zw * yv;
    yq = -xw * zv + zw * xv;
    zq = xw * yv - yw * xv;
/* Computing 2nd power */
    d__1 = xq;
/* Computing 2nd power */
    d__2 = yq;
/* Computing 2nd power */
    d__3 = zq;
    norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    if (norm < *epz) {
	*ipossi = 1;
	goto L1000;
    }

/*     compute orthogonal center of facet abc */

    denom = xu * xq + yu * yq + zu * zq;
    if (denom > -(*epz) && denom < *epz) {
	*ipossi = 1;
	goto L1000;
    }
    lambda = (xu * (xm - xn) + yu * (ym - yn) + zu * (zm - zn)) / denom;

    xe = xn + lambda * xq;
    ye = yn + lambda * yq;
    ze = zn + lambda * zq;

/*     compute edge vector t for edge ad */

    xl = (x[*a] + x[*d__]) / 2.;
    yl = (y[*a] + y[*d__]) / 2.;
    zl = (z__[*a] + z__[*d__]) / 2.;

    xt = x[*d__] - x[*a];
    yt = y[*d__] - y[*a];
    zt = z__[*d__] - z__[*a];
/* Computing 2nd power */
    d__1 = xt;
/* Computing 2nd power */
    d__2 = yt;
/* Computing 2nd power */
    d__3 = zt;
    norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    if (norm < *epz) {
	*ipossi = 1;
	goto L1000;
    }
    xt /= norm;
    yt /= norm;
    zt /= norm;

/*     compute orthogonal center of edge ad */

    if (! (*delaun)) {
	lambda = (w[*a] - w[*d__]) / norm / 2.;
	xl += lambda * xt;
	yl += lambda * yt;
	zl += lambda * zt;
    }

/*     compute orthogonal center of tetrahedron */

    denom = xt * xw + yt * yw + zt * zw;
    if (denom > -(*epz) && denom < *epz) {
	*ipossi = 1;
	goto L1000;
    }
    lambda = (xt * (xl - xe) + yt * (yl - ye) + zt * (zl - ze)) / denom;

    *xctr = xe + lambda * xw;
    *yctr = ye + lambda * yw;
    *zctr = ze + lambda * zw;

L1000:
    return 0;
} /* ctrad_ */

/* BISPHR */

/* This subroutine will compute the distance from a point */
/* (xctr,yctr,zctr) to the chordale plane between two points */
/* opvert and ivrt */

/* Subroutine */ int bisphr_(doublereal *x, doublereal *y, doublereal *z__, 
	doublereal *w, integer *opvert, integer *ivrt, doublereal *epz, 
	doublereal *xctr, doublereal *yctr, doublereal *zctr, doublereal *
	tdist, logical *delaun, integer *ipossi)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal xd, yd, zd, xm, ym, zm, xu, yu, zu, xu2, yu2, zu2, dif, dmax__,
	     norm, wambda;



/*     find midpoint of edge from opvert to ivrt */

    /* Parameter adjustments */
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    xm = (x[*opvert] + x[*ivrt]) / 2.;
    ym = (y[*opvert] + y[*ivrt]) / 2.;
    zm = (z__[*opvert] + z__[*ivrt]) / 2.;

/*     find vector from ivrt to opvert */

    xu = x[*opvert] - x[*ivrt];
    yu = y[*opvert] - y[*ivrt];
    zu = z__[*opvert] - z__[*ivrt];

/* Computing 2nd power */
    d__1 = xu;
/* Computing 2nd power */
    d__2 = yu;
/* Computing 2nd power */
    d__3 = zu;
    norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    if (norm < *epz) {
	*ipossi = 1;
	goto L1000;
    }
    xu2 = xu / norm;
    yu2 = yu / norm;
    zu2 = zu / norm;

/*     compute orthogonal center of edge ivrt-opvert */

    if (! (*delaun)) {
	wambda = (w[*ivrt] - w[*opvert]) / norm / 2.;
	xm += wambda * xu2;
	ym += wambda * yu2;
	zm += wambda * zu2;
    }

/*     compute distance */

    xd = *xctr - xm;
    yd = *yctr - ym;
    zd = *zctr - zm;
/* Computing 2nd power */
    d__1 = xd;
/* Computing 2nd power */
    d__2 = yd;
/* Computing 2nd power */
    d__3 = zd;
    dif = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    dmax__ = max(norm,dif);
    *tdist = (xd * xu + yd * yu + zd * zu) / dmax__;
    if (*tdist > -(*epz) && *tdist < *epz) {
	*ipossi = 1;
    }

L1000:
    return 0;
} /* bisphr_ */

/* IPSIGN */

/*     subroutine for determining position of point ifou with respect */
/*     to plane that contains points ifir, isec, ithi */
/*     if positive then ifou is on positive side of plane */
/*     if negative then ifou is on negative side of plane */
/*     if zero then ifou is in plane */

/* Subroutine */ int ipsign_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *ifir, integer *isec, integer *
	ithi, integer *ifou, integer *mhalf, integer *mfull, integer *isclp, 
	integer *ipout)
{
    integer io[30], iu[30], iv[30], iw[30], ix2[30], iy2[30], iz2[30], ix3[30]
	    , iy3[30], iz3[30], ix4[30], iy4[30], iz4[30], iko, ixf[30], iyf[
	    30], izf[30], iku, ikv, ikw, ikx2, iky2, ikz2, ikx3, iky3, ikz3, 
	    ikx4, iky4, ikz4, ikxf, ikyf, ikzf, isgo, isgu, isgv, isgw, ixfi2,
	     iyfi2, izfi2, ixfo2, iyfo2, izfo2, ixse2, iyse2, izse2, isgx2, 
	    ixth2, iyth2, izth2, isgy2, isgz2, isgx3, isgy3, isgz3, isgx4, 
	    isgy4, isgz4, isgxf, isgyf, isgzf, ixfiw, iyfiw, izfiw, ixfow, 
	    iyfow, izfow, ixsew, iysew, izsew, ixthw, iythw, izthw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), muldif_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), mulmul_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclp;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    ixfiw = x[*ifir];
    iyfiw = y[*ifir];
    izfiw = z__[*ifir];
    ixsew = x[*isec];
    iysew = y[*isec];
    izsew = z__[*isec];
    ixthw = x[*ithi];
    iythw = y[*ithi];
    izthw = z__[*ithi];
    ixfow = x[*ifou];
    iyfow = y[*ifou];
    izfow = z__[*ifou];

    ixfi2 = x2[*ifir];
    iyfi2 = y2[*ifir];
    izfi2 = z2[*ifir];
    ixse2 = x2[*isec];
    iyse2 = y2[*isec];
    izse2 = z2[*isec];
    ixth2 = x2[*ithi];
    iyth2 = y2[*ithi];
    izth2 = z2[*ithi];
    ixfo2 = x2[*ifou];
    iyfo2 = y2[*ifou];
    izfo2 = z2[*ifou];

    decmp2_(ixf, &isgxf, &ikxf, &ixfiw, &ixfi2, mhalf, mfull, &isclp[1]);
    decmp2_(iyf, &isgyf, &ikyf, &iyfiw, &iyfi2, mhalf, mfull, &isclp[1]);
    decmp2_(izf, &isgzf, &ikzf, &izfiw, &izfi2, mhalf, mfull, &isclp[1]);

    decmp2_(io, &isgo, &iko, &ixsew, &ixse2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix2, &isgo, &isgxf, &isgx2, &iko, &ikxf, &ikx2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iysew, &iyse2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy2, &isgo, &isgyf, &isgy2, &iko, &ikyf, &iky2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izsew, &izse2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz2, &isgo, &isgzf, &isgz2, &iko, &ikzf, &ikz2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &ixthw, &ixth2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix3, &isgo, &isgxf, &isgx3, &iko, &ikxf, &ikx3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iythw, &iyth2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy3, &isgo, &isgyf, &isgy3, &iko, &ikyf, &iky3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izthw, &izth2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz3, &isgo, &isgzf, &isgz3, &iko, &ikzf, &ikz3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &ixfow, &ixfo2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix4, &isgo, &isgxf, &isgx4, &iko, &ikxf, &ikx4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iyfow, &iyfo2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy4, &isgo, &isgyf, &isgy4, &iko, &ikyf, &iky4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izfow, &izfo2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz4, &isgo, &isgzf, &isgz4, &iko, &ikzf, &ikz4, &c__30, 
	    mhalf);

    mulmul_(iy2, iz3, iv, &isgy2, &isgz3, &isgv, &iky2, &ikz3, &ikv, &c__30, 
	    mhalf);
    mulmul_(iz2, iy3, iu, &isgz2, &isgy3, &isgu, &ikz2, &iky3, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, iw, &isgv, &isgu, &isgw, &ikv, &iku, &ikw, &c__30, mhalf);
    mulmul_(iw, ix4, io, &isgw, &isgx4, &isgo, &ikw, &ikx4, &iko, &c__30, 
	    mhalf);

    mulmul_(iz2, ix3, iv, &isgz2, &isgx3, &isgv, &ikz2, &ikx3, &ikv, &c__30, 
	    mhalf);
    mulmul_(ix2, iz3, iu, &isgx2, &isgz3, &isgu, &ikx2, &ikz3, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, iw, &isgv, &isgu, &isgw, &ikv, &iku, &ikw, &c__30, mhalf);
    mulmul_(iw, iy4, iu, &isgw, &isgy4, &isgu, &ikw, &iky4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(io, iu, iw, &isgo, &isgu, &isgw, &iko, &iku, &ikw, &c__30, mhalf);

    mulmul_(ix2, iy3, iv, &isgx2, &isgy3, &isgv, &ikx2, &iky3, &ikv, &c__30, 
	    mhalf);
    mulmul_(iy2, ix3, iu, &isgy2, &isgx3, &isgu, &iky2, &ikx3, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);
    mulmul_(io, iz4, iu, &isgo, &isgz4, &isgu, &iko, &ikz4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iw, iu, io, &isgw, &isgu, &isgo, &ikw, &iku, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* ipsign_ */

/* IPSIG1 */

/* Subroutine */ int ipsig1_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *xc, integer *yc, integer *zc, 
	integer *ifir, integer *isec, integer *ithi, integer *ifou, integer *
	ifif, integer *isix, integer *mhalf, integer *mfull, integer *isclp, 
	integer *ipout)
{
    integer io[30], iu[30], iv[30], ix2[30], iy2[30], iz2[30], ix3[30], iy3[
	    30], iz3[30], ix4[30], iy4[30], iz4[30], ix5[30], iy5[30], iz5[30]
	    , ix6[30], iy6[30], iz6[30], iko, ixf[30], iyf[30], izf[30], iku, 
	    ikv, ikx2, iky2, ikz2, ikx3, iky3, ikz3, ikx4, iky4, ikz4, ikx5, 
	    iky5, ikz5, ikx6, iky6, ikz6, ikxf, ikyf, ikzf, isgo, isgu, isgv, 
	    ixfi2, iyfi2, izfi2, ixfo2, iyfo2, izfo2, isgx2, isgy2, ixth2, 
	    iyth2, izth2, ixsi2, iysi2, izsi2, isgz2, isgx3, isgy3, isgz3, 
	    isgx4, isgy4, isgz4, isgx5, isgy5, isgz5, isgx6, isgy6, isgz6, 
	    isgxf, isgyf, isgzf, ixfiw, iyfiw, izfiw, ixfow, iyfow, izfow, 
	    ixfuw, ixthw, iythw, izthw, ixsiw, iysiw, izsiw, iyfuw, izfuw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    integer ixsuw, iysuw, izsuw;
    extern /* Subroutine */ int decomp_(integer *, integer *, integer *, 
	    integer *), muldif_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), mulmul_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    ixthw = x[*ithi];
    iythw = y[*ithi];
    izthw = z__[*ithi];
    ixfow = x[*ifou];
    iyfow = y[*ifou];
    izfow = z__[*ifou];
    ixfiw = x[*ifif];
    iyfiw = y[*ifif];
    izfiw = z__[*ifif];
    ixsiw = x[*isix];
    iysiw = y[*isix];
    izsiw = z__[*isix];

    ixth2 = x2[*ithi];
    iyth2 = y2[*ithi];
    izth2 = z2[*ithi];
    ixfo2 = x2[*ifou];
    iyfo2 = y2[*ifou];
    izfo2 = z2[*ifou];
    ixfi2 = x2[*ifif];
    iyfi2 = y2[*ifif];
    izfi2 = z2[*ifif];
    ixsi2 = x2[*isix];
    iysi2 = y2[*isix];
    izsi2 = z2[*isix];

    ixfuw = xc[*ifir];
    iyfuw = yc[*ifir];
    izfuw = zc[*ifir];
    ixsuw = xc[*isec];
    iysuw = yc[*isec];
    izsuw = zc[*isec];

    ikxf = 2;
    ikyf = 2;
    ikzf = 2;
    decomp_(ixf, &isgxf, &ixfuw, mhalf);
    decomp_(iyf, &isgyf, &iyfuw, mhalf);
    decomp_(izf, &isgzf, &izfuw, mhalf);

    decmp2_(ix3, &isgx3, &ikx3, &ixthw, &ixth2, mhalf, mfull, &isclp[1]);
    decmp2_(iy3, &isgy3, &iky3, &iythw, &iyth2, mhalf, mfull, &isclp[1]);
    decmp2_(iz3, &isgz3, &ikz3, &izthw, &izth2, mhalf, mfull, &isclp[1]);
    decmp2_(ix5, &isgx5, &ikx5, &ixfiw, &ixfi2, mhalf, mfull, &isclp[1]);
    decmp2_(iy5, &isgy5, &iky5, &iyfiw, &iyfi2, mhalf, mfull, &isclp[1]);
    decmp2_(iz5, &isgz5, &ikz5, &izfiw, &izfi2, mhalf, mfull, &isclp[1]);

    iko = 2;
    decomp_(io, &isgo, &ixsuw, mhalf);
    muldif_(io, ixf, ix2, &isgo, &isgxf, &isgx2, &iko, &ikxf, &ikx2, &c__30, 
	    mhalf);
    decomp_(io, &isgo, &iysuw, mhalf);
    muldif_(io, iyf, iy2, &isgo, &isgyf, &isgy2, &iko, &ikyf, &iky2, &c__30, 
	    mhalf);
    decomp_(io, &isgo, &izsuw, mhalf);
    muldif_(io, izf, iz2, &isgo, &isgzf, &isgz2, &iko, &ikzf, &ikz2, &c__30, 
	    mhalf);

    decmp2_(io, &isgo, &iko, &ixfow, &ixfo2, mhalf, mfull, &isclp[1]);
    muldif_(io, ix3, ix4, &isgo, &isgx3, &isgx4, &iko, &ikx3, &ikx4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iyfow, &iyfo2, mhalf, mfull, &isclp[1]);
    muldif_(io, iy3, iy4, &isgo, &isgy3, &isgy4, &iko, &iky3, &iky4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izfow, &izfo2, mhalf, mfull, &isclp[1]);
    muldif_(io, iz3, iz4, &isgo, &isgz3, &isgz4, &iko, &ikz3, &ikz4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &ixsiw, &ixsi2, mhalf, mfull, &isclp[1]);
    muldif_(io, ix5, ix6, &isgo, &isgx5, &isgx6, &iko, &ikx5, &ikx6, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iysiw, &iysi2, mhalf, mfull, &isclp[1]);
    muldif_(io, iy5, iy6, &isgo, &isgy5, &isgy6, &iko, &iky5, &iky6, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izsiw, &izsi2, mhalf, mfull, &isclp[1]);
    muldif_(io, iz5, iz6, &isgo, &isgz5, &isgz6, &iko, &ikz5, &ikz6, &c__30, 
	    mhalf);

    mulmul_(iy2, iz4, io, &isgy2, &isgz4, &isgo, &iky2, &ikz4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix6, iv, &isgo, &isgx6, &isgv, &iko, &ikx6, &ikv, &c__30, 
	    mhalf);

    mulmul_(iz2, ix4, io, &isgz2, &isgx4, &isgo, &ikz2, &ikx4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iy6, iu, &isgo, &isgy6, &isgu, &iko, &iky6, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix2, iy4, iv, &isgx2, &isgy4, &isgv, &ikx2, &iky4, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iz6, iu, &isgv, &isgz6, &isgu, &ikv, &ikz6, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iz2, iy4, io, &isgz2, &isgy4, &isgo, &ikz2, &iky4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix6, iu, &isgo, &isgx6, &isgu, &iko, &ikx6, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix2, iz4, iv, &isgx2, &isgz4, &isgv, &ikx2, &ikz4, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iy6, iu, &isgv, &isgy6, &isgu, &ikv, &iky6, &iku, &c__30, 
	    mhalf);
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iy2, ix4, io, &isgy2, &isgx4, &isgo, &iky2, &ikx4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iz6, iu, &isgo, &isgz6, &isgu, &iko, &ikz6, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* ipsig1_ */

/* IPSIG2 */

/* Subroutine */ int ipsig2_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *xc, integer *yc, integer *zc, 
	integer *ifir, integer *isec, integer *ithi, integer *ifou, integer *
	ifif, integer *isix, integer *mhalf, integer *mfull, integer *isclp, 
	integer *ipout)
{
    integer io[30], iu[30], iv[30], ix2[30], iy2[30], iz2[30], ix3[30], iy3[
	    30], iz3[30], ix4[30], iy4[30], iz4[30], ix5[30], iy5[30], iz5[30]
	    , ix6[30], iy6[30], iz6[30], iko, ixf[30], iyf[30], izf[30], iku, 
	    ikv, ikx2, iky2, ikz2, ikx3, iky3, ikz3, ikx4, iky4, ikz4, ikx5, 
	    iky5, ikz5, ikx6, iky6, ikz6, ikxf, ikyf, ikzf, isgo, isgu, isgv, 
	    ixfi2, iyfi2, izfi2, isgx2, isgy2, ixsi2, iysi2, izsi2, isgz2, 
	    isgx3, isgy3, isgz3, isgx4, isgy4, isgz4, isgx5, isgy5, isgz5, 
	    isgx6, isgy6, isgz6, isgxf, isgyf, isgzf, ixfiw, iyfiw, izfiw, 
	    ixfow, iyfow, izfow, ixfuw, ixthw, iythw, izthw, ixsiw, iysiw, 
	    izsiw, iyfuw, izfuw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    integer ixsuw, iysuw, izsuw;
    extern /* Subroutine */ int decomp_(integer *, integer *, integer *, 
	    integer *), muldif_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), mulmul_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    ixfiw = x[*ifif];
    iyfiw = y[*ifif];
    izfiw = z__[*ifif];
    ixsiw = x[*isix];
    iysiw = y[*isix];
    izsiw = z__[*isix];

    ixfi2 = x2[*ifif];
    iyfi2 = y2[*ifif];
    izfi2 = z2[*ifif];
    ixsi2 = x2[*isix];
    iysi2 = y2[*isix];
    izsi2 = z2[*isix];

    ixfuw = xc[*ifir];
    iyfuw = yc[*ifir];
    izfuw = zc[*ifir];
    ixsuw = xc[*isec];
    iysuw = yc[*isec];
    izsuw = zc[*isec];
    ixthw = xc[*ithi];
    iythw = yc[*ithi];
    izthw = zc[*ithi];
    ixfow = xc[*ifou];
    iyfow = yc[*ifou];
    izfow = zc[*ifou];

    ikxf = 2;
    ikyf = 2;
    ikzf = 2;
    ikx3 = 2;
    iky3 = 2;
    ikz3 = 2;
    decomp_(ixf, &isgxf, &ixfuw, mhalf);
    decomp_(iyf, &isgyf, &iyfuw, mhalf);
    decomp_(izf, &isgzf, &izfuw, mhalf);
    decomp_(ix3, &isgx3, &ixthw, mhalf);
    decomp_(iy3, &isgy3, &iythw, mhalf);
    decomp_(iz3, &isgz3, &izthw, mhalf);

    decmp2_(ix5, &isgx5, &ikx5, &ixfiw, &ixfi2, mhalf, mfull, &isclp[1]);
    decmp2_(iy5, &isgy5, &iky5, &iyfiw, &iyfi2, mhalf, mfull, &isclp[1]);
    decmp2_(iz5, &isgz5, &ikz5, &izfiw, &izfi2, mhalf, mfull, &isclp[1]);

    iko = 2;
    decomp_(io, &isgo, &ixsuw, mhalf);
    muldif_(io, ixf, ix2, &isgo, &isgxf, &isgx2, &iko, &ikxf, &ikx2, &c__30, 
	    mhalf);
    decomp_(io, &isgo, &iysuw, mhalf);
    muldif_(io, iyf, iy2, &isgo, &isgyf, &isgy2, &iko, &ikyf, &iky2, &c__30, 
	    mhalf);
    decomp_(io, &isgo, &izsuw, mhalf);
    muldif_(io, izf, iz2, &isgo, &isgzf, &isgz2, &iko, &ikzf, &ikz2, &c__30, 
	    mhalf);
    decomp_(io, &isgo, &ixfow, mhalf);
    muldif_(io, ix3, ix4, &isgo, &isgx3, &isgx4, &iko, &ikx3, &ikx4, &c__30, 
	    mhalf);
    decomp_(io, &isgo, &iyfow, mhalf);
    muldif_(io, iy3, iy4, &isgo, &isgy3, &isgy4, &iko, &iky3, &iky4, &c__30, 
	    mhalf);
    decomp_(io, &isgo, &izfow, mhalf);
    muldif_(io, iz3, iz4, &isgo, &isgz3, &isgz4, &iko, &ikz3, &ikz4, &c__30, 
	    mhalf);

    decmp2_(io, &isgo, &iko, &ixsiw, &ixsi2, mhalf, mfull, &isclp[1]);
    muldif_(io, ix5, ix6, &isgo, &isgx5, &isgx6, &iko, &ikx5, &ikx6, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iysiw, &iysi2, mhalf, mfull, &isclp[1]);
    muldif_(io, iy5, iy6, &isgo, &isgy5, &isgy6, &iko, &iky5, &iky6, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izsiw, &izsi2, mhalf, mfull, &isclp[1]);
    muldif_(io, iz5, iz6, &isgo, &isgz5, &isgz6, &iko, &ikz5, &ikz6, &c__30, 
	    mhalf);

    mulmul_(iy2, iz4, io, &isgy2, &isgz4, &isgo, &iky2, &ikz4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix6, iv, &isgo, &isgx6, &isgv, &iko, &ikx6, &ikv, &c__30, 
	    mhalf);

    mulmul_(iz2, ix4, io, &isgz2, &isgx4, &isgo, &ikz2, &ikx4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iy6, iu, &isgo, &isgy6, &isgu, &iko, &iky6, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix2, iy4, iv, &isgx2, &isgy4, &isgv, &ikx2, &iky4, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iz6, iu, &isgv, &isgz6, &isgu, &ikv, &ikz6, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iz2, iy4, io, &isgz2, &isgy4, &isgo, &ikz2, &iky4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix6, iu, &isgo, &isgx6, &isgu, &iko, &ikx6, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix2, iz4, iv, &isgx2, &isgz4, &isgv, &ikx2, &ikz4, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iy6, iu, &isgv, &isgy6, &isgu, &ikv, &iky6, &iku, &c__30, 
	    mhalf);
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iy2, ix4, io, &isgy2, &isgx4, &isgo, &iky2, &ikx4, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iz6, iu, &isgo, &isgz6, &isgu, &iko, &ikz6, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* ipsig2_ */

/* IPSIG3 */

/* Subroutine */ int ipsig3_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *xc, integer *yc, integer *zc, 
	integer *ifir, integer *isec, integer *ithi, integer *ifou, integer *
	mhalf, integer *mfull, integer *isclp, integer *ifn, integer *ipout)
{
    integer io[30], iu[30], iv[30], ix2[30], iy2[30], iz2[30], ix3[30], iy3[
	    30], iz3[30], ix4[30], iy4[30], iz4[30], iko, ixf[30], iyf[30], 
	    izf[30], iku, ikv, ikx2, iky2, ikz2, ikx3, iky3, ikz3, ikx4, iky4,
	     ikz4, ikxf, ikyf, ikzf, isgo, isgu, isgv, ixfi2, iyfi2, izfi2, 
	    ixse2, iyse2, izse2, isgx2, ixth2, iyth2, izth2, isgy2, isgz2, 
	    isgx3, isgy3, isgz3, isgx4, isgy4, isgz4, isgxf, isgyf, isgzf, 
	    ixfiw, iyfiw, izfiw, ixfow, iyfow, izfow, ixsew, iysew, izsew, 
	    ixthw, iythw, izthw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), decomp_(
	    integer *, integer *, integer *, integer *), muldif_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), mulmul_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    ixfiw = x[*ifir];
    iyfiw = y[*ifir];
    izfiw = z__[*ifir];
    ixsew = x[*isec];
    iysew = y[*isec];
    izsew = z__[*isec];
    ixthw = x[*ithi];
    iythw = y[*ithi];
    izthw = z__[*ithi];

    ixfi2 = x2[*ifir];
    iyfi2 = y2[*ifir];
    izfi2 = z2[*ifir];
    ixse2 = x2[*isec];
    iyse2 = y2[*isec];
    izse2 = z2[*isec];
    ixth2 = x2[*ithi];
    iyth2 = y2[*ithi];
    izth2 = z2[*ithi];

    ixfow = xc[*ifou];
    iyfow = yc[*ifou];
    izfow = zc[*ifou];

    decmp2_(ixf, &isgxf, &ikxf, &ixfiw, &ixfi2, mhalf, mfull, &isclp[1]);
    decmp2_(iyf, &isgyf, &ikyf, &iyfiw, &iyfi2, mhalf, mfull, &isclp[1]);
    decmp2_(izf, &isgzf, &ikzf, &izfiw, &izfi2, mhalf, mfull, &isclp[1]);

    ikx4 = 2;
    iky4 = 2;
    ikz4 = 2;
    decomp_(ix4, &isgx4, &ixfow, mhalf);
    decomp_(iy4, &isgy4, &iyfow, mhalf);
    decomp_(iz4, &isgz4, &izfow, mhalf);

    decmp2_(io, &isgo, &iko, &ixsew, &ixse2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix2, &isgo, &isgxf, &isgx2, &iko, &ikxf, &ikx2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iysew, &iyse2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy2, &isgo, &isgyf, &isgy2, &iko, &ikyf, &iky2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izsew, &izse2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz2, &isgo, &isgzf, &isgz2, &iko, &ikzf, &ikz2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &ixthw, &ixth2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix3, &isgo, &isgxf, &isgx3, &iko, &ikxf, &ikx3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iythw, &iyth2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy3, &isgo, &isgyf, &isgy3, &iko, &ikyf, &iky3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izthw, &izth2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz3, &isgo, &isgzf, &isgz3, &iko, &ikzf, &ikz3, &c__30, 
	    mhalf);

    mulmul_(iy4, iz2, io, &isgy4, &isgz2, &isgo, &iky4, &ikz2, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix3, iv, &isgo, &isgx3, &isgv, &iko, &ikx3, &ikv, &c__30, 
	    mhalf);

    mulmul_(iz4, ix2, io, &isgz4, &isgx2, &isgo, &ikz4, &ikx2, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iy3, iu, &isgo, &isgy3, &isgu, &iko, &iky3, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix4, iy2, iv, &isgx4, &isgy2, &isgv, &ikx4, &iky2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iz3, iu, &isgv, &isgz3, &isgu, &ikv, &ikz3, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iz4, iy2, io, &isgz4, &isgy2, &isgo, &ikz4, &iky2, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix3, iu, &isgo, &isgx3, &isgu, &iko, &ikx3, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix4, iz2, iv, &isgx4, &isgz2, &isgv, &ikx4, &ikz2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iy3, iu, &isgv, &isgy3, &isgu, &ikv, &iky3, &iku, &c__30, 
	    mhalf);
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iy4, ix2, io, &isgy4, &isgx2, &isgo, &iky4, &ikx2, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iz3, iu, &isgo, &isgz3, &isgu, &iko, &ikz3, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    *ipout = isgo;
    if (*ifn == 1) {
	*ipout = -isgo;
    }

    return 0;
} /* ipsig3_ */

/* IPSIG4 */

/* Subroutine */ int ipsig4_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *xc, integer *yc, integer *zc, 
	integer *ifir, integer *isec, integer *ithi, integer *ifou, integer *
	mhalf, integer *mfull, integer *isclp, integer *ipout)
{
    integer io[30], iu[30], iv[30], ix2[30], iy2[30], iz2[30], ix3[30], iy3[
	    30], iz3[30], ix4[30], iy4[30], iz4[30], iko, ixf[30], iyf[30], 
	    izf[30], iku, ikv, ikx2, iky2, ikz2, ikx3, iky3, ikz3, ikx4, iky4,
	     ikz4, ikxf, ikyf, ikzf, isgo, isgu, isgv, ixfo2, iyfo2, izfo2, 
	    isgx2, isgy2, ixth2, iyth2, izth2, isgz2, isgx3, isgy3, isgz3, 
	    isgx4, isgy4, isgz4, isgxf, isgyf, isgzf, ixfiw, iyfiw, izfiw, 
	    ixfow, iyfow, izfow, ixsew, iysew, izsew, ixthw, iythw, izthw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), decomp_(
	    integer *, integer *, integer *, integer *), muldif_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), mulmul_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    ixthw = x[*ithi];
    iythw = y[*ithi];
    izthw = z__[*ithi];
    ixfow = x[*ifou];
    iyfow = y[*ifou];
    izfow = z__[*ifou];

    ixth2 = x2[*ithi];
    iyth2 = y2[*ithi];
    izth2 = z2[*ithi];
    ixfo2 = x2[*ifou];
    iyfo2 = y2[*ifou];
    izfo2 = z2[*ifou];

    ixfiw = xc[*ifir];
    iyfiw = yc[*ifir];
    izfiw = zc[*ifir];
    ixsew = xc[*isec];
    iysew = yc[*isec];
    izsew = zc[*isec];

    decmp2_(ixf, &isgxf, &ikxf, &ixfow, &ixfo2, mhalf, mfull, &isclp[1]);
    decmp2_(iyf, &isgyf, &ikyf, &iyfow, &iyfo2, mhalf, mfull, &isclp[1]);
    decmp2_(izf, &isgzf, &ikzf, &izfow, &izfo2, mhalf, mfull, &isclp[1]);

    ikx2 = 2;
    iky2 = 2;
    ikz2 = 2;
    ikx3 = 2;
    iky3 = 2;
    ikz3 = 2;
    decomp_(ix2, &isgx2, &ixfiw, mhalf);
    decomp_(iy2, &isgy2, &iyfiw, mhalf);
    decomp_(iz2, &isgz2, &izfiw, mhalf);
    decomp_(ix3, &isgx3, &ixsew, mhalf);
    decomp_(iy3, &isgy3, &iysew, mhalf);
    decomp_(iz3, &isgz3, &izsew, mhalf);

    decmp2_(io, &isgo, &iko, &ixthw, &ixth2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix4, &isgo, &isgxf, &isgx4, &iko, &ikxf, &ikx4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iythw, &iyth2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy4, &isgo, &isgyf, &isgy4, &iko, &ikyf, &iky4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izthw, &izth2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz4, &isgo, &isgzf, &isgz4, &iko, &ikzf, &ikz4, &c__30, 
	    mhalf);

    mulmul_(iy2, iz3, io, &isgy2, &isgz3, &isgo, &iky2, &ikz3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix4, iv, &isgo, &isgx4, &isgv, &iko, &ikx4, &ikv, &c__30, 
	    mhalf);

    mulmul_(iz2, ix3, io, &isgz2, &isgx3, &isgo, &ikz2, &ikx3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iy4, iu, &isgo, &isgy4, &isgu, &iko, &iky4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix2, iy3, iv, &isgx2, &isgy3, &isgv, &ikx2, &iky3, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iz4, iu, &isgv, &isgz4, &isgu, &ikv, &ikz4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iz2, iy3, io, &isgz2, &isgy3, &isgo, &ikz2, &iky3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix4, iu, &isgo, &isgx4, &isgu, &iko, &ikx4, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix2, iz3, iv, &isgx2, &isgz3, &isgv, &ikx2, &ikz3, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iy4, iu, &isgv, &isgy4, &isgu, &ikv, &iky4, &iku, &c__30, 
	    mhalf);
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iy2, ix3, io, &isgy2, &isgx3, &isgo, &iky2, &ikx3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iz4, iu, &isgo, &isgz4, &isgu, &iko, &ikz4, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* ipsig4_ */

/* IPSIG6 */

/* Subroutine */ int ipsig6_(integer *x, integer *y, integer *z__, integer *
	x2, integer *y2, integer *z2, integer *xc, integer *yc, integer *zc, 
	integer *ifir, integer *isec, integer *ithi, integer *ifou, integer *
	mhalf, integer *mfull, integer *isclp, integer *ipout)
{
    integer io[30], iu[30], iv[30], ix2[30], iy2[30], iz2[30], ix3[30], iy3[
	    30], iz3[30], ix4[30], iy4[30], iz4[30], ix5[30], iy5[30], iz5[30]
	    , ix6[30], iy6[30], iz6[30], iko, ixf[30], iyf[30], izf[30], iku, 
	    ikv, ikx2, iky2, ikz2, ikx3, iky3, ikz3, ikx4, iky4, ikz4, ikx5, 
	    iky5, ikz5, ikx6, iky6, ikz6, ikxf, ikyf, ikzf, isgo, isgu, isgv, 
	    ixfi2, iyfi2, izfi2, ixfo2, iyfo2, izfo2, ixse2, iyse2, izse2, 
	    isgx2, ixth2, iyth2, izth2, isgy2, isgz2, isgx3, isgy3, isgz3, 
	    isgx4, isgy4, isgz4, isgx5, isgy5, isgz5, isgx6, isgy6, isgz6, 
	    isgxf, isgyf, isgzf, ixfiw, iyfiw, izfiw, ixfow, iyfow, izfow, 
	    ixsew, iysew, izsew, ixthw, iythw, izthw, ixfuw, iyfuw, izfuw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);
    integer ixsuw, iysuw, izsuw;
    extern /* Subroutine */ int decomp_(integer *, integer *, integer *, 
	    integer *), muldif_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), mulmul_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    --isclp;
    --zc;
    --yc;
    --xc;
    --z2;
    --y2;
    --x2;
    --z__;
    --y;
    --x;

    /* Function Body */
    ixfiw = x[*ifir];
    iyfiw = y[*ifir];
    izfiw = z__[*ifir];
    ixsew = x[*isec];
    iysew = y[*isec];
    izsew = z__[*isec];
    ixthw = x[*ithi];
    iythw = y[*ithi];
    izthw = z__[*ithi];
    ixfow = x[*ifou];
    iyfow = y[*ifou];
    izfow = z__[*ifou];

    ixfi2 = x2[*ifir];
    iyfi2 = y2[*ifir];
    izfi2 = z2[*ifir];
    ixse2 = x2[*isec];
    iyse2 = y2[*isec];
    izse2 = z2[*isec];
    ixth2 = x2[*ithi];
    iyth2 = y2[*ithi];
    izth2 = z2[*ithi];
    ixfo2 = x2[*ifou];
    iyfo2 = y2[*ifou];
    izfo2 = z2[*ifou];

    ixfuw = xc[*ifir];
    iyfuw = yc[*ifir];
    izfuw = zc[*ifir];
    ixsuw = xc[*isec];
    iysuw = yc[*isec];
    izsuw = zc[*isec];

    decmp2_(ixf, &isgxf, &ikxf, &ixfow, &ixfo2, mhalf, mfull, &isclp[1]);
    decmp2_(iyf, &isgyf, &ikyf, &iyfow, &iyfo2, mhalf, mfull, &isclp[1]);
    decmp2_(izf, &isgzf, &ikzf, &izfow, &izfo2, mhalf, mfull, &isclp[1]);

    ikx5 = 2;
    iky5 = 2;
    ikz5 = 2;
    ikx6 = 2;
    iky6 = 2;
    ikz6 = 2;
    decomp_(ix5, &isgx5, &ixfuw, mhalf);
    decomp_(iy5, &isgy5, &iyfuw, mhalf);
    decomp_(iz5, &isgz5, &izfuw, mhalf);
    decomp_(ix6, &isgx6, &ixsuw, mhalf);
    decomp_(iy6, &isgy6, &iysuw, mhalf);
    decomp_(iz6, &isgz6, &izsuw, mhalf);

    decmp2_(io, &isgo, &iko, &ixfiw, &ixfi2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix2, &isgo, &isgxf, &isgx2, &iko, &ikxf, &ikx2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iyfiw, &iyfi2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy2, &isgo, &isgyf, &isgy2, &iko, &ikyf, &iky2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izfiw, &izfi2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz2, &isgo, &isgzf, &isgz2, &iko, &ikzf, &ikz2, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &ixsew, &ixse2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix3, &isgo, &isgxf, &isgx3, &iko, &ikxf, &ikx3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iysew, &iyse2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy3, &isgo, &isgyf, &isgy3, &iko, &ikyf, &iky3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izsew, &izse2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz3, &isgo, &isgzf, &isgz3, &iko, &ikzf, &ikz3, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &ixthw, &ixth2, mhalf, mfull, &isclp[1]);
    muldif_(io, ixf, ix4, &isgo, &isgxf, &isgx4, &iko, &ikxf, &ikx4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &iythw, &iyth2, mhalf, mfull, &isclp[1]);
    muldif_(io, iyf, iy4, &isgo, &isgyf, &isgy4, &iko, &ikyf, &iky4, &c__30, 
	    mhalf);
    decmp2_(io, &isgo, &iko, &izthw, &izth2, mhalf, mfull, &isclp[1]);
    muldif_(io, izf, iz4, &isgo, &isgzf, &isgz4, &iko, &ikzf, &ikz4, &c__30, 
	    mhalf);

    mulmul_(iy5, iz3, io, &isgy5, &isgz3, &isgo, &iky5, &ikz3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix4, iv, &isgo, &isgx4, &isgv, &iko, &ikx4, &ikv, &c__30, 
	    mhalf);

    mulmul_(iz5, ix3, io, &isgz5, &isgx3, &isgo, &ikz5, &ikx3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iy4, iu, &isgo, &isgy4, &isgu, &iko, &iky4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix5, iy3, iv, &isgx5, &isgy3, &isgv, &ikx5, &iky3, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iz4, iu, &isgv, &isgz4, &isgu, &ikv, &ikz4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iz5, iy3, io, &isgz5, &isgy3, &isgo, &ikz5, &iky3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix4, iu, &isgo, &isgx4, &isgu, &iko, &ikx4, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix5, iz3, iv, &isgx5, &isgz3, &isgv, &ikx5, &ikz3, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iy4, iu, &isgv, &isgy4, &isgu, &ikv, &iky4, &iku, &c__30, 
	    mhalf);
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iy5, ix3, io, &isgy5, &isgx3, &isgo, &iky5, &ikx3, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iz4, iu, &isgo, &isgz4, &isgu, &iko, &ikz4, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(iy6, iz2, iv, &isgy6, &isgz2, &isgv, &iky6, &ikz2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, ix4, iu, &isgv, &isgx4, &isgu, &ikv, &ikx4, &iku, &c__30, 
	    mhalf);
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iz6, ix2, io, &isgz6, &isgx2, &isgo, &ikz6, &ikx2, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iy4, iu, &isgo, &isgy4, &isgu, &iko, &iky4, &iku, &c__30, 
	    mhalf);
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix6, iy2, iv, &isgx6, &isgy2, &isgv, &ikx6, &iky2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iz4, iu, &isgv, &isgz4, &isgu, &ikv, &ikz4, &iku, &c__30, 
	    mhalf);
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iz6, iy2, io, &isgz6, &isgy2, &isgo, &ikz6, &iky2, &iko, &c__30, 
	    mhalf);
    mulmul_(io, ix4, iu, &isgo, &isgx4, &isgu, &iko, &ikx4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    mulmul_(ix6, iz2, iv, &isgx6, &isgz2, &isgv, &ikx6, &ikz2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, iy4, iu, &isgv, &isgy4, &isgu, &ikv, &iky4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(io, iu, iv, &isgo, &isgu, &isgv, &iko, &iku, &ikv, &c__30, mhalf);

    mulmul_(iy6, ix2, io, &isgy6, &isgx2, &isgo, &iky6, &ikx2, &iko, &c__30, 
	    mhalf);
    mulmul_(io, iz4, iu, &isgo, &isgz4, &isgu, &iko, &ikz4, &iku, &c__30, 
	    mhalf);
    isgu = -isgu;
    muldif_(iv, iu, io, &isgv, &isgu, &isgo, &ikv, &iku, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* ipsig6_ */

/* DECMP2 */

/*     subroutine decmp2 */

/*     to decompose a regular or non-regular length integer */

/* Subroutine */ int decmp2_(integer *ia, integer *isga, integer *ika, 
	integer *iwi, integer *iwi2, integer *mhalf, integer *mfull, integer *
	isclp)
{
    integer io[30], iu[30], iko, iku, ikcl, isgo, isgu, isgcl;
    extern /* Subroutine */ int decomp_(integer *, integer *, integer *, 
	    integer *), muldif_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), mulmul_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);



    /* Parameter adjustments */
    --isclp;
    --ia;

    /* Function Body */
    decomp_(&ia[1], isga, iwi, mhalf);
    *ika = 2;
    if (*iwi2 != 0) {
	isgcl = 1;
	ikcl = 2;
	mulmul_(&ia[1], &isclp[1], iu, isga, &isgcl, &isgu, ika, &ikcl, &iku, 
		&c__30, mhalf);
	if (*iwi2 == *mfull) {
	    *iwi2 = 0;
	}
	decomp_(io, &isgo, iwi2, mhalf);
	isgo = -isgo;
	iko = 2;
	muldif_(iu, io, &ia[1], &isgu, &isgo, isga, &iku, &iko, ika, &c__30, 
		mhalf);
    }

    return 0;
} /* decmp2_ */

/* DECOMP */

/*     subroutine decomp */

/*     to decompose a regular length integer */

/*     iwi = isga*(ia(1) + ia(2) * mhalf) */

/*     iwi  is a regular length integer */
/*     isga is a sign integer (-1, 0, 1) */
/*     ia(1) and ia(2) are integers less than mhalf */

/* Subroutine */ int decomp_(integer *ia, integer *isga, integer *iwi, 
	integer *mhalf)
{
    integer ivi;



    /* Parameter adjustments */
    --ia;

    /* Function Body */
    if (*iwi > 0) {
	*isga = 1;
	ivi = *iwi;
    } else if (*iwi < 0) {
	*isga = -1;
	ivi = -(*iwi);
    } else {
	*isga = 0;
	ia[1] = 0;
	ia[2] = 0;
	return 0;
    }
    ia[2] = ivi / *mhalf;
    ia[1] = ivi - ia[2] * *mhalf;

    return 0;
} /* decomp_ */

/* MULMUL */

/*     subroutine mulmul */

/*     to perform a multiple precision integer multiplication */
/*     (for multiplying 2 or more integers) */

/*     io = ia * ib */

/*     ia represents a decomposed integer */
/*     ib represents a decomposed integer */
/*     io is the product of ia and ib in its decomposed form */

/* Subroutine */ int mulmul_(integer *ia, integer *ib, integer *io, integer *
	isga, integer *isgb, integer *isgo, integer *ika, integer *ikb, 
	integer *iko, integer *nkmax, integer *mhalf)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, j, k, ipr, ipt, iko1;



    /* Parameter adjustments */
    --io;
    --ib;
    --ia;

    /* Function Body */
    if (*isga == 0 || *isgb == 0) {
	*isgo = 0;
	*iko = 2;
	io[1] = 0;
	io[2] = 0;
	return 0;
    }

    *iko = *ika + *ikb;
    if (*iko > *nkmax) {
	s_stop("4710", (ftnlen)4);
    }

    if (*isga > 0) {
	if (*isgb > 0) {
	    *isgo = 1;
	} else {
	    *isgo = -1;
	}
    } else {
	if (*isgb > 0) {
	    *isgo = -1;
	} else {
	    *isgo = 1;
	}
    }

    iko1 = *iko - 1;
    ipr = 0;

    i__1 = iko1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipt = ipr;
	k = i__;
	i__2 = *ikb;
	for (j = 1; j <= i__2; ++j) {
	    if (k < 1) {
		goto L190;
	    }
	    if (k > *ika) {
		goto L150;
	    }
	    ipt += ia[k] * ib[j];
L150:
	    --k;
/* L180: */
	}
L190:
	ipr = ipt / *mhalf;
	io[i__] = ipt - ipr * *mhalf;
/* L200: */
    }

    io[*iko] = ipr;
    if (ipr >= *mhalf) {
	s_stop("4720", (ftnlen)4);
    }

    iko1 = *iko;
    i__1 = *ika + 1;
    for (i__ = iko1; i__ >= i__1; --i__) {
	if (io[i__] != 0) {
	    goto L400;
	}
	--(*iko);
/* L300: */
    }
L400:

    return 0;
} /* mulmul_ */

/* MULDIF */

/*     subroutine muldif */

/*     to perform a multiple precision integer subtraction */
/*     (for subtracting a decomposed product from another) */

/*     io = ia - ib */

/*     ia represents a decomposed regular length integer or the */
/*        decomposed product of two or more regular length integers */
/*     ib is similarly described */
/*     io is a decomposed integer which represents ia - ib */

/* Subroutine */ int muldif_(integer *ia, integer *ib, integer *io, integer *
	isga, integer *isgb, integer *isgo, integer *ika, integer *ikb, 
	integer *iko, integer *nkmax, integer *mhalf)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer i__, iko1, irel;



    /* Parameter adjustments */
    --io;
    --ib;
    --ia;

    /* Function Body */
    if (*isgb == 0) {
	if (*isga == 0) {
	    *isgo = 0;
	    *iko = 2;
	    io[1] = 0;
	    io[2] = 0;
	    return 0;
	}
	*isgo = *isga;
	*iko = *ika;
	i__1 = *iko;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io[i__] = ia[i__];
/* L100: */
	}
    } else if (*isga == 0) {
	*isgo = -(*isgb);
	*iko = *ikb;
	i__1 = *iko;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    io[i__] = ib[i__];
/* L200: */
	}
    } else {
	*iko = *ika;
	if (*ikb < *ika) {
	    i__1 = *ika;
	    for (i__ = *ikb + 1; i__ <= i__1; ++i__) {
		ib[i__] = 0;
/* L300: */
	    }
	} else if (*ika < *ikb) {
	    *iko = *ikb;
	    i__1 = *ikb;
	    for (i__ = *ika + 1; i__ <= i__1; ++i__) {
		ia[i__] = 0;
/* L400: */
	    }
	}
	if (*isga * *isgb > 0) {
	    irel = 0;
	    for (i__ = *iko; i__ >= 1; --i__) {
		if (ia[i__] > ib[i__]) {
		    irel = 1;
		    goto L600;
		} else if (ia[i__] < ib[i__]) {
		    irel = -1;
		    goto L600;
		}
/* L500: */
	    }
L600:
	    if (irel == 0) {
		*isgo = 0;
		i__1 = *iko;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    io[i__] = 0;
/* L700: */
		}
	    } else {
		*isgo = *isga * irel;
		io[1] = (ia[1] - ib[1]) * irel;
		i__1 = *iko;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    if (io[i__ - 1] < 0) {
			io[i__] = -1;
			io[i__ - 1] += *mhalf;
		    } else {
			io[i__] = 0;
		    }
		    io[i__] += (ia[i__] - ib[i__]) * irel;
/* L800: */
		}
		if (io[*iko] < 0) {
		    s_stop("4810", (ftnlen)4);
		}
	    }
	} else {
	    *isgo = *isga;
	    io[1] = ia[1] + ib[1];
	    i__1 = *iko;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		if (io[i__ - 1] >= *mhalf) {
		    io[i__] = 1;
		    io[i__ - 1] -= *mhalf;
		} else {
		    io[i__] = 0;
		}
		io[i__] = io[i__] + ia[i__] + ib[i__];
/* L900: */
	    }
	    if (io[*iko] >= *mhalf) {
		++(*iko);
		if (*iko > *nkmax) {
		    s_stop("4820", (ftnlen)4);
		}
		io[*iko] = 1;
		io[*iko - 1] -= *mhalf;
	    }
	}
    }

    if (*iko == 2) {
	goto L1400;
    }
    iko1 = *iko;
    for (i__ = iko1; i__ >= 3; --i__) {
	if (io[i__] != 0) {
	    goto L1400;
	}
	--(*iko);
/* L1300: */
    }
L1400:

    return 0;
} /* muldif_ */

/* IWSIGN */

/*     subroutine for determining sign of weight of point ifir minus */
/*     weight of point isec */

/* Subroutine */ int iwsign_(integer *w, integer *w2, integer *ifir, integer *
	isec, integer *mhalf, integer *mfull, integer *isclw, integer *ipout)
{
    integer iu[30], iw1[30], iw2[30], iku, ikw1, ikw2, isgu, iwfi2, iwse2, 
	    isgw1, isgw2, iwfiw, iwsew;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), muldif_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclw;
    --w2;
    --w;

    /* Function Body */
    iwfiw = w[*ifir];
    iwsew = w[*isec];

    iwfi2 = w2[*ifir];
    iwse2 = w2[*isec];

    decmp2_(iw1, &isgw1, &ikw1, &iwfiw, &iwfi2, mhalf, mfull, &isclw[1]);
    decmp2_(iw2, &isgw2, &ikw2, &iwsew, &iwse2, mhalf, mfull, &isclw[1]);
    muldif_(iw1, iw2, iu, &isgw1, &isgw2, &isgu, &ikw1, &ikw2, &iku, &c__30, 
	    mhalf);

    *ipout = isgu;

    return 0;
} /* iwsign_ */

/* IQSIGN */

/*     subroutine for determining position of point ifif with respect */
/*     to sphere determined by (weighted) points ifir, isec, ithi, ifou */
/*     if positive then ifif is outside the sphere */
/*     if negative then ifif is inside the sphere */
/*     if zero then ifif is in the surface of the sphere */

/* Subroutine */ int iqsign_(integer *x, integer *y, integer *z__, integer *w,
	 integer *x2, integer *y2, integer *z2, integer *w2, integer *ifir, 
	integer *isec, integer *ithi, integer *ifou, integer *ifif, integer *
	mhalf, integer *mfull, integer *isclp, integer *isclw, integer *isclr,
	 logical *delaun, integer *ipout)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer io[30], ip[30], iu[30], iv[30], iq2[30], iq3[30], iq4[30], iq5[30]
	    , iw2[30], ix2[30], iy2[30], iz2[30], ix3[30], iy3[30], iz3[30], 
	    ix4[30], iy4[30], iz4[30], ix5[30], iy5[30], iz5[30], iw3[30], 
	    iw4[30], iw5[30], iko, ikp, iwf[30], ixf[30], iyf[30], izf[30], 
	    iku, ikv, ikq2, ikq3, ixf2[30], iyf2[30], izf2[30], ikq4, ikw2, 
	    ikw3, ikw4, ikw5, ikq5, ikx2, iky2, ikz2, ikx3, iky3, ikz3, ikx4, 
	    iky4, ikz4, ikx5, iky5, ikz5, ikcl, ikwf, ikxf, ikyf, ikzf, isgo, 
	    isgp, isgu, isgv, iwfi2, ixfi2, iyfi2, izfi2, ikxf2, isgq2, iwfo2,
	     ixfo2, iyfo2, iwse2, izfo2, ixse2, iwfu2, iwth2, ixth2, iyth2, 
	    izth2, iyse2, izse2, ixfu2, iyfu2, izfu2, isgw2, isgw3, isgw4, 
	    isgw5, isgq3, isgq4, isgq5, ikyf2, ikzf2, isgx2, isgy2, isgz2, 
	    isgx3, isgy3, isgz3, isgx4, isgy4, isgz4, isgx5, isgy5, isgz5, 
	    isgcl, isgxf, isgyf, isgzf, isgwf, iwfiw, ixfiw, iyfiw, izfiw, 
	    iwfow, ixfow, iyfow, iwsew, izfow, iwfuw, iwthw, ixthw, iythw, 
	    izthw, ixsew, iysew, izsew, ixfuw, iyfuw, izfuw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), detrm3_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    integer isgxf2, isgyf2, isgzf2;
    extern /* Subroutine */ int muldif_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), frterm_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), mulmul_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    if (*delaun) {
	isgw2 = 0;
	isgw3 = 0;
	isgw4 = 0;
	isgw5 = 0;
    } else {
	iwfuw = w[*ifir];
	iwsew = w[*isec];
	iwthw = w[*ithi];
	iwfow = w[*ifou];
	iwfiw = w[*ifif];

	iwfu2 = w2[*ifir];
	iwse2 = w2[*isec];
	iwth2 = w2[*ithi];
	iwfo2 = w2[*ifou];
	iwfi2 = w2[*ifif];

	decmp2_(iwf, &isgwf, &ikwf, &iwfuw, &iwfu2, mhalf, mfull, &isclw[1]);
	isgcl = 1;
	ikcl = 2;
	decmp2_(io, &isgo, &iko, &iwsew, &iwse2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw2, &isgu, &isgcl, &isgw2, &iku, &ikcl, &ikw2,
		 &c__30, mhalf);
	decmp2_(io, &isgo, &iko, &iwthw, &iwth2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw3, &isgu, &isgcl, &isgw3, &iku, &ikcl, &ikw3,
		 &c__30, mhalf);
	decmp2_(io, &isgo, &iko, &iwfow, &iwfo2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw4, &isgu, &isgcl, &isgw4, &iku, &ikcl, &ikw4,
		 &c__30, mhalf);
	decmp2_(io, &isgo, &iko, &iwfiw, &iwfi2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw5, &isgu, &isgcl, &isgw5, &iku, &ikcl, &ikw5,
		 &c__30, mhalf);
    }

    ixfuw = x[*ifir];
    iyfuw = y[*ifir];
    izfuw = z__[*ifir];
    ixsew = x[*isec];
    iysew = y[*isec];
    izsew = z__[*isec];
    ixthw = x[*ithi];
    iythw = y[*ithi];
    izthw = z__[*ithi];
    ixfow = x[*ifou];
    iyfow = y[*ifou];
    izfow = z__[*ifou];
    ixfiw = x[*ifif];
    iyfiw = y[*ifif];
    izfiw = z__[*ifif];

    ixfu2 = x2[*ifir];
    iyfu2 = y2[*ifir];
    izfu2 = z2[*ifir];
    ixse2 = x2[*isec];
    iyse2 = y2[*isec];
    izse2 = z2[*isec];
    ixth2 = x2[*ithi];
    iyth2 = y2[*ithi];
    izth2 = z2[*ithi];
    ixfo2 = x2[*ifou];
    iyfo2 = y2[*ifou];
    izfo2 = z2[*ifou];
    ixfi2 = x2[*ifif];
    iyfi2 = y2[*ifif];
    izfi2 = z2[*ifif];

    decmp2_(ixf, &isgxf, &ikxf, &ixfuw, &ixfu2, mhalf, mfull, &isclp[1]);
    decmp2_(iyf, &isgyf, &ikyf, &iyfuw, &iyfu2, mhalf, mfull, &isclp[1]);
    decmp2_(izf, &isgzf, &ikzf, &izfuw, &izfu2, mhalf, mfull, &isclp[1]);
    mulmul_(ixf, ixf, ixf2, &isgxf, &isgxf, &isgxf2, &ikxf, &ikxf, &ikxf2, &
	    c__30, mhalf);
    mulmul_(iyf, iyf, iyf2, &isgyf, &isgyf, &isgyf2, &ikyf, &ikyf, &ikyf2, &
	    c__30, mhalf);
    mulmul_(izf, izf, izf2, &isgzf, &isgzf, &isgzf2, &ikzf, &ikzf, &ikzf2, &
	    c__30, mhalf);
    if (isgxf2 < 0 || isgyf2 < 0 || isgzf2 < 0) {
	s_stop("5105", (ftnlen)4);
    }

    frterm_(&ixsew, &iysew, &izsew, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw2, ix2, iy2, iz2, iq2, &isgw2, &isgx2, &
	    isgy2, &isgz2, &isgq2, &ikw2, &ikx2, &iky2, &ikz2, &ikq2, mhalf, 
	    mfull, &ixse2, &iyse2, &izse2, &isclp[1]);

    frterm_(&ixthw, &iythw, &izthw, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw3, ix3, iy3, iz3, iq3, &isgw3, &isgx3, &
	    isgy3, &isgz3, &isgq3, &ikw3, &ikx3, &iky3, &ikz3, &ikq3, mhalf, 
	    mfull, &ixth2, &iyth2, &izth2, &isclp[1]);

    frterm_(&ixfow, &iyfow, &izfow, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw4, ix4, iy4, iz4, iq4, &isgw4, &isgx4, &
	    isgy4, &isgz4, &isgq4, &ikw4, &ikx4, &iky4, &ikz4, &ikq4, mhalf, 
	    mfull, &ixfo2, &iyfo2, &izfo2, &isclp[1]);

    frterm_(&ixfiw, &iyfiw, &izfiw, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw5, ix5, iy5, iz5, iq5, &isgw5, &isgx5, &
	    isgy5, &isgz5, &isgq5, &ikw5, &ikx5, &iky5, &ikz5, &ikq5, mhalf, 
	    mfull, &ixfi2, &iyfi2, &izfi2, &isclp[1]);

    mulmul_(iq5, ix2, iv, &isgq5, &isgx2, &isgv, &ikq5, &ikx2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iq5, ix3, iu, &isgq5, &isgx3, &isgu, &ikq5, &ikx3, &iku, &c__30, 
	    mhalf);
    mulmul_(iq5, ix4, ip, &isgq5, &isgx4, &isgp, &ikq5, &ikx4, &ikp, &c__30, 
	    mhalf);
    detrm3_(iv, iy2, iz2, &isgv, &isgy2, &isgz2, iu, iy3, iz3, &isgu, &isgy3, 
	    &isgz3, ip, iy4, iz4, &isgp, &isgy4, &isgz4, &ikv, &iku, &ikp, &
	    iky2, &iky3, &iky4, &ikz2, &ikz3, &ikz4, io, &isgo, &iko, mhalf);

    detrm3_(iq2, iy2, iz2, &isgq2, &isgy2, &isgz2, iq3, iy3, iz3, &isgq3, &
	    isgy3, &isgz3, iq4, iy4, iz4, &isgq4, &isgy4, &isgz4, &ikq2, &
	    ikq3, &ikq4, &iky2, &iky3, &iky4, &ikz2, &ikz3, &ikz4, iu, &isgu, 
	    &iku, mhalf);
    mulmul_(iu, ix5, ip, &isgu, &isgx5, &isgp, &iku, &ikx5, &ikp, &c__30, 
	    mhalf);
    muldif_(io, ip, iv, &isgo, &isgp, &isgv, &iko, &ikp, &ikv, &c__30, mhalf);

    detrm3_(iq2, iz2, ix2, &isgq2, &isgz2, &isgx2, iq3, iz3, ix3, &isgq3, &
	    isgz3, &isgx3, iq4, iz4, ix4, &isgq4, &isgz4, &isgx4, &ikq2, &
	    ikq3, &ikq4, &ikz2, &ikz3, &ikz4, &ikx2, &ikx3, &ikx4, iu, &isgu, 
	    &iku, mhalf);
    mulmul_(iu, iy5, io, &isgu, &isgy5, &isgo, &iku, &iky5, &iko, &c__30, 
	    mhalf);
    muldif_(iv, io, ip, &isgv, &isgo, &isgp, &ikv, &iko, &ikp, &c__30, mhalf);

    detrm3_(iq2, ix2, iy2, &isgq2, &isgx2, &isgy2, iq3, ix3, iy3, &isgq3, &
	    isgx3, &isgy3, iq4, ix4, iy4, &isgq4, &isgx4, &isgy4, &ikq2, &
	    ikq3, &ikq4, &ikx2, &ikx3, &ikx4, &iky2, &iky3, &iky4, iu, &isgu, 
	    &iku, mhalf);
    mulmul_(iu, iz5, iv, &isgu, &isgz5, &isgv, &iku, &ikz5, &ikv, &c__30, 
	    mhalf);
    muldif_(ip, iv, io, &isgp, &isgv, &isgo, &ikp, &ikv, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* iqsign_ */

/* IQSIG1 */

/*     subroutine for determining position of point ifif with respect to */
/*     line segment determined by (weighted) points ifir, isec, assuming */
/*     ifir, isec and ifif are colinear */
/*     if positive then ifif is in the exterior of the line segment */
/*     if negative then ifif is in the interior of the line segment */
/*     if zero then ifif is one of the endpoints */

/* Subroutine */ int iqsig1_(integer *x, integer *y, integer *z__, integer *w,
	 integer *x2, integer *y2, integer *z2, integer *w2, integer *ifir, 
	integer *isec, integer *ifif, integer *mhalf, integer *mfull, integer 
	*isclp, integer *isclw, integer *isclr, logical *delaun, integer *
	ipout)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer io[30], iu[30], iv[30], iq2[30], iq5[30], iw2[30], ix2[30], iy2[
	    30], iz2[30], ix5[30], iy5[30], iz5[30], iw5[30], iko, iwf[30], 
	    ixf[30], iyf[30], izf[30], iku, ikv, ikq2, ixf2[30], iyf2[30], 
	    izf2[30], ikq5, ikw2, ikx2, iky2, ikw5, ikz2, ikx5, iky5, ikz5, 
	    ikcl, ikwf, ikxf, ikyf, ikzf, isgo, isgu, isgv, iwfi2, ixfi2, 
	    iyfi2, izfi2, ikxf2, isgq2, ikyf2, ikzf2, isgq5, iwse2, ixse2, 
	    iyse2, iwfu2, izse2, ixfu2, iyfu2, izfu2, isgw2, isgw5, isgx2, 
	    isgy2, isgz2, isgx5, isgy5, isgz5, isgcl, isgxf, isgyf, isgzf, 
	    isgwf, iwfiw, ixfiw, iyfiw, izfiw, iwsew, ixsew, iwfuw, iysew, 
	    izsew, ixfuw, iyfuw, izfuw;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), detrm1_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *);
    integer isgxf2, isgyf2, isgzf2;
    extern /* Subroutine */ int muldif_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), frterm_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), mulmul_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    if (*delaun) {
	isgw2 = 0;
	isgw5 = 0;
    } else {
	iwfuw = w[*ifir];
	iwsew = w[*isec];
	iwfiw = w[*ifif];

	iwfu2 = w2[*ifir];
	iwse2 = w2[*isec];
	iwfi2 = w2[*ifif];

	decmp2_(iwf, &isgwf, &ikwf, &iwfuw, &iwfu2, mhalf, mfull, &isclw[1]);
	isgcl = 1;
	ikcl = 2;
	decmp2_(io, &isgo, &iko, &iwsew, &iwse2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw2, &isgu, &isgcl, &isgw2, &iku, &ikcl, &ikw2,
		 &c__30, mhalf);
	decmp2_(io, &isgo, &iko, &iwfiw, &iwfi2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw5, &isgu, &isgcl, &isgw5, &iku, &ikcl, &ikw5,
		 &c__30, mhalf);
    }

    ixfuw = x[*ifir];
    iyfuw = y[*ifir];
    izfuw = z__[*ifir];
    ixsew = x[*isec];
    iysew = y[*isec];
    izsew = z__[*isec];
    ixfiw = x[*ifif];
    iyfiw = y[*ifif];
    izfiw = z__[*ifif];

    ixfu2 = x2[*ifir];
    iyfu2 = y2[*ifir];
    izfu2 = z2[*ifir];
    ixse2 = x2[*isec];
    iyse2 = y2[*isec];
    izse2 = z2[*isec];
    ixfi2 = x2[*ifif];
    iyfi2 = y2[*ifif];
    izfi2 = z2[*ifif];

    decmp2_(ixf, &isgxf, &ikxf, &ixfuw, &ixfu2, mhalf, mfull, &isclp[1]);
    decmp2_(iyf, &isgyf, &ikyf, &iyfuw, &iyfu2, mhalf, mfull, &isclp[1]);
    decmp2_(izf, &isgzf, &ikzf, &izfuw, &izfu2, mhalf, mfull, &isclp[1]);
    mulmul_(ixf, ixf, ixf2, &isgxf, &isgxf, &isgxf2, &ikxf, &ikxf, &ikxf2, &
	    c__30, mhalf);
    mulmul_(iyf, iyf, iyf2, &isgyf, &isgyf, &isgyf2, &ikyf, &ikyf, &ikyf2, &
	    c__30, mhalf);
    mulmul_(izf, izf, izf2, &isgzf, &isgzf, &isgzf2, &ikzf, &ikzf, &ikzf2, &
	    c__30, mhalf);
    if (isgxf2 < 0 || isgyf2 < 0 || isgzf2 < 0) {
	s_stop("5205", (ftnlen)4);
    }

    frterm_(&ixsew, &iysew, &izsew, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw2, ix2, iy2, iz2, iq2, &isgw2, &isgx2, &
	    isgy2, &isgz2, &isgq2, &ikw2, &ikx2, &iky2, &ikz2, &ikq2, mhalf, 
	    mfull, &ixse2, &iyse2, &izse2, &isclp[1]);

    frterm_(&ixfiw, &iyfiw, &izfiw, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw5, ix5, iy5, iz5, iq5, &isgw5, &isgx5, &
	    isgy5, &isgz5, &isgq5, &ikw5, &ikx5, &iky5, &ikz5, &ikq5, mhalf, 
	    mfull, &ixfi2, &iyfi2, &izfi2, &isclp[1]);

    detrm1_(iq5, ix2, iy2, iz2, ix2, iy2, iz2, &isgq5, &isgx2, &isgy2, &isgz2,
	     &isgx2, &isgy2, &isgz2, &ikq5, &ikx2, &iky2, &ikz2, &ikx2, &iky2,
	     &ikz2, iu, &isgu, &iku, mhalf);

    detrm1_(iq2, ix2, iy2, iz2, ix5, iy5, iz5, &isgq2, &isgx2, &isgy2, &isgz2,
	     &isgx5, &isgy5, &isgz5, &ikq2, &ikx2, &iky2, &ikz2, &ikx5, &iky5,
	     &ikz5, iv, &isgv, &ikv, mhalf);

    muldif_(iu, iv, io, &isgu, &isgv, &isgo, &iku, &ikv, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* iqsig1_ */

/* IQSIG2 */

/*     subroutine for determining position of point ifif with respect */
/*     to circle determined by (weighted) points ifir, isec, ithi, */
/*     if positive then ifif is outside the circle */
/*     if negative then ifif is inside the circle */
/*     if zero then ifif is in the circle */

/* Subroutine */ int iqsig2_(integer *x, integer *y, integer *z__, integer *w,
	 integer *x2, integer *y2, integer *z2, integer *w2, integer *ifir, 
	integer *isec, integer *ithi, integer *ifif, integer *mhalf, integer *
	mfull, integer *isclp, integer *isclw, integer *isclr, logical *
	delaun, integer *ipout)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    integer io[30], ip[30], iu[30], iv[30], iq2[30], iq3[30], iq5[30], iw2[30]
	    , ix2[30], iy2[30], iz2[30], ix3[30], iy3[30], iz3[30], ix5[30], 
	    iy5[30], iz5[30], iw3[30], iw5[30], iko, ikp, iwf[30], ixf[30], 
	    iyf[30], izf[30], iku, ikv, ikq2, ikq3, ixf2[30], iyf2[30], izf2[
	    30], ikq5, ikw2, ikw3, ikx2, ikw5, iky2, ikz2, ikx3, iky3, ikz3, 
	    ikx5, iky5, ikz5, ikcl, ikwf, ikxf, ikyf, ikzf, isgo, isgp, isgu, 
	    isgv, iwfi2, ixfi2, iyfi2, izfi2, ikxf2, isgq2, isgq3, ikyf2, 
	    isgq5, iwse2, ixse2, iyse2, iwfu2, iwth2, ixth2, iyth2, izth2, 
	    ixfu2, iyfu2, izfu2, izse2, isgw2, isgw3, isgw5, ikzf2, isgx2, 
	    isgy2, isgz2, isgx3, isgy3, isgz3, isgx5, isgy5, isgz5, isgcl, 
	    isgxf, isgyf, isgzf, isgwf, iwfiw, ixfiw, iyfiw, izfiw, iwsew, 
	    ixsew, iwfuw, iwthw, ixthw, iythw, izthw, ixfuw, iyfuw, izfuw, 
	    iysew, izsew;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), detrm0_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *), detrm2_(integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *),
	     detrm3_(integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);
    integer isgxf2, isgyf2, isgzf2;
    extern /* Subroutine */ int muldif_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), frterm_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *), mulmul_(integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclr;
    --isclw;
    --isclp;
    --w2;
    --z2;
    --y2;
    --x2;
    --w;
    --z__;
    --y;
    --x;

    /* Function Body */
    if (*delaun) {
	isgw2 = 0;
	isgw3 = 0;
	isgw5 = 0;
    } else {
	iwfuw = w[*ifir];
	iwsew = w[*isec];
	iwthw = w[*ithi];
	iwfiw = w[*ifif];

	iwfu2 = w2[*ifir];
	iwse2 = w2[*isec];
	iwth2 = w2[*ithi];
	iwfi2 = w2[*ifif];

	decmp2_(iwf, &isgwf, &ikwf, &iwfuw, &iwfu2, mhalf, mfull, &isclw[1]);
	isgcl = 1;
	ikcl = 2;
	decmp2_(io, &isgo, &iko, &iwsew, &iwse2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw2, &isgu, &isgcl, &isgw2, &iku, &ikcl, &ikw2,
		 &c__30, mhalf);
	decmp2_(io, &isgo, &iko, &iwthw, &iwth2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw3, &isgu, &isgcl, &isgw3, &iku, &ikcl, &ikw3,
		 &c__30, mhalf);
	decmp2_(io, &isgo, &iko, &iwfiw, &iwfi2, mhalf, mfull, &isclw[1]);
	muldif_(io, iwf, iu, &isgo, &isgwf, &isgu, &iko, &ikwf, &iku, &c__30, 
		mhalf);
	mulmul_(iu, &isclr[1], iw5, &isgu, &isgcl, &isgw5, &iku, &ikcl, &ikw5,
		 &c__30, mhalf);
    }

    ixfuw = x[*ifir];
    iyfuw = y[*ifir];
    izfuw = z__[*ifir];
    ixsew = x[*isec];
    iysew = y[*isec];
    izsew = z__[*isec];
    ixthw = x[*ithi];
    iythw = y[*ithi];
    izthw = z__[*ithi];
    ixfiw = x[*ifif];
    iyfiw = y[*ifif];
    izfiw = z__[*ifif];

    ixfu2 = x2[*ifir];
    iyfu2 = y2[*ifir];
    izfu2 = z2[*ifir];
    ixse2 = x2[*isec];
    iyse2 = y2[*isec];
    izse2 = z2[*isec];
    ixth2 = x2[*ithi];
    iyth2 = y2[*ithi];
    izth2 = z2[*ithi];
    ixfi2 = x2[*ifif];
    iyfi2 = y2[*ifif];
    izfi2 = z2[*ifif];

    decmp2_(ixf, &isgxf, &ikxf, &ixfuw, &ixfu2, mhalf, mfull, &isclp[1]);
    decmp2_(iyf, &isgyf, &ikyf, &iyfuw, &iyfu2, mhalf, mfull, &isclp[1]);
    decmp2_(izf, &isgzf, &ikzf, &izfuw, &izfu2, mhalf, mfull, &isclp[1]);
    mulmul_(ixf, ixf, ixf2, &isgxf, &isgxf, &isgxf2, &ikxf, &ikxf, &ikxf2, &
	    c__30, mhalf);
    mulmul_(iyf, iyf, iyf2, &isgyf, &isgyf, &isgyf2, &ikyf, &ikyf, &ikyf2, &
	    c__30, mhalf);
    mulmul_(izf, izf, izf2, &isgzf, &isgzf, &isgzf2, &ikzf, &ikzf, &ikzf2, &
	    c__30, mhalf);
    if (isgxf2 < 0 || isgyf2 < 0 || isgzf2 < 0) {
	s_stop("5305", (ftnlen)4);
    }

    frterm_(&ixsew, &iysew, &izsew, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw2, ix2, iy2, iz2, iq2, &isgw2, &isgx2, &
	    isgy2, &isgz2, &isgq2, &ikw2, &ikx2, &iky2, &ikz2, &ikq2, mhalf, 
	    mfull, &ixse2, &iyse2, &izse2, &isclp[1]);

    frterm_(&ixthw, &iythw, &izthw, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw3, ix3, iy3, iz3, iq3, &isgw3, &isgx3, &
	    isgy3, &isgz3, &isgq3, &ikw3, &ikx3, &iky3, &ikz3, &ikq3, mhalf, 
	    mfull, &ixth2, &iyth2, &izth2, &isclp[1]);

    frterm_(&ixfiw, &iyfiw, &izfiw, ixf, iyf, izf, &isgxf, &isgyf, &isgzf, &
	    ikxf, &ikyf, &ikzf, ixf2, iyf2, izf2, &isgxf2, &isgyf2, &isgzf2, &
	    ikxf2, &ikyf2, &ikzf2, iw5, ix5, iy5, iz5, iq5, &isgw5, &isgx5, &
	    isgy5, &isgz5, &isgq5, &ikw5, &ikx5, &iky5, &ikz5, &ikq5, mhalf, 
	    mfull, &ixfi2, &iyfi2, &izfi2, &isclp[1]);

    detrm0_(iq5, iy2, iz2, iy3, iz3, &isgq5, &isgy2, &isgz2, &isgy3, &isgz3, &
	    ikq5, &iky2, &ikz2, &iky3, &ikz3, iv, &isgv, &ikv, mhalf);

    detrm0_(iq5, iz2, ix2, iz3, ix3, &isgq5, &isgz2, &isgx2, &isgz3, &isgx3, &
	    ikq5, &ikz2, &ikx2, &ikz3, &ikx3, iu, &isgu, &iku, mhalf);

    detrm0_(iq5, ix2, iy2, ix3, iy3, &isgq5, &isgx2, &isgy2, &isgx3, &isgy3, &
	    ikq5, &ikx2, &iky2, &ikx3, &iky3, ip, &isgp, &ikp, mhalf);

    detrm3_(iv, ix2, ix3, &isgv, &isgx2, &isgx3, iu, iy2, iy3, &isgu, &isgy2, 
	    &isgy3, ip, iz2, iz3, &isgp, &isgz2, &isgz3, &ikv, &iku, &ikp, &
	    ikx2, &iky2, &ikz2, &ikx3, &iky3, &ikz3, io, &isgo, &iko, mhalf);

    detrm2_(iq2, ix2, iy2, iz2, &isgq2, &isgx2, &isgy2, &isgz2, iq3, ix3, iy3,
	     iz3, &isgq3, &isgx3, &isgy3, &isgz3, &ikq2, &ikx2, &iky2, &ikz2, 
	    &ikq3, &ikx3, &iky3, &ikz3, iu, &isgu, &iku, mhalf);
    mulmul_(iu, ix5, ip, &isgu, &isgx5, &isgp, &iku, &ikx5, &ikp, &c__30, 
	    mhalf);
    muldif_(io, ip, iv, &isgo, &isgp, &isgv, &iko, &ikp, &ikv, &c__30, mhalf);

    detrm2_(iq2, iy2, iz2, ix2, &isgq2, &isgy2, &isgz2, &isgx2, iq3, iy3, iz3,
	     ix3, &isgq3, &isgy3, &isgz3, &isgx3, &ikq2, &iky2, &ikz2, &ikx2, 
	    &ikq3, &iky3, &ikz3, &ikx3, iu, &isgu, &iku, mhalf);
    mulmul_(iu, iy5, io, &isgu, &isgy5, &isgo, &iku, &iky5, &iko, &c__30, 
	    mhalf);
    muldif_(iv, io, ip, &isgv, &isgo, &isgp, &ikv, &iko, &ikp, &c__30, mhalf);

    detrm2_(iq2, iz2, ix2, iy2, &isgq2, &isgz2, &isgx2, &isgy2, iq3, iz3, ix3,
	     iy3, &isgq3, &isgz3, &isgx3, &isgy3, &ikq2, &ikz2, &ikx2, &iky2, 
	    &ikq3, &ikz3, &ikx3, &iky3, iu, &isgu, &iku, mhalf);
    mulmul_(iu, iz5, iv, &isgu, &isgz5, &isgv, &iku, &ikz5, &ikv, &c__30, 
	    mhalf);
    muldif_(ip, iv, io, &isgp, &isgv, &isgo, &ikp, &ikv, &iko, &c__30, mhalf);

    *ipout = isgo;

    return 0;
} /* iqsig2_ */

/* FRTERM */

/* Subroutine */ int frterm_(integer *ixsew, integer *iysew, integer *izsew, 
	integer *ixf, integer *iyf, integer *izf, integer *isgxf, integer *
	isgyf, integer *isgzf, integer *ikxf, integer *ikyf, integer *ikzf, 
	integer *ixf2, integer *iyf2, integer *izf2, integer *isgxf2, integer 
	*isgyf2, integer *isgzf2, integer *ikxf2, integer *ikyf2, integer *
	ikzf2, integer *iw2, integer *ix2, integer *iy2, integer *iz2, 
	integer *iq2, integer *isgw2, integer *isgx2, integer *isgy2, integer 
	*isgz2, integer *isgq2, integer *ikw2, integer *ikx2, integer *iky2, 
	integer *ikz2, integer *ikq2, integer *mhalf, integer *mfull, integer 
	*ixse2, integer *iyse2, integer *izse2, integer *isclp)
{
    integer io[30], ip[30], iu[30], iv[30], iko, ikp, iku, ikv, isgo, isgp, 
	    isgu, isgv;
    extern /* Subroutine */ int decmp2_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), muldif_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *), mulmul_(
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *);



    /* Parameter adjustments */
    --isclp;
    --iq2;
    --iz2;
    --iy2;
    --ix2;
    --iw2;
    --izf2;
    --iyf2;
    --ixf2;
    --izf;
    --iyf;
    --ixf;

    /* Function Body */
    decmp2_(io, &isgo, &iko, ixsew, ixse2, mhalf, mfull, &isclp[1]);
    muldif_(io, &ixf[1], &ix2[1], &isgo, isgxf, isgx2, &iko, ikxf, ikx2, &
	    c__30, mhalf);
    mulmul_(io, io, iu, &isgo, &isgo, &isgu, &iko, &iko, &iku, &c__30, mhalf);
    muldif_(iu, &ixf2[1], iv, &isgu, isgxf2, &isgv, &iku, ikxf2, &ikv, &c__30,
	     mhalf);
    muldif_(iv, &iw2[1], ip, &isgv, isgw2, &isgp, &ikv, ikw2, &ikp, &c__30, 
	    mhalf);

    decmp2_(io, &isgo, &iko, iysew, iyse2, mhalf, mfull, &isclp[1]);
    muldif_(io, &iyf[1], &iy2[1], &isgo, isgyf, isgy2, &iko, ikyf, iky2, &
	    c__30, mhalf);
    mulmul_(io, io, iu, &isgo, &isgo, &isgu, &iko, &iko, &iku, &c__30, mhalf);
    muldif_(iu, &iyf2[1], iv, &isgu, isgyf2, &isgv, &iku, ikyf2, &ikv, &c__30,
	     mhalf);
    isgv = -isgv;
    muldif_(ip, iv, iu, &isgp, &isgv, &isgu, &ikp, &ikv, &iku, &c__30, mhalf);

    decmp2_(io, &isgo, &iko, izsew, izse2, mhalf, mfull, &isclp[1]);
    muldif_(io, &izf[1], &iz2[1], &isgo, isgzf, isgz2, &iko, ikzf, ikz2, &
	    c__30, mhalf);
    mulmul_(io, io, iv, &isgo, &isgo, &isgv, &iko, &iko, &ikv, &c__30, mhalf);
    muldif_(iv, &izf2[1], ip, &isgv, isgzf2, &isgp, &ikv, ikzf2, &ikp, &c__30,
	     mhalf);
    isgp = -isgp;
    muldif_(iu, ip, &iq2[1], &isgu, &isgp, isgq2, &iku, &ikp, ikq2, &c__30, 
	    mhalf);

    return 0;
} /* frterm_ */

/* DETRM0 */

/* Subroutine */ int detrm0_(integer *iq, integer *ix2, integer *iy2, integer 
	*ix3, integer *iy3, integer *isgq, integer *isgx2, integer *isgy2, 
	integer *isgx3, integer *isgy3, integer *ikq, integer *ikx2, integer *
	iky2, integer *ikx3, integer *iky3, integer *io, integer *isgo, 
	integer *iko, integer *mhalf)
{
    integer iu[30], iv[30], iw[30], iku, ikv, ikw, isgu, isgv, isgw;
    extern /* Subroutine */ int muldif_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), mulmul_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);



    /* Parameter adjustments */
    --io;
    --iy3;
    --ix3;
    --iy2;
    --ix2;
    --iq;

    /* Function Body */
    mulmul_(&iq[1], &ix2[1], iv, isgq, isgx2, &isgv, ikq, ikx2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, &iy3[1], iu, &isgv, isgy3, &isgu, &ikv, iky3, &iku, &c__30, 
	    mhalf);
    mulmul_(&iq[1], &iy2[1], iv, isgq, isgy2, &isgv, ikq, iky2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, &ix3[1], iw, &isgv, isgx3, &isgw, &ikv, ikx3, &ikw, &c__30, 
	    mhalf);
    muldif_(iu, iw, &io[1], &isgu, &isgw, isgo, &iku, &ikw, iko, &c__30, 
	    mhalf);

    return 0;
} /* detrm0_ */

/* DETRM1 */

/* Subroutine */ int detrm1_(integer *iq, integer *ix2, integer *iy2, integer 
	*iz2, integer *ix3, integer *iy3, integer *iz3, integer *isgq, 
	integer *isgx2, integer *isgy2, integer *isgz2, integer *isgx3, 
	integer *isgy3, integer *isgz3, integer *ikq, integer *ikx2, integer *
	iky2, integer *ikz2, integer *ikx3, integer *iky3, integer *ikz3, 
	integer *io, integer *isgo, integer *iko, integer *mhalf)
{
    integer iv[30], iw[30], ikv, ikw, isgv, isgw;
    extern /* Subroutine */ int muldif_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), mulmul_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);



    /* Parameter adjustments */
    --io;
    --iz3;
    --iy3;
    --ix3;
    --iz2;
    --iy2;
    --ix2;
    --iq;

    /* Function Body */
    mulmul_(&iq[1], &ix2[1], iv, isgq, isgx2, &isgv, ikq, ikx2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, &ix3[1], &io[1], &isgv, isgx3, isgo, &ikv, ikx3, iko, &c__30, 
	    mhalf);

    mulmul_(&iq[1], &iy2[1], iv, isgq, isgy2, &isgv, ikq, iky2, &ikv, &c__30, 
	    mhalf);
    mulmul_(iv, &iy3[1], iw, &isgv, isgy3, &isgw, &ikv, iky3, &ikw, &c__30, 
	    mhalf);
    isgw = -isgw;
    muldif_(&io[1], iw, iv, isgo, &isgw, &isgv, iko, &ikw, &ikv, &c__30, 
	    mhalf);

    mulmul_(&iq[1], &iz2[1], &io[1], isgq, isgz2, isgo, ikq, ikz2, iko, &
	    c__30, mhalf);
    mulmul_(&io[1], &iz3[1], iw, isgo, isgz3, &isgw, iko, ikz3, &ikw, &c__30, 
	    mhalf);
    isgw = -isgw;
    muldif_(iv, iw, &io[1], &isgv, &isgw, isgo, &ikv, &ikw, iko, &c__30, 
	    mhalf);

    return 0;
} /* detrm1_ */

/* DETRM2 */

/* Subroutine */ int detrm2_(integer *iq2, integer *ix2, integer *iy2, 
	integer *iz2, integer *isgq2, integer *isgx2, integer *isgy2, integer 
	*isgz2, integer *iq3, integer *ix3, integer *iy3, integer *iz3, 
	integer *isgq3, integer *isgx3, integer *isgy3, integer *isgz3, 
	integer *ikq2, integer *ikx2, integer *iky2, integer *ikz2, integer *
	ikq3, integer *ikx3, integer *iky3, integer *ikz3, integer *io, 
	integer *isgo, integer *iko, integer *mhalf)
{
    integer ip[30], ir[30], iv[30], iw[30], ikp, ikr, ikv, ikw, isgp, isgr, 
	    isgv, isgw;
    extern /* Subroutine */ int detrm0_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), muldif_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), mulmul_(integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);



    /* Parameter adjustments */
    --io;
    --iz3;
    --iy3;
    --ix3;
    --iq3;
    --iz2;
    --iy2;
    --ix2;
    --iq2;

    /* Function Body */
    detrm0_(&iq2[1], &ix2[1], &iy2[1], &ix3[1], &iy3[1], isgq2, isgx2, isgy2, 
	    isgx3, isgy3, ikq2, ikx2, iky2, ikx3, iky3, iv, &isgv, &ikv, 
	    mhalf);
    mulmul_(iv, &iy3[1], &io[1], &isgv, isgy3, isgo, &ikv, iky3, iko, &c__30, 
	    mhalf);

    detrm0_(&iq2[1], &iz2[1], &ix2[1], &iz3[1], &ix3[1], isgq2, isgz2, isgx2, 
	    isgz3, isgx3, ikq2, ikz2, ikx2, ikz3, ikx3, iv, &isgv, &ikv, 
	    mhalf);
    mulmul_(iv, &iz3[1], ip, &isgv, isgz3, &isgp, &ikv, ikz3, &ikp, &c__30, 
	    mhalf);

    muldif_(&io[1], ip, ir, isgo, &isgp, &isgr, iko, &ikp, &ikr, &c__30, 
	    mhalf);

    detrm0_(&iq3[1], &ix2[1], &iy2[1], &ix3[1], &iy3[1], isgq3, isgx2, isgy2, 
	    isgx3, isgy3, ikq3, ikx2, iky2, ikx3, iky3, iv, &isgv, &ikv, 
	    mhalf);
    mulmul_(iv, &iy2[1], &io[1], &isgv, isgy2, isgo, &ikv, iky2, iko, &c__30, 
	    mhalf);

    detrm0_(&iq3[1], &iz2[1], &ix2[1], &iz3[1], &ix3[1], isgq3, isgz2, isgx2, 
	    isgz3, isgx3, ikq3, ikz2, ikx2, ikz3, ikx3, iv, &isgv, &ikv, 
	    mhalf);
    mulmul_(iv, &iz2[1], ip, &isgv, isgz2, &isgp, &ikv, ikz2, &ikp, &c__30, 
	    mhalf);

    muldif_(&io[1], ip, iw, isgo, &isgp, &isgw, iko, &ikp, &ikw, &c__30, 
	    mhalf);
    muldif_(ir, iw, &io[1], &isgr, &isgw, isgo, &ikr, &ikw, iko, &c__30, 
	    mhalf);

    return 0;
} /* detrm2_ */

/* DETRM3 */

/* Subroutine */ int detrm3_(integer *ix2, integer *iy2, integer *iz2, 
	integer *isgx2, integer *isgy2, integer *isgz2, integer *ix3, integer 
	*iy3, integer *iz3, integer *isgx3, integer *isgy3, integer *isgz3, 
	integer *ix4, integer *iy4, integer *iz4, integer *isgx4, integer *
	isgy4, integer *isgz4, integer *ikx2, integer *ikx3, integer *ikx4, 
	integer *iky2, integer *iky3, integer *iky4, integer *ikz2, integer *
	ikz3, integer *ikz4, integer *io, integer *isgo, integer *iko, 
	integer *mhalf)
{
    integer iu[30], iv[30], iw[30], iku, ikv, ikw, isgu, isgv, isgw;
    extern /* Subroutine */ int muldif_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *), mulmul_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *);



    /* Parameter adjustments */
    --io;
    --iz4;
    --iy4;
    --ix4;
    --iz3;
    --iy3;
    --ix3;
    --iz2;
    --iy2;
    --ix2;

    /* Function Body */
    mulmul_(&ix3[1], &iy4[1], iv, isgx3, isgy4, &isgv, ikx3, iky4, &ikv, &
	    c__30, mhalf);
    mulmul_(&ix4[1], &iy3[1], iu, isgx4, isgy3, &isgu, ikx4, iky3, &iku, &
	    c__30, mhalf);
    muldif_(iv, iu, iw, &isgv, &isgu, &isgw, &ikv, &iku, &ikw, &c__30, mhalf);
    mulmul_(iw, &iz2[1], &io[1], &isgw, isgz2, isgo, &ikw, ikz2, iko, &c__30, 
	    mhalf);

    mulmul_(&ix2[1], &iy4[1], iv, isgx2, isgy4, &isgv, ikx2, iky4, &ikv, &
	    c__30, mhalf);
    mulmul_(&ix4[1], &iy2[1], iu, isgx4, isgy2, &isgu, ikx4, iky2, &iku, &
	    c__30, mhalf);
    muldif_(iv, iu, iw, &isgv, &isgu, &isgw, &ikv, &iku, &ikw, &c__30, mhalf);
    mulmul_(iw, &iz3[1], iu, &isgw, isgz3, &isgu, &ikw, ikz3, &iku, &c__30, 
	    mhalf);
    muldif_(&io[1], iu, iw, isgo, &isgu, &isgw, iko, &iku, &ikw, &c__30, 
	    mhalf);

    mulmul_(&ix3[1], &iy2[1], iv, isgx3, isgy2, &isgv, ikx3, iky2, &ikv, &
	    c__30, mhalf);
    mulmul_(&ix2[1], &iy3[1], iu, isgx2, isgy3, &isgu, ikx2, iky3, &iku, &
	    c__30, mhalf);
    muldif_(iv, iu, &io[1], &isgv, &isgu, isgo, &ikv, &iku, iko, &c__30, 
	    mhalf);
    mulmul_(&io[1], &iz4[1], iu, isgo, isgz4, &isgu, iko, ikz4, &iku, &c__30, 
	    mhalf);
    muldif_(iw, iu, &io[1], &isgw, &isgu, isgo, &ikw, &iku, iko, &c__30, 
	    mhalf);

    return 0;
} /* detrm3_ */

#ifdef __cplusplus
	}
#endif
