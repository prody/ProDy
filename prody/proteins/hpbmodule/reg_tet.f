!c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
!c
!c     Program regtet
!c
!c     This program will compute a Regular tetrahedralization for
!c     a set of points in 3-dimensional space.
!c
!c     A Regular tetrahedralization is the dual of a Power diagram.
!c     It is essentially a Delaunay tetrahedralization with weights.
!c
!c     If no weights are present the program will compute a Delaunay
!c     tetrahedralization. It is noted that if weights are present
!c     and degeneracies exist there may be points whose Power cells
!c     are not 3-dimensional.
!c
!c     The program is based on an algorithm for constructing Regular
!c     tetrahedralizations with incremental topological flipping.
!c
!c     Computations in this program for the purpose of computing
!c     the tetrahedralization are done in exact arithmetic whenever
!c     floating point arithmetic (done in double precision) does not
!c     seem appropriate.
!c
!c     Subroutine regtet is the driver routine of the program and is
!c     called from the main routine.
!c
!c     Documentation appears below in driver routine regtet.
!c
!c     Author: Javier Bernal
!c             National Institute of Standards and Technology (NIST)
!c
!c     Disclaimer:
!c
!c     This software was developed at the National Institute of Standards
!c     and Technology by employees of the Federal Government in the
!c     course of their official duties. Pursuant to title 17 Section 105
!c     of the United States Code this software is not subject to
!c     copyright protection and is in the public domain. This software is
!c     experimental. NIST assumes no responsibility whatsoever for its
!c     use by other parties, and makes no guarantees, expressed or
!c     implied, about its quality, reliability, or any other
!c     characteristic. We would appreciate acknowledgement if the
!c     software is used.
!c
!*MAIN
!c
!      program main
!c
!c     Setting of parameters:
!c
!c     1. Flipping history used for locating points:
!c
!c     integer nmax, nvmax, nhmax
!c     parameter (nmax=150000, nvmax=55*nmax, nhmax=1500)
!c
!c     2. Flipping history not used for locating points:
!c
!      integer nmax, nvmax, nhmax
!      parameter (nmax=150000, nvmax= 7*nmax, nhmax=1500)
!c
!      double precision x(nmax), y(nmax), z(nmax), w(nmax)
!      real v(nmax)
!      integer ix(nmax), iy(nmax), iz(nmax), iw(nmax)
!      integer ix2(nmax), iy2(nmax), iz2(nmax), iw2(nmax)
!      integer icon(8,nvmax), is(nmax), ifl(nvmax), io(nmax)
!      integer id(nvmax), ih(nhmax)
!      integer nv, nw, nt, nd, naddl, isu, jsu, ksu, nsu, icfig, iwfig
!      double precision wlenx, wleny, wlenz, wlenw, epz
!      logical delaun, pntoff, flphis, artfcl, random, reccor, redchk
!c
!      logical prompt, bigbox
!      character*1 answ
!      double precision xcor, ycor, zcor, wght
!      integer np, i, j, iric, irec
!      integer ideli, ipnti, iflpi, iarti, irani, ireci, iredi
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      write(*,*)'This program is for computing a Regular ',
!     *          'tetrahedralization for a set'
!      write(*,*)'of points in 3-d space, ',
!     *          'i. e. a Delaunay tetrahedralization with weights.'
!      write(*,*)' '
!      write(*,*)'If no weights are present a Delaunay ',
!     *          'tetrahedralization is computed.'
!      write(*,*)'Note that if weights are present and degeneracies ',
!     *          'exist there may'
!      write(*,*)'be points whose Power cells are ',
!     *          'not 3-dimensional.'
!      write(*,*)' '
!      write(*,*)'Computations in this program are done in exact ',
!     *          'arithmetic whenever'
!      write(*,*)'floating point arithmetic (done in double precision) ',
!     *          'does not seem'
!      write(*,*)'appropriate.'
!      write(*,*)' '
!      write(*,*)'Documentation about input/output variables, etc. ',
!     *          'appears in driver'
!      write(*,*)'routine regtet'
!c
!c----------------------------------------------------------------------
!c
!c     inquire about existence of prompts file
!c
!      inquire(file = 'prompts.data', exist = prompt)
!      write(*,*)' '
!      if(prompt)then
!         write(*,*)'User intervention not required: responses to ',
!     *             'would-be program generated'
!         write(*,*)'prompts will be read from file prompts.data ',
!     *             '(prompts will be suppressed).'
!         open (9, file = 'prompts.data')
!      else
!         write(*,*)'User intervention required: file named ',
!     *             'prompts.data does not exist, thus'
!         write(*,*)'program will generate prompts for which ',
!     *             'user must provide responses (in'
!         write(*,*)'order to avoid user intervention a file ',
!     *             'named prompts.data must exist that'
!         write(*,*)'contains the responses to the would-be ',
!     *             'prompts; prompts are then suppressed;'
!         write(*,*)'for this purpose the program will generate ',
!     *             'a file named prompts.tentative'
!         write(*,*)'that will contain the responses provided ',
!     *             'by the user during the current run,'
!         write(*,*)'and which can then be renamed prompts.data ',
!     *             ' by the user for future runs).'
!         open (10, file = 'prompts.tentative')
!      endif
!c
!c----------------------------------------------------------------------
!c
!   10 format (a1)
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'A Delaunay tetrahedralization (no weights) to be ',
!     *             'computed?(y/n)'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         delaun=.true.
!         write(*,*)'A Delaunay tetrahedralization will be computed.'
!      else
!         delaun=.false.
!         write(*,*)'A Regular tetrahedralization will be computed.'
!      endif
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'Any input points to be inactive during ',
!     *             'tetrahedralization computation?(y/n)'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         pntoff =.true.
!         write(*,*)'Some input points will be inactive during ',
!     *             'tetrahedralization computation.'
!      else
!         pntoff =.false.
!         write(*,*)'All input points will be active during ',
!     *             'tetrahedralization computation.'
!      endif
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'Tetrahedron list with flipping history to be ',
!     *             'used for locating points?(y/n)'
!         write(*,*)'(Otherwise shishkebab method will be used. Note ',
!     *             'that in order'
!         write(*,*)'to use flipping history, parameter nvmax in the ',
!     *             'program should'
!         write(*,*)'be set to about 55*nmax, otherwise to about ',
!     *             '7*nmax).'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         flphis=.true.
!         write(*,*)'Tetrahedron list with flipping history will be ',
!     *             'used.'
!      else
!         flphis=.false.
!         write(*,*)'Shishkebab method will be used.'
!      endif
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'Do you want output tetrahedron list to include ',
!     *             'both real and artificial'
!         write(*,*)'tetrahedra in the final tetrahedralization?(y/n)'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         artfcl=.true.
!         write(*,*)'Output tetrahedron list will include both real ',
!     *             'and artificial tetrahedra'
!         write(*,*)'in the final tetrahedralization together with ',
!     *             'flipping history tetrahedra'
!         write(*,*)'if flipping history was used for locating points.'
!      else
!         artfcl=.false.
!         write(*,*)'Output tetrahedron list will only include ',
!     *             'real tetrahedra in the final'
!         write(*,*)'tetrahedralization.'
!      endif
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'Do you want input points to be inserted in ',
!     *             'a random order?(y/n)'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         random = .true.
!         write(*,*)'Input points will be inserted in a random order.'
!         write(*,*)' '
!         if(prompt)then
!            read(9,*) isu, jsu, ksu, nsu
!            write(*,*)'The four seeds for randomizing are:'
!            write(*,*)isu, jsu, ksu, nsu
!         else
!            write(*,*)'Enter seeds for randomizing (4 integers):'
!            read(5,*) isu, jsu, ksu, nsu
!c
!c           isu = 521288629
!c           jsu = 362436069
!c           ksu = 16163801
!c           nsu = 131199299
!c
!            write(10,*) isu, jsu, ksu, nsu
!         endif
!      else
!         random = .false.
!         write(*,*)'Input points will be inserted in their ',
!     *             'current order.'
!      endif
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'Do you want points that define a rectangular ',
!     *             'regular grid on the surface'
!         write(*,*)'of a rectangular polyhedron that contains ',
!     *             'set of input points to become,'
!         write(*,*)'together with set of input points and ',
!     *             ' artificial points, the set for'
!         write(*,*)'which a tetrahedralization is to be computed?(y/n)'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         reccor=.true.
!         write(*,*)'Points that define a rectangular regular grid on ',
!     *             'the surface of a'
!         write(*,*)'rectangular polyhedron that contains set of input ',
!     *             'points will be part'
!         write(*,*)'of the set for which tetrahedralization is to be ',
!     *             'computed; dimensions'
!         write(*,*)'of polyhedron and choice of points that define ',
!     *             'grid are according to'
!         write(*,*)'specifications provided by user.'
!         write(*,*)' '
!         write(*,*)'If xmax, ymax, zmax are the maximum values of the ',
!     *             'x, y, and z coordinates'
!         write(*,*)'of the input points, respectively, and xmin, ',
!     *             'ymin, zmin are the minimum'
!         write(*,*)'values, then for positive numbers ',
!     *             'wlenx, wleny, wlenz (provided by user),'
!         write(*,*)'the eight vertices of the the polyhedron will be: '
!         write(*,*)'(xmin-wlenx, ymin-wleny, zmin-wlenz), ',
!     *             '(xmax+wlenx, ymin-wleny, zmin-wlenz),'
!         write(*,*)'(xmax+wlenx, ymax+wleny  zmin-wlenz), ',
!     *             '(xmin-wlenx, ymax+wleny, zmin-wlenz),'
!         write(*,*)'(xmin-wlenx, ymin-wleny, zmax+wlenz), ',
!     *             '(xmax+wlenx, ymin-wleny, zmax+wlenz),'
!         write(*,*)'(xmax+wlenx, ymax+wleny  zmax+wlenz), ',
!     *             '(xmin-wlenx, ymax+wleny, zmax+wlenz).'
!         write(*,*)' '
!         write(*,*)'For positive integer naddl (provided by user) ',
!     *             'for each facet of the'
!         write(*,*)'polyhedron a set of naddl x naddl points is ',
!     *             'generated by program;'
!         write(*,*)'this set defines a rectangular regular grid ',
!     *             'on the facet and contains'
!         write(*,*)'the four vertices of the facet; the points in ',
!     *             'the union of the six'
!         write(*,*)'sets thus generated define the rectangular ',
!     *             'grid on the surface of'
!         write(*,*)'the polyhedron; naddl can not be less than 2; ',
!     *             'if it equals 2 then'
!         write(*,*)'the grid is defined exactly by the 8 vertices ',
!     *             'of the polyhedron.'
!         if(.not.delaun) then
!            write(*,*)' '
!            write(*,*)'If wmin is the minimum value of the weights of ',
!     *                'the input points then'
!            write(*,*)'for a real number wlenw (provided ',
!     *                'by user) a weight equal to wmin - wlenw'
!            write(*,*)'is assigned by the program to ',
!     *                'each point in the rectangular grid on the'
!            write(*,*)'surface of the polyhedron.'
!         endif
!         write(*,*)' '
!         if(prompt)then
!            read(9,*) wlenx, wleny, wlenz
!            if(.not.delaun) read(9,*) wlenw
!            read(9,*) naddl
!            write(*,*)'The values of wlenx, wleny, wlenz are:'
!            write(*,*) wlenx, wleny, wlenz
!            if(.not.delaun) then
!               write(*,*)'The value of wlenw is:'
!               write(*,*) wlenw
!            endif
!            write(*,*)'The value of naddl is:'
!            write(*,*) naddl
!         else
!            write(*,*)'Enter wlenx, wleny, wlenz ',
!     *      '(3 positive real numbers):'
!            read(5,*) wlenx, wleny, wlenz
!            if(.not.delaun) then
!               write(*,*)'Enter wlenw (a real number):'
!               read(5,*) wlenw
!            endif
!            write(*,*)'Enter naddl (an integer greater than 1):'
!            read(5,*) naddl
!            write(10,*) wlenx, wleny, wlenz
!            if(.not.delaun) write(10,*) wlenw
!            write(10,*) naddl
!         endif
!      else
!         reccor=.false.
!         write(*,*)'Points that define a rectangular regular grid on ',
!     *             'the surface of a'
!         write(*,*)'rectangular polyhedron that ',
!     *             'contains set of input points will not'
!         write(*,*)'be part of the set for which a tetrahedralization ',
!     *             'is to be computed.'
!      endif
!c
!c----------------------------------------------------------------------
!c
!      if(delaun) go to 50
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'Do you want redundant points to be checked for ',
!     *             'redundancy after '
!         write(*,*)'final tetrahedralization has been computed?(y/n)'
!         write(*,*)'(doing so may require some additional time).'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         redchk=.true.
!         write(*,*)'Redundant points will be checked for redundancy.'
!      else
!         redchk=.false.
!         write(*,*)'Redundant points will not be checked for ',
!     *             'redundancy.'
!      endif
!   50 continue
!c
!c----------------------------------------------------------------------
!c
!      bigbox=.false.
!      if(.not.reccor) go to 55
!      write(*,*)' '
!      answ = ' '
!      if(prompt)then
!         read(9,10) answ
!      else
!         write(*,*)'Do you want points that define rectangular grid ',
!     *             'as described above'
!         write(*,*)'to be saved in a file?(y/n)'
!         read(5,10) answ
!         write(10,10) answ
!      endif
!      if(answ.eq.'y'.or.answ.eq.'Y') then
!         bigbox=.true.
!         write(*,*)'Rectangular grid will be saved in a file.'
!      else
!         write(*,*)'Rectangular grid will not be saved in a file.'
!      endif
!   55 continue
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      if(prompt)then
!         read(9,*) icfig
!         write(*,*)'The number of significant figures of decimal part ',
!     *             'of point coordinates is ', icfig
!      else
!         write(*,*)'Enter icfig, i. e. number of significant figures ',
!     *             'of decimal part'
!         write(*,*)'of point coordinates, -1 < icfig < 10 (total ',
!     *             'number of significant'
!         write(*,*)'figures should be at most 14 with at most 9 to ',
!     *             'either the left or'
!         write(*,*)'the right of the decimal point):'
!         read(5,*) icfig
!         write(10,*) icfig
!      endif
!      if(delaun) go to 60
!      write(*,*)' '
!      if(prompt)then
!         read(9,*) iwfig
!         write(*,*)'The number of significant figures of decimal part ',
!     *             'of point weights is ', iwfig
!      else
!         write(*,*)'Enter iwfig, i. e. number of significant figures ',
!     *             'of decimal part'
!         write(*,*)'of point weights, -1 < iwfig < 10, ',
!     *             '-1 < 2*icfig - iwfig < 10 (total'
!         write(*,*)'number of significant figures should ',
!     *             'be at most 14 with at most 9'
!         write(*,*)'to either the left or the right of the ',
!     *             'decimal point):'
!         read(5,*) iwfig
!         write(10,*) iwfig
!      endif
!   60 continue
!c
!c----------------------------------------------------------------------
!c
!c     open files
!c
!      open (11, file = 'pnts-wghts')
!      open (12, file = 'tetrahedra')
!c
!c----------------------------------------------------------------------
!c
!c     set tolerance
!c
!      epz = 0.001d0
!c
!c     read vertex data
!c
!      write(*,*)' '
!      write(*,*)'Reading vertex data ...'
!      nv = 0
!      if(delaun) go to 130
!c
!c     read vertex data with weights
!c
!  100 continue
!      read (11, *, end = 120) xcor, ycor, zcor, wght
!      nv = nv + 1
!      if (nv .gt. nmax) stop 10
!      x(nv) = xcor
!      y(nv) = ycor
!      z(nv) = zcor
!      w(nv) = wght
!      go to 100
!  120 continue
!      go to 140
!c
!c     read vertex data without weights
!c
!  130 continue
!      read (11, *, end = 140) xcor, ycor, zcor
!      nv = nv + 1
!      if (nv .gt. nmax) stop 20
!      x(nv) = xcor
!      y(nv) = ycor
!      z(nv) = zcor
!      w(nv) = 0.0d0
!      go to 130
!  140 continue
!c
!c     read off-on information about input points if there are any that
!c     are to be inactive during the tetrahedralization computation
!c
!      if(pntoff) then
!         open (13, file = 'points-off')
!         read (13,*) np
!         if(np .ne. nv) stop 30
!         read (13,150) (is(i), i = 1, np)
!      endif
!  150 format (40(1x,i1))
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      write(*,*)'Computation of tetrahedralization to begin: driver ',
!     *          'subroutine regtet'
!      write(*,*)'will be called from the ',
!     *          'main routine and points will be processed.'
!      write(*,*)' '
!      write(*,*)'Please wait ...'
!      write(*,*)' '
!c
!c----------------------------------------------------------------------
!c 
!c     call regtet to compute tetrahedralization
!c
!      call regtet(x, y, z, w, v, ix, iy, iz, iw, ix2, iy2, iz2, iw2,
!     *            icon, is, ifl, io, id, ih, nv, nw, nt, nd, nmax,
!     *            nvmax, nhmax, wlenx, wleny, wlenz, wlenw, naddl,
!     *            isu, jsu, ksu, nsu, icfig, iwfig, epz, delaun,
!     *            pntoff, flphis, artfcl, random, reccor, redchk)
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      write(*,*)'(Back to the main routine).'
!      write(*,*)' '
!      write(*,*)'Computation of tetrahedralization has been completed.'
!      write(*,*)' '
!      write(*,*)'Number of input vertices = ',nv
!      write(*,*)' '
!      write(*,*)'Number of output vertices = ',nw
!      write(*,*)' '
!      write(*,*)'Length of final tetrahedron list = ',nt
!      write(*,*)' '
!      write(*,*)'Number of real tetrahedra = ',nd
!      write(*,*)' '
!      write(*,*)'(The output vertices are the vertices of tetrahedra ',
!     *          'in the final'
!      write(*,*)'tetrahedron list).'
!c
!c----------------------------------------------------------------------
!c
!c     write tetrahedralization information in files
!c
!      ideli = 0
!      ipnti = 0
!      iflpi = 0
!      iarti = 0
!      irani = 0
!      ireci = 0
!      iredi = 0
!      if(delaun) ideli = 1
!      if(pntoff) ipnti = 1
!      if(flphis) iflpi = 1
!      if(artfcl) iarti = 1
!      if(random) irani = 1
!      if(reccor) ireci = 1
!      if(redchk) iredi = 1
!      write(*,*)' '
!      write(*,*)'Saving tetrahedralization data in output file ...'
!      if(nt.ne.0)then
!         write (12,480) ideli, ipnti, iflpi, iarti, irani, ireci, iredi
!         if(.not.delaun) then
!            write (12,*) nw, nt, icfig, iwfig
!         else
!            write (12,*) nw, nt, icfig
!         endif
!         write (12,500) ((icon(i,j), i = 1, 8), j = 1, nt)
!         write (12,550) (is(i), i = 1, nw)
!         if(reccor .and. .not.delaun) then
!            write (12,*) wlenx, wleny, wlenz, wlenw
!            write (12,*) naddl
!         elseif(reccor) then
!            write (12,*) wlenx, wleny, wlenz
!            write (12,*) naddl
!         endif
!      else
!         write (*,*)'warning: no real tetrahedra were created.'
!         write (*,*)' '
!         write (12,*)'warning: no real tetrahedra were created.'
!         write (12,*)' '
!      endif
!c
!  480 format (7(1x,i1))
!  500 format (8i10)
!  550 format (7i10)
!c
!      if(bigbox) then
!         open (14, file = 'bgbox-pnts')
!         iric = 1
!         if(artfcl) iric = 9
!         irec = iric + 6*(naddl**2) - 12*naddl + 7
!         if(delaun) then
!            do 580 i = iric, irec
!               write (14,*) x(i), y(i), z(i)
!  580       continue
!         else
!            do 590 i = iric, irec
!               write (14,*) x(i), y(i), z(i), w(i)
!  590       continue
!         endif
!      endif
!c
!c----------------------------------------------------------------------
!c
!c     write information in file in a readable form
!c
!c     open (23, file = 'exeinfor')
!c     write (23,*)' '
!c     write (23,*)'Number of vertices   =',nw
!c     write (23,*)'Number of tetrahedra =',nt
!c     write (23,*)' '
!c     if(nt.ne.0)then
!c        do 600 i = 1, nt
!c           write (23,800) i, (icon(j,i), j=1,8)
!c 600    continue
!c        do 650 i = 1, nw
!c           write (23,850) i, is(i)
!c 650    continue
!c     else
!c        write (23,*)'warning: no real tetrahedra were created.'
!c        write (23,*)' '
!c     endif
!c 800 format ('tetr: ', i4, ' info: ', 8(1x,i4))
!c 850 format ('vert: ', i4, ' is: ',i4)
!c
!c----------------------------------------------------------------------
!c
!      write(*,*)' '
!      stop
!      end
!*REGTET
c**********************************************************************
c
c     Driver subroutine of Fortran 77 program REGTET [3] for computing
c     a Regular tetrahedralization for a set of points in 3-dimensional
c     space.
c
c     A Regular tetrahedralization is the dual of a Power diagram.
c     It is essentially a Delaunay tetrahedralization with weights.
c     In the absence of weights the program simply computes a Delaunay
c     tetrahedralization. It is noted that if weights are present
c     and degeneracies exist there may be points whose Power cells
c     are not 3-dimensional.
c
c     Computations in this program for the purpose of computing
c     the tetrahedralization are done in exact arithmetic whenever
c     floating point arithmetic (done in double precision) does not
c     seem appropriate.
c
c     The program is based on an algorithm by Edelsbrunner and Shah [1]
c     for constructing Regular tetrahedralizations with incremental
c     topological flipping. At the start of the execution of the
c     program a Regular tetrahedralization for the vertices of an
c     artificial cube that contains the input points is constructed.
c     Throughout the execution the vertices of this cube (artificial
c     points) are treated in the proper lexicographical manner [2] so
c     that the final tetrahedralization is correct. The program has the
c     capability of maintaining at all times during the execution of
c     the program a list of all tetrahedra in the current and previous
c     tetrahedralizations. This list is in the form of a directed
c     acyclic graph that represents the history of the flips the
c     program has performed, and it is used by the program to identify
c     a tetrahedron in the current tetrahedralization that contains a
c     new input point. Finally, the program has the capability of
c     adding the input points in a random sequence.
c
c     [1] H. Edelsbrunner and N. R. Shah, Incremental topological
c         flipping works for regular triangulations, Algorithmica 15(3),
c         223-241 (1996).
c
c     [2] J. Bernal, Lexicographical manipulations for correctly
c         computing regular tetrahedralizations with incremental
c         topological flipping, NISTIR 6335 (1999).
c
c     [3] J. Bernal, REGTET: A program for computing regular
c         tetrahedralizations (long version), NISTIR 6786 (2001).
c
c         Author: Javier Bernal
c
c**********************************************************************
c
c     The following are examples of how parameters that are dimensions
c     of arrays can be defined in main routine:
c
c     1. Flipping history used for locating points:
c
c     integer nmax, nvmax, nhmax
c     parameter (nmax=150000, nvmax=55*nmax, nhmax=1500)
c
c     2. Flipping history not used for locating points:
c
c     integer nmax, nvmax, nhmax
c     parameter (nmax=150000, nvmax= 7*nmax, nhmax=1500)
c
c     Arrays and logical variables should be defined in main routine
c     as follows:
c
c     double precision x(nmax), y(nmax), z(nmax), w(nmax)
c     real v(nmax)
c     integer ix(nmax), iy(nmax), iz(nmax), iw(nmax)
c     integer ix2(nmax), iy2(nmax), iz2(nmax), iw2(nmax)
c     integer icon(8,nvmax), is(nmax), ifl(nvmax), io(nmax)
c     integer id(nvmax), ih(nhmax)
c     integer nv, nw, nt, nd, naddl, isu, jsu, ksu, nsu, icfig, iwfig
c     double precision wlenx, wleny, wlenz, wlenw, epz
c     logical delaun, pntoff, flphis, artfcl, random, reccor, redchk
c
c     Subroutine regtet should be called in main routine as follows:
c
c     call regtet(x, y, z, w, v, ix, iy, iz, iw, ix2, iy2, iz2, iw2,
c    *            icon, is, ifl, io, id, ih, nv, nw, nt, nd, nmax,
c    *            nvmax, nhmax, wlenx, wleny, wlenz, wlenw, naddl,
c    *            isu, jsu, ksu, nsu, icfig, iwfig, epz, delaun,
c    *            pntoff, flphis, artfcl, random, reccor, redchk)
c
c**********************************************************************
c
c     Description of variables in calling sequence
c     --------------------------------------------
c
c     For the purpose of describing these variables six sets of points
c     are first described. It should be noted that the descriptions of
c     some of these sets may refer to variables that have not yet been
c     described.
c
c     SI     Set of points on input, i. e. initial set of points
c            provided by the user. Input integer variable nv (described
c            below) equals the number of points in this set.
c
c     SP     Set of points that define a rectangular regular grid on
c            the surface of a rectangular polyhedron that contains SI
c            in its interior and which are internally generated by the
c            program if logical variable reccor (described below) is
c            .true.. This polyhedron is big enough to contain SI and
c            its size and points in SP are determined by the values
c            of variables wlenx, wleny, wlenz and naddl (described
c            below) provided by the user. If logical variable delaun
c            (described below) is .false. and if wmin is the minimum
c            weight of points in SI then a weight equal to wmin - wlenw
c            is assigned by the program to each point in SP where the
c            value of variable wlenw (described below) is provided by
c            the user.
c
c     SA     Set of eight points which are the vertices of an artificial
c            cube. These points are essentially at infinity and must be
c            treated by the program in the proper lexicographical
c            manner. These points are the artificial points.
c
c     SD     Set of points for which a Regular tetrahedralization is
c            desired. Set is defined according to the value of logical
c            variable reccor (described below) as follows:
c            if reccor is .false. then SD equals SI;
c            if reccor is .true. then SD equals the union of SP and SI.
c            If a point in SD has a weight as a member of either SP or
c            SI, it gets the same weight in SD.
c            The points in SD are the real points.
c
c     SU     Set of points for which the program actually computes a
c            Regular tetrahedralization.
c            Set equals the union of SA and SD.
c            If a point in SU is in either SP or SI and has a weight as
c            a member of either set, it gets the same weight in SU.
c            The final tetrahedralization computed by the program is
c            Regular for SU and contains a Regular tetrahedralization
c            for SD.
c
c     SO     Set of points on output. Points in this set are the
c            vertices of tetrahedra in output tetrahedron list. Output
c            integer variable nw (described below) equals the number of
c            points in this set. Set is defined according to the value
c            of logical variable artfcl (described below) as follows:
c            if artfcl is .false. then SO equals SD;
c            if artfcl is .true. then SO equals SU.
c            (As described above SD and SU are defined according to the
c            value of logical variable reccor (described below)).
c            If a point in SO is in either SP or SI and has a weight as
c            a member of either set, it gets the same weight in SO.
c
c            Depending on how SO is defined points in SO are ordered as
c            follows (here ng is the number of points in SP):
c
c            case 1: SO equals SI: for each i, i = 1, ..., nv, the
c            ith point in SI is the ith point in SO;
c            case 2: SO equals the union of SP and SI: for each i,
c            i = 1, ..., ng, the ith point in SP is the ith point in SO,
c            and for each i, i = 1, ..., nv, the ith point in SI is the
c            (i+ng)th point in SO;
c            case 3: SO equals the union of SA and SI: for each i,
c            i = 1, ..., 8, the ith point in SA is the ith point in SO,
c            and for each i, i = 1, ..., nv, the ith point in SI is the
c            (i+8)th point in SO;
c            case 4: SO equals the union of SA, SP, and SI: for each i,
c            i = 1, ..., 8, the ith point in SA is the ith point in SO,
c            for each i, i = 1, ..., ng, the ith point in SP is the
c            (i+8)th point in SO, and for each i, i = 1, ..., nv, the
c            ith point in SI is the (i+ng+8)th point in SO.
c
c            In what follows for each i, i = 1, ..., nw, the ith point
c            in SO is also refered to as point i.
c
c     The description of the variables follows:
c
c     delaun Input logical flag;
c            if .true. (no weights) then a Delaunay tetrahedralization
c            is to be computed;
c            if .false. then a Regular tetrahedralization is to be
c            computed.
c
c     pntoff Input logical flag;
c            if .true. then some input points will be inactive during
c            tetrahedralization computation;
c            if .false. then all input points will be active during
c            tetrahedralization computation.
c
c     flphis Input logical flag;
c            if .true. then tetrahedron list with flipping history will
c            be used for locating points;
c            if .false. then tetrahedron list with flipping history will
c            not be used for locating points; a shishkebab method that
c            locates points by checking tetrahedra in the current
c            tetrahedralization will be used.
c
c            The use of the tetrahedron list with flipping history
c            requires a lot more space than the alternative; parameter
c            nvmax (described below) must be set to about 55 times
c            the value of parameter nmax (described below) if flphis
c            equals .true.; otherwise it must be set to about 7 times
c            the value of parameter nmax.
c
c     artfcl Input logical flag;
c            if .true. and if logical variable flphis (described above)
c            is set to .true. then output tetrahedron list will include
c            the final tetrahedralization for SO (equal to SU) together
c            with the flipping history tetrahedra, i.e. tetrahedra that 
c            at some point during the execution of the program were part
c            of a tetrahedralization for SO but that are not in the
c            final tetrahedralization;
c            if .true. and if flphis equals .false. then output
c            tetrahedron list will just be the final tetrahedralization
c            for SO (equal to SU).
c            if .false. then output tetrahedron list will only include
c            the final tetrahedralization for SO (equal to SD)).
c
c     random Input logical flag;
c            if .true. then the points in SI are to be added by the
c            program in a random fashion;
c            if .false. then the points in SI are to be added by the
c            program in their original order.
c
c            If the points in SI are already randomized on input then
c            there is no need for the program to randomize them again
c            so that logical variable random should be set to .false..
c            If the points are in a nice order and the shishkebab
c            method is to be used for locating points then random
c            should be set to .false..
c
c     reccor Input logical flag;
c            if .true. then SD includes SP;
c            if .false. then SD does not include SP.
c
c            Including SP in SD is recommended if it is not desirable
c            to have points from SI in the boundary of the convex hull
c            of SD. However the final tetrahedralization for SO, even
c            though it is Regular and equals or contains a Regular
c            tetrahedralization for SD, will not necessarily contain a
c            Regular tetrahedralization for SI.
c
c     redchk Input logical flag;
c            used only if logical variable delaun (described above) is
c            set to .false. (weights are being used);
c            if .true. then redundant points are to be tested for
c            redundancy after final tetrahedralization has been
c            computed;
c            if .false. then they are not to be tested for redundancy.
c
c     nmax   Input integer variable;
c            must be defined in a parameter statement in main routine;
c            dimension of single arrays x, y, z, w, v, ix, iy, iz, iw,
c            ix2, iy2, iz2, iw2, is, io (all are described below).
c
c     nvmax  Input integer variable;
c            must be defined in a parameter statement in main routine;
c            second dimension of double array icon and dimension of
c            single arrays ifl, id (these arrays are described below);
c            nvmax should be set to 55*nmax if logical variable flphis
c            (described above) is set to .true.; to 7*nmax otherwise.
c
c     nhmax  Input integer variable;
c            must be defined in a parameter statement in main routine;
c            dimension of single array ih (described below).
c            
c     nv     Input integer variable that can not exceed nmax;
c            number of points or vertices in SI;
c            same value on output.
c
c     nw     Output integer variable that can not exceed nmax;
c            number of points or vertices in SO.
c
c     nt     Output integer variable that can not exceed nvmax;
c            number of tetrahedra in final tetrahedron list;
c            if logical variable artfcl is .false. then nt equals
c            the number of tetrahedra in the final tetrahedralization
c            for SO (equal to SD);
c            if artfcl is .true. and logical variable flphis is .true.
c            then nt equals the number of tetrahedra in the  final
c            tetrahedralization for SO (equal to SU) plus the
c            number of flipping history tetrahedra, i. e. tetrahedra
c            that at some point during the execution of the program were
c            part of a tetrahedralization for SO but that are not part
c            of the final tetrahedralization (a tetrahedron in a
c            previous tetrahedralization will not be in the final
c            tetrahedralization if at some time during the execution of
c            the program it has been eliminated and replaced by other
c            tetrahedra through a flip);
c            if artfcl is .true. and flphis is .false. then nt equals
c            the number of tetrahedra in the final tetrahedralization
c            for SO (equal to SU).
c
c     nd     Output integer variable that can not exceed nvmax;
c            number of real tetrahedra in final tetrahedron list;
c            i. e. number of tetrahedra in Regular tetrahedralization
c            for SD (the real points) which is contained in final
c            Regular tetrahedralization computed by the program for SU.
c
c     x      Input/output real*8 single array of dimension nmax;
c            on input for each i, i = 1, ..., nv, x(i) is the
c            x-coordinate of the ith point in SI;
c            with icfig as decribed below, on output for each i,
c            i = 1, ..., nw, if the ith point in SO is not in SA
c            then x(i) is the x-coordinate of the ith point in SO
c            rounded off so that its decimal part has icfig
c            significant figures; if it is in SA then x(i) is a
c            program generated value associated internally by the
c            program with the ith point in SO.
c
c     y      Input/output real*8 single array of dimension nmax;
c            on input for each i, i = 1, ..., nv, y(i) is the
c            y-coordinate of the ith point in SI;
c            with icfig as decribed below, on output for each i,
c            i = 1, ..., nw, if the ith point in SO is not in SA
c            then y(i) is the y-coordinate of the ith point in SO
c            rounded off so that its decimal part has icfig
c            significant figures; if it is in SA then y(i) is a
c            program generated value associated internally by the
c            program with the ith point in SO.
c
c     z      Input/output real*8 single array of dimension nmax;
c            on input for each i, i = 1, ..., nv, z(i) is the
c            z-coordinate of the ith point in SI;
c            with icfig as decribed below, on output for each i,
c            i = 1, ..., nw, if the ith point in SO is not in SA
c            then z(i) is the z-coordinate of the ith point in SO
c            rounded off so that its decimal part has icfig
c            significant figures; if it is in SA then z(i) is a
c            program generated value associated internally by the
c            program with the ith point in SO.
c
c     w      Input/output real*8 single array of dimension nmax;
c            on input for each i, i = 1, ..., nv, w(i) is the
c            weight of the ith point in SI;
c            with iwfig as decribed below, on output for each i,
c            i = 1, ..., nw, if the ith point in SO is not in SA
c            then w(i) is the weight of the ith point in SO
c            rounded off so that its decimal part has iwfig
c            significant figures; if it is in SA then w(i) is a
c            program generated value associated internally by the
c            program with the ith point in SO.
c
c     v      Real single array of dimension nmax;
c            internally used by program.
c
c     ix     Integer single array of dimension nmax;
c            internally used by the program.
c
c     iy     Integer single array of dimension nmax;
c            internally used by the program.
c
c     iz     Integer single array of dimension nmax;
c            internally used by the program.
c
c     iw     Integer single array of dimension nmax;
c            internally used by the program.
c
c     ix2    Integer single array of dimension nmax;
c            internally used by the program.
c
c     iy2    Integer single array of dimension nmax;
c            internally used by the program.
c
c     iz2    Integer single array of dimension nmax;
c            internally used by the program.
c
c     iw2    Integer single array of dimension nmax;
c            internally used by the program.
c
c     icon   Output integer double array of dimensions 8 and nvmax;
c            this is the tetrahedron list;
c            actually this is a list of 8 x nt integers, and it is a
c            list of tetrahedra in the sense that for each j,
c            j = 1, ..., nt, the 8 integers icon(i,j), i = 1, .., 8, are
c            associated with a tetrahedron, the jth tetrahedron or
c            tetrahedron j, as will be described below;
c            if logical variables artfcl and flphis are both .true.
c            then tetrahedra in this list are those in the final
c            tetrahedralization for SO (equal to SU) together with the
c            flipping history tetrahedra, i. e. tetrahedra that at some
c            point during the execution of the program were part of a
c            tetrahedralization for SO but that are not part  of the
c            final tetrahedralization;
c            if artfcl is .true. and flphis is .false. then tetrahedra
c            in this list are just those in the final tetrahedralization
c            for SO (equal to SU);
c            if artfcl is .false. then tetrahedra in this list are those
c            in the final tetrahedralization for SO (equal to SD);
c            for each j, j = 1, ..., nt, if icon(5,j) is negative (can
c            only happen if artfcl equals .true.) then tetrahedron j
c            is not in the final tetrahedralization for SO, it was in a
c            previous tetrahedralization (it was eliminated) and its
c            vertices are the points -icon(5,j), icon(6,j), icon(7,j),
c            and icon(8,j) in SO; in addition if flphis is .true. the
c            tetrahedra by which tetrahedron j was replaced through a
c            flip can be identified  as follows: for each i,
c            i = 1, ..., 4, if icon(i,j) is positive then tetrahedron
c            icon(i,j) is one of those tetrahedra;
c            for each j, j = 1, ..., nt, if icon(5,j) is positive then
c            tetrahedron j is in the final tetrahedralization for SO,
c            and its vertices are the points icon(5,j), icon(6,j),
c            icon(7,j), and icon(8,j) in SO; in addition the tetrahedra
c            in the final tetrahedralization that share a facet with
c            tetrahedron j can be identified as follows: for each i,
c            i = 1, ..., 4, if icon(i,j) is positive then tetrahedron
c            icon(i,j) is one of those tetrahedra;
c            for each j, j = 1, ..., nt, the vertices of tetrahedron j
c            are ordered as follows: when viewed from vertex icon(5,j)
c            (-icon(5,j) if icon(5,j) is negative) the other three
c            vertices icon(6,j), icon(7,j), icon(8,j) appear in this
c            order in a clockwise direction around the circle that
c            contains them;
c            for each j, j = 1, ..., nt, if tetrahedron j is in the
c            final tetrahedralization, i. e. icon(5,j) is positive,
c            then the tetrahedra in the final tetrahedralization that
c            share a facet with tetrahedron j are ordered as follows:
c            for each i, i = 1, ..., 4, if icon(i,j) is positive,
c            tetrahedron j shares with tetrahedron icon(i,j) the
c            facet of tetrahedron j that does not contain vertex
c            icon(i+4,j).
c
c     is     Input/output integer single array of dimension nmax;
c            on input if logical variable pntoff is .true. then for
c            each i, i = 1, ..., nv, if the value of is(i) equals 1
c            then the ith point in SI is to be active during the
c            tetrahedralization computation; if it equals 0 then the
c            ith point in SI is to be inactive during the computation;
c            on input if logical variable pntoff is .false. then
c            all points in SI are to be active during the
c            tetrahedralization computation and array is does not
c            have to be set to any values;
c            on output for each i, i = 1, ..., nw, the value of is(i)
c            indicates how the ith point in SO was processed by the
c            program as follows:
c            if is(i) is zero then point i was not considered as a
c            vertex for tetrahedralization;
c            if is(i) is positive then point i is part of the final
c            tetrahedralization for SO, i. e. there is at least one
c            tetrahedron in the final tetrahedralization with point i
c            as a vertex, and tetrahedron is(i) is one such tetrahedron
c            (actually if point i is in SA then is(i) is always
c            positive);
c            if is(i) is less than -8 then point i was found to be
c            redundant as the program was trying to insert it into the
c            current tetrahedralization because a point previously
c            inserted (point -is(i) in SO if artfcl is .true. (SO equals
c            SU), point -is(i)-8 in SO if artfcl is .false. (SO equals
c            SD)) was identical to it and either the weight of the
c            previously inserted point was larger or equal to the weight
c            of point i or there were no weights;
c            if is(i) equals -2 then point i had been inserted by the
c            program into the tetrahedralization but was found to be
c            redundant because another point was later inserted by the
c            program that was identical to point i and whose weight was
c            larger than that of point i (this case is not possible if
c            there are no weights);
c            if is(i) equals -3 then point i was found to be redundant
c            in the sense of a Regular tetrahedralization as the program
c            was trying to insert it into the current tetrahedralization
c            because of its weight as compared to the weights of the
c            vertices of the tetrahedron in the current
c            tetrahedralization that contains it even though it was not
c            identical to a previously inserted point (this case is not
c            possible if there are no weights);
c            if is(i) equals -4 then point i had been inserted by the
c            program into the tetrahedralization but was found to be
c            redundant in the sense of a Regular tetrahedralization
c            because of the weight of another point, not identical to
c            point i, that was later inserted by the program together
c            with the weights of three other previously inserted points
c            as compared to the weight of point i (this case is not
c            possible if there are no weights).
c
c     ifl    Integer single array of dimension nvmax;
c            internally used by program.
c
c     io     Integer single array of dimension nmax;
c            internally used by program.
c
c     id     Integer single array of dimension nvmax;
c            internally used by program.
c
c     ih     Integer single array of dimension nhmax;
c            internally used by program.
c
c     wlenx
c     wleny
c     wlenz  Input real*8 variables;
c            If reccor is .true. then these are three positive real
c            numbers provided by the user to be used by the program
c            to identify a rectangular polyhedron that contains SI in
c            its interior. This is the polyhedron whose surface will
c            contain the set SP. If xmax, ymax, zmax are the maximum
c            values of the x, y, and z coordinates of the points in SI,
c            respectively, and xmin, ymin, zmin are the minimum values,
c            then the eight vertices of the polyhedron will be:
c            (xmin-wlenx, ymin-wleny, zmin-wlenz),
c            (xmax+wlenx, ymin-wleny, zmin-wlenz),
c            (xmax+wlenx, ymax+wleny  zmin-wlenz),
c            (xmin-wlenx, ymax+wleny, zmin-wlenz),
c            (xmin-wlenx, ymin-wleny, zmax+wlenz),
c            (xmax+wlenx, ymin-wleny, zmax+wlenz),
c            (xmax+wlenx, ymax+wleny  zmax+wlenz),
c            (xmin-wlenx, ymax+wleny, zmax+wlenz).
c
c     wlenw  Input real*8 variable;
c            If reccor is .true. and delaun is .false. then this is a
c            real number provided  by the user to be used by the program
c            to determine a weight to be assigned to each point in SP.
c            If wmin is the minimum value of the weights of the points
c            in SI then this weight will be wmin - wlenw.
c
c     naddl  Input integer variable;
c            If reccor is .true. then this is a positive integer
c            greater than 1 provided by the user to be used by the
c            program to determine the set SP. The points in SP define
c            a rectangular regular grid on the surface of a rectangular
c            polyhedron (described above) that contains SI in its
c            interior. For each facet of the polyhedron a set of naddl
c            x naddl points is generated by the program that defines a
c            rectangular regular grid on the facet and that contains
c            the four vertices of the facet. SP is then the union of
c            the six sets thus generated (one per facet). It then
c            follows that the number of points in SP must be
c            6(naddl-2)(naddl-2)+12(naddl-2)+8 which reduces to 
c            6(naddl)(naddl)-12(naddl)+8.
c            It also follows that if naddl equals 2 then the points in
c            SP are exactly the 8 vertices of the polyhedron.
c
c     isu
c     jsu
c     ksu
c     nsu    Input integer variables;
c            If random is .true. then these are four integers provided
c            by the user to be used by program as seeds for identifying
c            random order in which points in SI are to be added;
c            they can be any four integers.
c
c     icfig  Input integer variable;
c            value is the number of significant figures of decimal part
c            of coordinates of input points; value should be nonnegative
c            and less than 10.
c
c     iwfig  Input integer variable;
c            value is the number of significant figures of decimal part
c            of weights (if any) of input points; value should be
c            nonnegative, less than 10, and not greater than twice the
c            value of variable icfig (described above).
c
c     epz    Input real*8 variable;
c            tolerance used by the program to switch from floating
c            point arithmetic to exact arithmetic by testing against
c            this tolerance whether certain quantities are too close
c            to zero; setting it equal to numbers such as 0.1, 0.01
c            has worked well so far.
c
c**********************************************************************
c
c     Examples of settings for logical variables delaun, pntoff,
c     flphis, artfcl, random, reccor, and redchk:
c
c     Delaunay tetrahedralization for set of randomized input points
c     is desired and nothing else (number of input points equals number
c     of output points, all input points are to be active during the
c     tetrahedralization computation, and tetrahedron list is exactly
c     a list of the tetrahedra in the final tetrahedralization for the
c     set of input points):
c
c     delaun = .true.
c     pntoff = .false.
c     flphis = .true.
c     artfcl = .false.
c     random = .false.
c     reccor = .false.
c
c     The same as above but a Regular tetrahedralization is desired,
c     the input points are not randomized, the flipping history
c     is to be used for locating points, and redundant points are not
c     to be tested for redundancy after the final tetrahedralization
c     has been computed:
c
c     delaun = .false.
c     pntoff = .false.
c     flphis = .true.
c     artfcl = .false.
c     random = .true.
c     reccor = .false.
c     redchk = .false.
c
c     The same as above but a Regular tetrahedralization is desired,
c     the input points are not randomized, the output tetrahedron
c     list is to include artificial tetrahedra information, and the
c     shishkebab method is to be used for locating points:
c
c     delaun = .false.
c     pntoff = .false.
c     flphis = .false.
c     artfcl = .true.
c     random = .true.
c     reccor = .false.
c     redchk = .false.
c
c**********************************************************************
c 
      subroutine regtet(x, y, z, w, v, ix, iy, iz, iw, ix2, iy2, iz2,
     *                  iw2, icon, is, ifl, io, id, ih, nv, nw, nt, nd,
     *                  nmax, nvmax, nhmax, wlenx, wleny, wlenz, wlenw,
     *                  naddl, isu, jsu, ksu, nsu, icfig, iwfig, epz,
     *                  delaun, pntoff, flphis, artfcl, random, reccor,
     *                  redchk)
c
      integer nmax, nvmax, nhmax
      double precision x(*), y(*), z(*), w(*)
      real v(*)
      integer ix(*), iy(*), iz(*), iw(*)
      integer ix2(*), iy2(*), iz2(*), iw2(*)
      integer icon(8,*), is(*), ifl(*), io(*)
      integer id(*), ih(*)
      integer nv, nw, nt, nd, naddl, isu, jsu, ksu, nsu, icfig, iwfig
      double precision wlenx, wleny, wlenz, wlenw, epz
      logical delaun, pntoff, flphis, artfcl, random, reccor, redchk
c
      integer isclp(2), isclw(2), isclr(2), tetra, tetru
      double precision xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax
      double precision xpre, ypre, zpre, wpre, xnow, ynow, znow, wnow
      integer mhalf, mfull, i, ihn, iftal, no, ipre, jpre, inow, jnow
      integer irec, irec1, nv1, nu, mxlook, j, npre, ni
      integer iredx, iconx, iorix, isphx
      integer iredp, ired0, ired1, ired2, ired3, ired4
c     integer nva
c
c     initialize Fortran 77 word lengths
c
      mhalf=32768
      mfull=1073741824
c
c     testing parameters and number of input points or vertices
c
      if (nv .lt. 1 .or. nv .gt. nmax .or. nvmax .lt. 12) then 
          write(*,*) nv, nmax, nvmax
          stop 110
      endif
c
c     initialize arrays ih, w, is and id 
c
      do 50 i = 1, nhmax
         ih(i) = 0
   50 continue
      if(delaun)then
         do 60 i = 1, nmax
            w(i) = 0.0d0
   60    continue
      endif
      if(.not.pntoff)then
         do 80 i = 1, nmax
            is(i) = 1
   80    continue
      endif
      if(.not.flphis)then
         ihn = 0
         iftal = 0
         do 100 i = 1, nvmax
            id(i) = 0
  100    continue
      endif
c
c     test variables associated with a possible rectangular polyhedron
c
      if(reccor)then
         if(wlenx.le.0.0d0 .or. wleny.le.0.0d0 .or. wlenz.le.0.0d0)
     *      stop 120
         if(naddl.lt.2) stop 130
      else
         wlenx = 0.0d0
         wleny = 0.0d0
         wlenz = 0.0d0
         wlenw = 0.0d0
         naddl = 0
      endif
c
c     calculating min and max
c
      xmax = x(1)
      xmin = x(1)
      ymax = y(1)
      ymin = y(1)
      zmax = z(1)
      zmin = z(1)
      wmax = w(1)
      wmin = w(1)
      do 150 no = 1, nv
         if (x(no) .gt. xmax) xmax = x(no)
         if (x(no) .lt. xmin) xmin = x(no)
         if (y(no) .gt. ymax) ymax = y(no)
         if (y(no) .lt. ymin) ymin = y(no)
         if (z(no) .gt. zmax) zmax = z(no)
         if (z(no) .lt. zmin) zmin = z(no)
         if (w(no) .gt .wmax) wmax = w(no)
         if (w(no) .lt. wmin) wmin = w(no)
  150 continue
c
c     if randomizing input data then obtain random order of integers
c     from 1 to nv and randomize data
c
      if(.not.random) go to 185
      call rdmord(v, io, nv, isu, jsu, ksu, nsu)
      do 180 no = 1, nv
         if(io(no).lt.0) go to 180
         ipre = no
         xpre = x(ipre)
         ypre = y(ipre)
         zpre = z(ipre)
         wpre = w(ipre)
         jpre = is(ipre)
  170    continue
         inow = io(ipre)
         io(ipre) = -inow
         xnow = x(inow)
         ynow = y(inow)
         znow = z(inow)
         wnow = w(inow)
         jnow = is(inow)
         x(inow) = xpre
         y(inow) = ypre
         z(inow) = zpre
         w(inow) = wpre
         is(inow) = jpre
         if(inow .eq. no) go to 180
         ipre = inow
         xpre = xnow
         ypre = ynow
         zpre = znow
         wpre = wnow
         jpre = jnow
         if(io(ipre) .lt. 0) stop 140
         go to 170
  180 continue
c
c     OPEN(22,FILE='ran.dat')
c     DO 182 I=1,NV
c        WRITE(22,*)X(I),Y(I),Z(I),W(I)
c 182 CONTINUE
c
c     shift data
c
  185 continue
      irec = 8
      if(reccor) irec = irec + 6*(naddl**2) - 12*naddl + 8
      irec1 = irec + 1
      nv1 = nv
      nv = nv + irec
      if(nv .gt. nmax) stop 150
      do 190 no = nv, irec1, -1
         nu = no - irec
         x(no) = x(nu)
         y(no) = y(nu)
         z(no) = z(nu)
         w(no) = w(nu)
         is(no) = is(nu)
         if(random)then
            if(io(nu).ge.0) stop 160
            io(nu) = -io(nu) + irec
         endif
  190 continue
c
c     initialize is for additional data
c
      do 200 i = 1, irec
         is(i) = 1
  200 continue
c
c     write(*,*)' '
c     write(*,*)'Entering poltri ...'
      call poltri(x, y, z, w, ix, iy, iz, iw, ix2, iy2, iz2, iw2, icon,
     *            is, ifl, id, ih, ihn, xmin, ymin, zmin, wmin, xmax,
     *            ymax, zmax, wmax, iftal, nv, nmax, nvmax, nhmax,
     *            wlenx, wleny, wlenz, wlenw, tetra, mxlook, irec,
     *            naddl, iredx, delaun, flphis, redchk, icfig, iwfig,
     *            mhalf, mfull, isclp, isclw, isclr, epz)
!     write(*,*)' '
!     write(*,*)'Checking tetrahedralization ...'
c     WRITE(*,*)' '
c     WRITE(*,*)'MAXLOOK=',MAXLOOK,' IHN=',IHN
c     write (*,*)' '
c     write (*,*)'Leaving poltri ...'
c
c     write(*,*)' '
c     write(*,*)'Entering consis ...'
      call consis(icon, is, ifl, nv, tetra)
c
c     write(*,*)' '
c     write(*,*)'Entering convex ...'
      call convex(icon, tetra, ifl, x, y, z, ix, iy, iz, ix2, iy2, iz2,
     *            iconx, mhalf, mfull, isclp, epz)
c
c     write(*,*)' '
c     write(*,*)'Entering orient ...'
      call orient(tetra, icon, ifl, x, y, z, ix, iy, iz, ix2, iy2, iz2,
     *            iorix, mhalf, mfull, isclp, epz)
c
c     write(*,*)' '
c     write(*,*)'Entering delchk ...'
      call delchk(tetra, icon, ifl, x, y, z, w, ix, iy, iz, iw, ix2,
     *            iy2, iz2, iw2, isphx, delaun, mhalf, mfull, isclp,
     *            isclw, isclr, epz)
c
c     checking for possible warnings
c
      if(iredx.ne.0) then
         write(*,*)' '
         write(*,*)'Warning: redundancy violations detected'
         write(*,*)'Number of violations = ',iredx
      endif
      if(iconx.ne.0) then
         write(*,*)' '
         write(*,*)'Warning: convexity violations detected'
         write(*,*)'Number of violations = ',iconx
      endif
      if(iorix.ne.0) then
         write(*,*)' '
         write(*,*)'Warning: tetrahedra orientation violations detected'
         write(*,*)'Number of violations = ',iorix
      endif
      if(isphx.ne.0) then
         write(*,*)' '
         write(*,*)'Warning: sphere criterion violations detected'
         write(*,*)'Number of violations = ',isphx
      endif
      if(iredx.ne.0 .or. iconx.ne.0 .or. iorix.ne.0 .or. isphx.ne.0)then
         write(*,*)' '
         write(*,*)'Increasing tolerance epz could improve situation'
      endif
c
c     readjust data structure for randomizing
c
      if(.not.random) go to 290
c     write(*,*)' '
c     write(*,*)'Readjusting data structure for randomizing ...'
      if(nv .gt. nvmax) stop 170
      nu = nv-irec
      do 230 no = 1, nu
         ifl(io(no)) = no + irec
  230 continue
      do 250 i = 1, tetra
         do 240 j = 5, 8
            if(icon(j,i).gt.irec) then
               icon(j,i) = ifl(icon(j,i))
            elseif(icon(j,i).lt.-irec) then
               icon(j,i) = -ifl(-icon(j,i))
            endif
  240    continue
  250 continue
      do 255 i = irec1, nv
         if(is(i).lt.-8) is(i) = -ifl(-is(i))
  255 continue
      do 260 i = irec1, nv
         ifl(i) = is(i)
  260 continue
      do 270 i = irec1, nv
         is(i) = ifl(io(i-irec))
  270 continue
      do 280 no = 1, nv1
         if(io(no).lt.0) go to 280
         nu = no + irec
         ipre = nu
         npre = no
         xpre = x(ipre)
         ypre = y(ipre)
         zpre = z(ipre)
         wpre = w(ipre)
  275    continue
         inow = io(npre)
         io(npre) = -inow
         if(inow .eq. nu) then
            x(ipre) = xpre
            y(ipre) = ypre
            z(ipre) = zpre
            w(ipre) = wpre
            go to 280
         endif
         x(ipre) = x(inow)
         y(ipre) = y(inow)
         z(ipre) = z(inow)
         w(ipre) = w(inow)
         ipre = inow
         npre = ipre - irec
         if(io(npre) .lt. 0) stop 180
         go to 275
  280 continue
c
c     write(*,*)' '
c     write(*,*)'Entering consis ...'
      call consis(icon, is, ifl, nv, tetra)
c
  290 continue
      nu=nv-8
      if(.not.artfcl) then
c        write(*,*)' '
c        write(*,*)'Entering revtet ...'
         call revtet(tetra, tetru, icon, nv, is, ifl, flphis)
         do 293 no = 1, nu
            x(no) = x(no+8)
            y(no) = y(no+8)
            z(no) = z(no+8)
            w(no) = w(no+8)
  293    continue
         if(tetru .eq. 0) go to 300
c
c        write(*,*)' '
c        write(*,*)'Entering consis ...'
         call consis(icon, is, ifl, nu, tetru)
      elseif(.not.flphis) then
         call ruvtet(tetra, tetru, icon, is, ifl)
         call consis(icon, is, ifl, nv, tetra)
      else
c 
c        count true tetrahedra
c
         tetru = 0
         do 295 i = 1, tetra
            if ((icon(6,i) .le. 8) .or. (icon(7,i) .le. 8) .or.
     *      (icon(8,i) .le. 8) .or. (icon(5,i) .le. 8)) goto 295
            tetru = tetru + 1
  295    continue
      endif
c
c     count redundant vertices
c
  300 continue
      nd = tetru
      iredp = 0
      ired0 = 0
      ired1 = 0
      ired2 = 0
      ired3 = 0
      ired4 = 0
      if(artfcl) then
         ni = 9
         nw = nv
         nt = tetra
      else
         ni = 1
         nw = nu
         nt = tetru
      endif
      do 400 i = ni, nw
         if(is(i) .gt. 0) then
            iredp = iredp + 1
         elseif(is(i) .eq.  0) then
            ired0 = ired0 + 1
         elseif(is(i) .lt. -8) then
            ired1 = ired1 + 1
         elseif(is(i) .eq. -2) then
            ired2 = ired2 + 1
         elseif(is(i) .eq. -3) then
            ired3 = ired3 + 1
         elseif(is(i) .eq. -4) then
            ired4 = ired4 + 1
         else
            stop 190
         endif
  400 continue
c
c     OPEN(23,FILE='unr.dat')
c     DO 500 I=1,NW
c        WRITE(23,*)X(I),Y(I),Z(I),W(I)
c 500 CONTINUE
c
c     nva = nv
      nv = nv1
c 
c     write info to screen
c
c     wtenv=float(tetra)/float(nva)
c     wtena=float(tetra)/float(iredp)
c     wtuna=float(tetru)/float(iredp)
c     write (*,*) ' '
c     write (*,*) 'Tetrahedralization data ...'
c     write (*,*) ' '
c     write (*,*) 'minimum weight = ',wmin
c     write (*,*) 'maximum weight = ',wmax
c     write (*,*) 'number of true vertices: ', nu
c     write (*,*) 'number of active vertices: ',iredp
c     write (*,*) 'maximum number of vertices parameter = ', nmax
c     write (*,*) 'maximum number of tetrahed parameter = ', nvmax
c     write (*,*) 'number of tetrahedra of all kinds: ', tetra
c     write (*,*) 'all tetrahedra-all vertices ratio: ',wtenv
c     write (*,*) 'number of true tetrahedra: ', tetru
c     write (*,*) ' all tetrahedra-active vertices ratio: ',wtena
c     write (*,*) 'true tetrahedra-active vertices ratio: ',wtuna
c     write (*,*) 'max levels gone down in hierarchy = ', mxlook
c     write (*,*) 'points active at the end of current run   = ',iredp
c     write (*,*) 'points inactive at the end of current run = ',ired0
c     write (*,*) 'points redundant by initial substitution  = ',ired1
c     write (*,*) 'points redundant by eventual substitution = ',ired2
c     write (*,*) 'points redundant by initial elimination   = ',ired3
c     write (*,*) 'points redundant by eventual elimination  = ',ired4
c
      return
      end
*POLTRI
c
c     This subroutine will obtain initial cube and will divide it
c     into 12 tetrahedra; insert points into tetrahedralization
c 
      subroutine poltri(x, y, z, w, ix, iy, iz, iw, ix2, iy2, iz2, iw2,
     *                  icon, is, ifl, id, ih, ihn, xmin, ymin, zmin,
     *                  wmin, xmax, ymax, zmax, wmax, iftal, nv, nmax,
     *                  nvmax, nhmax, wlenx, wleny, wlenz, wlenw, tetra,
     *                  mxlook, irec, naddl, iredx, delaun, flphis,
     *                  redchk, icsfig, iwsfig, mhalf, mfull, isclp,
     *                  isclw, isclr, epz)
c
      double precision x(*), y(*), z(*), w(*)
      integer ix(*), iy(*), iz(*), iw(*)
      integer ix2(*), iy2(*), iz2(*), iw2(*)
      integer icon(8,*), is(*), ifl(*), id(*), ih(*)
      double precision xc(8), yc(8), zc(8)
      integer ixc(8), iyc(8), izc(8)
      integer nmax, nvmax, nhmax
      integer ihn, nv, naddl, icsfig, iwsfig, mxlook, irec, iredx
      integer mhalf, mfull, ibsfig, itsfig
      double precision wlenx, wleny, wlenz, wlenw, epz, wbig, wmen, wman
      logical delaun, flphis, redchk
      integer isclp(*), isclw(*), isclr(*), tetra
      double precision xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax
      double precision xlen, ylen, zlen, wlen, xctr, yctr, zctr
      double precision xlan, ylan, zlan, wlan, xint, yint, zint
      double precision xcor, ycor, zcor
      double precision xmx, ymx, zmx, xmn, ymn, zmn
      double precision dscle, dscli, dfull, dfill, decml
      integer iftal, i, naddm, ng, j, i9, icsfi2, irsfig, isclu, isgcl
      integer i10, issin, k, iredu
      integer itchk, itcn1, itcn2, itcn3, itcn4
c
c     initialize
c
      mxlook = 0
c
c     test parameters
c
      if (nv .gt. nmax) stop 210
      if (nv .lt. 9) stop 220
      if (nvmax .lt. 12) stop 230
c
c     construct cube
c
      xlen=xmax-xmin
      ylen=ymax-ymin
      zlen=zmax-zmin
      wlan=(dmax1(xlen,ylen,zlen)/2.0d0) + dmax1(wlenx,wleny,wlenz)
      wlen=wlan + 4.0d0
      if(wlen.le.wlan) stop 235
c
      xctr=(xmax+xmin)/2.0d0
      yctr=(ymax+ymin)/2.0d0
      zctr=(zmax+zmin)/2.0d0
c     WRITE(*,*)'XYCTR=',XCTR,YCTR,ZCTR,' WLEN=',WLEN
c
c     identify cube corner direction vectors
c
      xc(1) = - 1.0d0
      yc(1) = - 1.0d0
      zc(1) = + 1.0d0
c
      xc(2) = - 1.0d0
      yc(2) = + 1.0d0
      zc(2) = + 1.0d0
c
      xc(3) = + 1.0d0
      yc(3) = + 1.0d0
      zc(3) = + 1.0d0
c
      xc(4) = + 1.0d0
      yc(4) = - 1.0d0
      zc(4) = + 1.0d0
c
      xc(5) = - 1.0d0
      yc(5) = - 1.0d0
      zc(5) = - 1.0d0
c
      xc(6) = - 1.0d0
      yc(6) = + 1.0d0
      zc(6) = - 1.0d0
c
      xc(7) = + 1.0d0
      yc(7) = + 1.0d0
      zc(7) = - 1.0d0
c
      xc(8) = + 1.0d0
      yc(8) = - 1.0d0
      zc(8) = - 1.0d0
c
c     compute corners of cube
c
      do 50 i=1,8
         x(i)=xctr+wlen*xc(i)
         y(i)=yctr+wlen*yc(i)
         z(i)=zctr+wlen*zc(i)
         if((x(i).ge.xmin.and.x(i).le.xmax).or.(y(i).ge.ymin.and.
     *       y(i).le.ymax).or.(z(i).ge.zmin.and.z(i).le.zmax)) stop 240
   50 continue
c
c     make coordinates of corners of cube into whole numbers
c
      dfull = dble(mfull)
      if(dabs(x(3)+1.0d0).ge.dfull) stop 242
      if(dabs(y(3)+1.0d0).ge.dfull) stop 243
      if(dabs(z(3)+1.0d0).ge.dfull) stop 244
      if(dabs(x(5)-1.0d0).ge.dfull) stop 245
      if(dabs(y(5)-1.0d0).ge.dfull) stop 246
      if(dabs(z(5)-1.0d0).ge.dfull) stop 247
c
      xmx = dble(idnint(x(3)+1.0d0))
      ymx = dble(idnint(y(3)+1.0d0))
      zmx = dble(idnint(z(3)+1.0d0))
      xmn = dble(idnint(x(5)-1.0d0))
      ymn = dble(idnint(y(5)-1.0d0))
      zmn = dble(idnint(z(5)-1.0d0))
c
      xlan = xmx - xmn
      ylan = ymx - ymn
      zlan = zmx - zmn
      wlan = dmax1(xlan,ylan,zlan) 
c
      x(1) = xmn
      y(1) = ymn
      z(1) = zmn + wlan
c
      x(2) = xmn
      y(2) = ymn + wlan
      z(2) = zmn + wlan
c
      x(3) = xmn + wlan
      y(3) = ymn + wlan
      z(3) = zmn + wlan
c
      x(4) = xmn + wlan
      y(4) = ymn
      z(4) = zmn + wlan
c
      x(5) = xmn
      y(5) = ymn
      z(5) = zmn
c
      x(6) = xmn
      y(6) = ymn + wlan
      z(6) = zmn
c
      x(7) = xmn + wlan
      y(7) = ymn + wlan
      z(7) = zmn
c
      x(8) = xmn + wlan
      y(8) = ymn
      z(8) = zmn
c
      do 55 i=1,8
         if((x(i).ge.xmin.and.x(i).le.xmax).or.(y(i).ge.ymin.and.
     *       y(i).le.ymax).or.(z(i).ge.zmin.and.z(i).le.zmax)) stop 250
   55 continue
c
      if(irec.eq.8) go to 155
c
c     compute corners of rectangular polyhedron
c
      x( 9) = xmin - wlenx
      y( 9) = ymin - wleny
      z( 9) = zmax + wlenz
c
      x(10) = xmin - wlenx
      y(10) = ymax + wleny
      z(10) = zmax + wlenz
c
      x(11) = xmax + wlenx
      y(11) = ymax + wleny
      z(11) = zmax + wlenz
c
      x(12) = xmax + wlenx
      y(12) = ymin - wleny
      z(12) = zmax + wlenz
c
      x(13) = xmin - wlenx
      y(13) = ymin - wleny
      z(13) = zmin - wlenz
c
      x(14) = xmin - wlenx
      y(14) = ymax + wleny
      z(14) = zmin - wlenz
c
      x(15) = xmax + wlenx
      y(15) = ymax + wleny
      z(15) = zmin - wlenz
c
      x(16) = xmax + wlenx
      y(16) = ymin - wleny
      z(16) = zmin - wlenz
c
      do 60 i=9,16
         if((x(i).ge.xmin.and.x(i).le.xmax).or.(y(i).ge.ymin.and.
     *       y(i).le.ymax).or.(z(i).ge.zmin.and.z(i).le.zmax)) stop 260
   60 continue
      if(x(1).ge.x( 9) .or. y(1).ge.y( 9) .or. z(1).le.z( 9) .or.
     *   x(2).ge.x(10) .or. y(2).le.y(10) .or. z(2).le.z(10) .or.
     *   x(3).le.x(11) .or. y(3).le.y(11) .or. z(3).le.z(11) .or.
     *   x(4).le.x(12) .or. y(4).ge.y(12) .or. z(4).le.z(12) .or.
     *   x(5).ge.x(13) .or. y(5).ge.y(13) .or. z(5).ge.z(13) .or.
     *   x(6).ge.x(14) .or. y(6).le.y(14) .or. z(6).ge.z(14) .or.
     *   x(7).le.x(15) .or. y(7).le.y(15) .or. z(7).ge.z(15) .or.
     *   x(8).le.x(16) .or. y(8).ge.y(16) .or. z(8).ge.z(16)) stop 270
c
      xmin = xmin - wlenx
      ymin = ymin - wleny
      zmin = zmin - wlenz
      xmax = xmax + wlenx
      ymax = ymax + wleny
      zmax = zmax + wlenz
c
      if(naddl.eq.2) go to 155
c
c     compute other points in grid on surface of polyhedron
c
      naddm = naddl-2
      xint = (xmax-xmin)/dble(naddl-1)
      yint = (ymax-ymin)/dble(naddl-1)
      zint = (zmax-zmin)/dble(naddl-1)
      ng = 16
c
      do 119 i = 1, naddm
         xcor = xmin + dble(i)*xint
         ng = ng + 4
         x(ng-3) = xcor
         y(ng-3) = ymin
         z(ng-3) = zmin
         x(ng-2) = xcor
         y(ng-2) = ymin
         z(ng-2) = zmax
         x(ng-1) = xcor
         y(ng-1) = ymax
         z(ng-1) = zmin
         x(ng) = xcor
         y(ng) = ymax
         z(ng) = zmax
  119 continue
c
      do 120 i = 1, naddm
         ycor = ymin + dble(i)*yint
         ng = ng + 4
         y(ng-3) = ycor
         z(ng-3) = zmin
         x(ng-3) = xmin
         y(ng-2) = ycor
         z(ng-2) = zmin
         x(ng-2) = xmax
         y(ng-1) = ycor
         z(ng-1) = zmax
         x(ng-1) = xmin
         y(ng) = ycor
         z(ng) = zmax
         x(ng) = xmax
  120 continue
c
      do 121 i = 1, naddm
         zcor = zmin + dble(i)*zint
         ng = ng + 4
         z(ng-3) = zcor
         x(ng-3) = xmin
         y(ng-3) = ymin
         z(ng-2) = zcor
         x(ng-2) = xmin
         y(ng-2) = ymax
         z(ng-1) = zcor
         x(ng-1) = xmax
         y(ng-1) = ymin
         z(ng) = zcor
         x(ng) = xmax
         y(ng) = ymax
  121 continue
c
      do 130 i = 1, naddm
         xcor = xmin + dble(i)*xint
         do 125 j = 1, naddm
            ycor = ymin + dble(j)*yint
            ng = ng + 2
            x(ng-1) = xcor
            y(ng-1) = ycor
            z(ng-1) = zmin
            x(ng) = xcor
            y(ng) = ycor
            z(ng) = zmax
  125    continue
  130 continue
c
      do 140 i = 1, naddm
         ycor = ymin + dble(i)*yint
         do 135 j = 1, naddm
            zcor = zmin + dble(j)*zint
            ng = ng + 2
            y(ng-1) = ycor
            z(ng-1) = zcor
            x(ng-1) = xmin
            y(ng) = ycor
            z(ng) = zcor
            x(ng) = xmax
  135    continue
  140 continue
c
      do 150 i = 1, naddm
         zcor = zmin + dble(i)*zint
         do 145 j = 1, naddm
            xcor = xmin + dble(j)*xint
            ng = ng + 2
            z(ng-1) = zcor
            x(ng-1) = xcor
            y(ng-1) = ymin
            z(ng) = zcor
            x(ng) = xcor
            y(ng) = ymax
  145    continue
  150 continue
c
      if(ng.ne.irec) stop 280
c
  155 continue
c
c     find i9
c
      do 157 i = 9, nv
         if(is(i).ne.0) go to 158
  157 continue
      stop 282
  158 continue
      i9 = i
c
c     divide cube into 12 tetrahedra
c
      icon(1,1) = 0
      icon(2,1) = 5
      icon(3,1) = 9
      icon(4,1) = 2
      icon(5,1) = i9
      icon(6,1) = 7
      icon(7,1) = 3
      icon(8,1) = 4
c
      icon(1,2) = 0
      icon(2,2) = 1
      icon(3,2) = 12
      icon(4,2) = 3
      icon(5,2) = i9
      icon(6,2) = 2
      icon(7,2) = 3
      icon(8,2) = 7
c
      icon(1,3) = 0
      icon(2,3) = 5
      icon(3,3) = 2
      icon(4,3) = 4
      icon(5,3) = i9
      icon(6,3) = 2
      icon(7,3) = 1
      icon(8,3) = 3
c
      icon(1,4) = 0
      icon(2,4) = 3
      icon(3,4) = 11
      icon(4,4) = 6
      icon(5,4) = i9
      icon(6,4) = 5
      icon(7,4) = 1
      icon(8,4) = 2
c
      icon(1,5) = 0
      icon(2,5) = 6
      icon(3,5) = 1
      icon(4,5) = 3
      icon(5,5) = i9
      icon(6,5) = 3
      icon(7,5) = 1
      icon(8,5) = 4
c
      icon(1,6) = 0
      icon(2,6) = 4
      icon(3,6) = 7
      icon(4,6) = 5
      icon(5,6) = i9
      icon(6,6) = 4
      icon(7,6) = 1
      icon(8,6) = 5
c
      icon(1,7) = 0
      icon(2,7) = 6
      icon(3,7) = 8
      icon(4,7) = 9
      icon(5,7) = i9
      icon(6,7) = 8
      icon(7,7) = 4
      icon(8,7) = 5
c
      icon(1,8) = 0
      icon(2,8) = 10
      icon(3,8) = 9
      icon(4,8) = 7
      icon(5,8) = i9
      icon(6,8) = 8
      icon(7,8) = 5
      icon(8,8) = 7
c
      icon(1,9) = 0
      icon(2,9) = 7
      icon(3,9) = 8
      icon(4,9) = 1
      icon(5,9) = i9
      icon(6,9) = 7
      icon(7,9) = 4
      icon(8,9) = 8
c
      icon(1,10) = 0
      icon(2,10) = 11
      icon(3,10) = 12
      icon(4,10) = 8
      icon(5,10) = i9
      icon(6,10) = 7
      icon(7,10) = 5
      icon(8,10) = 6
c
      icon(1,11) = 0
      icon(2,11) = 4
      icon(3,11) = 12
      icon(4,11) = 10
      icon(5,11) = i9
      icon(6,11) = 6
      icon(7,11) = 5
      icon(8,11) = 2
c
      icon(1,12) = 0
      icon(2,12) = 2
      icon(3,12) = 10
      icon(4,12) = 11
      icon(5,12) = i9
      icon(6,12) = 6
      icon(7,12) = 2
      icon(8,12) = 7
c
      tetra = 12
c
      wmen = wmin
      if(.not.delaun .and. irec.gt.8) then
         wmen = wmen - wlenw
         do 160 i=9,irec
            w(i) = wmen
  160    continue
      endif
      if(wmen.lt.wmin) wmin = wmen
      if(wmen.gt.wmax) wmax = wmen
      wman = wmin - 10.0d0
      do 170 i=1,8
         w(i) = wman
  170 continue
c
      is(1) = 5
      is(2) = 12
      is(3) = 1
      is(4) = 9
      is(5) = 10
      is(6) = 12
      is(7) = 12
      is(8) = 9
      is(i9) = 12
c
c     test # of significant figures of nondecimal part of coordinates
c
      wbig = 0.0d0
      if(wbig .lt. dabs(xmax)) wbig = dabs(xmax)
      if(wbig .lt. dabs(xmin)) wbig = dabs(xmin)
      if(wbig .lt. dabs(ymax)) wbig = dabs(ymax)
      if(wbig .lt. dabs(ymin)) wbig = dabs(ymin)
      if(wbig .lt. dabs(zmax)) wbig = dabs(zmax)
      if(wbig .lt. dabs(zmin)) wbig = dabs(zmin)
      wbig = wbig + epz
c     WRITE(*,*)'COORDINATES WBIG=',WBIG
      ibsfig = 0
  180 continue
      ibsfig = ibsfig+1
      wbig = wbig/10.0d0
      if(wbig .ge. 1.0d0) go to 180
      if(ibsfig.gt.9) then
         write(*,*)'Number of significant figures of largest ',
     *             'nondecimal part of'
         write(*,*)'a point coordinate appears to be greater than 9.'
         write(*,*)'Program is terminated.'
         stop 286
      endif
      itsfig = ibsfig + icsfig
c     WRITE(*,*)'ITSFIG=',ITSFIG,' IBSFIG=',IBSFIG,' ICSFIG=',ICSFIG
      if(itsfig.gt.14) then
         write(*,*)' '
         write(*,*)'For this execution of the program the largest ',
     *             'total number of'
         write(*,*)'significant figures ',
     *             'that a coordinate requires appears to be ',itsfig
         write(*,*)'Since the maximum allowed is 14, the number of ',
     *             'significant'
         write(*,*)'figures of the decimal part of the coordinates ',
     *             'for this run is '
         write(*,*)'decreased accordingly.'
         icsfig = 14 - ibsfig
         write(*,*)' '
         write(*,*)'Now icfig = ',icsfig
         write(*,*)' '
      endif
c
c     if a Regular tetrahedralization test # of significant figures
c     of nondecimal part of weights
c
      if(delaun) go to 200
      wbig = 0.0d0
      if(wbig .lt. dabs(wmax)) wbig = dabs(wmax)
      if(wbig .lt. dabs(wmin)) wbig = dabs(wmin)
      wbig = wbig + epz
c     WRITE(*,*)'COORDINATES WBIG=',WBIG
      ibsfig = 0
  190 continue
      ibsfig = ibsfig+1
      wbig = wbig/10.0d0
      if(wbig .ge. 1.0d0) go to 190
      if(ibsfig.gt.9) then
         write(*,*)'Number of significant figures of largest ',
     *             'nondecimal part of'
         write(*,*)'a weight appears to be greater than 9.'
         write(*,*)'Program is terminated.'
         stop 288
      endif
      itsfig = ibsfig + iwsfig
c     WRITE(*,*)'ITSFIG=',ITSFIG,' IBSFIG=',IBSFIG,' IWSFIG=',IWSFIG
      if(itsfig.gt.14) then
         write(*,*)' '
         write(*,*)'For this execution of the program the largest ',
     *             'total number of'
         write(*,*)'significant figures ',
     *             'that a weight requires appears to be ',itsfig
         write(*,*)'Since the maximum allowed is 14, the number of ',
     *             'significant'
         write(*,*)'figures of the decimal part of the weights for ',
     *             'this run is '
         write(*,*)'decreased accordingly.'
         iwsfig = 14 - ibsfig
         icsfi2 = 2*icsfig
         irsfig = icsfi2 - iwsfig
         if(irsfig.gt.9) then
            write(*,*)'In order to make this number compatible with ',
     *                'that of the'
            write(*,*)'coordinates, the latter is also decreased ',
     *                'accordingly.'
            icsfig = (iwsfig+9)/2
         endif
         write(*,*)' '
         write(*,*)'Now icfig = ',icsfig,' iwfig = ',iwsfig
         write(*,*)' '
      endif
  200 continue
c
c     test number of significant figures of decimal part of coordinates
c
      if(icsfig.lt.0 .or. icsfig.gt.9) then
         write(*,*)'Number of significant figures of decimal'
         write(*,*)'part of coordinates is out of range.'
         write(*,*)'Program is terminated.'
         stop 290
      endif
      isclu = 1
      dscle = 1.0d0
      if(icsfig.eq.0) go to 220
      do 210 i = 1, icsfig
         isclu = 10*isclu
         dscle = 10.0d0*dscle
  210 continue
  220 continue
      if(iabs(isclu).ge.mfull) stop 295
      call decomp(isclp, isgcl, isclu, mhalf)
      if(isgcl.ne.1) stop 297
c
c     test lengths of x, y, z-coordinates, shift and make them integers
c
      dfull = dble(mfull)
      dfill=dfull/dscle
      do 235 i = 1, nv
         ix2(i) = 0
         iy2(i) = 0
         iz2(i) = 0
         if(dabs(x(i)).lt.dfill) then
            ix(i) = idnint(dscle*x(i))
            if(iabs(ix(i)).lt.mfull) then
               x(i) = dble(ix(i))/dscle
               go to 225
            endif
         endif
         if(dabs(x(i)).ge.dfull) stop 305
         ix(i) = idint(x(i))
         if(iabs(ix(i)).ge.mfull) stop 310
         decml = (x(i) - dint(x(i)))*dscle
         if(dabs(decml).ge.dfull) stop 312
         ix2(i) = idnint(decml)
         if(iabs(ix2(i)).ge.mfull) stop 315
         if(iabs(ix2(i)).eq.0) then
            x(i) = dble(ix(i))
            ix2(i) = mfull
         else
            x(i) = dble(ix(i)) + (dble(ix2(i))/dscle)
         endif
  225    continue
         if(dabs(y(i)).lt.dfill) then
            iy(i) = idnint(dscle*y(i))
            if(iabs(iy(i)).lt.mfull) then
               y(i) = dble(iy(i))/dscle
               go to 230
            endif
         endif
         if(dabs(y(i)).ge.dfull) stop 320
         iy(i) = idint(y(i))
         if(iabs(iy(i)).ge.mfull) stop 325
         decml = (y(i) - dint(y(i)))*dscle
         if(dabs(decml).ge.dfull) stop 327
         iy2(i) = idnint(decml)
         if(iabs(iy2(i)).ge.mfull) stop 330
         if(iabs(iy2(i)).eq.0) then
            y(i) = dble(iy(i))
            iy2(i) = mfull
         else
            y(i) = dble(iy(i)) + (dble(iy2(i))/dscle)
         endif
  230    continue
         if(dabs(z(i)).lt.dfill) then
            iz(i) = idnint(dscle*z(i))
            if(iabs(iz(i)).lt.mfull) then
               z(i) = dble(iz(i))/dscle
               go to 235
            endif
         endif
         if(dabs(z(i)).ge.dfull) stop 335
         iz(i) = idint(z(i))
         if(iabs(iz(i)).ge.mfull) stop 340
         decml = (z(i) - dint(z(i)))*dscle
         if(dabs(decml).ge.dfull) stop 342
         iz2(i) = idnint(decml)
         if(iabs(iz2(i)).ge.mfull) stop 345
         if(iabs(iz2(i)).eq.0) then
            z(i) = dble(iz(i))
            iz2(i) = mfull
         else
            z(i) = dble(iz(i)) + (dble(iz2(i))/dscle)
         endif
  235 continue
c
c     if a Regular tetrahedralization test number of significant figures
c     of decimal part of weights, test lengths of weights, shift and
c     make them integers
c
      if(delaun) go to 300
      icsfi2 = 2*icsfig
      irsfig = icsfi2 - iwsfig
      if(iwsfig.lt.0.or.iwsfig.gt.9 .or. irsfig.lt.0.or.irsfig.gt.9)then
         write(*,*)'Either number of significant figures of decimal'
         write(*,*)'part of weights is out of range or it is not'
         write(*,*)'compatible with that of point coordinates.'
         write(*,*)'Program is terminated.'
         stop 350
      endif
      isclu = 1
      dscli = 1.0d0
      if(iwsfig.eq.0) go to 250
      do 240 i = 1, iwsfig
         isclu = 10*isclu
         dscli = 10.0d0*dscli
  240 continue
  250 continue
      if(iabs(isclu).ge.mfull) stop 360
      call decomp(isclw, isgcl, isclu, mhalf)
      if(isgcl.ne.1) stop 363
      dfill = dfull/dscli
      do 260 i = 1, nv
         iw2(i) = 0
         if(dabs(w(i)).lt.dfill) then
            iw(i) = idnint(dscli*w(i))
            if(iabs(iw(i)).lt.mfull) then
               w(i) = dble(iw(i))/dscli
               go to 260
            endif
         endif
         if(dabs(w(i)).ge.dfull) stop 365
         iw(i) = idint(w(i))
         if(iabs(iw(i)).ge.mfull) stop 370
         decml = (w(i) - dint(w(i)))*dscli
         if(dabs(decml).ge.dfull) stop 372
         iw2(i) = idnint(decml)
         if(iabs(iw2(i)).ge.mfull) stop 375
         if(iabs(iw2(i)).eq.0) then
            w(i) = dble(iw(i))
            iw2(i) = mfull
         else
            w(i) = dble(iw(i)) + (dble(iw2(i))/dscli)
         endif
  260 continue
      isclu = 1
      if(irsfig.eq.0) go to 290
      do 270 i = 1, irsfig
         isclu = 10*isclu
  270 continue
  290 continue
      if(iabs(isclu).ge.mfull) stop 385
      call decomp(isclr, isgcl, isclu, mhalf)
      if(isgcl.ne.1) stop 390
  300 continue
c
c     get cube corner directions in their integer form
c
      do 320 i = 1,8
         ixc(i) = idnint(xc(i))
         iyc(i) = idnint(yc(i))
         izc(i) = idnint(zc(i))
  320 continue
c
c     add all points to tetrahedralization
c
      i10 = i9+1
      if(nv .lt. i10) go to 400
c     write(*,*)' '
c     write(*,*)'Adding all points to tetrahedralization ...'
      ITCHK = 0
      ITCN1 = 0
      ITCN2 = 0
      ITCN3 = 0
      ITCN4 = 0
      issin = i9
      do 380 k = i10, nv
c         if(k.le.(k/1000)*1000)
c     *   write(*,*)'Number of points processed = ',k
         if(is(k).eq.0) go to 380
         call pntins(x, y, z, w, ix, iy, iz, iw, ix2, iy2, iz2, iw2,
     *               icon, is, ifl, id, ih, ihn, k, nvmax, nhmax, tetra,
     *               mxlook, ixc, iyc, izc, issin, iftal, delaun,
     *               flphis, mhalf, mfull, isclp, isclw, isclr, epz,
     *               ITCHK,ITCN1,ITCN2,ITCN3,ITCN4)
  380 continue
c     WRITE(*,*)' '
c     WRITE(*,*)'POINTINS ITCHK=',ITCHK
c     WCPNT = REAL(ITCHK)/REAL(NV-9)
c     WRITE(*,*)' '
c     WRITE(*,*)'DURING INSERTION AVERAGE NUMBER OF ',
c    *          'TETRAHEDRA CHECKED PER POINT WAS '
c     WRITE(*,*) WCPNT
c     WRITE(*,*)' '
c     WRITE(*,*)'ITCN1234=',ITCN1,ITCN2,ITCN3,ITCN4
  400 continue
c
c     test redundant points
c
      iredx = 0
      if(delaun) go to 1000
      if(nv.lt.10) go to 1000
      do 420 k = 9, nv
         if(is(k).le.0) go to 420
  420 continue
      if(redchk) then
         write(*,*)' '
         write(*,*)'Testing redundant points ...'
      endif
      ITCHK = 0
      ITCN1 = 0
      ITCN2 = 0
      ITCN3 = 0
      ITCN4 = 0
      iredu = 0
      do 500 k = 9, nv
         if(is(k).ge.0) go to 450
         iredu = iredu+1
         if(.not.redchk) go to 500
         call pntred(x, y, z, w, ix, iy, iz, iw, ix2, iy2, iz2, iw2,
     *               icon, is, id, k, tetra, ixc, iyc, izc,
     *               iredx, issin, iftal, delaun, flphis, mhalf,
     *               mfull, isclp, isclw, isclr, epz,
     *               ITCHK,ITCN1,ITCN2,ITCN3,ITCN4)
  450    continue
         if(is(k).le.0) go to 500
         issin = k
  500 continue
c     WRITE(*,*)' '
c     WRITE(*,*)'POINTRED ITCHK=',ITCHK
c     WRITE(*,*)' '
c     WRITE(*,*)'ITCN1234=',ITCN1,ITCN2,ITCN3,ITCN4
c     write(*,*)' '
c     write(*,*)'Number of redundant points = ',iredu
c
 1000 continue
      return
      end
*PNTINS
c 
c     This subroutine will find location of new point
c
c This routine also calls routine 'sphere' for the purpose of optimizing
c for the locally Regular property
c
      subroutine pntins(xi, yi, zi, wi, x, y, z, w, x2, y2, z2, w2,
     *                  icon, is, ifl, id, ih, ihn, k, nvmax, nhmax,
     *                  tetra, mxlook, xc, yc, zc, issin, iftal, delaun,
     *                  flphis, mhalf, mfull, isclp, isclw, isclr, epz,
     *                  ITCHK,ITCN1,ITCN2,ITCN3,ITCN4)
c
      double precision xi(*), yi(*), zi(*), wi(*)
      integer x(*), y(*), z(*), w(*)
      integer x2(*), y2(*), z2(*), w2(*)
      integer icon(8,*), is(*), ifl(*), id(*), ih(*)
      integer xc(*), yc(*), zc(*)
      integer ihn, k, nvmax, nhmax, mxlook, issin, iftal
      integer isclp(*), isclw(*), isclr(*), curr, side1, side2, tetra
      double precision epz
      logical delaun, flphis
      integer mhalf, mfull, itype, look, newtts, ired
      integer itchk, itcn1, itcn2, itcn3, itcn4
c
      if(flphis) then
         itype = -1
         look = 1
         call gettet(itype, k, xi, yi, zi, x, y, z, x2, y2, z2, icon,
     *               curr, side1, side2, xc, yc, zc, mhalf, mfull,
     *               isclp, epz, ITCHK)
c
   50    continue
         if (icon(5,curr) .lt. 0) then
            call lkdown(icon, curr, xi, yi, zi, x, y, z, x2, y2, z2,
     *                  itype, k, side1, side2, xc, yc, zc, mhalf,
     *                  mfull, isclp, epz, ITCHK)
            look = look + 1
            if(itype.eq.-1) stop 410
            goto 50
         endif
      else
         call shishk(xi, yi, zi, x, y, z, x2, y2, z2, is, icon, id,
     *               issin, k, side1, side2, curr, iftal, itype, tetra,
     *               xc, yc, zc, mhalf, mfull, isclp, epz, ITCHK)
         if(itype.eq.-1) stop 420
      endif
c
      if (itype .eq. 1) then
         call vrtins (k, w, w2, icon, nvmax, tetra, curr, is, id, iftal,
     *                side1, ifl, newtts, ired, delaun, flphis, mhalf,
     *                mfull, isclw)
         ITCN1 = ITCN1+1
      elseif (itype .eq. 2) then
         call intins (k, xi, yi, zi, wi, x, y, z, w, x2, y2, z2, w2,
     *                icon, ih, ihn, nhmax, nvmax, tetra, curr, is, ifl,
     *                newtts, ired, delaun, flphis, mhalf, mfull,
     *                isclp, isclw, isclr, epz)
         ITCN2 = ITCN2+1
      elseif (itype .eq. 3) then
         call edgins (k, x, y, z, w, x2, y2, z2, w2, icon, ih, ihn,
     *                nhmax, nvmax, tetra, curr, is, side1, side2,
     *                ifl, newtts, ired, delaun, flphis, mhalf, mfull,
     *                isclp, isclw, isclr)
         ITCN3 = ITCN3+1
      elseif (itype .eq. 4) then
         call sidins (k, x, y, z, w, x2, y2, z2, w2, icon, ih, ihn,
     *                nhmax, nvmax, tetra, curr, is, side1, ifl,
     *                newtts, ired, delaun, flphis, mhalf, mfull,
     *                isclp, isclw, isclr)
         ITCN4 = ITCN4+1
      else
         stop 430
      endif
      if(ired .eq. 1) go to 1000
      issin = k
c
c     optimize for Regular/Delaunay property
c
      call sphere(k, icon, ifl, ih, ihn, nhmax, newtts, xi, yi, zi, wi,
     *            x, y, z, w, x2, y2, z2, w2, tetra, is, nvmax, xc, yc,
     *            zc, delaun, flphis, mhalf, mfull, isclp, isclw, isclr,
     *            epz)
 1000 continue
c
      if (look .gt. mxlook) mxlook = look
c
      return
      end
*PNTRED
c 
c This subroutine will test a redundant point for redundancy
c
c
      subroutine pntred(xi, yi, zi, wi, x, y, z, w, x2, y2, z2, w2,
     *                  icon, is, id, k, tetra, xc, yc, zc,
     *                  idmax, issin, iftal, delaun, flphis, mhalf,
     *                  mfull, isclp, isclw, isclr, epz,
     *                  ITCHK,ITCN1,ITCN2,ITCN3,ITCN4)
c
      double precision xi(*), yi(*), zi(*), wi(*)
      integer x(*), y(*), z(*), w(*)
      integer x2(*), y2(*), z2(*), w2(*)
      integer icon(8,*), is(*), id(*)
      integer xc(*), yc(*), zc(*)
      integer k, idmax, issin, iftal, mhalf, mfull
      integer curr, side1, side2, tetra, a, b, c, d
      integer isclp(*), isclw(*), isclr(*)
      integer site1, site2, site3, site4, fndsit
      double precision epz, tdist, xctr, yctr, zctr
      logical delaun, flphis
      integer itype, itide, ipossi, i
      integer itchk, itcn1, itcn2, itcn3, itcn4
c
      if(flphis) then
         itype = -1
         call gettet(itype, k, xi, yi, zi, x, y, z, x2, y2, z2, icon,
     *               curr, side1, side2, xc, yc, zc, mhalf, mfull,
     *               isclp, epz, ITCHK)
c
   50    continue
         if (icon(5,curr) .lt. 0) then
            call lkdown(icon, curr, xi, yi, zi, x, y, z, x2, y2, z2,
     *                  itype, k, side1, side2, xc, yc, zc, mhalf,
     *                  mfull, isclp, epz, ITCHK)
            if(itype.eq.-1) stop 510
            goto 50
         endif
      else
         call shishk(xi, yi, zi, x, y, z, x2, y2, z2, is, icon, id,
     *               issin, k, side1, side2, curr, iftal, itype, tetra,
     *               xc, yc, zc, mhalf, mfull, isclp, epz, ITCHK)
         if(itype.eq.-1) stop 515
      endif
c
      if (itype .eq. 1) then
         ITCN1 = ITCN1+1
         site1 = icon(side1+4,curr)
c        itide = 1
c        if(w(k) .gt. w(site1)) itide =-1
         call iwsign(w, w2, site1, k, mhalf, mfull, isclw, itide)
      elseif (itype .eq. 2) then
         ITCN2 = ITCN2+1
         a = icon(5,curr)
         b = icon(6,curr)
         c = icon(7,curr)
         d = icon(8,curr)
         if(a.le.8 .or. b.le.8 .or. c.le.8 .or. d.le.8) stop 520
         call ctrad(xi, yi, zi, wi, xctr, yctr, zctr, a, b, c, d,
     *              epz, delaun, ipossi)
         if(ipossi.eq.1) go to 60
         call bisphr(xi, yi, zi, wi, k, a, epz, xctr, yctr, zctr, tdist,
     *               delaun, ipossi)
         if(ipossi.eq.1) go to 60
         itide = 1
         if(tdist.gt.0.0d0) itide = -1
         go to 1000
   60    continue
         call iqsign(x, y, z, w, x2, y2, z2, w2, a, b, c, d, k,
     *               mhalf, mfull, isclp, isclw, isclr, delaun, itide)
      elseif (itype .eq. 3) then
         ITCN3 = ITCN3+1
         fndsit = 0
         do 100 i = 5, 8
            if (fndsit .eq. 0) then
               if (i .eq. (side1+4) .or. i .eq. (side2+4)) then
                  goto 100
               else
                  site1 = icon(i,curr)
                  fndsit = 1
               endif
            else 
               if (i .eq. (side1+4) .or. i .eq. (side2+4)) then
                  goto 100
               else
                  site2 = icon(i,curr)
                  goto 150
               endif
            endif
  100    continue
         stop 530
  150    continue
c
         if(site1.le.8 .or. site2.le.8) stop 540
         call iqsig1(x, y, z, w, x2, y2, z2, w2, site1, site2, k,
     *               mhalf, mfull, isclp, isclw, isclr, delaun, itide)
      elseif (itype .eq. 4) then
         ITCN4 = ITCN4+1
         site1 = icon(side1+4,curr)
         call sitord(icon, site1, curr)
c
         site2 = icon(6,curr)
         site3 = icon(7,curr)
         site4 = icon(8,curr)
c
         if(site2.le.8 .or. site3.le.8 .or. site4.le.8) stop 550
         call iqsig2(x, y, z, w, x2, y2, z2, w2, site2, site3, site4, k,
     *               mhalf, mfull, isclp, isclw, isclr, delaun, itide)
      else
         stop 560
      endif
c
 1000 continue
      if(itide.lt.0) idmax = idmax+1
c
      return
      end
*GETTET
c
c This subroutine will test each of the 1st 12 tetrahedra to find
c where new point is located.  It'll do so by calling 'gette2'.
c
      subroutine gettet(itype, k, xi, yi, zi, x, y, z, x2, y2, z2, icon,
     *                  curr, side1, side2, xc, yc, zc, mhalf, mfull,
     *                  isclp, epz, ITCHK)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer icon(8,*), curr, a, b, c, d, side1, side2, flag
      integer isclp(*), mhalf, mfull
      double precision epz
      integer itype, k, itchk
c
      do 100 curr = 1, 12
         ITCHK = ITCHK+1
         a = iabs(icon(5, curr))
         b = icon(6, curr)
         c = icon(7, curr)
         d = icon(8, curr)
c
         flag = icon(5,curr)
         icon(5,curr)=a
c
         call vrtord(icon, curr, a, b, c, d)
         if(flag.lt.0) icon(5,curr)=-a
         if(a.le.8.or.b.gt.8.or.c.gt.8.or.d.gt.8) stop 610
c
         call gette2(a, b, c, d, xi, yi, zi, x, y, z, x2, y2, z2,
     *               itype, k, side1, side2, flag, xc, yc, zc,
     *               mhalf, mfull, isclp, epz)
         if (itype .ne. -1) goto 200
c
  100 continue
      stop 620
  200 continue
      return
      end
*LKDOWN
c
c This subroutine will traverse thru children of curr until point
c is found to be in one of them
c
      subroutine lkdown(icon, curr, xi, yi, zi, x, y, z, x2, y2, z2,
     *                  itype, k, side1, side2, xc, yc, zc, mhalf,
     *                  mfull, isclp, epz, ITCHK)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*)
      integer icon(8,*), xc(*), yc(*), zc(*)
      double precision epz
      integer isclp(*), itype, k, mhalf, mfull, i, newcur, itchk
      integer curr, side1, side2, a, b, c, d, flag
c
c     test children of current tetrahedron
c
      itype = -1
      do 100 i = 1, 4
         newcur = icon(i,curr)
         if (newcur .le. 0) goto 100
         ITCHK = ITCHK+1
c
         a = iabs(icon(5,newcur))
         b = icon(6,newcur)
         c = icon(7,newcur)
         d = icon(8,newcur)
c
         flag = icon(5,newcur)
         icon(5,newcur)=a
c
         call vrtord(icon, newcur, a, b, c, d)
         if(flag.lt.0) icon(5,newcur)=-a
         if(a.le.8) stop 710
c
         call gette2(a, b, c, d, xi, yi, zi, x, y, z, x2, y2, z2,
     *               itype, k, side1, side2, flag, xc, yc, zc,
     *               mhalf, mfull, isclp, epz)
         if (itype .eq. -1) goto 100
         curr = newcur
         goto 1000
  100 continue
c
 1000 continue
      return
      end
*GETTE2
c
c This subroutine will check for each tetra, if the point is equal to an
c existing vertex, inside (interior, edge, side), or outside curr tetra.
c
      subroutine gette2(a, b, c, d, xi, yi, zi, x, y, z, x2, y2, z2,
     *                  itype, k, side1, side2, flag, xc, yc, zc,
     *                  mhalf, mfull, isclp, epz)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*)
      integer xc(*), yc(*), zc(*)
      double precision epz
      integer isclp(*), itype, k, mhalf, mfull, ifn, ipout
      integer iside(4), a, b, c, d, side1, side2, flag
c
c     determine position of point k relative to facets of tetrahedron
c
      if(b.le.8 .or. c.le.8 .or. d.le.8) go to 100
c
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, b, d, c,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(1) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, a, c, d,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(2) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, a, d, b,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(3) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, a, b, c,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(4) = ipout
c
c     point k is not in tetrahedron
c
   50 continue
      if(iside(1).lt.0 .or. iside(2).lt.0 .or. iside(3).lt.0 .or.
     *   iside(4).lt.0) go to 1000
c
      if(flag.lt.0) then
         itype = 0
         go to 1000
      endif
c
c     point k is in the interior of tetrahedron
c
      if(iside(1).gt.0 .and. iside(2).gt.0 .and. iside(3).gt.0 .and.
     *   iside(4).gt.0) then
         itype = 2
         go to 1000
      endif
c
c     unacceptable situation
c
      if(iside(1).eq.0 .and. iside(2).eq.0 .and. iside(3).eq.0 .and.
     *   iside(4).eq.0) stop 805
c
c     point k is a vertex of tetrahedron
c
      if(iside(1).eq.0 .and. iside(2).eq.0 .and. iside(3).eq.0) then
         itype = 1
         side1 = 4
         go to 1000
      elseif(iside(1).eq.0 .and. iside(2).eq.0 .and. iside(4).eq.0) then
         itype = 1
         side1 = 3
         go to 1000
      elseif(iside(1).eq.0 .and. iside(3).eq.0 .and. iside(4).eq.0) then
         itype = 1
         side1 = 2
         go to 1000
      elseif(iside(2).eq.0 .and. iside(3).eq.0 .and. iside(4).eq.0) then
         itype = 1
         side1 = 1
         go to 1000
      endif
c
c     point k is in the interior of an edge of tetrahedron
c
      if (iside(1).eq.0 .and. iside(2).eq.0) then
         itype = 3 
         side1 = 1
         side2 = 2
         go to 1000
      elseif (iside(1).eq.0 .and. iside(3).eq.0) then
         itype = 3 
         side1 = 1
         side2 = 3
         go to 1000
      elseif (iside(1).eq.0 .and. iside(4).eq.0) then
         itype = 3 
         side1 = 1
         side2 = 4
         go to 1000
      elseif (iside(2).eq.0 .and. iside(3).eq.0) then
         itype = 3 
         side1 = 2
         side2 = 3
         go to 1000
      elseif (iside(2).eq.0 .and. iside(4).eq.0) then
         itype = 3 
         side1 = 2
         side2 = 4
         go to 1000
      elseif (iside(3).eq.0 .and. iside(4).eq.0) then
         itype = 3 
         side1 = 3
         side2 = 4
         go to 1000
      endif
c
c     point k is in the interior of a facet of tetrahedron
c
      itype = 4
      if (iside(1) .eq. 0) then
         side1 = 1
      elseif (iside(2) .eq. 0) then
         side1 = 2
      elseif (iside(3) .eq. 0) then
         side1 = 3
      elseif (iside(4) .eq. 0) then 
         side1 = 4
      else
         stop 807
      endif
      go to 1000
c
c     there is at least one artificial vertex
c
  100 continue
      if(b.le.8) then
         iside(1) = 1
         go to 120
      elseif(d.le.8.and.c.le.8) then
         call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *               d, c, k, b, mhalf, mfull, isclp, ipout)
         iside(1) = ipout
         if(iside(1).ne.0) go to 120
         call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *               d, c, k, b, mhalf, mfull, isclp, ipout)
         iside(1) = ipout
         go to 120
      elseif(d.le.8) then
         ifn = 0
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *               b, c, k, d, mhalf, mfull, isclp, ifn, ipout)
         iside(1) = ipout
      else
         ifn = 1
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *               b, d, k, c, mhalf, mfull, isclp, ifn, ipout)
         iside(1) = ipout
      endif
      if(iside(1).ne.0) go to 120
      call ipsign(x, y, z, x2, y2, z2, b, d, c, k, mhalf,
     *            mfull, isclp, ipout)
      iside(1) = ipout
  120 continue
c
      if(c.le.8.and.d.le.8) then
         call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *               c, d, k, a, mhalf, mfull, isclp, ipout)
         iside(2) = ipout
         if(iside(2).ne.0) go to 140
         call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *               c, d, k, a, mhalf, mfull, isclp, ipout)
         iside(2) = ipout
         go to 140
      elseif(c.le.8) then
         ifn = 0
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *               a, d, k, c, mhalf, mfull, isclp, ifn, ipout)
         iside(2) = ipout
      else
         ifn = 1
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *               a, c, k, d, mhalf, mfull, isclp, ifn, ipout)
         iside(2) = ipout
      endif
      if(iside(2).ne.0) go to 140
      call ipsign(x, y, z, x2, y2, z2, a, c, d, k, mhalf,
     *            mfull, isclp, ipout)
      iside(2) = ipout
  140 continue
c
      if(d.gt.8) then
         call ipsign(x, y, z, x2, y2, z2, a, d, b, k, mhalf,
     *               mfull, isclp, ipout)
         iside(3) = ipout
         go to 160
      elseif(b.le.8) then
         call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *               d, b, k, a, mhalf, mfull, isclp, ipout)
         iside(3) = ipout
         if(iside(3).ne.0) go to 160
         call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *               d, b, k, a, mhalf, mfull, isclp, ipout)
         iside(3) = ipout
         go to 160
      else
         ifn = 0
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *               a, b, k, d, mhalf, mfull, isclp, ifn, ipout)
         iside(3) = ipout
      endif
      if(iside(3).ne.0) go to 160
      call ipsign(x, y, z, x2, y2, z2, a, d, b, k, mhalf,
     *            mfull, isclp, ipout)
      iside(3) = ipout
  160 continue
c
      if(c.gt.8) then
         call ipsign(x, y, z, x2, y2, z2, a, b, c, k, mhalf,
     *               mfull, isclp, ipout)
         iside(4) = ipout
         go to 180
      elseif(b.le.8) then
         call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *               b, c, k, a, mhalf, mfull, isclp, ipout)
         iside(4) = ipout
         if(iside(4).ne.0) go to 180
         call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *               b, c, k, a, mhalf, mfull, isclp, ipout)
         iside(4) = ipout
         go to 180
      else
         ifn = 1
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *               a, b, k, c, mhalf, mfull, isclp, ifn, ipout)
         iside(4) = ipout
      endif
      if(iside(4).ne.0) go to 180
      call ipsign(x, y, z, x2, y2, z2, a, b, c, k, mhalf,
     *            mfull, isclp, ipout)
      iside(4) = ipout
  180 continue
      go to 50
c
 1000 continue
      return
      end
*SHISHK
c
c     shishkebab routines
c
c     subroutine shishk to -
c
c     move from a vertex in the tetrahedralization to a tetrahedron
c     that contains a point, and identify the type of location of
c     the point with respect to the tetrahedron
c     
      subroutine shishk(xi, yi, zi, x, y, z, x2, y2, z2, is, icon, id,
     *                  ileft, k, side1, side2, iscur, iftal, itype,
     *                  ivnxt, xc, yc, zc, mhalf, mfull, isclp, epz,
     *                  ITCHK)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*)
      integer x2(*), y2(*), z2(*)
      integer is(*), icon(8,*), id(*)
      double precision epz
      integer xc(*), yc(*), zc(*)
      integer isclp(*), ileft, k, iftal, itype, ivnxt, mhalf, mfull
      integer a, b, c, d, side1, side2, site0, site1
      integer iscur, isini, imist, isadj, ilift, islst, itchk, i
c
c     reinitialize array id if necessary
c
      if(iftal.gt.10000000) then
         iftal = 0
         do 50 i = 1, ivnxt
            id(i) = 0
   50    continue
      endif
c
      if(ileft .le. 8) stop 911
      a = ileft
  100 continue
c
c     find tetrahedron with point a as a vertex for which the ray with
c     origin point a and through point k intersects the facet of the
c     tetrahedron opposite to point a
c
      itype = 0
      iftal = iftal + 1
      iscur = is(a)
      if(iscur.le.0.or.iscur.gt.ivnxt) stop 912
      isini = iscur
c
c     reorder isini so that vertex a equals icon(5,isini)
c
      call sitord(icon, a, isini)
c
c     test current facet
c
  400 continue
      b = icon(6,iscur)
      c = icon(7,iscur)
      d = icon(8,iscur)
      id(iscur) = iftal
c
      ITCHK = ITCHK+1
      call fctest(xi, yi, zi, x, y, z, x2, y2, z2, xc, yc, zc, itype, k,
     *            imist, a, b, c, d, side1, side2, mhalf, mfull, isclp,
     *            epz)
      if(itype .gt. 0) go to 9000
      if(itype .eq. 0) go to 500
      if(itype .eq. -2) go to 1100
      if(itype .eq. -3) then
         a = imist
         go to 100
      elseif(itype .eq. -4) then
         site0 = a
         site1 = imist
         go to 2000
      else
         stop 913
      endif
c
c     obtain next tetrahedron with point a as a vertex
c
  500 continue
      isadj = iabs(icon(2,iscur))
      if(isadj.le.0.or.isadj.gt.ivnxt) stop 914
      if(id(isadj) .eq. iftal) go to 600
      ilift = icon(8,iscur)
      go to 900
  600 continue
      isadj = iabs(icon(3,iscur))
      if(isadj.le.0.or.isadj.gt.ivnxt) stop 915
      if(id(isadj) .eq. iftal) go to 700
      ilift = icon(6,iscur)
      go to 900
  700 continue
      isadj = iabs(icon(4,iscur))
      if(isadj.le.0.or.isadj.gt.ivnxt) stop 916
      if(iscur .eq. isini) go to 800
      if(iabs(icon(3,isadj)) .eq. iscur) then
          iscur = isadj
          go to 700
      elseif(iabs(icon(2,isadj)) .eq. iscur) then
          iscur = isadj
          go to 600
      elseif(iabs(icon(4,isadj)) .eq. iscur) then
          if(isadj .ne. isini) stop 917
          go to 1000
      else
          stop 918
      endif
  800 continue
      if(id(isadj) .eq. iftal) go to 1000
      ilift = icon(7,iscur)
c
c     reorder isadj so that a equals icon(5,isadj) and ilift
c     equals icon(6,isadj)
c
  900 continue
      call reordr(icon, a, ilift, isadj)
      iscur = isadj
      go to 400
c
c     can not find intersected tetrahedron
c
 1000 continue
      stop 919
c
c     obtain next tetrahedron along line segment as it crosses a facet
c
 1100 continue
      islst = iscur
      isadj = iabs(icon(1,iscur))
      if(isadj.le.0.or.isadj.gt.ivnxt) stop 921
      iscur = isadj
      if(iabs(icon(1,iscur)) .eq. islst) then
          ilift = icon(5,iscur)
      elseif(iabs(icon(2,iscur)) .eq. islst) then
          ilift = icon(6,iscur)
      elseif(iabs(icon(3,iscur)) .eq. islst) then
          ilift = icon(7,iscur)
      elseif(iabs(icon(4,iscur)) .eq. islst) then
          ilift = icon(8,iscur)
      else
          stop 922
      endif
c
c     obtain opposite facet of tetrahedron intersected by line
c     segment
c
      ITCHK = ITCHK+1
      call fcfind(xi, yi, zi, x, y, z, x2, y2, z2, xc, yc, zc, itype,
     *        ileft, k, ilift, imist, b, c, d, side1, side2, mhalf,
     *        mfull, isclp, epz)
      if(itype .gt. 0) then
         call reordr(icon, ilift, b, iscur)
         go to 9000
      elseif(itype .eq. -2) then
         call sitord(icon, imist, iscur)
         b = icon(6,iscur)
         c = icon(7,iscur)
         d = icon(8,iscur)
         go to 1100
      elseif(itype .eq. -3) then
         a = ilift
         go to 100
      elseif(itype .eq. -4) then
         if(imist.eq.b)then
            site0 = c
         elseif(imist.eq.c)then
            site0 = d
         else
            site0 = b
         endif
         site1 = imist
         go to 2000
      else
         stop 923
      endif
c
c     obtain next tetrahedron along line segment as it crosses an edge
c
 2000 continue
      call fcedge(x, y, z, x2, y2, z2, xc, yc, zc, itype, ileft, k,
     *            icon, iscur, imist, ivnxt, site0, site1, side1,
     *            side2, mhalf, mfull, isclp, ITCHK)
      if(itype .gt. 0) go to 9000
      if(itype .eq. -2) then
         call sitord(icon, imist, iscur)
         b = icon(6,iscur)
         c = icon(7,iscur)
         d = icon(8,iscur)
         go to 1100
      elseif(itype .eq. -3) then
         a = imist
         go to 100
      elseif(itype.eq.-4) then
         if(imist.eq.site1)then
            site0 = icon(7,iscur)
         elseif(imist.eq.site0)then
            site0 = site1
         endif
         site1 = imist
         go to 2000
      else
         stop 924
      endif
c
 9000 continue
      return
      end
*IRSIGN
c
c     subroutine to determine position of point site0 with respect
c     to the plane spanned by points site1, site2, site3
c
      subroutine irsign(xi, yi, zi, x, y, z, x2, y2, z2, site0, site1,
     *                  site2, site3, mhalf, mfull, isclp, epz, ipout)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*)
      double precision epz, dist
      integer isclp(*), mhalf, mfull, ipossi
      integer site0, site1, site2, site3, ipout
c
      call dstnce(xi, yi, zi, site1, site2, site3, epz, site0, dist,
     *            ipossi)
      if(ipossi.eq.0) then
         ipout = 1
         if(dist.lt.0.0d0) ipout = -1
      else
         call ipsign(x, y, z, x2, y2, z2, site1, site2, site3, site0,
     *               mhalf, mfull, isclp, ipout)
      endif
c
      return
      end
*FCTEST
c
c     This subroutine will test whether a ray with origin a vertex of
c     a tetrahedron intersects the facet opposite the vertex of the
c     tetrahedron and whether a point in the interior of the ray is
c     contained in the tetrahedron
c
      subroutine fctest(xi, yi, zi, x, y, z, x2, y2, z2, xc, yc, zc,
     *                  itype, k, imist, a, b, c, d, side1, side2,
     *                  mhalf, mfull, isclp, epz)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      double precision epz
      integer isclp(*), itype, k, imist, mhalf, mfull, iasign, ipout
      integer iside(4), a, b, c, d, side1, side2
c
c     determine whether ray with origin point a and through point k 
c     intersects facet of current tetrahedron opposite to point a
c
      if(b.le.8 .or. c.le.8 .or. d.le.8) go to 100
c
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, a, c, d,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(2) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, a, d, b,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(3) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, a, b, c,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(4) = ipout
c
      if(iside(2).lt.0 .or. iside(3).lt.0 .or. iside(4).lt.0) go to 1000
c
c     determine whether point k is in tetrahedron
c
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, b, d, c,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(1) = ipout
      if(iside(1).lt.0) go to 500
c
   50 continue
      call pntype(iside, itype, side1, side2)
      go to 1000
c
c     there is at least one artificial vertex
c
  100 continue
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, a, c, d, k,
     *            mhalf, mfull, isclp, iasign)
      iside(2) = iasign
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, a, d, b, k,
     *            mhalf, mfull, isclp, iasign)
      iside(3) = iasign
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, a, b, c, k,
     *            mhalf, mfull, isclp, iasign)
      iside(4) = iasign
c
      if(iside(2).lt.0 .or. iside(3).lt.0 .or. iside(4).lt.0) go to 1000
c
c     determine whether point k is in tetrahedron
c
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, b, d, c, k,
     *            mhalf, mfull, isclp, iasign)
      iside(1) = iasign
      if(iside(1).lt.0) go to 500
      go to 50
c
c     ray intersects facet but point k is not in tetrahedron
c
  500 continue
c
c     ray intersects interior of facet
c
      if(iside(2).gt.0 .and. iside(3).gt.0 .and. iside(4).gt.0) then
         itype = -2
         go to 1000
      endif
c
      if(iside(2).eq.0 .and. iside(3).eq.0 .and. iside(4).eq.0) stop 931
c
c     ray intersects a vertex of facet
c
      if (iside(2).eq.0 .and. iside(3).eq.0) then
         itype = -3
         imist = d
         go to 1000
      elseif (iside(2).eq.0 .and. iside(4).eq.0) then
         itype = -3
         imist = c
         go to 1000
      elseif (iside(3).eq.0 .and. iside(4).eq.0) then
         itype = -3
         imist = b
         go to 1000
      endif
c
c     ray intersects the interior of an edge of facet
c
      itype = -4
      if (iside(2) .eq. 0) then
         imist = c
      elseif (iside(3) .eq. 0) then
         imist = d
      elseif (iside(4) .eq. 0) then 
         imist = b
      else
         stop 932
      endif
c
 1000 continue
      return
      end
*FCFIND
c
c     This subroutine tests whether a point on a ray that intersects
c     the interior of a facet of a tetrahedron is in the tetrahedron
c     and if not finds other facet of the tetrahedron intersected by
c     the ray
c
      subroutine fcfind(xi, yi, zi, x, y, z, x2, y2, z2, xc, yc, zc,
     *                  itype, ileft, k, ilift, imist, b, c, d,
     *                  side1, side2, mhalf, mfull, isclp, epz)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      double precision epz
      integer isclp(*), mhalf, mfull
      integer itype, ileft, k, ilift, imist, iasign, ipout
      integer idut1, idut2, idut3, idot1, idot2, idot3
      integer iside(4), b, c, d, side1, side2
c
c     determine whether point k is in tetrahedron
c
      if(b.le.8 .or. c.le.8 .or. d.le.8 .or. ilift.le.8) go to 100
c
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ilift, d, c,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(2) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ilift, c, b,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(3) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ilift, b, d,
     *            mhalf, mfull, isclp, epz, ipout)
      iside(4) = ipout
c
      if(iside(2).lt.0 .or. iside(3).lt.0 .or. iside(4).lt.0) go to 200
c
   50 continue
      iside(1) = 1
      call pntype(iside, itype, side1, side2)
      go to 1000
c
c     there is at least one artificial vertex
c
  100 continue
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ilift, d, c, k,
     *            mhalf, mfull, isclp, iasign)
      iside(2) = iasign
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ilift, b, d, k,
     *            mhalf, mfull, isclp, iasign)
      iside(4) = iasign
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ilift, c, b, k,
     *            mhalf, mfull, isclp, iasign)
      iside(3) = iasign
c
      if(iside(2).lt.0 .or. iside(3).lt.0 .or. iside(4).lt.0) go to 200
      go to 50
c
c     k is not in tetrahedron
c
c     determine position of ilift with repect to current situation
c
  200 continue
      if(b.le.8 .or. c.le.8 .or. d.le.8 .or. ilift.le.8) go to 300
c
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, ileft, ilift, c,
     *            b, mhalf, mfull, isclp, epz, idut1)
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, ileft, ilift, d,
     *            c, mhalf, mfull, isclp, epz, idut2)
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, ileft, ilift, b,
     *            d, mhalf, mfull, isclp, epz, idut3)
c
      if(idut1.le.0 .or. idut2.le.0 .or. idut3.le.0) go to 700
      go to 400
c
c     there is at least one artificial vertex
c
  300 continue
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ilift, c, b, ileft,
     *            mhalf, mfull, isclp, idut1)
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ilift, d, c, ileft,
     *            mhalf, mfull, isclp, idut2)
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ilift, b, d, ileft,
     *            mhalf, mfull, isclp, idut3)
c
      if(idut1.le.0 .or. idut2.le.0 .or. idut3.le.0) go to 700
c
c     ilift, ileft, b, d, c, form a strictly convex set
c
  400 continue
      if(b.le.8 .or. c.le.8 .or. d.le.8 .or. ilift.le.8) go to 500
c
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ileft, b,
     *            ilift, mhalf, mfull, isclp, epz, idot1)
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ileft, c,
     *            ilift, mhalf, mfull, isclp, epz, idot2)
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ileft, d,
     *            ilift, mhalf, mfull, isclp, epz, idot3)
      go to 600
c
  500 continue
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ileft, b, ilift, k,
     *            mhalf, mfull, isclp, idot1)
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ileft, c, ilift, k,
     *            mhalf, mfull, isclp, idot2)
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ileft, d, ilift, k,
     *            mhalf, mfull, isclp, idot3)
c
  600 continue
      itype = -2
      if(idot1 .lt. 0 .and. idot2 .gt. 0) then
          imist = d
      elseif(idot2 .lt. 0 .and. idot3 .gt. 0) then
          imist = b
      elseif(idot3 .lt. 0 .and. idot1 .gt. 0) then
          imist = c
      elseif(idot1 .eq. 0 .and. idot2 .eq. 0 .and. idot3 .eq.0) then
          itype = -3
      elseif(idot1 .eq. 0) then
          itype = -4
          imist = b
      elseif(idot2 .eq. 0) then
          itype = -4
          imist = c
      elseif(idot3 .eq. 0) then
          itype = -4
          imist = d
      else
          stop 951
      endif
      go to 1000
c
  700 continue
      if(idut1.le.0 .and. idut2.le.0 .and. idut3.le.0) stop 952
      itype = -2
      if(idut1.le.0 .and. idut2.le.0)then
         imist = c
      elseif(idut2.le.0 .and. idut3.le.0)then
         imist = d
      elseif(idut3.le.0 .and. idut1.le.0)then
         imist = b
      elseif(idut1.le.0)then
         if(d.gt.8 .and. ilift.gt.8) then
            call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ileft, d,
     *                  ilift, mhalf, mfull, isclp, epz, idot3)
         else
            call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ileft, d,
     *                  ilift, k, mhalf, mfull, isclp, idot3)
         endif
         if(idot3.gt.0)then
            imist = b
         elseif(idot3.lt.0)then
            imist = c
         else
            itype = -4
            imist = d
         endif
      elseif(idut2.le.0)then
         if(b.gt.8 .and. ilift.gt.8) then
            call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ileft, b,
     *                  ilift, mhalf, mfull, isclp, epz, idot1)
         else
            call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ileft, b,
     *                  ilift, k, mhalf, mfull, isclp, idot1)
         endif
         if(idot1.gt.0)then
            imist = c
         elseif(idot1.lt.0)then
            imist = d
         else
            itype = -4
            imist = b
         endif
      else
         if(c.gt.8 .and. ilift.gt.8) then
            call irsign(xi, yi, zi, x, y, z, x2, y2, z2, k, ileft, c,
     *                  ilift, mhalf, mfull, isclp, epz, idot2)
         else
            call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ileft, c,
     *                  ilift, k, mhalf, mfull, isclp, idot2)
         endif
         if(idot2.gt.0)then
            imist = d
         elseif(idot2.lt.0)then
            imist = b
         else
            itype = -4
            imist = c
         endif
      endif
c
 1000 continue
      return
      end
*FCEDGE
c
c     This subroutine will test whether a ray through an edge of a
c     tetrahedron intersects either of the facets of the tetrahedron
c     opposite the edge and whether a point in the interior of the
c     ray is contained in the tetrahedron
c
      subroutine fcedge(x, y, z, x2, y2, z2, xc, yc, zc, itype, ileft,
     *                  k, icon, iscur, imist, ivnxt, site0, site1,
     *                  side1, side2, mhalf, mfull, isclp, ITCHK)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer icon(8,*)
      integer isclp(*), mhalf, mfull
      integer itype, ileft, k, iscur, imist, ivnxt
      integer iside(4), site0, site1, site2, site3, side1, side2
      integer isnow, idut, iasign, idot0, itchk, ipout
c
c     find intersecting facet
c
      call reordr(icon, site0, site1, iscur)
      site2 = icon(7,iscur)
      site0 = icon(8,iscur)
      isnow = iabs(icon(1,iscur))
c
  300 continue
      ITCHK = ITCHK+1
      if(isnow.le.0.or.isnow.gt.ivnxt) stop 961
      if(isnow.eq.iscur) stop 963
      call reordr(icon, site0, site1, isnow)
      site3 = icon(8,isnow)
      if(site1.gt.8 .and. site2.gt.8 .and. site3.gt.8) then
         call ipsign(x, y, z, x2, y2, z2, site1, site3, site2, k,
     *               mhalf, mfull, isclp, idut)
      else
         call artsig(x, y, z, x2, y2, z2, xc, yc, zc, site1, site3,
     *               site2, k, mhalf, mfull, isclp, idut)
      endif
      if(idut.ge.0) go to 400
      site0 = site3
      isnow = iabs(icon(1,isnow))
      go to 300
c
  400 continue
      iscur = isnow
c
c     determine whether point k is in tetrahedron
c
      if(site0.le.8 .or. site1.le.8 .or. site2.le.8 .or. site3.le.8)
     *   go to 500
c
      call ipsign(x, y, z, x2, y2, z2, site1, site0, site3, k,
     *            mhalf, mfull, isclp, ipout)
      iside(3) = ipout
      call ipsign(x, y, z, x2, y2, z2, site2, site3, site0, k,
     *            mhalf, mfull, isclp, ipout)
      iside(2) = ipout
c
      if(iside(2).lt.0 .or. iside(3).lt.0) go to 600
c
  450 continue
      iside(1) = idut
      iside(4) = 1
      call pntype(iside, itype, side1, side2)
      go to 1000
c
c     there is at least one artificial vertex
c
  500 continue
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, site1, site0, site3,
     *            k, mhalf, mfull, isclp, iasign)
      iside(3) = iasign
      call artsig(x, y, z, x2, y2, z2, xc, yc, zc, site2, site3, site0,
     *            k, mhalf, mfull, isclp, iasign)
      iside(2) = iasign
c
      if(iside(2).lt.0 .or. iside(3).lt.0) go to 600
      go to 450
c
c     k is not in tetrahedron but ray intersects one of the facets
c     of the tetrahedron opposite the edge
c
  600 continue
      if(site0.gt.8 .and. site3.gt.8) then
         call ipsign(x, y, z, x2, y2, z2, ileft, site0, site3, k,
     *               mhalf, mfull, isclp, idot0)
      else
         call artsig(x, y, z, x2, y2, z2, xc, yc, zc, ileft, site0,
     *               site3, k, mhalf, mfull, isclp, idot0)
      endif
      if(idot0.gt.0)then
         if(idut.gt.0) then
            itype = -2
            imist = site1
         else
            itype = -4
            imist = site2
         endif
      elseif(idot0.lt.0)then
         if(idut.gt.0) then
            itype = -2
            imist = site2
         else
            itype = -4
            imist = site1
         endif
      else
         if(idut.gt.0) then
            itype = -4
            imist = site0
         else
            itype = -3
            imist = site3
         endif
      endif
c
 1000 continue
      return
      end
*PNTYPE
c
c     This subroutine determines point type with respect to a
c     tetrahedron that contains a point
c
      subroutine pntype(iside, itype, side1, side2)
c
      integer iside(*), itype, side1, side2
c
c     point is in the interior of tetrahedron
c
      if(iside(1).gt.0 .and. iside(2).gt.0 .and. iside(3).gt.0 .and.
     *   iside(4).gt.0) then
         itype = 2
         go to 1000
      endif
c
c     unacceptable situation
c
      if(iside(1).eq.0 .and. iside(2).eq.0 .and. iside(3).eq.0 .and.
     *   iside(4).eq.0) stop 971
c
c     point is a vertex of tetrahedron
c
      if(iside(1).eq.0 .and. iside(2).eq.0 .and. iside(3).eq.0) then
         itype = 1
         side1 = 4
         go to 1000
      elseif(iside(1).eq.0 .and. iside(2).eq.0 .and. iside(4).eq.0) then
         itype = 1
         side1 = 3
         go to 1000
      elseif(iside(1).eq.0 .and. iside(3).eq.0 .and. iside(4).eq.0) then
         itype = 1
         side1 = 2
         go to 1000
      elseif(iside(2).eq.0 .and. iside(3).eq.0 .and. iside(4).eq.0) then
         itype = 1
         side1 = 1
         go to 1000
      endif
c
c     point is in the interior of an edge of tetrahedron
c
      if (iside(1).eq.0 .and. iside(2).eq.0) then
         itype = 3 
         side1 = 1
         side2 = 2
         go to 1000
      elseif (iside(1).eq.0 .and. iside(3).eq.0) then
         itype = 3 
         side1 = 1
         side2 = 3
         go to 1000
      elseif (iside(1).eq.0 .and. iside(4).eq.0) then
         itype = 3 
         side1 = 1
         side2 = 4
         go to 1000
      elseif (iside(2).eq.0 .and. iside(3).eq.0) then
         itype = 3 
         side1 = 2
         side2 = 3
         go to 1000
      elseif (iside(2).eq.0 .and. iside(4).eq.0) then
         itype = 3 
         side1 = 2
         side2 = 4
         go to 1000
      elseif (iside(3).eq.0 .and. iside(4).eq.0) then
         itype = 3 
         side1 = 3
         side2 = 4
         go to 1000
      endif
c
c     point is in the interior of a facet of tetrahedron
c
      itype = 4
      if (iside(1) .eq. 0) then
         side1 = 1
      elseif (iside(2) .eq. 0) then
         side1 = 2
      elseif (iside(3) .eq. 0) then
         side1 = 3
      elseif (iside(4) .eq. 0) then 
         side1 = 4
      else
         stop 972
      endif
c
 1000 continue
      return
      end
*ARTSIG
c
c     This subroutine determines the position of a point with respect
c     to a plane when artificial points may be involved
c
      subroutine artsig(x, y, z, x2, y2, z2, xc, yc, zc, ib, id, ic, k,
     *                  mhalf, mfull, isclp, iasign)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer isclp(*), mhalf, mfull
      integer b, d, c, ib, id, ic, k, iasign, ifn
c
      call vrtarr(ib, id, ic, b, d, c)
c
      if(d.gt.8 .and. c.gt.8) then
         call ipsign(x, y, z, x2, y2, z2, b, d, c, k, mhalf,
     *               mfull, isclp, iasign)
         go to 100
      elseif(b.le.8) then
         iasign = 1
         go to 100
      elseif(d.le.8.and.c.le.8) then
         call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc, d, c, k, b,
     *               mhalf, mfull, isclp, iasign)
         if(iasign.ne.0) go to 100
         call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc, d, c, k, b,
     *               mhalf, mfull, isclp, iasign)
         go to 100
      elseif(d.le.8) then
         ifn = 0
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc, b, c, k, d,
     *               mhalf, mfull, isclp, ifn, iasign)
      else
         ifn = 1
         call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc, b, d, k, c,
     *               mhalf, mfull, isclp, ifn, iasign)
      endif
      if(iasign.ne.0) go to 100
      call ipsign(x, y, z, x2, y2, z2, b, d, c, k, mhalf,
     *            mfull, isclp, iasign)
c
  100 continue
      return
      end
*VRTINS
c
c This subroutine will insert a point located at a vertex of the current
c tetrahedron
c
      subroutine vrtins (k, w, w2, icon, nvmax, tetra, curr, is, id,
     *                   iftal, side, ifl, newtts, ired, delaun,
     *                   flphis, mhalf, mfull, isclw)
c
      integer w(*), w2(*), icon(8,*), is(*), ifl(*), id(*)
      integer k, nvmax, iftal, newtts, ired, mhalf, mfull, itide
      integer tetra, curr, side, site1, site2, isclw(*)
      logical delaun, flphis
      integer isini, isfrt, islst, iscur, indx, i, j, isadj, ilift
c
      site1 = icon(side+4,curr)
      if(site1.le.8) stop 1010
c
      if(delaun) go to 50
      call iwsign(w, w2, k, site1, mhalf, mfull, isclw, itide)
      if(itide .gt. 0) go to 100
   50 continue
      is(k) = -site1
      ired = 1
      go to 2000
c
  100 continue
      if(.not.flphis) go to 1000
      ired = 0
      is(site1) = -2
      call sitord(icon, site1, curr)
      isini = curr
      tetra = tetra + 1
      if(tetra .gt. nvmax) stop 1020
      isfrt = tetra
      icon(5,isini) = -tetra
      icon(5,tetra) = isini
      islst = isini
      iscur = icon(2,islst)
      if(iscur.le.0) stop 1030
      site2 = icon(8,islst)
      is(k) = tetra
      is(site2) = tetra
      is(icon(6,islst)) = tetra
      is(icon(7,islst)) = tetra
      newtts = 1
      ifl(1) = tetra
c
  200 continue
      call reordr(icon, site1, site2, iscur)
      if(icon(4,iscur).ne.islst) stop 1040
      tetra = tetra + 1
      if(tetra .gt. nvmax) stop 1050
      icon(5,iscur) = -tetra
      icon(5,tetra) = iscur
      islst = iscur
      indx = 2
      iscur = icon(2,islst)
      if(iscur.le.0) stop 1060
      site2 = icon(8,islst)
      is(site2) = tetra
      newtts = newtts + 1
      ifl(newtts) = tetra
      if(icon(5,iscur).gt.0) go to 200
c
  500 continue
      if(indx .eq. 2) then
         indx = 3
         iscur = icon(3,islst)
         if(iscur.le.0) stop 1070
         site2 = icon(6,islst)
         if(icon(5,iscur).gt.0) go to 200
         go to 500
      elseif(indx .eq. 3) then
         if(islst .ne. isini) then
            iscur = islst
            islst = icon(4,iscur)
            if(icon(2,islst) .eq. iscur) then
               indx = 2
            elseif(icon(3,islst) .eq. iscur) then
               indx = 3
            elseif(icon(4,islst) .eq. iscur) then
               indx = 4
            else
               stop 1080
            endif
            go to 500
         else
            indx = 4
            iscur = icon(4,islst)
            if(iscur.le.0) stop 1090
            site2 = icon(7,islst)
            if(icon(5,iscur).gt.0) go to 200
            go to 500
         endif
      endif
      if(islst .ne. isini) stop 1110
c
      do 800 i = isfrt, tetra
         iscur = icon(5,i)
         do 600 j = 2, 4
            icon(j,i)=-icon(5,icon(j,iscur))
  600    continue
         do 700 j = 6, 8
            icon(j,i) = icon(j,iscur)
  700    continue
  800 continue
c
      do 900 i = isfrt, tetra
         iscur = icon(5,i)
         icon(5,i) = k
         icon(5,iscur) = -site1
         isadj = icon(1,iscur)
         icon(1,iscur) = i
         icon(2,iscur) = 0
         icon(3,iscur) = 0
         icon(4,iscur) = 0
         icon(1,i) = isadj
         if(isadj .eq. 0) go to 900
         do 840 j = 1, 4
            if(icon(j,isadj) .eq. iscur) go to 860
  840    continue
         stop 1120
  860    continue
         icon(j,isadj) = i
  900 continue
      go to 2000
c
 1000 continue
      ired = 0
      iftal = iftal + 1
      iscur = is(site1)
      if(iscur.le.0.or.iscur.gt.tetra) stop 1130
      isini = iscur
      call sitord(icon, site1, isini)
      is(site1) = -2
      is(k) = isini
      newtts = 0
c
 1400 continue
      newtts = newtts+1
      ifl(newtts) = iscur
      id(iscur) = iftal
      icon(5,iscur) = k
c
      isadj = icon(2,iscur)
      if(isadj.le.0.or.isadj.gt.tetra) stop 1140
      if(id(isadj).eq.iftal) go to 1600
      ilift = icon(8,iscur)
      go to 1900
 1600 continue
      isadj = icon(3,iscur)
      if(isadj.le.0.or.isadj.gt.tetra) stop 1150
      if(id(isadj).eq.iftal) go to 1700
      ilift = icon(6,iscur)
      go to 1900
 1700 continue
      isadj = icon(4,iscur)
      if(isadj.le.0.or.isadj.gt.tetra) stop 1160
      if(iscur .eq. isini) go to 1800
      if(icon(3,isadj) .eq. iscur) then
          iscur = isadj
          go to 1700
      elseif(icon(2,isadj) .eq. iscur) then
          iscur = isadj
          go to 1600
      elseif(icon(4,isadj) .eq. iscur) then
          if(isadj .ne. isini) stop 1170
          go to 2000
      else
          stop 1180
      endif
 1800 continue
      if(id(isadj).eq.iftal) go to 2000
      ilift = icon(7,iscur)
c
 1900 continue
      call reordr(icon, site1, ilift, isadj)
      iscur = isadj
      go to 1400
c
 2000 continue
      return
      end
*INTINS
c
c This subroutine will insert a point located in the interior of the
c current tetrahedron.  Four new tetra will be created.
c
      subroutine intins(k, xi, yi, zi, wi, x, y, z, w, x2, y2, z2, w2,
     *                  icon, ih, ihn, nhmax, nvmax, tetra, curr, is,
     *                  ifl, newtts, ired, delaun, flphis, mhalf,
     *                  mfull, isclp, isclw, isclr, epz)
c
      double precision xi(*), yi(*), zi(*), wi(*)
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer icon(8, *), is(*), ifl(*), ih(*)
      double precision epz, tdist, xctr, yctr, zctr
      integer ihn, nhmax, nvmax, newtts, ired, mhalf, mfull
      integer isclp(*), isclw(*), isclr(*)
      integer tetra, curr, adj, a, b, c, d, newi(4)
      logical delaun, flphis
      integer k, ipossi, itide, new1, new2, new3, new4, i, j
c
      a = icon(5,curr)
      b = icon(6,curr)
      c = icon(7,curr)
      d = icon(8,curr)
c
      if(delaun) go to 30
      if(a.le.8 .or. b.le.8 .or. c.le.8 .or. d.le.8) go to 30
      call ctrad(xi, yi, zi, wi, xctr, yctr, zctr, a, b, c, d,
     *           epz, delaun, ipossi)
      if(ipossi.eq.1) go to 20
      call bisphr(xi, yi, zi, wi, k, a, epz, xctr, yctr, zctr, tdist,
     *            delaun, ipossi)
      if(ipossi.eq.1) go to 20
      itide = 1
      if(tdist.gt.0.0d0) itide = -1
      go to 25
   20 continue
      call iqsign(x, y, z, w, x2, y2, z2, w2, a, b, c, d, k,
     *            mhalf, mfull, isclp, isclw, isclr, delaun, itide)
   25 continue
      if(itide.le.0) go to 30
      is(k) = -3
      ired = 1
      go to 1000
c
   30 continue
      ired = 0
c
      if(flphis) then
         new1 = tetra + 1
         new2 = tetra + 2
         new3 = tetra + 3
         new4 = tetra + 4
         tetra = new4
         if (tetra .gt. nvmax) stop 1210
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            tetra = new1
            if(tetra .gt. nvmax) stop 1220
         else
            new1 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new2 = tetra+1
            tetra = new2
            if(tetra .gt. nvmax) stop 1230
         else
            new2 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new3 = tetra+1
            tetra = new3
            if(tetra .gt. nvmax) stop 1240
         else
            new3 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new4 = tetra+1
            tetra = new4
            if(tetra .gt. nvmax) stop 1250
         else
            new4 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      ifl(1) = new1
      ifl(2) = new2
      ifl(3) = new3
      ifl(4) = new4
c
      newtts = 4
c
      do 50 i = 1, 8
         icon(i,new1) = icon(i,curr)
         icon(i,new2) = icon(i,curr)
         icon(i,new3) = icon(i,curr)
         icon(i,new4) = icon(i,curr)
   50 continue
c
c     update new1
c
      icon(2,new1) = new2
      icon(3,new1) = new3
      icon(4,new1) = new4
      icon(5,new1) = k
c
c     update new2
c
      icon(1,new2) = new1
      icon(3,new2) = new3
      icon(4,new2) = new4
      icon(6,new2) = k
      call sitord (icon, k, new2)
c
c     update new3
c
      icon(1,new3) = new1
      icon(2,new3) = new2
      icon(4,new3) = new4
      icon(7,new3) = k
      call sitord (icon, k, new3)
c
c     update new4
c
      icon(1,new4) = new1
      icon(2,new4) = new2
      icon(3,new4) = new3
      icon(8,new4) = k
      call sitord (icon, k, new4)
c      
c     update is(*)
c
      is(k) = new1
c
      is(a) = new2
      is(b) = new1
      is(c) = new1
      is(d) = new1
c
c     update neighboring tetra
c
      newi(1) = new1
      newi(2) = new2
      newi(3) = new3
      newi(4) = new4
      do 200 i = 1, 4
         adj = icon(i,curr)
         if (adj .ne. 0) then
            do 175 j = 1, 4
               if (icon(j,adj) .eq. curr) then
                  icon(j,adj) = newi(i)
                  goto 200
               endif
  175       continue
            stop 1260
         endif
  200 continue
c
c     flag current tetra to denote that it has children
c
      if(flphis) then
         icon(5,curr) = - icon(5,curr)
         icon(1,curr) = new1  
         icon(2,curr) = new2
         icon(3,curr) = new3
         icon(4,curr) = new4
      else
         icon(5,curr) = - icon(5,curr)
         ihn = ihn+1
         if(ihn.gt.nhmax) stop 1270
         ih(ihn) = curr
      endif
c
 1000 continue
      return
      end
*SIDINS
c
c This subroutine will insert a point which is on side
c of curr (and adj). Six new tetra will be created.
c 
      subroutine sidins(k, x, y, z, w, x2, y2, z2, w2, icon, ih, ihn,
     *                   nhmax, nvmax, tetra, curr, is, side, ifl,
     *                   newtts, ired, delaun, flphis, mhalf, mfull,
     *                   isclp, isclw, isclr)
c
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer icon(8, *), is(*), ifl(*), ih(*)
      integer tetra, curr, adj, side, site1, site2, site3, site4, temp
      integer ihn, nhmax, nvmax, newtts, ired, mhalf, mfull
      integer isclp(*), isclw(*), isclr(*), newi(6)
      logical delaun, flphis
      integer k, itide, new1, new2, new3, new4, new5, new6, i, j
c
      adj = icon(side, curr)
      if (adj .eq. 0) stop 1310
c
c     rearrange curr
c
      site1 = icon(side+4,curr)
      call sitord(icon, site1, curr)
c
      site2 = icon(6,curr)
      site3 = icon(7,curr)
      site4 = icon(8,curr)
c
      if(delaun) go to 30
      if(site2.le.8 .or. site3.le.8 .or. site4.le.8) go to 30
      call iqsig2(x, y, z, w, x2, y2, z2, w2, site2, site3, site4,
     *            k, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
      if(itide .le. 0) go to 30
      is(k) = -3
      ired = 1
      go to 1000
c
   30 continue
      ired = 0
c
      if(flphis) then
         new1 = tetra + 1
         new2 = tetra + 2
         new3 = tetra + 3
         new4 = tetra + 4
         new5 = tetra + 5
         new6 = tetra + 6
         tetra = new6
         if(tetra .gt. nvmax) stop 1320
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            tetra = new1
            if(tetra .gt. nvmax) stop 1330
         else
            new1 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new2 = tetra+1
            tetra = new2
            if(tetra .gt. nvmax) stop 1340
         else
            new2 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new3 = tetra+1
            tetra = new3
            if(tetra .gt. nvmax) stop 1350
         else
            new3 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new4 = tetra+1
            tetra = new4
            if(tetra .gt. nvmax) stop 1360
         else
            new4 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new5 = tetra+1
            tetra = new5
            if(tetra .gt. nvmax) stop 1370
         else
            new5 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new6 = tetra+1
            tetra = new6
            if(tetra .gt. nvmax) stop 1380
         else
            new6 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      ifl(1) = new1
      ifl(2) = new2
      ifl(3) = new3
      ifl(4) = new4
      ifl(5) = new5
      ifl(6) = new6
c
      newtts = 6
c
c     create new1, new2, new3
c
      do 50 i = 1, 8
         icon(i,new1) = icon(i,curr)
         icon(i,new2) = icon(i,curr)
         icon(i,new3) = icon(i,curr)
   50 continue  
c
      icon(1,new1) = new4
      icon(3,new1) = new2
      icon(4,new1) = new3
      icon(6,new1) = k
      call sitord (icon, k, new1)
c
      icon(1,new2) = new6
      icon(2,new2) = new1
      icon(4,new2) = new3
      icon(7,new2) = k
      call sitord (icon, k, new2)
c
      icon(1,new3) = new5
      icon(2,new3) = new1
      icon(3,new3) = new2
      icon(8,new3) = k
      call sitord (icon, k, new3)
c
c     update is(*)
c
      is(k) = new1
c
      is(site1) = new1
      is(site2) = new2
      is(site3) = new1
      is(site4) = new1
c
c     update curr's neighbors
c
      newi(1) = new1
      newi(2) = new2
      newi(3) = new3
      do 200 i = 2, 4
         temp = icon(i,curr)
         if (temp .ne. 0) then
            do 150 j = 1, 4
               if (icon(j,temp) .eq. curr) then
                  icon(j,temp) = newi(i-1)
                  goto 200
               endif
  150       continue
            stop 1390
         endif
  200 continue
c
c     flag curr to show its children
c
      if(flphis) then
         icon(5,curr) = - icon(5,curr)
         icon(1,curr) = new1 
         icon(2,curr) = new2
         icon(3,curr) = new3
         icon(4,curr) = -adj
      else
         icon(5,curr) = - icon(5,curr)
         ihn = ihn+1
         if(ihn.gt.nhmax) stop 1410
         ih(ihn) = curr
      endif
c
c     update 2nd tet (adj)
c
      do 300 i = 1, 4
         if (icon(i,adj) .eq. curr) then
            site1 = icon(i+4,adj)
            goto 325
         endif
  300 continue
      stop 1420
  325 continue
c
      call reordr (icon, site1, site2, adj)
c
c     create new4, new5, new6
c
      do 350 i = 1, 8
         icon(i,new4) = icon(i,adj)
         icon(i,new5) = icon(i,adj)
         icon(i,new6) = icon(i,adj)
  350 continue  
c
      icon(1,new4) = new1
      icon(3,new4) = new5
      icon(4,new4) = new6
      icon(6,new4) = k
      call sitord (icon, k, new4)
c
      icon(1,new5) = new3
      icon(2,new5) = new4
      icon(4,new5) = new6
      icon(7,new5) = k
      call sitord (icon, k, new5)
c
      icon(1,new6) = new2
      icon(2,new6) = new4
      icon(3,new6) = new5
      icon(8,new6) = k
      call sitord (icon, k, new6)
c
c     update is(*)
c
      is(site1) = new4
c
c     update adj's neighbors
c
      newi(4) = new4
      newi(5) = new5
      newi(6) = new6
      do 500 i = 2, 4
         temp = icon(i,adj)
         if (temp .ne. 0) then
            do 450 j = 1, 4
               if (icon(j,temp) .eq. adj) then
                  icon(j,temp) = newi(i+2) 
                  goto 500
               endif
  450       continue
            stop 1430
         endif
  500 continue
c
c     flag adj to show its children
c
      if(flphis) then
         icon(5,adj) = - icon(5,adj)
         icon(1,adj) = new4 
         icon(2,adj) = new5
         icon(3,adj) = new6
         icon(4,adj) = -curr
      else
         icon(5,adj) = - icon(5,adj)
         ihn = ihn+1
         if(ihn.gt.nhmax) stop 1440
         ih(ihn) = adj
      endif
c
 1000 continue
      return
      end
*EDGINS
c
c This subroutine will insert point on edge of curr tetra.
c
      subroutine edgins(k, x, y, z, w, x2, y2, z2, w2, icon, ih, ihn,
     *                  nhmax, nvmax, tetra, curr, is, side1, side2,
     *                  ifl, newtts, ired, delaun, flphis,
     *                  mhalf, mfull, isclp, isclw, isclr)
c
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer icon(8,*), is(*), ifl(*), ih(*)
      integer ihn, nhmax, nvmax, newtts, ired
      integer isclp(*), isclw(*), isclr(*), mhalf, mfull
      integer tetra, curr, side1, side2, fndsit
      integer site1, site2, prvtet, adj, firtet
      logical delaun, flphis
      integer k, i, itide, new1, new2, nel1, nel2, now1, now2, j
c
c     find endpoints of edge
c
      fndsit = 0
      do 100 i = 5, 8
         if (fndsit .eq. 0) then
            if (i .eq. (side1+4) .or. i .eq. (side2+4)) then
               goto 100
            else
               site1 = icon(i,curr)
               fndsit = 1
            endif
         else 
            if (i .eq. (side1+4) .or. i .eq. (side2+4)) then
               goto 100
            else
               site2 = icon(i,curr)
               goto 150
            endif
         endif
  100 continue
      stop 1510
  150 continue
c
      if(delaun) go to 160
      if(site1.le.8 .or. site2.le.8) go to 160
      call iqsig1(x, y, z, w, x2, y2, z2, w2, site1, site2, k,
     *            mhalf, mfull, isclp, isclw, isclr, delaun, itide)
      if(itide .le. 0) go to 160
      is(k) = -3
      ired = 1
      go to 1000
c
c     order vertices of tetrahedra around edge
c
  160 continue
      firtet = curr
c
  163 continue
      call reordr(icon, site1, site2, curr)
      curr = icon(3,curr)
      if(curr.ne.firtet) go to 163
c
      ired = 0
      if(flphis) then
         new1 = tetra + 1
         new2 = tetra + 2
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            new2 = tetra+2
         elseif(ihn.eq.1) then
            new1 = ih(1)
            new2 = tetra+1
         else
            new1 = ih(ihn)
            new2 = ih(ihn-1)
         endif
      endif
      nel1 = new1
      nel2 = new2
      newtts = 0
c
      is(k) = new1
      is(site1) = new2
      is(site2) = new1
c
c     create 2 new tetra  
c
  175 continue
      if(flphis) then
         now1 = tetra + 1
         now2 = tetra + 2
         tetra = now2
         if(tetra .gt. nvmax) stop 1520
      else
         if(ihn.eq.0) then
            now1 = tetra+1
            tetra = now1
            if(tetra .gt. nvmax) stop 1530
         else
            now1 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            now2 = tetra+1
            tetra = now2
            if(tetra .gt. nvmax) stop 1540
         else
            now2 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      newtts = newtts + 1
      ifl(newtts) = now1 
      newtts = newtts + 1
      ifl(newtts) = now2 
c
      do 180 i = 1, 8
         icon(i,now1) = icon(i,curr)
         icon(i,now2) = icon(i,curr)
  180 continue
c
      icon(2,now1) = now2
      icon(5,now1) = k
      icon(1,now2) = now1
      icon(6,now2) = k
c
      icon(3,nel1) = now1
      icon(3,nel2) = now2
      icon(4,now1) = nel1
      icon(4,now2) = nel2
      call sitord (icon, k, now1)
      call sitord (icon, k, now2)
c
c     update is(*)
c
      is(icon(7,curr)) = now1
c
c     update neighbors of curr
c
      do 300 i = 1, 2
         adj = icon(i,curr)
         if (adj .ne. 0) then
            do 250 j = 1, 4
               if (icon(j,adj) .eq. curr) then
                  if (i .eq. 1) then
                     icon(j,adj) = now1
                  else
                     icon(j,adj) = now2
                  endif
                  goto 300
               endif
  250       continue
            stop 1550
         endif
  300 continue
c
      prvtet = curr
      curr = icon(3,curr)
c
c     show children of old tetra
c
      if(flphis) then
         icon(5,prvtet) = - icon(5,prvtet)
         icon(1,prvtet) = now1
         icon(2,prvtet) = now2
         icon(3,prvtet) = -curr
         icon(4,prvtet) = -icon(4,prvtet)
      else
         icon(5,prvtet) = - icon(5,prvtet)
         ihn = ihn+1
         if(ihn.gt.nhmax) stop 1560
         ih(ihn) = prvtet
      endif
c
c     go to next tetrahedron until we're back at firtet
c
      if (curr .ne. firtet) then
         nel1 = now1
         nel2 = now2
         go to 175
      else
         icon(4,new1) = now1
         icon(2,new2) = now2
         icon(3,now1) = new1
         icon(3,now2) = new2
      endif
c
 1000 continue
      return
      end
*SPHERE
c
c This subroutine will optimize locally at point k for the
c Regular/Delaunay property
c
      subroutine sphere(k, icon, ifl, ih, ihn, nhmax, newtts, xi, yi,
     *                  zi, wi, x, y, z, w, x2, y2, z2, w2, tetra, is,
     *                  nvmax, xc, yc, zc, delaun, flphis, mhalf, mfull,
     *                  isclp, isclw, isclr, epz)
c
      double precision xi(*), yi(*), zi(*), wi(*)
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer icon(8,*), ifl(*), is(*), ih(*)
      double precision epz, tdist, xctr, yctr, zctr
      integer ihn, nhmax, newtts, nvmax, mhalf, mfull
      integer xc(*), yc(*), zc(*)
      integer adj, opvert, tetra, a, b, c, d, oddsid
      integer isclp(*), isclw(*), isclr(*), iside(4), sidist(4)
      logical delaun, flphis
      integer k, i, now, ipossi, itide, j, ib, ic, id, ifn, istt
      integer isodd, isite, ipout
c
      i = 1
  100 continue
      if (i .gt. newtts) goto 1500
c
      now = ifl(i)
      if (icon(5,now) .lt. 0) goto 1000
      a = icon(5,now)
      if(a.ne.k) stop 1610
c
c     look at adj tet 
c
      adj = icon(1,now)
      if (adj .eq. 0) goto 1000
      if(icon(5,adj).eq.k) stop 1620
c
      b = icon(6,now)
      c = icon(7,now)
      d = icon(8,now)
      if(b.lt.c) b=c
      if(b.lt.d) b=d
      if(b.le.8) stop 1630
      call reordr(icon,a,b,now)
      c=icon(7,now)
      d=icon(8,now)
c
c     reorder adj
c
      call reordr (icon, b, c, adj)
c
      if (icon(7,adj) .ne. d) stop 1640
      if (icon(4,adj) .ne. now) stop 1650
c
      opvert = icon(8,adj)
      if(opvert.le.8.or.c.le.8.or.d.le.8) go to 300
c
c     test whether now and adj form a Regular configuration
c
      call ctrad(xi, yi, zi, wi, xctr, yctr, zctr, b, c, d, opvert,
     *           epz, delaun, ipossi)
      if(ipossi.eq.1) go to 150
      call bisphr(xi, yi, zi, wi, a, b, epz, xctr, yctr, zctr, tdist,
     *            delaun, ipossi)
      if(ipossi.eq.1) go to 150
      if(tdist.le.0.0d0) go to 1000
      go to 170
  150 continue
      call iqsign(x, y, z, w, x2, y2, z2, w2, a, b, c, d, opvert,
     *            mhalf, mfull, isclp, isclw, isclr, delaun, itide)
      if(itide .ge. 0) go to 1000
c
c     compute distances from opvert to facets of now
c
  170 continue
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, a, c, opvert,
     *            d, mhalf, mfull, isclp, epz, ipout)
      iside(2) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, a, b, d,
     *            opvert, mhalf, mfull, isclp, epz, ipout)
      iside(3) = ipout
      call irsign(xi, yi, zi, x, y, z, x2, y2, z2, a, b, opvert,
     *            c, mhalf, mfull, isclp, epz, ipout)
      iside(4) = ipout
c
c     set sidist array
c     
      do 200 j = 2, 4
         if(iside(j) .gt. 0) then
            sidist(j) = 0
         elseif(iside(j) .lt. 0) then
            sidist(j) = -1
         else
            sidist(j) = 1
         endif
  200 continue
c
c     flip according to type of flip
c
      if ((sidist(2). eq. 0) .and. (sidist(3). eq .0) .and.
     *(sidist(4). eq. 0)) then
         call flip23 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, nvmax)
         go to 1000
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.0)) then 
         oddsid = 2
         call flip32 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.-1) .and.
     *(sidist(4).eq.0)) then
         oddsid = 3
         call flip32 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.-1)) then
         oddsid = 4
         call flip32 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.0)) then
         oddsid = 2
         call flip22 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.1) .and.
     *(sidist(4).eq.0)) then
         oddsid = 3
         call flip22 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.1)) then
         oddsid = 4
         call flip22 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif(delaun) then
         write(*,*)'Warning: Delaunay sphere violation'
         go to 1000
      endif
c
      if ((sidist(2).eq.0) .and. (sidist(3).eq.-1) .and.
     *(sidist(4).eq.-1)) then 
         oddsid = 2
         call flip41 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.-1)) then
         oddsid = 3
         call flip41 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.-1) .and.
     *(sidist(4).eq.0)) then
         oddsid = 4
         call flip41 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.1) .and.
     *(sidist(4).eq.1)) then
         oddsid = 2
         call flip21 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.1)) then
         oddsid = 3
         call flip21 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.1) .and.
     *(sidist(4).eq.0)) then
         oddsid = 4
         call flip21 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.-1) .and.
     *(sidist(4).eq.1)) then 
         oddsid = 2
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.-1)) then
         oddsid = 3
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.1) .and.
     *(sidist(4).eq.0)) then
         oddsid = 4
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.1) .and.
     *(sidist(4).eq.-1)) then
         oddsid = -2
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.0) .and.
     *(sidist(4).eq.1)) then
         oddsid = -3
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.-1) .and.
     *(sidist(4).eq.0)) then
         oddsid = -4
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      else
         write(*,*)'Warning: Regular sphere violation'
      endif
      go to 1000
c 
  300 continue
      if(opvert.le.8.and.c.gt.8.and.d.gt.8) go to 1000
c
c     determine signs of distances from opvert to facets of
c     tetrahedron now
c
      if(opvert.gt.8)then
         if(c.le.8.and.d.le.8)then
            call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  d, c, a, opvert, mhalf, mfull, isclp, ipout)
            iside(2) = ipout
            if(iside(2).ne.0) go to 310
            call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  d, c, a, opvert, mhalf, mfull, isclp, ipout)
            iside(2) = ipout
            go to 310
         endif
         call vrtarr(c,d,opvert,ib,ic,id)
         if(ic.gt.8)then
            ifn = 0
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  ib, ic, a, id, mhalf, mfull, isclp, ifn, ipout)
            iside(2) = ipout
         else
            ifn = 1
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  ib, id, a, ic, mhalf, mfull, isclp, ifn, ipout)
            iside(2) = ipout
         endif
         if(iside(2).ne.0) go to 310
         call ipsign(x, y, z, x2, y2, z2, ib, id, ic, a, mhalf,
     *               mfull, isclp, ipout)
         iside(2) = ipout
  310    continue
c
         call vrtarr(b,opvert,d,ib,ic,id)
         if(d.gt.8)then
            call ipsign(x, y, z, x2, y2, z2, ib, id, ic, a, mhalf,
     *                  mfull, isclp, ipout)
            iside(3) = ipout
            go to 320
         endif
         if(ic.gt.8)then
            ifn = 0
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  ib, ic, a, id, mhalf, mfull, isclp, ifn, ipout)
            iside(3) = ipout
         else
            ifn = 1
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  ib, id, a, ic, mhalf, mfull, isclp, ifn, ipout)
            iside(3) = ipout
         endif
         if(iside(3).ne.0) go to 320
         call ipsign(x, y, z, x2, y2, z2, ib, id, ic, a, mhalf,
     *               mfull, isclp, ipout)
         iside(3) = ipout
  320    continue
c
         call vrtarr(b,c,opvert,ib,ic,id)
         if(c.gt.8)then
            call ipsign(x, y, z, x2, y2, z2, ib, id, ic, a, mhalf,
     *                  mfull, isclp, ipout)
            iside(4) = ipout
            go to 330
         endif
         if(ic.gt.8)then
            ifn = 0
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  ib, ic, a, id, mhalf, mfull, isclp, ifn, ipout)
            iside(4) = ipout
         else
            ifn = 1
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  ib, id, a, ic, mhalf, mfull, isclp, ifn, ipout)
            iside(4) = ipout
         endif
         if(iside(4).ne.0) go to 330
         call ipsign(x, y, z, x2, y2, z2, ib, id, ic, a, mhalf,
     *               mfull, isclp, ipout)
         iside(4) = ipout
  330    continue
c
      else
         if(c.le.8.and.d.le.8) then
            iside(2) = 1
         elseif(c.gt.8) then
            call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  opvert, d, a, c, mhalf, mfull, isclp, ipout)
            iside(2) = ipout
            if(iside(2).ne.0) go to 340
            call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  opvert, d, a, c, mhalf, mfull, isclp, ipout)
            iside(2) = ipout
         else
            call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  c, opvert, a, d, mhalf, mfull, isclp, ipout)
            iside(2) = ipout
            if(iside(2).ne.0) go to 340
            call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  c, opvert, a, d, mhalf, mfull, isclp, ipout)
            iside(2) = ipout
         endif
  340    continue
c
         if(d.gt.8)then
            ifn = 1
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                 b, d, a, opvert, mhalf, mfull, isclp, ifn, ipout)
            iside(3) = ipout
         else
            call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  d, opvert, a, b, mhalf, mfull, isclp, ipout)
            iside(3) = ipout
            if(iside(3).ne.0) go to 350
            call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  d, opvert, a, b, mhalf, mfull, isclp, ipout)
            iside(3) = ipout
            go to 350
         endif
         if(iside(3).ne.0) go to 350
         call ipsign(x, y, z, x2, y2, z2, b, d, opvert, a, mhalf,
     *               mfull, isclp, ipout)
         iside(3) = ipout
  350    continue
c
         if(c.gt.8)then
            ifn = 0
            call ipsig3(x, y, z, x2, y2, z2, xc, yc, zc,
     *                 b, c, a, opvert, mhalf, mfull, isclp, ifn, ipout)
            iside(4) = ipout
         else
            call ipsig4(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  opvert, c, a, b, mhalf, mfull, isclp, ipout)
            iside(4) = ipout
            if(iside(4).ne.0) go to 360
            call ipsig6(x, y, z, x2, y2, z2, xc, yc, zc,
     *                  opvert, c, a, b, mhalf, mfull, isclp, ipout)
            iside(4) = ipout
            go to 360
         endif
         if(iside(4).ne.0) go to 360
         call ipsign(x, y, z, x2, y2, z2, b, opvert, c, a, mhalf,
     *               mfull, isclp, ipout)
         iside(4) = ipout
  360    continue
      endif
c
c     set sidist array
c     
      do 400 j = 2, 4
         if(iside(j) .gt. 0) then
            sidist(j) = 0
         elseif(iside(j) .lt. 0) then
            sidist(j) = -1
         else
            sidist(j) = 1
         endif
  400 continue
c
c     flip according to type of flip if possible
c
      if ((sidist(2). eq. 0) .and. (sidist(3). eq .0) .and.
     *    (sidist(4). eq. 0)) then
         if(opvert.gt.8) go to 420
         if(c.le.8.and.d.le.8)then
            call ipsig2(x, y, z, x2, y2, z2, xc, yc, zc, opvert, d,
     *                  opvert, c, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 420
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         elseif(c.le.8)then
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, c, opvert,
     *                  b, d, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 420
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         else
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, d, opvert,
     *                  c, b, a, c, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 420
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         endif
         if(itide .ge. 0) go to 1000
  420    continue
         call flip23 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, nvmax)
         go to 1000
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.0)) then 
         if(c.le.8.and.d.le.8) stop 1660
         if(opvert.gt.8) go to 440
         if(c.le.8)then
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, c, opvert,
     *                  b, d, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 440
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         else
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, d, opvert,
     *                  c, b, a, c, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 440
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         endif
         if(itide .ge. 0) go to 1000
  440    continue
         oddsid = 2
         call flip32 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.-1) .and.
     *        (sidist(4).eq.0)) then
         if(d.gt.8) go to 1000
         if(opvert.gt.8.and.c.gt.8) go to 460
         if(opvert.gt.8)then
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, d, c, b,
     *                  opvert, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 460
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         elseif(c.le.8)then
            call ipsig2(x, y, z, x2, y2, z2, xc, yc, zc, opvert, d,
     *                  opvert, c, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 460
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         else
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, d, opvert,
     *                  c, b, a, c, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 460
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         endif
         if(itide .ge. 0) go to 1000
  460    continue
         oddsid = 3
         call flip32 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.-1)) then
         if(c.gt.8) go to 1000
         if(opvert.gt.8.and.d.gt.8) go to 480
         if(opvert.gt.8)then
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, d, c, b,
     *                  opvert, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 480
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         elseif(d.le.8)then
            call ipsig2(x, y, z, x2, y2, z2, xc, yc, zc, opvert, d,
     *                  opvert, c, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 480
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         else
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, c, opvert,
     *                  b, d, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 480
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         endif
         if(itide .ge. 0) go to 1000
  480    continue
         oddsid = 4
         call flip32 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.0)) then
         if(c.le.8.and.d.le.8) stop 1670
         if(opvert.gt.8) go to 500
         if(c.le.8)then
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, c, opvert,
     *                  b, d, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 500
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         else
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, d, opvert,
     *                  c, b, a, c, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 500
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         endif
         if(itide .ge. 0) go to 1000
  500    continue
         oddsid = 2
         call flip22 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.1) .and.
     *        (sidist(4).eq.0)) then
         if(opvert.gt.8.and.d.gt.8)then
            call iqsig2(x, y, z, w, x2, y2, z2, w2, b, d, opvert, a,
     *           mhalf, mfull, isclp, isclw, isclr, delaun, itide)
            if(itide.ge.0) go to 1000
            go to 520
         endif
         if(opvert.gt.8) go to 520
         if(d.gt.8) go to 1000
         if(c.gt.8)then
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, d, opvert,
     *                  c, b, a, c, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 520
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         else
            call ipsig2(x, y, z, x2, y2, z2, xc, yc, zc, opvert, d,
     *                  opvert, c, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 520
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         endif
         if(itide .ge. 0) go to 1000
  520    continue
         oddsid = 3
         call flip22 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.1)) then
         if(opvert.gt.8.and.c.gt.8)then
            call iqsig2(x, y, z, w, x2, y2, z2, w2, b, opvert, c, a,
     *           mhalf, mfull, isclp, isclw, isclr, delaun, itide)
            if(itide.ge.0) go to 1000
            go to 540
         endif
         if(opvert.gt.8) go to 540
         if(c.gt.8) go to 1000
         if(d.gt.8)then
            call ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, c, opvert,
     *                  b, d, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 540
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         else
            call ipsig2(x, y, z, x2, y2, z2, xc, yc, zc, opvert, d,
     *                  opvert, c, a, b, mhalf, mfull, isclp, istt)
            if(istt.gt.0) go to 540
            if(istt.lt.0) go to 1000
            call iqsign(x, y, z, w, x2, y2, z2, w2, opvert, d, c, b,
     *           a, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         endif
         if(itide .ge. 0) go to 1000
  540    continue
         oddsid = 4
         call flip22 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
         go to 1000
      elseif(delaun) then
         go to 1000
      endif
c
      if ((sidist(2).eq.0) .and. (sidist(3).eq.-1) .and.
     *    (sidist(4).eq.-1)) then 
         oddsid = 2
         go to 900
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.-1)) then
         oddsid = 3
         go to 900
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.-1) .and.
     *        (sidist(4).eq.0)) then
         oddsid = 4
         go to 900
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.1) .and.
     *        (sidist(4).eq.1)) then
         oddsid = 2
         if(opvert.le.8)go to 900
         isodd = icon(6,now)
         if(isodd.le.8) stop 1680
         call iqsig1(x, y, z, w, x2, y2, z2, w2, isodd, opvert, k,
     *               mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         if(itide .ge. 0) go to 1000
         call flip21 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.1)) then
         oddsid = 3
         if(opvert.le.8)go to 900
         isodd = icon(7,now)
         if(isodd.le.8) stop 1690
         call iqsig1(x, y, z, w, x2, y2, z2, w2, isodd, opvert, k,
     *               mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         if(itide .ge. 0) go to 1000
         call flip21 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.1) .and.
     *        (sidist(4).eq.0)) then
         oddsid = 4
         if(opvert.le.8)go to 900
         isodd = icon(8,now)
         if(isodd.le.8) stop 1710
         call iqsig1(x, y, z, w, x2, y2, z2, w2, isodd, opvert, k,
     *               mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         if(itide .ge. 0) go to 1000
         call flip21 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.-1) .and.
     *        (sidist(4).eq.1)) then 
         oddsid = 2
         isite = icon(7,now)
         if(opvert.le.8 .or. isite.le.8)go to 900
         isodd = icon(6,now)
         if(isodd.le.8) stop 1720
         call iqsig2(x, y, z, w, x2, y2, z2, w2, isodd, isite, opvert,
     *        k, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         if(itide .ge. 0) go to 1000
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.-1)) then
         oddsid = 3
         isite = icon(8,now)
         if(opvert.le.8 .or. isite.le.8)go to 900
         stop 1730
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.1) .and.
     *        (sidist(4).eq.0)) then
         oddsid = 4
         isite = icon(6,now)
         if(opvert.le.8 .or. isite.le.8)go to 900
         isodd = icon(8,now)
         if(isodd.le.8) stop 1740
         call iqsig2(x, y, z, w, x2, y2, z2, w2, isodd, isite, opvert,
     *        k, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         if(itide .ge. 0) go to 1000
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.0) .and. (sidist(3).eq.1) .and.
     *        (sidist(4).eq.-1)) then
         oddsid = -2
         isite = icon(8,now)
         if(opvert.le.8 .or. isite.le.8)go to 800
         isodd = icon(6,now)
         if(isodd.le.8) stop 1750
         call iqsig2(x, y, z, w, x2, y2, z2, w2, isodd, isite, opvert,
     *        k, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         if(itide .ge. 0) go to 1000
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.-1) .and. (sidist(3).eq.0) .and.
     *        (sidist(4).eq.1)) then
         oddsid = -3
         isite = icon(6,now)
         if(opvert.le.8 .or. isite.le.8)go to 800
         isodd = icon(7,now)
         if(isodd.le.8) stop 1760
         call iqsig2(x, y, z, w, x2, y2, z2, w2, isodd, isite, opvert,
     *        k, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
         if(itide .ge. 0) go to 1000
         call flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax, flphis,
     *                tetra, newtts, is, oddsid, nvmax)
      elseif ((sidist(2).eq.1) .and. (sidist(3).eq.-1) .and.
     *(sidist(4).eq.0)) then
         oddsid = -4
         isite = icon(7,now)
         if(opvert.le.8 .or. isite.le.8) go to 800
         stop 1770
      endif
      go to 1000
c
  800 continue
      oddsid=-oddsid
  900 continue
      isite = icon(oddsid+4,now)
      if(isite.le.8) stop 1780
c
 1000 continue
      i = i + 1
      goto 100
c
 1500 continue
c
      return
      end
*FLIP23
c
c     This subroutine will perform a 2->3 flip.
c
      subroutine flip23 (icon, k, now, adj, ifl, ih, ihn, nhmax,
     *                   flphis, tetra, newtts, is, nvmax)
c
      integer icon(8,*), adj, ifl(*), tetra, is(*), ih(*)
      integer k, now, ihn, nhmax, newtts, nvmax
      logical flphis
      integer new1, new2, new3, i, next, j
c
      if(flphis) then
         new1 = tetra + 1
         new2 = tetra + 2
         new3 = tetra + 3
         tetra = new3
         if(tetra.gt.nvmax) stop 1810
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            tetra = new1
            if(tetra .gt. nvmax) stop 1820
         else
            new1 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new2 = tetra+1
            tetra = new2
            if(tetra .gt. nvmax) stop 1830
         else
            new2 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new3 = tetra+1
            tetra = new3
            if(tetra .gt. nvmax) stop 1840
         else
            new3 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      newtts = newtts + 3
      if(newtts.gt.nvmax) stop 1845
      ifl(newtts-2) = new1
      ifl(newtts-1) = new2
      ifl(newtts) = new3
c
c     create new1
c
      icon(1,new1) = icon(1,adj) 
      icon(2,new1) = new2
      icon(3,new1) = new3
      icon(4,new1) = icon(2,now)
      icon(5,new1) = k
      icon(6,new1) = icon(6,adj)
      icon(7,new1) = icon(7,adj)
      icon(8,new1) = icon(8,adj)
c
c     create new2
c
      icon(1,new2) = icon(2,adj)
      icon(2,new2) = new3
      icon(3,new2) = new1
      icon(4,new2) = icon(3,now)
      icon(5,new2) = k
      icon(6,new2) = icon(7,adj)
      icon(7,new2) = icon(5,adj)
      icon(8,new2) = icon(8,adj)
c
c     create new3
c
      icon(1,new3) = icon(3,adj)
      icon(2,new3) = new2
      icon(3,new3) = icon(4,now)
      icon(4,new3) = new1
      icon(5,new3) = k
      icon(6,new3) = icon(6,adj)
      icon(7,new3) = icon(8,adj)
      icon(8,new3) = icon(5,adj)
c
c     update neighboring tetrahedra
c
      do 100 i = 2, 4
         next = icon(i,now)
         if (next .eq. 0) goto 100
         do 50 j = 1,4
            if (icon(j,next) .eq. now) goto 60
   50    continue
         stop 1850
   60    continue
         if (i .eq. 2) then
            icon(j,next) = new1
         elseif (i .eq. 3) then
            icon(j,next) = new2
         else
            icon(j,next) = new3
         endif
  100 continue
c
      do 200 i = 1, 3
         next = icon(i,adj)
         if (next .eq. 0) goto 200
         do 150 j = 1,4
            if (icon(j,next) .eq. adj) goto 160
  150    continue
         stop 1860
  160    continue
         if (i .eq. 1) then
            icon(j,next) = new1
         elseif (i .eq. 2) then
            icon(j,next) = new2
         else
            icon(j,next) = new3
         endif
  200 continue
c
c     update is(*)
c
      is(icon(5,now)) = new3
      is(icon(6,now)) = new3
      is(icon(7,now)) = new3
      is(icon(8,now)) = new1
      is(icon(8,adj)) = new3
c
c     mark 2 old tetra to show children
c
      if(flphis) then
         icon(5,now) = -icon(5,now)
         icon(1,now) = new1 
         icon(2,now) = new2
         icon(3,now) = new3
         icon(4,now) = 0
c
         icon(5,adj) = -icon(5,adj)
         icon(1,adj) = new1 
         icon(2,adj) = new2
         icon(3,adj) = new3
         icon(4,adj) = 0
      else
         icon(5,now) = - icon(5,now)
         icon(5,adj) = - icon(5,adj)
         ihn = ihn+2
         if(ihn.gt.nhmax) stop 1870
         ih(ihn) = now
         ih(ihn-1) = adj
      endif
c
      return
      end
*FLIP32
c
c This subroutine will perform a 3->2 flip.  As a result, 3 new
c tetra will be created.
c
      subroutine flip32 (icon, k, now, adj, ifl, ih, ihn, nhmax,
     *                   flphis, tetra, newtts, is, negsid, nvmax)
c
      integer icon(8,*), adj, ifl(*), tetra, is(*), ih(*), site1, site2
      integer k, now, ihn, nhmax, newtts, negsid, nvmax
      logical flphis
      integer nxtnow, nxtadj, i, new1, new2, neigh, j
c
c     reorder now
c
      site1 = icon(negsid+4, now)
      call reordr (icon, k, site1, now)
c
c     check if now & adj have same neighbor, reorder adj
c
      nxtnow = icon(2,now)
      if(icon(5,nxtnow).ne.k) stop 1910
c
      call sitord (icon, site1, adj)
      nxtadj = icon(1,adj)
      if (nxtnow .ne. nxtadj) goto 1000
c
      do 210 i = 1, 4
         if (icon(i,adj) .eq. now) then
            site2 = icon(i+4,adj)
            goto 215
         endif
  210 continue
      stop 1920
  215 continue
      call reordr (icon, site1, site2, adj)
c
c     reorder nxtnow
c
      call reordr (icon, k, site2, nxtnow)
c
      if(flphis) then
         new1 = tetra + 1
         new2 = tetra + 2
         tetra = new2
         if (tetra .gt. nvmax) stop 1930
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            tetra = new1
            if(tetra .gt. nvmax) stop 1940
         else
            new1 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new2 = tetra+1
            tetra = new2
            if(tetra .gt. nvmax) stop 1950
         else
            new2 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      newtts = newtts + 2
      if(newtts.gt.nvmax) stop 1955
      ifl(newtts-1) = new1
      ifl(newtts) = new2
c
c     create new1
c
      icon(1,new1) = icon(4,adj)
      icon(2,new1) = icon(3,nxtnow)
      icon(3,new1) = new2
      icon(4,new1) = icon(4,now)
      icon(5,new1) = k
      icon(6,new1) = icon(6,now)
      icon(7,new1) = icon(7,now)
      icon(8,new1) = icon(6,adj)
c
c     create new2
c
      icon(1,new2) = icon(3,adj)
      icon(2,new2) = icon(4,nxtnow)
      icon(3,new2) = icon(3,now)
      icon(4,new2) = new1
      icon(5,new2) = k 
      icon(6,new2) = icon(6,now)
      icon(7,new2) = icon(6,adj)
      icon(8,new2) = icon(8,now)
c
c     update neighboring tetrahedra
c
      do 400 i = 3, 4
         neigh = icon(i,now)
         if (neigh .eq. 0) goto 400 
         do 350 j = 1, 4
            if (icon(j,neigh) .eq. now) goto 375
  350    continue
         stop 1960
  375    continue
         if (i .eq. 3) then
            icon(j,neigh) = new2
         else
            icon(j,neigh) = new1
         endif
  400 continue
c
      do 500 i = 3, 4
         neigh = icon(i,adj)
         if (neigh .eq. 0) goto 500
         do 450 j = 1, 4
            if (icon(j,neigh) .eq. adj) goto 475
  450    continue
         stop 1970
  475    continue
         if (i .eq. 3) then
            icon(j,neigh) = new2
         else
            icon(j,neigh) = new1
         endif
  500 continue
c
      do 600 i = 3, 4
         neigh = icon(i,nxtnow)
         if (neigh .eq. 0) goto 600
         do 550 j = 1, 4
            if (icon(j,neigh) .eq. nxtnow) goto 575
  550    continue
         stop 1980
  575    continue
         if (i .eq. 3) then
            icon(j,neigh) = new1
         else
            icon(j,neigh) = new2
         endif
  600 continue
c
c     update is(*)
c
      is(icon(5,now)) = new1
      is(icon(6,now)) = new1
      is(icon(7,now)) = new1
      is(icon(8,now)) = new2
      is(icon(6,adj)) = new1
c
c     show children of adj, now, nxtnow
c
      if(flphis) then
         icon(5,now) = -icon(5,now)
         icon(1,now) = new1
         icon(2,now) = new2
         icon(3,now) = 0
         icon(4,now) = 0
c
         icon(5,adj) = -icon(5,adj)
         icon(1,adj) = new1
         icon(2,adj) = new2
         icon(3,adj) = 0
         icon(4,adj) = 0
c
         icon(5,nxtnow) = -icon(5,nxtnow)
         icon(1,nxtnow) = new1
         icon(2,nxtnow) = new2
         icon(3,nxtnow) = 0 
         icon(4,nxtnow) = 0
      else
         icon(5,now) = - icon(5,now)
         icon(5,adj) = - icon(5,adj)
         icon(5,nxtnow) = -icon(5,nxtnow)
         ihn = ihn+3
         if(ihn.gt.nhmax) stop 1990
         ih(ihn) = now
         ih(ihn-1) = adj
         ih(ihn-2) = nxtnow
      endif
c
 1000 continue
      return
      end
*FLIP22
c
c This subroutine will perform a 2->2 flip.  Four new tetra will be
c created.
c
      subroutine flip22 (icon, k, now, adj, ifl, ih, ihn, nhmax,
     *                   flphis, tetra, newtts, is, zersid, nvmax)
c
      integer icon(8,*), adj, ifl(*), tetra, is(*), ih(*), zersid
      integer k, now, ihn, nhmax, newtts, nvmax, site1, site2
      logical flphis
      integer nxtnow, nxtadj, i, new1, new2, new3, new4, neigh, j
c
c     reorder now 
c
      site1 = icon(zersid+4, now)
      call reordr (icon, k, site1, now)
c
c     define nxtnow
c
      nxtnow = icon(2,now)
      if(icon(5,nxtnow).ne.k) stop 2010
c
c     reorder adj
c
      call sitord (icon, site1, adj)
c
c     define nxtadj
c
      nxtadj = icon(1,adj)
c
c     are nxtnow and nxtadj neighbors?
c
      do 5 i = 1, 4
         if (icon(i,nxtnow) .eq. nxtadj) goto 6
    5 continue
      goto 2000
    6 continue
c
      if(flphis) then
         new1 = tetra + 1
         new2 = tetra + 2
         new3 = tetra + 3
         new4 = tetra + 4
         tetra = new4
         if (tetra .gt. nvmax) stop 2020
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            tetra = new1
            if(tetra .gt. nvmax) stop 2030
         else
            new1 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new2 = tetra+1
            tetra = new2
            if(tetra .gt. nvmax) stop 2040
         else
            new2 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new3 = tetra+1
            tetra = new3
            if(tetra .gt. nvmax) stop 2050
         else
            new3 = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            new4 = tetra+1
            tetra = new4
            if(tetra .gt. nvmax) stop 2060
         else
            new4 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      newtts = newtts+4
      if(newtts.gt.nvmax) stop 2065
      ifl(newtts-3) = new1
      ifl(newtts-2) = new2
      ifl(newtts-1) = new3
      ifl(newtts) = new4
c
c     reorder adj, nxtnow, nxtadj
c
      do 10 i = 2, 4
         if (icon(i,adj) .eq. now) then
            site2 = icon(i+4,adj)
            goto 15
         endif
   10 continue
      stop 2070
   15 continue
      call reordr (icon, site1, site2, adj)
c
      do 20 i = 2, 4
         if (icon(i,nxtnow) .eq. now) then
            site1 = icon(i+4,nxtnow)
            goto 25
         endif
  20  continue
      stop 2080
  25  continue
      call reordr (icon, k, site1, nxtnow)
c
      call reordr (icon, site1, site2, nxtadj)
c
c     create new1
c
      icon(1,new1) = icon(4,adj)
      icon(2,new1) = new3
      icon(3,new1) = new2
      icon(4,new1) = icon(4,now)
      icon(5,new1) = k
      icon(6,new1) = icon(6,now)
      icon(7,new1) = icon(7,now)
      icon(8,new1) = icon(6,adj)
c
c     create new2
c
      icon(1,new2) = icon(3,adj)
      icon(2,new2) = new4
      icon(3,new2) = icon(3,now)
      icon(4,new2) = new1
      icon(5,new2) = k
      icon(6,new2) = icon(6,now)
      icon(7,new2) = icon(6,adj)
      icon(8,new2) = icon(8,now)
c
c     create new3
c
      icon(1,new3) = icon(3,nxtadj)
      icon(2,new3) = new4
      icon(3,new3) = new1
      icon(4,new3) = icon(3,nxtnow)
      icon(5,new3) = k
      icon(6,new3) = icon(7,now)
      icon(7,new3) = icon(6,nxtnow)
      icon(8,new3) = icon(6,adj)
c
c     create new4
c
      icon(1,new4) = icon(4,nxtadj)
      icon(2,new4) = new2
      icon(3,new4) = new3
      icon(4,new4) = icon(4,nxtnow)
      icon(5,new4) = k
      icon(6,new4) = icon(6,nxtnow)
      icon(7,new4) = icon(8,now)
      icon(8,new4) = icon(6,adj)
c
c     update is(*)
c
      is(icon(5,now)) = new1
      is(icon(6,now)) = new1
      is(icon(7,now)) = new1
      is(icon(8,now)) = new2
      is(icon(6,adj)) = new1
      is(icon(6,nxtnow)) = new3
c
c     update neighbors of now
c
      do 100 i = 3, 4
         neigh = icon(i,now)
         if (neigh .eq. 0) goto 100
         do 50 j = 1, 4
            if (icon(j,neigh) .eq. now) goto 75
   50    continue
         stop 2090
   75    continue
         if (i .eq. 3) then
            icon(j,neigh) = new2
         else
            icon(j,neigh) = new1
         endif
  100 continue
c
c     update neighbors of adj
c
      do 400 i = 3, 4
         neigh = icon(i,adj)
         if (neigh .eq. 0) goto 400
         do 350 j = 1, 4
            if (icon(j,neigh) .eq. adj) goto 375
  350    continue
         stop 2110
  375    continue
         if (i .eq. 3) then
            icon(j,neigh) = new2
         else
            icon(j,neigh) = new1
         endif
  400 continue
c
c     update neighbors of nxtnow
c
      do 600 i = 3, 4
         neigh = icon(i,nxtnow)
         if (neigh .eq. 0) goto 600
         do 575 j = 1, 4
            if (icon(j,neigh) .eq. nxtnow) goto 590
  575    continue
         stop 2120
  590    continue
         if (i .eq. 3) then
            icon(j,neigh) = new3
         else
            icon(j,neigh) = new4
         endif
  600    continue   
c
c     update neighbors of nxtadj
c
      do 900 i = 3, 4
         neigh = icon(i,nxtadj)
         if (neigh. eq. 0) goto 900
         do 875 j = 1, 4
            if (icon(j,neigh) .eq. nxtadj) goto 890
  875    continue
         stop 2130
  890    continue
         if (i .eq. 3) then
            icon(j,neigh) = new3
         else
            icon(j,neigh) = new4
         endif
  900 continue
c
c     show children of old tetra
c
      if(flphis) then
         icon(5,now) = -icon(5,now)
         icon(1,now) = new1
         icon(2,now) = new2
         icon(3,now) = -nxtnow
         icon(4,now) = 0
c
         icon(5,adj) = -icon(5,adj)
         icon(1,adj) = new1
         icon(2,adj) = new2
         icon(3,adj) = -nxtadj
         icon(4,adj) = 0
c
         icon(5,nxtnow) = -icon(5,nxtnow)
         icon(1,nxtnow) = new3 
         icon(2,nxtnow) = new4
         icon(3,nxtnow) = -now
         icon(4,nxtnow) = 0
c
         icon(5,nxtadj) = -icon(5,nxtadj)
         icon(1,nxtadj) = new3 
         icon(2,nxtadj) = new4
         icon(3,nxtadj) = -adj
         icon(4,nxtadj) = 0
      else
         icon(5,now) = - icon(5,now)
         icon(5,adj) = - icon(5,adj)
         icon(5,nxtnow) = -icon(5,nxtnow)
         icon(5,nxtadj) = -icon(5,nxtadj)
         ihn = ihn+4
         if(ihn.gt.nhmax) stop 2140
         ih(ihn) = now
         ih(ihn-1) = adj
         ih(ihn-2) = nxtnow
         ih(ihn-3) = nxtadj
      endif
c
 2000 continue
      return
      end
*FLIP41
c
c This subroutine will perform a 4->1 flip. 1 tetrahedron will be
c created from 4 tetrahedra
c
      subroutine flip41 (icon, k, now, adj, ifl, ih, ihn, nhmax,
     *                   flphis, tetra, newtts, is, zersid, nvmax)
c
      integer icon(8,*), ifl(*), is(*), ih(*)
      integer k, now, ihn, nhmax, newtts, nvmax
      integer adj, tetra, zersid, site1, site2
      logical flphis
      integer nxtnow, nxtadj, new1, neigh, j
c
c     reorder now
c
      site1 = icon(zersid+4,now)
      call reordr(icon, k, site1, now)
      site2 = icon(7,now)
c
c     reorder adj
c
      call reordr(icon, site1, site2, adj)
c
c     define nxtnow and nxtadj
c
      nxtnow = icon(4,now)
      nxtadj = icon(3,now)
c
c     do now, adj, nxtnow, nxtadj form a tetrahedron?
c
      if(icon(3,adj).ne.nxtnow .or. icon(2,adj).ne.nxtadj) go to 2000
c
c     flip
c
      if(flphis) then
         new1 = tetra + 1
         tetra = new1
         if(tetra .gt. nvmax) stop 2210
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            tetra = new1
            if(tetra .gt. nvmax) stop 2220
         else
            new1 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      newtts = newtts + 1
      if(newtts.gt.nvmax) stop 2225
      ifl(newtts) = new1
c
c     reorder nxtnow and nxtadj
c
      call reordr(icon, k, site2, nxtnow)
      call reordr(icon, k, site1, nxtadj)
      if(icon(2,nxtnow) .ne. nxtadj) stop 2230
c
c     create tetra
c
      icon(1,new1) = icon(1,adj)
      icon(2,new1) = icon(2,nxtadj)
      icon(3,new1) = icon(3,nxtnow)
      icon(4,new1) = icon(2,now)
      icon(5,new1) = k
      icon(6,new1) = site2
      icon(7,new1) = icon(8,now)
      icon(8,new1) = icon(8,adj)
c
c     update is(*)
c
      is(site1) = -4
      is(icon(5,new1)) = new1 
      is(icon(6,new1)) = new1
      is(icon(7,new1)) = new1
      is(icon(8,new1)) = new1
c
c     update neighbor of now
c
      neigh = icon(2,now)
      if(neigh .eq. 0) go to 200
      do 140 j = 1, 4
         if(icon(j,neigh) .eq. now) go to 160
  140 continue
      stop 2240
  160 continue
      icon(j,neigh) = new1
  200 continue
c
c     update neighbor of adj
c
      neigh = icon(1,adj)
      if(neigh .eq. 0) go to 300
      do 240 j = 1, 4
         if(icon(j,neigh) .eq. adj) go to 260
  240 continue
      stop 2250
  260 continue
      icon(j,neigh) = new1
  300 continue
c
c     update neighbor of nxtnow
c
      neigh = icon(3,nxtnow)
      if(neigh .eq. 0) go to 400
      do 340 j = 1, 4
         if(icon(j,neigh) .eq. nxtnow) go to 360
  340 continue
      stop 2260
  360 continue
      icon(j,neigh) = new1
  400 continue
c
c     update neighbor of nxtadj
c
      neigh = icon(2,nxtadj)
      if(neigh .eq. 0) go to 500
      do 440 j = 1, 4
         if(icon(j,neigh) .eq. nxtadj) go to 460
  440 continue
      stop 2270
  460 continue
      icon(j,neigh) = new1
  500 continue
c
c     show children of old tetrahedra
c
      if(flphis) then
         icon(5,now) = -icon(5,now)
         icon(1,now) = new1
         icon(2,now) = 0
         icon(3,now) = 0
         icon(4,now) = 0
c
         icon(5,adj) = -icon(5,adj)
         icon(1,adj) = new1
         icon(2,adj) = 0
         icon(3,adj) = 0
         icon(4,adj) = 0
c
         icon(5,nxtnow) = -icon(5,nxtnow)
         icon(1,nxtnow) = new1
         icon(2,nxtnow) = 0
         icon(3,nxtnow) = 0
         icon(4,nxtnow) = 0
c
         icon(5,nxtadj) = -icon(5,nxtadj)
         icon(1,nxtadj) = new1
         icon(2,nxtadj) = 0
         icon(3,nxtadj) = 0
         icon(4,nxtadj) = 0
      else
         icon(5,now) = - icon(5,now)
         icon(5,adj) = - icon(5,adj)
         icon(5,nxtnow) = -icon(5,nxtnow)
         icon(5,nxtadj) = -icon(5,nxtadj)
         ihn = ihn+4
         if(ihn.gt.nhmax) stop 2280
         ih(ihn) = now
         ih(ihn-1) = adj
         ih(ihn-2) = nxtnow
         ih(ihn-3) = nxtadj
      endif
c
 2000 continue
      return
      end
*FLIP21
c
c This subroutine will perform 2->1 flips. 1 tetrahedron will
c be created from 2 tetrahedra for each flip
c
      subroutine flip21 (icon, k, now, adj, ifl, ih, ihn, nhmax,
     *                   flphis, tetra, newtts, is, zersid, nvmax)
c
      integer icon(8,*), ifl(*), is(*), ih(*)
      integer k, now, ihn, nhmax, newtts, nvmax
      integer adj, tetra, zersid, site1, site2, site3, site4
      logical flphis
      integer nxtnow, nxtadj, initet, new1, nel1, neigh, j
      integer lstnow, lstadj
c
c     reorder now
c
      site1 = icon(zersid+4,now)
      call reordr(icon, k, site1, now)
      site2 = icon(7,now)
      site3 = icon(8,now)
c
c     reordr adj
c
      call reordr(icon, site1, site2, adj)
c
      nxtnow = icon(3,now)
      nxtadj = icon(2,adj)
      if(nxtnow.eq.nxtadj) stop 2310
c
c     go around edge to test for flipping
c
  100 continue
      call reordr(icon, k, site1, nxtnow)
      call reordr(icon, site1, site3, nxtadj)
      if(icon(1,nxtnow).ne.nxtadj) go to 2000
      site3 = icon(8,nxtnow)
      nxtnow = icon(3,nxtnow)
      nxtadj = icon(2,nxtadj)
      if(nxtnow .eq. now) go to 200
      if(nxtadj .eq. adj) stop 2320
      go to 100
c
c     flip
c
  200 continue
c
      if(nxtadj .ne. adj) stop 2330
      if(flphis) then
         initet = tetra+1
      else
         if(ihn.eq.0) then
            initet = tetra+1
         else
            initet = ih(ihn)
         endif
      endif
      new1 = initet
      site4 = icon(8,adj)
c
c     go around edge for creating new tetrahedra
c
  300 continue
      nel1 = new1
      if(flphis) then
         new1 = tetra + 1
         tetra = new1
         if(tetra .gt. nvmax) stop 2340
      else
         if(ihn.eq.0) then
            new1 = tetra+1
            tetra = new1
            if(tetra .gt. nvmax) stop 2350
         else
            new1 = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      newtts = newtts + 1
      if(newtts.gt.nvmax) stop 2355
      ifl(newtts) = new1
c
c     create tetra
c
      icon(1,new1) = icon(1,nxtadj)
      icon(2,new1) = icon(2,nxtnow)
      icon(3,nel1) = new1
      icon(4,new1) = nel1
      icon(5,new1) = k
      icon(6,new1) = site4
      icon(7,new1) = icon(7,nxtnow)
      icon(8,new1) = icon(8,nxtnow)
c
c     update is(*)
c
      is(icon(8,new1)) = new1
c
c     update neighbor of nxtnow
c
      neigh = icon(2,nxtnow)
      if(neigh .eq. 0) go to 400
      do 340 j = 1, 4
         if(icon(j,neigh) .eq. nxtnow) go to 360
  340 continue
      stop 2360
  360 continue
      icon(j,neigh) = new1
  400 continue
c
c     update neighbor of nxtadj
c
      neigh = icon(1,nxtadj)
      if(neigh .eq. 0) go to 500
      do 440 j = 1, 4
         if(icon(j,neigh) .eq. nxtadj) go to 460
  440 continue
      stop 2370
  460 continue
      icon(j,neigh) = new1
  500 continue
c
c     show children of nxtnow, nxtadj
c
      lstnow = nxtnow
      lstadj = nxtadj
      nxtnow = icon(3,lstnow)
      nxtadj = icon(2,lstadj)
c
      if(flphis) then
         icon(5,lstnow) = -icon(5,lstnow)
         icon(1,lstnow) = new1
         icon(2,lstnow) = 0
         icon(3,lstnow) = -nxtnow
         icon(4,lstnow) = -icon(4,lstnow)
c
         icon(5,lstadj) = -icon(5,lstadj)
         icon(1,lstadj) = new1
         icon(2,lstadj) = -nxtadj
         icon(3,lstadj) = -icon(3,lstadj)
         icon(4,lstadj) = 0
      else
         icon(5,lstnow) = -icon(5,lstnow)
         icon(5,lstadj) = -icon(5,lstadj)
         ihn = ihn+2
         if(ihn.gt.nhmax) stop 2380
         ih(ihn) = lstnow
         ih(ihn-1) = lstadj
      endif
c
      if(nxtnow .ne. now) go to 300
      icon(3,new1) = initet
      icon(4,initet) = new1
      is(k) = new1
      is(site4) = new1
      is(site1) = -4
c
 2000 continue
      return
      end
*FLIP31
c
c This subroutine will perform 3->1 flips. 1 tetrahedron will
c be created from 3 tetrahedra for each flip
c
      subroutine flip31 (icon, k, now, adj, ifl, ih, ihn, nhmax,
     *                   flphis, tetra, newtts, is, zersid, nvmax)
c
      integer icon(8,*), ifl(*), is(*), ih(*)
      integer k, now, ihn, nhmax, newtts, nvmax
      integer adj, tetra, zersid
      integer site1, site2, site3, site4, site5
      logical flphis
      integer nside, nrt, nrd, nxtnow, nxtadj, nxtnrd
      integer nitet, nitat, neigh, j
c
      nside = iabs(zersid)
c
c     reorder now
c
      site1 = icon(nside+4,now)
      call reordr(icon, k, site1, now)
c
      if(zersid .gt. 0) go to 10
c
c     define nrt
c
      nrt = icon(4,now)
c
c     reordr nrt
c
      call reordr(icon, k, site1, nrt)
      if(icon(1,nrt) .ne. adj) go to 2000
c
c     define nrd, redefine now
c
      nrd = now
      now = nrt
      go to 20
c
c     define nrd
c
   10 continue
      nrd = icon(3,now)
c
c     reorder nrd
c
      call reordr(icon, k, site1, nrd)
      if(icon(1,nrd) .ne. adj) go to 2000
c
   20 continue
c
c     reorder adj
c
      site2 = icon(7,now)
      call reordr(icon, site1, site2, adj)
c
c     define nxtnow, nxtadj, nxtnrd, and reorder
c
      nxtnow = icon(4,now)
      nxtadj = icon(3,adj)
      nxtnrd = icon(3,nrd)
      call reordr(icon, k, site1, nxtnow)
      call reordr(icon, site1, site2, nxtadj)
      call reordr(icon, k, site1, nxtnrd)
      site5 = icon(7,nxtnow)
      if(icon(8,nxtadj) .ne. site5 .or.
     *   icon(8,nxtnrd) .ne. site5) go to 2000
c
c     flip
c
      site3 = icon(8,now)
      site4 = icon(8,adj)
      if(flphis) then
         nitet = tetra + 1
         nitat = tetra + 2
         tetra = nitat
         if(tetra .gt. nvmax) stop 2410
      else
         if(ihn.eq.0) then
            nitet = tetra+1
            tetra = nitet
            if(tetra .gt. nvmax) stop 2420
         else
            nitet = ih(ihn)
            ihn = ihn-1
         endif
         if(ihn.eq.0) then
            nitat = tetra+1
            tetra = nitat
            if(tetra .gt. nvmax) stop 2430
         else
            nitat = ih(ihn)
            ihn = ihn-1
         endif
      endif
c
      newtts = newtts + 2
      if(newtts.gt.nvmax) stop 2435
      ifl(newtts-1) = nitet
      ifl(newtts) = nitat
c
c     create new tetrahedra
c
      icon(1,nitet) = icon(1,adj)
      icon(2,nitet) = icon(2,nrd)
      icon(3,nitet) = nitat
      icon(4,nitet) = icon(2,now)
      icon(5,nitet) = k
      icon(6,nitet) = site2
      icon(7,nitet) = site3
      icon(8,nitet) = site4
c
      icon(1,nitat) = icon(1,nxtadj)
      icon(2,nitat) = icon(2,nxtnrd)
      icon(3,nitat) = icon(2,nxtnow)
      icon(4,nitat) = nitet
      icon(5,nitat) = k
      icon(6,nitat) = site2
      icon(7,nitat) = site4
      icon(8,nitat) = site5
c
c     update is(*)
c
      is(k) = nitat
      is(site1) = -4
      is(site2) = nitat
      is(site3) = nitet
      is(site4) = nitat
      is(site5) = nitat
c
c     update neighbors of adj, nrd, now
c
      neigh = icon(2,now)
      if(neigh .eq. 0) go to 200
      do 140 j = 1, 4
         if(icon(j,neigh) .eq. now) go to 160
  140 continue
      stop 2440
  160 continue
      icon(j,neigh) = nitet
  200 continue
c
      neigh = icon(1,adj)
      if(neigh .eq. 0) go to 300
      do 240 j = 1, 4
         if(icon(j,neigh) .eq. adj) go to 260
  240 continue
      stop 2450
  260 continue
      icon(j,neigh) = nitet
  300 continue
c
      neigh = icon(2,nrd)
      if(neigh .eq. 0) go to 330
      do 310 j = 1, 4
         if(icon(j,neigh) .eq. nrd) go to 320
  310 continue
      stop 2460
  320 continue
      icon(j,neigh) = nitet
  330 continue
c
c     update neighbors of nxtnow, nxtadj, nxtnrd
c
      neigh = icon(2,nxtnow)
      if(neigh .eq. 0) go to 400
      do 340 j = 1, 4
         if(icon(j,neigh) .eq. nxtnow) go to 360
  340 continue
      stop 2470
  360 continue
      icon(j,neigh) = nitat
  400 continue
c
      neigh = icon(1,nxtadj)
      if(neigh .eq. 0) go to 500
      do 440 j = 1, 4
         if(icon(j,neigh) .eq. nxtadj) go to 460
  440 continue
      stop 2480
  460 continue
      icon(j,neigh) = nitat
  500 continue
c
      neigh = icon(2,nxtnrd)
      if(neigh .eq. 0) go to 600
      do 540 j = 1, 4
         if(icon(j,neigh) .eq. nxtnrd) go to 560
  540 continue
      stop 2490
  560 continue
      icon(j,neigh) = nitat
  600 continue
c
c     show children of old tetrahedra
c
      if(flphis) then
         icon(5,now) = -icon(5,now)
         icon(1,now) = nitet
         icon(2,now) = -nxtnow
         icon(3,now) = 0
         icon(4,now) = 0
c
         icon(5,adj) = -icon(5,adj)
         icon(1,adj) = nitet
         icon(2,adj) = -nxtadj
         icon(3,adj) = 0
         icon(4,adj) = 0
c
         icon(5,nrd) = -icon(5,nrd)
         icon(1,nrd) = nitet
         icon(2,nrd) = -nxtnrd
         icon(3,nrd) = 0
         icon(4,nrd) = 0
c
         icon(5,nxtnow) = -icon(5,nxtnow)
         icon(1,nxtnow) = nitat
         icon(2,nxtnow) = -now
         icon(3,nxtnow) = 0
         icon(4,nxtnow) = 0
c
         icon(5,nxtadj) = -icon(5,nxtadj)
         icon(1,nxtadj) = nitat
         icon(2,nxtadj) = -adj
         icon(3,nxtadj) = 0
         icon(4,nxtadj) = 0
c
         icon(5,nxtnrd) = -icon(5,nxtnrd)
         icon(1,nxtnrd) = nitat
         icon(2,nxtnrd) = -nrd
         icon(3,nxtnrd) = 0
         icon(4,nxtnrd) = 0
      else
         icon(5,now) = -icon(5,now)
         icon(5,adj) = -icon(5,adj)
         icon(5,nrd) = -icon(5,nrd)
         icon(5,nxtnow) = -icon(5,nxtnow)
         icon(5,nxtadj) = -icon(5,nxtadj)
         icon(5,nxtnrd) = -icon(5,nxtnrd)
         ihn = ihn+6
         if(ihn.gt.nhmax) stop 2495
         ih(ihn) = now
         ih(ihn-1) = adj
         ih(ihn-2) = nrd
         ih(ihn-3) = nxtnow
         ih(ihn-4) = nxtadj
         ih(ihn-5) = nxtnrd
      endif
c
 2000 continue
      return
      end
*CONVEX
c
c This subroutine will classify all tetra in array ifl, where: 
c 	ifl(curr) = -1 -> tetra has children
c 	ifl(curr) =  0 -> tetra is outside convex hull
c 	ifl(curr) =  1 -> tetra is inside convex hull
c Then, a verification of the surface's convexity will be run. 
c
      subroutine convex (icon, tetra, ifl, xi, yi, zi, x, y, z, x2, y2,
     *                   z2, idmin, mhalf, mfull, isclp, epz)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*)
      integer icon(8,*), tetra, ifl(*), curr, adj
      double precision epz
      integer isclp(*), mhalf, mfull
      integer idmin, i, j, nop, now, iside
      integer ikon(8,1), site1, site2, site3, a, b, c, vert
c
      idmin = 0
c
c     classify all tetra 
c
      do 100 curr = 1, tetra
         if (icon(5,curr) .lt. 0) then
            ifl(curr) = -1
         elseif ((icon(5,curr) .le. 8) .or. (icon(6,curr) .le. 8) .or.
     *   (icon(7,curr) .le. 8) .or. (icon(8,curr) .le. 8)) then
            ifl(curr) = 0
         else
            ifl(curr) = 1
         endif
  100 continue
c
c     take all tetra s.t. ifl(tetra) = 1 and check convexity
c
      do 300 curr = 1, tetra
         if (ifl(curr) .ne. 1) goto 300
         do 200 i = 1, 4
            adj = icon(i,curr)
            if (adj .eq. 0) stop 2510
            if (ifl(adj) .ne. 0) goto 200
            do 150 j = 1, 8
               ikon(j,1) = icon(j,curr)
  150       continue
            site1 = icon(i+4,curr) 
            call sitord (ikon, site1, 1)
            a = ikon(6,1)
            b = ikon(7,1)
            c = ikon(8,1)
            site1 = a
            site2 = b
            site3 = c
            call reordr (ikon, site1, site2, 1)
            do 175 j = 1, 3
               nop = curr
               now = ikon(4,1)
  160          continue
               if (now. eq. curr .or. now .eq. 0) stop 2520
               call reordr (icon, site1, site2, now)
               if (ifl(now) .ne. 0) then
                  nop = now
                  now = icon(4,nop)
                  goto 160
               else
                  if(nop .eq. curr) go to 170
                  vert = icon(8,now)
                  call irsign(xi, yi, zi, x, y, z, x2, y2, z2, vert,
     *                        site1, site2, site3, mhalf, mfull,
     *                        isclp, epz, iside)
                  if(iside .gt. 0) idmin = idmin+1
               endif
c
  170          continue
               if (j .eq. 1) then
                  site1 = b
                  site2 = c
                  site3 = a
                  call reordr (ikon, site1, site2, 1)
               endif
               if (j .eq. 2) then
                  site1 = c
                  site2 = a
                  site3 = b
                  call reordr (ikon, site1, site2, 1)
               endif
  175       continue
  200    continue
  300 continue
c
      return
      end
*DELCHK
c
c This subroutine will test how well the Regular/Delaunay property is
c satisfied by the tetrahedra inside convex hull of point set
c
      subroutine delchk(tetra, icon, ifl, xi, yi, zi, wi, x, y, z, w,
     *                  x2, y2, z2, w2, idmax, delaun, mhalf, mfull,
     *                  isclp, isclw, isclr, epz)
c
      double precision xi(*), yi(*), zi(*), wi(*)
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer icon(8,*), ifl(*)
      double precision epz, tdist, xctr, yctr, zctr
      integer idmax, mhalf, mfull, i, j, k, isite, ipossi, itide
      integer isclp(*), isclw(*), isclr(*)
      integer ikon(8,1), site1, site2, site3, tetra
      integer opvert, a, b, c, d, adj
      logical delaun
c
c     initialize
c
      idmax = 0
c
c     test all tetra with ifl(tetra)=1 for the Regular/Delaunay property
c
      do 200 i = 1, tetra
         if(ifl(i) .ne. 1) go to 200
         a = icon(5,i)
         b = icon(6,i)
         c = icon(7,i)
         d = icon(8,i)
         do 100 j = 1, 4
            adj = icon(j,i)
            if (adj .eq. 0) stop 2610
            if(ifl(adj).ne.1) go to 100
            if(adj.gt.i) go to 100
            do 50 k = 1, 8
               ikon(k,1) = icon(k,i)
   50       continue
            isite = icon(j+4,i)
            call sitord (ikon, isite, 1)
c
            site1 = ikon(6,1)
            site2 = ikon(7,1)
            site3 = ikon(8,1)
            call reordr (icon, site1, site2, adj)
            if (icon(7,adj) .ne. site3) stop 2620
            if (icon(4,adj) .ne. i) stop 2630
c
            opvert = icon(8,adj)
            call ctrad(xi, yi, zi, wi, xctr, yctr, zctr, site1, site2,
     *                 site3, opvert, epz, delaun, ipossi)
            if(ipossi.eq.1) go to 70
            call bisphr(xi, yi, zi, wi, isite, site1, epz, xctr, yctr,
     *                  zctr, tdist, delaun, ipossi)
            if(ipossi.eq.1) go to 70
            if(tdist.le.0.0d0) go to 100
            go to 80
   70       continue
            call iqsign(x, y, z, w, x2, y2, z2, w2, a, b, c, d,
     *      opvert, mhalf, mfull, isclp, isclw, isclr, delaun, itide)
            if(itide.ge.0) go to 100
   80       continue
            idmax = idmax+1
  100    continue 
  200 continue
c
      return
      end
*REVTET
c
c This subroutine will compress data structure in order to save space
c by eliminating artificial and discarded tetrahedra
c
      subroutine revtet(tetra, tetru, icon, nv, is, ifl, flphis)
c
      integer tetra, tetru, icon(8,*), is(*), ifl(*)
      integer nv, i, j, ii, ielm
      logical flphis
c
c     identify true tetrahedra
c
      tetru = 0
      ielm = 0
      do 100 i = 1, tetra
         if ((icon(5,i).le.0) .or. (icon(6,i).le.0) .or.
     *       (icon(7,i).le.0) .or. (icon(8,i).le.0)) then
            ielm = ielm + 1
            ifl(i) = 0
         elseif ((icon(5,i).le.8) .or. (icon(6,i).le.8) .or.
     *           (icon(7,i).le.8) .or. (icon(8,i).le.8)) then
            ifl(i) = 0
         else
            tetru=tetru+1
            ifl(i) = 1
         endif
  100 continue
      if(tetru.eq.0) go to 1000
c
c     zero out nonexistent tetrahedra in icon
c
      do 300 i=1,tetra
         if(ifl(i).eq.0) go to 300
         do 200 j=1,4
            if(icon(j,i).le.0.or.icon(j,i).gt.tetra) stop 2710
            if(ifl(icon(j,i)).eq.0) icon(j,i)=0
  200    continue
  300 continue
c
c     compress icon
c
      ii=0
      do 500 i=1,tetra
         if(ifl(i).eq.0) go to 500
         ii=ii+1
         ifl(i)=ii
         do 400 j=1,8
            icon(j,ii)=icon(j,i)
  400    continue
  500 continue
c
c     update icon for tetrahedra and is for vertices
c
      do 550 i=9,nv
         if(is(i).gt.0)then
            is(i-8)=1
         else
            is(i-8)=is(i)
         endif
  550 continue
      do 800 i=1,tetru
         do 600 j=1,4
            if(icon(j,i).eq.0) go to 600
            icon(j,i)=ifl(icon(j,i))
  600    continue
         do 700 j=5,8
            icon(j,i)=icon(j,i)-8
            if(is(icon(j,i)).le.0) stop 2720
            is(icon(j,i))=i
  700    continue
  800 continue
c
 1000 continue
      if(.not.flphis) tetra = tetra - ielm
c
      return
      end
*RUVTET
c
c This subroutine will compress data structure in order to save space
c by eliminating discarded tetrahedra while keeping artificial ones
c
      subroutine ruvtet(tetra, tetru, icon, is, ifl)
c
      integer tetra, tetru, icon(8,*), is(*), ifl(*)
      integer i, j, ii, ielm
c
c     identify true tetrahedra
c
      tetru = 0
      ielm = 0
      do 100 i = 1, tetra
         if ((icon(5,i).le.0) .or. (icon(6,i).le.0) .or.
     *       (icon(7,i).le.0) .or. (icon(8,i).le.0)) then
            ielm = ielm + 1
            ifl(i) = 0
         elseif ((icon(5,i).le.8) .or. (icon(6,i).le.8) .or.
     *           (icon(7,i).le.8) .or. (icon(8,i).le.8)) then
            ifl(i) = 1
         else
            tetru=tetru+1
            ifl(i) = 1
         endif
  100 continue
c
c     compress icon
c
      ii=0
      do 500 i=1,tetra
         if(ifl(i).eq.0) go to 500
         ii=ii+1
         ifl(i)=ii
         do 400 j=1,8
            icon(j,ii)=icon(j,i)
  400    continue
  500 continue
c
      tetra = tetra - ielm
c
c     update icon for tetrahedra and is for vertices
c
      do 800 i=1,tetra
         do 600 j=1,4
            if(icon(j,i).eq.0) go to 600
            icon(j,i)=ifl(icon(j,i))
  600    continue
         do 700 j=5,8
            if(is(icon(j,i)).le.0) stop 2730
            is(icon(j,i))=i
  700    continue
  800 continue
c
      return
      end
*CONSIS
c
c     subroutine consis to -
c
c     test consistency of diagram
c
c     May 1, 1989
c
      subroutine consis(icon, is, ifl, n, ivnxt)
c
      integer icon(8,*), is(*), ifl(*), ikon(8,1)
      integer site0, site1, site2, site3, n, ivnxt
      integer i, iscur, isone, islst, isini, indx
c
c     test initial tetrahedron for each site
c
      do 50 i = 1, n
          iscur = is(i)
          if (iscur .le. 0) goto 50
          if(icon(5,iscur) .ne. i .and. icon(6,iscur) .ne. i .and.
     *       icon(7,iscur) .ne. i .and. icon(8,iscur) .ne. i) stop 2810
   50 continue
c
c     initialize
c
      isone = 1
      do 60 i = 1, n
         if(is(i) .gt. 0) go to 80
   60 continue
      stop 2820
   80 continue
      islst = is(i)
      isini = islst
c
      do 100 i = 1, ivnxt
          ifl(i) = 0
  100 continue
c
      ifl(isini) = 1
      indx = 1
      iscur = icon(1,isini)
      if(iscur.eq.0) go to 500
      site0 = icon(5,isini)
      site1 = icon(6,isini)
      site2 = icon(7,isini)
      site3 = icon(8,isini)
c
c     reorder iscur relative to site1 and site2, and test
c
  200 continue
      if(site0.eq.site1 .or. site0.eq.site2 .or. site0.eq.site3 .or.
     *   site1.eq.site2 .or. site1.eq.site3 .or. site2.eq.site3)
     *   stop 2830
      call reordr(icon, site1, site2, iscur)
      if(icon(7,iscur) .ne. site3) stop 2840
      if(icon(4,iscur) .ne. islst) stop 2850
      if(icon(8,iscur) .eq. site0) stop 2855
      ifl(iscur) = 1
c
c     obtain next tetrahedron
c
      islst = iscur
      indx = 1
      iscur = icon(1,islst)
      if(iscur.eq.0) go to 500
      site0 = icon(5,islst)
      site1 = icon(6,islst)
      site2 = icon(7,islst)
      site3 = icon(8,islst)
      if(ifl(iscur) .ne. 1) go to 200
c
c     reorder iscur relative to site1 and site2, and test
c
  300 continue
      if(site0.eq.site1 .or. site0.eq.site2 .or. site0.eq.site3 .or.
     *   site1.eq.site2 .or. site1.eq.site3 .or. site2.eq.site3)
     *   stop 2860
      do 400 i = 1, 8
          ikon(i,1) = icon(i,iscur)
  400 continue
      call reordr(ikon, site1, site2, isone)
      if(ikon(7,1) .ne. site3) stop 2865
      if(ikon(4,1) .ne. islst) stop 2870
      if(ikon(8,1) .eq. site0) stop 2875
c
c     obtain next tetrahedron
c
  500 continue
      if(indx.eq.1) then
          indx = 2
          iscur = icon(2,islst)
          if(iscur.eq.0) go to 500
          site0 = icon(6,islst)
          site1 = icon(5,islst)
          site2 = icon(8,islst)
          site3 = icon(7,islst)
          if(ifl(iscur) .ne. 1) go to 200
          go to 300
      elseif(indx.eq.2) then
          indx = 3
          iscur = icon(3,islst)
          if(iscur.eq.0) go to 500
          site0 = icon(7,islst)
          site1 = icon(5,islst)
          site2 = icon(6,islst)
          site3 = icon(8,islst)
          if(ifl(iscur) .ne. 1) go to 200
          go to 300
      elseif(indx.eq.3) then
          if(islst .ne. isini) then
              iscur = islst
              islst = icon(4,iscur)
              if(islst. eq. 0) stop 2880
              if(icon(1,islst) .eq. iscur) then
                  indx = 1
              elseif(icon(2,islst) .eq. iscur) then
                  indx = 2
              elseif(icon(3,islst) .eq. iscur) then
                  indx = 3
              elseif(icon(4,islst) .eq. iscur) then
                  indx = 4
              else
                  stop 2885
              endif
              go to 500
          else
              indx = 4
              iscur = icon(4,islst)
              if(iscur.eq.0) go to 500
              site0 = icon(8,islst)
              site1 = icon(5,islst)
              site2 = icon(7,islst)
              site3 = icon(6,islst)
              if(ifl(iscur) .ne. 1) go to 200
              go to 300
          endif
      endif
      if(islst .ne. isini) stop 2890
c
c     write (*,*) ' '
c     write (*,*) '**************************************'
c     write (*,*) 'consistency check satisfied'
c     write (*,*) '**************************************'
c     write (*,*) ' '
c
      return
      end
*ORIENT
c
c     This subroutine will test the orientation of the tetrahedra
c
      subroutine orient(tetra, icon, ifl, xi, yi, zi, x, y, z, x2, y2,
     *                  z2, idmin, mhalf, mfull, isclp, epz)
c
      double precision xi(*), yi(*), zi(*)
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*)
      double precision epz
      integer tetra, icon(8,*), ifl(*), a, b, c, d, idmin
      integer isclp(*), mhalf, mfull
      integer i, iside
c
c     test all tetrahedra with ifl equal to 1
c
      idmin = 0
      do 200 i=1,tetra
         if(ifl(i).ne.1) go to 200
         a=icon(5,i)
         b=icon(6,i)
         c=icon(7,i)
         d=icon(8,i)
         call irsign(xi, yi, zi, x, y, z, x2, y2, z2, d, a, b, c,
     *               mhalf, mfull, isclp, epz, iside)
         if(iside .le. 0) idmin = idmin+1
  200 continue
c
      return
      end
*REORDR
c
c     subroutine reordr to -
c
c     reorder icon(i,iscur), i = 1, ..., 8, so that site1 equals
c     icon(5,iscur) and site2 equals icon(6,iscur)
c
c     July 22, 1988
c
      subroutine reordr(icon, site1, site2, iscur)
c
      integer icon(8,*), site1, site2, iscur, itemp
c
      if(icon(5,iscur) .eq. site1) go to 200
      if(icon(6,iscur) .eq. site1) then
          itemp = icon(1,iscur)
          icon(1,iscur) = icon(2,iscur)
          icon(2,iscur) = icon(4,iscur)
          icon(4,iscur) = itemp
          itemp = icon(5,iscur)
          icon(5,iscur) = icon(6,iscur)
          icon(6,iscur) = icon(8,iscur)
          icon(8,iscur) = itemp
      elseif(icon(7,iscur) .eq. site1) then
          itemp = icon(1,iscur)
          icon(1,iscur) = icon(3,iscur)
          icon(3,iscur) = icon(2,iscur)
          icon(2,iscur) = itemp
          itemp = icon(5,iscur)
          icon(5,iscur) = icon(7,iscur)
          icon(7,iscur) = icon(6,iscur)
          icon(6,iscur) = itemp
      elseif(icon(8,iscur) .eq. site1) then
          itemp = icon(1,iscur)
          icon(1,iscur) = icon(4,iscur)
          icon(4,iscur) = icon(3,iscur)
          icon(3,iscur) = itemp
          itemp = icon(5,iscur)
          icon(5,iscur) = icon(8,iscur)
          icon(8,iscur) = icon(7,iscur)
          icon(7,iscur) = itemp
      else
          stop 2910
      endif
  200 continue
c
      if(icon(6,iscur) .eq. site2) go to 300
      if(icon(7,iscur) .eq. site2) then
          itemp = icon(2,iscur)
          icon(2,iscur) = icon(3,iscur)
          icon(3,iscur) = icon(4,iscur)
          icon(4,iscur) = itemp
          itemp = icon(6,iscur)
          icon(6,iscur) = icon(7,iscur)
          icon(7,iscur) = icon(8,iscur)
          icon(8,iscur) = itemp
      elseif(icon(8,iscur) .eq. site2) then
          itemp = icon(2,iscur)
          icon(2,iscur) = icon(4,iscur)
          icon(4,iscur) = icon(3,iscur)
          icon(3,iscur) = itemp
          itemp = icon(6,iscur)
          icon(6,iscur) = icon(8,iscur)
          icon(8,iscur) = icon(7,iscur)
          icon(7,iscur) = itemp
      else
          stop 2920
      endif
  300 continue
c
      return
      end
*SITORD
c
c     subroutine sitord to -
c
c     reorder icon(i,iscur), i = 1, ..., 8, so that site1 equals
c     icon(5,iscur) 
c
c     July 22, 1988
c
      subroutine sitord(icon, site1, iscur)
c
      integer icon(8,*), site1, iscur, itemp
c
      if(icon(5,iscur) .eq. site1) go to 200
      if(icon(6,iscur) .eq. site1) then
          itemp = icon(1,iscur)
          icon(1,iscur) = icon(2,iscur)
          icon(2,iscur) = icon(4,iscur)
          icon(4,iscur) = itemp
          itemp = icon(5,iscur)
          icon(5,iscur) = icon(6,iscur)
          icon(6,iscur) = icon(8,iscur)
          icon(8,iscur) = itemp
      elseif(icon(7,iscur) .eq. site1) then
          itemp = icon(1,iscur)
          icon(1,iscur) = icon(3,iscur)
          icon(3,iscur) = icon(2,iscur)
          icon(2,iscur) = itemp
          itemp = icon(5,iscur)
          icon(5,iscur) = icon(7,iscur)
          icon(7,iscur) = icon(6,iscur)
          icon(6,iscur) = itemp
      elseif(icon(8,iscur) .eq. site1) then
          itemp = icon(1,iscur)
          icon(1,iscur) = icon(4,iscur)
          icon(4,iscur) = icon(3,iscur)
          icon(3,iscur) = itemp
          itemp = icon(5,iscur)
          icon(5,iscur) = icon(8,iscur)
          icon(8,iscur) = icon(7,iscur)
          icon(7,iscur) = itemp
      else
          stop 3010
      endif
  200 continue
      return
      end
*VRTORD
c
c This routine will order vertices a, b, c, d of a tetrahedron
c so that a>b, b>c, b>d. Data structure is updated
c
      subroutine vrtord(icon, curr, a, b, c, d)
c
      integer icon(8,*), a, b, c, d, curr, it
c
      if(a.lt.b)then
         it=a
         a=b
         b=it
      endif
      if(a.lt.c)then
         it=a
         a=c
         c=it
      endif
      if(a.lt.d)then
         it=a
         a=d
         d=it
      endif
      if(b.lt.c) b=c
      if(b.lt.d) b=d
      call reordr(icon,a,b,curr)
      c = icon(7,curr)
      d = icon(8,curr)
      if(b.gt.a.or.c.gt.b.or.d.gt.b) stop 3110
c
      return
      end
*VRTARR
c
c This routine will arrange vertices b, c, d of a tetrahedron
c so that b>c, b>d. Data structure is not updated
c
      subroutine vrtarr(i2,i3,i4,b,c,d)
c
      integer b, c, d, i2, i3, i4, ix
c
      b=i2
      c=i3
      d=i4
      ix = max0(b,c,d)
      if(b .eq. ix) go to 100
      if(c .eq. ix) then
         ix = b
         b = c
         c = d
         d = ix
      else
         ix = b
         b = d
         d = c
         c = ix
      endif
  100 continue
      if(c.gt.b .or. d.gt.b) stop 3210
c
      return
      end
*RDMORD
c
c     subroutine to reorder randomly integers from 1 to n
c
c
      subroutine rdmord(wr, io, n, isu, jsu, ksu, nsu)
c
      real wr(*), rumi
      integer io(*), n, isu, jsu, ksu, nsu, i
      integer iu, ju, ku, nu
c
c     initialize
c
c     isu = 521288629
c     jsu = 362436069
c     ksu = 16163801
c     nsu = 131199299
c
      call mzrans(isu,jsu,ksu,nsu,iu,ju,ku,nu)
c
c     get numbers
c
      do 10 i = 1, n
          call umi(rumi,iu,ju,ku,nu)
          wr(i) = rumi
          io(i) = i
   10 continue
c
c     sort in increasing order thus obtaining random order
c     of integers from 1 to n
c
c     OPEN(21,FILE='umi.dat')
c     WRITE(21,*)(WR(I),I=1,N)
c
      call trsort(wr, io, n)
c
      return
      end
*TRSORT
c
c     subroutine trsort to -
c
c     sort an array of real numbers in increasing order
c     in O(k log k) time
c
c     January 15, 1988
c
      subroutine trsort(var, ji, klen)
c
      real var(*)
      integer ji(*), klen, iast, k, i, jj, iis
c
c     create initial tree in decreasing order
c
      iast = klen
      do 100 k = 1, klen
          i = k
   50     continue
c
c     check if current node is as small as father
c
          if(i .eq. 1) go to 100
          if(var(ji(i)) .le. var(ji(i/2))) go to 100
          jj = ji(i)
          ji(i) = ji(i/2)
          ji(i/2) = jj
          i = i/2
          go to 50
  100 continue
      if(iast .eq. 1) go to 160
c
c     sort by shrinking tree: last element is moved to the
c     first position
c
  102 continue
      i = 1
      jj = ji(1)
      ji(1) = ji(iast)
      ji(iast) = jj
      if(iast .eq. 2) go to 160
      iast = iast - 1
c
  105 continue
      iis = 2*i
c
c     check which sons exist
c
      if(iis - iast) 110, 140, 150
c
c     both sons exist
c
  110 continue
      if(var(ji(i)) .lt. var(ji(iis))) go to 120
      if(var(ji(i)) .ge. var(ji(iis+1))) go to 150
      go to 125
c
c     check which son to be switched
c
  120 continue
      if(var(ji(iis)) .ge. var(ji(iis+1))) go to 130
c
c     adjust to switch with right son
c
  125 continue
      iis = iis + 1
c
c     switch
c
  130 continue
      jj = ji(i)
      ji(i) = ji(iis)
      ji(iis) = jj
      i = iis
      go to 105
c
c     only left son exists
c
  140 continue
      if(var(ji(i)) .ge. var(ji(iis))) go to 150
      jj = ji(i)
      ji(i) = ji(iis)
      ji(iis) = jj
c
c     no more switching needed
c
  150 continue
      go to 102
c
c     sorting is finished
c
  160 continue
c
      return
      end
*UMI
c
c     subroutine umi to -
c
c     generate random numbers
c
      subroutine umi(rumi,i,j,k,n)
c
      real rumi
      integer mzrn,i,j,k,n
c
      call mzran(i,j,k,n,mzrn)
      rumi = .5 + .2328306e-9*mzrn
c
      return
      end
*MZRAN
c
c     subroutine mzran to -
c
c     do computations in order to generate random numbers
c
      subroutine mzran(i,j,k,n,mzrn)
c
      integer i, j, k, n, mzrn
c
      mzrn = i - k
      if(mzrn .lt. 0) mzrn = mzrn + 2147483579
      i = j
      j = k
      k = mzrn
      n = 69069*n + 1013904243
      mzrn = mzrn + n
c
      return
      end
*MZRANS
c
c     subroutine mzrans to -
c
c     initialize in order to generate random numbers
c
      subroutine mzrans(is,js,ks,ns,i,j,k,n)
c
c     save is,js,ks,ns
c     data is,js,ks,ns/521288629,362436069,16163801,1131199299/
c
      integer i, j, k, n
      integer is, js, ks, ns
c
      i = 1+iabs(is)
      j = 1+iabs(js)
      k = 1+iabs(ks)
      n = ns
c
      return
      end
*DSTNCE
c
c This subroutine will compute the distance from a point to a facet of
c a tetrahedron.
c
      subroutine dstnce(x, y, z, p, q, r, epz, k, dist, ipossi)
c
      integer p, q, r, k, ipossi
      double precision x(*), y(*), z(*)
      double precision epz, dist
      double precision xvec1, yvec1, zvec1, xvec2, yvec2, zvec2
      double precision xvec3, yvec3, zvec3, dst1, dst2, dst3
      double precision dotx, doty, dotz, dmax, dlun
      double precision xvecp, yvecp, zvecp, dstp, dlen
c
      ipossi = 0
      xvec1 = x(q) - x(p)
      yvec1 = y(q) - y(p)
      zvec1 = z(q) - z(p)
      xvec2 = x(r) - x(p)
      yvec2 = y(r) - y(p)
      zvec2 = z(r) - z(p)
      xvec3 = x(q) - x(r)
      yvec3 = y(q) - y(r)
      zvec3 = z(q) - z(r)
      dst1=dsqrt(xvec1**2+yvec1**2+zvec1**2)
      dst2=dsqrt(xvec2**2+yvec2**2+zvec2**2)
      dst3=dsqrt(xvec3**2+yvec3**2+zvec3**2)
      if(dst1.lt.epz .or. dst2.lt.epz .or. dst3.lt.epz) then
         ipossi = 1
         go to 1000
      endif
      dmax = dmax1(dst1,dst2,dst3)
c
      dotx = yvec1 * zvec2 - yvec2 * zvec1
      doty = - xvec1 * zvec2 + xvec2 * zvec1
      dotz = xvec1 * yvec2 - xvec2 * yvec1
      dlen = dsqrt (dotx**2 + doty**2 + dotz**2)
      if(dlen.lt.epz .or. dlen/dmax.lt.epz)then
         ipossi = 1
         go to 1000
      endif
c
      xvecp = x(k) - x(p)
      yvecp = y(k) - y(p)
      zvecp = z(k) - z(p)
      dstp=dsqrt(xvecp**2+yvecp**2+zvecp**2)
      if(dstp.lt.epz) then
         ipossi = 1
         go to 1000
      endif
c
      dlun=dstp*dmax
      dlun=dmax1(dlen,dlun)
      dist=(xvecp*dotx+yvecp*doty+zvecp*dotz)/dlun
      if(dist.gt.-epz .and. dist.lt.epz)then
         ipossi = 1
      endif
c
 1000 continue
      return
      end
*CTRAD
c
c This subroutine will compute the orthogonal center of a tetrahedron
c
      subroutine ctrad(x, y, z, w, xctr, yctr, zctr, a, b, c, d,
     *                 epz, delaun, ipossi)
c
      double precision x(*), y(*), z(*), w(*)
      double precision epz, xctr, yctr, zctr
      double precision xm, ym, zm, xn, yn, zn, xu, yu, zu, xv, yv, zv
      double precision xw, yw, zw, xq, yq, zq, xe, ye, ze, xl, yl, zl
      double precision xt, yt, zt
      double precision norm, lambda, normu, normv, denom, dmax
      integer a, b, c, d, ipossi
      logical delaun
c
c     initialize
c
      ipossi = 0
c
c     find midpoints of edges ac and ab
c
      xm = (x(a) + x(c)) / 2.0d0
      ym = (y(a) + y(c)) / 2.0d0
      zm = (z(a) + z(c)) / 2.0d0
c
      xn = (x(a) + x(b)) / 2.0d0
      yn = (y(a) + y(b)) / 2.0d0
      zn = (z(a) + z(b)) / 2.0d0
c
c     compute edge vectors u and v for edges ac and ab
c
      xu = x(c) - x(a)
      yu = y(c) - y(a)
      zu = z(c) - z(a)
c
      xv = x(b) - x(a)
      yv = y(b) - y(a)
      zv = z(b) - z(a)
c
c     compute lengths of u and v
c
      normu = dsqrt(xu**2 + yu**2 + zu**2)
      normv = dsqrt(xv**2 + yv**2 + zv**2)
      if(normu.lt.epz .or. normv.lt.epz) then
         ipossi = 1
         go to 1000
      endif
      dmax = dmax1(normu,normv)
c
c     find perpendicular to facet abc of tetrahedron
c
      xw = yu * zv - zu * yv
      yw = -xu * zv + zu * xv
      zw = xu * yv - yu * xv
c
c     test whether edges ac, ab are colinear
c
      norm = dsqrt(xw**2 + yw**2 + zw**2)/dmax
      if(norm .lt. epz)then
         ipossi = 1
         go to 1000
      endif
      xw = xw/normu
      yw = yw/normu
      zw = zw/normu
c
c     normalize u and v
c
      xu = xu / normu
      yu = yu / normu
      zu = zu / normu
      xv = xv / normv
      yv = yv / normv
      zv = zv / normv
c
c     compute orthogonal center of edge ac
c
      if(.not.delaun)then
         lambda = ((w(a)-w(c))/normu)/2.0d0
         xm = xm + lambda*xu
         ym = ym + lambda*yu
         zm = zm + lambda*zu
      endif
c
c     compute orthogonal center of edge ab
c
      if(.not.delaun)then
         lambda = ((w(a)-w(b))/normv)/2.0d0
         xn = xn + lambda*xv
         yn = yn + lambda*yv
         zn = zn + lambda*zv
      endif
c
c     find perpendicular to edge v in plane that contains facet abc
c
      xq = yw * zv - zw * yv
      yq = -xw * zv + zw * xv
      zq = xw * yv - yw * xv
      norm = dsqrt(xq**2 + yq**2 + zq**2)
      if(norm.lt.epz) then
         ipossi = 1
         go to 1000
      endif
c
c     compute orthogonal center of facet abc
c
      denom = xu*xq + yu*yq + zu*zq
      if(denom .gt. -epz .and. denom .lt. epz) then
         ipossi = 1
         go to 1000
      endif
      lambda = (xu*(xm-xn) + yu*(ym-yn) + zu*(zm-zn)) / denom
c
      xe = xn + lambda*xq
      ye = yn + lambda*yq
      ze = zn + lambda*zq
c
c     compute edge vector t for edge ad
c
      xl = (x(a) + x(d)) / 2.0d0
      yl = (y(a) + y(d)) / 2.0d0
      zl = (z(a) + z(d)) / 2.0d0
c
      xt = x(d) - x(a)
      yt = y(d) - y(a)
      zt = z(d) - z(a)
      norm = dsqrt(xt**2 + yt**2 + zt**2)
      if(norm .lt. epz) then
         ipossi = 1
         go to 1000
      endif
      xt = xt / norm
      yt = yt / norm
      zt = zt / norm
c
c     compute orthogonal center of edge ad
c
      if(.not.delaun)then
         lambda = ((w(a)-w(d))/norm)/2.0d0
         xl = xl + lambda*xt
         yl = yl + lambda*yt
         zl = zl + lambda*zt
      endif
c
c     compute orthogonal center of tetrahedron
c
      denom = xt*xw + yt*yw + zt*zw
      if(denom .gt. -epz .and. denom .lt. epz) then
         ipossi = 1
         go to 1000
      endif
      lambda = (xt*(xl-xe) + yt*(yl-ye) + zt*(zl-ze)) / denom
c
      xctr = xe + lambda*xw
      yctr = ye + lambda*yw
      zctr = ze + lambda*zw
c 
 1000 continue
      return
      end
*BISPHR
c
c This subroutine will compute the distance from a point
c (xctr,yctr,zctr) to the chordale plane between two points
c opvert and ivrt
c
      subroutine bisphr(x, y, z, w, opvert, ivrt, epz,
     *                  xctr, yctr, zctr, tdist, delaun, ipossi)
c
      double precision x(*), y(*), z(*), w(*), norm
      integer opvert, ivrt, ipossi
      double precision epz, tdist, wambda, xctr, yctr, zctr, dif
      double precision xm, ym, zm, xu, yu, zu, xd, yd, zd, dmax
      double precision xu2, yu2, zu2
      logical delaun
c
c     find midpoint of edge from opvert to ivrt
c
      xm = (x(opvert) + x(ivrt)) / 2.0d0
      ym = (y(opvert) + y(ivrt)) / 2.0d0
      zm = (z(opvert) + z(ivrt)) / 2.0d0
c
c     find vector from ivrt to opvert
c
      xu = x(opvert) - x(ivrt)
      yu = y(opvert) - y(ivrt)
      zu = z(opvert) - z(ivrt)
c 
      norm = dsqrt(xu**2 + yu**2 + zu**2)
      if(norm .lt. epz) then
        ipossi = 1
        go to 1000
      endif
      xu2 = xu/norm
      yu2 = yu/norm
      zu2 = zu/norm
c
c     compute orthogonal center of edge ivrt-opvert
c
      if(.not.delaun)then
         wambda = ((w(ivrt)-w(opvert))/norm)/2.0d0
         xm = xm + wambda*xu2
         ym = ym + wambda*yu2
         zm = zm + wambda*zu2
      endif
c
c     compute distance
c
      xd = xctr - xm
      yd = yctr - ym
      zd = zctr - zm
      dif = dsqrt(xd**2 + yd**2 + zd**2)
      dmax = dmax1(norm,dif)
      tdist = (xd*xu + yd*yu + zd*zu) / dmax
      if(tdist.gt. -epz .and. tdist.lt. epz) then
         ipossi = 1
      endif
c
 1000 continue
      return
      end
*IPSIGN
c
c     subroutine for determining position of point ifou with respect
c     to plane that contains points ifir, isec, ithi
c     if positive then ifou is on positive side of plane
c     if negative then ifou is on negative side of plane
c     if zero then ifou is in plane
c
      subroutine ipsign(x, y, z, x2, y2, z2, ifir, isec, ithi,
     *                  ifou, mhalf, mfull, isclp, ipout)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*)
      integer ifir, isec, ithi, ifou
      integer isclp(*), mhalf, mfull, nkmax, ipout
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax), iw(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ix4(nkmax), iy4(nkmax), iz4(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixfiw, iyfiw, izfiw, ixsew, iysew, izsew
      integer ixthw, iythw, izthw, ixfow, iyfow, izfow
      integer ixfi2, iyfi2, izfi2, ixse2, iyse2, izse2
      integer ixth2, iyth2, izth2, ixfo2, iyfo2, izfo2
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx4, isgy4, isgz4, ikx4, iky4, ikz4
      integer isgo, isgu, isgv, isgw, iko, iku, ikv, ikw
c
      ixfiw = x(ifir)
      iyfiw = y(ifir)
      izfiw = z(ifir)
      ixsew = x(isec)
      iysew = y(isec)
      izsew = z(isec)
      ixthw = x(ithi)
      iythw = y(ithi)
      izthw = z(ithi)
      ixfow = x(ifou)
      iyfow = y(ifou)
      izfow = z(ifou)
c
      ixfi2 = x2(ifir)
      iyfi2 = y2(ifir)
      izfi2 = z2(ifir)
      ixse2 = x2(isec)
      iyse2 = y2(isec)
      izse2 = z2(isec)
      ixth2 = x2(ithi)
      iyth2 = y2(ithi)
      izth2 = z2(ithi)
      ixfo2 = x2(ifou)
      iyfo2 = y2(ifou)
      izfo2 = z2(ifou)
c
      call decmp2(ixf, isgxf, ikxf, ixfiw, ixfi2, mhalf, mfull, isclp)
      call decmp2(iyf, isgyf, ikyf, iyfiw, iyfi2, mhalf, mfull, isclp)
      call decmp2(izf, isgzf, ikzf, izfiw, izfi2, mhalf, mfull, isclp)
c
      call decmp2(io, isgo, iko, ixsew, ixse2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix2, isgo, isgxf, isgx2, iko, ikxf, ikx2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iysew, iyse2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy2, isgo, isgyf, isgy2, iko, ikyf, iky2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izsew, izse2, mhalf, mfull, isclp)
      call muldif(io, izf, iz2, isgo, isgzf, isgz2, iko, ikzf, ikz2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, ixthw, ixth2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix3, isgo, isgxf, isgx3, iko, ikxf, ikx3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iythw, iyth2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy3, isgo, isgyf, isgy3, iko, ikyf, iky3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izthw, izth2, mhalf, mfull, isclp)
      call muldif(io, izf, iz3, isgo, isgzf, isgz3, iko, ikzf, ikz3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, ixfow, ixfo2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix4, isgo, isgxf, isgx4, iko, ikxf, ikx4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iyfow, iyfo2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy4, isgo, isgyf, isgy4, iko, ikyf, iky4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izfow, izfo2, mhalf, mfull, isclp)
      call muldif(io, izf, iz4, isgo, isgzf, isgz4, iko, ikzf, ikz4,
     *              nkmax, mhalf)
c
      call mulmul(iy2, iz3, iv, isgy2, isgz3, isgv, iky2, ikz3, ikv,
     *              nkmax, mhalf)
      call mulmul(iz2, iy3, iu, isgz2, isgy3, isgu, ikz2, iky3, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, iw, isgv, isgu, isgw, ikv, iku, ikw,
     *              nkmax, mhalf)
      call mulmul(iw, ix4, io, isgw, isgx4, isgo, ikw, ikx4, iko,
     *              nkmax, mhalf)
c
      call mulmul(iz2, ix3, iv, isgz2, isgx3, isgv, ikz2, ikx3, ikv,
     *              nkmax, mhalf)
      call mulmul(ix2, iz3, iu, isgx2, isgz3, isgu, ikx2, ikz3, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, iw, isgv, isgu, isgw, ikv, iku, ikw,
     *              nkmax, mhalf)
      call mulmul(iw, iy4, iu, isgw, isgy4, isgu, ikw, iky4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(io, iu, iw, isgo, isgu, isgw, iko, iku, ikw,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iy3, iv, isgx2, isgy3, isgv, ikx2, iky3, ikv,
     *              nkmax, mhalf)
      call mulmul(iy2, ix3, iu, isgy2, isgx3, isgu, iky2, ikx3, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz4, iu, isgo, isgz4, isgu, iko, ikz4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iw, iu, io, isgw, isgu, isgo, ikw, iku, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*IPSIG1
c
      subroutine ipsig1(x, y, z, x2, y2, z2, xc, yc, zc, ifir,
     *                  isec, ithi, ifou, ifif, isix, mhalf,
     *                  mfull, isclp, ipout)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer ifir, isec, ithi, ifou, ifif, isix
      integer isclp(*), mhalf, mfull, nkmax, ipout
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ix4(nkmax), iy4(nkmax), iz4(nkmax)
      integer ix5(nkmax), iy5(nkmax), iz5(nkmax)
      integer ix6(nkmax), iy6(nkmax), iz6(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixthw, iythw, izthw, ixfow, iyfow, izfow
      integer ixfiw, iyfiw, izfiw, ixsiw, iysiw, izsiw
      integer ixfuw, iyfuw, izfuw, ixsuw, iysuw, izsuw
      integer ixth2, iyth2, izth2, ixfo2, iyfo2, izfo2
      integer ixfi2, iyfi2, izfi2, ixsi2, iysi2, izsi2
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx4, isgy4, isgz4, ikx4, iky4, ikz4
      integer isgx5, isgy5, isgz5, ikx5, iky5, ikz5
      integer isgx6, isgy6, isgz6, ikx6, iky6, ikz6
      integer isgo, isgu, isgv, iko, iku, ikv
c
      ixthw = x(ithi)
      iythw = y(ithi)
      izthw = z(ithi)
      ixfow = x(ifou)
      iyfow = y(ifou)
      izfow = z(ifou)
      ixfiw = x(ifif)
      iyfiw = y(ifif)
      izfiw = z(ifif)
      ixsiw = x(isix)
      iysiw = y(isix)
      izsiw = z(isix)
c
      ixth2 = x2(ithi)
      iyth2 = y2(ithi)
      izth2 = z2(ithi)
      ixfo2 = x2(ifou)
      iyfo2 = y2(ifou)
      izfo2 = z2(ifou)
      ixfi2 = x2(ifif)
      iyfi2 = y2(ifif)
      izfi2 = z2(ifif)
      ixsi2 = x2(isix)
      iysi2 = y2(isix)
      izsi2 = z2(isix)
c
      ixfuw = xc(ifir)
      iyfuw = yc(ifir)
      izfuw = zc(ifir)
      ixsuw = xc(isec)
      iysuw = yc(isec)
      izsuw = zc(isec)
c
      ikxf = 2
      ikyf = 2
      ikzf = 2
      call decomp(ixf, isgxf, ixfuw, mhalf)
      call decomp(iyf, isgyf, iyfuw, mhalf)
      call decomp(izf, isgzf, izfuw, mhalf)
c
      call decmp2(ix3, isgx3, ikx3, ixthw, ixth2, mhalf, mfull, isclp)
      call decmp2(iy3, isgy3, iky3, iythw, iyth2, mhalf, mfull, isclp)
      call decmp2(iz3, isgz3, ikz3, izthw, izth2, mhalf, mfull, isclp)
      call decmp2(ix5, isgx5, ikx5, ixfiw, ixfi2, mhalf, mfull, isclp)
      call decmp2(iy5, isgy5, iky5, iyfiw, iyfi2, mhalf, mfull, isclp)
      call decmp2(iz5, isgz5, ikz5, izfiw, izfi2, mhalf, mfull, isclp)
c
      iko = 2
      call decomp(io, isgo, ixsuw, mhalf)
      call muldif(io, ixf, ix2, isgo, isgxf, isgx2, iko, ikxf, ikx2,
     *              nkmax, mhalf)
      call decomp(io, isgo, iysuw, mhalf)
      call muldif(io, iyf, iy2, isgo, isgyf, isgy2, iko, ikyf, iky2,
     *              nkmax, mhalf)
      call decomp(io, isgo, izsuw, mhalf)
      call muldif(io, izf, iz2, isgo, isgzf, isgz2, iko, ikzf, ikz2,
     *              nkmax, mhalf)
c
      call decmp2(io, isgo, iko, ixfow, ixfo2, mhalf, mfull, isclp)
      call muldif(io, ix3, ix4, isgo, isgx3, isgx4, iko, ikx3, ikx4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iyfow, iyfo2, mhalf, mfull, isclp)
      call muldif(io, iy3, iy4, isgo, isgy3, isgy4, iko, iky3, iky4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izfow, izfo2, mhalf, mfull, isclp)
      call muldif(io, iz3, iz4, isgo, isgz3, isgz4, iko, ikz3, ikz4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, ixsiw, ixsi2, mhalf, mfull, isclp)
      call muldif(io, ix5, ix6, isgo, isgx5, isgx6, iko, ikx5, ikx6,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iysiw, iysi2, mhalf, mfull, isclp)
      call muldif(io, iy5, iy6, isgo, isgy5, isgy6, iko, iky5, iky6,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izsiw, izsi2, mhalf, mfull, isclp)
      call muldif(io, iz5, iz6, isgo, isgz5, isgz6, iko, ikz5, ikz6,
     *              nkmax, mhalf)
c
      call mulmul(iy2, iz4, io, isgy2, isgz4, isgo, iky2, ikz4, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix6, iv, isgo, isgx6, isgv, iko, ikx6, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz2, ix4, io, isgz2, isgx4, isgo, ikz2, ikx4, iko,
     *              nkmax, mhalf)
      call mulmul(io, iy6, iu, isgo, isgy6, isgu, iko, iky6, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iy4, iv, isgx2, isgy4, isgv, ikx2, iky4, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iz6, iu, isgv, isgz6, isgu, ikv, ikz6, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz2, iy4, io, isgz2, isgy4, isgo, ikz2, iky4, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix6, iu, isgo, isgx6, isgu, iko, ikx6, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iz4, iv, isgx2, isgz4, isgv, ikx2, ikz4, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy6, iu, isgv, isgy6, isgu, ikv, iky6, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iy2, ix4, io, isgy2, isgx4, isgo, iky2, ikx4, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz6, iu, isgo, isgz6, isgu, iko, ikz6, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*IPSIG2
c
      subroutine ipsig2(x, y, z, x2, y2, z2, xc, yc, zc, ifir,
     *                  isec, ithi, ifou, ifif, isix, mhalf,
     *                  mfull, isclp, ipout)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer ifir, isec, ithi, ifou, ifif, isix
      integer isclp(*), mhalf, mfull, nkmax, ipout
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ix4(nkmax), iy4(nkmax), iz4(nkmax)
      integer ix5(nkmax), iy5(nkmax), iz5(nkmax)
      integer ix6(nkmax), iy6(nkmax), iz6(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixthw, iythw, izthw, ixfow, iyfow, izfow
      integer ixfiw, iyfiw, izfiw, ixsiw, iysiw, izsiw
      integer ixfuw, iyfuw, izfuw, ixsuw, iysuw, izsuw
      integer ixfi2, iyfi2, izfi2, ixsi2, iysi2, izsi2
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx4, isgy4, isgz4, ikx4, iky4, ikz4
      integer isgx5, isgy5, isgz5, ikx5, iky5, ikz5
      integer isgx6, isgy6, isgz6, ikx6, iky6, ikz6
      integer isgo, isgu, isgv, iko, iku, ikv
c
      ixfiw = x(ifif)
      iyfiw = y(ifif)
      izfiw = z(ifif)
      ixsiw = x(isix)
      iysiw = y(isix)
      izsiw = z(isix)
c
      ixfi2 = x2(ifif)
      iyfi2 = y2(ifif)
      izfi2 = z2(ifif)
      ixsi2 = x2(isix)
      iysi2 = y2(isix)
      izsi2 = z2(isix)
c
      ixfuw = xc(ifir)
      iyfuw = yc(ifir)
      izfuw = zc(ifir)
      ixsuw = xc(isec)
      iysuw = yc(isec)
      izsuw = zc(isec)
      ixthw = xc(ithi)
      iythw = yc(ithi)
      izthw = zc(ithi)
      ixfow = xc(ifou)
      iyfow = yc(ifou)
      izfow = zc(ifou)
c
      ikxf = 2
      ikyf = 2
      ikzf = 2
      ikx3 = 2
      iky3 = 2
      ikz3 = 2
      call decomp(ixf, isgxf, ixfuw, mhalf)
      call decomp(iyf, isgyf, iyfuw, mhalf)
      call decomp(izf, isgzf, izfuw, mhalf)
      call decomp(ix3, isgx3, ixthw, mhalf)
      call decomp(iy3, isgy3, iythw, mhalf)
      call decomp(iz3, isgz3, izthw, mhalf)
c
      call decmp2(ix5, isgx5, ikx5, ixfiw, ixfi2, mhalf, mfull, isclp)
      call decmp2(iy5, isgy5, iky5, iyfiw, iyfi2, mhalf, mfull, isclp)
      call decmp2(iz5, isgz5, ikz5, izfiw, izfi2, mhalf, mfull, isclp)
c
      iko = 2
      call decomp(io, isgo, ixsuw, mhalf)
      call muldif(io, ixf, ix2, isgo, isgxf, isgx2, iko, ikxf, ikx2,
     *              nkmax, mhalf)
      call decomp(io, isgo, iysuw, mhalf)
      call muldif(io, iyf, iy2, isgo, isgyf, isgy2, iko, ikyf, iky2,
     *              nkmax, mhalf)
      call decomp(io, isgo, izsuw, mhalf)
      call muldif(io, izf, iz2, isgo, isgzf, isgz2, iko, ikzf, ikz2,
     *              nkmax, mhalf)
      call decomp(io, isgo, ixfow, mhalf)
      call muldif(io, ix3, ix4, isgo, isgx3, isgx4, iko, ikx3, ikx4,
     *              nkmax, mhalf)
      call decomp(io, isgo, iyfow, mhalf)
      call muldif(io, iy3, iy4, isgo, isgy3, isgy4, iko, iky3, iky4,
     *              nkmax, mhalf)
      call decomp(io, isgo, izfow, mhalf)
      call muldif(io, iz3, iz4, isgo, isgz3, isgz4, iko, ikz3, ikz4,
     *              nkmax, mhalf)
c
      call decmp2(io, isgo, iko, ixsiw, ixsi2, mhalf, mfull, isclp)
      call muldif(io, ix5, ix6, isgo, isgx5, isgx6, iko, ikx5, ikx6,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iysiw, iysi2, mhalf, mfull, isclp)
      call muldif(io, iy5, iy6, isgo, isgy5, isgy6, iko, iky5, iky6,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izsiw, izsi2, mhalf, mfull, isclp)
      call muldif(io, iz5, iz6, isgo, isgz5, isgz6, iko, ikz5, ikz6,
     *              nkmax, mhalf)
c
      call mulmul(iy2, iz4, io, isgy2, isgz4, isgo, iky2, ikz4, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix6, iv, isgo, isgx6, isgv, iko, ikx6, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz2, ix4, io, isgz2, isgx4, isgo, ikz2, ikx4, iko,
     *              nkmax, mhalf)
      call mulmul(io, iy6, iu, isgo, isgy6, isgu, iko, iky6, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iy4, iv, isgx2, isgy4, isgv, ikx2, iky4, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iz6, iu, isgv, isgz6, isgu, ikv, ikz6, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz2, iy4, io, isgz2, isgy4, isgo, ikz2, iky4, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix6, iu, isgo, isgx6, isgu, iko, ikx6, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iz4, iv, isgx2, isgz4, isgv, ikx2, ikz4, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy6, iu, isgv, isgy6, isgu, ikv, iky6, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iy2, ix4, io, isgy2, isgx4, isgo, iky2, ikx4, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz6, iu, isgo, isgz6, isgu, iko, ikz6, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*IPSIG3
c
      subroutine ipsig3(x, y, z, x2, y2, z2, xc, yc, zc, ifir, isec,
     *                  ithi, ifou, mhalf, mfull, isclp, ifn, ipout)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer ifir, isec, ithi, ifou
      integer isclp(*), mhalf, mfull, ifn, nkmax, ipout
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ix4(nkmax), iy4(nkmax), iz4(nkmax)
      integer ixfiw, iyfiw, izfiw, ixsew, iysew, izsew
      integer ixthw, iythw, izthw, ixfow, iyfow, izfow
      integer ixfi2, iyfi2, izfi2, ixse2, iyse2, izse2
      integer ixth2, iyth2, izth2
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx4, isgy4, isgz4, ikx4, iky4, ikz4
      integer isgo, isgu, isgv, iko, iku, ikv
c
      ixfiw = x(ifir)
      iyfiw = y(ifir)
      izfiw = z(ifir)
      ixsew = x(isec)
      iysew = y(isec)
      izsew = z(isec)
      ixthw = x(ithi)
      iythw = y(ithi)
      izthw = z(ithi)
c
      ixfi2 = x2(ifir)
      iyfi2 = y2(ifir)
      izfi2 = z2(ifir)
      ixse2 = x2(isec)
      iyse2 = y2(isec)
      izse2 = z2(isec)
      ixth2 = x2(ithi)
      iyth2 = y2(ithi)
      izth2 = z2(ithi)
c
      ixfow = xc(ifou)
      iyfow = yc(ifou)
      izfow = zc(ifou)
c
      call decmp2(ixf, isgxf, ikxf, ixfiw, ixfi2, mhalf, mfull, isclp)
      call decmp2(iyf, isgyf, ikyf, iyfiw, iyfi2, mhalf, mfull, isclp)
      call decmp2(izf, isgzf, ikzf, izfiw, izfi2, mhalf, mfull, isclp)
c
      ikx4 = 2
      iky4 = 2
      ikz4 = 2
      call decomp(ix4, isgx4, ixfow, mhalf)
      call decomp(iy4, isgy4, iyfow, mhalf)
      call decomp(iz4, isgz4, izfow, mhalf)
c
      call decmp2(io, isgo, iko, ixsew, ixse2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix2, isgo, isgxf, isgx2, iko, ikxf, ikx2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iysew, iyse2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy2, isgo, isgyf, isgy2, iko, ikyf, iky2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izsew, izse2, mhalf, mfull, isclp)
      call muldif(io, izf, iz2, isgo, isgzf, isgz2, iko, ikzf, ikz2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, ixthw, ixth2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix3, isgo, isgxf, isgx3, iko, ikxf, ikx3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iythw, iyth2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy3, isgo, isgyf, isgy3, iko, ikyf, iky3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izthw, izth2, mhalf, mfull, isclp)
      call muldif(io, izf, iz3, isgo, isgzf, isgz3, iko, ikzf, ikz3,
     *              nkmax, mhalf)
c
      call mulmul(iy4, iz2, io, isgy4, isgz2, isgo, iky4, ikz2, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix3, iv, isgo, isgx3, isgv, iko, ikx3, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz4, ix2, io, isgz4, isgx2, isgo, ikz4, ikx2, iko,
     *              nkmax, mhalf)
      call mulmul(io, iy3, iu, isgo, isgy3, isgu, iko, iky3, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix4, iy2, iv, isgx4, isgy2, isgv, ikx4, iky2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iz3, iu, isgv, isgz3, isgu, ikv, ikz3, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz4, iy2, io, isgz4, isgy2, isgo, ikz4, iky2, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix3, iu, isgo, isgx3, isgu, iko, ikx3, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix4, iz2, iv, isgx4, isgz2, isgv, ikx4, ikz2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy3, iu, isgv, isgy3, isgu, ikv, iky3, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iy4, ix2, io, isgy4, isgx2, isgo, iky4, ikx2, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz3, iu, isgo, isgz3, isgu, iko, ikz3, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
      if(ifn.eq.1) ipout = -isgo
c
      return
      end
*IPSIG4
c
      subroutine ipsig4(x, y, z, x2, y2, z2, xc, yc, zc, ifir, isec,
     *                  ithi, ifou, mhalf, mfull, isclp, ipout)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer ifir, isec, ithi, ifou
      integer isclp(*), mhalf, mfull, nkmax, ipout
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ix4(nkmax), iy4(nkmax), iz4(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixthw, iythw, izthw, ixfow, iyfow, izfow
      integer ixfiw, iyfiw, izfiw, ixsew, iysew, izsew
      integer ixth2, iyth2, izth2, ixfo2, iyfo2, izfo2
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx4, isgy4, isgz4, ikx4, iky4, ikz4
      integer isgo, isgu, isgv, iko, iku, ikv
c
      ixthw = x(ithi)
      iythw = y(ithi)
      izthw = z(ithi)
      ixfow = x(ifou)
      iyfow = y(ifou)
      izfow = z(ifou)
c
      ixth2 = x2(ithi)
      iyth2 = y2(ithi)
      izth2 = z2(ithi)
      ixfo2 = x2(ifou)
      iyfo2 = y2(ifou)
      izfo2 = z2(ifou)
c
      ixfiw = xc(ifir)
      iyfiw = yc(ifir)
      izfiw = zc(ifir)
      ixsew = xc(isec)
      iysew = yc(isec)
      izsew = zc(isec)
c
      call decmp2(ixf, isgxf, ikxf, ixfow, ixfo2, mhalf, mfull, isclp)
      call decmp2(iyf, isgyf, ikyf, iyfow, iyfo2, mhalf, mfull, isclp)
      call decmp2(izf, isgzf, ikzf, izfow, izfo2, mhalf, mfull, isclp)
c
      ikx2 = 2
      iky2 = 2
      ikz2 = 2
      ikx3 = 2
      iky3 = 2
      ikz3 = 2
      call decomp(ix2, isgx2, ixfiw, mhalf)
      call decomp(iy2, isgy2, iyfiw, mhalf)
      call decomp(iz2, isgz2, izfiw, mhalf)
      call decomp(ix3, isgx3, ixsew, mhalf)
      call decomp(iy3, isgy3, iysew, mhalf)
      call decomp(iz3, isgz3, izsew, mhalf)
c
      call decmp2(io, isgo, iko, ixthw, ixth2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix4, isgo, isgxf, isgx4, iko, ikxf, ikx4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iythw, iyth2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy4, isgo, isgyf, isgy4, iko, ikyf, iky4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izthw, izth2, mhalf, mfull, isclp)
      call muldif(io, izf, iz4, isgo, isgzf, isgz4, iko, ikzf, ikz4,
     *              nkmax, mhalf)
c
      call mulmul(iy2, iz3, io, isgy2, isgz3, isgo, iky2, ikz3, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix4, iv, isgo, isgx4, isgv, iko, ikx4, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz2, ix3, io, isgz2, isgx3, isgo, ikz2, ikx3, iko,
     *              nkmax, mhalf)
      call mulmul(io, iy4, iu, isgo, isgy4, isgu, iko, iky4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iy3, iv, isgx2, isgy3, isgv, ikx2, iky3, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iz4, iu, isgv, isgz4, isgu, ikv, ikz4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz2, iy3, io, isgz2, isgy3, isgo, ikz2, iky3, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix4, iu, isgo, isgx4, isgu, iko, ikx4, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iz3, iv, isgx2, isgz3, isgv, ikx2, ikz3, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy4, iu, isgv, isgy4, isgu, ikv, iky4, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iy2, ix3, io, isgy2, isgx3, isgo, iky2, ikx3, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz4, iu, isgo, isgz4, isgu, iko, ikz4, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*IPSIG6
c
      subroutine ipsig6(x, y, z, x2, y2, z2, xc, yc, zc, ifir, isec,
     *                  ithi, ifou, mhalf, mfull, isclp, ipout)
c
      integer x(*), y(*), z(*), x2(*), y2(*), z2(*), xc(*), yc(*), zc(*)
      integer ifir, isec, ithi, ifou
      integer isclp(*), mhalf, mfull, nkmax, ipout
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ix4(nkmax), iy4(nkmax), iz4(nkmax)
      integer ix5(nkmax), iy5(nkmax), iz5(nkmax)
      integer ix6(nkmax), iy6(nkmax), iz6(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixthw, iythw, izthw, ixfow, iyfow, izfow
      integer ixfiw, iyfiw, izfiw, ixsew, iysew, izsew
      integer ixfuw, iyfuw, izfuw, ixsuw, iysuw, izsuw
      integer ixth2, iyth2, izth2, ixfo2, iyfo2, izfo2
      integer ixfi2, iyfi2, izfi2, ixse2, iyse2, izse2
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx4, isgy4, isgz4, ikx4, iky4, ikz4
      integer isgx5, isgy5, isgz5, ikx5, iky5, ikz5
      integer isgx6, isgy6, isgz6, ikx6, iky6, ikz6
      integer isgo, isgu, isgv, iko, iku, ikv
c
      ixfiw = x(ifir)
      iyfiw = y(ifir)
      izfiw = z(ifir)
      ixsew = x(isec)
      iysew = y(isec)
      izsew = z(isec)
      ixthw = x(ithi)
      iythw = y(ithi)
      izthw = z(ithi)
      ixfow = x(ifou)
      iyfow = y(ifou)
      izfow = z(ifou)
c
      ixfi2 = x2(ifir)
      iyfi2 = y2(ifir)
      izfi2 = z2(ifir)
      ixse2 = x2(isec)
      iyse2 = y2(isec)
      izse2 = z2(isec)
      ixth2 = x2(ithi)
      iyth2 = y2(ithi)
      izth2 = z2(ithi)
      ixfo2 = x2(ifou)
      iyfo2 = y2(ifou)
      izfo2 = z2(ifou)
c
      ixfuw = xc(ifir)
      iyfuw = yc(ifir)
      izfuw = zc(ifir)
      ixsuw = xc(isec)
      iysuw = yc(isec)
      izsuw = zc(isec)
c
      call decmp2(ixf, isgxf, ikxf, ixfow, ixfo2, mhalf, mfull, isclp)
      call decmp2(iyf, isgyf, ikyf, iyfow, iyfo2, mhalf, mfull, isclp)
      call decmp2(izf, isgzf, ikzf, izfow, izfo2, mhalf, mfull, isclp)
c
      ikx5 = 2
      iky5 = 2
      ikz5 = 2
      ikx6 = 2
      iky6 = 2
      ikz6 = 2
      call decomp(ix5, isgx5, ixfuw, mhalf)
      call decomp(iy5, isgy5, iyfuw, mhalf)
      call decomp(iz5, isgz5, izfuw, mhalf)
      call decomp(ix6, isgx6, ixsuw, mhalf)
      call decomp(iy6, isgy6, iysuw, mhalf)
      call decomp(iz6, isgz6, izsuw, mhalf)
c
      call decmp2(io, isgo, iko, ixfiw, ixfi2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix2, isgo, isgxf, isgx2, iko, ikxf, ikx2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iyfiw, iyfi2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy2, isgo, isgyf, isgy2, iko, ikyf, iky2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izfiw, izfi2, mhalf, mfull, isclp)
      call muldif(io, izf, iz2, isgo, isgzf, isgz2, iko, ikzf, ikz2,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, ixsew, ixse2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix3, isgo, isgxf, isgx3, iko, ikxf, ikx3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iysew, iyse2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy3, isgo, isgyf, isgy3, iko, ikyf, iky3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izsew, izse2, mhalf, mfull, isclp)
      call muldif(io, izf, iz3, isgo, isgzf, isgz3, iko, ikzf, ikz3,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, ixthw, ixth2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix4, isgo, isgxf, isgx4, iko, ikxf, ikx4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, iythw, iyth2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy4, isgo, isgyf, isgy4, iko, ikyf, iky4,
     *              nkmax, mhalf)
      call decmp2(io, isgo, iko, izthw, izth2, mhalf, mfull, isclp)
      call muldif(io, izf, iz4, isgo, isgzf, isgz4, iko, ikzf, ikz4,
     *              nkmax, mhalf)
c
      call mulmul(iy5, iz3, io, isgy5, isgz3, isgo, iky5, ikz3, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix4, iv, isgo, isgx4, isgv, iko, ikx4, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz5, ix3, io, isgz5, isgx3, isgo, ikz5, ikx3, iko,
     *              nkmax, mhalf)
      call mulmul(io, iy4, iu, isgo, isgy4, isgu, iko, iky4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix5, iy3, iv, isgx5, isgy3, isgv, ikx5, iky3, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iz4, iu, isgv, isgz4, isgu, ikv, ikz4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz5, iy3, io, isgz5, isgy3, isgo, ikz5, iky3, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix4, iu, isgo, isgx4, isgu, iko, ikx4, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix5, iz3, iv, isgx5, isgz3, isgv, ikx5, ikz3, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy4, iu, isgv, isgy4, isgu, ikv, iky4, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iy5, ix3, io, isgy5, isgx3, isgo, iky5, ikx3, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz4, iu, isgo, isgz4, isgu, iko, ikz4, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(iy6, iz2, iv, isgy6, isgz2, isgv, iky6, ikz2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, ix4, iu, isgv, isgx4, isgu, ikv, ikx4, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz6, ix2, io, isgz6, isgx2, isgo, ikz6, ikx2, iko,
     *              nkmax, mhalf)
      call mulmul(io, iy4, iu, isgo, isgy4, isgu, iko, iky4, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix6, iy2, iv, isgx6, isgy2, isgv, ikx6, iky2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iz4, iu, isgv, isgz4, isgu, ikv, ikz4, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iz6, iy2, io, isgz6, isgy2, isgo, ikz6, iky2, iko,
     *              nkmax, mhalf)
      call mulmul(io, ix4, iu, isgo, isgx4, isgu, iko, ikx4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix6, iz2, iv, isgx6, isgz2, isgv, ikx6, ikz2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy4, iu, isgv, isgy4, isgu, ikv, iky4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(io, iu, iv, isgo, isgu, isgv, iko, iku, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iy6, ix2, io, isgy6, isgx2, isgo, iky6, ikx2, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz4, iu, isgo, isgz4, isgu, iko, ikz4, iku,
     *              nkmax, mhalf)
      isgu =-isgu
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*DECMP2
c
c     subroutine decmp2 
c
c     to decompose a regular or non-regular length integer
c
      subroutine decmp2(ia, isga, ika, iwi, iwi2, mhalf, mfull, isclp)
c
      integer ia(*), isga, ika, iwi, iwi2, mhalf, mfull, isclp(*)
      integer nkmax
      parameter (nkmax = 30)
      integer iu(nkmax), io(nkmax), isgu, isgo, iku, iko, isgcl, ikcl
c
      call decomp(ia, isga, iwi, mhalf)
      ika = 2
      if(iwi2.ne.0) then
         isgcl = 1
         ikcl = 2
         call mulmul(ia, isclp, iu, isga, isgcl, isgu, ika, ikcl,
     *                 iku, nkmax, mhalf)
         if(iwi2.eq.mfull) iwi2 = 0
         call decomp(io, isgo, iwi2, mhalf)
         isgo = -isgo
         iko = 2
         call muldif(iu, io, ia, isgu, isgo, isga, iku, iko, ika,
     *                 nkmax, mhalf)
      endif
c
      return
      end
*DECOMP
c
c     subroutine decomp
c
c     to decompose a regular length integer
c
c     iwi = isga*(ia(1) + ia(2) * mhalf)
c
c     iwi  is a regular length integer
c     isga is a sign integer (-1, 0, 1)
c     ia(1) and ia(2) are integers less than mhalf
c
      subroutine decomp(ia, isga, iwi, mhalf)
c
      integer ia(*), isga, iwi, mhalf, ivi
c
      if(iwi.gt.0) then
         isga = 1
         ivi = iwi
      elseif(iwi.lt.0) then
         isga =-1
         ivi = -iwi
      else
         isga = 0
         ia(1) = 0
         ia(2) = 0
         return
      endif
      ia(2) = ivi/mhalf
      ia(1) = ivi - ia(2)*mhalf
c
      return
      end
*MULMUL
c
c     subroutine mulmul
c
c     to perform a multiple precision integer multiplication
c     (for multiplying 2 or more integers)
c
c     io = ia * ib
c
c     ia represents a decomposed integer
c     ib represents a decomposed integer
c     io is the product of ia and ib in its decomposed form
c
      subroutine mulmul(ia, ib, io, isga, isgb, isgo, ika, ikb, iko,
     *                    nkmax, mhalf)
c
      integer ia(*), ib(*), io(*)
      integer isga, isgb, isgo, ika, ikb, iko, nkmax, mhalf
      integer i, ipt, ipr, iko1, k, j
c
      if(isga.eq.0.or.isgb.eq.0)then
         isgo=0
         iko = 2
         io(1) = 0
         io(2) = 0
         return
      endif
c
      iko = ika + ikb
      if(iko.gt.nkmax) stop 4710
c
      if(isga.gt.0)then
         if(isgb.gt.0)then
            isgo = 1
         else
            isgo =-1
         endif
      else
         if(isgb.gt.0)then
            isgo =-1
         else
            isgo = 1
         endif
      endif
c
      iko1 = iko - 1
      ipr = 0
c
      do 200 i = 1, iko1
         ipt = ipr
         k = i
         do 180 j = 1, ikb
            if(k .lt. 1) go to 190
            if(k .gt. ika) go to 150
            ipt = ipt + ia(k)*ib(j)
  150       continue
            k = k - 1
  180    continue
  190    continue
         ipr = ipt/mhalf
         io(i) = ipt - ipr*mhalf
  200 continue
c
      io(iko) = ipr
      if(ipr.ge.mhalf) stop 4720
c
      iko1 = iko
      do 300 i = iko1, ika+1, -1
         if(io(i) .ne. 0) go to 400
         iko = iko - 1
  300 continue
  400 continue
c
      return
      end
*MULDIF
c
c     subroutine muldif
c
c     to perform a multiple precision integer subtraction
c     (for subtracting a decomposed product from another)
c
c     io = ia - ib
c
c     ia represents a decomposed regular length integer or the 
c        decomposed product of two or more regular length integers
c     ib is similarly described
c     io is a decomposed integer which represents ia - ib
c
      subroutine muldif(ia, ib, io, isga, isgb, isgo, ika, ikb, iko,
     *                    nkmax, mhalf)
c
      integer ia(*), ib(*), io(*)
      integer isga, isgb, isgo, ika, ikb, iko, nkmax, mhalf
      integer i, iko1, irel
c
      if(isgb.eq.0)then
         if(isga.eq.0)then
            isgo=0
            iko = 2
            io(1) = 0
            io(2) = 0
            return
         endif
         isgo = isga
         iko = ika
         do 100 i=1,iko
            io(i) = ia(i)
  100    continue
      elseif(isga.eq.0)then
         isgo =-isgb
         iko = ikb
         do 200 i=1,iko
            io(i) = ib(i)
  200    continue
      else
         iko = ika
         if(ikb.lt.ika) then
            do 300 i=ikb+1,ika
               ib(i) = 0
  300       continue
         elseif(ika.lt.ikb) then
            iko = ikb
            do 400 i=ika+1,ikb
               ia(i) = 0
  400       continue
         endif
         if(isga*isgb.gt.0)then
            irel = 0
            do 500 i = iko, 1, -1
               if(ia(i).gt.ib(i))then
                  irel = 1
                  go to 600
               elseif(ia(i).lt.ib(i))then
                  irel = -1
                  go to 600
               endif
  500       continue
  600       continue
            if(irel.eq.0)then
               isgo = 0
               do 700 i=1,iko
                  io(i) = 0
  700          continue
            else
               isgo=isga*irel
               io(1) = (ia(1)-ib(1))*irel
               do 800 i=2,iko
                  if(io(i-1).lt.0) then
                     io(i) =-1
                     io(i-1) = io(i-1) + mhalf
                  else
                     io(i) = 0
                  endif
                  io(i) = io(i) + (ia(i)-ib(i))*irel
  800          continue
               if(io(iko).lt.0) stop 4810
            endif
         else
            isgo=isga
            io(1) = ia(1)+ib(1)
            do 900 i=2,iko
               if(io(i-1).ge.mhalf) then
                  io(i) = 1
                  io(i-1) = io(i-1) - mhalf
               else
                  io(i) = 0
               endif
               io(i) = io(i) + ia(i)+ib(i)
  900       continue
            if(io(iko).ge.mhalf) then
               iko = iko+1
               if(iko.gt.nkmax) stop 4820
               io(iko) = 1
               io(iko-1) = io(iko-1) - mhalf
            endif
         endif
      endif
c
      if(iko .eq. 2) go to 1400
      iko1 = iko
      do 1300 i = iko1, 3, -1
         if(io(i) .ne. 0) go to 1400
         iko = iko - 1
 1300 continue
 1400 continue
c
      return
      end
*IWSIGN
c
c     subroutine for determining sign of weight of point ifir minus
c     weight of point isec
c
      subroutine iwsign(w, w2, ifir, isec, mhalf, mfull, isclw, ipout)
c
      integer w(*), w2(*)
      integer ifir, isec, mhalf, mfull, nkmax, ipout
      parameter (nkmax = 30)
      integer isclw(*), iu(nkmax), iw1(nkmax), iw2(nkmax)
      integer iwfiw, iwsew, iwfi2, iwse2
      integer isgw1, isgw2, ikw1, ikw2, isgu, iku
c
      iwfiw = w(ifir)
      iwsew = w(isec)
c
      iwfi2 = w2(ifir)
      iwse2 = w2(isec)
c
      call decmp2(iw1,isgw1,ikw1, iwfiw,iwfi2, mhalf, mfull, isclw)
      call decmp2(iw2,isgw2,ikw2, iwsew,iwse2, mhalf, mfull, isclw)
      call muldif(iw1, iw2, iu, isgw1, isgw2, isgu, ikw1, ikw2, iku,
     *              nkmax, mhalf)
c
      ipout = isgu
c
      return
      end
*IQSIGN
c
c     subroutine for determining position of point ifif with respect
c     to sphere determined by (weighted) points ifir, isec, ithi, ifou
c     if positive then ifif is outside the sphere
c     if negative then ifif is inside the sphere
c     if zero then ifif is in the surface of the sphere
c
      subroutine iqsign(x, y, z, w, x2, y2, z2, w2, ifir, isec,
     *                  ithi, ifou, ifif, mhalf, mfull, isclp,
     *                  isclw, isclr, delaun, ipout)
c
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer ifir, isec, ithi, ifou, ifif, mhalf, mfull, nkmax, ipout
      integer isclp(*), isclw(*), isclr(*)
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax), ip(nkmax)
      integer iq2(nkmax), iq3(nkmax), iq4(nkmax), iq5(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ix4(nkmax), iy4(nkmax), iz4(nkmax)
      integer ix5(nkmax), iy5(nkmax), iz5(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixf2(nkmax), iyf2(nkmax), izf2(nkmax)
      integer iwf(nkmax), iw2(nkmax), iw3(nkmax), iw4(nkmax), iw5(nkmax)
      logical delaun
      integer iwfuw, iwsew, iwthw, iwfow, iwfiw
      integer ixthw, iythw, izthw, ixfow, iyfow, izfow
      integer ixfiw, iyfiw, izfiw, ixsew, iysew, izsew
      integer ixfuw, iyfuw, izfuw
      integer iwfu2, iwse2, iwth2, iwfo2, iwfi2
      integer ixth2, iyth2, izth2, ixfo2, iyfo2, izfo2
      integer ixfi2, iyfi2, izfi2, ixse2, iyse2, izse2
      integer ixfu2, iyfu2, izfu2
      integer isgw2, isgw3, isgw4, isgw5, ikw2, ikw3, ikw4, ikw5
      integer isgq2, isgq3, isgq4, isgq5, ikq2, ikq3, ikq4, ikq5
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgxf2, isgyf2, isgzf2, ikxf2, ikyf2, ikzf2
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx4, isgy4, isgz4, ikx4, iky4, ikz4
      integer isgx5, isgy5, isgz5, ikx5, iky5, ikz5
      integer isgo, isgu, isgv, isgp, iko, iku, ikv, ikp
      integer isgwf, isgcl, ikwf, ikcl
c
      if(delaun) then
         isgw2 = 0
         isgw3 = 0
         isgw4 = 0
         isgw5 = 0
      else
         iwfuw = w(ifir)
         iwsew = w(isec)
         iwthw = w(ithi)
         iwfow = w(ifou)
         iwfiw = w(ifif)
c
         iwfu2 = w2(ifir)
         iwse2 = w2(isec)
         iwth2 = w2(ithi)
         iwfo2 = w2(ifou)
         iwfi2 = w2(ifif)
c
         call decmp2(iwf,isgwf,ikwf, iwfuw,iwfu2, mhalf, mfull, isclw)
         isgcl = 1
         ikcl = 2
         call decmp2(io, isgo, iko, iwsew, iwse2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw2, isgu, isgcl, isgw2, iku, ikcl,
     *                 ikw2, nkmax, mhalf)
         call decmp2(io, isgo, iko, iwthw, iwth2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw3, isgu, isgcl, isgw3, iku, ikcl,
     *                 ikw3, nkmax, mhalf)
         call decmp2(io, isgo, iko, iwfow, iwfo2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw4, isgu, isgcl, isgw4, iku, ikcl,
     *                 ikw4, nkmax, mhalf)
         call decmp2(io, isgo, iko, iwfiw, iwfi2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw5, isgu, isgcl, isgw5, iku, ikcl,
     *                 ikw5, nkmax, mhalf)
      endif
c
      ixfuw = x(ifir)
      iyfuw = y(ifir)
      izfuw = z(ifir)
      ixsew = x(isec)
      iysew = y(isec)
      izsew = z(isec)
      ixthw = x(ithi)
      iythw = y(ithi)
      izthw = z(ithi)
      ixfow = x(ifou)
      iyfow = y(ifou)
      izfow = z(ifou)
      ixfiw = x(ifif)
      iyfiw = y(ifif)
      izfiw = z(ifif)
c
      ixfu2 = x2(ifir)
      iyfu2 = y2(ifir)
      izfu2 = z2(ifir)
      ixse2 = x2(isec)
      iyse2 = y2(isec)
      izse2 = z2(isec)
      ixth2 = x2(ithi)
      iyth2 = y2(ithi)
      izth2 = z2(ithi)
      ixfo2 = x2(ifou)
      iyfo2 = y2(ifou)
      izfo2 = z2(ifou)
      ixfi2 = x2(ifif)
      iyfi2 = y2(ifif)
      izfi2 = z2(ifif)
c
      call decmp2(ixf, isgxf, ikxf, ixfuw, ixfu2, mhalf, mfull, isclp)
      call decmp2(iyf, isgyf, ikyf, iyfuw, iyfu2, mhalf, mfull, isclp)
      call decmp2(izf, isgzf, ikzf, izfuw, izfu2, mhalf, mfull, isclp)
      call mulmul(ixf, ixf, ixf2, isgxf, isgxf, isgxf2, ikxf, ikxf,
     *              ikxf2, nkmax, mhalf)
      call mulmul(iyf, iyf, iyf2, isgyf, isgyf, isgyf2, ikyf, ikyf,
     *              ikyf2, nkmax, mhalf)
      call mulmul(izf, izf, izf2, isgzf, isgzf, isgzf2, ikzf, ikzf,
     *              ikzf2, nkmax, mhalf)
      if(isgxf2.lt.0 .or. isgyf2.lt.0 .or. isgzf2.lt.0) stop 5105
c
      call frterm(ixsew, iysew, izsew, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw2, ix2, iy2, iz2, iq2, isgw2, isgx2,
     *            isgy2, isgz2, isgq2, ikw2, ikx2, iky2, ikz2, ikq2,
     *            mhalf, mfull, ixse2, iyse2, izse2, isclp)
c
      call frterm(ixthw, iythw, izthw, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw3, ix3, iy3, iz3, iq3, isgw3, isgx3,
     *            isgy3, isgz3, isgq3, ikw3, ikx3, iky3, ikz3, ikq3,
     *            mhalf, mfull, ixth2, iyth2, izth2, isclp)
c
      call frterm(ixfow, iyfow, izfow, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw4, ix4, iy4, iz4, iq4, isgw4, isgx4,
     *            isgy4, isgz4, isgq4, ikw4, ikx4, iky4, ikz4, ikq4,
     *            mhalf, mfull, ixfo2, iyfo2, izfo2, isclp)
c
      call frterm(ixfiw, iyfiw, izfiw, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw5, ix5, iy5, iz5, iq5, isgw5, isgx5,
     *            isgy5, isgz5, isgq5, ikw5, ikx5, iky5, ikz5, ikq5,
     *            mhalf, mfull, ixfi2, iyfi2, izfi2, isclp)
c
      call mulmul(iq5, ix2, iv, isgq5, isgx2, isgv, ikq5, ikx2, ikv,
     *              nkmax, mhalf)
      call mulmul(iq5, ix3, iu, isgq5, isgx3, isgu, ikq5, ikx3, iku,
     *              nkmax, mhalf)
      call mulmul(iq5, ix4, ip, isgq5, isgx4, isgp, ikq5, ikx4, ikp,
     *              nkmax, mhalf)
      call detrm3(iv, iy2, iz2, isgv, isgy2, isgz2,
     *            iu, iy3, iz3, isgu, isgy3, isgz3,
     *            ip, iy4, iz4, isgp, isgy4, isgz4,
     *            ikv, iku, ikp, iky2, iky3, iky4,
     *            ikz2, ikz3, ikz4, io, isgo, iko, mhalf)
c
      call detrm3(iq2, iy2, iz2, isgq2, isgy2, isgz2,
     *            iq3, iy3, iz3, isgq3, isgy3, isgz3,
     *            iq4, iy4, iz4, isgq4, isgy4, isgz4,
     *            ikq2, ikq3, ikq4, iky2, iky3, iky4,
     *            ikz2, ikz3, ikz4, iu, isgu, iku, mhalf)
      call mulmul(iu, ix5, ip, isgu, isgx5, isgp, iku, ikx5, ikp,
     *              nkmax, mhalf)
      call muldif(io, ip, iv, isgo, isgp, isgv, iko, ikp, ikv,
     *              nkmax, mhalf)
c
      call detrm3(iq2, iz2, ix2, isgq2, isgz2, isgx2,
     *            iq3, iz3, ix3, isgq3, isgz3, isgx3,
     *            iq4, iz4, ix4, isgq4, isgz4, isgx4,
     *            ikq2, ikq3, ikq4, ikz2, ikz3, ikz4,
     *            ikx2, ikx3, ikx4, iu, isgu, iku, mhalf)
      call mulmul(iu, iy5, io, isgu, isgy5, isgo, iku, iky5, iko,
     *              nkmax, mhalf)
      call muldif(iv, io, ip, isgv, isgo, isgp, ikv, iko, ikp,
     *              nkmax, mhalf)
c
      call detrm3(iq2, ix2, iy2, isgq2, isgx2, isgy2,
     *            iq3, ix3, iy3, isgq3, isgx3, isgy3,
     *            iq4, ix4, iy4, isgq4, isgx4, isgy4,
     *            ikq2, ikq3, ikq4, ikx2, ikx3, ikx4,
     *            iky2, iky3, iky4, iu, isgu, iku, mhalf)
      call mulmul(iu, iz5, iv, isgu, isgz5, isgv, iku, ikz5, ikv,
     *              nkmax, mhalf)
      call muldif(ip, iv, io, isgp, isgv, isgo, ikp, ikv, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*IQSIG1
c
c     subroutine for determining position of point ifif with respect to
c     line segment determined by (weighted) points ifir, isec, assuming
c     ifir, isec and ifif are colinear
c     if positive then ifif is in the exterior of the line segment
c     if negative then ifif is in the interior of the line segment
c     if zero then ifif is one of the endpoints
c
      subroutine iqsig1(x, y, z, w, x2, y2, z2, w2, ifir, isec, ifif,
     *                 mhalf, mfull, isclp, isclw, isclr, delaun, ipout)
c
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer ifir, isec, ifif, mhalf, mfull, nkmax, ipout
      integer isclp(*), isclw(*), isclr(*)
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax)
      integer iq2(nkmax), iq5(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix5(nkmax), iy5(nkmax), iz5(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixf2(nkmax), iyf2(nkmax), izf2(nkmax)
      integer iwf(nkmax), iw2(nkmax), iw5(nkmax)
      logical delaun
      integer iwfuw, iwsew, iwfiw
      integer ixfiw, iyfiw, izfiw, ixsew, iysew, izsew
      integer ixfuw, iyfuw, izfuw
      integer iwfu2, iwse2, iwfi2
      integer ixfi2, iyfi2, izfi2, ixse2, iyse2, izse2
      integer ixfu2, iyfu2, izfu2
      integer isgw2, isgw5, ikw2, ikw5
      integer isgq2, isgq5, ikq2, ikq5
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgxf2, isgyf2, isgzf2, ikxf2, ikyf2, ikzf2
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx5, isgy5, isgz5, ikx5, iky5, ikz5
      integer isgo, isgu, isgv, iko, iku, ikv
      integer isgwf, isgcl, ikwf, ikcl
c
      if(delaun) then
         isgw2 = 0
         isgw5 = 0
      else
         iwfuw = w(ifir)
         iwsew = w(isec)
         iwfiw = w(ifif)
c
         iwfu2 = w2(ifir)
         iwse2 = w2(isec)
         iwfi2 = w2(ifif)
c
         call decmp2(iwf,isgwf,ikwf, iwfuw,iwfu2, mhalf, mfull, isclw)
         isgcl = 1
         ikcl = 2
         call decmp2(io, isgo, iko, iwsew, iwse2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw2, isgu, isgcl, isgw2, iku, ikcl,
     *                 ikw2, nkmax, mhalf)
         call decmp2(io, isgo, iko, iwfiw, iwfi2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw5, isgu, isgcl, isgw5, iku, ikcl,
     *                 ikw5, nkmax, mhalf)
      endif
c
      ixfuw = x(ifir)
      iyfuw = y(ifir)
      izfuw = z(ifir)
      ixsew = x(isec)
      iysew = y(isec)
      izsew = z(isec)
      ixfiw = x(ifif)
      iyfiw = y(ifif)
      izfiw = z(ifif)
c
      ixfu2 = x2(ifir)
      iyfu2 = y2(ifir)
      izfu2 = z2(ifir)
      ixse2 = x2(isec)
      iyse2 = y2(isec)
      izse2 = z2(isec)
      ixfi2 = x2(ifif)
      iyfi2 = y2(ifif)
      izfi2 = z2(ifif)
c
      call decmp2(ixf, isgxf, ikxf, ixfuw, ixfu2, mhalf, mfull, isclp)
      call decmp2(iyf, isgyf, ikyf, iyfuw, iyfu2, mhalf, mfull, isclp)
      call decmp2(izf, isgzf, ikzf, izfuw, izfu2, mhalf, mfull, isclp)
      call mulmul(ixf, ixf, ixf2, isgxf, isgxf, isgxf2, ikxf, ikxf,
     *              ikxf2, nkmax, mhalf)
      call mulmul(iyf, iyf, iyf2, isgyf, isgyf, isgyf2, ikyf, ikyf,
     *              ikyf2, nkmax, mhalf)
      call mulmul(izf, izf, izf2, isgzf, isgzf, isgzf2, ikzf, ikzf,
     *              ikzf2, nkmax, mhalf)
      if(isgxf2.lt.0 .or. isgyf2.lt.0 .or. isgzf2.lt.0) stop 5205
c
      call frterm(ixsew, iysew, izsew, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw2, ix2, iy2, iz2, iq2, isgw2, isgx2,
     *            isgy2, isgz2, isgq2, ikw2, ikx2, iky2, ikz2, ikq2,
     *            mhalf, mfull, ixse2, iyse2, izse2, isclp)
c
      call frterm(ixfiw, iyfiw, izfiw, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw5, ix5, iy5, iz5, iq5, isgw5, isgx5,
     *            isgy5, isgz5, isgq5, ikw5, ikx5, iky5, ikz5, ikq5,
     *            mhalf, mfull, ixfi2, iyfi2, izfi2, isclp)
c
      call detrm1(iq5, ix2, iy2, iz2, ix2, iy2, iz2, isgq5, isgx2,
     *            isgy2, isgz2, isgx2, isgy2, isgz2, ikq5, ikx2, iky2,
     *            ikz2, ikx2, iky2, ikz2, iu, isgu, iku, mhalf)
c
      call detrm1(iq2, ix2, iy2, iz2, ix5, iy5, iz5, isgq2, isgx2,
     *            isgy2, isgz2, isgx5, isgy5, isgz5, ikq2, ikx2, iky2,
     *            ikz2, ikx5, iky5, ikz5, iv, isgv, ikv, mhalf)
c
      call muldif(iu, iv, io, isgu, isgv, isgo, iku, ikv, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*IQSIG2
c
c     subroutine for determining position of point ifif with respect
c     to circle determined by (weighted) points ifir, isec, ithi,
c     if positive then ifif is outside the circle
c     if negative then ifif is inside the circle
c     if zero then ifif is in the circle
c
      subroutine iqsig2(x, y, z, w, x2, y2, z2, w2, ifir, isec, ithi,
     *           ifif, mhalf, mfull, isclp, isclw, isclr, delaun, ipout)
c
      integer x(*), y(*), z(*), w(*), x2(*), y2(*), z2(*), w2(*)
      integer ifir, isec, ithi, ifif, mhalf, mfull, nkmax, ipout
      integer isclp(*), isclw(*), isclr(*)
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax), ip(nkmax)
      integer iq2(nkmax), iq3(nkmax), iq5(nkmax)
      integer ix2(nkmax), iy2(nkmax), iz2(nkmax)
      integer ix3(nkmax), iy3(nkmax), iz3(nkmax)
      integer ix5(nkmax), iy5(nkmax), iz5(nkmax)
      integer ixf(nkmax), iyf(nkmax), izf(nkmax)
      integer ixf2(nkmax), iyf2(nkmax), izf2(nkmax)
      integer iwf(nkmax), iw2(nkmax), iw3(nkmax), iw5(nkmax)
      logical delaun
      integer iwfuw, iwsew, iwthw, iwfiw
      integer ixthw, iythw, izthw, ixfiw, iyfiw, izfiw
      integer ixfuw, iyfuw, izfuw, ixsew, iysew, izsew
      integer iwfu2, iwse2, iwth2, iwfi2
      integer ixth2, iyth2, izth2, ixfi2, iyfi2, izfi2
      integer ixfu2, iyfu2, izfu2, ixse2, iyse2, izse2
      integer isgw2, isgw3, isgw5, ikw2, ikw3, ikw5
      integer isgq2, isgq3, isgq5, ikq2, ikq3, ikq5
      integer isgxf, isgyf, isgzf, ikxf, ikyf, ikzf
      integer isgxf2, isgyf2, isgzf2, ikxf2, ikyf2, ikzf2
      integer isgx2, isgy2, isgz2, ikx2, iky2, ikz2
      integer isgx3, isgy3, isgz3, ikx3, iky3, ikz3
      integer isgx5, isgy5, isgz5, ikx5, iky5, ikz5
      integer isgo, isgu, isgv, isgp, iko, iku, ikv, ikp
      integer isgwf, isgcl, ikwf, ikcl
c
      if(delaun) then
         isgw2 = 0
         isgw3 = 0
         isgw5 = 0
      else
         iwfuw = w(ifir)
         iwsew = w(isec)
         iwthw = w(ithi)
         iwfiw = w(ifif)
c
         iwfu2 = w2(ifir)
         iwse2 = w2(isec)
         iwth2 = w2(ithi)
         iwfi2 = w2(ifif)
c
         call decmp2(iwf,isgwf,ikwf, iwfuw,iwfu2, mhalf, mfull, isclw)
         isgcl = 1
         ikcl = 2
         call decmp2(io, isgo, iko, iwsew, iwse2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw2, isgu, isgcl, isgw2, iku, ikcl,
     *                 ikw2, nkmax, mhalf)
         call decmp2(io, isgo, iko, iwthw, iwth2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw3, isgu, isgcl, isgw3, iku, ikcl,
     *                 ikw3, nkmax, mhalf)
         call decmp2(io, isgo, iko, iwfiw, iwfi2, mhalf, mfull, isclw)
         call muldif(io, iwf, iu, isgo, isgwf, isgu, iko, ikwf, iku,
     *                 nkmax, mhalf)
         call mulmul(iu, isclr, iw5, isgu, isgcl, isgw5, iku, ikcl,
     *                 ikw5, nkmax, mhalf)
      endif
c
      ixfuw = x(ifir)
      iyfuw = y(ifir)
      izfuw = z(ifir)
      ixsew = x(isec)
      iysew = y(isec)
      izsew = z(isec)
      ixthw = x(ithi)
      iythw = y(ithi)
      izthw = z(ithi)
      ixfiw = x(ifif)
      iyfiw = y(ifif)
      izfiw = z(ifif)
c
      ixfu2 = x2(ifir)
      iyfu2 = y2(ifir)
      izfu2 = z2(ifir)
      ixse2 = x2(isec)
      iyse2 = y2(isec)
      izse2 = z2(isec)
      ixth2 = x2(ithi)
      iyth2 = y2(ithi)
      izth2 = z2(ithi)
      ixfi2 = x2(ifif)
      iyfi2 = y2(ifif)
      izfi2 = z2(ifif)
c
      call decmp2(ixf, isgxf, ikxf, ixfuw, ixfu2, mhalf, mfull, isclp)
      call decmp2(iyf, isgyf, ikyf, iyfuw, iyfu2, mhalf, mfull, isclp)
      call decmp2(izf, isgzf, ikzf, izfuw, izfu2, mhalf, mfull, isclp)
      call mulmul(ixf, ixf, ixf2, isgxf, isgxf, isgxf2, ikxf, ikxf,
     *              ikxf2, nkmax, mhalf)
      call mulmul(iyf, iyf, iyf2, isgyf, isgyf, isgyf2, ikyf, ikyf,
     *              ikyf2, nkmax, mhalf)
      call mulmul(izf, izf, izf2, isgzf, isgzf, isgzf2, ikzf, ikzf,
     *              ikzf2, nkmax, mhalf)
      if(isgxf2.lt.0 .or. isgyf2.lt.0 .or. isgzf2.lt.0) stop 5305
c
      call frterm(ixsew, iysew, izsew, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw2, ix2, iy2, iz2, iq2, isgw2, isgx2,
     *            isgy2, isgz2, isgq2, ikw2, ikx2, iky2, ikz2, ikq2,
     *            mhalf, mfull, ixse2, iyse2, izse2, isclp)
c
      call frterm(ixthw, iythw, izthw, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw3, ix3, iy3, iz3, iq3, isgw3, isgx3,
     *            isgy3, isgz3, isgq3, ikw3, ikx3, iky3, ikz3, ikq3,
     *            mhalf, mfull, ixth2, iyth2, izth2, isclp)
c
      call frterm(ixfiw, iyfiw, izfiw, ixf, iyf, izf, isgxf, isgyf,
     *            isgzf, ikxf, ikyf, ikzf, ixf2, iyf2, izf2,
     *            isgxf2, isgyf2, isgzf2, ikxf2,
     *            ikyf2, ikzf2, iw5, ix5, iy5, iz5, iq5, isgw5, isgx5,
     *            isgy5, isgz5, isgq5, ikw5, ikx5, iky5, ikz5, ikq5,
     *            mhalf, mfull, ixfi2, iyfi2, izfi2, isclp)
c
      call detrm0(iq5, iy2, iz2, iy3, iz3, isgq5, isgy2, isgz2,
     *            isgy3, isgz3, ikq5, iky2, ikz2, iky3, ikz3,
     *            iv, isgv, ikv, mhalf)
c
      call detrm0(iq5, iz2, ix2, iz3, ix3, isgq5, isgz2, isgx2,
     *            isgz3, isgx3, ikq5, ikz2, ikx2, ikz3, ikx3,
     *            iu, isgu, iku, mhalf)
c
      call detrm0(iq5, ix2, iy2, ix3, iy3, isgq5, isgx2, isgy2,
     *            isgx3, isgy3, ikq5, ikx2, iky2, ikx3, iky3,
     *            ip, isgp, ikp, mhalf)
c
      call detrm3(iv, ix2, ix3, isgv, isgx2, isgx3,
     *            iu, iy2, iy3, isgu, isgy2, isgy3,
     *            ip, iz2, iz3, isgp, isgz2, isgz3,
     *            ikv, iku, ikp, ikx2, iky2, ikz2,
     *            ikx3, iky3, ikz3, io, isgo, iko, mhalf)
c
      call detrm2(iq2, ix2, iy2, iz2, isgq2, isgx2, isgy2, isgz2,
     *            iq3, ix3, iy3, iz3, isgq3, isgx3, isgy3, isgz3,
     *            ikq2, ikx2, iky2, ikz2, ikq3, ikx3, iky3, ikz3,
     *            iu, isgu, iku, mhalf)
      call mulmul(iu, ix5, ip, isgu, isgx5, isgp, iku, ikx5, ikp,
     *              nkmax, mhalf)
      call muldif(io, ip, iv, isgo, isgp, isgv, iko, ikp, ikv,
     *              nkmax, mhalf)
c
      call detrm2(iq2, iy2, iz2, ix2, isgq2, isgy2, isgz2, isgx2,
     *            iq3, iy3, iz3, ix3, isgq3, isgy3, isgz3, isgx3,
     *            ikq2, iky2, ikz2, ikx2, ikq3, iky3, ikz3, ikx3,
     *            iu, isgu, iku, mhalf)
      call mulmul(iu, iy5, io, isgu, isgy5, isgo, iku, iky5, iko,
     *              nkmax, mhalf)
      call muldif(iv, io, ip, isgv, isgo, isgp, ikv, iko, ikp,
     *              nkmax, mhalf)
c
      call detrm2(iq2, iz2, ix2, iy2, isgq2, isgz2, isgx2, isgy2,
     *            iq3, iz3, ix3, iy3, isgq3, isgz3, isgx3, isgy3,
     *            ikq2, ikz2, ikx2, iky2, ikq3, ikz3, ikx3, iky3,
     *            iu, isgu, iku, mhalf)
      call mulmul(iu, iz5, iv, isgu, isgz5, isgv, iku, ikz5, ikv,
     *              nkmax, mhalf)
      call muldif(ip, iv, io, isgp, isgv, isgo, ikp, ikv, iko,
     *              nkmax, mhalf)
c
      ipout = isgo
c
      return
      end
*FRTERM
c
      subroutine frterm(ixsew, iysew, izsew, ixf, iyf, izf, isgxf,
     *                  isgyf, isgzf, ikxf, ikyf, ikzf, ixf2, iyf2,
     *                  izf2, isgxf2, isgyf2, isgzf2, ikxf2,
     *                  ikyf2, ikzf2, iw2, ix2, iy2, iz2,
     *                  iq2, isgw2, isgx2, isgy2, isgz2, isgq2, ikw2,
     *                  ikx2, iky2, ikz2, ikq2, mhalf, mfull,
     *                  ixse2, iyse2, izse2, isclp)
c
      integer ixsew, iysew, izsew, isgxf, isgyf, isgzf
      integer isgxf2, isgyf2, isgzf2
      integer ikxf, ikyf, ikzf, ikxf2, ikyf2, ikzf2
      integer isgw2, isgx2, isgy2, isgz2, isgq2
      integer ikw2, ikx2, iky2, ikz2, ikq2, mhalf, mfull
      integer ixse2, iyse2, izse2, isclp(*)
      integer isgo, isgu, isgv, isgp, iko, iku, ikv, ikp
      integer ixf(*), iyf(*), izf(*), ixf2(*), iyf2(*), izf2(*)
      integer iw2(*), ix2(*), iy2(*), iz2(*), iq2(*)
      integer nkmax
      parameter (nkmax = 30)
      integer io(nkmax), iu(nkmax), iv(nkmax), ip(nkmax)
c
      call decmp2(io, isgo, iko, ixsew, ixse2, mhalf, mfull, isclp)
      call muldif(io, ixf, ix2, isgo, isgxf, isgx2, iko, ikxf, ikx2,
     *              nkmax, mhalf)
      call mulmul(io, io, iu, isgo, isgo, isgu, iko, iko, iku,
     *              nkmax, mhalf)
      call muldif(iu, ixf2, iv, isgu, isgxf2, isgv, iku, ikxf2, ikv,
     *              nkmax, mhalf)
      call muldif(iv, iw2, ip, isgv, isgw2, isgp, ikv, ikw2, ikp,
     *              nkmax, mhalf)
c
      call decmp2(io, isgo, iko, iysew, iyse2, mhalf, mfull, isclp)
      call muldif(io, iyf, iy2, isgo, isgyf, isgy2, iko, ikyf, iky2,
     *              nkmax, mhalf)
      call mulmul(io, io, iu, isgo, isgo, isgu, iko, iko, iku,
     *              nkmax, mhalf)
      call muldif(iu, iyf2, iv, isgu, isgyf2, isgv, iku, ikyf2, ikv,
     *              nkmax, mhalf)
      isgv=-isgv
      call muldif(ip, iv, iu, isgp, isgv, isgu, ikp, ikv, iku,
     *              nkmax, mhalf)
c
      call decmp2(io, isgo, iko, izsew, izse2, mhalf, mfull, isclp)
      call muldif(io, izf, iz2, isgo, isgzf, isgz2, iko, ikzf, ikz2,
     *              nkmax, mhalf)
      call mulmul(io, io, iv, isgo, isgo, isgv, iko, iko, ikv,
     *              nkmax, mhalf)
      call muldif(iv, izf2, ip, isgv, isgzf2, isgp, ikv, ikzf2, ikp,
     *              nkmax, mhalf)
      isgp=-isgp
      call muldif(iu, ip, iq2, isgu, isgp, isgq2, iku, ikp, ikq2,
     *              nkmax, mhalf)
c
      return
      end
*DETRM0
c
      subroutine detrm0(iq, ix2, iy2, ix3, iy3, isgq, isgx2, isgy2,
     *                  isgx3, isgy3, ikq, ikx2, iky2, ikx3, iky3,
     *                  io, isgo, iko, mhalf)
c
      integer nkmax
      parameter (nkmax = 30)
      integer io(*), iu(nkmax), iv(nkmax), iw(nkmax)
      integer iq(*), ix2(*), iy2(*), ix3(*), iy3(*)
      integer isgq, isgx2, isgy2, isgx3, isgy3
      integer ikq, ikx2, iky2, ikx3, iky3, isgo, iko, mhalf
      integer isgu, isgv, isgw, iku, ikv, ikw
c
      call mulmul(iq, ix2, iv, isgq, isgx2, isgv, ikq, ikx2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy3, iu, isgv, isgy3, isgu, ikv, iky3, iku,
     *              nkmax, mhalf)
      call mulmul(iq, iy2, iv, isgq, isgy2, isgv, ikq, iky2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, ix3, iw, isgv, isgx3, isgw, ikv, ikx3, ikw,
     *              nkmax, mhalf)
      call muldif(iu, iw, io, isgu, isgw, isgo, iku, ikw, iko,
     *              nkmax, mhalf)
c
      return
      end
*DETRM1
c
      subroutine detrm1(iq, ix2, iy2, iz2, ix3, iy3, iz3,
     *                  isgq, isgx2, isgy2, isgz2, isgx3, isgy3, isgz3,
     *                  ikq, ikx2, iky2, ikz2, ikx3, iky3, ikz3,
     *                  io, isgo, iko, mhalf)
c
      integer nkmax
      parameter (nkmax = 30)
      integer io(*), iv(nkmax), iw(nkmax)
      integer iq(*), ix2(*), iy2(*), iz2(*), ix3(*), iy3(*), iz3(*)
      integer isgq, isgx2, isgy2, isgz2, isgx3, isgy3, isgz3
      integer ikq, ikx2, iky2, ikz2, ikx3, iky3, ikz3
      integer isgo, isgv, isgw, iko, ikv, ikw, mhalf
c
      call mulmul(iq, ix2, iv, isgq, isgx2, isgv, ikq, ikx2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, ix3, io, isgv, isgx3, isgo, ikv, ikx3, iko,
     *              nkmax, mhalf)
c
      call mulmul(iq, iy2, iv, isgq, isgy2, isgv, ikq, iky2, ikv,
     *              nkmax, mhalf)
      call mulmul(iv, iy3, iw, isgv, isgy3, isgw, ikv, iky3, ikw,
     *              nkmax, mhalf)
      isgw =-isgw
      call muldif(io, iw, iv, isgo, isgw, isgv, iko, ikw, ikv,
     *              nkmax, mhalf)
c
      call mulmul(iq, iz2, io, isgq, isgz2, isgo, ikq, ikz2, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz3, iw, isgo, isgz3, isgw, iko, ikz3, ikw,
     *              nkmax, mhalf)
      isgw =-isgw
      call muldif(iv, iw, io, isgv, isgw, isgo, ikv, ikw, iko,
     *              nkmax, mhalf)
c
      return
      end
*DETRM2
c
      subroutine detrm2(iq2, ix2, iy2, iz2, isgq2, isgx2, isgy2, isgz2,
     *                  iq3, ix3, iy3, iz3, isgq3, isgx3, isgy3, isgz3,
     *                  ikq2, ikx2, iky2, ikz2, ikq3, ikx3, iky3, ikz3,
     *                  io, isgo, iko, mhalf)
c
      integer nkmax
      parameter (nkmax = 30)
      integer io(*), ip(nkmax), ir(nkmax), iv(nkmax), iw(nkmax)
      integer iq2(*), iq3(*)
      integer ix2(*), iy2(*), iz2(*), ix3(*), iy3(*), iz3(*)
      integer isgq2, isgx2, isgy2, isgz2, isgq3, isgx3, isgy3, isgz3
      integer ikq2, ikx2, iky2, ikz2, ikq3, ikx3, iky3, ikz3
      integer isgo, isgp, isgr, isgv, isgw, iko, ikp, ikr, ikv, ikw
      integer mhalf
c
      call detrm0(iq2, ix2, iy2, ix3, iy3, isgq2, isgx2, isgy2,
     *            isgx3, isgy3, ikq2, ikx2, iky2, ikx3, iky3,
     *            iv, isgv, ikv, mhalf)
      call mulmul(iv, iy3, io, isgv, isgy3, isgo, ikv, iky3, iko,
     *              nkmax, mhalf)
c
      call detrm0(iq2, iz2, ix2, iz3, ix3, isgq2, isgz2, isgx2,
     *            isgz3, isgx3, ikq2, ikz2, ikx2, ikz3, ikx3,
     *            iv, isgv, ikv, mhalf)
      call mulmul(iv, iz3, ip, isgv, isgz3, isgp, ikv, ikz3, ikp,
     *              nkmax, mhalf)
c
      call muldif(io, ip, ir, isgo, isgp, isgr, iko, ikp, ikr,
     *              nkmax, mhalf)
c
      call detrm0(iq3, ix2, iy2, ix3, iy3, isgq3, isgx2, isgy2,
     *            isgx3, isgy3, ikq3, ikx2, iky2, ikx3, iky3,
     *            iv, isgv, ikv, mhalf)
      call mulmul(iv, iy2, io, isgv, isgy2, isgo, ikv, iky2, iko,
     *              nkmax, mhalf)
c
      call detrm0(iq3, iz2, ix2, iz3, ix3, isgq3, isgz2, isgx2,
     *            isgz3, isgx3, ikq3, ikz2, ikx2, ikz3, ikx3,
     *            iv, isgv, ikv, mhalf)
      call mulmul(iv, iz2, ip, isgv, isgz2, isgp, ikv, ikz2, ikp,
     *              nkmax, mhalf)
c
      call muldif(io, ip, iw, isgo, isgp, isgw, iko, ikp, ikw,
     *              nkmax, mhalf)
      call muldif(ir, iw, io, isgr, isgw, isgo, ikr, ikw, iko,
     *              nkmax, mhalf)
c
      return
      end
*DETRM3
c
      subroutine detrm3(ix2, iy2, iz2, isgx2, isgy2, isgz2,
     *                  ix3, iy3, iz3, isgx3, isgy3, isgz3,
     *                  ix4, iy4, iz4, isgx4, isgy4, isgz4,
     *                  ikx2, ikx3, ikx4, iky2, iky3, iky4,
     *                  ikz2, ikz3, ikz4, io, isgo, iko, mhalf)
c
      integer nkmax
      parameter (nkmax = 30)
      integer io(*), iu(nkmax), iv(nkmax), iw(nkmax)
      integer ix2(*), iy2(*), iz2(*), ix3(*), iy3(*), iz3(*)
      integer ix4(*), iy4(*), iz4(*)
      integer isgx2, isgy2, isgz2, isgx3, isgy3, isgz3
      integer isgx4, isgy4, isgz4, ikx2, ikx3, ikx4, iky2, iky3, iky4
      integer ikz2, ikz3, ikz4, isgo, iko, mhalf
      integer isgu, isgv, isgw, iku, ikv, ikw
c
      call mulmul(ix3, iy4, iv, isgx3, isgy4, isgv, ikx3, iky4, ikv,
     *              nkmax, mhalf)
      call mulmul(ix4, iy3, iu, isgx4, isgy3, isgu, ikx4, iky3, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, iw, isgv, isgu, isgw, ikv, iku, ikw,
     *              nkmax, mhalf)
      call mulmul(iw, iz2, io, isgw, isgz2, isgo, ikw, ikz2, iko,
     *              nkmax, mhalf)
c
      call mulmul(ix2, iy4, iv, isgx2, isgy4, isgv, ikx2, iky4, ikv,
     *              nkmax, mhalf)
      call mulmul(ix4, iy2, iu, isgx4, isgy2, isgu, ikx4, iky2, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, iw, isgv, isgu, isgw, ikv, iku, ikw,
     *              nkmax, mhalf)
      call mulmul(iw, iz3, iu, isgw, isgz3, isgu, ikw, ikz3, iku,
     *              nkmax, mhalf)
      call muldif(io, iu, iw, isgo, isgu, isgw, iko, iku, ikw,
     *              nkmax, mhalf)
c
      call mulmul(ix3, iy2, iv, isgx3, isgy2, isgv, ikx3, iky2, ikv,
     *              nkmax, mhalf)
      call mulmul(ix2, iy3, iu, isgx2, isgy3, isgu, ikx2, iky3, iku,
     *              nkmax, mhalf)
      call muldif(iv, iu, io, isgv, isgu, isgo, ikv, iku, iko,
     *              nkmax, mhalf)
      call mulmul(io, iz4, iu, isgo, isgz4, isgu, iko, ikz4, iku,
     *              nkmax, mhalf)
      call muldif(iw, iu, io, isgw, isgu, isgo, ikw, iku, iko,
     *              nkmax, mhalf)
c
      return
      end
