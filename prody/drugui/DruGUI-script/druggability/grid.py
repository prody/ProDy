# Druggability: Python Package and VMD GUI for Druggability Index Analysis
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""This module defines classes for parsing, smoothing and writing grid files.

In Druggability calculations, grid files contain probe atom counts (or 
occupancy). This module defines a base class and file format specific 
derived classes. These classes are direcly used by the 
:mod:`druggability.probe` module.

Classes:

* :class:`Grid`
* :class:`OpenDX`
* :class:`XPLOR`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'
__version__ = '0.5.2'



"""
File Type Documentation
-----------------------

XPLOR Crystallographic map files
--------------------------------
This information is from XPLOR documentation:
    http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html

The X-PLOR program is able to write electron density map files in either a 
binary format or in an ASCII format. The binary format is more compact and may 
be read more quickly than the ASCII format but has the disadvantage that it may
not be readable when transferred between different kinds of computer. The ASCII
formatted file is written if the FORMatted keyword is set TRUE (the default); 
setting FORMatted=FALSE will cause a binary map file to be written.

Map header

The X-PLOR map file begins with an eight-line header.
1.	   Line 1
An empty line written by the `/ ` FORTRAN format descriptor in the formatted 
map file.
2.	   Lines 2- 5
Title information written as character strings. These lines are written as 
80-character strings in the formatted file map.
3.	   Line 6
A series of nine integers NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX. The 
values NA, NB and NC indicate the total number of grid points along the a,b, 
and c cell edges. The items AMIN, AMAX, BMIN, BMAX, CMIN, CMAX indicate the 
starting and stopping grid points along each cell edge in the portion of the 
map that is written. In the formatted map file this line is written using the 
FORTRAN format statement (9I8).
4.	   Line 7
A series of six double-precision items corresponding to the crystal cell 
dimensions a, b, c, alpha, beta, gamma. In the formatted map file these items 
are written using the FORTRAN format statement (6E12.5).
5.	   Line 8
A three-letter character string which always reads `ZXY'.

Density array
Following the map header, the density matrix is then written section-by-section
with c moving slowest (in z-sections). Each section of the density map is 
preceded by a section number.
Thus, for the formatted map file each section of the density map is written 
using FORTRAN statements of the type

WRITE(OUNIT,'(I8)') KSECT

WRITE(OUNIT,'(6E12.5)') ((SECTON(I,J),I=1,ISECT),J=1,JSECT)

and the resulting map is written with six pixels on each line. The binary 
format is identical except the format statements are missing, so that each 
line that is written contains the entire length of map along the `fast' 
(a-axis) direction.


Map footer

Two lines follow the density array.

1.	   Line 1
The integer `-9999' is always written. For the formatted map file, The FORTRAN 
format statement (I8) is used to write this value.

2.	   Line 2
Two double-precision items corresponding to the average electron density and 
the standard deviation of the map. For the formatted map file these items are 
written using the FORTRAN format statement (2(E12.4,1X)).


OpenDX scalar data
------------------
This information is from APBS website:
    <http://www.poissonboltzmann.org/file-formats/mesh-and-data-formats/
    opendx-scalar-data>

We output most discretized scalar data (e.g., potential, accessibility, 
etc.) from APBS in the data format used by the OpenDX software package. 
The OpenDX data format is very flexible; the following sections describe 
the application of this format for APBS multigrid and finite element 
datasets.
The multigrid data format has the following form::
   
   # Comments
   object 1 class gridpositions counts nx ny nz
   origin xmin ymin zmin
   delta hx 0.0 0.0
   delta 0.0 hy 0.0 
   delta 0.0 0.0 hz
   object 2 class gridconnections counts nx ny nz
   object 3 class array type double rank 0 items n data follows
   u(0,0,0) u(0,0,1) u(0,0,2)
   ...
   u(0,0,nz-3) u(0,0,nz-2) u(0,0,nz-1)
   u(0,1,0) u(0,1,1) u(0,1,2)
   ...
   u(0,1,nz-3) u(0,1,nz-2) u(0,1,nz-1)
   ...
   u(0,ny-1,nz-3) u(0,ny-1,nz-2) u(0,ny-1,nz-1)
   u(1,0,0) u(1,0,1) u(1,0,2)
   ...
   attribute "dep" string "positions"
   object "regular positions regular connections" class field
   component "positions" value 1
   component "connections" value 2
   component "data" value 3
   
The variables in this format have been shown in bold and include
Comments
Any number of comment lines, each line starting with the "#" symbol
nx ny nz
The number of grid points in the x-, y-, and z-directions.
xmin ymin zmin
The coordinates of the grid lower corner.
hx hy hz
The grid spacings in the x-, y-, and z-directions.
n
The total number of grid points; n = nx * ny * nz
u(*,*,*)
The data values, ordered with the z-index increasing most quickly, 
followed by the y-index, and then the x-index.

"""

import os.path

import numpy as np

from druggability.exceptions import GridError

__all__ = ['Grid']

class Grid(object):

    """Base class for manipulating grid data.
    
    An important grid attribute is :attr:`state`. After the grid file is 
    parsed, this is set to *original*. This attribute is compared or changed
    when two grids are added or a grid :meth:`smooth` is called.

    .. attribute:: filename

       Name of the grid file.
       
    .. attribute:: format
    
       Grid file format.
    
    .. attribute:: name

       Name of the instance, which is deduced from the filename. Name is used
       as a prefix to output file names.

    .. attribute:: array 

       3-dimentional Numpy array that holds grid data.

    .. attribute:: offset

       Offset of the origin of the grid from origin of the Cartesian coordinate 
       space. 

    .. attribute:: shape

       Shape of the data array.

    .. attribute:: spacing

       Resolution of the grid.
    
    .. attribute:: state

       State of the grid. *original* and *smooth* are predefined states.
    
    """
    
    
    def __init__(self, filename=None):
        """Instantiate grid from a file.
        
        A grid maybe initialized without a filename. If a filename is provided,
        file format specific derivatives of this class parse the file at
        initialization. To delay this, one can initialize a grid without a 
        filename and then use :meth:`parse` method to parse the grid file. 
        
        """
        
        self.name = None
        self.filename = None
        self.array = None
        self.state = None
        self.offset = None
        self.spacing = None
        self.shape = None
        
    def __repr__(self):
        return 'Grid {0:s} (file: {1:s}) in {2:s} state'.format(self.name, 
                                                    self.filename, self.state)

    def __add__(self, other):
        """Compare grid attributes to determine if grids can be added.
        
        File type specific Grid derivatives need to implement method for 
        adding grid data array.
        
        Before adding data arrays of the grids, four attributes are compared. 
        These are :attr:`state`, :attr:`shape`, :attr:`spacing`, and 
        :attr:`offset`. All of these attributes must have the same values.  
        
        """
        if self.state != other.state:
            raise GridError('both grids must be at the same state')
        elif self.shape != other.shape:
            raise GridError('both grids must have the same shape')
        elif all(self.spacing != other.spacing):
            raise GridError('both grids must have the same spacing')
        elif all(self.offset != other.offset):
            raise GridError('both grids must have the same offset')
        else:
            return True


    def _smooth(self):
        """Smooth grid array by averaging over neighboring grid elements."""        
        self.array = (self.array[0:-2, 0:-2, 0:-2] +
                      self.array[0:-2, 0:-2, 1:-1] +
                      self.array[0:-2, 0:-2, 2:  ] +
                      self.array[0:-2, 1:-1, 0:-2] +
                      self.array[0:-2, 1:-1, 1:-1] +
                      self.array[0:-2, 1:-1, 2:  ] +
                      self.array[0:-2, 2:,   0:-2] +
                      self.array[0:-2, 2:,   1:-1] +
                      self.array[0:-2, 2:,   2:  ] +
                      self.array[1:-1, 0:-2, 0:-2] +
                      self.array[1:-1, 0:-2, 1:-1] +
                      self.array[1:-1, 0:-2, 2:  ] +
                      self.array[1:-1, 1:-1, 0:-2] +
                      self.array[1:-1, 1:-1, 1:-1] +
                      self.array[1:-1, 1:-1, 2:  ] +
                      self.array[1:-1, 2:,   0:-2] +
                      self.array[1:-1, 2:,   1:-1] +
                      self.array[1:-1, 2:,   2:  ] +
                      self.array[2:,   0:-2, 0:-2] +
                      self.array[2:,   0:-2, 1:-1] +
                      self.array[2:,   0:-2, 2:  ] +
                      self.array[2:,   1:-1, 0:-2] +
                      self.array[2:,   1:-1, 1:-1] +
                      self.array[2:,   1:-1, 2:  ] +
                      self.array[2:,   2:,   0:-2] +
                      self.array[2:,   2:,   1:-1] +
                      self.array[2:,   2:,   2:  ]) / 27.0
        self.state = 'smooth'
        
    def parse(self, filename):
        """File type specific method for parsing grid data."""
        pass

    def write(self, filename=None):
        """File type specific method for writing grid data."""
        pass

    def smooth(self):
        """File type specific method for averaging grid data.
        
        Smoothing is performed by assigning a grid element the value found by 
        averaging values of itself and its neighboring grid elements. After 
        this operation grid becomes smaller by two elements along each
        direction.

        Calling this method changes the :attr:`state` attribute of the grid 
        from *original* to *smooth*.

        """
        pass

        
class XPLOR(Grid):
    
    """A class to manipulate XPLOR formatted contour file.
    
    This class is tested using grid files outputed by ptraj from `AmberTools
    <http://ambermd.org/>`_.
    
    """
    
    
    def __init__(self, filename=None):
        """Instantiation arguments are passed to :class:`Grid`"""
        Grid.__init__(self, filename)
        self.format = 'Xplor'
        self._size = None
        self._ignored_lines = None
        self._first = None
        self._last = None
        self._angles = None
        
        if filename:
            self.parse(filename)

        
    def __add__(self, other):
        if not Grid.__add__(self, other):
            raise GridError('{0:s} and {1:s} cannot be added'
                            .format(self, other))
        grid = XPLOR(None)
        grid.name = '(' + self.name + ' + '+ other.name + ')'
        grid.array = self.array + other.array
        grid.state = self.state
        
        grid.offset = self.offset
        grid.spacing = self.spacing
        grid.shape = self.shape
        
        grid._ignored_lines = self._ignored_lines
        grid._size = self._size
        grid._first = self._first
        grid._last = self._last
        grid._angles = self._angles
        return grid

    def smooth(self):
        """Smooth grid and change grid attributes and state."""
        self._smooth()
        self.shape = np.shape(self.array)
        self.offset += self.spacing

        self._size -= self.spacing * 2
        self._first += 1
        self._last -= 1
    

    def parse(self, filename):
        """Parse grid data from file."""
        if not os.path.isfile(filename):
            raise IOError('{0:s} not found'.format(filename))
        self.filename = filename
        self.name = os.path.splitext(os.path.split(filename)[1])[0]

        xplor_file = open(filename)
        xplor = xplor_file.read()
        xplor_file.close()
        
        xyz = xplor.index('ZYX')
        lines = xplor[:xyz].split('\n')
        
        self._ignored_lines = lines[:3] #Hold lines that're not used
        line = lines[3].split() #Parse indexing related numbers
        self.shape = np.array((int(line[0]), int(line[3]), int(line[6])), 
                              'int64')
        self._first = np.array((int(line[1]), int(line[4]), int(line[7])), 
                              'int64')
        self._last = np.array((int(line[2]), int(line[5]), int(line[8])), 
                              'int64')
        
        line = lines[4].split()
        self._size = np.array((float(line[0]), float(line[1]), float(line[2])), 
                             'd')
        self.spacing = self._size / self.shape
        self._angles = (float(line[3]), float(line[4]), float(line[5]))
        
        self.offset = (self._first - 0.5) * self.spacing

        array = np.fromstring(xplor[xyz+4:], dtype='d', sep=' ')
        self.array = np.zeros(self.shape, 'd')
        yxshape = (self.shape[1], self.shape[0])
        length = yxshape[0] * yxshape[1] + 1  
        for k in range(self.shape[2]):
            self.array[:, :, k] = array[length*k+1:length*(k+1)].reshape(
                                                                    yxshape).T
        self.state = 'original'

    def write(self, filename=None):
        """Write grid data into a file.
        
        If a filename is not provided, gridname_state.xplor will be used.
            
        """
        if filename is None:
            filename = os.path.splitext(self.filename)
            filename = filename[0] + '_' + self.state + filename[1]

        xplor_file = open(filename, 'w')
        for line in self._ignored_lines:
            xplor_file.write(line+'\n')
        xplor_file.write(('{0[0]:8d}{1[0]:8d}{2[0]:8d}'
                          '{0[1]:8d}{1[1]:8d}{2[1]:8d}'
                          '{0[2]:8d}{1[2]:8d}{2[2]:8d}\n').format(
                                          self.shape, self._first, self._last))
        xplor_file.write(('{0[0]:12.3f}{0[1]:12.3f}{0[2]:12.3f}'
                          '{1[0]:12.3f}{1[1]:12.3f}{1[2]:12.3f}\n').format(
                                                     self._size, self._angles))
        xplor_file.write('ZYX\n')
        
        format_ = ''
        for i in range(self.shape[0]):
            if i != 0 and i % 6 == 0:
                format_ += '\n'
            format_ += '{0['+str(i)+']:12.5f}'
        else:
            if i % 6 != 0:
                format_ += '\n'

        for k in range(self.shape[2]):
            xplor_file.write('{0:8d}\n'.format(k + self._first[2]))
            for j in range(self.shape[1]):
                xplor_file.write(format_.format(self.array[:, j, k]))
        xplor_file.close()
        return filename
    
class OpenDX(Grid):

    """A class to manipulate OpenDX scalar data files.
    
    This classes is tested using grid files outputed by volmap in `VMD
    <http://www.ks.uiuc.edu/Research/vmd/current/ug/>`_.
    
    Additional relevant information on this file format may be found here 
    `APBS website <http://www.poissonboltzmann.org/file-formats/
    mesh-and-data-formats/opendx-scalar-data>`_
    
    """


    def __init__(self, filename=None):
        """Instantiation arguments are passed to :class:`Grid`"""
        Grid.__init__(self, filename)
        self.format = 'OpenDX'
        self._origin = None
        self._comments = None
        
        if filename:
            self.parse(filename)

    def __add__(self, other):
        if not Grid.__add__(self, other):
            raise GridError('{0:s} and {1:s} cannot be added'
                            .format(self, other))
        grid = OpenDX(None)
        grid.name = '(' + self.name + ' + '+ other.name + ')'
        grid.array = self.array + other.array
        grid.state = self.state
        
        grid.offset = self.offset
        grid.spacing = self.spacing
        grid.shape = self.shape

        grid._comments = self._comments
        grid._origin = self._origin
        return grid


    def parse(self, filename):
        """Parse grid data from file."""

        if not os.path.isfile(filename):
            raise IOError('{0:s} not found'.format(filename))
        self.filename = filename
        self.name = os.path.splitext(os.path.split(filename)[1])[0]

        opendx_file = open(filename)
        opendx = opendx_file.read()
        opendx_file.close()
        
        lindex = opendx.index('data follows') + len('data follows') + 1
        lines = opendx[:lindex].split('\n')
        self._comments = []
        self.spacing = np.zeros(3, 'd')
        for line in lines:
            if line.startswith('#'):
                self._comments.append(line)
            elif line.startswith('object 1'):
                items = line.split()
                self.shape = (int(items[-3]), int(items[-2]), int(items[-1]))
            elif line.startswith('origin'):
                items = line.split()
                self._origin = np.array(items[1:], 'd')
            elif line.startswith('delta'):
                items = line.split()
                self.spacing += np.array(items[1:], 'd')
        rindex = opendx.rindex('object')
        self._comments.append(opendx[rindex:].strip())
        self.offset = self._origin - self.spacing / 2
        self.array = np.fromstring(opendx[lindex:rindex], dtype='d', sep=' '
                       ).reshape((self.shape[0], self.shape[1], self.shape[2]))
      
        self.state = 'original'


    def smooth(self):
        """Smooth grid and change grid attributes and state."""
        self._smooth()
        self.shape = self.array.shape
        self.offset += self.spacing
        
        self._origin += self.spacing
        
    def write(self, filename=None):
        """Write grid data into a file.
        
        If a filename is not provided, gridname_state.dx will be used.

        """
        if filename is None:
            filename = os.path.splitext(self.filename)
            filename = filename[0] + '_' + self.state + filename[1]

        opendx = open(filename, 'w')
        opendx.write('{0:s} modified by Druggability\n'
                     .format(self._comments[0]))
        opendx.write('object 1 class gridpositions counts {0[0]:d} {0[1]:d} '
                     '{0[2]:d}\n'.format(self.shape))
        opendx.write('origin {0[0]:.9g} {0[1]:.9g} {0[2]:.9g}\n'
                     .format(self._origin))
        opendx.write('delta {0:.9g} 0 0\n'.format(self.spacing[0]))
        opendx.write('delta 0 {0:.9g} 0\n'.format(self.spacing[1]))
        opendx.write('delta 0 0 {0:.9g}\n'.format(self.spacing[2]))
        opendx.write('object 2 class gridconnections counts {0[0]:d} {0[1]:d} '
                     '{0[2]:d}\n'.format(self.shape))
        length = self.shape[0]*self.shape[1]*self.shape[2]
        opendx.write('object 3 class array type double rank 0 items {0:d} data'
                     ' follows\n'.format(length))
        
        array = self.array.flatten()
        string = ''
        times = length / 9
        n_times = int(times)
        for i in range(n_times):
            string += ('{0[0]:.9g} {0[1]:.9g} {0[2]:.9g}\n'
                       '{0[3]:.9g} {0[4]:.9g} {0[5]:.9g}\n'
                       '{0[6]:.9g} {0[7]:.9g} {0[8]:.9g}\n'
                      ).format(
                       array[i*9:i*9+9])
        length = length - n_times * 9
        n_times = int(length / 3)
        for i in range(n_times):
            string += '{0[0]:.9g} {0[1]:.9g} {0[2]:.9g}\n'.format(
                        array[i*3:i*3+3])
        length = length - n_times * 3
        if length == 2:
            string += '{0[0]:.9g} {0[1]:.9g}\n'.format(array[-2:])            
        elif length == 1:
            string += '{0:.9g}\n'.format(array[-1])
        opendx.write(string)
        opendx.write('\n{0:s}\n'.format(self._comments[1]))
        opendx.close()
        return filename
