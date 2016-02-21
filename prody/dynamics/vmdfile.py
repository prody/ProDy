"""This module defines input and output functions for NMD format.

.. _nmd-format:

NMD Format
-------------------------------------------------------------------------------

Description
^^^^^^^^^^^

NMD files (extension :file:`.nmd`) are plain text files that contain at
least normal mode and system coordinate data.

NMD files can be visualized using :ref:`nmwiz`.  ProDy functions
:func:`.writeNMD` and :func:`.parseNMD` can be used to read and write NMD
files.

Data fields
^^^^^^^^^^^

Data fields in bold face are required. All data arrays and lists must be in a
single line and items must be separated by one or more space characters.

**coordinates**: system coordinates as a list of decimal numbers
  Coordinate array is the most important line in an NMD file. All mode array
  lengths must match the length of the coordinate array. Also, number of atoms
  in the system is deduced from the length of the coordinate array.

::

  coordinates 27.552 4.354 23.629 24.179 4.807 21.907 ...

**mode**: normal mode array as a list of decimal numbers
  Optionally, mode index and a scaling factor may be provided
  in the same line as a mode array. Both of these must precede the mode array.
  Providing a scaling factor enables relative scaling of the mode arrows and
  the amplitude of the fluctuations in animations. For NMA, scaling factors
  may be chosen to be the square-root of the inverse-eigenvalue associated
  with the mode. Analogously, for PCA data, scaling factor would be the
  square-root of the eigenvalue.

  If a mode line contains numbers preceding the mode array, they are evaluated
  based on their type. If an integer is encountered, it is considered the mode
  index. If a decimal number is encountered, it is considered the scaling
  factor. Scaling factor may be the square-root of the inverse eigenvalue
  if data is from an elastic network model, or the square-root of the
  eigenvalue if data is from an essential dynamics (or principal component)
  analysis.

  For example, all of the following lines are valid. The first line contains
  mode index and scaling factor. Second and third lines contain mode index or
  scaling factor. Last line contains only the mode array.

::

  mode 1 2.37    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 1    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 2.37    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 0.039 0.009 0.058 0.038 -0.011 0.052 0.043  ...

*name*: name of the model

The length of all following data fields must be equal to the number of atoms in
the system. NMWiz uses such data when writing a temporary PDB files for
loading coordinate data into VMD.

*atomnames*: list of atom names
  If not provided, all atom names are set to "CA".

*resnames*: list of residue names
  If not provided, all residue names are set to "GLY".

*chainids*: list of chain identifiers
  If not provided, all chain identifiers are set to "A".

*resids*: list of residue numbers
  If not provided, residue numbers are started from 1 and incremented by one
  for each atom.

*bfactors*: list of experimental beta-factors
  If not provided, all beta-factors are set to zero.
  Beta-factors can be used to color the protein representation.

NMD files may contain additional lines. Only lines that start with one of the
above field names are evaluated by NMWiz.


Autoload Trick
^^^^^^^^^^^^^^

By adding a special line in an NMD file, file content can be automatically
loaded into VMD at startup. The first line calls a NMWiz function to load the
file itself (:file:`xyzeros.nmd`).

::

  nmwiz_load xyzeros.nmd
  coordinates 0 0 0 0 0 0  ...
  mode 0.039 0.009 0.058 0.038 -0.011 0.052 ...
  mode -0.045 -0.096 -0.009 -0.040 -0.076 -0.010 ...
  mode 0.007 -0.044 0.080 0.015 -0.037 0.062 ...


In this case, VMD must be started from the command line by typing
:program:`vmd -e xyzeros.nmd`."""


__all__ = ['writeVMDstiffness']

import os
from os.path import abspath, join, split, splitext

import numpy as np

from prody import LOGGER, SETTINGS, PY3K
from prody.atomic import AtomGroup
from prody.utilities import openFile, isExecutable, which, PLATFORM, addext
from prody.proteins import writePDB

from .nma import NMA
from .anm import ANM
from .gnm import GNM, ZERO
from .pca import PCA
from .mode import Vector, Mode
from .modeset import ModeSet
from .nmdfile import viewNMDinVMD, pathVMD, getVMDpath, setVMDpath

def writeVMDstiffness(filename, model, pdb, indices, k_range):
    """Return *filename* that contains *modes* and *atoms* data in NMD format
    described in :ref:`nmd-format`.  :file:`.nmd` extension is appended to
    filename, if it does not have an extension.

    .. note::
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`.Vector` instance is given, it will be normalized
          before it is written. It's length before normalization will be
          written as the scaling factor of the vector."""

    try:
        coords = (pdb._getCoords() if hasattr(pdb, '_getCoords') else
                pdb.getCoords())
    except AttributeError:
        try:
            checkCoords(pdb)
        except TypeError:
            raise TypeError('pdb must be a Numpy array or an object '
                            'with `getCoords` method')
    
    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance')
    elif not model.is3d():
        raise TypeError('model must be a 3-dimensional NMA instance')
    elif len(model) == 0:
        raise ValueError('model must have normal modes calculated')
    elif model.getStiffness() is None:
        raise ValueError('model must have stiffness matrix calculated')
    elif len(indices)==0:
        raise ValueError('indices cannot be an empty array')

    index = np.zeros(2)
    if len(indices)==1:
        index[1]=indices[0]
        index[0]=indices[0]
    elif len(indices)==2:
        index[1]=indices[1]
        index[0]=indices[0]

    out = openFile(addext(filename, '.tcl'), 'w')
    writePDB(filename + '.pdb', pdb)

    LOGGER.info('Creating VMD file.')
    
    out.write('display rendermode GLSL \n')
    out.write('display projection orthographic\n')
    out.write('color Display Background white\n')
    out.write('display shadows on\n')
    out.write('display depthcue off\n')
    out.write('axes location off\n')
    out.write('stage location off\n')
    out.write('light 0 on\n')
    out.write('light 1 on\n')
    out.write('light 2 off\n')
    out.write('light 3 on\n')
    out.write('mol addrep 0\n')
    out.write('display resetview\n')
    out.write('mol new {./'+str(filename)+'.pdb} type {pdb} first 0 last -1 step 1 waitfor 1\n')
    out.write('mol modselect 0 0 protein\n')
    out.write('mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0\n')
    out.write('mol modcolor 0 0 Structure\n')
    out.write('mol color Structure\n')
    out.write('mol representation NewCartoon 0.300000 10.000000 4.100000 0\n')
    out.write('mol selection protein\n')
    out.write('mol material Opaque\n')

    colors = ['blue', 'red', 'gray', 'orange','yellow', 'tan','silver', 'green', \
    'white', 'pink', 'cyan', 'purple', 'lime', 'mauve', 'ochre', 'iceblue', 'black', \
    'yellow2','yellow3','green2','green3','cyan2','cyan3','blue2','blue3','violet', \
    'violet2','magenta','magenta2','red2','red3','orange2','orange3']*50
    
    color_nr = 1 # starting from red color in VMD
    ResCounter = []
    for r in xrange(indices[0], indices[1]):
        baza_col = [] # Value of Kij is here for each residue
        nr_baza_col = [] # Resid of aa are here
        out.write("draw color "+str(colors[color_nr])+"\n")
            
        for nr_i, i in enumerate(model.getStiffness()[r]):
            if k_range[0] < float(i) < k_range[1]:
                baza_col.append(i)
                nr_baza_col.append(nr_i)
                coords_all_ca = pdb.select('calpha')
                resid_r = str(pdb.getResnames()[r])+str(r)
                resid_r2 = str(pdb.getResnames()[nr_i])+str(nr_i)
                    
                if len(baza_col) == 0: # if base is empty then it will not change the color
                    color_nr = 0
                else:
                    out.write("draw line "+'{'+str(coords[r])[1:-1]+'} {'+str(coords[nr_i])[1:-1]+'} width 3 style solid \n')
                    ResCounter.append(len(baza_col))
                        
            else: pass
        
            if len(baza_col) != 0:
                out.write('mol addrep 0\n')
                out.write('mol modselect '+str(color_nr+1)+' 0 protein and name CA and resid '+ str(r)+' '+str(nr_baza_col)[1:-1].replace(',','')+'\n')
                out.write('mol modcolor '+str(color_nr+1)+' 0 ColorID '+str(color_nr)+'\n')
                out.write('mol modstyle '+str(color_nr+1)+' 0 VDW 0.600000 12.000000\n')
                out.write('mol color ColorID '+str(color_nr)+'\n')
                out.write('mol representation VDW 1.000000 12.000000 \n')
                out.write('mol selection protein and name CA and resid '+ str(r)+' '+str(nr_baza_col)[1:-1].replace(',','')+'\n')
                out.write('mol material Opaque \n')
                color_nr = color_nr + 1
                
    out.write('mol addrep 0\n')
    out.close()
        
    if len(ResCounter) > 0:
        return out
    elif len(ResCounter) == 0:
        LOGGER.info('There is no residue pair in this Kij range.')
        return 'None'   
