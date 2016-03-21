# -*- coding: utf-8 -*-
"""This module defines TCL file for VMD program."""

__all__ = ['writeVMDstiffness', 'writeDeformProfile']

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

def writeVMDstiffness(model, pdb, indices, k_range, filename='vmd_out', \
                            selstr='protein and name CA', loadToVMD=True):
   
    """Return three *filename* files: (1) PDB file with coordinates. 
    (2) TCL file containing vmd commands for loading PDB file with accurate 	
    vmd representation. Pair of residues with selected *k_range* of 
    effective spring constant are shown in VMD respresentation with 
    solid line between them.
    If more than one residue will be selected in *indices*, different pair 
    for each residue will be colored in the different colors.    
    (3) TXT file contains pair of residues with effective spring constant in
    selected range *k_range*.    

    The effective spring constant calculation using ``buildSM`` method
    from :class:`.ANM`.
    
    .. note::
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`.Vector` instance is given, it will be normalized
          before it is written. It's length before normalization will be
          written as the scaling factor of the vector.
          

    :arg model: this is an 3-dimensional NMA instance from a :class:`.ANM
        calculations
    :type model: :class:`.ANM`
    :arg pdb: a coordinate set or an object with ``getCoords`` method
    :type pdb: :class:`numpy.ndarray`. 
    :arg indices: amino acid number.
    :type indices: ``[int, int]`` or ``[int]`` for one amino acid    
    :arg k_range: effective force constant value.
    :type k_range: int or float, ``[int, int]``
    
    By default files are saved as *filename* and loaded to VMD program and 
    *selstr* is a selection from :class:`.Select`
     """

    try:
        coords_sel = pdb.select(selstr)
        resnum_list = coords_sel.getResnums()  
        coords = (coords_sel._getCoords() if hasattr(coords_sel, '_getCoords') else
                coords_sel.getCoords())
    except AttributeError:
        try:
            checkCoords(coords_sel)
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

    if len(indices)==1:
        indices0=indices[0]-resnum_list[0]
        indices1=indices[0]-resnum_list[0]
    elif len(indices)==2:
        indices0=indices[0]-resnum_list[0]
        indices1=indices[1]-resnum_list[0]

    out = openFile(addext(filename, '.tcl'), 'w')
    out_txt = openFile(addext(filename,'.txt'), 'w')
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
    for r in xrange(indices0, indices1+1):
        baza_col = [] # Value of Kij is here for each residue
        nr_baza_col = [] # Resid of aa are here
        out.write("draw color "+str(colors[color_nr])+"\n")
            
        for nr_i, i in enumerate(model.getStiffness()[r]):
            if k_range[0] < float(i) < k_range[1]:
                baza_col.append(i)
                nr_baza_col.append(nr_i+resnum_list[0])
                resid_r = str(coords_sel.getResnames()[r])+str(r+resnum_list[0])
                resid_r2 = str(coords_sel.getResnames()[nr_i])+str(nr_i+resnum_list[0])
                    
                if len(baza_col) == 0: # if base is empty then it will not change the color
                    color_nr = 0
                else:
                    out.write("draw line "+'{'+str(coords[r])[1:-1]+'} {'+\
                       str(coords[nr_i])[1:-1]+'} width 3 style solid \n')
                    out_txt.write(str(resid_r)+'\t'+resid_r2+'\t'+str(i)+'\n')
                    ResCounter.append(len(baza_col))
                        
            else: pass
        
        if len(baza_col) != 0:
            out.write('mol addrep 0\n')
            out.write('mol modselect '+str(color_nr+1)+' 0 protein and name CA and resid '+ \
                       str(r+resnum_list[0])+' '+str(nr_baza_col)[1:-1].replace(',','')+'\n')
            out.write('mol modcolor '+str(color_nr+1)+' 0 ColorID '+str(color_nr)+'\n')
            out.write('mol modstyle '+str(color_nr+1)+' 0 VDW 0.600000 12.000000\n')
            out.write('mol color ColorID '+str(color_nr)+'\n')
            out.write('mol representation VDW 1.000000 12.000000 \n')
            out.write('mol selection protein and name CA and resid '+ \
            str(r+resnum_list[0])+' '+str(nr_baza_col)[1:-1].replace(',','')+'\n')
            
            out.write('mol material Opaque \n')
            color_nr = color_nr + 1
                
    out.write('mol addrep 0\n')
    out.close()
    out_txt.close()

    if (loadToVMD == True):
        from prody import pathVMD
        LOGGER.info('File will be loaded to VMD program.')
        os.system(pathVMD()+" -e "+str(filename)+".tcl")
                                        
    if len(ResCounter) > 0:
        return out
    elif len(ResCounter) == 0:
        LOGGER.info('There is no residue pair in this Kij range.')
        return 'None'   


def writeDeformProfile(model, pdb, filename='dp_out', selstr='protein and name CA',\
                                            pdb_selstr='protein', loadToVMD=True):

    """Calculate deformability (plasticity) profile of molecule based on mechanical
    stiffness matrix (see [EB08]_).

    :arg model: this is an 3-dimensional NMA instance from a :class:`.ANM
        calculations
    :type model: :class:`.ANM`
    :arg pdb: a coordinate set or an object with ``getCoords`` method
    :type pdb: :class:`numpy.ndarray`    
    
    Note: selection can be done usig ``selstr`` and ``pdb_selstr``. ``selstr`` define
    ``model`` selection (used for building :class:`.ANM` model) and ``pdb_selstr`` 
    will be used in VMD program for visualization. 
    
    By default files are saved as *filename* and loaded to VMD program. To change 
    it use ``loadToVMD=False``.
     
    Mean value of mechanical stiffness for molecule can be found in occupancy column
    in PDB file.
    """
    
    pdb = pdb.select(pdb_selstr)
    coords = pdb.select(selstr)
    meanSiff = np.mean(model.getStiffness(), axis=0)
    
    out_mean = open(filename+'_mean.txt','w')   # mean value of Kij for each residue
    for nr_i, i in enumerate(meanSiff):
        out_mean.write("{} {}\n".format(nr_i, i))
    out_mean.close()
    
    from collections import Counter
    aa_counter = Counter(pdb.getResindices()) 
    
    meanStiff_all = []        
    for i in range(coords.numAtoms()):
         meanStiff_all.extend(aa_counter.values()[i]*[round(meanSiff[i], 2)])
        
    from prody.proteins import writePDB    
    kw = {'occupancy': meanStiff_all}
    writePDB(filename, pdb, **kw)                
    LOGGER.info('PDB file with deformability profile has been saved.')
    LOGGER.info('Creating TCL file.')
    out_tcl = open(filename+'.tcl','w')
    out_tcl.write('menu files off \nmenu files on\ndisplay resetview \nmol addrep 0 \ndisplay resetview \n')
    out_tcl.write('mol new {./'+filename+'.pdb} type {pdb} first 0 last -1 step 1 waitfor 1 \n')
    out_tcl.write('animate style Loop \nmenu files off \ndisplay projection Orthographic \n')
    out_tcl.write('display depthcue off \ndisplay rendermode GLSL \naxes location Off \nmenu color off \n')
    out_tcl.write('menu color on \ncolor Display Background white \nmenu color off \nmenu graphics off \n')
    out_tcl.write('menu graphics on \nmol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0 \n')
    out_tcl.write('mol modmaterial 0 0 Diffuse \nmol modcolor 0 0 Occupancy \n')
    out_tcl.close()

    if (loadToVMD == True):
        from prody import pathVMD
        LOGGER.info('File will be loaded to VMD program.')
        os.system(pathVMD()+" -e "+str(filename)+".tcl")

