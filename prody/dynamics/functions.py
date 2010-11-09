# ProDy: A Python Package for Protein Structural Dynamics Analysis
# 
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>

""":mod:`functions` module defines supporting functions for dynamics modules."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import os.path
import time
import os
import gzip

import numpy as np

import prody
from prody import ProDyLogger as LOGGER
from .nma import ZERO
from .nma import *

VMDPATH = '/usr/local/bin/vmd'

__all__ = ['getANM', 'getGNM', 'writeNMD', 
           'viewNMDinVMD', 'setVMDpath', 'getVMDpath',
           'getProjection', 'getSumOfWeights',
           'getOverlap', 'reduceModel', 'printOverlapTable',
           'writeModes', 'writeArray',
           'getSqFlucts', 'getCrossCorrelations',
           'getSubspaceOverlap', 'getCumulativeOverlap',
           'writeOverlapTable', 'getCovariance'
           ]

def getVMDpath():
    """Returns path to the VMD executable."""
    return VMDPATH

def setVMDpath(path):
    """Set the path to VMD executable."""
    if not path.startswith('vmd') and not os.path.isfile(path):
        LOGGER.warning('{0:s} may not be a valid path for VMD executable.')
    VMDPATH = path

def writeNMD(filename, modes, atoms):
    """Writes an NMD file for given *modes* and includes applicable data from 
    *atoms*.
    
    Returns *filename*, if file is successfully written.
    
    This function skips modes with zero eigenvalues.    
    """
    if not isinstance(modes, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(modes)))
    if modes.getNumOfAtoms() != atoms.getNumOfAtoms():
        raise Exception('number of atoms do not match')
    out = open(filename, 'w')
    
    #out.write('#!{0:s} -e\n'.format(VMDPATH))
    out.write('nmwiz_load {0:s}\n'.format(os.path.abspath(filename)))
    out.write('name {0:s}\n'.format(modes.getName()))
    try:
        coords = atoms.getCoordinates()
    except:
        raise RuntimeError('coordinates could not be retrived from atoms instance')
    if coords is None:
        raise RuntimeError('coordinates could not be retrived from atoms instance')
    
    try:
        data = atoms.getAtomNames()
        if data is not None:
            out.write('atomnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResidueNames()
        if data is not None:
            out.write('resnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResidueNumbers()
        if data is not None:
            out.write('resids {0:s}\n'.format(' '.join(data.astype('|S5'))))
    except:
        pass
    try:
        data = atoms.getChainIdentifiers()
        if data is not None:
            out.write('chainids {0:s}\n'.format(' '.join(data)))
    except:
        pass
    
    out.write('coordinates {0:s}\n'.format(' '.join(['{0:.3f}'.format(x) for x in coords.flatten()])))
    
    count = 0
    if isinstance(modes, Vector):
        out.write('mode {0:s}\n'.format(' '.join(['{0:.3f}'.format(x) for x in modes.getArray()])))
        count += 1
    else:
        for mode in modes:
            if mode.getEigenvalue() < ZERO:
                continue
            out.write('mode {0:d} {1:.2f} {2:s}\n'.format(mode.getIndex()+1, mode.getVariance()**0.5, ' '.join(['{0:.3f}'.format(x) for x in mode.getArray()])))
            count += 1
    if count == 0:
        LOGGER.warning('No normal mode data was written. Given modes might have 0 eigenvalues.')
    out.close() 
    return filename  

def viewNMDinVMD(filename):
    """Start VMD in the current Python session and load NMD data."""
    os.system('{0:s} -e {1:s}'.format(VMDPATH, os.path.abspath(filename)))
    
def getANM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, zeros=False):
    """Returns an ANM instance for given PDB identifier or atom data.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`prody.proteins.AtomGroup`, :class:`prody.proteins.Selection`,  
    or :class:`prody.proteins.Chain` instance.  
    
    """
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        name = ag.getName()
    elif isinstance(pdb, (prody.AtomGroup, prody.Selection, prody.Chain, prody.AtomMap)):
        ag = pdb
        if isinstance(pdb, prody.AtomGroup):
            name = ag.getName()
        else: 
            name = ag.getAtomGroup().getName()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'.format(type(pdb)))
    anm = prody.ANM(name + ' ANM')
    anm.buildHessian(ag.select(selstr), cutoff, gamma)
    anm.calcModes(n_modes)
    return anm

def getGNM(pdb, selstr='all', cutoff=15., gamma=1., n_modes=20, zeros=False):
    """Returns an GNM instance for given PDB identifier or atom data.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`prody.proteins.AtomGroup`, :class:`prody.proteins.Selection`,  
    or :class:`prody.proteins.Chain` instance.  
    
    """
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        name = ag.getName()
    elif isinstance(pdb, (prody.AtomGroup, prody.Selection, prody.Chain, prody.AtomMap)):
        ag = pdb
        if isinstance(pdb, prody.AtomGroup):
            name = ag.getName()
        else: 
            name = ag.getAtomGroup().getName()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'.format(type(pdb)))
    gnm = prody.GNM(name + ' GNM')
    gnm.buildKirchhoff(ag.select(selstr), cutoff, gamma)
    gnm.calcModes(n_modes)
    return gnm

def getProjection(ensemble, modes):
    """Returns projection of conformational deviations onto given modes.

    For K conformations and M modes, a (K,M) matrix is returned.     
                   
    """
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    if not modes.is3d(): 
        raise ValueError('modes must be 3-dimensional')
    if ensemble.getNumOfAtoms() != modes.getNumOfAtoms():
        raise ValueError('number of atoms are not the same')
    deviations = ensemble.getDeviations()
    if deviations.ndim == 3:
        deviations = deviations.reshape((deviations.shape[0], 
                                         deviations.shape[1] * 3))
    else:
        deviations = deviations.reshape((1, deviations.shape[0] * 3))
    return np.dot(deviations, modes.getArray())


def getSumOfWeights(ensemble):
    """Returns sum of weights from an ensemble.
    
    Weights are summed for each atom over conformations in the ensemble.
    Size of the plotted array will be equal to the number of atoms.
    
    When analyzing an ensemble of X-ray structures, this function can be used 
    to see how many times a residue is resolved.
    
    """
    
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble instance')
    
    weights = ensemble.getWeights()
    
    if weights is None:
        return None
    
    return weights.sum(0)
    
    
def getOverlap(rows, cols):
    """Returns overlap (or correlation) between two sets of modes (*rows* and *cols*).
    
    Returns a matrix whose rows correspond to modes passed as 
    *rows* argument, and columns correspond to those passed as *cols* 
    argument.
    
    """
    if not isinstance(rows, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('rows must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(rows)))
    if not isinstance(rows, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('cols must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(cols)))
    
    if rows.getNumOfDegOfFreedom() != cols.getNumOfDegOfFreedom(): 
        raise ValueError('number of defrees of freedom of rows and cols must be the same')
        
    return np.dot(rows.getArray().T, cols.getArray())

def printOverlapTable(rows, cols):
    """Print table of overlaps (correlations) between two sets of modes.
    
    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the printed table.

    This function may be used to take a quick look into mode correspondences 
    between two models.

    """
    print getOverlapTable(rows, cols)

def writeOverlapTable(filename, rows, cols):
    """Write table of overlaps (correlations) between two sets of modes to a file.

    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the overlap table.
    
    See also :func:`printOverlapTable`.

    """
    out = open(filename, 'w')
    out.write(getOverlapTable(rows, cols))
    out.close()
    return filename
    
def getOverlapTable(rows, cols):
    overlap = getOverlap(rows, cols)
    rname = rows.getName()
    cname = cols.getName()
    if isinstance(rows, Mode):
        rids = [rows.getIndex()]
    elif isinstance(rows, NMA): 
        rids = range(len(rows))        
    elif isinstance(rows, ModeSet): 
        rids = rows.getIndices()
    else:
        rids = [1]
    rlen = len(rids)
    if isinstance(cols, Mode):
        cids = [cols.getIndex()]
    elif isinstance(cols, NMA): 
        cids = range(len(cols))    
    elif isinstance(cols, ModeSet):
        cids = cols.getIndices()
    else:
        cids = [1]
    clen = len(cids)

    table = 'Overlap/Correlation Table\n'
    table += ' '*(len(rname)+5) + cname.center(clen*7)+'\n'  
    table += ' '*(len(rname)+5)
    for j in range(clen):
        table += ('#{0}'.format(cids[j]+1)).center(7)
    table += '\n'
    for i in  range(rlen):
        table += rname + (' #{0}'.format(rids[i]+1)).ljust(5)
        for j in  range(clen):
            if overlap[i, j] < 0: 
                minplus = '-'
            else: 
                minplus = '+'
            table += (minplus+'{0:-.2f}').format(abs(overlap[i, j])).center(7)
        table += '\n'
    return table

def reduceModel(model, atoms, selstr):
    """Reduce dynamics model to a subset of *atoms* matching a selection *selstr*.

    Returns a tuple containing reduced model and atom selection.
    
    This function behaves depending on the type of the model.
    
    :arg model: dynamics model
    :type model: :class:`ANM`, :class:`GNM`, or :class:`PCA`
    :arg atoms: atoms that were used to build the model
    :arg selstr: a selection string specifying subset of atoms  

    For ANM and GNM:    
       This function implements [KH00]_. Selected atoms constitute the system 
       and the rest is the environment.
    
    For PCA:
       This function simply takes the sub-covariance matrix for the selected
       atoms.
       
    
    """
    if prody.la is None:
        prody.importScipyLinalg()

    #LOGGER.warning('Implementation of this function is not finalized. Use it with caution.')
    if not isinstance(model, prody.NMA):
        raise TypeError('model must be an NMA instance, not {0:s}'.format(type(model)))
    if not isinstance(atoms, (prody.AtomGroup, prody.AtomSubset, prody.AtomMap)):
        raise TypeError('atoms type is not valid')
    if len(atoms) <= 1:
        raise TypeError('atoms must contain more than 1 atoms')

    if isinstance(model, prody.GNM):
        matrix = model._kirchhoff
    elif isinstance(model, prody.ANM):
        matrix = model._hessian
    elif isinstance(model, prody.PCA):
        matrix = model._cov
    else:
        raise TypeError('model does not have a valid type derived from prody.NMA')

    if matrix is None:
        raise ValueError('model matrix (Hessian/Kirchhoff/Covariance) is not built')

    if isinstance(atoms, prody.AtomGroup):
        indices = np.arange(atoms.getNumOfAtoms())
    else:        
        indices = atoms.getIndices()
    
    selection = atoms.select(selstr)
    if len(selection) == 0:
        LOGGER.warning('selection has 0 atoms')
        return None
    sel = selection.getIndices()
    if len(indices) == len(sel):
        LOGGER.warning('selection results in same number of atoms, model is not reduced')
        return None
    ndim = 1
    if model._is3d:
        ndim = 3
    system = [] 
    other = []
    index = 0
    for i in indices:
        if i in sel:
            system.extend(range(index*ndim, (index+1)*ndim))
        else:
            other.extend(range(index*ndim, (index+1)*ndim))
        index += 1
    ss = matrix[system,:][:,system]
    if isinstance(model, prody.PCA):
        eda = prody.PCA('Reduced '+model.getName())
        eda.setCovariance(ss)
        return eda, selection
    so = matrix[system,:][:,other]
    os = matrix[other,:][:,system]
    oo = matrix[other,:][:,other]
    matrix = ss - np.dot(so, np.dot(prody.la.inv(oo), os))
    
    if isinstance(model, prody.GNM):
        gnm = prody.GNM('Reduced '+model.getName())
        gnm.setKirchhoff(matrix)
        return gnm, selection
    elif isinstance(model, prody.ANM):
        anm = prody.ANM('Reduced '+model.getName())
        anm.setHessian(matrix)
        return anm, selection
    elif isinstance(model, prody.PCA):
        eda = prody.PCA('Reduced '+model.getName())
        eda.setCovariance(matrix)
        return eda, selection

def writeModes(filename, modes, format='g', sep=' ', compressed=False):
    """Write *modes* (eigenvectors) into a plain text file with name *filename*.
    
    See also :func:`writeArray`.
    
    """
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    return writeArray(filename, modes.getArray(), format=format, sep=sep)
    
def writeArray(filename, array, format='g', sep=' ', compressed=False):
    """Write 1-d or 2-d array data into a delimited text file.
    
    If *filename* is ``None``, output will be returned as a formatted string. 
    Default *format* argument is "g", which writes enough many significant 
    digits automatically.  Default seperator (*sep*) argument is white 
    space " ". Optionally, a gzip *compressed* file may be outputted.
    
    If array is written into a file, *filename* will be returned upon
    successful writing. 
     
    """
    if not array.ndim in (1, 2):
        raise ValueError('array must be 1- or 2-dimensional')
    if array.ndim == 2:
        length = array.shape[1]
    else:
        length = 1
        
    if length == 1:
        line = '{' + '0:{0:s}'.format(format) + '}'
    else:
        line = '{' + '0[{0:d}]:{1:s}'.format(0, format) + '}'
    for j in range(1, length):
        line += sep + '{' + '0[{0:d}]:{1:s}'.format(j, format) + '}' 
    line += '\n'
    output = ''
    for row in array:
        output += line.format(row)
    if compressed:
        out = gzip.open(filename + '.gz', 'w')
    else:
        out = open(filename, 'w')
    out.write(output)
    out.close()
    return filename



def getSqFlucts(modes):
    """Returns sum of square-fluctuations for given set of normal *modes*."""
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if isinstance(modes, Mode):
        modes = [modes]
    square_fluctuations = np.zeros(modes.getNumOfAtoms()) 
    for mode in modes.modes:
        square_fluctuations += mode.getSqFlucts()
    return square_fluctuations
 
def getCrossCorrelations(modes, n_cpu=1):
    """Returns cross-correlations matrix.
    
    For a 3-d model, cross-correlations matrix is an NxN matrix, where N is the 
    number of atoms. Each element of this matrix is the trace of the 
    submatrix corresponding to a pair of atoms.
    
    Covariance matrix may be calculated using all modes or a subset of modes
    of an NMA instance.

    For large systems, calculation of cross-correlations matrix may be time 
    consuming. Optionally, multiple processors may be employed to perform
    calculations by passing ``n_cpu=2`` or more. 

    :arg n_cpu: number of CPUs to use 
    :type n_cpu: int, default is 1
    
    """
    if not isinstance(n_cpu, int):
        raise NMAError('n_cpu must be an integer')
    elif n_cpu < 1:
        raise NMAError('n_cpu must be equal to or greater than 1')
        
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
        
    if modes.is3d():
        model = modes
        if isinstance(modes, (Mode, ModeSet)):
            model = modes._model
            if isinstance(modes, (Mode)):
                indices = [modes.getIndex()]
                n_modes = 1
            else:
                indices = modes.getIndices()
                n_modes = len(modes)
        else:
            n_modes = len(modes)
            indices = np.arange(n_modes)
        array = model._array
        n_atoms = model._n_atoms
        variances = model._vars
        if n_cpu == 1:
            arvar = (array[:, indices]*variances[indices]).T.reshape((n_modes,
                                                                   n_atoms, 3))
            array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
            covariance = np.tensordot(array.transpose(2, 0, 1),
                                      arvar.transpose(0, 2, 1),
                                      axes=([0, 1], [1, 0]))
        else:
            import multiprocessing
            n_cpu = min(multiprocessing.cpu_count(), n_cpu)
            queue = multiprocessing.Queue()
            size = n_modes / n_cpu
            for i in range(n_cpu):
                if n_cpu - i == 1:
                    indices = modes.indices[i*size:]
                else:
                    indices = modes.indices[i*size:(i+1)*size]
                process = multiprocessing.Process(target=_cross_correlations, 
                              args=(queue, n_atoms, array, variances, indices))
                process.start()
            while queue.qsize() < n_cpu:
                time.sleep(0.05)
            covariance = queue.get()
            while queue.qsize() > 0:
                covariance += queue.get()
    else:
        covariance = modes.getCovariance()
    diag = np.power(covariance.diagonal(), 0.5)
    return covariance / np.outer(diag, diag)

def _crossCorrelations(queue, n_atoms, array, variances, indices):
    """Calculate covariance-matrix for a subset of modes."""
    n_modes = len(indices)
    arvar = (array[:, indices] * variances[indices]).T.reshape((n_modes,
                                                                n_atoms, 3))
    array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
    covariance = np.tensordot(array.transpose(2, 0, 1),
                              arvar.transpose(0, 2, 1),
                              axes=([0, 1], [1, 0]))
    queue.put(covariance)


def getCumulativeOverlap(modes1, modes2):
    """Returns cumulative overlap of modes in *modes2* with those in *modes1*.
    
    Returns an array. Elements of the array correspond to modes passed as 
    *modes1* argument. 
        
    """
    overlap = getOverlap(modes1, modes2)
    cumov = np.sqrt(np.power(overlap, 2).sum(axis=overlap.ndim-1))
    return cumov


def getSubspaceOverlap(modes1, modes2):
    """Returns subspace overlap between two sets of modes (*modes1* and *modes2*).
    
    Also known as the root mean square inner product (RMSIP) of essential 
    subspaces [AA99]_.
    
    This function returns a single number.
        
    """
    modes1 = get_dict(modes1)
    modes2 = get_dict(modes2)
    overlap = getOverlap(modes1, modes2)
    rmsip = np.sqrt(np.power(overlap, 2).sum() /
                               len(adict['modes']))
    return rmsip

def getCovariance(modes):
    """Calculate covariance matrix from given modes and return it."""
    return modes.getCovariance()
