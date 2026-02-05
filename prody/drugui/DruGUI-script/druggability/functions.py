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

"""Auxiliary functions for Druggability Index Analysis class.

Functions:
    
* :func:`pickler`
* :func:`get_histr`
* :func:`single_linkage`
* :func:`convert_dG_to_Kd`
* :func:`convert_Kd_to_dG`
* :func:`format_Kd`
* :func:`get_pdb_lines`
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'
__version__ = '0.5.2'


import numpy as np

import pickle
import gzip
import os.path

from .abase import ABase

__all__ = ['pickler']

GAS_CONSTANT = 0.0019858775                        # kcal/K/mol


def get_histr(array, bins=10, **kwargs):
    """Return histogram as a string for printing purposes.
    
    This function is making use of Numpy :func:`histogram` function.
    
    :arg array: data array
    :type array: array_like
    
    :arg bins: number of bins
    :type bins: int, default is 10
    
    :arg format: number format string
    :type format: str, default is {0:.2f}
    
    :arg orientation: horizontal ("h") or vertical ("v")
    :type orientation: str, default is "h"
    
    :arg line: line symbol
    :type line: str, defaults are "|" for vertical and "-" for horizontal
    
    :arg marker: line end symbol
    :type marker: str, default is "o"
    
    :arg label: axis label
    :type label: str

    :arg title: histogram title
    :type title: str

    """
    if isinstance(array, np.ndarray):
        if array.ndim > 1:
            raise ValueError('array must be 1-dimensional')
    counts, ranges = np.histogram(array, bins=bins)

    maxcount = max(counts)
    format = kwargs.get('format', '{0:.2f}')
    r_width = max([len(format.format(x)) for x in ranges])

    histr = ''
    marker = kwargs.get('marker', 'o')
    if kwargs.get('orientation', 'h')[0] == 'h': 
        counts = counts[-1::-1]
        ranges = ranges[-1::-1]
        if 'title' in kwargs:     
            histr += kwargs['title'].center(maxcount + r_width + 4) + '\n'
        if 'label' in kwargs:
            histr += kwargs['label'] + '\n'
        line = kwargs.get('line', '-')
        histr += ' ' * r_width + ' #' + '-' * maxcount + '-#\n'
        for i, count in enumerate(counts):
            histr += format.format(ranges[i]).strip().rjust(r_width) + ' |'
            histr += line * (count - 1) + marker * (count > 0)
            histr += ' ' * (maxcount - count + 1) + '|\n'
            
        histr += ' ' * r_width + ' #' + '-' * maxcount + '-#\n'
        histr += ' ' * r_width + ' 0'
        for i in range(0, maxcount+1, 5):
            if i > 0:
                histr += str(i).rjust(5)
    else:
        line = kwargs.get('line', '|') 
        r_width += 1
        clen = len(str(maxcount))
        length = clen + 3 + r_width * bins + 2
        if kwargs.has_key('title'):     
            histr += kwargs['title'].center(length) + '\n'

        histr += '#' * clen + ' --' + '-' * r_width * bins + '--' + '\n'
        
        for i in range(maxcount, 0, -1):
            histr += str(i).rjust(clen) + ' | '
            for j in range(bins):
                if counts[j] > i:
                    histr += line.center(r_width)
                elif counts[j] == i:
                    histr += marker.center(r_width)
                else:
                    histr += ' ' * r_width
            histr += ' |' + '\n'
        histr += ' ' * clen + ' --' + '-' * r_width * bins + '--' + '\n'
        histr += ' ' * clen + '   ' + ''.join(
                    [format.format(x).strip().center(r_width) 
                        for x in ranges[:-1]]
                    )  + '\n'
        if kwargs.has_key('label'):
            histr += kwargs['label'].center(length) + '\n'
    
    return histr

def pickler(filename, obj=None, **kwargs):
    """cPickle/uncPickle an object to/from a gzipped file with the given name.
    
    If the unpickled object has "set_logger" method, it will be called.
    
    This function is defined to complement the 
    :meth:`druggability.abase.ABase.pickle`
    method. It is recommended that ProDy objects are unpickled using this 
    function.
    
    """
    if obj is None:
        if not os.path.isfile(filename):
            raise IOError('{0:s} is not found.'.format(filename))        
        if filename.endswith('gz'):
            gzfile = gzip.open(filename)
        else:
            gzfile = open(filename)
        obj = pickle.load(gzfile)
        gzfile.close()
        
        if isinstance(obj, ABase):
            fileloc = os.path.split(os.path.join(os.getcwd(), filename))[0]
            if obj.workdir != fileloc:
                obj.workdir = fileloc
            
            obj.set_logger(**kwargs)
            if obj.__dict__.has_key('name'):
                obj.logger.info('{0:s} has been reloaded.'.format(obj.name))
            else:
                print('{0:s} has been unpickled.'.format(filename))
            
            obj.set_workdir(obj.workdir)
        else:
            print('{0:s} has been unpickled.'.format(filename))
    else:
        try:
            dumps = pickle.dumps(obj)
            gzfile = gzip.open(filename, 'w')
            gzfile.write(dumps)
            gzfile.close()
        except pickle.UnpickleableError:
            raise TypeError('Object is unpickleable. If it has a '
                            '"pickle" method, try using it.')
        obj = None
    return obj

def single_linkage(dist, cutoff):
    """Return clusters determined by single-linkage agglomerative clustering.
    
    For speedup purposes, distance matrix may contain square-distances and
    cutoff distance may be cutoff-square. 

    :arg dist: distance matrix, strictly lower triangular or a symmetric matrix
               with 0s along the diagonal are acceptable forms.
    :type dist: np.ndarray
    :arg cutoff: distance within which two cluster will be merged to make a new 
                 cluster
    :type cutoff: float
        
    :return: an array of cluster assignments, i.e. items in the same cluster 
             will have save cluster number
    :rtype: np.ndarray
        
         
    """
    if not isinstance(dist, np.ndarray):
        raise TypeError('dist must be of type numpy.ndarray')
    elif dist.ndim != 2:
        raise ValueError('dist must be a 2-dimensional array')
    elif dist.shape[0] != dist.shape[1]:
        raise ValueError('dist must be a square matrix')
    
    size = dist.shape[0]
    clust = [set([i]) for i in range(size)]
    for i in range(size-1):        
        which = (dist[i+1:, i] <= cutoff).nonzero()[0] + i+1
        new = clust[i]
        for j in which:
            new = new.union(clust[j])
        for j in new:
            clust[j] = new
    clusters = - np.ones(size, 'int64')
    iclust = 0
    for i in range(size):
        if clusters[i] == -1:
            clusters[list(clust[i])] = iclust
            iclust += 1
    return clusters

def convert_dG_to_Kd(dG, temperature=300.0):
    """Return molar affinity for given binding free energy (kcal/mol)."""
    return np.exp(dG / GAS_CONSTANT / temperature)

def format_Kd(Kd):
    """Return formatted Kd."""
    if Kd < 1e-14:
        return Kd * 1e15, 'fM'
    elif Kd < 1e-11:
        return Kd * 1e12, 'pM'
    elif Kd < 1e-7:
        return Kd * 1e9, 'nM'
    elif Kd < 1e-4:
        return Kd * 1e6, 'uM'
    elif Kd < 1e-1:
        return Kd * 1e3, 'mM'
    else:
        return 0, '--'

def convert_Kd_to_dG(Kd, temperature=300.0):
    """Return binding free energy (kcal/mol) for given molar affinity."""
    return GAS_CONSTANT * temperature * np.log(Kd)        
    
def get_pdb_lines(coordinates, **kwargs):
    """Return PDB lines for a given set of coordinates.
    
    :arg coordinates: Coordinate array, with shape (N,3). N is number of atoms.
    :type coordinates: numpy.ndarray
        
    :keyword atomnames: List of atom names. 
    :type atomnames: array_like, default is ["CA", "CA", ...]

    :keyword resnames: List of residue names.
    :type resnames: array_like, default is ["GLY", "GLY", ...]
    
    :keyword chainids: List of chain ids.
    :type chainids: array_like, default is ["A", "A", ...]
        
    :keyword resids: List of residue numbers. 
    :type resids: array_like, default is [1, 2, 3, ...]
    
    :keyword occupancy: Atomic occupancy values.
    :type occupancy: array_like, default is [1, 1, ...]
    
    :keyword bfactor: Atomic bfactor values.
    :type bfactor: array_like, default is [0, 0, ...]
    
    :keyword connect: PDB conect lines, as a list of tuples of atom indices.
    :type connect: array_like
    
    :keyword nonzero: Boolean option to skip atoms with zero occupancy
    :type nonzer: bool, default is True
        
    """
    if isinstance(coordinates, np.ndarray):
        if coordinates.ndim == 2:
            n_atoms = coordinates.shape[0]
        else:
            raise ValueError('coordinates must be 2-d. Given array is {0:d}-d.'
                            .format(coordinates.ndim))
    else:
        raise TypeError('coordinates must be of type numpy.ndarray.')
    
    atomnames = kwargs.get('atomnames', None)
    if atomnames is None: 
        atomnames = ['CA'] * n_atoms
    elif not isinstance(atomnames, (np.ndarray, tuple, list)):
        raise TypeError('atomnames must be an array like object')
    elif len(atomnames) != n_atoms:
        raise ValueError('length of atomnames must match number of atoms')

    resnames = kwargs.get('resnames', None)
    if resnames is None: 
        resnames = ['GLY'] * n_atoms
    elif not isinstance(resnames, (np.ndarray, tuple, list)):
        raise TypeError('resnames must be an array like object')
    elif len(resnames) != n_atoms:
        raise ValueError('length of resnames must match number of atoms')

    chainids = kwargs.get('chainids', None)
    if chainids is None: 
        chainids = ['A'] * n_atoms
    elif not isinstance(chainids, (np.ndarray, tuple, list)):
        raise TypeError('chainids must be an array like object')
    elif len(chainids) != n_atoms:
        raise ValueError('length of chainids must match number of atoms')

    resids = kwargs.get('resids', None)
    if resids is None: 
        resids = range(1, n_atoms+1)
    elif not isinstance(resids, (np.ndarray, tuple, list)):
        raise TypeError('resids must be an array like object')
    elif len(resids) != n_atoms:
        raise ValueError('length of resids must match number of atoms')

    occupancy = kwargs.get('occupancy', None)
    if occupancy is None:
        occupancy = np.ones(n_atoms)
    elif isinstance(occupancy, np.ndarray):
            bfactor = occupancy.flatten()
    elif not isinstance(occupancy, (tuple, list)):
        raise TypeError('occupancy must be an array like object')
    elif len(occupancy) != n_atoms:
        raise ValueError('length of occupancy must match number of atoms')
            
    bfactor = kwargs.get('bfactor', None)
    if bfactor is None:
        bfactor = np.zeros(n_atoms)
    elif isinstance(bfactor, np.ndarray):
            bfactor = bfactor.flatten()
    elif not isinstance(bfactor, (tuple, list)):
        raise TypeError('bfactor must be an array like object')
    elif len(bfactor) != n_atoms:
        raise ValueError('length of bfactor must match number of atoms')

    nonzero = kwargs.get('nonzero', False)
    if not isinstance(nonzero, bool):
        raise TypeError('nonzero must be of type bool')

    pdb = ''
    for i, xyz in enumerate(coordinates):
        if nonzero and not occupancy[i]:
            continue
        pdb += ('{0:6s}{1:5d}{2:2s}{3:3s}{4:1s}{5:4s}{6:1s}{7:4d}' +
           '{8:4s}{9:8.3f}{10:8.3f}{11:8.3f}{12:6.2f}{13:6.2f}' +
           '{14:6s}{15:4s}{16:2s}\n').format('ATOM  ', i+1, '  ', 
            atomnames[i].ljust(3), ' ', resnames[i].ljust(4), 
            chainids[i], int(resids[i]), '    ', float(xyz[0]), float(xyz[1]), float(xyz[2]), 
            float(occupancy[i]), float(bfactor[i]), '      ', '    ',  
            atomnames[i][0].rjust(2))

    conect = kwargs.get('connect', None)    
    if conect is not None:
        if not isinstance(nonzero, list):
            raise TypeError('connect must be a list')
        for bond in conect:
            if not isinstance(bond, (tuple, list)):
                raise TypeError('items of connect must be array_like')
            elif len(bond) > 1:
                raise ValueError('items of connect must have length 2')
            pdb += 'CONECT{0:5d}{1:5d}\n'.format(bond[0], bond[1])
    return pdb
