# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
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

"""This module defines a class and a function for Gaussian network model
(GNM) calculations.

""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import time
from types import FunctionType 

import numpy as np

from prody.proteins import parsePDB, Atomic, AtomGroup
from prody.measure import getKDTree
from prody.tools import checkCoords, importLA

from nma import NMA
from gamma import Gamma

__all__ = ['GNM', 'calcGNM']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

ZERO = 1e-6

class GNMBase(NMA):

    """Class for Gaussian Network Model analysis of proteins."""

    def __init__(self, name='Unknown'):
        NMA.__init__(self, name)
        self._is3d = False
        self._cutoff = None
        self._kirchhoff = None
        self._gamma = None
        
    def __repr__(self):
        return '<GNM: {0:s} ({1:d} modes, {2:d} nodes)>'.format(
                                    self._title, self.__len__(), self._n_atoms)
    
    def __str__(self):
        return 'GNM {0:s}'.format(self._title)
    
    def _reset(self):
        NMA._reset(self)
        self._cutoff = None
        self._gamma = None
        self._kirchhoff = None
        self._is3d = False
    
    def getCutoff(self):
        """Return cutoff distance."""
        
        return self._cutoff
    
    def getGamma(self):
        """Return spring constant (or the gamma function or :class:`Gamma`
        instance)."""
        
        return self._gamma

    def getKirchhoff(self):
        """Return a copy of the Kirchhoff matrix."""
        
        if self._kirchhoff is None: return None
        return self._kirchhoff.copy()

    def _getKirchhoff(self):
        """Return the Kirchhoff matrix."""
        
        return self._kirchhoff

def checkENMParameters(cutoff, gamma):
    """Check type and values of *cutoff* and *gamma*."""

    if not isinstance(cutoff, (float, int)):
        raise TypeError('cutoff must be a float or an integer')
    elif cutoff < 4:
        raise ValueError('cutoff must be greater or equal to 4')
    if isinstance(gamma, Gamma):
        gamma_func = gamma.gamma
    elif isinstance(gamma, FunctionType):
        gamma_func = gamma
    else:
        if not isinstance(gamma, (float, int)):
            raise TypeError('gamma must be a float, an integer, derived '
                             'from Gamma, or a function')
        elif gamma <= 0:
            raise ValueError('gamma must be greater than 0')
        gamma = float(gamma)
        gamma_func = lambda dist2, i, j: gamma 
    return cutoff, gamma, gamma_func 

class GNM(GNMBase):
    
    """A class for Gaussian Network Model (GNM) analysis of proteins 
    ([IB97]_, [TH97]_).
    
    |example| See example :ref:`gnm`.
    
    """
    
    def setKirchhoff(self, kirchhoff):
        """Set Kirchhoff matrix."""
        
        if not isinstance(kirchhoff, np.ndarray):
            raise TypeError('kirchhoff must be a Numpy array')
        elif not kirchhoff.ndim == 2 or \
                 kirchhoff.shape[0] != kirchhoff.shape[1]:
            raise ValueError('kirchhoff must be a square matrix')
        elif kirchhoff.dtype != float:
            try:
                kirchhoff = kirchhoff.astype(float)
            except:
                raise ValueError('kirchhoff.dtype must be float')
                
        self._reset()
        self._kirchhoff = kirchhoff
        self._n_atoms = kirchhoff.shape[0]
        self._dof = kirchhoff.shape[0]
    
    def buildKirchhoff(self, coords, cutoff=10., gamma=1., **kwargs):
        """Build Kirchhoff matrix for given coordinate set.
        
        :arg coords: a coordinate set or anything with getCoordinates method
        :type coords: :class:`~numpy.ndarray` or :class:`~prody.atomic.Atomic`
        
        :arg cutoff: cutoff distance (Å) for pairwise interactions
            default is 10.0 Å, , minimum is 4.0 Å
        :type cutoff: float
        
        :arg gamma: spring constant, default is 1.0
        :type gamma: float
        
        :arg sparse: Elect to use sparse matrices. Default is ``False``. If 
            Scipy is not found, :class:`ImportError` is raised.
        :type sparse: bool
        
        .. versionchanged:: 0.6
            Instances of :class:`Gamma` classes and custom functions are
            accepted as *gamma* argument.        

        .. versionchanged:: 0.7.3
           When Scipy is available, user can select to use sparse matrices for
           efficient usage of memory at the cost of computation speed."""
        
        slow = kwargs.get('slow', False)
        try:
            from KDTree import KDTree
        except ImportError:
            KDTree = False
        if not slow and not KDTree: 
            LOGGER.info('Using a slower method for building the Kirchhoff '
                         'matrix.')
        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoords()
            except AttributeError:
                raise TypeError('coords must be a Numpy array or must have '
                                'getCoordinates attribute')
        coords = checkCoords(coords, 'coords')
        cutoff, g, gamma = checkENMParameters(cutoff, gamma)
        self._reset()
        self._cutoff = cutoff
        self._gamma = g
                    
        n_atoms = coords.shape[0]
        start = time.time()
        if kwargs.get('sparse', False):
            try:
                from scipy import sparse as scipy_sparse
            except ImportError:    
                raise ImportError('failed to import scipy.sparse, which  is '
                                  'required for sparse matrix calculations')
            kirchhoff = scipy_sparse.lil_matrix((n_atoms, n_atoms))
        else:
            kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        
        if not slow and KDTree:
            kdtree = getKDTree(coords) 
            kdtree.all_search(cutoff)
            radii = kdtree.all_get_radii()
            r = 0
            for i, j in kdtree.all_get_indices():
                g = gamma(radii[r]**2, i, j)
                kirchhoff[i, j] = -g
                kirchhoff[j, i] = -g
                kirchhoff[i, i] = kirchhoff[i, i] + g
                kirchhoff[j, j] = kirchhoff[j, j] + g
                r += 1
        else:
            cutoff2 = cutoff * cutoff
            for i in range(n_atoms):
                xyz_i = coords[i, :]
                for j in range(i+1, n_atoms):
                    i2j = coords[j, :] - xyz_i
                    dist2 = np.dot(i2j, i2j)
                    if dist2 > cutoff2:
                        continue             
                    g = gamma(dist2, i, j)
                    kirchhoff[i, j] = -g
                    kirchhoff[j, i] = -g
                    kirchhoff[i, i] = kirchhoff[i, i] + g
                    kirchhoff[j, j] = kirchhoff[j, j] + g
            
        LOGGER.debug('Kirchhoff was built in {0:.2f}s.'
                     .format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._n_atoms = n_atoms
        self._dof = n_atoms
        
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh` 
        function to diagonalize the Kirchhoff matrix. When Scipy is not found, 
        :func:`numpy.linalg.eigh` is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
              If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        """
        
        if self._kirchhoff is None:
            raise ValueError('Kirchhoff matrix is not built or set')
        assert n_modes is None or isinstance(n_modes, int) and n_modes > 0, \
            'n_modes must be a positive integer'
        assert isinstance(zeros, bool), 'zeros must be a boolean'
        assert isinstance(turbo, bool), 'turbo must be a boolean'
        linalg = importLA()
        start = time.time()
        shift = 0
        if linalg.__package__.startswith('scipy'):
            if n_modes is None:
                eigvals = None
                n_modes = self._dof 
            else:
                if n_modes >= self._dof:
                    eigvals = None
                    n_modes = self._dof
                else: 
                    eigvals = (0, n_modes + shift)
            if eigvals: 
                turbo = False
            if isinstance(self._kirchhoff, np.ndarray):            
                values, vectors = linalg.eigh(self._kirchhoff, turbo=turbo, 
                                              eigvals=eigvals)
            else:
                try:
                    from scipy.sparse import linalg as scipy_sparse_la
                except ImportError:    
                    raise ImportError('failed to import scipy.sparse.linalg, '
                                      'which is required for sparse matrix '
                                      'decomposition')
                try:
                    values, vectors = scipy_sparse_la.eigsh(
                                self._kirchhoff, k=n_modes + 1, which='SA')
                except:
                    values, vectors = scipy_sparse_la.eigen_symmetric(
                                self._kirchhoff, k=n_modes + 1, which='SA')                
        else:
            values, vectors = linalg.eigh(self._kirchhoff)
        n_zeros = sum(values < ZERO)
        if n_zeros < 1: 
            LOGGER.warning('Less than 1 zero eigenvalues are calculated.')
            shift = n_zeros - 1
        elif n_zeros > 1: 
            LOGGER.warning('More than 1 zero eigenvalues are calculated.')
            shift = n_zeros - 1
        if zeros:
            shift = -1
        self._eigvals = values[1+shift:]
        self._vars = 1 / self._eigvals
        self._trace = self._vars.sum()
        self._array = vectors[:, 1+shift:]
        self._n_modes = len(self._eigvals)
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                          ''.format(self._n_modes, time.time()-start))


def calcGNM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, 
            zeros=False):
    """Return a :class:`GNM` instance and atoms used for the calculations.
    By default only alpha carbons are considered, but selection string helps 
    selecting a subset of it.  *pdb* can be :class:`~prody.atomic.Atomic` 
    instance.  
    
    .. versionchanged:: 0.6
       Returns also the :class:`~prody.atomic.Selection` instance."""
    
    if isinstance(pdb, str):
        ag = parsePDB(pdb)
        title = ag.getTitle()
    elif isinstance(pdb, Atomic):
        ag = pdb
        if isinstance(pdb, AtomGroup):
            title = ag.getTitle()
        else: 
            title = ag.getAtomGroup().getTitle()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'
                        .format(type(pdb)))
    gnm = GNM(title)
    sel = ag.select(selstr)
    gnm.buildKirchhoff(sel, cutoff, gamma)
    gnm.calcModes(n_modes)
    return gnm, sel
