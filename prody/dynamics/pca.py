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

"""This module defines classes for principal component analysis (PCA) and 
essential dynamics analysis (EDA) calculations.

""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import time

import numpy as np

from prody.atomic import Atomic
from prody.ensemble import Ensemble, PDBEnsemble
from prody.trajectory import TrajBase
from prody.tools import importLA

from nma import NMA

__all__ = ['PCA', 'EDA']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

class PCA(NMA):
    
    """A class for Principal Component Analysis (PCA) of conformational 
    ensembles.
    
    |example| See examples in :ref:`pca`.
    """

    def __init__(self, name='Unknown'):
        NMA.__init__(self, name)
    
    def __repr__(self):
        return '<PCA: {0:s} ({1:d} modes, {2:d} atoms)>'.format(
                self._title, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'PCA {0:s}'.format(self._title)
    
    def setCovariance(self, covariance):
        """Set covariance matrix."""
        
        if not isinstance(covariance, np.ndarray):
            raise TypeError('covariance must be an ndarray')
        elif not (covariance.ndim == 2 and covariance.shape[0] == covariance.shape[1]):
            raise TypeError('covariance must be square matrix')
        self._reset()
        self._cov = covariance
        self._dof = covariance.shape[0]
        self._n_atoms = self._dof / 3
        self._trace = self._cov.trace()

    def buildCovariance(self, coordsets):
        """Build a weighted covariance matrix for *coordsets*.  *coordsets* 
        argument may be an instance of one of the following:
          
        * :class:`~.Atomic`
        * :class:`~.Ensemble`
        * :class:`~.PDBEnsemble`
        * :class:`~.Trajectory`
        * :class:`numpy.ndarray`

       :class:`~numpy.ndarray` and :class:`~.TrajBase` 
       instances are accepted as *coordsets* argument.

        A NumPy array passed as *coordsets* argument must have the shape 
        (n_coordsets, n_atoms, 3).
        
        When a trajectory instance is passed as *coordsets* argument, 
        covariance will be built by aligning frames to the reference 
        coordinates.  Reference coordinates will be considered as the 
        average of coordinate sets in the trajectory. 
         
        .. note::        
           If *coordsets* is a :class:`~prody.ensemble.PDBEnsemble` instance,
           coordinates are treated specially.  Let's say **C**\_ij is the 
           super element of the covariance matrix that corresponds to atoms 
           *i* and *j*.  This super element is divided by number of coordinate
           sets (PDB structures) in which both of these atoms are observed 
           together."""
        
        if not isinstance(coordsets, (Ensemble, Atomic, TrajBase, np.ndarray)):
            raise TypeError('coordsets must be an Ensemble, Atomic, Numpy '
                            'array instance')
        LOGGER.timeit()
        weights = None
        if isinstance(coordsets, np.ndarray): 
            if coordsets.ndim != 3 or coordsets.shape[2] != 3 or \
                coordsets.dtype not in (np.float32, float):
                raise ValueError('coordsets is not a valid coordinate array')
        elif isinstance(coordsets, Atomic):
            coordsets = coordsets._getCoordsets()
        elif isinstance(coordsets, Ensemble):
            if isinstance(coordsets, PDBEnsemble):
                weights = coordsets.getWeights() > 0
            coordsets = coordsets._getCoordsets()
        
        if isinstance(coordsets, TrajBase):
            nfi = coordsets.getNextIndex()
            coordsets.reset()
            n_atoms = coordsets.numSelected()
            dof = n_atoms * 3
            cov = np.zeros((dof, dof))
            mean = coordsets._getCoords().flatten()
            n_confs = 0
            n_frames = len(coordsets)
            LOGGER.info('Covariance will be calculated using {0:d} frames.'
                            .format(n_frames))
            coordsum = np.zeros(dof)
            LOGGER.progress('Calculating covariance', n_frames)
            for frame in coordsets:
                frame.superpose()
                coords = frame._getCoords().flatten()
                coordsum += coords
                cov += np.outer(coords, coords)
                n_confs += 1
                LOGGER.update(n_confs)
            LOGGER.clear()
            cov /= n_confs
            coordsum /= n_confs
            cov -= np.outer(coordsum, coordsum)
            coordsets.goto(nfi)
            self._cov = cov
        else:
            n_confs = coordsets.shape[0]
            if n_confs < 3:
                raise ValueError('coordsets must have more than 3 coordinate '
                                 'sets')
            n_atoms = coordsets.shape[1]
            if n_atoms < 3:
                raise ValueError('coordsets must have more than 3 atoms')
            dof = n_atoms * 3
            LOGGER.info('Covariance is calculated using {0:d} coordinate sets.'
                            .format(len(coordsets)))
            if weights is None:
                if coordsets.dtype == float:
                    self._cov = np.cov(coordsets.reshape((n_confs, dof)).T, 
                                       bias=1)
                else:
                    cov = np.zeros((dof, dof))
                    coordsets = coordsets.reshape((n_confs, dof))
                    mean = coordsets.mean(0)
                    LOGGER.progress('Building covariance', n_confs)
                    for i, coords in enumerate(
                                            coordsets.reshape((n_confs, dof))):
                        deviations = coords - mean
                        cov += np.outer(deviations, deviations)
                        LOGGER.update(n_confs)
                    LOGGER.clear()
                    cov /= n_confs 
                    self._cov = cov
            else:
                # PDB ensemble case
                mean = np.zeros((n_atoms, 3))
                for i, coords in enumerate(coordsets):
                    mean += coords * weights[i]
                mean /= weights.sum(0)
                d_xyz = ((coordsets - mean) * weights).reshape((n_confs, dof))
                divide_by = weights.astype(float).repeat(3, 
                                                axis=2).reshape((n_confs, dof))
                self._cov = np.dot(d_xyz.T, d_xyz) / np.dot(divide_by.T, 
                                                            divide_by)
        self._trace = self._cov.trace()
        self._dof = dof
        self._n_atoms = n_atoms
        LOGGER.timing('Covariance matrix was calculated in %2fs.')
        
    def calcModes(self, n_modes=20, turbo=True):
        """Calculate principal (or essential) modes.  This method uses 
        :func:`scipy.linalg.eigh` function to diagonalize covariance matrix. 
        When Scipy is not found, :func:`numpy.linalg.eigh` is used.
        
        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        """
        
        linalg = importLA()
        if self._cov is None:
            raise ProDyException('covariance matrix is not built or set')
        start = time.time()
        dof = self._dof
        if linalg.__package__.startswith('scipy'):        
            if n_modes is None:
                eigvals = None
                n_modes = dof
            else:
                n_modes = int(n_modes)
                if n_modes >= self._dof:
                    eigvals = None
                    n_modes = dof
                else:
                    eigvals = (dof - n_modes, dof - 1)
            values, vectors = linalg.eigh(self._cov, turbo=turbo, 
                                          eigvals=eigvals)
        else:
            values, vectors = linalg.eigh(self._cov)
        # Order by descending SV
        revert = range(len(values)-1, -1, -1)
        values = values[revert]
        vectors = vectors[:, revert]
        which = values > 1e-8
        self._eigvals = values[which]
        self._array = vectors[:, which]
        self._vars = self._eigvals
        self._n_modes = len(self._eigvals)
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                     .format(self._n_modes, time.time()-start))

    def performSVD(self, coordsets):
        """Calculate principal modes using singular value decomposition (SVD).
        *coordsets* argument may be a :class:`~prody.atomic.Atomic`, 
        :class:`~prody.ensemble.Ensemble`, or :class:`numpy.ndarray` instance.
        If *coordsets* is a numpy array it must have the shape 
        ``(n_coordsets, n_atoms, 3)``.  :class:`numpy.ndarray` instances are 
        accepted as *coordsets* argument.
        
        This is a considerably faster way of performing PCA calculations 
        compared to eigenvalue decomposition of covariance matrix, but is
        an approximate method when heterogeneous datasets are analyzed. 
        Covariance method should be preferred over this one for analysis of 
        ensembles with missing atomic data.  See :ref:`pca-xray-calculations`
        example for comparison of results from SVD and covariance methods."""

        linalg = importLA()

        start = time.time()
        if not isinstance(coordsets, (Ensemble, Atomic, np.ndarray)):
            raise TypeError('coordsets must be an Ensemble, Atomic, Numpy '
                            'array instance')
        if isinstance(coordsets, np.ndarray):
            if coordsets.ndim != 3 or coordsets.shape[2] != 3 or \
                coordsets.dtype not in (np.float32, float):
                raise ValueError('coordsets is not a valid coordinate array')
            deviations = coordsets - coordsets.mean(0)
        else:
            if isinstance(coordsets, Ensemble):
                deviations = coordsets.getDeviations()
            elif isinstance(coordsets, Atomic):
                deviations = coordsets._getCoordsets() - \
                             coordsets._getCoords()

        n_confs = deviations.shape[0]
        if n_confs < 3:
            raise ValueError('coordsets must have more than 3 coordinate sets')
        n_atoms = deviations.shape[1]
        if n_atoms < 3:
            raise ValueError('coordsets must have more than 3 atoms')

        dof = n_atoms * 3        
        deviations = deviations.reshape((n_confs, dof)).T

        vectors, values, self._temp = linalg.svd(deviations, 
                                                 full_matrices=False)
        values = (values ** 2) / n_confs
        self._dof = dof
        self._n_atoms = n_atoms
        which = values > 1e-18
        self._eigvals = values[which]
        self._array = vectors[:, which]
        self._vars = self._eigvals
        self._trace = self._vars.sum()
        self._n_modes = len(self._eigvals)
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                         .format(self._n_modes, time.time()-start))
        
    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add *eigenvector* and *eigenvalue* pair to :class:`NMA` instance.
        If *eigenvalue* is not given, it will be set to 1.  Eigenvalue is also 
        set as the variance."""

        NMA.addEigenpair(self, eigenvector, eigenvalue)
        self._vars = self._eigvals.copy()


    def setEigens(self, vectors, values=None):
        """Set eigenvectors and eigenvalues.
        
        :arg vectors: eigenvectors
        :type vectors: numpy.ndarray
        
        :arg values: Eigenvalues. When ``None`` is passed (default value), 
            all eigenvalues will be set to ``1``.
        :type values: numpy.ndarray
        
        For M modes and N atoms, *vectors* must have shape ``(3*N, M)``
        and values must have shape ``(M,)``.  Eigenvalues are also set as the 
        variances."""
        
        NMA.setEigens(self, vectors, values)
        self._vars = self._eigvals.copy()

class EDA(PCA):
    
    """A class for Essential Dynamics Analysis (EDA) [AA93]_.
    
    |example| See examples in :ref:`eda`."""

    def __repr__(self):
        return '<EDA: {0:s} ({1:d} modes, {2:d} atoms)>'.format(
                self._title, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'EDA {0:s}'.format(self._title)



    
