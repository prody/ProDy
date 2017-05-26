# -*- coding: utf-8 -*-
"""This module defines classes for principal component analysis (PCA) and
essential dynamics analysis (EDA) calculations."""

import time

import numpy as np

from prody import LOGGER, PY2K
from prody.atomic import Atomic
from prody.ensemble import Ensemble, PDBEnsemble
from prody.trajectory import TrajBase
from prody.utilities import importLA

from .nma import NMA

if PY2K:
    range = xrange

__all__ = ['PCA', 'EDA']


class PCA(NMA):

    """A class for Principal Component Analysis (PCA) of conformational
    ensembles. See examples in :ref:`pca`."""

    def __init__(self, name='Unknown'):

        NMA.__init__(self, name)

    def setCovariance(self, covariance):
        """Set covariance matrix."""

        if not isinstance(covariance, np.ndarray):
            raise TypeError('covariance must be an ndarray')
        elif not (covariance.ndim == 2 and
                  covariance.shape[0] == covariance.shape[1]):
            raise TypeError('covariance must be square matrix')
        self._reset()
        self._cov = covariance
        self._dof = covariance.shape[0]
        self._n_atoms = self._dof / 3
        self._trace = self._cov.trace()

    def buildCovariance(self, coordsets, **kwargs):
        """Build a covariance matrix for *coordsets* using mean coordinates
        as the reference.  *coordsets* argument may be one of the following:

        * :class:`.Atomic`
        * :class:`.Ensemble`
        * :class:`.TrajBase`
        * :class:`numpy.ndarray` with shape ``(n_csets, n_atoms, 3)``

        For ensemble and trajectory objects, ``update_coords=True`` argument
        can be used to set the mean coordinates as the coordinates of the
        object.

        When *coordsets* is a trajectory object, such as :class:`.DCDFile`,
        covariance will be built by superposing frames onto the reference
        coordinate set (see :meth:`.Frame.superpose`).  If frames are already
        aligned, use ``aligned=True`` argument to skip this step.


        .. note::
           If *coordsets* is a :class:`.PDBEnsemble` instance, coordinates are
           treated specially.  Let's say **C**\_ij is the element of the
           covariance matrix that corresponds to atoms *i* and *j*.  This
           super element is divided by number of coordinate sets (PDB models or
           structures) in which both of these atoms are observed together."""

        if not isinstance(coordsets, (Ensemble, Atomic, TrajBase, np.ndarray)):
            raise TypeError('coordsets must be an Ensemble, Atomic, Numpy '
                            'array instance')
        LOGGER.timeit('_prody_pca')
        mean = None
        weights = None
        ensemble = None
        if isinstance(coordsets, np.ndarray):
            if (coordsets.ndim != 3 or coordsets.shape[2] != 3 or
                    coordsets.dtype not in (np.float32, float)):
                raise ValueError('coordsets is not a valid coordinate array')
        elif isinstance(coordsets, Atomic):
            coordsets = coordsets._getCoordsets()
        elif isinstance(coordsets, Ensemble):
            ensemble = coordsets
            if isinstance(coordsets, PDBEnsemble):
                weights = coordsets.getWeights() > 0
            coordsets = coordsets._getCoordsets()

        update_coords = bool(kwargs.get('update_coords', False))

        if isinstance(coordsets, TrajBase):
            nfi = coordsets.nextIndex()
            coordsets.reset()
            n_atoms = coordsets.numSelected()
            dof = n_atoms * 3
            cov = np.zeros((dof, dof))
            #mean = coordsets._getCoords().flatten()
            n_confs = 0
            n_frames = len(coordsets)
            LOGGER.info('Covariance will be calculated using {0} frames.'
                        .format(n_frames))
            coordsum = np.zeros(dof)
            LOGGER.progress('Building covariance', n_frames, '_prody_pca')
            align = not kwargs.get('aligned', False)
            for frame in coordsets:
                if align:
                    frame.superpose()
                coords = frame._getCoords().flatten()
                coordsum += coords
                cov += np.outer(coords, coords)
                n_confs += 1
                LOGGER.update(n_confs, '_prody_pca')
            LOGGER.clear()
            cov /= n_confs
            coordsum /= n_confs
            mean = coordsum
            cov -= np.outer(coordsum, coordsum)
            coordsets.goto(nfi)
            self._cov = cov
            if update_coords:
                coordsets.setCoords(mean.reshape((n_atoms, 3)))
        else:
            n_confs = coordsets.shape[0]
            if n_confs < 3:
                raise ValueError('coordsets must have more than 3 coordinate '
                                 'sets')
            n_atoms = coordsets.shape[1]
            if n_atoms < 3:
                raise ValueError('coordsets must have more than 3 atoms')
            dof = n_atoms * 3
            LOGGER.info('Covariance is calculated using {0} coordinate sets.'
                        .format(len(coordsets)))
            s = (n_confs, dof)
            if weights is None:
                if coordsets.dtype == float:
                    self._cov = np.cov(coordsets.reshape((n_confs, dof)).T,
                                       bias=1)
                else:
                    cov = np.zeros((dof, dof))
                    coordsets = coordsets.reshape((n_confs, dof))
                    mean = coordsets.mean(0)
                    LOGGER.progress('Building covariance', n_confs,
                                    '_prody_pca')
                    for i, coords in enumerate(coordsets.reshape(s)):
                        deviations = coords - mean
                        cov += np.outer(deviations, deviations)
                        LOGGER.update(n_confs, '_prody_pca')
                    LOGGER.clear()
                    cov /= n_confs
                    self._cov = cov
            else:
                # PDB ensemble case
                mean = np.zeros((n_atoms, 3))
                for i, coords in enumerate(coordsets):
                    mean += coords * weights[i]
                mean /= weights.sum(0)
                d_xyz = ((coordsets - mean) * weights).reshape(s)
                divide_by = weights.astype(float).repeat(3, axis=2).reshape(s)
                self._cov = np.dot(d_xyz.T, d_xyz) / np.dot(divide_by.T,
                                                            divide_by)
            if update_coords and ensemble is not None:
                if mean is None:
                    mean = coordsets.mean(0)
                ensemble.setCoords(mean)

        self._trace = self._cov.trace()
        self._dof = dof
        self._n_atoms = n_atoms
        LOGGER.report('Covariance matrix calculated in %2fs.', '_prody_pca')

    def calcModes(self, n_modes=20, turbo=True):
        """Calculate principal (or essential) modes.  This method uses
        :func:`scipy.linalg.eigh`, or :func:`numpy.linalg.eigh`, function
        to diagonalize the covariance matrix.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate,
            default is 20,
            if **None** or 'all' is given, all modes will be calculated
        :type n_modes: int

        :arg turbo: when available, use a memory intensive but faster way to
            calculate modes, default is **True**
        :type turbo: bool"""

        linalg = importLA()
        if self._cov is None:
            raise ValueError('covariance matrix is not built or set')
        start = time.time()
        dof = self._dof
        if linalg.__package__.startswith('scipy'):
            if n_modes is None or 'all':
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
            if n_modes is not None:
                LOGGER.info('Scipy is not found, all modes are calculated.')
            values, vectors = linalg.eigh(self._cov)
        # Order by descending SV
        revert = list(range(len(values)-1, -1, -1))
        values = values[revert]
        vectors = vectors[:, revert]
        which = values > 1e-8
        self._eigvals = values[which]
        self._array = vectors[:, which]
        self._vars = self._eigvals
        self._n_modes = len(self._eigvals)
        LOGGER.debug('{0} modes were calculated in {1:.2f}s.'
                     .format(self._n_modes, time.time()-start))

    def performSVD(self, coordsets):
        """Calculate principal modes using singular value decomposition (SVD).
        *coordsets* argument may be a :class:`.Atomic`, :class:`.Ensemble`,
        or :class:`numpy.ndarray` instance.  If *coordsets* is a numpy array,
        its shape must be ``(n_csets, n_atoms, 3)``.  Note that coordinate
        sets must be aligned prior to SVD calculations.

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
            if (coordsets.ndim != 3 or coordsets.shape[2] != 3 or
                    coordsets.dtype not in (np.float32, float)):
                raise ValueError('coordsets is not a valid coordinate array')
            deviations = coordsets - coordsets.mean(0)
        else:
            if isinstance(coordsets, Ensemble):
                deviations = coordsets.getDeviations()
            elif isinstance(coordsets, Atomic):
                deviations = (coordsets._getCoordsets() -
                              coordsets._getCoords())

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
        LOGGER.debug('{0} modes were calculated in {1:.2f}s.'
                     .format(self._n_modes, time.time()-start))

    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add eigen *vector* and eigen *value* pair(s) to the instance.
        If eigen *value* is omitted, it will be set to 1.  Eigenvalues
        are set as variances."""

        NMA.addEigenpair(self, eigenvector, eigenvalue)
        self._vars = self._eigvals

    def setEigens(self, vectors, values=None):
        """Set eigen *vectors* and eigen *values*.  If eigen *values* are
        omitted, they will be set to 1.  Eigenvalues are set as variances."""

        NMA.setEigens(self, vectors, values)
        self._vars = self._eigvals


class EDA(PCA):

    """A class for Essential Dynamics Analysis (EDA) [AA93]_.
    See examples in :ref:`eda`.

    .. [AA93] Amadei A, Linssen AB, Berendsen HJ. Essential dynamics of
       proteins. *Proteins* **1993** 17(4):412-25."""

    pass
