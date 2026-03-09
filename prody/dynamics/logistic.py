# -*- coding: utf-8 -*-
"""This module defines a class for logistic regression classification calculations."""

import time

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic
from prody.ensemble import Ensemble
from prody.trajectory import DCDFile
from prody.utilities import isListLike

from .nma import NMA

__all__ = ['LRA']


class LRA(NMA):

    """A class for logistic regression classification of conformational
    ensembles. See examples in :ref:`pca`."""

    def __init__(self, name='Unknown'):
        NMA.__init__(self, name)

    def calcModes(self, coordsets, labels, lasso=True, **kwargs):
        """Calculate logistic regression classification modes between classes.  
        This method uses :class:`sklearn.linear_model.LogisticRegression`
        on coordsets with class labels.

        *coordsets* argument may be one of the following:

        * :class:`.Atomic`
        * :class:`.Ensemble`
        * :class:`.TrajBase`
        * :class:`numpy.ndarray` with shape ``(n_csets, n_atoms, 3)``

        :arg labels: a set of labels for discriminating classes
        :type labels: :class:`~numpy.ndarray`

        :arg lasso: whether to use lasso regression
            This sets penalty='l1', solver='saga', max_iter=50000 if these are not set
            Default **True**
        :type lasso: bool

        :arg n_shuffles: number of random shuffles of labels to assess variability
        :type n_shuffles: int

        Other kwargs for the LogisticRegression class can also be used
        """
        try:
            from sklearn.linear_model import LogisticRegression
        except ImportError:
            raise ImportError("Please install scikit-learn to use LogisticRegression")

        start = time.time()
        self._clear()

        if isinstance(coordsets, np.ndarray):
            if (coordsets.ndim != 3 or coordsets.shape[2] != 3 or
                    coordsets.dtype not in (np.float32, float)):
                raise ValueError('coordsets is not a valid coordinate array')
            self._coordsets = coordsets
        elif isinstance(coordsets, (Atomic, Ensemble)):
            self._coordsets = coordsets._getCoordsets()
        elif isinstance(coordsets, DCDFile):
            self._coordsets = coordsets.getCoordsets()
        else:
            raise TypeError('coordsets should be Atomic, Ensemble or numpy.ndarray, not {0}'
                            .format(type(coordsets)))

        nconfs = self._coordsets.shape[0]
        if not isListLike(labels):
            raise TypeError('labels must be either a list or a numpy.ndarray, not {0}'
                            .format(type(labels)))
        if not isinstance(labels, np.ndarray):
            labels = np.array(labels)
        if labels.ndim != 1 or len(labels) != nconfs:
            raise ValueError('labels should have same number as conformers')
        
        self._n_atoms = self._coordsets.shape[1]
        
        self._coordsets = self._coordsets.reshape(nconfs, -1)
        self._labels = labels

        quiet = kwargs.pop('quiet', False)

        self._n_shuffles = kwargs.pop('n_shuffles', 0)

        if lasso:
            if 'penalty' not in kwargs:
                kwargs['penalty'] ='l1'
            else:
                LOGGER.warn('using provided penalty kwarg instead of l1 from lasso')

            if 'solver' not in kwargs:
                kwargs['solver'] ='saga'
            else:
                LOGGER.warn('using provided solver kwarg instead of saga for lasso')

            if 'max_iter' not in kwargs:
                kwargs['max_iter'] = 50000
            else:
                LOGGER.warn('using provided max_iter kwarg instead of 50000 that works for lasso')

        # The general multi-class case has one array per class
        self._n_modes = len(set(labels))

        if self._n_modes == 2:
            # binary classification just has one
            self._n_modes = 1

        self._lra = LogisticRegression(**kwargs)
        self._projection = self._lra.fit(self._coordsets, self._labels)
        self._array = self._lra.coef_.T/np.linalg.norm(self._lra.coef_)
        self._eigvals = np.ones(self._n_modes)
        self._vars = np.ones(self._n_modes)

        if not quiet:
            if self._n_modes > 1:
                LOGGER.debug('{0} modes were calculated in {1:.2f}s.'
                        .format(self._n_modes, time.time()-start))
            else:
                LOGGER.debug('{0} mode was calculated in {1:.2f}s.'
                        .format(self._n_modes, time.time()-start))

            if self._n_shuffles > 0:
                if self._n_modes > 1 and self._n_shuffles > 1:
                    LOGGER.debug('Calculating {0} modes for {1} shuffles.'
                        .format(self._n_modes, self._n_shuffles))
                elif self._n_modes > 1 and self._n_shuffles == 1:
                    LOGGER.debug('Calculating {0} modes for 1 shuffle.'
                        .format(self._n_modes))
                else:
                    LOGGER.debug('Calculating 1 mode for 1 shuffle.')

        self._shuffled_lras = [LRA('shuffle '+str(n)) for n in range(self._n_shuffles)]
        self._coordsets_reshaped = self._coordsets.reshape(self._coordsets.shape[0], self._n_atoms, -1)

        n = 0
        while n < self._n_shuffles:
            labelsNew = self._labels.copy()
            # use random generator with None, 
            # then fresh, unpredictable entropy will be pulled from the OS
            rng = np.random.default_rng()
            rng.shuffle(labelsNew) # in place

            self._shuffled_lras[n].calcModes(self._coordsets_reshaped,
                                             labelsNew, quiet=True)
            
            if np.allclose(abs(np.dot(self._shuffled_lras[n].getEigvecs()[0],
                                      self.getEigvecs()[0])), 
                           1):
                # LDA has flipped direction as labels match or are exactly flipped
                continue
            
            n += 1

        if self._n_shuffles > 0 and not quiet:
            if self._n_modes > 1:
                LOGGER.debug('{0} modes were calculated with {1} shuffles in {2:.2f}s.'
                        .format(self._n_modes, self._n_shuffles, time.time()-start))
            else:
                LOGGER.debug('{0} mode was calculated with {1} shuffles in {2:.2f}s.'
                        .format(self._n_modes, self._n_shuffles, time.time()-start))

    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add eigen *vector* and eigen *value* pair(s) to the instance.
        If eigen *value* is omitted, it will be set to 1.  Eigenvalues
        are set as variances."""

        NMA.addEigenpair(self, eigenvector, eigenvalue)
        self._vars = self._eigvals

    def setEigens(self, vectors, values=None):
        """Set eigen *vectors* and eigen *values*.  If eigen *values* are
        omitted, they will be set to 1.  Eigenvalues are set as variances."""

        self._clear()
        NMA.setEigens(self, vectors, values)
        self._vars = self._eigvals

    def getShuffledModes(self):
        return self._shuffled_lras.copy()
    
    def getShuffledEigvecs(self):
        return np.array([lda.getEigvecs() for lda in self._shuffled_lras])

    def getShuffledPercentile(self, percentile, take_abs=True):
        shuffles = self.getShuffledEigvecs()
        if take_abs:
            shuffles = abs(shuffles)
        return np.percentile(shuffles, percentile)
