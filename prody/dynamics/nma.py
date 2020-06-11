# -*- coding: utf-8 -*-
"""This module defines a class handling normal mode analysis data."""

import numpy as np

from .mode import Mode
from .modeset import ModeSet

from prody import PY2K, LOGGER
from prody.utilities import isListLike

if PY2K:
    range = xrange

__all__ = ['NMA', 'MaskedNMA']

class NMA(object):

    """A class for handling Normal Mode Analysis (NMA) data."""

    def __init__(self, title='Unknown'):

        self._title = str(title).strip()
        self._n_modes = 0
        self._cov = None
        self._n_atoms = 0
        self._dof = 0
        self._array = None      # modes/eigenvectors
        self._eigvals = None
        self._vars = None       # evs for PCA, inverse evs for ENM
        self._trace = None
        self._is3d = True       # is set to false for GNM
        self._indices = None

    def __len__(self):

        return self._n_modes

    def __getitem__(self, index):
        """A list or tuple of integers can be used for indexing."""

        if self._n_modes == 0:
            raise ValueError('{0} modes are not calculated, use '
                             'calcModes() method'.format(str(self)))
        if isinstance(index, slice):
            if (index.stop is not None and index.stop > self.numModes()) or (index.start is not None and index.start > self.numModes()):
                LOGGER.warn('The selection index contains a higher number than the total mode number ({0})'.format(self.numModes()))
            indices = np.arange(*index.indices(len(self)))
            if len(indices) > 0:
                return ModeSet(self, indices)
        elif isinstance(index, (list, tuple, np.ndarray)):
            if len(index) == 1:
                return self._getMode(index[0])
            return ModeSet(self, index)
        try:
            index = int(index)
        except Exception:
            raise IndexError('indices must be int, slice, list, or tuple')
        else:
            return self._getMode(index)

    def __iter__(self):

        for i in range(self._n_modes):
            yield self[i]

    def __repr__(self):

        return '<{0}: {1} ({2} modes; {3} atoms)>'.format(
                self.__class__.__name__, self._title, self._n_modes,
                self._n_atoms)

    def __str__(self):

        return self.__class__.__name__ + ' ' + self._title

    def _reset(self):
        self._n_modes = 0
        self._cov = None

        self._n_atoms = 0
        self._dof = 0

        self._array = None
        self._eigvals = None
        self._vars = None
        self._trace = None

        self._is3d = True

    def _clear(self):
        pass

    def _getMode(self, index):

        if self._n_modes == 0:
            raise ValueError('{0} modes are not calculated, use '
                             'calcModes() method'.format(str(self)))
        if index >= self._n_modes or index < -self._n_modes:
            raise IndexError('{0} contains {1} modes, try 0 <= index < '
                             '{1}'.format(str(self), self._n_modes))
        if index < 0:
            index += self._n_modes
        return Mode(self, index)

    def _getTrace(self):
        """Returns trace, and emit a warning message if trace is calculated
        using eigenvalues of a subset of variances (eigenvalues or inverse
        eigenvalues)."""

        trace = self._trace
        if trace is None:
            if self._vars is None:
                raise ValueError('variances are not set or calculated')
            trace = self._vars.sum()
            diff = self._dof - self._n_modes
            if self._is3d and diff > 6:
                diff = True
            elif diff > 1:
                diff = True
            else:
                diff = False
            if diff:
                from prody import LOGGER
                LOGGER.warn('Total variance for {0} is calculated using '
                            '{1} available modes out of {2} possible.'
                            .format(str(self), self._n_modes, self._dof))
        return trace

    def getModel(self):
        """Returns self."""

        return self

    def is3d(self):
        """Returns **True** if model is 3-dimensional."""

        return self._is3d

    def numAtoms(self):
        """Returns number of atoms."""

        return self._n_atoms

    def numModes(self):
        """Returns number of modes in the instance (not necessarily maximum
        number of possible modes)."""

        return self._n_modes

    def numDOF(self):
        """Returns number of degrees of freedom."""

        return self._dof

    def numEntries(self):
        """Returns number of entries in one eigenvector."""
        
        arr = self._getArray()
        if arr is None:
            return 0
        return arr.shape[0]

    def getTitle(self):
        """Returns title of the model."""

        return self._title

    def setTitle(self, title):
        """Set title of the model."""

        self._title = str(title)

    def getEigvals(self):
        """Returns eigenvalues.  For :class:`.PCA` and :class:`.EDA` models
        built using coordinate data in Å, unit of eigenvalues is |A2|.  For
        :class:`.ANM`, :class:`.GNM`, and :class:`.RTB`, on the other hand,
        eigenvalues are in arbitrary or relative units but they correlate with
        stiffness of the motion along associated eigenvector."""

        if self._eigvals is None: return None
        return self._eigvals.copy()

    def getVariances(self):
        """Returns variances.  For :class:`.PCA` and :class:`.EDA` models
        built using coordinate data in Å, unit of variance is |A2|.  For
        :class:`.ANM`, :class:`.GNM`, and :class:`.RTB`, on the other hand,
        variance is the inverse of the eigenvalue, so it has arbitrary or
        relative units."""

        if self._vars is None: return None
        return self._vars.copy()

    def getArray(self):
        """Returns a copy of eigenvectors array."""

        if self._array is None: return None
        return self._getArray().copy()

    getEigvecs = getArray

    def _getArray(self):
        """Returns eigenvectors array."""

        if self._array is None: return None
        return self._array

    def getCovariance(self):
        """Returns covariance matrix.  If covariance matrix is not set or yet
        calculated, it will be calculated using available modes."""

        if self._cov is None:
            array = self.getArray()
            if array is None:
                return None
            self._cov = np.dot(array, np.dot(np.diag(self._vars), array.T))
        return self._cov

    def calcModes(self):
        """"""

    def addEigenpair(self, vector, value=None):
        """Add eigen *vector* and eigen *value* pair(s) to the instance.
        If eigen *value* is omitted, it will be set to 1.  Inverse
        eigenvalues are set as variances."""

        try:
            ndim, shape = vector.ndim, vector.shape
        except AttributeError:
            raise TypeError('vector must be a numpy array')

        if ndim > 2:
            raise ValueError('vector must be 1- or 2-dimensional')
        elif ndim == 2 and shape[0] < shape[1]:
            raise ValueError('eigenvectors must correspond to vector columns')
        else:
            vector = vector.reshape((shape[0], 1))

        eigval = value

        if eigval is None:
            if ndim == 1:
                eigval = np.ones(1)
            else:
                eigval = np.ones(shape[2])
        else:
            try:
                ndim2, shape2 = eigval.ndim, eigval.shape
            except AttributeError:
                try:
                    eigval = np.array([value], float)
                except Exception:
                    raise TypeError('value must be a number or array')
            else:
                if value.ndim > 1:
                    raise ValueError('value must be a 1-dimensional array')
                elif value.shape[0] != value.shape[0]:
                    raise ValueError('number of eigenvectors and eigenvalues '
                                     'must match')

        if self._array is None:
            return self.setEigens(vector, value)

        if vector.shape[0] != self._array.shape[0]:
            raise ValueError('shape of vector do not match shape of '
                             'existing eigenvectors')
        self._array = np.concatenate((self._array, vector), 1)
        self._eigvals = np.concatenate((self._eigvals, value))
        self._n_modes += shape[1]
        self._vars = 1 / self._eigvals

    def setEigens(self, vectors, values=None):
        """Set eigen *vectors* and eigen *values*.  If eigen *values* are
        omitted, they will be set to 1.  Inverse eigenvalues are set as
        variances."""

        try:
            ndim, shape = vectors.ndim, vectors.shape
        except AttributeError:
            raise TypeError('vectors must be a numpy array')

        if ndim != 2:
            raise ValueError('vectors must be a 2-dimensional array')
        else:
            dof = shape[0]
            if self._is3d:
                n_atoms = dof // 3
            else:
                n_atoms = dof
            if self._n_atoms > 0 and n_atoms != self._n_atoms:
                LOGGER.warn('the number of atoms of the model will be changed to {0}'
                            .format(n_atoms))
            n_modes = vectors.shape[1]
        if values is not None:
            if not isinstance(vectors, np.ndarray):
                raise TypeError('values must be a numpy.ndarray, not {0}'
                                .format(type(vectors)))
            elif values.ndim != 1:
                raise ValueError('values must be a 1-dimensional array')
            else:
                if values.shape[0] != vectors.shape[1]:
                    raise ValueError('number of vectors and values do not match')
        else:
            values = np.ones(n_modes)

        self._array = vectors
        self._eigvals = values
        self._dof = dof
        self._n_atoms = n_atoms
        self._n_modes = n_modes
        self._vars = 1 / values

        self._clear()

    def getIndices(self):
      if self._indices is not None:
        return self._indices
      else:
        return np.arange(self.numModes())

class MaskedNMA(NMA):
    def __init__(self, name='Unknown', mask=False, masked=True):
        if isinstance(name, str):
            super(MaskedNMA, self).__init__(name)
        else:
            model = name
            name = model.getTitle()
            super(MaskedNMA, self).__init__(name)
            self.__dict__.update(model.__dict__)
            
        self.mask = False
        self.masked = masked
        self._maskedarray = None

        if not np.isscalar(mask):
            self.mask = np.array(mask)

    def _isOriginal(self):
        return self.masked or np.isscalar(self.mask)

    def __repr__(self):
        if self._isOriginal():
            n_dummies = 0
        else:
            n_dummies = len(self.mask) - self._n_atoms

        return ('<{0}: {1} ({2} modes; {3} nodes + {4} dummies)>'
                .format(self.__class__.__name__, self._title, self.__len__(), 
                        self._n_atoms, n_dummies))

    def numAtoms(self):
        """Returns number of atoms."""

        if self._isOriginal():
            return self._n_atoms
        else:
            return len(self.mask)

    def _extend(self, arr, axis=None, defval=0):
        mask = self.mask#.copy()
        if self.is3d():
            mask = np.repeat(mask, 3)

        n_true = np.sum(mask)
        N = len(mask)

        if axis is None:
            axes = [i for i in range(arr.ndim)]
        elif not isListLike(axis):
            axes = [axis]
        else:
            axes = axis

        shape = np.array(arr.shape)
        shape[axes] = N

        whole_array = np.empty(shape, dtype=arr.dtype)
        whole_array.fill(defval)

        I = [np.arange(s) for s in shape]
        J = [np.arange(s) for s in arr.shape]

        for ax in axes:
            I[ax] = mask
            J[ax] = np.arange(n_true)

        whole_array[np.ix_(*I)] = arr[np.ix_(*J)]

        return whole_array

    def _getArray(self):
        """Returns eigenvectors array. The function returns 
        a copy of the array if *masked* is **True**."""

        if self._array is None: return None

        if self._isOriginal():
            array = self._array
        else:
            if self._maskedarray is None:
                array = self._maskedarray = self._extend(self._array, axis=0)
            else:
                array = self._maskedarray

        return array

    def setNumAtoms(self, n_atoms):
        """Fixes the tail of the model. If *n_atoms* is greater than the original size 
        (number of nodes), then extra hidden nodes will be added to the model, and if 
        not, the model will be trimmed so that the total number of nodes, including the 
        hidden ones, will be equal to the *n_atoms*. Note that if *masked* is **True**, 
        the *n_atoms* should be the number of *visible* nodes.

        :arg n_atoms: number of atoms of the model 
        :type n_atoms: int
        """
        def _fixLength(vector, length, filled_value=0, axis=0):
            shape = vector.shape
            dim = len(shape)

            if shape[axis] < length:
                dl = length - shape[0]
                pad_width = []
                for i in range(dim):
                    if i == axis:
                        pad_width.append((0, dl))
                    else:
                        pad_width.append((0, 0))
                vector = np.pad(vector, pad_width, 'constant', constant_values=(filled_value,))
            elif shape[axis] > length:
                vector = vector[:length]
            return vector

        trimmed = self.masked
        if trimmed:
            trimmed_length = self.numAtoms()
            self.masked = False
            extra_length = self.numAtoms() - trimmed_length
            n_atoms += extra_length

        if np.isscalar(self.mask):
            self.mask = np.ones(self.numAtoms(), dtype=bool)

        self.mask = _fixLength(self.mask, n_atoms, False)
        self.masked = trimmed
        self._n_atoms = np.sum(self.mask, dtype=int)
        return

    def setEigens(self, vectors, values=None):
        mask = self.mask
        if not self.masked:
            if not np.isscalar(mask):
                if self.is3d():
                    mask = np.repeat(mask, 3)
                try:
                    vectors = vectors[mask, :]
                except IndexError:
                    raise IndexError('size mismatch between vectors (%d) and mask (%d).'
                                     'Try set masked to False or reset mask'%(vectors.shape[0], len(mask)))
        self._maskedarray = None
        super(MaskedNMA, self).setEigens(vectors, values)

    def calcModes(self):
        self._maskedarray = None
        super(MaskedNMA, self).calcModes()

    def _reset(self):
        super(MaskedNMA, self)._reset()
        self._maskedarray = None
