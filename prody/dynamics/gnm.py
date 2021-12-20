# -*- coding: utf-8 -*-
"""This module defines a class and a function for Gaussian network model
(GNM) calculations."""

import time
from types import FunctionType

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB
from prody.kdtree import KDTree
from prody.utilities import importLA, checkCoords, solveEig, ZERO

from .nma import NMA, MaskedNMA
from .gamma import Gamma

__all__ = ['GNM', 'MaskedGNM', 'calcGNM']

class GNMBase(NMA):

    """Class for Gaussian Network Model analysis of proteins."""

    def __init__(self, name='Unknown'):

        super(GNMBase, self).__init__(name)
        self._is3d = False
        self._cutoff = None
        self._kirchhoff = None
        self._gamma = None

    def __repr__(self):

        return ('<{0}: {1} ({2} modes; {3} nodes)>'
                .format(self.__class__.__name__, self._title, self.__len__(),
                        self._n_atoms))

    def __str__(self):

        return self.__class__.__name__ + ' ' + self._title

    def _reset(self):

        super(GNMBase, self)._reset()
        self._cutoff = None
        self._gamma = None
        self._kirchhoff = None
        self._is3d = False

    def _clear(self):
        self._trace = None
        self._cov = None
        
    def getCutoff(self):
        """Returns cutoff distance."""

        return self._cutoff

    def getGamma(self):
        """Returns spring constant (or the gamma function or :class:`Gamma`
        instance)."""

        return self._gamma

    def getKirchhoff(self):
        """Returns a copy of the Kirchhoff matrix."""

        if self._kirchhoff is None:
            return None
        return self._getKirchhoff().copy()

    def _getKirchhoff(self):
        """Returns the Kirchhoff matrix."""

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

    See example :ref:`gnm`.

    .. [IB97] Bahar I, Atilgan AR, Erman B. Direct evaluation of thermal
       fluctuations in protein using a single parameter harmonic potential.
       *Folding & Design* **1997** 2:173-181.

    .. [TH97] Haliloglu T, Bahar I, Erman B. Gaussian dynamics of folded
       proteins. *Phys. Rev. Lett.* **1997** 79:3090-3093."""

    def __init__(self, name='Unknown'):
        super(GNM, self).__init__(name)

    def setKirchhoff(self, kirchhoff):
        """Set Kirchhoff matrix."""

        if not isinstance(kirchhoff, np.ndarray):
            raise TypeError('kirchhoff must be a Numpy array')
        elif (not kirchhoff.ndim == 2 or
              kirchhoff.shape[0] != kirchhoff.shape[1]):
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

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray` or :class:`.Atomic`

        :arg cutoff: cutoff distance (Å) for pairwise interactions
            default is 10.0 Å, , minimum is 4.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float

        :arg sparse: elect to use sparse matrices, default is **False**. If
            Scipy is not found, :class:`ImportError` is raised.
        :type sparse: bool

        :arg kdtree: elect to use KDTree for building Kirchhoff matrix faster,
            default is **True**
        :type kdtree: bool


        Instances of :class:`Gamma` classes and custom functions are
        accepted as *gamma* argument.

        When Scipy is available, user can select to use sparse matrices for
        efficient usage of memory at the cost of computation speed."""

        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')

        cutoff, g, gamma = checkENMParameters(cutoff, gamma)
        self._reset()
        self._cutoff = cutoff
        self._gamma = g

        n_atoms = coords.shape[0]
        start = time.time()
        sparse = kwargs.get('sparse', False)
        if sparse:
            try:
                from scipy import sparse as scipy_sparse
            except ImportError:
                raise ImportError('failed to import scipy.sparse, which  is '
                                  'required for sparse matrix calculations')
            kirchhoff = scipy_sparse.lil_matrix((n_atoms, n_atoms))
        else:
            kirchhoff = np.zeros((n_atoms, n_atoms), 'd')

        if kwargs.get('kdtree', True):
            kdtree = KDTree(coords)
            kdtree.search(cutoff)
            dist2 = kdtree.getDistances() ** 2
            r = 0
            for i, j in kdtree.getIndices():
                g = gamma(dist2[r], i, j)
                kirchhoff[i, j] = -g
                kirchhoff[j, i] = -g
                kirchhoff[i, i] = kirchhoff[i, i] + g
                kirchhoff[j, j] = kirchhoff[j, j] + g
                r += 1
        else:
            LOGGER.info('Using slower method for building the Kirchhoff.')
            cutoff2 = cutoff * cutoff
            mul = np.multiply
            for i in range(n_atoms):
                xyz_i = coords[i, :]
                i_p1 = i+1
                i2j = coords[i_p1:, :] - xyz_i
                mul(i2j, i2j, i2j)
                for j, dist2 in enumerate(i2j.sum(1)):
                    if dist2 > cutoff2:
                        continue
                    j += i_p1
                    g = gamma(dist2, i, j)
                    kirchhoff[i, j] = -g
                    kirchhoff[j, i] = -g
                    kirchhoff[i, i] = kirchhoff[i, i] + g
                    kirchhoff[j, j] = kirchhoff[j, j] + g

        if sparse:
            kirchhoff = kirchhoff.tocsr()

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
              If **None** or ``'all'`` is given, all modes will be calculated.
        :type n_modes: int or None, default is 20

        :arg zeros: If **True**, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is **True**

        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is **True**

        """

        if self._kirchhoff is None:
            raise ValueError('Kirchhoff matrix is not built or set')
        if str(n_modes).lower() == 'all':
            n_modes = None
        assert n_modes is None or isinstance(n_modes, int) and n_modes > 0, \
            'n_modes must be a positive integer'
        assert isinstance(zeros, bool), 'zeros must be a boolean'
        assert isinstance(turbo, bool), 'turbo must be a boolean'
        self._clear()
        LOGGER.timeit('_gnm_calc_modes')
        values, vectors, vars = solveEig(self._kirchhoff, n_modes=n_modes, zeros=zeros, 
                                         turbo=turbo, expct_n_zeros=1)

        self._eigvals = values
        self._array = vectors
        self._vars = vars
        self._trace = self._vars.sum()
        self._n_modes = len(self._eigvals)
        if self._n_modes > 1:
            LOGGER.report('{0} modes were calculated in %.2fs.'
                        .format(self._n_modes), label='_gnm_calc_modes')
        else:
            LOGGER.report('{0} mode was calculated in %.2fs.'
                        .format(self._n_modes), label='_gnm_calc_modes')

    def getNormDistFluct(self, coords):
        """Normalized distance fluctuation
        """
            
        model = self.getModel()
        LOGGER.info('Number of chains: {0}, chains: {1}.'
                     .format(len(list(set(coords.getChids()))), \
                                 list(set(coords.getChids()))))

        try:
            #coords = coords.select('protein and name CA')
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                                'with `getCoords` method')
        
        if not isinstance(model, NMA):
            LOGGER.info('Calculating new model')
            model = GNM('prot analysis')
            model.buildKirchhoff(coords)
            model.calcModes() 
            
        LA = importLA()
        n_atoms = model.numAtoms()
        LOGGER.timeit('_ndf')
    
        from .analysis import calcCrossCorr
        # <dRi, dRi>, <dRj, dRj> = 1
        crossC = 2-2*calcCrossCorr(model)
        r_ij = np.zeros((n_atoms,n_atoms,3))

        for i in range(n_atoms):
            for j in range(i+1,n_atoms):
                r_ij[i][j] = coords[j,:] - coords[i,:]
                r_ij[j][i] = r_ij[i][j]
                
        r_ij_n = LA.norm(r_ij, axis=2)

        #with np.errstate(divide='ignore'):
        r_ij_n[np.diag_indices_from(r_ij_n)] = ZERO  # div by 0
        crossC = abs(crossC)
        normdistfluct = np.divide(np.sqrt(crossC), r_ij_n)
        LOGGER.report('NDF calculated in %.2lfs.', label='_ndf')
        normdistfluct[np.diag_indices_from(normdistfluct)] = 0  # div by 0
        return normdistfluct

    def setEigens(self, vectors, values=None):
        self._clear()
        super(GNM, self).setEigens(vectors, values)

class MaskedGNM(GNM, MaskedNMA):
    def __init__(self, name='Unknown', mask=False, masked=True):
        GNM.__init__(self, name)
        MaskedNMA.__init__(self, name, mask, masked)

    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        self._maskedarray = None
        super(MaskedGNM, self).calcModes(n_modes, zeros, turbo)

    def _reset(self):
        super(MaskedGNM, self)._reset()
        self._maskedarray = None

    def setKirchhoff(self, kirchhoff):
        """Set Kirchhoff matrix."""

        if not isinstance(kirchhoff, np.ndarray):
            raise TypeError('kirchhoff must be a Numpy array')
        elif (not kirchhoff.ndim == 2 or
              kirchhoff.shape[0] != kirchhoff.shape[1]):
            raise ValueError('kirchhoff must be a square matrix')
        elif kirchhoff.dtype != float:
            try:
                kirchhoff = kirchhoff.astype(float)
            except:
                raise ValueError('kirchhoff.dtype must be float')
        
        mask = self.mask
        if not self.masked:
            if not np.isscalar(mask):
                if self.is3d():
                    mask = np.repeat(mask, 3)
                try:
                    kirchhoff = kirchhoff[mask, :][:, mask]
                except IndexError:
                    raise IndexError('size mismatch between Kirchhoff (%d) and mask (%d).'
                                     'Try set masked to False or reset mask'%(len(kirchhoff), len(mask)))

        super(MaskedGNM, self).setKirchhoff(kirchhoff)

    def _getKirchhoff(self):
        """Returns the Kirchhoff matrix."""

        kirchhoff = self._kirchhoff

        if kirchhoff is None: return None

        if not self._isOriginal():
            kirchhoff = self._extend(kirchhoff, axis=None)

        return kirchhoff


def calcGNM(pdb, selstr='calpha', cutoff=10, gamma=1., n_modes=20,
            zeros=False):
    """Returns a :class:`GNM` instance and atoms used for the calculations.
    By default only alpha carbons are considered, but selection string helps
    selecting a subset of it.  *pdb* can be :class:`.Atomic` instance."""

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
        raise TypeError('pdb must be an atom container, not {0}'
                        .format(type(pdb)))
    gnm = GNM(title)
    sel = ag.select(selstr)
    gnm.buildKirchhoff(sel, cutoff, gamma)
    gnm.calcModes(n_modes, zeros)
    return gnm, sel
