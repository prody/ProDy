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
from prody.utilities import importLA, checkCoords

from .nma import NMA
from .gamma import Gamma

__all__ = ['GNM', 'calcGNM']

ZERO = 1e-6


class GNMBase(NMA):

    """Class for Gaussian Network Model analysis of proteins."""

    def __init__(self, name='Unknown'):

        super(GNMBase, self).__init__(name)
        self._is3d = False
        self._cutoff = None
        self._kirchhoff = None
        self._gamma = None
        self._hinges = None

    def __repr__(self):

        return ('<{0}: {1} ({2} modes; {3} nodes)>'
                .format(self.__class__.__name__, self._title, self.__len__(),
                        self._n_atoms))

    def __str__(self):

        return self.__class__.__name__ + ' ' + self._title

    def _reset(self):

        NMA._reset(self)
        self._cutoff = None
        self._gamma = None
        self._kirchhoff = None
        self._is3d = False

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
        return self._kirchhoff.copy()

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
        if kwargs.get('sparse', False):
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

        LOGGER.debug('Kirchhoff was built in {0:.2f}s.'
                     .format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._n_atoms = n_atoms
        self._dof = n_atoms

    def calcModes(self, n_modes=20, zeros=False, turbo=True, hinges=True):
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

        :arg hinges: Identify hinge sites after modes are computed.
        :type hinges: bool, default is ``True``
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
                    values, vectors = (
                        scipy_sparse_la.eigsh(self._kirchhoff,
                                              k=n_modes + 1, which='SA'))
                except:
                    values, vectors = (
                        scipy_sparse_la.eigen_symmetric(self._kirchhoff,
                                                        k=n_modes + 1,
                                                        which='SA'))
        else:
            if n_modes is not None:
                LOGGER.info('Scipy is not found, all modes are calculated.')
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
        if hinges:
            self.calcHinges()
        LOGGER.debug('{0} modes were calculated in {1:.2f}s.'
                     .format(self._n_modes, time.time()-start))

    def calcHinges(self):
        if self._array is None:
            raise ValueError('Modes are not calculated.')
        # obtain the eigenvectors
        V = self._array
        (m, n) = V.shape
        hinges = []
        for i in range(n):
            v = V[:,i]
            # obtain the signs of eigenvector
            s = np.insert(np.sign(v), 0, 0)
            # obtain the relative magnitude of eigenvector
            mag = np.insert(np.sign(np.diff(np.abs(v))), 0, 0)
            # obtain the cross-overs
            torf = np.diff(s)!=0
            indices = np.where(torf)[0]
            # find which side is more close to zero
            for i in range(len(indices)):
                idx = indices[i]
                if mag[idx] > 0:
                    indices[i] -= 1
            hinges.append(indices)
        self._hinges = np.array(hinges)
        return self._hinges

    def getHinges(self, modeIndex=None):
        """Get residue index of hinge sites given mode indices.

        :arg modeIndex: indices of modes. This parameter can be a scalar, a list, 
            or logical indices.
        :type modeIndex: int or list, default is ``None``
        """
        if self._hinges is None:
            raise ValueError('Hinges are not calculated.')
        if modeIndex is None:
            hinges = self._hinges
        else:
            hinges = self._hinges[modeIndex]
        if hinges.dtype is np.dtype('O'):
            hingelist = [j for i in hinges for j in i]
        else:
            hingelist = [i for i in hinges]
        return list(set(hingelist))

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
            
        linalg = importLA()
        n_atoms = model.numAtoms()
        n_modes = model.numModes()
        LOGGER.timeit('_ndf')
    
        from .analysis import calcCrossCorr
        from numpy import linalg as LA
        # <dRi, dRi>, <dRj, dRj> = 1
        crossC = 2-2*calcCrossCorr(model)
        r_ij = np.zeros((n_atoms,n_atoms,3))

        for i in range(n_atoms):
           for j in range(i+1,n_atoms):
               r_ij[i][j] = coords[j,:] - coords[i,:]
               r_ij[j][i] = r_ij[i][j]
               r_ij_n = LA.norm(r_ij, axis=2)

        #with np.errstate(divide='ignore'):
        r_ij_n[np.diag_indices_from(r_ij_n)] = 1e-5  # div by 0
        crossC=abs(crossC)
        normdistfluct = np.divide(np.sqrt(crossC),r_ij_n)
        LOGGER.report('NDF calculated in %.2lfs.', label='_ndf')
        normdistfluct[np.diag_indices_from(normdistfluct)] = 0  # div by 0
        return normdistfluct


def calcGNM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20,
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
