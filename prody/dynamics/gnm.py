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

__all__ = ['GNM', 'calcGNM', 'TrimedGNM', 'test']

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

    def buildAffinity(self, coords=None, cutoff=7.0):

        if self._kirchhoff is None:
            try: 
                self.buildKirchhoff(coords, cutoff)
            except:
                raise TypeError('You should provide the coordinates of the proteins')

        if not isinstance(self._kirchhoff, np.ndarray):
            raise TypeError('kirchhoff must be a Numpy array')
        elif (not self._kirchhoff.ndim == 2 or
              self._kirchhoff.shape[0] != self._kirchhoff.shape[1]):
            raise ValueError('kirchhoff must be a square matrix')
        elif self._kirchhoff.dtype != float:
            try:
                self.kirchhoff = self.kirchhoff.astype(float)
            except:
                raise ValueError('kirchhoff.dtype must be float')

        from scipy import sparse

        self._diagonal = np.diag(self._kirchhoff)
        
        self._affinity = sparse.spdiags(self._diagonal, 0, len(self._diagonal), 
            len(self._diagonal)).toarray() - self._kirchhoff

    def calcHitTime(self, method='Z'):

        if self._affinity is None:
            raise ValueError('Affinity matrix is not built or set')

        start = time.time()
        linalg = importLA()
        if method == 'Z':

            D = self._diagonal
            A = self._affinity

            st = D / sum(D)

            P = np.dot(np.diag(D**(-1)), A)

            W = np.ones((len(st),1)) * st.T

            Z = linalg.pinv(np.eye(P.shape[0], P.shape[1]) - P + W)

            H = np.ones((len(st),1)) * np.diag(Z).T - Z
            H = H / W
            H = H.T

        elif method == 'K':

            K = self._kirchhoff
            D = self._diagonal

            K_inv = linalg.pinv(K)
            sum_D = sum(D)

            T1 = (sum_D * np.ones((len(D),1)) * diag(K_inv)).T

            T2 = sum_D * K_inv
            T3_i = np.dot((np.ones((len(D),1)) * D), K_inv)

            H = T1 - T2 + T3_i - T3_i.T

        self._hitTime = H
        self._commuteTime = H + H.T


        LOGGER.debug('Hitting and commute time are calculated in  {0:.2f}s.'
                     .format(time.time()-start))    

    def getAffinity(self):
        """Returns a copy of the Kirchhoff matrix."""

        if self._affinity is None:
            return None
        return self._affinity.copy()

    def _getAffinity(self):
        """Returns the Kirchhoff matrix."""

        return self._affinity

    def getDiagonal(self):
        """Returns a copy of the Kirchhoff matrix."""

        if self._diagonal is None:
            return None
        return self._diagonal.copy()

    def _getDiagonal(self):
        """Returns the Kirchhoff matrix."""

        return self._diagonal

    def getHitTime(self):
        """Returns a copy of the hit time matrix."""

        if self._hitTime is None:
            return None
        return self._hitTime.copy()

    def _getHitTime(self):
        """Returns the hit time matrix."""

        return self._getHitTime

    def getCommuteTime(self):
        """Returns a copy of the Kirchhoff matrix."""

        if self._commuteTime is None:
            return None
        return self._commuteTime.copy()

    def _getCommuteTime(self):
        """Returns the Kirchhoff matrix."""

        return self._commuteTime    


    def calcModes(self, n_modes=20, zeros=False, turbo=True, hinges=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh`
        function to diagonalize the Kirchhoff matrix. When Scipy is not found,
        :func:`numpy.linalg.eigh` is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
              If ``None`` or 'all' is given, all modes will be calculated.
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
        if str(n_modes).lower() is 'all':
            n_modes = None
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
            s = np.sign(v)
            # obtain the relative magnitude of eigenvector
            mag = np.sign(np.diff(np.abs(v)))
            # obtain the cross-overs
            torf = np.diff(s)!=0
            indices = np.where(torf)[0]
            # find which side is more close to zero
            for i in range(len(indices)):
                idx = indices[i]
                if mag[idx] < 0:
                    indices[i] += 1
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
            LOGGER.info('Warning: hinges are not calculated, thus null is returned. '
                        'Please call GNM.calcHinges() to calculate the hinge sites first.')
            return None
        if modeIndex is None:
            hinges = self._hinges
        else:
            hinges = self._hinges[modeIndex]
        if hinges.dtype is np.dtype('O'):
            hingelist = [j for i in hinges for j in i]
        else:
            hingelist = [i for i in hinges]
        return sorted(set(hingelist))
    
    def numHinges(self, modeIndex=None):
        return len(self.getHinges(modeIndex=modeIndex))

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
            zeros=False, hinges=True):
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
    gnm.calcModes(n_modes, zeros, hinges=hinges)
    return gnm, sel

class TrimedGNM(GNM):
    def __init__(self, name='Unknown', mask=False, useTrimed=True):
        super(TrimedGNM, self).__init__(name)
        self.mask = False
        self.useTrimed = useTrimed

        if not np.isscalar(mask):
            self.mask = np.array(mask)

    def numAtoms(self):
        """Returns number of atoms."""

        if self.useTrimed or np.isscalar(self.mask):
            return self._n_atoms
        else:
            return len(self.mask)

    def getArray(self):
        """Returns a copy of eigenvectors array."""

        if self._array is None: return None

        array = self._array.copy()

        if self.useTrimed or np.isscalar(self.mask):
            return array

        mask = self.mask.copy()
        N = len(mask)
        n, m = array.shape
        whole_array = np.zeros((N,m))
        mask = np.expand_dims(mask, axis=1)
        mask = mask.repeat(m, axis=1)
        whole_array[mask] = array.flatten()
        return whole_array

    getEigvecs = getArray

    def _getArray(self):
        """Returns eigenvectors array. The function returns 
        a copy of the array if useTrimed is ``True``."""

        if self._array is None: return None

        if self.useTrimed or np.isscalar(self.mask):
            return self._array
        else:
            return self.getArray()

    def fixTail(self, length):
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

        if np.isscalar(self.mask):
            self.mask = ones(self.numAtoms(), dtype=bool)

        self.mask = _fixLength(self.mask, length, False)
        return

def test():
    
    from prody import parsePDB
    pdb = parsePDB('1z83',subset='ca',chain='A')

    gnm = GNM()
    gnm.buildAffinity(pdb)
    gnm.calcHitTime()

    hitTime = gnm.getHitTime()
    commuteTime = gnm.getCommuteTime()

    return hitTime, commuteTime