# -*- coding: utf-8 -*-
"""This module defines a class and a function for anisotropic network model
(ANM) calculations."""

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB
from prody.utilities import importLA, checkCoords
from prody.kdtree import KDTree

from .nma import NMA
from .gnm import GNMBase, ZERO, checkENMParameters

__all__ = ['ANM', 'calcANM']

class ANMBase(NMA):

    def __init__(self, name='Unknown'):

        super(ANMBase, self).__init__(name)
        self._is3d = True
        self._cutoff = None
        self._gamma = None
        self._hessian = None
        self._stiffness = None 

    def _reset(self):

        GNMBase._reset(self)
        self._hessian = None
        self._stiffness = None
        self._is3d = True

    def getHessian(self):
        """Returns a copy of the Hessian matrix."""

        if self._hessian is None:
            return None
        return self._hessian.copy()

    def _getHessian(self):
        """Returns the Hessian matrix."""

        return self._hessian

    def _getStiffness(self):
        """ Return the Stiffness matrix."""

        return self._stiffness

    def getStiffness(self):
        """ Return a copy of Stiffness matrix."""

        if self._stiffness is None:
            return None
        return self._stiffness.copy()
    
    def getStiffnessRange(self):
        """ Return the range of effective spring constant."""
        if self._stiffness is None:
            return None
        else:
            return np.array(np.min(self._stiffness[np.nonzero(self._stiffness)]), \
                                               np.amax(self._stiffness))

    def getMechStiffStatistic(self, rangeK, minAA=0, AA='all'):
        """Returns number of effective spring constant with set range of
        amino acids of protein structure.
        ``AA`` can be a list with a range of analysed amino acids as:
        [first_aa, last_aa, first_aa2, last_aa2],
        minAA - eliminate amino acids that are within 20aa and
        ``rangeK`` is a list [minK, maxK]"""
        
        model = self.getModel()
        if AA == 'all':
            sm = model.getStiffness()
        elif type(AA) == int:
            sm = model.getStiffness()[0: AA, (-1)*AA-1:-1]
        elif type(AA) == list and len(AA) == 1:
            sm = model.getStiffness()[0: AA, (-1)*AA-1:-1]
        elif type(AA) == list and len(AA) == 4:
            sm = model.getStiffness()[AA[0]:AA[1],AA[2]:AA[3]]
        if minAA > 0:
            sm2 = sm[minAA:-1,0:-1-minAA]  # matrix without close contacts
            sm3 = np.tril(sm2, k=-1)
            #sort_sm2 = np.sort((np.tril(sm2, k=-1)1).flatten())
            a = np.where(np.logical_and(sm3>rangeK[0], sm3<rangeK[1]))
        if minAA == 0:
            sm2 = np.tril(sm, k=-1)
            a = np.where(np.logical_and(sm2>rangeK[0], sm2<rangeK[1]))
        return len(a[0])
        

    def getStiffnessRangeSel(self, value, minAA=20, AA='all'):
        """Returns minimum or maximum value of sping constant from 
        mechanical stiffness calculations for residues that are within 
        more than ``min_aa`` from each other. ``Value`` should be 'minK'
        or 'maxK'. It alow to avoid residues near each other. 
        ``AA`` is a number of residues from both terminus (N and C)
        of protein strcuture, it can be ``all`` or int value (than first 
        and last ``AA`` residues will be analysed. 
        With ``minAA=0`` it can be used to search the highest/lowest 
        values of interactions between N-C terminus if protein structure
        has a shear, zipper or SD1-disconnected mechanical clamp 
        -it is common in FnIII/Ig like domains and determines the maximum 
        unfolding force in AFM or SMD method."""
        
        model = self.getModel()
        if AA == 'all':
            sm = model.getStiffness()
        elif type(AA) == int:
            sm = model.getStiffness()[0: AA, (-1)*AA-1:-1]
        minK = np.min(sm[np.nonzero(sm)]) 
        maxK = np.amax(sm)
        
        if value == 'minK':
            indices = np.where(sm == np.min(sm[np.nonzero(sm)]))
        elif value == 'maxK':
            indices = np.where(sm == np.amax(sm))
        try:
            residue_diff = abs(indices[0][0]-indices[0][1])
        except: residue_diff = abs(indices[0]-indices[1])
        sm_mod = sm
        checking=True
        while checking:
            if residue_diff < minAA:
                sm_mod[indices[0][0],indices[0][1]] = 0
                sm_mod[indices[1][0],indices[1][1]] = 0
                if value == 'minK':
                    mK = np.min(sm_mod[np.nonzero(sm_mod)])
                    indices = np.where(sm_mod == np.min(sm_mod[np.nonzero(sm_mod)]))
                elif value == 'maxK':
                    mK = np.amax(sm_mod)
                    indices = np.where(sm_mod == np.amax(sm_mod))
                try:
                    residue_diff = abs(indices[0][0]-indices[0][1])
                except: residue_diff = abs(indices[0]-indices[1])
            else: 
                if value == 'minK':
                    mK = np.min(sm_mod[np.nonzero(sm_mod)])
                    indices = np.where(sm_mod == np.min(sm_mod[np.nonzero(sm_mod)]))
                elif value == 'maxK':
                    mK = np.amax(sm_mod)
                    indices = np.where(sm_mod == np.amax(sm_mod))
                checking=False
                if len(indices[0]) == 2:
                    return mK, list(indices[0])
                else:
                    return mK, list(indices[0])+list(indices[1])

    
    def setHessian(self, hessian):
        """Set Hessian matrix.  A symmetric matrix is expected, i.e. not a
        lower- or upper-triangular matrix."""

        if not isinstance(hessian, np.ndarray):
            raise TypeError('hessian must be a Numpy array')
        elif hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
            raise ValueError('hessian must be square matrix')
        elif hessian.shape[0] % 3:
            raise ValueError('hessian.shape must be (3*n_atoms,3*n_atoms)')
        elif hessian.dtype != float:
            try:
                hessian = hessian.astype(float)
            except:
                raise ValueError('hessian.dtype must be float')
        self._reset()
        self._hessian = hessian
        self._dof = hessian.shape[0]
        self._n_atoms = self._dof / 3
        
    def buildHessian(self, coords, cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å, minimum is 4.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float, :class:`Gamma`

        :arg sparse: elect to use sparse matrices, default is **False**. If
            Scipy is not found, :class:`ImportError` is raised.
        :type sparse: bool

        :arg kdtree: elect to use KDTree for building Hessian matrix,
            default is **False** since KDTree method is slower
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

        dof = n_atoms * 3
        LOGGER.timeit('_anm_hessian')

        if kwargs.get('sparse', False):
            try:
                from scipy import sparse as scipy_sparse
            except ImportError:
                raise ImportError('failed to import scipy.sparse, which  is '
                                  'required for sparse matrix calculations')
            kirchhoff = scipy_sparse.lil_matrix((n_atoms, n_atoms))
            hessian = scipy_sparse.lil_matrix((dof, dof))
        else:
            kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
            hessian = np.zeros((dof, dof), float)

        if kwargs.get('kdtree', False):
            LOGGER.info('Using KDTree for building the Hessian.')
            kdtree = KDTree(coords)
            kdtree.search(cutoff)
            for i, j in kdtree.getIndices():
                i2j = coords[j] - coords[i]
                dist2 = np.dot(i2j, i2j)
                g = gamma(dist2, i, j)
                super_element = np.outer(i2j, i2j) * (- g / dist2)
                res_i3 = i*3
                res_i33 = res_i3+3
                res_j3 = j*3
                res_j33 = res_j3+3
                hessian[res_i3:res_i33, res_j3:res_j33] = super_element
                hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                hessian[res_i3:res_i33, res_i3:res_i33] = \
                    hessian[res_i3:res_i33, res_i3:res_i33] - super_element
                hessian[res_j3:res_j33, res_j3:res_j33] = \
                    hessian[res_j3:res_j33, res_j3:res_j33] - super_element
                kirchhoff[i, j] = -g
                kirchhoff[j, i] = -g
                kirchhoff[i, i] = kirchhoff[i, i] - g
                kirchhoff[j, j] = kirchhoff[j, j] - g
        else:
            cutoff2 = cutoff * cutoff
            for i in range(n_atoms):
                res_i3 = i*3
                res_i33 = res_i3+3
                i_p1 = i+1
                i2j_all = coords[i_p1:, :] - coords[i]
                for j, dist2 in enumerate((i2j_all ** 2).sum(1)):
                    if dist2 > cutoff2:
                        continue
                    i2j = i2j_all[j]
                    j += i_p1
                    g = gamma(dist2, i, j)
                    res_j3 = j*3
                    res_j33 = res_j3+3
                    super_element = np.outer(i2j, i2j) * (- g / dist2)
                    hessian[res_i3:res_i33, res_j3:res_j33] = super_element
                    hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                    hessian[res_i3:res_i33, res_i3:res_i33] = \
                        hessian[res_i3:res_i33, res_i3:res_i33] - super_element
                    hessian[res_j3:res_j33, res_j3:res_j33] = \
                        hessian[res_j3:res_j33, res_j3:res_j33] - super_element
                    kirchhoff[i, j] = -g
                    kirchhoff[j, i] = -g
                    kirchhoff[i, i] = kirchhoff[i, i] - g
                    kirchhoff[j, j] = kirchhoff[j, j] - g
        LOGGER.report('Hessian was built in %.2fs.', label='_anm_hessian')
        self._kirchhoff = kirchhoff
        self._hessian = hessian
        self._n_atoms = n_atoms
        self._dof = dof

    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh`
        function to diagonalize the Hessian matrix. When Scipy is not found,
        :func:`numpy.linalg.eigh` is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
            If ``None`` or 'all' is given, all modes will be calculated.
        :type n_modes: int or None, default is 20

        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``

        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        """

        if self._hessian is None:
            raise ValueError('Hessian matrix is not built or set')
        if str(n_modes) is 'all':
            n_modes = None
        assert n_modes is None or isinstance(n_modes, int) and n_modes > 0, \
            'n_modes must be a positive integer'
        assert isinstance(zeros, bool), 'zeros must be a boolean'
        assert isinstance(turbo, bool), 'turbo must be a boolean'
        linalg = importLA()
        LOGGER.timeit('_anm_calc_modes')
        shift = 5
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
            if isinstance(self._hessian, np.ndarray):
                values, vectors = linalg.eigh(self._hessian, turbo=turbo,
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
                        scipy_sparse_la.eigsh(self._hessian, k=n_modes+6,
                                              which='SA'))
                except:
                    values, vectors = (
                        scipy_sparse_la.eigen_symmetric(self._hessian,
                                                        k=n_modes+6,
                                                        which='SA'))

        else:
            if n_modes is not None:
                LOGGER.info('Scipy is not found, all modes are calculated.')
            values, vectors = np.linalg.eigh(self._hessian)
        n_zeros = sum(values < ZERO)

        if n_zeros < 6:
            LOGGER.warning('Less than 6 zero eigenvalues are calculated.')
            shift = n_zeros - 1
        elif n_zeros > 6:
            LOGGER.warning('More than 6 zero eigenvalues are calculated.')
            shift = n_zeros - 1
        if zeros:
            shift = -1
        if n_zeros > n_modes:
            self._eigvals = values[1+shift:]
        else:
            self._eigvals = values[1+shift:]
        self._vars = 1 / self._eigvals
        self._trace = self._vars.sum()
        
        if shift:
            self._array = vectors[:, 1+shift:].copy()
        else:
            self._array = vectors
        self._n_modes = len(self._eigvals)
        LOGGER.report('{0} modes were calculated in %.2fs.'
                     .format(self._n_modes), label='_anm_calc_modes')

    def buildMechStiff(self, coords, n_modes=None, kbt=1.):

        """Calculate stiffness matrix calculated using :class:`.ANM` instance. 
        Method described in [EB08]_. 
    
        .. [EB08] Eyal E., Bahar I. Toward a Molecular Understanding of 
            the Anisotropic Response of Proteins to External Forces:
            Insights from Elastic Network Models. *Biophys J* **2008** 94:3424-34355. 
    
        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`.
        :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
            If ``None`` is given, all modes will be calculated (3x number of atoms).
        :type n_modes: int or ``None``, default is 20.
        
        Author: Mustafa Tekpinar & Karolina Mikulska-Ruminska & Cihan Kaya
        """

        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')
        n_atoms = natoms = self._n_atoms
        n_modes = 3 * n_atoms

        self.calcModes(n_modes=None, zeros=True)
        
        LOGGER.timeit('_sm')
        eigvecs = (np.transpose(self._array)).flatten()
        eigvals = np.transpose(self._eigvals)
        natoms = n_atoms

        sm = np.zeros((n_atoms, n_atoms), np.double)
        from .smtools import calcSM
        LOGGER.info('Calculating stiffness matrix.')

        calcSM(coords, sm, eigvecs, eigvals,
                natoms, n_modes, float(kbt))

        LOGGER.report('Stiffness matrix calculated in %.2lfs.', label='_sm')

        self._stiffness = sm
        
        LOGGER.info('The range of effective force constant is: {0} to {1}.'
                                   .format(np.min(sm[np.nonzero(sm)]), np.amax(sm)))
        

class ANM(ANMBase, GNMBase):

    """Class for Anisotropic Network Model (ANM) analysis of proteins
    ([PD00]_, [ARA01]_).

    See a usage example in :ref:`anm`.

    .. [PD00] Doruker P, Atilgan AR, Bahar I. Dynamics of proteins predicted by
       molecular dynamics simulations and analytical approaches: Application to
       a-amylase inhibitor. *Proteins* **2000** 40:512-524.

    .. [ARA01] Atilgan AR, Durrell SR, Jernigan RL, Demirel MC, Keskin O,
       Bahar I. Anisotropy of fluctuation dynamics of proteins with an
       elastic network model. *Biophys. J.* **2001** 80:505-515."""


    def __init__(self, name='Unknown'):

        super(ANM, self).__init__(name)


def calcANM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20,
            zeros=False):
    """Returns an :class:`ANM` instance and atoms used for the calculations.
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
        raise TypeError('pdb must be an atomic class, not {0}'
                        .format(type(pdb)))
    anm = ANM(title)
    sel = ag.select(selstr)
    anm.buildHessian(sel, cutoff, gamma)
    anm.calcModes(n_modes, zeros)
    return anm, sel
