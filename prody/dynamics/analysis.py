# -*- coding: utf-8 -*-
"""This module defines functions for calculating physical properties from normal
modes."""

import time

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic
from prody.ensemble import Ensemble, Conformation
from prody.trajectory import TrajBase
from prody.utilities import importLA, checkCoords, div0
from numpy import sqrt, arange, log, polyfit, array

from .nma import NMA
from .modeset import ModeSet
from .mode import VectorBase, Mode, Vector
from .gnm import GNMBase

__all__ = ['calcCollectivity', 'calcCovariance', 'calcCrossCorr',
           'calcFractVariance', 'calcSqFlucts', 'calcRMSFlucts', 'calcTempFactors',
           'calcProjection', 'calcCrossProjection',
           'calcSpecDimension', 'calcPairDeformationDist',
           'calcDistFlucts', 'calcHinges', 'calcHitTime',
           'calcAnisousFromModel', 'calcScipionScore', 'calcHemnmaScore']
           #'calcEntropyTransfer', 'calcOverallNetEntropyTransfer']

def calcCollectivity(mode, masses=None, is3d=None):
    """Returns collectivity of the mode.  This function implements collectivity
    as defined in equation 5 of [BR95]_.  If *masses* are provided, they will
    be incorporated in the calculation.  Otherwise, atoms are assumed to have
    uniform masses.

    .. [BR95] Bruschweiler R. Collective protein dynamics and nuclear
       spin relaxation. *J Chem Phys* **1995** 102:3396-3403.

    :arg mode: mode(s) or vector(s)
    :type mode: :class:`.Mode`, :class:`.Vector`, :class:`.ModeSet`,
        :class:`.NMA`, :class:`~numpy.ndarray`

    :arg masses: atomic masses
    :type masses: :class:`numpy.ndarray`
    
    :arg is3d: whether mode is 3d. Default is **None** which means determine 
        the value based on ``mode.is3d()``.
    :type is3d: bool
    """

    if isinstance(mode, np.ndarray):
        V = mode
        ndim = V.ndim
        shape = V.shape

        if is3d is None:
            is3d = False

        if ndim == 0:
            raise ValueError('mode cannot be an empty array')
        elif ndim == 1:
            V = V[:, np.newaxis]

        n = shape[0]
        if is3d:
            n_atoms = n // 3
        else:
            n_atoms = n
    else:
        V, W, is3d_, n_atoms = _getModeProperties(mode)
        if is3d is None:
            is3d = is3d_
    
    colls = []

    def log0(a):
        return log(a + np.finfo(float).eps)

    for v in V.T:
        if is3d:
            u2in = (v ** 2)
            u2in_Nx3 = np.reshape(u2in, (n_atoms, 3))
            u2in = u2in_Nx3.sum(axis=1)
        else:
            u2in = (v ** 2)
        if masses is not None:
            if len(masses) != n_atoms:
                raise ValueError('length of masses must be equal to number of atoms')
            u2in = u2in / masses
        u2in = u2in * (1 / u2in.sum() ** 0.5)
        coll = np.exp(-(u2in * log0(u2in)).sum()) / n_atoms
        colls.append(coll)
    
    if len(colls) == 1:
        return coll
    else:
        return np.array(colls)

def calcSpecDimension(mode):

    """
    :arg mode: mode or vector
    :type mode: :class:`.Mode` or :class:`.Vector`

    """
    # if not isinstance(mode, Mode):
    #     raise TypeError('mode must be a Mode instance')
    
    length = mode.shape[0]
    numbers = arange(2,length+1)
    ds,p=polyfit(log(sqrt(mode[0:int(length*0.25)])),log(numbers[0:int(length*0.25)]),1)
    
    return ds

def calcFracDimension(mode):
    """
    :arg mode: mode or vector
    :type mode: mode or vector """




def calcFractVariance(mode):
    """Returns fraction of variance explained by the *mode*.  Fraction of
    variance is the ratio of the variance along a mode to the trace of the
    covariance matrix of the model."""

    if isinstance(mode, Mode):
        var = mode.getVariance()
        trace = mode.getModel()._getTrace()
    elif isinstance(mode, (ModeSet, NMA)):
        var = mode.getVariances()
        if isinstance(mode, ModeSet):
            trace = mode.getModel()._getTrace()
        else:
            trace = mode._getTrace()
    else:
        raise TypeError('mode must be a Mode instance')
    if trace is None:
        raise ValueError('modes are not calculated')

    return var / trace


def calcProjection(ensemble, modes, rmsd=True, norm=False):
    """Returns projection of conformational deviations onto given modes.
    *ensemble* coordinates are used to calculate the deviations that are
    projected onto *modes*.  For K conformations and M modes, a (K,M)
    matrix is returned.

    :arg ensemble: an ensemble, trajectory or a conformation for which
        deviation(s) will be projected, or a deformation vector
    :type ensemble: :class:`.Ensemble`, :class:`.Conformation`,
        :class:`.Vector`, :class:`.Trajectory`
        
    :arg modes: up to three normal modes
    :type modes: :class:`.Mode`, :class:`.ModeSet`, :class:`.NMA`

    By default, root-mean-square deviation (RMSD) along the normal mode is
    calculated. To calculate the raw projection pass ``rmsd=False``.

    By default, the projection is not normalized. If you would like it to be,
    pass ``norm=True``.

    :class:`.Vector` instances are accepted as *ensemble* argument to allow
    for projecting a deformation vector onto normal modes."""

    if not isinstance(ensemble, (Ensemble, Conformation, Vector, TrajBase)):
        raise TypeError('ensemble must be Ensemble, Conformation, Vector, '
                        'or a TrajBase, not {0}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, VectorBase)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0}'
                        .format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be 3-dimensional')
    if isinstance(ensemble, Vector):
        n_atoms = ensemble.numAtoms()
    else:
        n_atoms = ensemble.numSelected()
    if n_atoms != modes.numAtoms():
        raise ValueError('number of atoms are not the same')
    if isinstance(ensemble, Vector):
        if not ensemble.is3d():
            raise ValueError('ensemble must be a 3d vector instance')
        deviations = ensemble._getArray()
    elif isinstance(ensemble, (Ensemble, Conformation)):
        deviations = ensemble.getDeviations()
    else:
        nfi = ensemble.nextIndex()
        ensemble.goto(0)
        deviations = np.array([frame.getDeviations() for frame in ensemble])
        ensemble.goto(nfi)
    if deviations.ndim == 3:
        deviations = deviations.reshape((deviations.shape[0],
                                         deviations.shape[1] * 3))
    elif deviations.ndim == 2:
        deviations = deviations.reshape((1, deviations.shape[0] * 3))
    else:
        deviations = deviations.reshape((1, deviations.shape[0]))
    la = importLA()
    if norm:
        N = la.norm(deviations)
        if N != 0:
            deviations = deviations / N
    projection = np.dot(deviations, modes._getArray())
    if rmsd:
        projection = (1 / (n_atoms ** 0.5)) * projection
    return projection


def calcCrossProjection(ensemble, mode1, mode2, scale=None, **kwargs):
    """Returns projection of conformational deviations onto modes from
    different models.

    :arg ensemble: ensemble for which deviations will be projected
    :type ensemble: :class:`.Ensemble`

    :arg mode1: normal mode to project conformations onto
    :type mode1: :class:`.Mode`, :class:`.Vector`

    :arg mode2: normal mode to project conformations onto
    :type mode2: :class:`.Mode`, :class:`.Vector`

    :arg scale: scale width of the projection onto mode1 (``x``) or mode2(``y``),
        an optimized scaling factor (scalar) will be calculated by default 
        or a value of scalar can be passed.
        
    This function uses calcProjection and its arguments can be 
    passed to it as keyword arguments.
    By default, this function applies RMSD scaling and normalisation. 
    These can be turned off with ``rmsd=False`` and ``norm=False``."""

    if not isinstance(ensemble, (Ensemble, Conformation, Vector, TrajBase)):
        raise TypeError('ensemble must be Ensemble, Conformation, Vector, '
                        'or a Trajectory, not {0}'.format(type(ensemble)))
    if not isinstance(mode1, VectorBase):
        raise TypeError('mode1 must be a Mode instance, not {0}'
                        .format(type(mode1)))
    if not mode1.is3d():
        raise ValueError('mode1 must be 3-dimensional')
    if not isinstance(mode2, VectorBase):
        raise TypeError('mode2 must be a Mode instance, not {0}'
                        .format(type(mode2)))
    if not mode2.is3d():
        raise ValueError('mode2 must be 3-dimensional')

    if scale is not None:
        assert isinstance(scale, str), 'scale must be a string'
        scale = scale.lower()
        assert scale in ('x', 'y'), 'scale must be x or y'

    xcoords = calcProjection(ensemble, mode1, kwargs.get('rmsd', True), kwargs.get('norm', True))
    ycoords = calcProjection(ensemble, mode2, kwargs.pop('rmsd', True), kwargs.pop('norm', True))
    if scale:
        scalar = kwargs.get('scalar', None)
        if scalar:
            assert isinstance(scalar, (float, int)), 'scalar must be a number'
        else:
            scalar = ((ycoords.max() - ycoords.min()) /
                      (xcoords.max() - xcoords.min())
                      ) * np.sign(np.dot(xcoords, ycoords))
            if scale == 'x':
                LOGGER.info('Projection onto {0} is scaled by {1:.2f}'
                            .format(mode1, scalar))
            else:
                scalar = 1 / scalar
                LOGGER.info('Projection onto {0} is scaled by {1:.2f}'
                            .format(mode2, scalar))

        if scale == 'x':
            xcoords = xcoords * scalar
        else:
            ycoords = ycoords * scalar

    return xcoords, ycoords

def _getModeProperties(modes):
    V = []; W = []; is3d = None; n_atoms = 0
    if isinstance(modes, VectorBase):
        V = modes._getArray()
        if isinstance(modes, Mode):
            W = modes.getVariance()
        else:
            W = 1.
        V = np.asarray([V]).T
        W = np.asarray([[W]])
        is3d = modes.is3d()
        n_atoms = modes.numAtoms()
    elif isinstance(modes, (NMA, ModeSet)):
        V = modes._getArray()
        W = np.diag(modes.getVariances())
        is3d = modes.is3d()
        n_atoms = modes.numAtoms()
    elif isinstance(modes, list):
        for mode in modes:
            if not isinstance(mode, VectorBase):
                raise TypeError('modes can be a list of VectorBase instances, '
                                'not {0}'.format(type(mode)))
            V.append(mode._getArray())
            if isinstance(mode, Mode):
                W.append(modes.getVariance())
            else:
                W.append(1.)
            if is3d is None:
                is3d = mode.is3d()
                n_atoms = mode.numAtoms()
            else:
                if is3d != mode.is3d():
                    raise ValueError('modes must be either all from ANM or GNM')
                if n_atoms != mode.numAtoms():
                    raise ValueError('each mode in the list must have the same number of atoms')
        V = np.array(V).T
        W = np.diag(W)
    else:
        raise TypeError('modes must be a Mode, NMA, ModeSet instance, '
                        'or a list of Mode instances, not {0}'.format(type(modes)))
    return V, W, is3d, n_atoms

def calcSqFlucts(modes):
    """Returns sum of square-fluctuations for given set of normal *modes*.
    Square fluctuations for a single mode is obtained by multiplying the
    square of the mode array with the variance (:meth:`.Mode.getVariance`)
    along the mode.  For :class:`.PCA` and :class:`.EDA` models built using
    coordinate data in Å, unit of square-fluctuations is |A2|, for
    :class:`.ANM` and :class:`.GNM`, on the other hand, it is arbitrary or
    relative units."""

    V = []; W = []; is3d = None; n_atoms = 0
    V, W, is3d, n_atoms = _getModeProperties(modes)

    sq_flucts = np.dot(V * V, W).sum(axis=1)

    if is3d:
        sq_flucts_Nx3 = np.reshape(sq_flucts, (n_atoms, 3))
        sq_flucts = sq_flucts_Nx3.sum(axis=1)
    return sq_flucts

def calcRMSFlucts(modes):
    """Returns root mean square fluctuation(s) (RMSF) for given set of normal *modes*.
    This is calculated just by doing the square root of the square fluctuations """
    sq_flucts = calcSqFlucts(modes)
    if len(np.where(sq_flucts<0)[0]) != 0:
        raise ValueError("Square Fluctuation should not contain negative values, please check input modes")

    return sq_flucts ** 0.5

def calcCrossCorr(modes, n_cpu=1, norm=True):
    """Returns cross-correlations matrix.  For a 3-d model, cross-correlations
    matrix is an NxN matrix, where N is the number of atoms.  Each element of
    this matrix is the trace of the submatrix corresponding to a pair of atoms.
    Cross-correlations matrix may be calculated using all modes or a subset of modes
    of an NMA instance.  For large systems, calculation of cross-correlations
    matrix may be time consuming.  Optionally, multiple processors may be
    employed to perform calculations by passing ``n_cpu=2`` or more."""

    if not isinstance(n_cpu, int):
        raise TypeError('n_cpu must be an integer')
    elif n_cpu < 1:
        raise ValueError('n_cpu must be equal to or greater than 1')

    if not isinstance(modes, (Mode, Vector, NMA, ModeSet)):
        if isinstance(modes, list):
            try:
                is3d = modes[0].is3d()
            except:
                raise TypeError('modes must be a list of Mode or Vector instances, '
                            'not {0}'.format(type(modes)))
        else:
            raise TypeError('modes must be a Mode, Vector, NMA, or ModeSet instance, '
                            'not {0}'.format(type(modes)))
    else:
        is3d = modes.is3d()

    if is3d:
        model = modes
        if isinstance(modes, (Mode, ModeSet)):
            model = modes._model
            if isinstance(modes, (Mode)):
                indices = [modes.getIndex()]
                n_modes = 1
            else:
                indices = modes.getIndices()
                n_modes = len(modes)
        elif isinstance(modes, Vector):
                indices = [0]
                n_modes = 1
        else:
            n_modes = len(modes)
            indices = np.arange(n_modes)
            
        array = model._getArray()
        n_atoms = model._n_atoms

        if not isinstance(modes, Vector):
            variances = model._vars
        else:
            array = array.reshape(-1, 1)
            variances = np.ones(1)

        if n_cpu == 1:
            s = (n_modes, n_atoms, 3)
            arvar = (array[:, indices]*variances[indices]).T.reshape(s)
            array = array[:, indices].T.reshape(s)
            covariance = np.tensordot(array.transpose(2, 0, 1),
                                      arvar.transpose(0, 2, 1),
                                      axes=([0, 1], [1, 0]))
        else:
            import multiprocessing
            n_cpu = min(multiprocessing.cpu_count(), n_cpu)
            queue = multiprocessing.Queue()
            size = n_modes / n_cpu
            for i in range(n_cpu):
                if n_cpu - i == 1:
                    indices = modes.indices[i*size:]
                else:
                    indices = modes.indices[i*size:(i+1)*size]
                process = multiprocessing.Process(
                    target=_crossCorrelations,
                    args=(queue, n_atoms, array, variances, indices))
                process.start()
            while queue.qsize() < n_cpu:
                time.sleep(0.05)
            covariance = queue.get()
            while queue.qsize() > 0:
                covariance += queue.get()
    else:
        covariance = calcCovariance(modes)
    if norm:
        diag = np.power(covariance.diagonal(), 0.5)
        D = np.outer(diag, diag)
        covariance = div0(covariance, D)
    return covariance


def _crossCorrelations(queue, n_atoms, array, variances, indices):
    """Calculate covariance-matrix for a subset of modes."""

    n_modes = len(indices)
    arvar = (array[:, indices] * variances[indices]).T.reshape((n_modes,
                                                                n_atoms, 3))
    array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
    covariance = np.tensordot(array.transpose(2, 0, 1),
                              arvar.transpose(0, 2, 1),
                              axes=([0, 1], [1, 0]))
    queue.put(covariance)

def calcDistFlucts(modes, n_cpu=1, norm=True):
    """Returns the matrix of distance fluctuations (i.e. an NxN matrix
    where N is the number of residues, of MSFs in the inter-residue distances)
    computed from the cross-correlation matrix (see Eq. 12.E.1 in [IB18]_). 
    The arguments are the same as in :meth:`.calcCrossCorr`.

    .. [IB18] Dill K, Jernigan RL, Bahar I. Protein Actions: Principles and
       Modeling. *Garland Science* **2017**. """

    cc = calcCrossCorr(modes, n_cpu=n_cpu, norm=norm)
    cc_diag = np.diag(cc).reshape(-1,1)
    distFluct = cc_diag.T + cc_diag -2.*cc
    return distFluct

def calcTempFactors(modes, atoms):
    """Returns temperature (β) factors calculated using *modes* from a
    :class:`.ANM` or :class:`.GNM` instance scaled according to the 
    experimental B-factors from *atoms*."""

    model = modes.getModel()
    if not isinstance(model, GNMBase):
        raise TypeError('modes must come from GNM or ANM')
    if model.numAtoms() != atoms.numAtoms():
        raise ValueError('modes and atoms must have same number of nodes')
    sqf = calcSqFlucts(modes)
    expBetas = atoms.getBetas()
    # add warning message if experimental B-factors are zeros or meaningless (e.g., having same values)?
    if expBetas.max() < 0.5 or expBetas.std() < 0.5:
        LOGGER.warning('Experimental B-factors are quite small or meaningless. The calculated B-factors may be incorrect.')
    return sqf * (expBetas.sum() / sqf.sum())


def calcCovariance(modes):
    """Returns covariance matrix calculated for given *modes*.
    This is 3Nx3N for 3-d models and NxN (equivalent to cross-correlations) 
    for 1-d models such as GNM."""

    if isinstance(modes, NMA):
        return modes.getCovariance()
    else:
        V, W, _, _ = _getModeProperties(modes)
        return np.dot(V, np.dot(W, V.T))


def calcPairDeformationDist(model, coords, ind1, ind2, kbt=1.):                                       
    """Returns distribution of the deformations in the distance contributed by each mode 
    for selected pair of residues *ind1* *ind2* using *model* from a :class:`.ANM`.
    Method described in [EB08]_ equation (10) and figure (2).     
    
    .. [EB08] Eyal E., Bahar I. Toward a Molecular Understanding of 
        the Anisotropic Response of Proteins to External Forces:
        Insights from Elastic Network Models. *Biophys J* **2008** 94:3424-34355. 
    
    :arg model: this is an 3-dimensional :class:`NMA` instance from a :class:`.ANM`
        calculations.
    :type model: :class:`.ANM`  

    :arg coords: a coordinate set or an object with :meth:`getCoords` method.
        Recommended: ``coords = parsePDB('pdbfile').select('protein and name CA')``.
    :type coords: :class:`~numpy.ndarray`.

    :arg ind1: first residue number.
    :type ind1: int 
    
    :arg ind2: second residue number.
    :type ind2: int 
    """

    try:
        resnum_list = coords.getResnums()
        resnam_list = coords.getResnames()
        coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                coords.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be a Numpy array or an object '
                            'with `getCoords` method')
    
    if not isinstance(model, NMA):
        raise TypeError('model must be a NMA instance')
    elif not model.is3d():
        raise TypeError('model must be a 3-dimensional NMA instance')
    elif len(model) == 0:
        raise ValueError('model must have normal modes calculated')
    
    linalg = importLA()
    n_atoms = model.numAtoms()
    n_modes = model.numModes()
    LOGGER.timeit('_pairdef')

    r_ij = np.zeros((n_atoms,n_atoms,3))
    r_ij_norm = np.zeros((n_atoms,n_atoms,3))

    for i in range(n_atoms):
        for j in range(i+1,n_atoms):
            r_ij[i][j] = coords[j,:] - coords[i,:]
            r_ij[j][i] = r_ij[i][j]
            r_ij_norm[i][j] = r_ij[i][j]/linalg.norm(r_ij[i][j])
            r_ij_norm[j][i] = r_ij_norm[i][j]

    eigvecs = model.getEigvecs()
    eigvals = model.getEigvals()
    
    D_pair_k = []
    mode_nr = []
    ind1 = ind1 - resnum_list[0]
    ind2 = ind2 - resnum_list[0]

    for m in range(6,n_modes):
        U_ij_k = [(eigvecs[ind1*3][m] - eigvecs[ind2*3][m]), (eigvecs[ind1*3+1][m] \
            - eigvecs[ind2*3+1][m]), (eigvecs[ind1*3+2][m] - eigvecs[ind2*3+2][m])] 
        D_ij_k = abs(sqrt(kbt/eigvals[m])*(np.vdot(r_ij_norm[ind1][ind2], U_ij_k)))  
        D_pair_k.append(D_ij_k)
        mode_nr.append(m)

    LOGGER.report('Deformation was calculated in %.2lfs.', label='_pairdef')
    
    return mode_nr, D_pair_k

def calcHinges(modes, atoms=None, flag=False):
    """Returns the hinge sites identified using normal modes.     

    :arg modes: normal modes of which will be used to identify hinge sites
    :type modes: :class:`.GNM`  
    
    :arg atoms: an Atomic object on which to map hinges. The output will then be a selection. 
    :type atoms: :class:`.Atomic`

    :arg flag: whether return flag or index array. Default is **False**
    :type flag: bool

    """

    def identify(v):
        # obtain the signs of eigenvector
        s = np.sign(v)
        # obtain the relative magnitude of eigenvector
        mag = np.sign(np.diff(np.abs(v)))
        # obtain the cross-overs
        torf = np.diff(s)!=0
        torf = np.append(torf, [False], axis=0)
        # find which side is more close to zero
        for j, m in enumerate(mag):
            if torf[j] and m < 0:
                torf[j+1] = True
                torf[j] = False

        return torf

    if modes.is3d():
        raise ValueError('3D models are not supported.')

    # obtain the eigenvectors
    V = modes.getArray()
    if V.ndim == 1:
        hinges = identify(V)
    elif V.ndim == 2:
        _, n = V.shape
        hinges = []
        for i in range(n):
            v = V[:, i]
            torf = identify(v)
            hinges.append(torf)

        hinges = np.stack(hinges).T
    else:
        raise TypeError('wrong dimension of the array: %d'%V.ndim)
    
    if not flag:
        hinge_list = np.where(hinges)[0]
        if atoms is not None:
            if isinstance(atoms, Atomic):
                return atoms[hinge_list]
            else:
                raise TypeError('atoms should be an Atomic object')
        return sorted(set(hinge_list))
    return hinges

def calcHitTime(model, method='standard'):
    """Returns the hit and commute times between pairs of nodes calculated 
    based on a :class:`.NMA` object. 

    .. [CB95] Chennubhotla C., Bahar I. Signal Propagation in Proteins and Relation
    to Equilibrium Fluctuations. *PLoS Comput Biol* **2007** 3(9).

    :arg model: model to be used to calculate hit times
    :type model: :class:`.NMA`  

    :arg method: method to be used to calculate hit times. Available options are 
        ``"standard"`` or ``"kirchhoff"``. Default is ``"standard"``
    :type method: str

    :returns: (:class:`~numpy.ndarray`, :class:`~numpy.ndarray`)
    """

    try:
        K = model.getKirchhoff()
    except AttributeError:
        raise TypeError('model must be an NMA instance')

    if K is None:
        raise ValueError('model not built')
    
    method = method.lower()

    D = np.diag(K)
    A = np.diag(D) - K

    start = time.time()
    linalg = importLA()
    if method == 'standard':
        st = D / sum(D)

        P = np.dot(np.diag(D**(-1)), A)
        W = np.ones((len(st), 1)) * st.T
        Z = linalg.pinv(np.eye(P.shape[0], P.shape[1]) - P + W)

        H = np.ones((len(st), 1)) * np.diag(Z).T - Z
        H = H / W
        H = H.T

    elif method == 'kirchhoff':
        K_inv = linalg.pinv(K)
        sum_D = sum(D)

        T1 = (sum_D * np.ones((len(D),1)) * np.diag(K_inv)).T

        T2 = sum_D * K_inv
        T3_i = np.dot((np.ones((len(D),1)) * D), K_inv)

        H = T1 - T2 + T3_i - T3_i.T

    C = H + H.T

    LOGGER.debug('Hit and commute times are calculated in  {0:.2f}s.'
                 .format(time.time()-start)) 
    return H, C


def calcAnisousFromModel(model, ):
    """Returns a Nx6 matrix containing anisotropic B factors (ANISOU lines)
    from a covariance matrix calculated from **model**.

    :arg model: 3D model from which to calculate covariance matrix
    :type model: :class:`.ANM`, :class:`.PCA`

    .. ipython:: python

       from prody import *
       protein = parsePDB('1ejg')
       anm, calphas = calcANM(protein)
       adp_matrix = calcAnisousFromModel(anm)"""

    if not isinstance(model, (NMA, Mode)) or not model.is3d():
        raise TypeError('model must be of type ANM, PCA or Mode, not {0}'
                        .format(type(model)))

    cov = calcCovariance(model)
    n_atoms = model.numAtoms()
    
    submatrices = [cov[i*3:(i+1)*3, i*3:(i+1)*3] for i in range(n_atoms)]

    anisou = np.zeros((n_atoms, 6))
    for index, submatrix in enumerate(submatrices):
        anisou[index, 0] = submatrix[0, 0]
        anisou[index, 1] = submatrix[1, 1]
        anisou[index, 2] = submatrix[2, 2]
        anisou[index, 3] = submatrix[0, 1]
        anisou[index, 4] = submatrix[0, 2]
        anisou[index, 5] = submatrix[1, 2]
    return anisou


def calcScipionScore(modes):
    """Calculate the score from hybrid electron microscopy normal mode analysis (HEMNMA) 
    [CS14]_ as implemented in the Scipion continuousflex plugin [MH20]_. This score 
    prioritises modes as a function of mode number and collectivity order.

    .. [CS14] Sorzano COS, de la Rosa-Trevín JM, Tama F, Jonić S.
       Hybrid Electron Microscopy Normal Mode Analysis graphical interface and protocol.
       *J Struct Biol* **2014** 188:134-41.

    .. [MH20] Harastani M, Sorzano COS, Jonić S. 
       Hybrid Electron Microscopy Normal Mode Analysis with Scipion.
       *Protein Sci* **2020** 29:223-236.

    :arg modes: mode(s) or vector(s)
    :type modes: :class:`.Mode`, :class:`.Vector`, :class:`.ModeSet`, :class:`.NMA`
    """
    n_modes = modes.numModes()
    
    if n_modes > 1:
        collectivityList = list(calcCollectivity(modes))
    else:
        collectivityList = [calcCollectivity(modes)]

    idxSorted = [i[0] for i in sorted(enumerate(collectivityList),
                                      key=lambda x: x[1],
                                      reverse=True)]

    score = np.zeros(n_modes)
    modeNum = list(range(n_modes))

    for i in range(n_modes):
        score[idxSorted[i]] = idxSorted[i] + modeNum[i] + 2  

    score = score / (2.0 * n_modes) 

    return score

calcHemnmaScore = calcScipionScore
