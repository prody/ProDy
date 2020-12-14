# -*- coding: utf-8 -*-
"""This module defines a class and a function for explicit membrane GNM calculations."""

import numpy as np

from prody import LOGGER, PY3K
from prody.atomic import Atomic, AtomGroup
from prody.utilities import importLA, checkCoords, copy
from numpy import sqrt, zeros, array, ceil, dot

from .gnm import GNM, checkENMParameters
from .editing import _reduceModel

LA = importLA()
inv = LA.inv
pinv = LA.pinv
norm = LA.norm

__all__ = ['exGNM']

class exGNM(GNM):

    """Class for explicit GNM (exGNM) method ([FT00]_).
    Optional arguments build a membrane lattice permit analysis of membrane
     effect on elastic network models in *exGNM* method described in [TL12]_.

    .. [TL12] Lezon TR, Bahar I, Constraints Imposed by the Membrane
       Selectively Guide the Alternating Access Dynamics of the Glutamate
       Transporter GltPh

    """

    def __init__(self, name='Unknown'):

        super(exGNM, self).__init__(name)
        self._membrane = None
        self._combined = None

    def buildMembrane(self, coords, **kwargs):
        """Build Kirchhoff matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg membrane_high: the maximum z coordinate of the membrane. Default is **13.0**
        :type membrane_high: float

        :arg membrane_low: the minimum z coordinate of the membrane. Default is **-13.0**
        :type membrane_low: float

        :arg R: radius of all membrane in x-y direction. Default is **80**
        :type R: float

        :arg Ri: inner radius of the membrane in x-y direction if it needs to be hollow. 
                 Default is **0**, which is not hollow
        :type Ri: float

        :arg r: radius of each membrane node. Default is *.1**
        :type r: float
        
        :arg lat: lattice type which could be **FCC** (face-centered-cubic, default), 
                  **SC** (simple cubic), **SH** (simple hexagonal)
        :type lat: str

        :arg exr: exclusive radius of each protein node. Default is **5.0**
        :type exr: float

        :arg hull: whether use convex hull to determine the protein's interior. 
                   Turn it off if protein is multimer. Default is **True**
        :type hull: bool

        :arg center: whether transform the structure to the origin (only x- and y-axis). 
                     Default is **True**
        :type center: bool
        """
        
        atoms = coords

        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')

        self._n_atoms = int(coords.shape[0])

        LOGGER.timeit('_membrane')

        depth = kwargs.pop('depth', None)
        h = depth / 2 if depth is not None else None
            
        h = kwargs.pop('h', h)
        if h is not None:
            h = float(h)
            hu = h
            hl = -h
        else:
            hu = kwargs.pop('membrane_high', 13.0)
            hu = kwargs.pop('high', hu)
            hu = float(hu)
            
            hl = kwargs.pop('membrane_low', -13.0)
            hl = kwargs.pop('low', hl)
            hl = float(hl)

        R = float(kwargs.pop('R', 80.))
        Ri = float(kwargs.pop('Ri', 0.))
        r = float(kwargs.pop('r', 3.1))
        lat = str(kwargs.pop('lat', 'FCC'))
        exr = float(kwargs.pop('exr', 5.))
        use_hull = kwargs.pop('hull', True)
        centering = kwargs.pop('center', True)
        
        V = assign_lpvs(lat)

        if centering:
            c0 = coords.mean(axis=0)
            c0[-1] = 0.
            coords -= c0
        # determine transmembrane part
        torf = np.logical_and(coords[:, -1] < hu, coords[:, -1] > hl)
        transmembrane = coords[torf, :]

        if not np.any(torf):
            raise ValueError('No region was identified as membrane. Please use a structure from opm/ppm.')

        if use_hull:
            from scipy.spatial import ConvexHull
            hull = ConvexHull(transmembrane)
        else:
            hull = transmembrane

        ## determine the bound for ijk
        imax = (R + V[0,2] * (hu - hl)/2.)/r
        jmax = (R + V[1,2] * (hu - hl)/2.)/r
        kmax = (R + V[2,2] * (hu - hl)/2.)/r    

        imax = int(ceil(imax))
        jmax = int(ceil(jmax))
        kmax = int(ceil(kmax))

        membrane = []
        atm = 0
        for i in range(-imax, imax):
            for j in range(-jmax, jmax):
                for k in range(-kmax, kmax):
                    c = array([i, j, k])
                    xyz = 2.*r*dot(c, V)
                    
                    if xyz[2]>hl and xyz[2]<hu and \
                       xyz[0]>-R and xyz[0]<R and \
                       xyz[1]>-R and xyz[1]<R:
                        dd = norm(xyz[:2])
                        if dd < R and dd > Ri:
                            if checkClash(xyz, hull, radius=exr):
                                membrane.append(xyz)
                                atm = atm + 1 

        membrane = array(membrane)

        if len(membrane) == 0:
            self._membrane = None
            LOGGER.warn('no membrane is built. The protein should be transformed to the correct origin as in OPM')
            return coords
        else:
            self._membrane = AtomGroup(title="Membrane")
            self._membrane.setCoords(membrane)
            self._membrane.setResnums(range(atm))
            self._membrane.setResnames(["NE1" for i in range(atm)])
            self._membrane.setChids(["Q" for i in range(atm)])
            self._membrane.setElements(["Q1" for i in range(atm)])
            self._membrane.setNames(["Q1" for i in range(atm)])
            LOGGER.report('Membrane was built in %2.fs.', label='_membrane')

            coords = self._combineMembraneProtein(atoms)
            return coords

    def buildKirchhoff(self, coords, cutoff=15., gamma=1., **kwargs):
        """Build Kirchhoff matrix for given coordinate set. 
        **kwargs** are passed to :method:`.buildMembrane`.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float
        """

        atoms = coords
        turbo = kwargs.pop('turbo', True)

        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')

        n_atoms = int(coords.shape[0])

        if self._membrane is None:
            coords = self.buildMembrane(atoms, **kwargs)
        else:
            coords = self._combined.getCoords()

        system = zeros(coords.shape[0], dtype=bool)
        system[:n_atoms] = True

        LOGGER.timeit('_exgnm')

        if turbo:
            self._kirchhoff = buildReducedKirchhoff(coords, system, cutoff, gamma, **kwargs)
        else:
            super(exGNM, self).buildKirchhoff(coords, cutoff, gamma, **kwargs)
            system = np.repeat(system, 3)
            self._kirchhoff = _reduceModel(self._kirchhoff, system)

        LOGGER.report('Kirchhoff was built in %.2fs.', label='_exgnm')
        self._dof = self._kirchhoff.shape[0]
        self._n_atoms = n_atoms
    
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh`
        function to diagonalize the Kirchhoff matrix. When Scipy is not found,
        :func:`numpy.linalg.eigh` is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
            If **None** is given, all modes will be calculated.
        :type n_modes: int or None, default is 20

        :arg zeros: If **True**, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is **True**

        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is **True**
        """

        super(exGNM, self).calcModes(n_modes, zeros, turbo)

    def getMembrane(self):
        """Returns a copy of the membrane coordinates."""

        return self._membrane

    def getCombined(self):
        """Returns a copy of the combined atoms or coordinates."""

        return self._combined

    def _combineMembraneProtein(self, coords):
        if self._membrane is not None:
            if isinstance(coords, Atomic):
                self._combined = coords.copy() + self._membrane
                coords = self._combined.getCoords()
            else:
                self._combined = coords = np.concatenate((coords, self._membrane.getCoords()), axis=0)
        else:
            self._combined = coords = copy(coords)
        return coords
    

def assign_lpvs(lat):
    """ Given lattice type return 3 lattice primitive vectors"""
    lpv = zeros((3,3))
    if lat=='FCC':
        lpv[0,1]=1./sqrt(2)
        lpv[0,2]=1./sqrt(2)
        lpv[1,0]=1./sqrt(2)
        lpv[1,2]=1./sqrt(2)
        lpv[2,0]=1./sqrt(2)
        lpv[2,1]=1./sqrt(2)
    elif lat=='SC':
        lpv[0,0]=1
        lpv[1,1]=1
        lpv[2,2]=1
    elif lat=='SH':
        lpv[0,0]=1./2
        lpv[0,1]=-sqrt(3)/2
        lpv[1,0]=1./2
        lpv[1,1]=sqrt(3)/2
        lpv[2,2]=1.
    return lpv

def checkClash(node, hull, radius=5.):
    """ Check there is a clash between given coordinate and all pdb coordinates.
    **False** for clashing and **True** for not clashing."""

    if isinstance(hull, np.ndarray):
        H = hull
        ishull = False
    else:
        H = hull.points
        ishull = True

    lb = H.min(axis=0) - radius
    ub = H.max(axis=0) + radius
    if np.all(node > ub) or np.all(node < lb):
        return True
    
    if ishull:
        in_hull = all([dot(eq[:-1], node) + eq[-1] <= 0 
                    for eq in hull.equations])

        if in_hull:
            return False

    for coord in H:
        if norm(node-coord)<radius:
            return False
    return True

def peelr(coords, system, r0=20., dr=20.):
    n_sys_atoms = int(system.sum())
    n_atoms = len(system)
    labels = np.zeros(n_atoms, dtype=int)

    # identify system beads
    sys_coords = coords[system, :2]
    sys_norms = norm(sys_coords, axis=1)
    sys_r = max(sys_norms)
    r0 += sys_r

    # label environment beads
    env_coords = coords[~system, :2]
    env_norms = norm(env_coords, axis=1)
    L = (env_norms - r0) // dr + 1
    L = np.clip(L, 0, None) + 1
    labels[n_sys_atoms:] = L
    
    uniq_labels = np.unique(labels)
    if len(uniq_labels) >= 3:
        uniq_labels.sort()
        lbl_last = uniq_labels[-1]
        lbl_2nd_last = uniq_labels[-2]

        n_last = np.sum(labels == lbl_last)
        n_2nd_last = np.sum(labels == lbl_2nd_last)

        if n_last < 0.2 * n_2nd_last:
            LOGGER.debug('edge nodes detected (%d/%d)'%(n_2nd_last, n_last))
            labels[labels == lbl_last] = lbl_2nd_last

    if len(uniq_labels) >= 3:
        uniq_labels.sort()
        lbl_first = uniq_labels[1]
        lbl_2nd = uniq_labels[2]

        n_first = np.sum(labels == lbl_first)
        n_2nd = np.sum(labels == lbl_2nd)

        if n_first < 0.2 * n_2nd:
            LOGGER.debug('inner nodes detected (%d/%d)'%(n_2nd, n_first))
            labels[labels == lbl_first] = lbl_2nd
    if not any(uniq_labels == 1):
        LOGGER.debug('no layer inside the system')
        for i in range(len(labels)):
            if labels[i] > 1:
                labels[i] -= 1

    uniq_labels = np.unique(labels)
    for i, label in enumerate(uniq_labels):
        labels[labels==label] = i

    return labels

def buildReducedKirchhoff(coords, system, cutoff=15., gamma=1.0, **kwargs):
    
    r0 = kwargs.pop('r0', 20.)
    dr = kwargs.pop('dr', 20.)
    labels = peelr(coords, system, r0, dr)
    LOGGER.debug('layers: ' + str(np.unique(labels)))

    G = calcKirchhoffRecursion(coords, labels, 0, cutoff=cutoff, gamma=gamma, **kwargs)
    return G

def calcKirchhoffRecursion(coords, layers, layer, cutoff=15., gamma=1.0, **kwargs):
    if layer == 0:
        LOGGER.debug('max layer: %d'%max(layers))
    LOGGER.debug('layer: %d'%layer)
    Gss, Gse = buildLayerKirchhoff(coords, layers, layer, cutoff=cutoff, gamma=gamma, **kwargs)

    if Gse is None: # last layer, Gee=Gss
        G = Gss
    else:
        Gee = calcKirchhoffRecursion(coords, layers, layer+1, cutoff=cutoff, gamma=gamma, **kwargs)
        Cee = inv(Gee)
        #G = Gss - Gse.dot(Cee.dot(Gse.T))
        #G = Gss - Gse @ Cee @ Gse.T
        if PY3K:
            G = Gss - Gse.__matmul__(Cee).__matmul__(Gse.T)
        else:
            G = Gss - Gse.dot(Cee.dot(Gse.T))
    LOGGER.debug('layer: %d finished'%layer)
    return G

def buildLayerKirchhoff(coords, layers, layer, cutoff=15., gamma=1.0, **kwargs):
    gmem = kwargs.pop('gamma_memb', gamma)

    torf_inner = layers == (layer-1)
    torf_sys = layers == layer
    torf_env = layers == (layer+1)
    
    coords_inner = coords[torf_inner, :] # inner layer coords
    coords_sys = coords[torf_sys, :] # sys coords
    coords_env = coords[torf_env, :] # env coords
    
    n_sys_atoms = len(coords_sys)
    n_env_atoms = len(coords_env)
    n_inner_atoms = len(coords_inner)
    
    Gss = np.zeros((n_sys_atoms, n_sys_atoms))
    Gse = np.zeros((n_sys_atoms, n_env_atoms)) if n_env_atoms else None
    
    cutoff2 = cutoff * cutoff
    for i in range(n_sys_atoms):
        coordi = coords_sys[i, :]
        I = slice(i, (i+1))

        # sys-sys
        for j in range(i+1, n_sys_atoms):
            coordj = coords_sys[j, :]
            J = slice(j, (j+1))

            v = coordi - coordj
            dist2 = np.inner(v, v)
            if dist2 < cutoff2:
                g = gamma if layer == 0 else gmem
                Gss[I, J] = Gss[J, I] = -g
                Gss[I, I] += g
                Gss[J, J] += g

        # sys-env
        for k in range(n_env_atoms):
            coordk = coords_env[k, :]
            K = slice(k, (k+1))

            v = coordi - coordk
            dist2 = np.inner(v, v)
            if dist2 < cutoff2:
                g = gmem 
                Gse[I, K] = -g
                Gss[I, I] += g

        # sys-inner
        for k in range(n_inner_atoms):
            coordk = coords_inner[k, :]
            K = slice(k, (k+1))

            v = coordi - coordk
            dist2 = np.inner(v, v)
            if dist2 < cutoff2:
                g = gmem
                Gss[I, I] += g
    return Gss, Gse
