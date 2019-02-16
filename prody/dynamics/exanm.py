# -*- coding: utf-8 -*-
"""This module defines a class and a function for explicit membrane ANM calculations."""

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB, writePDB
from prody.utilities import importLA, checkCoords, copy
from prody.kdtree import KDTree
from numpy import sqrt, zeros, linalg, min, max, mean, array, ceil, dot
from numpy.linalg import norm, inv

from .anm import ANMBase, calcANM, ANM
from .gnm import checkENMParameters
from .editing import reduceModel

__all__ = ['exANM']

class Increment(object):

    def __init__(self, s=0):

        self._i = s

    def __call__(self, i=1):

        self._i += i
        return self._i


class exANM(ANM):

    """Class for explicit ANM (exANM) method ([FT00]_).
    Optional arguments build a membrane lattice permit analysis of membrane
     effect on elastic network models in *exANM* method described in [TL12]_.

    .. [TL12] Lezon TR, Bahar I, Constraints Imposed by the Membrane
       Selectively Guide the Alternating Access Dynamics of the Glutamate
       Transporter GltPh

    """

    def __init__(self, name='Unknown'):

        super(exANM, self).__init__(name)
        self._membrane = None
        self._combined = None

    def buildMembrane(self, coords, **kwargs):
        """Build Hessian matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg membrane_hi: the maximum z coordinate of the membrane. Default is **13.0**
        :type membrane_hi: float

        :arg membrane_lo: the minimum z coordinate of the membrane. Default is **-13.0**
        :type membrane_lo: float

        :arg R: radius of all membrane in x-y direction. Default is **80**
        :type R: float

        :arg r: radius of each membrane node. Default is **3.1**
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

        self._n_atoms = natoms = int(coords.shape[0])

        LOGGER.timeit('_membrane')

        h = kwargs.pop('h', None)
        if h is not None:
            h = float(h)
            hu = h
            hl = -h
        else:
            hu = float(kwargs.pop('membrane_hi', 13.0))
            hl = float(kwargs.pop('membrane_lo', -13.0))
        R = float(kwargs.pop('R', 80.))
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
                        if dd<R:
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

    def buildHessian(self, coords, cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set. 
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

        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')

        self._n_atoms = natoms = int(coords.shape[0])

        if self._membrane is None:
            coords = self.buildMembrane(atoms, **kwargs)
        else:
            coords = self._combined.getCoords()

        LOGGER.timeit('_exanm')

        total_natoms = int(coords.shape[0])
        self._hessian = np.zeros((natoms*3, natoms*3), float)
        total_hessian = np.zeros((total_natoms*3, total_natoms*3), float)
        cutoff, g, gamma = checkENMParameters(cutoff, gamma)
        cutoff2 = cutoff * cutoff
        for i in range(total_natoms):
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
                total_hessian[res_i3:res_i33, res_j3:res_j33] = super_element
                total_hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                total_hessian[res_i3:res_i33, res_i3:res_i33] = total_hessian[res_i3:res_i33, res_i3:res_i33] - super_element
                total_hessian[res_j3:res_j33, res_j3:res_j33] = total_hessian[res_j3:res_j33, res_j3:res_j33] - super_element

        ss = total_hessian[:natoms*3, :natoms*3]
        so = total_hessian[:natoms*3, natoms*3:]
        os = total_hessian[natoms*3:,:natoms*3]
        oo = total_hessian[natoms*3:, natoms*3:]
        self._hessian = ss - np.dot(so, np.dot(linalg.inv(oo), os))
        LOGGER.report('Hessian was built in %.2fs.', label='_exanm')
        self._dof = self._hessian.shape[0]
    
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh`
        function to diagonalize the Hessian matrix. When Scipy is not found,
        :func:`numpy.linalg.eigh` is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
            If **None** is given, all modes will be calculated.
        :type n_modes: int or None, default is 20

        :arg zeros: If **True**, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is **True**

        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is **True**
        """

        super(exANM, self).calcModes(n_modes, zeros, turbo)

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
        if linalg.norm(node-coord)<radius:
            return False
    return True



def test2(pdb='2nwl-mem.pdb'):
    from prody import parsePDB
    structure = parsePDB(pdb, subset='ca')
    exanm = exANM('2nwl')
    exanm.buildHessian(structure)
    exanm.calcModes()
    return exanm
