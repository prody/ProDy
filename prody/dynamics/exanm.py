# -*- coding: utf-8 -*-
"""This module defines a class and a function for explicit membrane ANM calculations."""

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB
from prody.utilities import importLA, checkCoords
from prody.kdtree import KDTree
from numpy import sqrt, zeros, linalg, min, max

from .anm import ANMBase, calcANM
from .gnm import checkENMParameters
from .editing import reduceModel

__all__ = ['exANM','writeMembranePDB']

class Increment(object):

    def __init__(self, s=0):

        self._i = s

    def __call__(self, i=1):

        self._i += i
        return self._i


class exANM(ANMBase):

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

    def buildHessian(self, coords, cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float

        :arg membrane_hi: the maximum z coordinate of the pdb default is 13.0
        :type membrane_hi: float

        :arg membrane_lo: the minimum z coordinate of the pdb default is -13.0
        :type membrane_lo: float

        :arg R: radius of all membrane in x-y direction default is 80. 
        :type R: float

        :arg r: radius of individual barrel-type membrane protein default is 2.5.
        :type 
        
        :arg lat: lattice type which could be FCC(face-centered-cubic)(default), 
        SC(simple cubic), SH(simple hexagonal)
        :type lat: str
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

        self._n_atoms = natoms = int(coords.shape[0])

        center = mean(coords, axis=0)
        coords = coords - center

        pxlo = min(np.append(coords[:,0],10000))
        pxhi = max(np.append(coords[:,0],-10000))
        pylo = min(np.append(coords[:,1],10000))
        pyhi = max(np.append(coords[:,1],-10000))
        pzlo = min(np.append(coords[:,2],10000))
        pzhi = max(np.append(coords[:,2],-10000))

        membrane_hi = float(kwargs.get('membrane_hi', 13.0))
        membrane_lo = float(kwargs.get('membrane_lo', -13.0))
        R = float(kwargs.get('R', 80))
        r = float(kwargs.get('r', 2.5))
        lat = str(kwargs.get('lat', 'FCC'))
        lpv = assign_lpvs(lat)

        imax = R/r
        jmax = R/r
        kmax = membrane_hi/r
        kmin = membrane_lo/r

        self._membrane = membrane = zeros((1,3))

        LOGGER.timeit('_membrane')
        coords = pdb.getCoords()
        membrane = zeros((1,3))
        membrane1 = zeros((1,3))
        atm = 0
        a = 0
        b = 0
        scale = 5
        lim = range(-int(scale/2), int(scale/2)+1)
        for i in range(-int(imax),int(imax+1)):
            for j in range(-int(jmax),int(jmax+1)):
                for k in range(int(kmin),int(kmax+1)):
                    X = zeros((1,3))
                    X[0,0]=-R+(2*i+1)*r
                    X[0,1]=-R+(2*j+1)*r
                    X[0,2]=-membrane_lo+(2*k+1)*r
                    dd=0
                    for p in range(3):
                        dd += X[0,p] ** 2
                    if dd<R**2 and X[0,2]>membrane_lo and X[0,2]<membrane_hi:
                        if X[0,2]-6<membrane_lo or k==int(kmin):
                            if X[0,0]>pxlo and X[0,0]<pxhi and X[0,1]>pylo and X[0,1]<pyhi and X[0,2]>pzlo and X[0,2]<pzhi:
                                if checkClash(X, coords[:natoms,:], radius=5):
                                    if atm ==0:
                                        membrane = X
                                    else:
                                        membrane = np.append(membrane, X, axis=0)        
                                    atm = atm + 1
                                    coords = np.append(coords, X, axis=0)
                            else:
                                if atm==0:
                                    membrane = X
                                else:
                                    membrane = np.append(membrane, X, axis=0)
                                atm = atm + 1
                                coords = np.append(coords, X, axis=0)
                        elif X[0,2]+6>membrane_hi or k==int(kmax):
                            if X[0,0]>pxlo and X[0,0]<pxhi and X[0,1]>pylo and X[0,1]<pyhi and X[0,2]>pzlo and X[0,2]<pzhi:
                                if checkClash(X, coords[:natoms,:], radius=5):
                                    if atm ==0:
                                        membrane = X
                                    else:
                                        membrane = np.append(membrane, X, axis=0)        
                                    atm = atm + 1
                                    coords = np.append(coords, X, axis=0)
                            else:
                                if atm==0:
                                    membrane = X
                                else:
                                    membrane = np.append(membrane, X, axis=0)
                                atm = atm + 1
                                coords = np.append(coords, X, axis=0)
                        else:
                            T1 = zeros((scale,3))
                            for l in range(scale):
                                T1[l,0]=X[0,0]
                                T1[l,1]=X[0,1]
                                T1[l,2]=X[0,2]-(1.*(lim[l])/(scale-1))*r
                            if X[0,0]>pxlo and X[0,0]<pxhi and X[0,1]>pylo and X[0,1]<pyhi and X[0,2]>pzlo and X[0,2]<pzhi:
                                if checkClash(X, coords[:natoms,:], radius=5):
                                    if atm ==0:
                                        membrane1 = X
                                    else:
                                        membrane1 = np.append(membrane1, X, axis=0)        
                                    atm = atm + 1
                                    coords = np.append(coords, X, axis=0)
                            else:
                                if atm==0:
                                    membrane1 = X
                                else:
                                    membrane1 = np.append(membrane1, X, axis=0)
                                atm = atm + 1
                                coords = np.append(coords, X, axis=0)
                            for l in range(scale):
                                if T1[l,0]>pxlo and T1[l,0]<pxhi and T1[l,1]>pylo and T1[l,1]<pyhi and T1[l,2]>pzlo and T1[l,2]<pzhi:
                                    if checkClash(T1[l,:], coords[:natoms,:], radius=5):
                                        membrane1 = np.append(membrane1, T1[l,:].reshape(1,3), axis=0)        
                                        atm = atm + 1
                                        coords = np.append(coords, T1[l,:].reshape(1,3), axis=0)
                                else:
                                    membrane1 = np.append(membrane1, T1[l,:].reshape(1,3), axis=0)
                                    atm = atm + 1
                                    coords = np.append(coords, T1[l,:].reshape(1,3), axis=0)

        length = membrane.shape[0]
        f = open(filename, 'w')
        for i in range(length):
            f.write('ATOM%7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n' % (i,i,membrane[i,0],membrane[i,1],membrane[i,2]))
        f.close()     
        length = membrane1.shape[0]
        f = open(filename2, 'w')
        for i in range(length):
            f.write('ATOM%7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n' % (i,i,membrane1[i,0],membrane1[i,1],membrane1[i,2]))
        f.close()                     

        self._membrane = membrane 
        LOGGER.report('Membrane was built in %2.fs.', label='_membrane')
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
        so = total_hessian[:natoms*3, natoms*3+1:]
        os = total_hessian[natoms*3+1:,:natoms*3]
        oo = total_hessian[natoms*3+1:, natoms*3+1:]
        self._hessian = ss - np.dot(so, np.dot(linalg.inv(oo), os))
        LOGGER.report('Hessian was built in %.2fs.', label='_exanm')
    self._dof = self._hessian.shape[0]
    
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh`
        function to diagonalize the Hessian matrix. When Scipy is not found,
        :func:`numpy.linalg.eigh` is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
            If ``None`` is given, all modes will be calculated.
        :type n_modes: int or None, default is 20

        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``

        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        """

        super(exANM, self).calcModes(n_modes, zeros, turbo)

    def getMembrane(self):
        """Return a copy of the membrane coordinates."""

        if self._membrane is not None:
            return self._membrane.copy()

    def _getMembrane(self):

        return self._membrane
    

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

def checkClash(coordinates, pdb_coords, radius):
    """ Check there is a clash between given coordinate and all pdb coordinates."""
    for i in range(pdb_coords.shape[0]):
        if linalg.norm(coordinates-pdb_coords[i])<radius:
            return False
    return True

def writeMembranePDB(filename, membrane):
    """ Given membrane coordinates it will write a pdb file with membrane coordinates. 
        :arg filename: filename for the pdb file. 
        :type filename: str

        :arg membrane: membrane coordinates or the membrane structure. 
        :type membrane: nd.array
        """
    length = membrane.shape[0]
    f = open(filename, 'w')
    for i in range(length):
        f.write('ATOM%7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n' % (i,i,membrane[i,0],membrane[i,1],membrane[i,2]))
    f.close()

def test2(pdb='2nwl-mem.pdb'):
    from prody import parsePDB
    structure = parsePDB(pdb, subset='ca')
    exanm = exANM('2nwl')
    exanm.buildHessian(structure)
    exanm.calcModes()
    return exanm
