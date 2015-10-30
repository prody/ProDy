# -*- coding: utf-8 -*-
"""This module defines a class and a function for rotating translating blocks
(RTB) calculations."""

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB
from prody.utilities import importLA, checkCoords

from .anm import ANMBase
from .gnm import GNMBase, ZERO, checkENMParameters
from numpy import eye, arccos, zeros, linalg, dot, tan, sqrt, pi


__all__ = ['bbENM']

class Increment(object):

    def __init__(self, s=0):

        self._i = s

    def __call__(self, i=1):

        self._i += i
        return self._i


class bbENM(ANMBase):

    """Class for bond bending correction for anisotropic network model ([AS12]_).
    Additional argument is the b parameter which is the bond bending constant in the model. 
    .. [AS12] Srivastava A, Halevi RB, Veksler A, Granek, R. Tensorial elastic network model 
    for protein dynamics: Integration of the anisotropic model with bond-bending and twist
    elasticities.
       *Proteins* **2012** 80:2692-2700.

    """

    def __init__(self, name='Unknown'):

        super(bbENM, self).__init__(name)


    def buildHessian(self, coords, B=1., cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg b: bond-bending constant, default is 1.0
        :type b: float

        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float

        :type scale: float

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

        LOGGER.timeit('_bbenm')
        self._n_atoms = natoms = int(coords.shape[0])

        self._hessian = hessian = np.zeros((3*natoms, 3*natoms), float)
        self._dof = 3*natoms - 6
<<<<<<< HEAD
        #vectorized version of all bonds.
        b = zeros((natoms,natoms,3))
        for i in range(natoms-1):
            for j in range(i+1, natoms):
                b[i,j]=coords[i,:]-coords[j,:]
                b[j,i]=-b[i,j]

        #kirchoff matrix calculation
        kirchoff = eye(natoms)
        for i in range(natoms-1):
            for j in range(i+1,natoms):
                kirchoff[i,j]=sqrt(linalg.norm(b[i,j])) < 15

        #theta zero calculation
        theta= zeros((natoms, natoms, natoms))
        for i in  range(natoms-1):
            for j in range(i+1,natoms):
                for k in range(i+1,natoms):
                    if j!=k:
                        if kirchoff[i,j] and kirchoff[j,k]:
                            theta[i,j,k]=arccos(dot(b[i,j],b[j,k])/linalg.norm(b[i,j])/linalg.norm(b[j,k]))
                            theta[k,j,i]=theta[i,j,k]
        # anm hessian calculation 
=======
        cutoffd = cutoff
        gammad = gamma
        
        #anm hessian calculation 
>>>>>>> prody/master
        cutoff, g, gamma = checkENMParameters(cutoff, gamma)
        cutoff2 = cutoff * cutoff
        for i in range(natoms):
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

<<<<<<< HEAD
        #hessian updates
        for i in range(natoms-1):
            for j in range(i+1, natoms):
                for k in range(i+1, natoms):
                    if j!= k:
                        if kirchoff[i,j] and kirchoff[j,k]:
                            denom = linalg.norm(b[i,j])**2*linalg.norm(b[j,k])**2
                            if theta[i,j,k] < 1e-5 and theta[i,j,k]> -1e-5:
                                #i,i update
                                hessian[3*i, 3*i]+=0.5*(b[j,k][1]**2+b[j,k][2]**2)/denom
                                hessian[3*i, 3*i+1]+=0.5*(-2*b[j,k][0]*b[j,k][1])/denom
                                hessian[3*i, 3*i+2]+=0.5*(-2*b[j,k][0]*b[j,k][2])/denom
                                hessian[3*i+1, 3*i+1]+=0.5*(b[j,k][0]**2+b[j,k][2]**2)/denom
                                hessian[3*i+2, 3*i+2]+=0.5*(b[j,k][0]**2+b[j,k][1]**2)/denom
                                hessian[3*i+1, 3*i+2]+=0.5*(-2*b[j,k][1]*b[j,k][2])/denom
                                #j,j update
                                hessian[3*j, 3*j]+=0.5*(b[i,j][1]**2+b[i,j][2]**2)/denom
                                hessian[3*j, 3*j+1]+=0.5*(-2*b[i,j][0]*b[i,j][1])/denom
                                hessian[3*j, 3*j+2]+=0.5*(-2*b[i,j][0]*b[i,j][2])/denom
                                hessian[3*j+1, 3*j+1]+=0.5*(b[i,j][0]**2+b[i,j][2]**2)/denom
                                hessian[3*j+2, 3*j+2]+=0.5*(b[i,j][0]**2+b[i,j][1]**2)/denom
                                hessian[3*j+1, 3*j+2]+=0.5*(-2*b[i,j][1]*b[i,j][2])/denom
                                #k,k update
                                hessian[3*k, 3*k]+=0.5*(b[i,j][1]**2-2*b[i,j][1]*b[j,k][1]+b[j,k][1]**2 \
                                    +b[i,j][2]**2-2*b[i,j][2]*b[j,k][2]+b[j,k][2]**2)/denom
                                hessian[3*k, 3*k+1]+=0.5*(2*b[i,j][1]*b[j,k][0]-b[i,j][0]*b[i,j][1] \
                                    +b[i,j][0]*b[j,k][1]-b[j,k][0]*b[j,k][1])/denom
                                hessian[3*k, 3*k+2]+=0.5*(2*b[i,j][2]*b[j,k][0]-b[i,j][0]*b[j,k][2] \
                                    +b[i,j][0]*b[j,k][2]-b[j,k][0]*b[j,k][2])/denom
                                hessian[3*k+1, 3*k+1]+=0.5*(b[i,j][0]**2-2*b[i,j][0]*b[j,k][0] +b[j,k][0]**2 \
                                    +b[i,j][0]**2-2*b[i,j][2]*b[j,k][2]+b[j,k][2]**2)/denom
                                hessian[3*k+2, 3*k+2]+=0.5*(b[i,j][0]**2-2*b[i,j][0]*b[j,k][0] +b[j,k][0]**2 \
                                    +b[i,j][1]**2-2*b[i,j][1]*b[j,k][1]+b[j,k][1]**2)/denom
                                hessian[3*k+1, 3*k+2]+=(b[i,j][2]*b[j,k][1]-b[i,j][1]*b[i,j][2]+b[i,j][1]*b[j,k][2] \
                                    -b[j,k][1]*b[j,k][2])/denom
                                #i,j update
                                hessian[3*i, 3*j]-=(b[i,j][1]*b[j,k][1]+b[i,j][2]*b[j,k][2])/denom
                                hessian[3*i+1, 3*j+1]-=(b[i,j][0]*b[j,k][0]+b[i,j][2]*b[j,k][2])/denom
                                hessian[3*i+2, 3*j+2]-=(b[i,j][0]*b[j,k][0]+b[i,j][1]*b[j,k][1])/denom
                                hessian[3*i, 3*j+1]+=(b[i,j][0]*b[j,k][1])/denom
                                hessian[3*i, 3*j+2]+=(b[i,j][0]*b[j,k][2])/denom
                                hessian[3*i+1, 3*j]+=(b[i,j][1]*b[j,k][0])/denom
                                hessian[3*i+1, 3*j+2]+=(b[i,j][1]*b[j,k][2])/denom
                                hessian[3*i+2, 3*j]+=(b[i,j][2]*b[j,k][0])/denom
                                hessian[3*i+2, 3*j+1]+=(b[i,j][2]*b[j,k][1])/denom
                                #i,k update
                                hessian[3*i, 3*k]+=(b[i,j][1]*b[j,k][1]-b[j,k][1]**2+b[i,j][2]*b[j,k][2]-b[j,k][2]**2)/denom
                                hessian[3*i+1, 3*k+1]+=(b[i,j][0]*b[j,k][0]-b[j,k][0]**2+b[i,j][2]*b[j,k][2]-b[j,k][2]**2)/denom
                                hessian[3*i+2, 3*k+2]+=(b[i,j][0]*b[j,k][0]-b[j,k][0]**2+b[i,j][1]*b[j,k][1]-b[j,k][1]**2)/denom
                                hessian[3*i, 3*k+1]+=(b[j,k][0]*b[j,k][1]-b[i,j][0]*b[j,k][1])/denom
                                hessian[3*i, 3*k+2]+=(b[j,k][0]*b[j,k][2]-b[i,j][0]*b[j,k][2])/denom
                                hessian[3*i+1, 3*k+2]+=(b[j,k][1]*b[j,k][2]-b[i,j][1]*b[j,k][2])/denom
                                hessian[3*i+1, 3*k]+=(b[j,k][0]*b[j,k][1]-b[i,j][1]*b[j,k][0])/denom
                                hessian[3*i+2, 3*k]+=(b[j,k][0]*b[j,k][2]-b[i,j][2]*b[j,k][0])/denom
                                hessian[3*i+2, 3*k+1]+=(b[j,k][1]*b[j,k][2]-b[i,j][2]*b[j,k][1])/denom
                                #j,k update
                                hessian[3*j, 3*k]+=(b[i,j][1]*b[j,k][1]+b[i,j][2]*b[j,k][2]-b[i,j][1]**2-b[i,j][2]**2)/denom
                                hessian[3*j+1, 3*k+1]+=(b[i,j][0]*b[j,k][0]+b[i,j][2]*b[j,k][2]-b[i,j][0]**2-b[i,j][2]**2)/denom
                                hessian[3*j+2, 3*k+2]+=(b[i,j][0]*b[j,k][0]+b[i,j][1]*b[j,k][1]-b[i,j][0]**2-b[i,j][1]**2)/denom
                                hessian[3*j, 3*k+1]+=(b[i,j][0]*b[i,j][1]-b[i,j][1]*b[j,k][0])/denom
                                hessian[3*j, 3*k+2]+=(b[i,j][0]*b[i,j][2]-b[i,j][2]*b[j,k][0])/denom
                                hessian[3*j+1, 3*k+2]+=(b[i,j][1]*b[i,j][2]-b[i,j][1]*b[j,k][2])/denom
                                hessian[3*j+1, 3*k]+=(b[i,j][0]*b[i,j][1]-b[i,j][0]*b[j,k][1])/denom
                                hessian[3*j+2, 3*k]+=(b[i,j][0]*b[i,j][2]-b[i,j][0]*b[j,k][2])/denom
                                hessian[3*j+2, 3*k+1]+=(b[i,j][1]*b[i,j][2]-b[i,j][1]*b[j,k][2])/denom
                            elif theta[i,j,k] - pi/2 < 1e-5 and theta[i,j,k] -pi/2 > -1e-5:
                                #ii update
                                hessian[3*i, 3*i]+=0.5*(b[j,k][0]**2)/denom
                                hessian[3*i, 3*i+1]+=0.5*(2*b[j,k][0]*b[j,k][1])/denom
                                hessian[3*i, 3*i+2]+=0.5*(2*b[j,k][0]*b[j,k][2])/denom
                                hessian[3*i+1, 3*i+1]+=0.5*(b[j,k][1]**2)/denom
                                hessian[3*i+2, 3*i+2]+=0.5*(b[j,k][2]**2)/denom
                                hessian[3*i+1, 3*i+2]+=0.5*(2*b[j,k][1]*b[j,k][2])/denom
                                #j,j update
                                hessian[3*j, 3*j]+=0.5*(b[i,j][0]**2)/denom
                                hessian[3*j, 3*j+1]+=0.5*(2*b[i,j][0]*b[i,j][1])/denom
                                hessian[3*j, 3*j+2]+=0.5*(2*b[i,j][0]*b[i,j][2])/denom
                                hessian[3*j+1, 3*j+1]+=0.5*(b[i,j][1]**2)/denom
                                hessian[3*j+2, 3*j+2]+=0.5*(b[i,j][2]**2)/denom
                                hessian[3*j+1, 3*j+2]+=0.5*(2*b[i,j][1]*b[i,j][2])/denom
                                #k,k update
                                hessian[3*k, 3*k]+=0.5*(b[i,j][0]**2+2*b[i,j][0]*b[j,k][0]+b[j,k][0]**2)/denom
                                hessian[3*k+1, 3*k+1]+=0.5*(b[i,j][1]**2+2*b[i,j][1]*b[j,k][1]+b[j,k][1]**2)/denom
                                hessian[3*k+2, 3*k+2]+=0.5*(b[i,j][2]**2+2*b[i,j][2]*b[j,k][2]+b[j,k][2]**2)/denom
                                hessian[3*k, 3*k+1]+=(b[i,j][0]*b[i,j][1]+b[i,j][1]*b[j,k][0] \
                                    +b[i,j][0]*b[j,k][1]+b[j,k][0]*b[j,k][1])/denom
                                hessian[3*k, 3*k+2]+=(b[i,j][0]*b[i,j][2]+b[i,j][2]*b[j,k][0] \
                                    +b[i,j][0]*b[j,k][2]+b[j,k][0]*b[j,k][2])/denom
                                hessian[3*k+1, 3*k+2]+=(b[i,j][1]*b[i,j][2]+b[i,j][2]*b[j,k][1] \
                                    +b[i,j][1]*b[j,k][2]+b[j,k][1]*b[j,k][2])/denom
                                
                                #i,j update
                                hessian[3*i, 3*j]+=(b[i,j][0]*b[j,k][0])/denom
                                hessian[3*i+1, 3*j+1]+=(b[i,j][1]*b[j,k][1])/denom
                                hessian[3*i+2, 3*j+2]+=(b[i,j][2]*b[j,k][2])/denom
                                hessian[3*i, 3*j+1]+=(b[i,j][1]*b[j,k][0])/denom
                                hessian[3*i, 3*j+2]+=(b[i,j][2]*b[j,k][0])/denom
                                hessian[3*i+1, 3*j]+=(b[i,j][0]*b[j,k][1])/denom
                                hessian[3*i+1, 3*j+2]+=(b[i,j][2]*b[j,k][1])/denom
                                hessian[3*i+2, 3*j]+=(b[i,j][0]*b[j,k][2])/denom
                                hessian[3*i+2, 3*j+1]+=(b[i,j][1]*b[j,k][2])/denom
                                #i,k update
                                hessian[3*i, 3*k]-=(b[i,j][0]*b[j,k][0]+b[j,k][0]**2)/denom
                                hessian[3*i+1, 3*k+1]-=(b[i,j][1]*b[j,k][1]+b[j,k][1]**2)/denom
                                hessian[3*i+2, 3*k+2]-=(b[i,j][2]*b[j,k][2]+b[j,k][2]**2)/denom
                                hessian[3*i, 3*k+1]-=(b[i,j][1]*b[j,k][0]+b[j,k][0]*b[j,k][1])/denom
                                hessian[3*i, 3*k+2]-=(b[i,j][2]*b[j,k][1]+b[j,k][1]*b[j,k][2])/denom
                                hessian[3*i+1, 3*k+2]-=(b[i,j][1]*b[j,k][0]+b[j,k][0]*b[j,k][1])/denom
                                hessian[3*i+1, 3*k]-=(b[i,j][0]*b[j,k][1]+b[j,k][0]*b[j,k][1])/denom
                                hessian[3*i+2, 3*k]-=(b[i,j][0]*b[j,k][2]+b[j,k][0]*b[j,k][2])/denom
                                hessian[3*i+2, 3*k+1]-=(b[i,j][1]*b[j,k][2]+b[j,k][1]*b[j,k][2])/denom
                                
                                #j,k update
                                hessian[3*j, 3*k]-=(b[j,k][0]*b[i,j][0]+b[i,j][0]**2)/denom
                                hessian[3*j+1, 3*k+1]-=(b[j,k][1]*b[i,j][1]+b[i,j][1]**2)/denom
                                hessian[3*j+2, 3*k+2]-=(b[j,k][2]*b[i,j][2]+b[i,j][2]**2)/denom
                                hessian[3*j, 3*k+1]-=(b[j,k][1]*b[i,j][0]+b[i,j][0]*b[i,j][1])/denom
                                hessian[3*j, 3*k+2]-=(b[j,k][2]*b[i,j][0]+b[i,j][0]*b[i,j][2])/denom
                                hessian[3*j+1, 3*k+2]-=(b[j,k][1]*b[i,j][2]+b[i,j][2]*b[i,j][1])/denom
                                hessian[3*j+1, 3*k]-=(b[j,k][1]*b[i,j][0]+b[i,j][0]*b[i,j][1])/denom
                                hessian[3*j+2, 3*k]-=(b[j,k][2]*b[i,j][0]+b[i,j][0]*b[i,j][2])/denom
                                hessian[3*j+2, 3*k+1]-=(b[j,k][2]*b[i,j][1]+b[i,j][1]*b[i,j][2])/denom
                            else:
                                cot=1/tan(theta[i,j,k])
                                #ii update
                                hessian[3*i, 3*i]+=0.5*cot*(b[i,j][0]**2/linalg.norm(b[i,j])**4+b[j,k][0]**2/ \
                                    dot(b[i,j],b[j,k])**2-2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*i, 3*i+1]+=0.5*cot*(2*b[i,j][0]*b[i,j][1]/linalg.norm(b[i,j])**4+2*b[j,k][0]*b[j,k][1]/ \
                                    dot(b[i,j],b[j,k])**2-2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2 \
                                    -2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*i, 3*i+2]+=0.5*cot*(2*b[i,j][0]*b[i,j][2]/linalg.norm(b[i,j])**4+2*b[j,k][0]*b[j,k][2]/ \
                                    dot(b[i,j],b[j,k])**2-2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2 \
                                    -2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*i+1, 3*i+1]+=0.5*cot*(b[i,j][1]**2/linalg.norm(b[i,j])**4+b[j,k][1]**2/ \
                                    dot(b[i,j],b[j,k])**2-2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*i+2, 3*i+2]+=0.5*cot*(b[i,j][2]**2/linalg.norm(b[i,j])**4+b[j,k][2]**2/ \
                                    dot(b[i,j],b[j,k])**2-2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*i+1, 3*i+2]+=0.5*cot*(2*b[i,j][1]*b[i,j][2]/linalg.norm(b[i,j])**4+2*b[j,k][1]*b[j,k][2]/ \
                                    dot(b[i,j],b[j,k])**2-2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2 \
                                    -2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)

                                ##jj update
                                hessian[3*j, 3*j]+=0.5*cot*(b[i,j][0]**2/dot(b[i,j],b[j,k])**2+b[j,k][0]**2/ \
                                    linalg.norm(b[j,k])**4-2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*j+1, 3*j+1]+=0.5*cot*(b[i,j][1]**2/dot(b[i,j],b[j,k])**2+b[j,k][1]**2/ \
                                    linalg.norm(b[j,k])**4-2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*j+2, 3*j+2]+=0.5*cot*(b[i,j][2]**2/dot(b[i,j],b[j,k])**2+b[j,k][2]**2/ \
                                    linalg.norm(b[j,k])**4-2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2)
                                hessian[3*j, 3*j+1]+=0.5*cot*(2*b[i,j][0]*b[i,j][1]/dot(b[i,j],b[j,k])**2+2*b[j,k][0]*b[j,k][1]/ \
                                    linalg.norm(b[j,k])**4-2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2 \
                                    -2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*j, 3*j+2]+=0.5*cot*(2*b[i,j][0]*b[i,j][2]/dot(b[i,j],b[j,k])**2+2*b[j,k][0]*b[j,k][2]/ \
                                    linalg.norm(b[j,k])**4-2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2 \
                                    -2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*j+1, 3*j+2]+=0.5*cot*(2*b[i,j][1]*b[i,j][2]/dot(b[i,j],b[j,k])**2+2*b[j,k][1]*b[j,k][2]/ \
                                    linalg.norm(b[j,k])**4-2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2 \
                                    -2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)

                                #k,k update
                                hessian[3*k, 3*k]+=0.5*cot*(b[i,j][0]**2/linalg.norm(b[i,j])**4+ \
                                    b[i,j][0]**2/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    b[j,k][0]**2/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    b[j,k][0]**2/linalg.norm(b[j,k])**4+ \
                                    2*b[i,j][0]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][0]**2/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*k+1, 3*k+1]+=0.5*cot*(b[i,j][1]**2/linalg.norm(b[i,j])**4+ \
                                    b[i,j][1]**2/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    b[j,k][1]**2/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][1]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    b[j,k][1]**2/linalg.norm(b[j,k])**4+ \
                                    2*b[i,j][1]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][1]**2/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*k+2, 3*k+2]+=0.5*cot*(b[i,j][2]**2/linalg.norm(b[i,j])**4+ \
                                    b[i,j][2]**2/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    b[j,k][2]**2/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][2]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    b[j,k][2]**2/linalg.norm(b[j,k])**4+ \
                                    2*b[i,j][2]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][2]**2/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*k, 3*k+1]+=0.5*cot*(2*b[i,j][0]*b[i,j][1]/linalg.norm(b[i,j])**4+ \
                                    2*b[i,j][0]*b[i,j][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[j,k][0]*b[j,k][1]/dot(b[i,j],b[j,k])**2- \
                                    4*b[i,j][0]*b[i,j][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[j,k][0]*b[j,k][1]/linalg.norm(b[j,k])**4+ \
                                    2*b[i,j][1]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][0]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    4*b[j,k][0]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*k, 3*k+2]+=0.5*cot*(2*b[i,j][0]*b[i,j][2]/linalg.norm(b[i,j])**4+ \
                                    2*b[i,j][0]*b[i,j][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[j,k][0]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    4*b[i,j][0]*b[i,j][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[j,k][0]*b[j,k][2]/linalg.norm(b[j,k])**4+ \
                                    2*b[i,j][2]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][0]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    4*b[j,k][0]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*k+1, 3*k+2]+=0.5*cot*(2*b[i,j][1]*b[i,j][2]/linalg.norm(b[i,j])**4+ \
                                    2*b[i,j][1]*b[i,j][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[j,k][1]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    4*b[i,j][1]*b[i,j][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[j,k][1]*b[j,k][2]/linalg.norm(b[j,k])**4+ \
                                    2*b[i,j][2]*b[j,k][1]/linalg.norm(b[i,j])/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][1]*b[j,k][2]/linalg.norm(b[i,j])/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2- \
                                    4*b[j,k][1]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                                                        
                                #i,j update
                                hessian[3*i, 3*j]+=0.5*cot*(2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][0]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][0]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+1, 3*j+1]+=0.5*cot*(2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][1]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][1]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][1]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+2, 3*j+2]+=0.5*cot*(2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][2]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][2]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][2]*b[j,k][2]/dot(b[i,j],b[j,k])**2/linalg.norm(b[j,k])**2)
                                hessian[3*i, 3*j+1]+=0.5*cot*(2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]*b[i,j][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][0]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][0]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i, 3*j+2]+=0.5*cot*(2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]*b[i,j][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][0]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][0]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+1, 3*j+2]+=0.5*cot*(2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][1]*b[i,j][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][1]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][1]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)   
                                hessian[3*i+1, 3*j]+=0.5*cot*(2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][1]*b[i,j][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][1]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][1]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+2, 3*j]+=0.5*cot*(2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][2]*b[i,j][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][2]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][2]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+2, 3*j+1]+=0.5*cot*(2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][2]*b[i,j][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2+ \
                                    2*b[i,j][2]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2- \
                                    2*b[j,k][2]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                
                                #i,k update
                                hessian[3*i, 3*k]+=0.5*cot*(-2*b[i,j][0]**2/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][0]**2/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]**2/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    4*b[i,j][0]*b[j,k][0]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][0]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+1, 3*k+1]+=0.5*cot*(-2*b[i,j][1]**2/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][1]**2/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]**2/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    4*b[i,j][1]*b[j,k][1]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][1]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][1]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+2, 3*k+2]+=0.5*cot*(-2*b[i,j][2]**2/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][2]**2/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]**2/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    4*b[i,j][2]*b[j,k][2]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][2]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][2]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i, 3*k+1]+=0.5*cot*(-2*b[i,j][0]*b[i,j][1]/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][0]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]*b[i,j][1]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][1]*b[j,k][0]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][0]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][0]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i, 3*k+2]+=0.5*cot*(-2*b[i,j][0]*b[i,j][2]/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][0]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]*b[i,j][2]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][2]*b[j,k][0]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][0]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][0]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+1, 3*k+2]+=0.5*cot*(-2*b[i,j][1]*b[i,j][2]/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][1]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]*b[i,j][2]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][2]*b[j,k][1]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][1]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][1]*b[j,k][2]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+1, 3*k]+=0.5*cot*(-2*b[i,j][1]*b[i,j][0]/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][1]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]*b[i,j][0]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][0]*b[j,k][1]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][1]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][1]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+2, 3*k]+=0.5*cot*(-2*b[i,j][2]*b[i,j][0]/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][2]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]*b[i,j][0]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][0]*b[j,k][2]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][2]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][2]*b[j,k][0]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)
                                hessian[3*i+2, 3*k+1]+=0.5*cot*(-2*b[i,j][2]*b[i,j][1]/linalg.norm(b[i,j])**4- \
                                    2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[j,k][2]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]*b[i,j][1]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][1]*b[j,k][2]/linalg.norm(b[i,j])**2/dot(b[i,j],b[j,k])- \
                                    2*b[i,j][2]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[j,k][2]*b[j,k][1]/dot(b[i,j],b[j,k])/linalg.norm(b[j,k])**2)

                                #j,k update
                                hessian[3*j, 3*k]+=0.5*cot*(-2*b[i,j][0]**2/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][0]**2/linalg.norm(b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    4*b[i,j][0]*b[j,k][0]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][0]**2/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j+1, 3*k+1]+=0.5*cot*(-2*b[i,j][1]**2/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][1]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][1]**2/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][1]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    4*b[i,j][1]*b[j,k][1]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][1]**2/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j+2, 3*k+2]+=0.5*cot*(-2*b[i,j][2]**2/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][2]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]**2/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][2]**2/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][2]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    4*b[i,j][2]*b[j,k][2]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][2]**2/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j, 3*k+1]+=0.5*cot*(-2*b[i,j][0]*b[i,j][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]*b[i,j][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][0]*b[j,k][1]/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][1]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][1]*b[j,k][0]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][0]*b[j,k][1]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][0]*b[j,k][1]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j, 3*k+2]+=0.5*cot*(-2*b[i,j][0]*b[i,j][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][0]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][0]*b[i,j][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][0]*b[j,k][2]/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][2]*b[j,k][0]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][2]*b[j,k][0]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][0]*b[j,k][2]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][0]*b[j,k][2]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j+1, 3*k+2]+=0.5*cot*(-2*b[i,j][1]*b[i,j][2]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][1]*b[j,k][2]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]*b[i,j][2]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][1]*b[j,k][2]/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][2]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][2]*b[j,k][1]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][1]*b[j,k][2]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][1]*b[j,k][2]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j+1, 3*k]+=0.5*cot*(-2*b[i,j][1]*b[i,j][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][1]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][1]*b[i,j][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][1]*b[j,k][0]/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][0]*b[j,k][1]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][0]*b[j,k][1]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][1]*b[j,k][0]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][1]*b[j,k][0]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j+2, 3*k]+=0.5*cot*(-2*b[i,j][2]*b[i,j][0]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][2]*b[j,k][0]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]*b[i,j][0]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][2]*b[j,k][0]/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][0]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][0]*b[j,k][2]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][2]*b[j,k][0]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][2]*b[j,k][0]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))
                                hessian[3*j+2, 3*k+1]+=0.5*cot*(-2*b[i,j][2]*b[i,j][1]/dot(b[i,j],b[j,k])**2- \
                                    2*b[i,j][2]*b[j,k][1]/dot(b[i,j],b[j,k])**2+ \
                                    2*b[i,j][2]*b[i,j][1]/dot(b[i,j],b[j,k])/linalg.norm(b[i,j])**2- \
                                    2*b[j,k][2]*b[j,k][1]/linalg.norm(b[j,k])**4- \
                                    2*b[i,j][1]*b[j,k][2]/linalg.norm(b[i,j])**2/linalg.norm(b[j,k])**2+ \
                                    2*b[i,j][1]*b[j,k][2]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[i,j][2]*b[j,k][1]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k])+ \
                                    2*b[j,k][2]*b[j,k][1]/linalg.norm(b[j,k])**2/dot(b[i,j],b[j,k]))

        for i in range(hessian.shape[0]-1):
            for j in range(i+1, hessian.shape[0]):
                hessian[j,i]=hessian[i,j]

        LOGGER.report('Hessian was built in %.2fs.', label='_rtb')
=======
        # hessian updates
        from .bbenmtools import buildhessian

        buildhessian(coords, hessian, natoms, 
                     float(cutoffd), float(gammad),)

        LOGGER.report('Hessian was built in %.2fs.', label='_bbenm')
>>>>>>> prody/master

    # def calcModes(self, n_modes=20, zeros=False, turbo=True):
    #     """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh`
    #     function to diagonalize the Hessian matrix. When Scipy is not found,
    #     :func:`numpy.linalg.eigh` is used.

    #     :arg n_modes: number of non-zero eigenvalues/vectors to calculate.
    #         If ``None`` is given, all modes will be calculated.
    #     :type n_modes: int or None, default is 20

    #     :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
    #     :type zeros: bool, default is ``False``

    #     :arg turbo: Use a memory intensive, but faster way to calculate modes.
    #     :type turbo: bool, default is ``True``
    #     """

    #     super(RTB, self).calcModes(n_modes, zeros, turbo)

    #     self._array = np.dot(self._project, self._array)


def test(pdb='1p38'):

    from prody import parsePDB
    from numpy import zeros


    pdb = parsePDB(pdb, subset='ca')
    bbenm = bbENM('1p38')
    bbenm.buildHessian(pdb)
    return bbenm

