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
        
        # anm hessian calculation 
        cutoff, gamma, gamma_func = checkENMParameters(cutoff, gamma)
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
                g = gamma_func(dist2, i, j)
                res_j3 = j*3
                res_j33 = res_j3+3
                super_element = np.outer(i2j, i2j) * (- g / dist2)
                hessian[res_i3:res_i33, res_j3:res_j33] = super_element
                hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                hessian[res_i3:res_i33, res_i3:res_i33] = \
                    hessian[res_i3:res_i33, res_i3:res_i33] - super_element
                hessian[res_j3:res_j33, res_j3:res_j33] = \
                    hessian[res_j3:res_j33, res_j3:res_j33] - super_element

        # hessian updates
        from .bbenmtools import buildhessian

        buildhessian(coords, hessian, natoms, 
                     float(cutoff), float(gamma),)

        LOGGER.report('Hessian was built in %.2fs.', label='_bbenm')

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

         super(bbENM, self).calcModes(n_modes, zeros, turbo)

    #     self._array = np.dot(self._project, self._array)


def test(pdb='2ci2'):

    from prody import parsePDB
    from numpy import zeros


    pdb = parsePDB(pdb, subset='ca')
    bbenm = bbENM('2ci2')
    bbenm.buildHessian(pdb, cutoff=7.)
    bbenm.calcModes(n_modes = None)
    return bbenm

