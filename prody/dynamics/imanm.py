# -*- coding: utf-8 -*-
"""This module defines a class and a function for explicit membrane ANM calculations."""

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.utilities import importLA, checkCoords, copy
from numpy import sqrt, zeros, ones, array, ceil, dot

from .anm import ANMBase
from .rtb import RTB
from .editing import reduceModel

LA = importLA()
inv = LA.inv
norm = LA.norm

__all__ = ['imANM']

class imANM(RTB):

    """Class for implicit ANM (imANM) method ([FT00]_).

    .. [TL12] Lezon TR, Bahar I, Constraints Imposed by the Membrane
       Selectively Guide the Alternating Access Dynamics of the Glutamate
       Transporter GltPh

    """

    def __init__(self, name='Unknown'):

        super(imANM, self).__init__(name)

    def buildHessian(self, coords, blocks, cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg blocks: a list or array of block identifiers
        :type blocks: list, :class:`numpy.ndarray`

        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float

        :arg scale: scaling factor for force constant along Z-direction,
            default is 16.0
        :type scale: float

        :arg membrane_low: minimum z-coordinate at which membrane scaling
            is applied
            default is 1.0
        :type membrane_low: float

        :arg membrane_high: maximum z-coordinate at which membrane scaling
            is applied.  If membrane_high < membrane_low, scaling will be 
            applied to the entire structure
            default is -1.0
        :type membrane_high: float
        """

        scale = kwargs.pop('scale', 16.0)
        depth = kwargs.pop('depth', None)
        h = depth / 2 if depth is not None else None
            
        h = kwargs.pop('h', h)
        if h is not None:
            h = float(h)
            hu = h
            hl = -h
        else:
            hu = kwargs.pop('membrane_hi', 13.0)
            hu = kwargs.pop('high', hu)
            hu = float(hu)
            
            hl = kwargs.pop('membrane_lo', -13.0)
            hl = kwargs.pop('low', hl)
            hl = float(hl)

        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')

        ANMBase.buildHessian(self, coords, cutoff=cutoff, gamma=gamma, **kwargs)

        ## Scale horizontal spring constants ##
        natm = self._n_atoms
        H = self._hessian

        s = sqrt(sqrt(scale))
        S0 = array([[s*s, s*s, s],
                    [s*s, s*s, s],
                    [s  , s,   1]], dtype=float)

        super_element = lambda i, j: H[i*3:(i+1)*3, j*3:(j+1)*3]
        scaler = lambda coords: S0 if coords[2] < hu and coords[2] > hl else ones((3, 3), dtype=float)

        for i in range(natm):
            Si = scaler(coords[i])
            for j in range(i+1, natm):
                Sj = scaler(coords[j])
                S = Si * Sj
                H33 = super_element(i, j)

                if np.any(H33 != 0) and np.any(S != 1):
                    H0 = H33.copy()
                    H33 *= S

                    D33 = super_element(i, i)
                    D33 += H0 - H33

        self.calcProjection(coords, blocks, **kwargs)
    
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

        super(imANM, self).calcModes(n_modes, zeros, turbo)


def test(pdb='2nwl-mem.pdb', blk='2nwl.blk'):

    from prody import parsePDB
    from numpy import zeros, linalg

    pdb = parsePDB(pdb, subset='ca')
    pdb.setData('block', zeros(len(pdb), int))
    with open(blk) as inp:
        for line in inp:
            if line.startswith('BLOCK'):
                _, b, n1, c1, r1, n2, c2, r2 = line.split()
                sel = pdb.select('chain {} and resnum {} to {}'
                                 .format(c1, r1, r2))              
                if sel:
                    sel.setData('block', int(b))
    pdb.setBetas(pdb.getData('block'))
    coords = pdb.getCoords() 
    blocks = pdb.getBetas()
    from prody import writePDB
    writePDB('pdb2gb1_truncated.pdb', pdb)
    anm = imANM('2nwl')
    anm.buildHessian(coords, blocks, scale=64)
    #anm.calcModes()
    return anm
