# -*- coding: utf-8 -*-
"""This module defines a class and a function for rotating translating blocks
(RTB) calculations."""

import numpy as np

from prody import LOGGER
from prody.utilities import checkCoords

from .anm import ANMBase

__all__ = ['RTB']

class Increment(object):

    def __init__(self, s=0):

        self._i = s

    def __call__(self, i=1):

        self._i += i
        return self._i


class RTB(ANMBase):

    """Class for Rotations and Translations of Blocks (RTB) method ([FT00]_).

    .. [FT00] Tama F, Gadea FJ, Marques O, Sanejouand YH. Building-block
       approach for determining low-frequency normal modes of macromolecules.
       *Proteins* **2000** 41:1-7.

    """

    def __init__(self, name='Unknown'):

        super(RTB, self).__init__(name)
        self._project = None


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

        super(RTB, self).buildHessian(coords, cutoff=cutoff, gamma=gamma, **kwargs)

        self.calcProjection(coords, blocks, **kwargs)


    def calcProjection(self, coords, blocks, **kwargs):
        natoms = self._n_atoms

        if natoms != len(blocks):
            raise ValueError('len(blocks) must match number of atoms')

        LOGGER.timeit('_rtb')
        from collections import defaultdict
        i = Increment()
        d = defaultdict(i)
        blocks = np.array([d[b] for b in blocks], dtype='int32')

        try:
            from collections import Counter
        except ImportError:
            counter = defaultdict(int)
            for b in blocks:
                counter[b] += 1
        else:
            counter = Counter(blocks)

        nblocks = len(counter)
        maxsize = 1
        nones = 0
        while counter:
            _, size = counter.popitem()
            if size == 1:
                nones += 1
            if size > maxsize:
                maxsize = size
        LOGGER.info('System has {0} blocks, the largest with {1} of {2} units.'
                    .format(nblocks, maxsize, natoms))
        nb6 = nblocks * 6 - nones * 3

        coords = coords.T.astype(float, order='C')

        hessian = self._hessian
        self._project = project = np.zeros((natoms * 3, nb6), float)

        from .rtbtools import calc_projection

        calc_projection(coords, blocks, project, natoms, nblocks, nb6, maxsize)

        self._hessian = project.T.dot(hessian).dot(project)
        self._dof = self._hessian.shape[0]
        LOGGER.report('Block Hessian and projection matrix were calculated in %.2fs.', label='_rtb')


    def getProjection(self):
        """Returns a copy of the projection matrix."""

        if self._project is not None:
            return self._project.copy()

    def _getProjection(self):

        return self._project

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
        if n_modes is None:
            n_modes = self._dof
        super(RTB, self).calcModes(n_modes, zeros, turbo)
        self._array = np.dot(self._project, self._array)
