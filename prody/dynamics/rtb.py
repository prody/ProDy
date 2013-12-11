# -*- coding: utf-8 -*-
"""This module defines a class and a function for rotating translating blocks
(RTB) calculations."""

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.proteins import parsePDB
from prody.utilities import importLA, checkCoords

from .anm import ANM

__all__ = ['RTB']


class RTB(ANM):

    """Class for Rotating Translating Blocks (RTB) analysis of proteins."""

    def buildHessian(self, coords, blocks, cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set.

        :arg coords: a coordinate set or an object with ``getCoords`` method
        :type coords: :class:`numpy.ndarray`

        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å, minimum is 4.0 Å
        :type cutoff: float

        :arg gamma: spring constant, default is 1.0
        :type gamma: float, :class:`Gamma`"""


        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                      coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with `getCoords` method')

        LOGGER.timeit('_rtb')
        natoms = coords.shape[0]
        if (natoms,) != blocks.shape:
            raise ValueError('blocks.shape must be (natoms,)')

        nblocks = len(set(blocks))
        nb6 = nblocks * 6

        coords = coords.T.copy()

        self._hessian = hessian = np.zeros((nb6, nb6), float)
        self._project = project = np.zeros((natoms * 3, nb6), float)

        from rtbtools import buildhessian
        buildhessian(coords, blocks, hessian, project,
                     natoms, nblocks, float(cutoff), float(gamma))

        LOGGER.report('Hessian was built in %.2fs.', label='_rtb')


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

        super(TRB, self).calcModes(n_modes, zeros, turbo)

        # do the projection here
