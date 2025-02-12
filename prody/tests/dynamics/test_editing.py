"""This module contains unit tests for :mod:`~prody.KDTree` module."""

from numpy.testing import assert_array_equal, assert_equal

from numpy import concatenate, array

from prody.dynamics import calcGNM, calcANM
from prody.dynamics import (extendModel, extendMode, extendVector,
                            sliceModel, sliceMode, sliceVector,
                            interpolateModel)

from prody.tests import unittest
from prody.tests.datafiles import parseDatafile

from prody import _PY3K, LOGGER

LOGGER.verbosity = 'none'


if not _PY3K:
    range = xrange

ATOMS = parseDatafile('1ubi').protein.copy()
ANM, NODES = calcANM(ATOMS)
GNM = calcGNM(ATOMS)[0]

EXT1D = concatenate([[i] * len(res)
                     for i, res in enumerate(ATOMS.iterResidues())])

EXT3D = concatenate([list(range(i*3, (i+1)*3)) * len(res)
                     for i, res in enumerate(ATOMS.iterResidues())])

SLC1D = concatenate([[i] for i in ATOMS.ca.getIndices()])

SLC3D = concatenate([list(range(i*3, (i+1)*3))
                     for i in ATOMS.ca.getIndices()])

ANM_AA, NODES_AA = calcANM(ATOMS, selstr='all', cutoff=8)
GNM_AA = calcGNM(ATOMS, selstr='all', cutoff=4)[0]


class TestExtending(unittest.TestCase):

    def testModel3d(self):

        ext = extendModel(ANM, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), ANM._getArray()[EXT3D, :])

    def testModel1d(self):

        ext = extendModel(GNM, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), GNM._getArray()[EXT1D, :])

    def testMode3d(self):

        mode = ANM[0]
        ext = extendMode(mode, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), mode._getArray()[EXT3D] *
                                            mode.getVariance()**0.5)

    def testMode1d(self):

        mode = GNM[0]
        ext = extendMode(mode, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), mode._getArray()[EXT1D] *
                                            mode.getVariance()**0.5)

    def testVector3d(self):

        vector = ANM[0] * 1
        ext = extendVector(vector, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), vector._getArray()[EXT3D])

    def testVector1d(self):

        vector = GNM[0] * 1
        ext = extendVector(vector, NODES, ATOMS)[0]
        assert_array_equal(ext._getArray(), vector._getArray()[EXT1D])


class TestSlicing(unittest.TestCase):

    def testModel3d(self):

        slc = sliceModel(ANM_AA, NODES_AA, NODES)[0]
        assert_array_equal(slc._getArray(), ANM_AA._getArray()[SLC3D, :])

    def testModel1d(self):

        slc = sliceModel(GNM_AA, NODES_AA, NODES)[0]
        assert_array_equal(slc._getArray(), GNM_AA._getArray()[SLC1D, :])

    def testMode3d(self):

        mode = ANM_AA[0]
        slc = sliceMode(mode, NODES_AA, NODES)[0]
        assert_array_equal(slc._getArray(), mode._getArray()[SLC3D] *
                                            mode.getVariance()**0.5)

    def testMode1d(self):

        mode = GNM_AA[0]
        slc = sliceMode(mode, NODES_AA, NODES)[0]
        assert_array_equal(slc._getArray(), mode._getArray()[SLC1D] *
                                            mode.getVariance()**0.5)

    def testVector3d(self):

        vector = ANM_AA[0] * 1
        slc = sliceVector(vector, NODES_AA, NODES)[0]
        assert_array_equal(slc._getArray(), vector._getArray()[SLC3D])

    def testVector1d(self):

        vector = GNM_AA[0] * 1
        slc = sliceVector(vector, NODES_AA, NODES)[0]
        assert_array_equal(slc._getArray(), vector._getArray()[SLC1D])


class TestInterpolating(unittest.TestCase):

    def testModel3d(self):

        interp, atoms = interpolateModel(ANM, NODES, ATOMS)
        assert_equal(interp._getArray().shape, ANM._getArray()[EXT3D, :].shape)
        assert_equal(atoms.numAtoms(), NODES_AA.numAtoms())
