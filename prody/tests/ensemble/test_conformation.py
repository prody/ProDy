"""This module contains unit tests for :mod:`~prody.ensemble`."""

from prody.tests import TestCase

from numpy import ones
from numpy.testing import assert_equal

from . import ATOMS, CONF, COORDS, CONFW, ENSEMBLE
from . import PDBCONF, PDBENSEMBLE, WEIGHTS, WEIGHTS_BOOL


class TestConformation(TestCase):


    def testCoordinates(self):

        assert_equal(CONF.getCoords(), COORDS,
                     'failed to get coordinates for conformation')

    def testWeights(self):

        weights = CONFW.getWeights()
        self.assertEqual(weights.ndim, 2,
                        'wrong ndim for weights of Conformation')
        self.assertTupleEqual(weights.shape,
                              (CONFW.numAtoms(), 1),
                              'wrong shape for weights of Conformation')
        assert_equal(weights, ones((ATOMS.numAtoms(), 1), float),
                     'failed to set weights for Conformation')

    def testCoordinatesForAll(self):

        for i, conf in enumerate(ENSEMBLE):
            assert_equal(conf.getCoords(), ATOMS.getCoordsets(i),
                         'failed to get coordinates for Conformation')

    def testGetIndex(self):
        """Test get index function."""

        for i, conf in enumerate(ENSEMBLE):
            self.assertEqual(conf.getIndex(), i,
                             'failed to get correct index')

    def testGetNumAtoms(self):
        """Test get index function."""

        for i, conf in enumerate(ENSEMBLE):
            self.assertEqual(conf.numAtoms(), ATOMS.numAtoms(),
                             'failed to get correct number of atoms')


class TestPDBConformation(TestCase):


    def testCoordinates(self):

        assert_equal(PDBCONF.getCoords(), COORDS,
                     'failed to get correct coordinates')

    def testWeights(self):

        weights = PDBCONF.getWeights()
        self.assertEqual(weights.ndim, 2,
                        'wrong ndim for weights')
        self.assertTupleEqual(weights.shape, (PDBCONF.numAtoms(), 1),
                              'wrong shape for weights')
        assert_equal(PDBCONF.getWeights(), WEIGHTS[0],
                     'failed to set weights')

    def testWeightsForAll(self):

        for i, conf in enumerate(PDBENSEMBLE):
            assert_equal(conf.getWeights(), WEIGHTS[i],
                         'failed to set weights')


    def testCoordinatesForAll(self):

        for i, conf in enumerate(PDBENSEMBLE):
            assert_equal(conf.getCoords()[WEIGHTS_BOOL[i]],
                         ATOMS.getCoordsets(i)[WEIGHTS_BOOL[i]],
                         'failed to get coordinates')

    def testGetIndex(self):

        for i, conf in enumerate(PDBENSEMBLE):
            self.assertEqual(conf.getIndex(), i,
                             'failed to get correct index')

    def testGetNumAtoms(self):

        for i, conf in enumerate(PDBENSEMBLE):
            self.assertEqual(conf.numAtoms(), ATOMS.numAtoms(),
                             'failed to get correct number of atoms')
