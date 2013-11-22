"""This module contains unit tests for :mod:`~prody.ensemble`."""

from prody.tests import TestCase

from numpy import arange
from numpy.testing import assert_equal

from . import ATOMS, PDBENSEMBLE, COORDS, WEIGHTS_BOOL, ENSEMBLE, WEIGHTS

class TestPDBEnsemble(TestCase):

    def testGetCoordinates(self):

        assert_equal(PDBENSEMBLE.getCoords(), COORDS,
                     'failed to set reference coordinates for PDBEnsemble')

    def testGetCoordsets(self):

        assert_equal(PDBENSEMBLE.getCoordsets()[WEIGHTS_BOOL],
                     ATOMS.getCoordsets()[WEIGHTS_BOOL],
                     'failed to add coordinate sets for PDBEnsemble')

    def testGetWeights(self):

        self.assertEqual(PDBENSEMBLE.getWeights().ndim, 3,
                        'wrong ndim for weights of PDBEnsemble')
        self.assertTupleEqual(PDBENSEMBLE.getWeights().shape,
                              (PDBENSEMBLE.numCoordsets(),
                               PDBENSEMBLE.numAtoms(), 1),
                               'wrong shape for weights of PDBEnsemble')
        assert_equal(PDBENSEMBLE.getWeights(), WEIGHTS,
                     'failed to get correct weights')


    def testSlicingCopy(self):

        SLICE = PDBENSEMBLE[:]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE.getCoords(),
                     'slicing copy failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE.getCoordsets(),
                     'slicing copy failed to add coordinate sets')

    def testSlicing(self):

        SLICE = PDBENSEMBLE[:2]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE.getCoordsets([0,1]),
                     'slicing failed to add coordinate sets')

    def testSlicingList(self):

        SLICE = PDBENSEMBLE[[0,2]]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE.getCoordsets([0,2]),
                     'slicing failed to add coordinate sets')

    def testSlicingWeights(self):

        SLICE = PDBENSEMBLE[:2]
        assert_equal(SLICE.getWeights(), PDBENSEMBLE.getWeights()[:2],
                     'slicing failed to set weights')

    def testIterCoordsets(self):

        for i, xyz in enumerate(ENSEMBLE.iterCoordsets()):
            assert_equal(xyz[WEIGHTS_BOOL[i]],
                         ATOMS.getCoordsets(i)[WEIGHTS_BOOL[i]],
                         'failed iterate coordinate sets')

    def testGetNumAtoms(self):

        self.assertEqual(PDBENSEMBLE.numAtoms(), ATOMS.numAtoms(),
                         'failed to get correct number of atoms')

    def testGetNumCsets(self):

        self.assertEqual(PDBENSEMBLE.numCoordsets(),
                         ATOMS.numCoordsets(),
                         'failed to get correct number of coordinate sets')

    def testDelCoordsetMiddle(self):

        ensemble = PDBENSEMBLE[:]
        ensemble.delCoordset(1)
        assert_equal(ensemble.getCoordsets()[WEIGHTS_BOOL[[0,2]]],
                     ATOMS.getCoordsets([0,2])[WEIGHTS_BOOL[[0,2]]],
                    'failed to delete middle coordinate set')

    def testDelCoordsetAll(self):
        """Test consequences of deleting all coordinate sets."""

        ensemble = PDBENSEMBLE[:]
        ensemble.delCoordset(arange(len(PDBENSEMBLE)))
        self.assertIsNone(ensemble.getCoordsets(),
                        'failed to delete all coordinate sets')
        self.assertIsNone(ensemble.getWeights(), 'failed to delete weights '
                          'with all coordinate sets')
        assert_equal(ensemble.getCoords(), COORDS,
                     'failed to delete all coordinate sets')


    def testConcatenation(self):
        """Test concatenation of PDB ensembles."""

        ensemble = PDBENSEMBLE + PDBENSEMBLE
        assert_equal(ensemble.getCoordsets(arange(3)),
                     PDBENSEMBLE.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(arange(3,6)),
                     PDBENSEMBLE.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')
        assert_equal(ensemble.getWeights()[arange(3)],
                     PDBENSEMBLE.getWeights(),
                     'concatenation failed')
        assert_equal(ensemble.getWeights()[arange(3,6)],
                     PDBENSEMBLE.getWeights(),
                     'concatenation failed')
