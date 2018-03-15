"""This module contains unit tests for :mod:`~prody.ensemble`."""

from prody.tests import TestCase

from numpy import arange
from numpy.testing import assert_equal

from . import ATOMS, PDBENSEMBLE, PDBENSEMBLE2, COORDS, WEIGHTS_BOOL, ENSEMBLE, WEIGHTS

class TestPDBEnsemble(TestCase):

    def testGetCoordinates(self):

        assert_equal(PDBENSEMBLE.getCoords(), COORDS,
                     'failed to set reference coordinates for PDBEnsemble')
        assert_equal(PDBENSEMBLE2.getCoords(), COORDS,
                     'failed to set reference coordinates for PDBEnsemble')

    def testGetCoordsets(self):

        assert_equal(PDBENSEMBLE.getCoordsets()[WEIGHTS_BOOL],
                     ATOMS.getCoordsets()[WEIGHTS_BOOL],
                     'failed to add coordinate sets for PDBEnsemble')
        assert_equal(PDBENSEMBLE2.getCoordsets()[WEIGHTS_BOOL],
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
        
        self.assertEqual(PDBENSEMBLE2.getWeights().ndim, 3,
                        'wrong ndim for weights of PDBEnsemble')
        self.assertTupleEqual(PDBENSEMBLE2.getWeights().shape,
                              (PDBENSEMBLE2.numCoordsets(),
                               PDBENSEMBLE2.numAtoms(), 1),
                               'wrong shape for weights of PDBEnsemble')
        assert_equal(PDBENSEMBLE2.getWeights(), WEIGHTS,
                     'failed to get correct weights')


    def testSlicingCopy(self):

        SLICE = PDBENSEMBLE2[:]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE2.getCoords(),
                     'slicing copy failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE2.getCoordsets(),
                     'slicing copy failed to add coordinate sets')

    def testSlicing(self):

        SLICE = PDBENSEMBLE2[:2]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE2.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE2.getCoordsets([0,1]),
                     'slicing failed to add coordinate sets')

    def testSlicingList(self):

        SLICE = PDBENSEMBLE2[[0,2]]
        assert_equal(SLICE.getCoords(), PDBENSEMBLE2.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), PDBENSEMBLE2.getCoordsets([0,2]),
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

        ensemble = PDBENSEMBLE2[:]
        ensemble.delCoordset(1)
        assert_equal(ensemble.getCoordsets()[WEIGHTS_BOOL[[0,2]]],
                     ATOMS.getCoordsets([0,2])[WEIGHTS_BOOL[[0,2]]],
                    'failed to delete middle coordinate set')

    def testDelCoordsetAll(self):
        """Test consequences of deleting all coordinate sets."""

        ensemble = PDBENSEMBLE2[:]
        ensemble.delCoordset(arange(len(PDBENSEMBLE2)))
        self.assertIsNone(ensemble.getCoordsets(),
                        'failed to delete all coordinate sets')
        self.assertIsNone(ensemble.getWeights(), 'failed to delete weights '
                          'with all coordinate sets')
        assert_equal(ensemble.getCoords(), COORDS,
                     'failed to delete all coordinate sets')


    def testConcatenation(self):
        """Test concatenation of PDB ensembles."""

        ensemble = PDBENSEMBLE + PDBENSEMBLE2
        assert_equal(ensemble.getCoordsets(arange(3)),
                     PDBENSEMBLE.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(arange(3,6)),
                     PDBENSEMBLE2.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')
        assert_equal(ensemble.getWeights()[arange(3)],
                     PDBENSEMBLE.getWeights(),
                     'concatenation failed')
        assert_equal(ensemble.getWeights()[arange(3,6)],
                     PDBENSEMBLE2.getWeights(),
                     'concatenation failed')

        msa = PDBENSEMBLE2.getMSA()
        msa2 = ensemble.getMSA(arange(3,6))

        assert_equal(msa, msa2, 'associated MSA concatenation failed')


