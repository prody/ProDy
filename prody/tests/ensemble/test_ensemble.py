"""This module contains unit tests for :mod:`~prody.ensemble`."""

import os.path
from prody.tests import TestCase

import numpy as np
from numpy import arange
from numpy.testing import assert_equal, assert_allclose

from . import ATOMS, COORDS, ENSEMBLE, ENSEMBLEW
from . import ENSEMBLE_RMSD, ENSEMBLE_SUPERPOSE
from . import ATOL, RTOL

class TestEnsemble(TestCase):

    def testGetCoordinates(self):
        """Test correctness of reference coordinates."""

        assert_equal(ENSEMBLE.getCoords(), COORDS,
                     'failed to get correct coordinates')


    def testGetCoordsets(self):
        """Test correctness of all coordinates."""

        assert_equal(ATOMS.getCoordsets(), ENSEMBLE.getCoordsets(),
                     'failed to add coordinate sets for Ensemble')


    def testGetWeights(self):
        """Test correctness of all weights."""

        assert_equal(ENSEMBLEW.getWeights().ndim, 2,
                     'failed to get correct weights ndim')
        assert_equal(ENSEMBLEW.getWeights().shape,
                     (ENSEMBLEW.numAtoms(), 1),
                    'failed to get correct weights shape')
        assert_equal(ENSEMBLEW.getWeights(),
                     np.ones((ENSEMBLEW.numAtoms(), 1), float),
                     'failed to get expected weights')

    def testSlicingCopy(self):
        """Test making a copy by slicing operation."""

        SLICE = ENSEMBLE[:]
        assert_equal(SLICE.getCoords(), ENSEMBLE.getCoords(),
                     'slicing copy failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), ENSEMBLE.getCoordsets(),
                     'slicing copy failed to add coordinate sets')

    def testSlicing(self):
        """Test slicing operation."""

        SLICE = ENSEMBLE[:2]
        assert_equal(SLICE.getCoords(), ENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), ENSEMBLE.getCoordsets([0,1]),
                     'slicing failed to add coordinate sets')

    def testSlicingList(self):
        """Test slicing with a list."""

        SLICE = ENSEMBLE[[0,2]]
        assert_equal(SLICE.getCoords(), ENSEMBLE.getCoords(),
                     'slicing failed to set reference coordinates')
        assert_equal(SLICE.getCoordsets(), ENSEMBLE.getCoordsets([0,2]),
                     'slicing failed to add coordinate sets for Ensemble')


    def testSlicingWeights(self):
        """Test slicing operation with weights."""

        SLICE = ENSEMBLEW[:2]
        assert_equal(SLICE.getWeights(), ENSEMBLEW.getWeights(),
                     'slicing failed to set weights')

    def testIterCoordsets(self):
        """Test coordinate iteration."""

        for i, xyz in enumerate(ENSEMBLE.iterCoordsets()):
            assert_equal(xyz, ATOMS.getCoordsets(i),
                         'failed yield correct coordinates')

    def testGetNumAtoms(self):

        self.assertEqual(ENSEMBLE.numAtoms(), ATOMS.numAtoms(),
                         'failed to get correct number of atoms')

    def testGetNumCsets(self):

        self.assertEqual(ENSEMBLE.numCoordsets(), ATOMS.numCoordsets(),
                         'failed to get correct number of coordinate sets')

    def testGetRMSDs(self):

        assert_allclose(ENSEMBLE.getRMSDs(), ENSEMBLE_RMSD,
                        rtol=0, atol=1e-3,
                        err_msg='failed to calculate RMSDs sets')

    def testSuperpose(self):

        ensemble = ENSEMBLE[:]
        ensemble.superpose()
        assert_allclose(ensemble.getRMSDs(), ENSEMBLE_SUPERPOSE,
                        rtol=0, atol=1e-3,
                        err_msg='failed to superpose coordinate sets')

    def testGetRMSDsWeights(self):

        assert_allclose(ENSEMBLEW.getRMSDs(), ENSEMBLE_RMSD,
                        rtol=0, atol=1e-3,
                        err_msg='failed to calculate RMSDs')

    def testSuperposeWeights(self):

        ensemble = ENSEMBLEW[:]
        ensemble.superpose()
        assert_allclose(ensemble.getRMSDs(), ENSEMBLE_SUPERPOSE,
                        rtol=0, atol=1e-3,
                        err_msg='failed to superpose coordinate sets')

    def testDelCoordsetMiddle(self):

        ensemble = ENSEMBLE[:]
        ensemble.delCoordset(1)
        assert_equal(ensemble.getCoordsets(), ATOMS.getCoordsets([0,2]),
                     'failed to delete middle coordinate set for Ensemble')

    def testDelCoordsetAll(self):

        ensemble = ENSEMBLE[:]
        ensemble.delCoordset(arange(len(ENSEMBLE)))
        self.assertIsNone(ensemble.getCoordsets(),
                        'failed to delete all coordinate sets for Ensemble')
        assert_equal(ensemble.getCoords(), COORDS,
                     'failed when deleting all coordinate sets')


    def testConcatenation(self):
        """Test concatenation of ensembles without weights."""

        ensemble = ENSEMBLE + ENSEMBLE
        assert_equal(ensemble.getCoordsets(arange(3)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(arange(3,6)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')

    def testConcatenationWeights(self):
        """Test concatenation of ensembles with weights."""

        ensemble = ENSEMBLEW + ENSEMBLEW
        assert_equal(ensemble.getCoordsets(arange(3)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(arange(3,6)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')
        assert_equal(ensemble.getWeights(), ENSEMBLEW.getWeights(),
                     'concatenation failed')

    def testConcatenationNoweightsWeights(self):
        """Test concatenation of ensembles without and with weights."""

        ensemble = ENSEMBLE + ENSEMBLEW
        assert_equal(ensemble.getCoordsets(arange(3)), ATOMS.getCoordsets(),
                     'concatenation failed')
        assert_equal(ensemble.getCoordsets(arange(3,6)), ATOMS.getCoordsets(),
                    'concatenation failed')
        assert_equal(ensemble.getCoords(), COORDS,
                     'concatenation failed')
        self.assertIsNone(ensemble.getWeights(), 'concatenation failed')

    def testConcatenationWeightsNoweights(self):
        """Test concatenation of ensembles with and without weights."""

        ensemble = ENSEMBLEW + ENSEMBLE
        assert_equal(ensemble.getCoordsets(arange(3)), ATOMS.getCoordsets(),
                     'failed at concatenation for Ensemble')
        assert_equal(ensemble.getCoordsets(arange(3,6)), ATOMS.getCoordsets(),
                     'failed at concatenation for Ensemble')
        assert_equal(ensemble.getCoords(), COORDS,
                     'failed at concatenation for Ensemble')
        assert_equal(ensemble.getWeights(), ENSEMBLEW.getWeights(),
                     'failed at concatenation for Ensemble')
