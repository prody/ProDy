"""This module contains unit tests for :mod:`~prody.ensemble`."""

from prody.tests import TestCase

from numpy.testing import assert_equal

from prody import calcOccupancies, PDBEnsemble
from . import PDBENSEMBLE, WEIGHTS, ENSEMBLE


class TestCalcOccupancies(TestCase):

    def testResults(self):
        assert_equal(calcOccupancies(PDBENSEMBLE), WEIGHTS.sum(0).flatten(),
                     'calcOccupancies failed')

    def testInvalidType(self):

        self.assertRaises(TypeError, calcOccupancies, ENSEMBLE)

    def testWeightsNone(self):

        self.assertRaises(ValueError, calcOccupancies, PDBEnsemble())

