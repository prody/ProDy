"""This module contains unit tests for :mod:`~prody.ensemble`."""

from prody.tests import TestCase

from numpy.testing import assert_equal

from prody import calcOccupancies, trimPDBEnsemble, PDBEnsemble
from . import PDBENSEMBLE, WEIGHTS, ENSEMBLE, ATOMS, PDBENSEMBLEA


class TestCalcOccupancies(TestCase):

    def testResults(self):
        assert_equal(calcOccupancies(PDBENSEMBLE), WEIGHTS.sum(0).flatten(),
                     'calcOccupancies failed')

    def testInvalidType(self):

        self.assertRaises(TypeError, calcOccupancies, ENSEMBLE)

    def testWeightsNone(self):

        self.assertRaises(ValueError, calcOccupancies, PDBEnsemble())

class TestTrimPDBEnsemble(TestCase):

    def testResults(self):
        ensemble = PDBENSEMBLEA[:]
        weights = ensemble.getWeights()
        weights[:, 0, :] = 0.
        ensemble.setWeights(weights)
        ensemble_trimed = trimPDBEnsemble(ensemble, occupancy=0.9)
        assert calcOccupancies(ensemble_trimed).min() > 0.9, \
                     'soft trimPDBEnsemble failed'
        
        ensemble_trimed2 = trimPDBEnsemble(ensemble, occupancy=0.9, hard=True)
        assert calcOccupancies(ensemble_trimed2).min() > 0.9, \
                     'hard trimPDBEnsemble failed'

        atoms = ensemble_trimed.getAtoms()
        atoms2 = ATOMS[ensemble_trimed._indices]
        assert_equal(atoms.getIndices(), atoms2.getIndices(), 
                    'soft trimPDBEnsemble returns a wrong result')

        msa1 = ensemble_trimed.getMSA()
        msa2 = ensemble_trimed2.getMSA()
        assert_equal(msa1.getArray(), msa2.getArray(), 
                    'soft trimPDBEnsemble returns a wrong result')

