"""This module contains unit tests for :mod:`~prody.ensemble`."""

from prody.tests import TestCase

from numpy.testing import assert_equal

from prody import (calcOccupancies, trimPDBEnsemble, PDBEnsemble, 
                   parsePDB, buildPDBEnsemble, bestMatch, sameChid, 
                   sameChainPos)
from prody.tests.datafiles import *
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

class TestBuildPDBEnsemble(TestCase):
    
    def testResults(self):
        
        ags = parsePDB([DATA_FILES['3hsy']['path'], 
                        DATA_FILES['3o21']['path'], 
                        DATA_FILES['3p3w']['path']], 
                       subset='ca')
        
        ens1 = buildPDBEnsemble(ags, match_func=bestMatch)
        assert_equal(ens1.numConfs(), 5, 
            'buildPDBEnsemble with bestMatch on full structures did not extract all AMPAR dimers')
        
        ens2 = buildPDBEnsemble(ags, match_func=sameChid)
        assert_equal(ens2.numConfs(), 2, 
            'buildPDBEnsemble with sameChid on full structures did not just take AB dimers')
        
        ens3 = buildPDBEnsemble(ags, match_func=sameChainPos)
        assert_equal(ens3.numConfs(), 2, 
            'buildPDBEnsemble with sameChainPos on full structures did not just take AB dimers')        
   
        
        biomols = parsePDB([DATA_FILES["3hsy"]["path"],
                            DATA_FILES["3o21"]["path"],
                            DATA_FILES["3p3w"]["path"]], 
                        subset="ca",
                        biomol=True, extend_biomol=True)
        
        ens4 = buildPDBEnsemble(ags, match_func=bestMatch)
        assert_equal(ens4.numConfs(), 5, 
            'buildPDBEnsemble with bestMatch on biomols did not include all AMPAR dimers')    

        ens5 = buildPDBEnsemble(biomols, match_func=sameChid)
        assert_equal(ens5.numConfs(), 2, 
            'buildPDBEnsemble with sameChid on biomols did not just take AB dimer')
        
        ens6 = buildPDBEnsemble(biomols, match_func=bestMatch)
        assert_equal(ens6.numConfs(), 5, 
            'buildPDBEnsemble with sameChainPos on biomols did not include all AMPAR dimers')        
        
