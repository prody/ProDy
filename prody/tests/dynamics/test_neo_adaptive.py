"""This module contains unit tests for :mod:`~prody.dynamics.neo_adaptive`."""

import numpy as np
from numpy.testing import assert_array_equal

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

# adenylate kinase chain A, open (4AKE) and closed (1AKE), matching Calpha atoms
OPEN = parsePDB(pathDatafile('pdb4ake_fixed'), chain='A').select('calpha')
CLOSED = parsePDB(pathDatafile('pdb1ake'), chain='A').select('calpha')
N_STEPS = 4          # short runs: these tests check the guarantees of every step, not convergence


def bond_lengths(X):
    return np.linalg.norm(X[1:] - X[:-1], axis=1)


def closest_nonbonded(X, sep=3):
    i, j = np.triu_indices(len(X), k=sep)
    return np.linalg.norm(X[i] - X[j], axis=1).min()


def final_rmsd(ens):
    X = ens.getCoordsets()[-1]
    t = calcTransformation(X, CLOSED.getCoords())
    return calcRMSD(applyTransformation(t, X.copy()), CLOSED.getCoords())


class TestNeoAdaptiveANM(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ens = calcNeoAdaptiveANM(OPEN, CLOSED, N_STEPS, model='genanm')

    def testReturnsEnsembleWithProgress(self):
        self.assertIsInstance(self.ens, Ensemble)
        self.assertGreater(self.ens.numConfs(), 1)
        start = calcRMSD(*[applyTransformation(calcTransformation(OPEN.getCoords(), CLOSED.getCoords()),
                                               OPEN.getCoords().copy()), CLOSED.getCoords()])
        self.assertLess(final_rmsd(self.ens), start)

    def testEveryStepKeepsTheBackbone(self):
        """No accepted step stretches a virtual bond beyond the default guard (6.0 A)."""
        for X in self.ens.getCoordsets():
            self.assertLessEqual(bond_lengths(X).max(), 6.0 + 1e-6)

    def testEveryStepAvoidsClashes(self):
        """No accepted step brings non-bonded residues closer than the steric limit
        (3.0 A), except pairs already that close in an endpoint."""
        floor = min(3.0, closest_nonbonded(OPEN.getCoords()), closest_nonbonded(CLOSED.getCoords()))
        for X in self.ens.getCoordsets():
            self.assertGreaterEqual(closest_nonbonded(X), floor - 1e-6)

    def testEveryStepKeepsChainGeometry(self):
        """No frame shortens a virtual bond, or brings residues i and i+2 together, beyond the native-aware floors
        (defaults 3.4 A and 4.4 A, or the native value where that is smaller)."""
        from prody.dynamics.neo_adaptive import build_chain_floors
        bf, df = build_chain_floors([OPEN.getCoords(), CLOSED.getCoords()])
        for X in self.ens.getCoordsets():
            self.assertLessEqual((bf - bond_lengths(X)).max(), 0.05 + 1e-6)
            self.assertLessEqual((df - np.linalg.norm(X[2:] - X[:-2], axis=1)).max(), 0.05 + 1e-6)

    def testChainFloorsHonourNativeGeometry(self):
        """A bond shorter than the default in a native structure keeps its native length as the floor."""
        from prody.dynamics.neo_adaptive import build_chain_floors
        bf, df = build_chain_floors([OPEN.getCoords(), CLOSED.getCoords()], bond_lo=3.4, d13=4.4)
        native = np.minimum(bond_lengths(OPEN.getCoords()), bond_lengths(CLOSED.getCoords()))
        np.testing.assert_allclose(bf, np.minimum(3.4, native))
        self.assertTrue((bf <= 3.4).all() and (df <= 4.4).all())
        off = build_chain_floors([OPEN.getCoords()], bond_lo=0, d13=0)
        self.assertEqual(off, (None, None))

    def testEveryStepKeepsBondsBelowCeiling(self):
        """No frame stretches a virtual bond beyond its native-aware ceiling (4.2 A, or the native length where
        longer) by more than the 0.05 A tolerance."""
        from prody.dynamics.neo_adaptive import build_bond_ceiling
        ceiling = build_bond_ceiling([OPEN.getCoords(), CLOSED.getCoords()], 4.2)
        for X in self.ens.getCoordsets():
            self.assertLessEqual((bond_lengths(X) - ceiling).max(), 0.05 + 1e-6)

    def testChainGeometryCheck(self):
        """The chain test flags a bond above its ceiling and a bond below its floor, and passes the native chain."""
        from prody.dynamics.neo_adaptive import build_chain_floors, build_bond_ceiling, checkChainGeometry
        X = OPEN.getCoords().copy()
        bf, df = build_chain_floors([X]); ceiling = build_bond_ceiling([X], 4.2)
        self.assertFalse(checkChainGeometry(X, bf, df, bond_ceiling=ceiling)[0])
        Y = X.copy(); Y[10:] += 1.2 * (Y[10] - Y[9]) / np.linalg.norm(Y[10] - Y[9])     # stretch bond 9-10
        self.assertTrue(checkChainGeometry(Y, bf, df, bond_ceiling=ceiling)[0])
        Z = X.copy(); Z[10:] -= 0.8 * (Z[10] - Z[9]) / np.linalg.norm(Z[10] - Z[9])     # compress bond 9-10
        self.assertTrue(checkChainGeometry(Z, bf, None)[0])
        self.assertIsNone(build_bond_ceiling([X], 0))

    def testExemptPairsKeepNativeDistance(self):
        """A pair exempt from the steric test (already close in a native structure) may not come closer than its
        native distance: the chain test flags it, and the floors come from the closer native."""
        from prody.dynamics.neo_adaptive import build_clash_exempt, build_exempt_floors, checkChainGeometry
        X = OPEN.getCoords().copy(); Y = CLOSED.getCoords()
        ex = build_clash_exempt([X, Y], 4.5, 3); self.assertTrue(len(ex) > 0)
        P, fl = build_exempt_floors([X, Y], ex)
        native = np.minimum(np.linalg.norm(X[P[:, 0]] - X[P[:, 1]], axis=1), np.linalg.norm(Y[P[:, 0]] - Y[P[:, 1]], axis=1))
        np.testing.assert_allclose(fl, native)
        self.assertFalse(checkChainGeometry(X, exempt_floor=(P, fl))[0])
        i, j = P[0]; Z = X.copy(); Z[j] += 0.8 * (X[i] - X[j]) / np.linalg.norm(X[i] - X[j])    # push one exempt pair together
        self.assertTrue(checkChainGeometry(Z, exempt_floor=(P, fl))[0])
        self.assertIsNone(build_exempt_floors([X], set()))

    def testDeterministic(self):
        again = calcNeoAdaptiveANM(OPEN, CLOSED, N_STEPS, model='genanm')
        assert_array_equal(again.getCoordsets(), self.ens.getCoordsets())

    def testComponentsCanBeSwitchedOff(self):
        """Every NeoAdaptiveANM component can be turned off with one keyword."""
        ens = calcNeoAdaptiveANM(OPEN, CLOSED, N_STEPS, model='anm', soft_filter=0, clash_min=0,
                                 geom_bond_hi=0, geom_clash=0, contact_on=0)
        self.assertGreater(ens.numConfs(), 1)

    def testInvalidInput(self):
        self.assertRaises(TypeError, calcNeoAdaptiveANM, OPEN, CLOSED, 2.5)
        self.assertRaises(ValueError, calcNeoAdaptiveANM, OPEN, CLOSED, N_STEPS, mode=7)


if __name__ == '__main__':
    unittest.main()
