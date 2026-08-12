"""Regression tests for three bugs found in prody.dynamics.analysis:

1. getGlobalHinges silently drops hinges from every chain except the last
   because the accumulation block is mis-indented outside the per-chain
   loop, and also slices vecs with an exclusive upper bound (vecs[fst:lst])
   even though `chains` stores inclusive (fst, lst) pairs, dropping the
   last atom of every chain from consideration.
2. getHinges raises IndexError when an eigenvector segment has no sign
   change at all (regs ends up empty and `regs[0]` blows up).
3. calcCollectivity normalizes by dividing by sqrt(sum) instead of sum when
   masses are supplied, which breaks the probability-distribution
   requirement of the Bruschweiler 1995 (eq. 5) collectivity definition.
"""

import unittest
from unittest.mock import patch

import numpy as np
from prody import *
from prody.dynamics.analysis import getHinges, getGlobalHinges, calcCollectivity


def _make_two_chain_atoms(n1, n2):
    """Tiny synthetic 2-chain AtomGroup (CA-only) for hinge tests."""
    n = n1 + n2
    ag = AtomGroup('synthetic')
    ag.setCoords(np.zeros((n, 3)))
    ag.setChids(['A'] * n1 + ['B'] * n2)
    ag.setResnums(list(range(1, n1 + 1)) + list(range(1, n2 + 1)))
    ag.setNames(['CA'] * n)
    ag.setResnames(['ALA'] * n)
    return ag


class TestGetHingesEmptyRegions(unittest.TestCase):
    """Bug 2: getHinges must not crash when there is no sign change."""

    def test_no_sign_change_returns_empty_list(self):
        v = np.array([1.0, 1.0, 1.0, 1.0])
        # Should return [] (no hinges), not raise IndexError.
        hinges = getHinges(v)
        self.assertEqual(hinges, [])

    def test_normal_case_still_works(self):
        # Sanity check the fix doesn't break the normal path.
        v = np.array([-1, -1, -1, 1, 1, 1.0])
        self.assertEqual(getHinges(v), [2])


class TestGetGlobalHingesChainAccumulation(unittest.TestCase):
    """Bug 1a: hinges from all chains but the last are silently dropped."""

    def test_all_chains_contribute_hinges(self):
        ag = _make_two_chain_atoms(6, 5)

        # Chain A (global idx 0-5): crossing well inside -> local hinge 2
        v_a = np.array([-1, -1, -1, 1, 1, 1.0])
        # Chain B (global idx 6-10): crossing well inside -> local hinge 1
        v_b = np.array([1, 1, -1, -1, -1.0])
        vecs = np.concatenate([v_a, v_b])[:, np.newaxis]

        with patch('prody.dynamics.analysis._getModeVecs', return_value=vecs):
            hinges = getGlobalHinges(object(), atoms=ag, threshold=15,
                                      space=None, trim=False)

        # Expect a hinge from BOTH chains: global 2 (chain A) and 7 (chain B).
        self.assertEqual(hinges, [[2, 7]])


class TestGetGlobalHingesOffByOne(unittest.TestCase):
    """Bug 1b: vecs[fst:lst] should be vecs[fst:lst+1] (inclusive chains)."""

    def test_hinge_at_last_residue_of_chain_is_detected(self):
        # Single chain so the accumulation-indentation bug can't interfere
        # (loop runs exactly once) -- isolates the slicing off-by-one.
        ag = _make_two_chain_atoms(6, 0)  # single 6-atom chain 'A'
        # crossings at local (1,2) and (4,5); the second one needs the very
        # last atom of the chain (global index 5) to be included in the slice.
        v = np.array([-1, -1, 1, 1, 1, -1.0])
        vecs = v[:, np.newaxis]

        with patch('prody.dynamics.analysis._getModeVecs', return_value=vecs):
            hinges = getGlobalHinges(object(), atoms=ag, threshold=15,
                                      space=None, trim=False)

        self.assertEqual(hinges, [[1, 4]])


class TestCalcCollectivityMasses(unittest.TestCase):
    """Bug 3: masses path must divide by sum(u2in), not sqrt(sum(u2in))."""

    def test_masses_path_matches_bruschweiler_definition(self):
        v = np.array([0.6, 0.8])
        masses = np.array([1.0, 4.0])

        u2 = v ** 2 / masses
        p = u2 / u2.sum()
        self.assertAlmostEqual(p.sum(), 1.0, places=12)
        expected = np.exp(-(p * np.log(p)).sum()) / len(v)

        actual = calcCollectivity(v, masses=masses, is3d=False)
        self.assertAlmostEqual(actual, expected, places=8)

    def test_default_path_without_masses_is_unaffected(self):
        # Unit-norm vector: sum(u2in) == sqrt(sum(u2in)) == 1, so this path
        # is identical before and after the fix. Guards against regressions.
        v = np.array([0.6, 0.8])
        u2 = v ** 2
        p = u2 / u2.sum()
        expected = np.exp(-(p * np.log(p)).sum()) / len(v)

        actual = calcCollectivity(v, is3d=False)
        self.assertAlmostEqual(actual, expected, places=8)


if __name__ == '__main__':
    unittest.main(verbosity=2)
