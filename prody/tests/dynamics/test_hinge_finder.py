import unittest
from unittest.mock import patch

import numpy as np
from prody import *

# TODO: Import your function here
# from your_module_name import getGlobalHinges


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


# --- THE UNIT TEST SUITE ---
class TestHingesWithRealGNM(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Runs once before all tests. 
        Downloads 1AKE and calculates modes to save time.
        """
        print("\n[Setup] Downloading 1AKE and building GNM...")
        
        # 1. Parse PDB (Downloads if not found)
        cls.pdb = parsePDB("1AKE", subset='ca')
        
        # 2. Build GNM
        cls.gnm = GNM('1AKE')
        cls.gnm.buildKirchhoff(cls.pdb)
        cls.gnm.calcModes(n_modes='all')
        
        print("[Setup] GNM built successfully.")

    def test_output_structure(self):
        """Verify the output is a list of lists of integers."""
        # Use explicit n_modes to keep it simple
        hinges = getGlobalHinges(self.gnm, n_modes=2)
        
        self.assertIsInstance(hinges, list, "Output must be a list")
        self.assertEqual(len(hinges), 2, "Must return results for exactly 2 modes")
        self.assertIsInstance(hinges[0], list, "Each mode result must be a list")
        
        # Check if contents are integers (if any hinges found)
        if len(hinges[0]) > 0:
            self.assertIsInstance(hinges[0][0], (int, np.integer), "Hinge indices must be integers")

    def test_auto_mode_selection(self):
        """Verify that passing n_modes=None triggers auto-selection."""
        # The function defaults to 33% cumulative variance
        hinges = getGlobalHinges(self.gnm, n_modes=None)
        
        num_modes_selected = len(hinges)
        
        # Calculate expected number manually
        fv = calcFractVariance(self.gnm)
        cum_var = np.cumsum(fv)
        expected_modes = np.argmax(cum_var >= 0.33) + 1
        
        print(f"\n[Test Auto] Auto-selected {num_modes_selected} modes.")
        self.assertEqual(num_modes_selected, expected_modes, 
                         f"Should select {expected_modes} modes for 33% variance")

    def test_trim_functionality(self):
        """Verify that 'trim=True' removes hinges at the very ends."""
        # 1AKE has ~214 residues. 
        # trim = length // 20 = 214 // 20 = 10 residues cut from each end.
        
        # First, run WITHOUT trim
        hinges_untrimmed = getGlobalHinges(self.gnm, n_modes=1, trim=False)[0]
        
        # Then, run WITH trim
        hinges_trimmed = getGlobalHinges(self.gnm, n_modes=1, trim=True)[0]
        
        # Calculate the forbidden zones
        n_atoms = self.gnm.numAtoms()
        trim_zone = n_atoms // 20
        
        # Ensure no hinges in trimmed zones
        for h in hinges_trimmed:
            self.assertTrue(trim_zone <= h < (n_atoms - trim_zone),
                            f"Hinge at {h} should have been trimmed (Limit: {trim_zone} to {n_atoms-trim_zone})")
            
        # Verify we didn't lose Valid hinges from the middle
        # (This assumes the un-trimmed version had middle hinges, which 1AKE Mode 1 usually does)
        self.assertGreater(len(hinges_trimmed), 0, "1AKE Mode 1 should have valid central hinges")


class TestGetHingesEmptyRegions(unittest.TestCase):
    """getHinges must not crash when there is no sign change.

    ``regs`` ends up empty and the old code indexed ``regs[0]``.
    """

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
    """Hinges from all chains but the last were silently dropped.

    The accumulation block sat outside the per-chain loop.
    """

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
    """``vecs[fst:lst]`` should be ``vecs[fst:lst+1]``.

    ``chains`` stores inclusive (fst, lst) pairs, so the exclusive slice
    dropped the last atom of every chain from consideration.
    """

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


if __name__ == '__main__':
    # Verbosity=2 shows individual test status (OK/FAIL)
    unittest.main(verbosity=2)
