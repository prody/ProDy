import unittest
import numpy as np
from prody import *

# TODO: Import your function here
# from your_module_name import getGlobalHinges

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


if __name__ == '__main__':
    # Verbosity=2 shows individual test status (OK/FAIL)
    unittest.main(verbosity=2)
