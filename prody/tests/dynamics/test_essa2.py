import unittest
import os
import shutil
import tempfile
import numpy as np
from prody import parsePDB, AtomGroup
from prody import ESSA2 

class TestESSA2(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Create a minimal dummy PDB file for testing."""
        cls.test_dir = tempfile.mkdtemp()
        cls.pdb_path = os.path.join(cls.test_dir, "test.pdb")
        # Minimal PDB content: 3 Alanine residues (CA, N, C, O, CB)
        pdb_content = (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.009   1.427   0.000  1.00  0.00           C\n"
            "ATOM      4  O   ALA A   1       1.291   2.394   0.000  1.00  0.00           O\n"
            "ATOM      5  CB  ALA A   1       2.000  -0.766   1.200  1.00  0.00           C\n"
            "ATOM      6  N   ALA A   2       3.328   1.545   0.000  1.00  0.00           N\n"
            "ATOM      7  CA  ALA A   2       4.032   2.845   0.000  1.00  0.00           C\n"
            "ATOM      8  C   ALA A   2       5.523   2.600   0.000  1.00  0.00           C\n"
            "ATOM      9  O   ALA A   2       6.333   3.528   0.000  1.00  0.00           O\n"
            "ATOM     10  CB  ALA A   2       3.650   3.642   1.250  1.00  0.00           C\n"
            "ATOM     11  N   ALA A   3       5.860   1.319   0.000  1.00  0.00           N\n"
            "ATOM     12  CA  ALA A   3       7.247   0.908   0.000  1.00  0.00           C\n"
            "ATOM     13  C   ALA A   3       7.550  -0.584   0.000  1.00  0.00           C\n"
            "ATOM     14  O   ALA A   3       6.657  -1.432   0.000  1.00  0.00           O\n"
            "ATOM     15  CB  ALA A   3       7.800   1.600   1.250  1.00  0.00           C\n"
            "TER\n"
            "END"
        )
        with open(cls.pdb_path, "w") as f:
            f.write(pdb_content)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.test_dir)

    def setUp(self):
        self.essa = ESSA2()

    def test_set_system(self):
        """Test initialization and GNM calculation."""
        self.essa.setSystem(self.pdb_path, num_modes=2)
        self.assertTrue(self.essa._system_ready)
        self.assertEqual(len(self.essa._ca), 3)
        self.assertIsNotNone(self.essa._ref_msfs)
        self.assertEqual(len(self.essa._ref_msfs), 3)

    def test_scan_residues(self):
        """Test the mutation scanning process."""
        self.essa.setSystem(self.pdb_path, num_modes=2)
        # Scan only the first residue to save time
        msf_matrix = self.essa.scanResidues(residue_indices=[1], num_rotamers=1)
        
        self.assertEqual(msf_matrix.shape, (1, 3))
        self.assertIsNotNone(self.essa._msf_matrix)

    def test_zscore_calculation(self):
        """Test Z-score generation from scanned MSFs."""
        self.essa.setSystem(self.pdb_path, num_modes=2)
        # Mocking a small MSF matrix manually for logic check
        self.essa._msf_matrix = np.array([[0.1, 0.2, 0.1]])
        self.essa._ref_msfs = np.array([0.5, 0.5, 0.5])
        
        zscores = self.essa.zscores()
        self.assertEqual(len(zscores), 1)
        # In this 1-row case, std is 0, logic returns avg_diff
        self.assertAlmostEqual(zscores[0], (0.4 + 0.3 + 0.4) / 3.0)

    def test_io_operations(self):
        """Test saving and loading Z-scores."""
        self.essa.setSystem(self.pdb_path, num_modes=2)
        self.essa._zscores = np.array([1.0, 2.0, 3.0])
        
        save_path = os.path.join(self.test_dir, "test_zs.npy")
        self.essa.saveZscores(save_path)
        
        loaded_zs = self.essa.loadZscores(save_path)
        np.testing.assert_array_equal(loaded_zs, self.essa._zscores)

    def test_ligand_residue_identification(self):
        """Test that ligand selection logic correctly maps indices."""
        # Using the same residues as 'ligand' by mocking proximity
        # Note: In a real PDB, this would look for atoms near a specified ligand.
        self.essa.setSystem(self.pdb_path, lig='A 1', lig_dist=10.0)
        self.assertIn('A1', self.essa._ligres_idx)

if __name__ == '__main__':
    unittest.main()
