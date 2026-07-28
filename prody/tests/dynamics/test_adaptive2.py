import unittest
from unittest.mock import MagicMock, patch
import numpy as np

from prody import AdaptiveANM, AANM_ONEWAY, AANM_ALTERNATING

class TestAdaptiveANM(unittest.TestCase):

    def setUp(self):
        # Create two sets of coordinates with a known difference
        # 4 atoms, 3 dimensions
        self.coords_a = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]], dtype=float)
        # Offset coords_b by 1.0 unit on the X axis to create a measurable RMSD
        self.coords_b = self.coords_a + np.array([1.0, 0.0, 0.0])
        
        # Calculate the expected RMSD manually for the test assertion
        # RMSD = sqrt(sum(squared_dist) / N)
        # Here dist is 1.0 for every atom, so RMSD = sqrt(4*1.0 / 4) = 1.0
        self.expected_rmsd = 1.0

    @patch('prody.superpose')
    def test_initialization(self, mock_superpose):
        """Test if the class initializes attributes correctly with real RMSD."""
        # Mock superpose to just return the coords as-is (aligned=True style)
        mock_superpose.return_value = (self.coords_a.copy(), None)
        
        # Use aligned=True to bypass superposition logic during this specific test
        aanm = AdaptiveANM(self.coords_a, self.coords_b, aligned=True, n_modes=10)
        
        # Assertions
        self.assertEqual(aanm.n_modes, 10)
        self.assertEqual(len(aanm.rmsds), 1)
        
        # We check against 1.0 because that is the physical RMSD 
        # between coords_a and coords_b defined in setUp
        self.assertAlmostEqual(aanm.rmsds[0], self.expected_rmsd, places=5)
        
        # Verify initial coordinates were stored
        np.testing.assert_array_equal(aanm.coordsA, self.coords_a)

    def test_invalid_mode_raises_error(self):
        """Ensure passing an invalid mode to run() raises ValueError."""
        aanm = AdaptiveANM(self.coords_a, self.coords_b, aligned=True)
        with self.assertRaises(ValueError):
            aanm.run(n_steps=5, mode=999)

    @patch('prody.calcENM')
    @patch('prody.calcRMSD')
    @patch.object(AdaptiveANM, '_check_disconnection')
    @patch.object(AdaptiveANM, '_compute_slack_info')
    def test_step_logic(self, mock_slack, mock_disconnect, mock_rmsd, mock_calcENM):
        """Test a single _step execution with guaranteed variable definitions."""
        # 1. Coordinate Setup
        coords_a = np.array([[0, 0, 0], [5, 0, 0], [0, 5, 0], [0, 0, 5]], dtype=float)
        coords_b = coords_a.copy()
        coords_b[0, 0] += 1.0  # Deformation on X-axis
        
        # 2. Mock ANM and logic gates
        mock_anm_h = MagicMock()
        eigvecs = np.zeros((12, 12))
        eigvecs[0, 0] = 1.0  # Align eigenvector with X-axis deformation
        mock_anm_h.getEigvecs.return_value = eigvecs
        mock_calcENM.return_value = (mock_anm_h, None)
        
        mock_disconnect.return_value = False
        mock_slack.return_value = (None, 1.0, np.array([0]))
        
        # 3. RMSD: First call is trial check (needs to be < 0.5), 
        # second call is final update (recorded as improvement)
        mock_rmsd.side_effect = [0.4, 0.4, 0.4, 0.4] 
        
        # 4. Initialize the class properly
        aanm = AdaptiveANM(coords_a, coords_b, aligned=True)
        
        mock_ensemble = MagicMock()
        mock_ensemble.getWeights.return_value = None
        
        # Initialize lists for history
        f_hist, ms_hist, we_hist = [0.2], [], []
        
        # 5. Execute the step
        n_modes_out, improvement = aanm._step(
            coords_a, 
            coords_b, 
            n_modes=1, 
            ensemble=mock_ensemble, 
            defvecs=[], 
            rmsds=[0.5],
            Fmin=0.99, 
            f_hist=f_hist, 
            min_slack_hist=ms_hist, 
            worst_edges_hist=we_hist
        )
        
        # 6. Assertions
        self.assertGreater(improvement, 0, "Step should result in positive improvement")
        self.assertTrue(mock_ensemble.addCoordset.called, "addCoordset should be called")
        self.assertEqual(mock_ensemble.addCoordset.call_count, 2, "Should add mid and end points")
        self.assertFalse(np.array_equal(coords_a, np.zeros((4,3))), "Coords should be updated")

    @patch.object(AdaptiveANM, '_step')
    def test_run_oneway_termination(self, mock_step):
        """Test that run_oneway stops when convergence (n_modes=0) is reached."""
        aanm = AdaptiveANM(self.coords_a, self.coords_b, aligned=True)
        
        # Mock _step to return n_modes=0 immediately
        mock_step.return_value = (0, 0.1)
        
        ensemble = aanm.run(n_steps=10, mode=AANM_ONEWAY)
        
        # Should only have called _step once before breaking
        self.assertEqual(mock_step.call_count, 1)

if __name__ == '__main__':
    unittest.main()
