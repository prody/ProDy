#!/usr/bin/env python3
import unittest
import os
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open

# 1. Import the module under test
from prody.drugui import no_gui

# We mock the internal 'probe' module entirely to avoid any imports from DruGUI-script
mock_probe = MagicMock()
mock_probe.PROBETYPES = {'core': ['BENZ', 'IPRO']}
mock_probe._LOADED = True

# Save the original exists to allow checking our actual temp directory
original_exists = os.path.exists

def isolated_exists_mock(path):
    """
    Returns True for any DruGUI-related file check.
    Returns real disk status for the temp protein files.
    """
    path_str = str(path).lower()
    
    # Keyword list based on no_gui.py file validation logic
    drugui_files = [
        'probedata.txt', 'probev2.dat', 'top_all36_cgenff', 
        'par_all36_cgenff', '.psf', '.pdb', '.top', '.prm'
    ]
    
    # If checking our actual temp pdb/psf, use the real disk
    if 'protein.' in path_str:
        return original_exists(path)
        
    # Otherwise, fake it for DruGUI
    if any(k in path_str for k in drugui_files):
        return True
        
    return original_exists(path)

class TestDruguiIsolated(unittest.TestCase):
    def setUp(self):
        # Create a real temp location for the "user" protein files
        self.temp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.temp_dir, ignore_errors=True)

        self.pdb = str(Path(self.temp_dir) / "protein.pdb")
        self.psf = str(Path(self.temp_dir) / "protein.psf")
        Path(self.pdb).touch()
        Path(self.psf).touch()

        # START PATCHES
        # Patch the probe module directly in the no_gui namespace
        patch('prody.drugui.no_gui.druggability.probe', mock_probe, create=True).start()
        
        # Patch file existence to bypass disk requirements
        patch('os.path.exists', side_effect=isolated_exists_mock).start()

        # Mock open() to provide content for probedata.txt
        # and prevent writing real files to protected areas
        self.m_open = mock_open(read_data="IPRO 1.0\nBENZ 1.0\n")
        patch('builtins.open', self.m_open).start()

        self.addCleanup(patch.stopall)

    @patch("prody.drugui.no_gui.subprocess.run")
    @patch("prody.drugui.no_gui.os.chdir")
    def test_missing_vmd_executable_raises(self, mock_chdir, mock_run):
        """Tests that empty VMD path triggers ValueError even if files 'exist'"""
        with self.assertRaises(ValueError) as cm:
            no_gui.drugui_prepare(
                pdb=self.pdb, psf=self.psf, vmd="", 
                outdir_location=self.temp_dir, prefix="test"
            )
        self.assertIn("VMD executable", str(cm.exception))

    @patch("prody.drugui.no_gui.subprocess.run")
    @patch("prody.drugui.no_gui.os.chdir")
    def test_bad_probe_percentage_raises(self, mock_chdir, mock_run):
        """Tests percentage validation logic by bypassing file and VMD hurdles"""
        with self.assertRaises(ValueError) as cm:
            no_gui.drugui_prepare(
                pdb=self.pdb, 
                psf=self.psf, 
                vmd="/fake/vmd/path", # Mocked path, doesn't need to exist
                Probes=True,
                probes_and_percent={"IPRO": 50, "BENZ": 10}, # Sums to 60 (Invalid)
                outdir_location=self.temp_dir, 
                prefix="testbad"
            )
        self.assertIn("percentages must sum up to 100", str(cm.exception).lower())

if __name__ == "__main__":
    unittest.main()
