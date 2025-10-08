"""This modeule contains unit tests for :mod:'prody.drugui'."""

from prody import *
from prody import LOGGER
from prody.tests import unittest
from prody.tests.datafiles import parseDatafile, pathDatafile, DATA_FILES
from prody.utilities import which
from prody.drugui import drugui_prepare, drugui_analysis


LOGGER.verbosity = 'none'

DATA = DATA_FILES['drugui']

class TestDruGUI(unittest.TestCase):

    def setUp(self):
        """Set PDB and PSF file data"""

        self.pdb = DATA_FILES['mdm2_pdb']
        self.psf = DATA_FILES['mdm2_psf']
        self.outdir_location = DATA_FILES['drugui_outdir']
        
    def testUsualCase(self):
        """Test the out come of preparing and analyzing a simple druggability simulation"""
        vmd_executable = which('vmd_MACOSXARM64')

        if vmd_executable == None:
            vmd_executable = which('vmd_LINUXAMD64')

        drugui_prepare(pdb = f"{self.pdb['file']}", psf = f"{self.psf['file']}", prefix = "md", Probes = True, 
                       probe_and_percent = {"IPRO": 16, "IMID": 14, "ACTT": 14, "ACAM": 14, "IPAM": 14, "IBTN": 14, "BENZ": 14},
                       solvent_padding = 15, boundary_padding = 2.4, neutralize = True, lipids = False, write_conf = True, 
                       n_sims = 4, sim_length = 40, outdir_location = f"{self.outdir_location['outdir_location']}", constrain = "heavy",
                       vmd = vmd_executable)