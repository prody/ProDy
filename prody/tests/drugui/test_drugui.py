"""This modeule contains unit tests for :mod:'prody.drugui'."""

from prody import *
from prody import LOGGER
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *
from prody.utilities import which
from prody.drugui import drugui_prepare, drugui_analysis, drugui_evaluate
import os


LOGGER.verbosity = 'none'


class TestDruGUI(unittest.TestCase):

    def setUp(self):
        """Set PDB and PSF file data"""

        self.pdb = {'file': pathDatafile('mdm2.pdb')}
        self.psf = {'file': pathDatafile('mdm2.psf')}
        self.pdb1 = {'file': pathDatafile('md.pdb')}
        self.psf1 = {'file': pathDatafile('md.psf')}
        self.dcd = {'file': pathDatafile('final.dcd')}
        self.pdb2 = {'file': pathDatafile('ligand.pdb')}
        
    def testUsualCase(self):
        """Test the out come of preparing and analyzing a simple druggability simulation"""
        vmd_executable = which('vmd_MACOSXARM64')

        if vmd_executable == None:
            vmd_executable = which('vmd_LINUXAMD64')

        outdir = os.path.join(TEMPDIR, "drugui_test_output")
        os.makedirs(outdir, exist_ok=True)

        drugui_prepare(pdb = self.pdb['file'], psf = self.psf['file'], prefix = "md", Probes = True, 
                       probe_and_percent = {"IPRO": 16, "IMID": 14, "ACTT": 14, "ACAM": 14, "IPAM": 14, "IBTN": 14, "BENZ": 14},
                       solvent_padding = 15, boundary_padding = 2.4, neutralize = True, lipids = False, write_conf = True, 
                       n_sims = 4, sim_length = 40, outdir_location = outdir, 
                       constrain = "heavy", vmd = vmd_executable)
        
        md_pdb = os.path.join(outdir, "md/simulation_run1/md.pdb")  
       
        structure = parsePDB(f'{md_pdb}', long_resname = True)

        self.assertIsNotNone(structure.select('protein'), 'DruGUI failed to include protein')
        self.assertIsNotNone(structure.select('water'), 'DruGUI failed to include water molecules')
        self.assertIsNotNone(structure.select('chain I'), 'DruGUI failed to include ions')
        self.assertIsNotNone(structure.select('segment PROB'), 'DruGUI failed to include probe molecules')

        drugui_analysis(pdb = self.pdb1['file'], psf = self.psf1['file'], dcds= self.dcd['file'], prefix = 'dg',
                        outdir_location = outdir, selection = "noh and protein", grid_spacing = 0.5, contact_distance = 4.0, 
                        align = 'calpha', temperature = 300., merge_radius = 5.6, n_frames = 1, n_probes = 7, delta_g = -1.0, 
                        min_n_probes = 6, low_affinity = 10, max_charge = 2, n_solutions = 3, n_charged = 3, 
                        probes = ['IPAM', 'IPRO', 'ACTT', 'IBTN', 'ACAM', 'IMID', 'BENZ'])
        
        dg_protein_heavyatom = os.path.join(outdir, "dg_protein_heavyatoms.pdb")
        
        heavyatom = parsePDB(f'{dg_protein_heavyatom}')

        self.assertIsNotNone(heavyatom.select('protein'), "DruGUI failed to produce protein's heavyatoms")

        probe_dcds = ['dg_IPRO.dcd', 'dg_IMID.dcd', 'dg_ACTT.dcd', 'dg_ACAM.dcd', 'dg_IPAM.dcd', 'dg_IBTN.dcd', 'dg_BENZ.dcd']
        dcd_files = [os.path.join(outdir, name) for name in probe_dcds]

        for dcd in dcd_files:
            probe = parseDCD(f'{dcd}')
            self.assertIsNotNone(probe, 'DruGUI failed to produce probe dcd files')

        solution = os.path.join(outdir, "dg_site_1.pdb")
        solution_pdb = parsePDB(f'{solution}')
        self.assertIsNotNone(solution_pdb, "DruGUI failed to produce druggable sites")


        dso = os.path.join(outdir, "dg.dso.gz")
        drugui_evaluate(pdb =self.pdb2['file'], outdir_location =outdir, dso = dso, radius=1.5, delta_g=-0.5)

        ligand = os.path.join(outdir, "dg_ligand.pdb")
        ligand_pdb = parsePDB(f'{ligand}')
        self.assertIsNotNone(ligand_pdb, "DruGUI failed to evaluate the inhibitor")
        
if __name__ == '__main__':
    unittest.main()
