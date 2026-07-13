"""This module contains unit tests for :mod:`~prody.proteins`."""

from collections import OrderedDict
import os

import numpy as np
from numpy.testing import *

from prody.utilities import importDec
dec = importDec()

from prody import *
from prody import LOGGER
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

class TestParseMMCIF(unittest.TestCase):

    def setUp(self):
        """Set MMCIF file data and parse the MMCIF file."""
        self.multi = DATA_FILES['multi_model_cif']
        self.no_pdb = DATA_FILES['long_chid_cif']
        self.biomols = DATA_FILES['biomols_cif']
        self.big_biomols = DATA_FILES['big_biomols_cif']
        self.chimerax = DATA_FILES['chimerax_cif']
        self.pymol = DATA_FILES['pymol_cif']
        self.wrapped = DATA_FILES['boltz_wrapped_line_cif']
        self.scipion = DATA_FILES['scipion_cif']

        self.altlocs = DATA_FILES['cif_6flr']
        self.his_selstr = 'resname HIS and chain B and resnum 234 and name CA'

    def testUsualCase(self):
        """Test the outcome of a simple parsing scenario."""

        ag = parseDatafile(self.multi['file'])

        self.assertIsInstance(ag, prody.AtomGroup,
            'parseMMCIF failed to return an AtomGroup instance')

        self.assertEqual(ag.numAtoms(), self.multi['atoms'],
            'parseMMCIF failed to parse correct number of atoms')

        self.assertEqual(ag.numCoordsets(), self.multi['models'],
            'parseMMCIF failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(ag.getTitle(),
             os.path.splitext(self.multi['file'])[0],
            'failed to set AtomGroup title based on filename')

    def testPDBArgument(self):
        """Test outcome of invalid *pdb* arguments."""

        self.assertRaises(IOError, parseMMCIF, self.multi['file'] + '.gz')
        self.assertRaises(TypeError, parseMMCIF, None)

    def testModelArgument(self):
        """Test outcome of valid and invalid *model* arguments."""

        path = pathDatafile(self.multi['file'])
        self.assertRaises(TypeError, parseMMCIF, path, model='0')
        self.assertRaises(ValueError, parseMMCIF, path, model=-1)
        self.assertRaises(proteins.MMCIFParseError, parseMMCIF, path,
                          model=self.multi['models']+1)
        self.assertIsNone(parseMMCIF(path, model=0),
            'parseMMCIF failed to parse no coordinate sets')

        self.assertEqual(parseMMCIF(path, model=1).numCoordsets(), 1,
            'parseMMCIF failed to parse the first coordinate set')

        self.assertEqual(parseMMCIF(path, model=2).numCoordsets(), 1,
            'parseMMCIF failed to parse the 2nd coordinate set')

        self.assertEqual(parseMMCIF(path, model=1).numAtoms(),
                        self.multi['atoms'],
                        'parseMMCIF failed to parse the 1st coordinate set')

        self.assertEqual(parseMMCIF(path, model=2).numAtoms(),
                        self.multi['atoms'],
                        'parseMMCIF failed to parse the 2nd coordinate set')
            
        self.assertEqual(parseMMCIF(path, 
                                    model=self.multi['models']).numCoordsets(), 
                        1, 'parseMMCIF failed to parse the last coordinate set')

    def testTitleArgument(self):
        """Test outcome of *title* argument."""

        path = pathDatafile(self.multi['file'])
        title = 'small protein'
        self.assertEqual(parseMMCIF(path, title=title).getTitle(),
             title, 'parseMMCIF failed to set user given title')

        name = 1999
        self.assertEqual(parseMMCIF(path, title=name).getTitle(),
             str(name), 'parseMMCIF failed to set user given non-string name')

    def testChainArgument(self):
        """Test outcome of valid and invalid *chain* arguments."""

        path = pathDatafile(self.multi['file'])
        self.assertRaises(TypeError, parseMMCIF, path, chain=['A'])
        self.assertRaises(ValueError, parseMMCIF, path, chain='')
        self.assertIsNone(parseMMCIF(path, chain='$'))
        self.assertEqual(parseMMCIF(path, chain='A').numAtoms(),
                        self.multi['chainA_atoms'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'when chain is specified')

    def testLongChainArgument(self):
        """Test outcome of valid and invalid *segment* arguments."""

        path = pathDatafile(self.no_pdb['file'])
        self.assertRaises(TypeError, parseMMCIF, path, segment=['SX0'])
        self.assertRaises(ValueError, parseMMCIF, path, segment='')
        self.assertIsNone(parseMMCIF(path, segment='$'))
        self.assertEqual(parseMMCIF(path, segment='SX0').numAtoms(),
                        self.no_pdb['segment_SX0_atoms'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'when segment SX0 is specified')
        
    def testUniteChainsArgument(self):
        """Test outcome of valid and invalid *segment* arguments."""

        path = pathDatafile(self.biomols['file'])
        self.assertEqual(parseMMCIF(path, chain='A').numAtoms(),
                        self.biomols['chainA_atoms_alone'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'when chain A is specified with unite_chain default (False)')
        self.assertEqual(parseMMCIF(path, chain='A', unite_chains=True).numAtoms(),
                        self.biomols['chainA_atoms_united'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'when chain A is specified with unite_chain True')
        self.assertEqual(parseMMCIF(path, chain='A', header=True)[0].numAtoms(),
                        self.biomols['chainA_atoms_alone'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'when chain A is specified with unite_chain default (False) '
                        'with header True')
        self.assertEqual(parseMMCIF(path, chain='A', header=True, unite_chains=True)[0].numAtoms(),
                        self.biomols['chainA_atoms_united'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'when chain A is specified with unite_chain True '
                        'with header True')

    def testUniteChainsAndBiomolArguments(self):
        """Test outcome of valid and invalid *segment* arguments."""

        path = pathDatafile(self.biomols['file'])

        bm_united = parseMMCIF(path, biomol=True, unite_chains=True)
        self.assertEqual(bm_united[0].numAtoms(),
                        self.biomols['bm0_atoms'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'with biomol True and unite_chain True')
        self.assertEqual([b.numChains() for b in bm_united], 
                        self.biomols['bm_chains_united'],
                        'parseMMCIF failed to parse correct numbers of chains '
                        'with biomol True and unite_chain True')

        bm_non_united = parseMMCIF(path, biomol=True)
        self.assertEqual(bm_non_united[0].numAtoms(),
                        self.biomols['bm0_atoms'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'with biomol True unite_chain default (False)')
        self.assertEqual([b.numChains() for b in bm_non_united], 
                        self.biomols['bm_chains_alone'],
                        'parseMMCIF failed to parse correct numbers of chains '
                        'with biomol True and unite_chain default (False)')

        bm_header = parseMMCIF(path, biomol=True, header=True, unite_chains=True)[0]
        self.assertEqual(bm_header[0].numAtoms(),
                        self.biomols['bm0_atoms'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'with biomol True and unite_chain True')
        self.assertEqual([b.numChains() for b in bm_header], 
                        self.biomols['bm_chains_united'],
                        'parseMMCIF failed to parse correct numbers of chains '
                        'with biomol True and unite_chain True')

    def testBiomolOperationRange(self):
        """Test outcome of valid and invalid *segment* arguments."""

        path = pathDatafile(self.big_biomols['file'])

        non_bm = parseMMCIF(path, biomol=False)
        self.assertEqual(non_bm.numAtoms(),
                         self.big_biomols['atoms'],
                         'parseMMCIF failed to parse correct number of atoms '
                         'for 7cth with biomol False')
        self.assertEqual(non_bm.numChains(),
                         self.big_biomols['num_chains'],
                        'parseMMCIF failed to parse correct numbers of chains '
                        'for 7cth with biomol False')

        bm_header = parseMMCIF(path, biomol=True, header=True)[0]
        self.assertEqual(bm_header[0].numAtoms(),
                         self.big_biomols['bm0_atoms'],
                         'parseMMCIF failed to parse correct number of atoms '
                         'for 7cth with biomol True')
        self.assertEqual(bm_header[0].numChains(),
                         self.big_biomols['bm0_chains'],
                         'parseMMCIF failed to parse correct number of chains '
                         'for 7cth with biomol True')

    def testWrappedLines(self):
        """Test that we can handle wrapped lines, as generatd by Boltz and Chai"""
        path = pathDatafile(self.wrapped['file'])

        prot = parseMMCIF(path)
        self.assertEqual(prot.numAtoms(),
                         self.wrapped['atoms'],
                         'parseMMCIF failed to parse correct number of atoms '
                         f'for {self.wrapped["file"]}')
        self.assertEqual(prot.numChains(),
                         self.wrapped['num_chains'],
                        'parseMMCIF failed to parse correct numbers of chains '
                         f'for {self.wrapped["file"]}')       

    def testNonpolyAuthSeqIdDot(self):
        """Test that non-polymer atoms can fall back when auth_seq_id is not numeric."""

        path = os.path.join(TEMPDIR, 'openfold3_auth_seq_id_dot.cif')
        with open(path, 'w') as out:
            out.write("""data_openfold3
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A CA ALA A 1 1 ? 1.0 2.0 3.0 1.00 10.00 1
ATOM 2 N N . ALA A N ALA A 1 1 ? 1.1 2.1 3.1 1.00 11.00 1
HETATM 3 C C1 . LIG B C1 LIG B . . ? 4.0 5.0 6.0 1.00 12.00 1
HETATM 4 O O1 . LIG B O1 LIG B . . ? 4.1 5.1 6.1 1.00 13.00 1
#
loop_
_pdbx_nonpoly_scheme.asym_id
_pdbx_nonpoly_scheme.entity_id
_pdbx_nonpoly_scheme.mon_id
_pdbx_nonpoly_scheme.ndb_seq_num
_pdbx_nonpoly_scheme.pdb_seq_num
_pdbx_nonpoly_scheme.auth_seq_num
_pdbx_nonpoly_scheme.pdb_mon_id
_pdbx_nonpoly_scheme.auth_mon_id
_pdbx_nonpoly_scheme.pdb_strand_id
_pdbx_nonpoly_scheme.pdb_ins_code
B 2 LIG 1 1 1 LIG LIG B ?
#
""")

        ag = parseMMCIF(path)
        self.assertEqual(ag.numAtoms(), 4,
                        'parseMMCIF failed to parse atoms when non-polymer auth_seq_id is "."')

        ligand = ag.select('chain B and resname LIG')
        self.assertEqual(ligand.numAtoms(), 2,
                        'parseMMCIF failed to keep non-polymer atoms with auth_seq_id "."')
        assert_equal(ligand.getResnums(), np.array([1, 1]),
                     err_msg='parseMMCIF failed to use non-polymer scheme numbering when auth_seq_id is "."')
        
    def testChimeraxCIFBiomolArguments(self):
        """Test outcome of valid and invalid *segment* arguments."""

        path = pathDatafile(self.chimerax['file'])

        bm_non_united = parseMMCIF(path, biomol=True)
        self.assertEqual(bm_non_united[0].numAtoms(),
                        self.chimerax['bm0_atoms'],
                        'parseMMCIF failed to parse correct number of atoms '
                        'with biomol True for the chimerax mmcif of 1ake')

    def testSubsetArgument(self):
        """Test outcome of valid and invalid *subset* arguments."""

        path = pathDatafile(self.multi['file'])
        self.assertRaises(TypeError, parseMMCIF, path, subset=['A'])
        self.assertEqual(parseMMCIF(path, subset='ca').numAtoms(),
                        self.multi['ca_atoms'],
                        'failed to parse correct number of "ca" atoms')
        self.assertEqual(parseMMCIF(path, subset='bb').numAtoms(),
                        self.multi['bb_atoms'],
                        'failed to parse correct number of "bb" atoms')

    def testAgArgument(self):
        """Test outcome of valid and invalid *ag* arguments."""

        path = pathDatafile(self.multi['file'])
        self.assertRaises(TypeError, parseMMCIF, path, ag='AtomGroup')
        ag = prody.AtomGroup('One atom')
        ag.setCoords(np.array([[0., 0., 0.]]))
        self.assertRaises(ValueError, parseMMCIF, path, ag=ag)
        ag = prody.AtomGroup('Test')
        self.assertEqual(parseMMCIF(path, ag=ag).numAtoms(),
            self.multi['atoms'],
            'parseMMCIF failed to parse correct number of atoms')

    def testAgArgMultiModel(self):
        """Test number of coordinate sets when using *ag* arguments."""

        path = pathDatafile(self.multi['file'])
        ag = parseMMCIF(path)
        coords = ag.getCoordsets()
        ncsets = ag.numCoordsets()
        ag = parseMMCIF(path, ag=ag)
        self.assertEqual(ag.numCoordsets(), ncsets*2,
            'parseMMCIF failed to append coordinate sets to given ag')
        assert_equal(coords, ag.getCoordsets(np.arange(ncsets, ncsets*2)))

    def testUnobsHeaderArgument(self):
        """Test outcome of valid and invalid *subset* arguments."""

        path = pathDatafile(self.biomols['file'])
        header = parseCIFHeader(path)
        self.assertIsInstance(header, dict,
            'parseCIFHeader failed to return an dict instance')
        self.assertTrue('unobserved' in header,
                        'parseCIFHeader failed to return a header containing '
                        'unobserved')     

        unobs_header = parseCIFHeader(path, 'unobserved')
        self.assertIsInstance(unobs_header, OrderedDict,
            'parseCIFHeader failed to return an OrderedDict instance when '
            'providing a key')
        self.assertEqual(len(unobs_header), 
                        self.biomols['num_chains_united'],
                        'failed to parse unobserved for correct number of chains')
        self.assertEqual(unobs_header['B'][0][:40], 
                        self.biomols['unobs_B_start'],
                        'failed to parse unobserved alignment correctly')

    def testAltlocAllToMultiAtoms(self):
        """Test number of coordinate sets and atoms with altloc='all'."""

        path = pathDatafile(self.altlocs['file'])

        ag = parsePDB(path, altloc="all")
        self.assertEqual(ag.numAtoms(), self.altlocs['atoms_altloc'],
            'parsePDB failed to parse correct number of atoms with altloc "all"')
        self.assertEqual(ag.numCoordsets(), 1,
            'parsePDB failed to parse correct number of coordsets (1) with altloc "all"')

        hisB234 = ag.select(self.his_selstr)
        self.assertEqual(hisB234.numAtoms(), self.altlocs['num_altlocs'],
            'parsePDB failed to parse correct number of His B234 CA atoms (2) with altloc "all"')

        self.assertEqual(hisB234.getAnisous().shape, (self.altlocs['num_altlocs'], 6),
            'parsePDB failed to have right shape for His B234 CA atoms getAnisous (2, 6) with altloc "all"')

        assert_allclose(hisB234.getAnisous()[0], self.altlocs['anisousA'][0],
            err_msg='parsePDB failed to have right His B234 CA atoms getAnisous A with altloc "all"')

        assert_allclose(hisB234.getAnisous()[1], self.altlocs['anisousB'][0],
            err_msg='parsePDB failed to have right His B234 CA atoms getAnisous B with altloc "all"')

    def testAltlocAllPymol(self):
        """Test number of coordinate sets and atoms for PyMOL CIF file with altloc='all'."""

        path = pathDatafile(self.pymol['file'])

        ag = parsePDB(path, altloc="all")
        self.assertEqual(ag.numAtoms(), self.pymol['atoms'],
            'parsePDB failed to parse correct number of atoms from pymol cif with altloc "all"')
        self.assertEqual(ag.numCoordsets(), 1,
            'parsePDB failed to parse correct number of coordsets (1) from pymol cif with altloc "all"')
        assert_allclose(ag.getCoords()[-1], self.pymol['last_coords'],
            err_msg='parsePDB failed to parse correct last coords from pymol cif with altloc "all"')

    def testAltlocAPymol(self):
        """Test number of coordinate sets and atoms for PyMOL CIF file with altloc='all'."""

        path = pathDatafile(self.pymol['file'])

        ag = parsePDB(path, altloc="A")
        self.assertEqual(ag.numAtoms(), self.pymol['atoms'],
            'parsePDB failed to parse correct number of atoms from pymol cif with altloc "all"')
        self.assertEqual(ag.numCoordsets(), 1,
            'parsePDB failed to parse correct number of coordsets (1) from pymol cif with altloc "all"')
        assert_allclose(ag.getCoords()[-1], self.pymol['last_coords'],
            err_msg='parsePDB failed to parse correct last coords from pymol cif with altloc "all"')

    def testAltlocNoneToLessAtoms(self):
        """Test number of coordinate sets and atoms with altloc=None."""

        path = pathDatafile(self.altlocs['file'])

        ag = parsePDB(path, altloc=None)
        self.assertEqual(ag.numAtoms(), self.altlocs['atoms_single'],
            'parsePDB failed to parse correct number of atoms with altloc None')
        self.assertEqual(ag.numCoordsets(), 1,
            'parsePDB failed to parse correct number of coordsets (1) with altloc None')

        hisB234 = ag.select(self.his_selstr)
        self.assertEqual(hisB234.numAtoms(), 1,
            'parsePDB failed to parse correct number of His B234 CA atoms (1) with altloc None')

        self.assertEqual(hisB234.getAnisous().shape, (1, 6),
            'parsePDB failed to have right shape for His B234 CA atoms getAnisous (1, 6) with altloc None')

        assert_allclose(hisB234.getAnisous(), self.altlocs['anisousA'],
            err_msg='parsePDB failed to have right His B234 CA atoms getAnisous A with altloc "all"')

    def testAltlocAllMultiModels(self):
        """Test number of coordinate sets and atoms for multi-model case with altloc='all'."""

        path = pathDatafile(self.multi['file'])

        ag = parsePDB(path, altloc="all")
        self.assertEqual(ag.numAtoms(), self.multi['atoms'],
            'parsePDB failed to parse correct number of atoms for multi-model with altloc "all"')
        self.assertEqual(ag.numCoordsets(), self.multi['models'],
            'parsePDB failed to parse correct number of coordsets ({0}) with altloc "all"'.format(self.multi['models']))

    def testStopScipion(self):
        """Test number of coordinate sets and atoms for PyMOL CIF file with altloc='all'."""

        path = pathDatafile(self.scipion['file'])

        ag = parsePDB(path)
        self.assertEqual(ag.numAtoms(), self.scipion['atoms'],
            'parsePDB failed to parse correct number of atoms from pymol cif with altloc "all"')
        self.assertEqual(ag.numCoordsets(), 1,
            'parsePDB failed to parse correct number of coordsets (1) from pymol cif with altloc "all"')
        assert_allclose(ag.getCoords()[-1], self.scipion['last_coords'],
            err_msg='parsePDB failed to parse correct last coords from pymol cif with altloc "all"')
