"""This module contains unit tests for :mod:`~prody.proteins`."""

from collections import OrderedDict
import os

import numpy as np
from numpy.testing import *
try:
    import numpy.testing.decorators as dec
except ImportError:
    from numpy.testing import dec

from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'

class TestParseMMCIF(unittest.TestCase):

    def setUp(self):
        """Set MMCIF file data and parse the MMCIF file."""
        self.multi = DATA_FILES['multi_model_cif']
        self.no_pdb = DATA_FILES['long_chid_cif']
        self.biomols = DATA_FILES['biomols_cif']

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
                        'when chain A is specified with unite_chain default (False)')
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

