#!/usr/bin/python
# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module contains unit tests for :mod:`~prody.proteins`."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os

import unittest
import numpy as np
from numpy.testing import *

from prody import *
from prody.tools import *
from test_datafiles import *

prody.setVerbosity('none')

class TestFetchPDB(unittest.TestCase):
    
    """Test :func:`~prody.proteins.fetchPDB` function."""
    
    @dec.slow
    def setUp(self):
        """Instantiate a list for storing downloaded file names."""
        
        self.filenames = []
    
    @dec.slow
    def testFetchingSingleFile(self):
        """Test the outcome of fetching a single PDB file."""
        
        fn = fetchPDB('1p38', folder=TEMPDIR, copy=True)
        self.assertTrue(os.path.isfile(fn), 
            'fetching a single PDB file failed')
        self.filenames.append(fn)
        
    @dec.slow
    def testFetchingMultipleFiles(self): 
        """Test the outcome of fetching multiple PDB files."""
        
        fns = fetchPDB(['1p38', '1r39'], folder=TEMPDIR, copy=True)
        self.assertIsInstance(fns, list, 
            'failed to return a list of filenames')
        self.assertTrue(all([os.path.isfile(fn) for fn in fns]),
            'fetching multiple PDB files failed')
        self.filenames.extend(fns)
    
    @dec.slow
    def testCompressedArgument(self):
        """Test decompressing fetched PDB files."""
        
        fns = fetchPDB(['1p38', '1r39'], folder=TEMPDIR, compressed=False)
        self.assertTrue(all([os.path.isfile(fn) for fn in fns]),
            'fetching decompressed PDB files failed')
        self.assertTrue(all([os.path.splitext(fn)[1] != '.gz' for fn in fns]),
            'decompressing PDB files failed')
        self.filenames.extend(fns)

    
    def testInvalidPDBIdentifier(self):
        """Test outcome of passing invalid PDB identifiers."""
        
        self.assertIsNone(fetchPDB('XXXXX', folder=TEMPDIR),
            'failed to return None for an invalid PDB identifier')
            
        self.assertFalse(all(fetchPDB(
                    ['XXXXX', '654654', '-/-*/+', ''], folder=TEMPDIR)),
            'failed to return None for invalid PDB identifiers')
    @dec.slow    
    def tearDown(self):
        """Remove downloaded files from disk."""
        
        for fn in self.filenames:
            if os.path.isfile(fn):
                os.remove(fn)
        
class TestParsePDB(unittest.TestCase):
    
    def setUp(self):
        """Set PDB file data and parse the PDB file."""
        
        self.pdb = DATA_FILES['multi_model_truncated']
        self.one = DATA_FILES['oneatom']
        self.ca = DATA_FILES['1ubi']
         
    def testUsualCase(self):
        """Test the outcome of a simple parsing scenario."""

        ag = parseDatafile(self.pdb['file'])
        
        self.assertIsInstance(ag, prody.AtomGroup,
            'parsePDB failed to return an AtomGroup instance')
    
        self.assertEqual(ag.numAtoms(), self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms')
    
        self.assertEqual(ag.numCoordsets(), self.pdb['models'],
            'parsePDB failed to parse correct number of coordinate sets '
            '(models)')

        self.assertEqual(ag.getTitle(), 
             os.path.splitext(self.pdb['file'])[0],
            'failed to set AtomGroup title based on filename')

    def testPDBArgument(self):
        """Test outcome of invalid *pdb* arguments."""
        
        self.assertRaises(IOError, parsePDB, self.pdb['file'] + '.gz')
        self.assertRaises(TypeError, parsePDB, None)

    def testModelArgument(self):
        """Test outcome of valid and invalid *model* arguments."""
        
        path = getDatafilePath(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, model='0')
        self.assertRaises(ValueError, parsePDB, path, model=-1)
        self.assertRaises(proteins.PDBParseError, parsePDB, path, 
                          model=self.pdb['models']+1)
        self.assertIsNone(parsePDB(path, model=0),
            'parsePDB failed to parse no coordinate sets')
        self.assertEqual(parsePDB(path, model=1).numCoordsets(), 1,
            'parsePDB failed to parse the first coordinate set')
        self.assertEqual(parsePDB(path, model=self.pdb['models'])
            .numCoordsets(), 1,
            'parsePDB failed to parse the last coordinate set')

    def testTitleArgument(self):
        """Test outcome of *title* argument."""
        
        path = getDatafilePath(self.pdb['file'])
        title = 'small protein'    
        self.assertEqual(parsePDB(path, title=title).getTitle(), 
             title, 'parsePDB failed to set user given title')

        name = 1999
        self.assertEqual(parsePDB(path, title=title).getTitle(), 
             str(title), 'parsePDB failed to set user given non-string name')
            
    def testChainArgument(self):
        """Test outcome of valid and invalid *chain* arguments."""
        
        path = getDatafilePath(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, chain=['A'])
        self.assertRaises(ValueError, parsePDB, path, chain='')
        self.assertIsNone(parsePDB(path, chain='$'))
        self.assertEqual(parsePDB(path, chain='A')
            .numAtoms(), self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms when chain is '
            'specified')

    def testSubsetArgument(self):
        """Test outcome of valid and invalid *subset* arguments."""

        path = getDatafilePath(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, subset=['A'])
        self.assertRaises(ValueError, parsePDB, path, subset='')
        self.assertEqual(parsePDB(path, subset='ca').numAtoms(), 10,
                        'failed to parse correct number of "ca" atoms')
        self.assertEqual(parsePDB(path, subset='bb').numAtoms(), 40,
                        'failed to parse correct number of "bb" atoms')

    def testAgArgument(self):
        """Test outcome of valid and invalid *ag* arguments."""

        path = getDatafilePath(self.pdb['file'])
        self.assertRaises(TypeError, parsePDB, path, ag='AtomGroup')
        ag = prody.AtomGroup('One atom')
        ag.setCoords(prody.np.array([[0, 0, 0]]))
        self.assertRaises(ValueError, parsePDB, path, ag=ag)
        ag = prody.AtomGroup('Test')
        self.assertEqual(parsePDB(path, ag=ag).numAtoms(), 
            self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms')
    
    def testAgArgMultiModel(self):
        """Test number of coordinate sets when using *ag* arguments."""
        
        path = getDatafilePath(self.pdb['file'])
        ag = parsePDB(path)
        coords = ag.getCoordsets()
        ncsets = ag.numCoordsets()
        ag = parsePDB(path, ag=ag)
        self.assertEqual(ag.numCoordsets(), ncsets*2,
            'parsePDB failed to append coordinate sets to given ag')
        assert_equal(coords, ag.getCoordsets(np.arange(ncsets, ncsets*2)))
        
    def testAgArgSingleModel(self):
        """Test number of coordinate sets when using *ag* arguments."""
        
        path = getDatafilePath(self.ca['file'])
        ag = parsePDB(path)
        coords = ag.getCoordsets()
        ncsets = ag.numCoordsets()
        ag = parsePDB(path, ag=ag, subset='bb')
        self.assertEqual(ag.numCoordsets(), ncsets*2,
            'parsePDB failed to append coordinate sets to given ag')
        assert_equal(coords, ag.getCoordsets(np.arange(ncsets, ncsets*2)))

    
    def testBiomolArgument(self):
        
        self.assertRaises(proteins.PDBParseError, parsePDB, self.one['path'], 
                          biomol=True)


    def testSecondaryArgument(self):

        self.assertRaises(proteins.PDBParseError, parsePDB, self.one['path'], 
                          secondary=True)

class TestWritePDB(unittest.TestCase):
    
    @dec.slow
    def setUp(self):
        """Set PDB file data and parse the PDB file."""
        
        self.pdb = DATA_FILES['multi_model_truncated']
        self.ag = parsePDB(self.pdb['path'])
        self.tmp = os.path.join(TEMPDIR, 'test.pdb')

    msg = 'user does not have write access to temp dir {0:s}'.format(TEMPDIR) 
    
    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testParsingOutput(self):
        """Test if parsing output is the same as parsing original file."""
        
        out = writePDB(self.tmp, self.ag)
        self.assertEqual(self.tmp, out,
            'writePDB failed to return correct output filename')
        self.assertTrue(os.path.isfile(out),
            'writePDB failed to write output')
        out = parsePDB(out)
        self.assertEqual(self.ag.numAtoms(), out.numAtoms(),
            'writePDB failed to write correct number of atoms')                
        self.assertEqual(self.ag.numCoordsets(), out.numCoordsets(),
            'writePDB failed to write correct number of atoms')
            
    @dec.slow
    @unittest.skipUnless(os.access(TEMPDIR, os.W_OK), msg)
    def testModelArgument(self):
        """Test valid and invalid model arguments and if specified model
        is correctly written."""
        
        self.assertRaises(TypeError, writePDB, self.tmp, self.ag, model='s')
        self.assertRaises(ValueError, writePDB, self.tmp, self.ag, model=-1)
        self.assertRaises(ValueError, writePDB, self.tmp, self.ag, model=0)
        for i in range(self.ag.numCoordsets()):
            out = parsePDB(writePDB(self.tmp, self.ag, model=i+1))
            self.assertEqual(out.numCoordsets(), 1,
                'writePDB failed to write correct number of models')
            self.assertTrue(np.all(out.getCoords() == self.ag.getCoordsets(i)),
                'writePDB failed to write coordinates correctly')
                
    @dec.slow
    def tearDown(self):
        """Remove test file."""
        
        if os.path.isfile(self.tmp):
            os.remove(self.tmp)


class TestParsePDBHeaderOnly(unittest.TestCase):
    
    def setUp(self):
        self.header = parsePDB(getDatafilePath('pdb2k39_truncated.pdb'), 
                               header=True, model=0)

    def testHeaderType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')
        
    def testHeaderContent(self):
        self.assertEqual(self.header.get('classification'), 
            'SIGNALING PROTEIN',
            'failed to get expected value for classification from header')
        self.assertEqual(self.header.get('experiment'), 'SOLUTION NMR',
            'failed to get expected value for experiment from header')
        self.assertEqual(self.header.get('deposition_date'), '25-APR-08',
            'failed to get expected value for deposition_date from header')
        self.assertEqual(self.header.get('identifier'), '2K39',
            'failed to get expected value for identifier from header')
        self.assertEqual(self.header.get('title'), 
            'RECOGNITION DYNAMICS UP TO MICROSECONDS REVEALED FROM RDC '
            'DERIVED UBIQUITIN ENSEMBLE IN SOLUTION',
            'failed to get expected value for title from header dictionary')

    def tearDown(self):
        
        self.header = None


class TestParsePDBHeaderAndAllModels(unittest.TestCase):

    def setUp(self):
        self.atomgroup, self.header = \
            parsePDB(getDatafilePath('pdb2k39_truncated.pdb'), header=True)

    def testAtomGroupType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')
        self.assertIsInstance(self.atomgroup, prody.AtomGroup,
            'atom group type is incorrect')
        
    def testAtomGroupContent(self):
        
        self.assertEqual(self.atomgroup.numAtoms(), 167,
            'incorrect number of atoms')
        self.assertEqual(self.atomgroup.numCoordsets(), 3,
            'incorrect number of coordinate sets (models)')

    def tearDown(self):
        
        self.header = None
        self.atomgroup = None


class TestParsePDBAltloc(unittest.TestCase):
    
    def setUp(self):
        
        self.pdbfile = getDatafilePath('pdb1ejg.pdb')
    
    def testAltlocNone(self):
        
        self.assertEqual(len(parsePDB(self.pdbfile)), 637,
            'failed to parse unspecified alternate locations correctly')
    
    def testAltlocA(self):
        
        self.assertEqual(len(parsePDB(self.pdbfile, altloc='A')), 637,
            'failed to parse alternate locations A correctly')
        
    def testAltlocB(self):
        
        self.assertEqual(len(parsePDB(self.pdbfile, altloc='B')), 634,
            'failed to parse alternate locations B correctly')

    def testAltlocC(self):
        
        self.assertEqual(len(parsePDB(self.pdbfile, altloc='C')), 496,
            'failed to parse alternate locations C correctly')



class TestParsePDBBiomolecule(unittest.TestCase):
    pass

class TestParsePDBSecondary(unittest.TestCase):
    pass

class TestParsePSF(unittest.TestCase):
    pass

class TestParsePSFandPDB(unittest.TestCase):
    pass
 

class TestDSSPFunctions(unittest.TestCase):
    
    @dec.slow
    def setUp(self):
        """Setup the testing framework."""

        self.pdbs = [DATA_FILES['dssp']]
    
    @dec.slow
    @unittest.skipIf(which('dssp') is None, 'dssp is not found')
    def testDSSPBridgePartners(self):
        """Check if the DSSP bridge-partners were correctly parsed and 
        assigned."""

        for pdb in self.pdbs:
            prot_ag = parseDatafile(pdb['file'], folder=TEMPDIR)
            dssp = execDSSP(getDatafilePath(pdb['file']), outputdir=TEMPDIR, 
                            stderr=False)
            parseDSSP(dssp, prot_ag, parseall=True)
    
            # Map a dssp_resnum to its Residue object.
            dssp_dict = {}

            for chain in prot_ag.select("protein").getHierView():
                for res in chain:
                    dssp_resnum = res.getData("dssp_resnum")[0]
                    dssp_dict[dssp_resnum] = res

            for res in dssp_dict.itervalues():
                bp1 = res.getData("dssp_bp1")[0]
                bp2 = res.getData("dssp_bp2")[0]

                if bp1 != 0:
                    msg_ = "BP1 (dssp_resnum: %d) of %s is missing" % \
                        (bp1, str(res))
                    self.assertTrue(dssp_dict.has_key(bp1), msg=msg_)

                if bp2 != 0:
                    msg_ = "BP2 (dssp_resnum: %d) of %s is missing" % \
                        (bp2, str(res))
                    self.assertTrue(dssp_dict.has_key(bp2), msg=msg_)


if __name__ == '__main__':
    unittest.main()
