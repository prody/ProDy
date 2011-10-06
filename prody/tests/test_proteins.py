# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

"""This module contains unit tests for :mod:`~prody.proteins` module.
Data files used in tests are truncated PDB files, e.g. most of atoms and/or
models and/or header sections are removed for having a compact installation
package that contains test modules and files as well.  
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import os
import os.path
import sys
import unittest
import tempfile
import prody
from prody.proteins import *
from prody import proteins

TEMPDIR = tempfile.gettempdir()
prody.changeVerbosity('none')

PDB_FILES = {
    'multi_model_truncated': {
        'pdb': '2k39',
        'path': 'data/pdb2k39_m1to3_r1to10.pdb',
        'atoms': 167,
        'models': 3
    },
    'dssp': {
        'pdb': '1r19',
        'path': 'data/pdb1r19_dssp.pdb',
        'atoms': 8216,
        'models': 1
    },
}

class TestFetchPDB(unittest.TestCase):
    
    """Test :func:`~prody.proteins.fetchPDB` function."""
    
    def setUp(self):
        """Instantiate a list for storing downloaded file names."""
        
        self.filenames = []
    
    def testFetchingSingleFile(self):
        """Test the outcome of fetching a single PDB file."""
        
        fn = fetchPDB('1p38', folder=TEMPDIR, copy=True)
        self.assertTrue(os.path.isfile(fn), 
            'fetching a single PDB file failed')
        self.filenames.append(fn)
        
    def testFetchingMultipleFiles(self): 
        """Test the outcome of fetching multiple PDB files."""
        
        fns = fetchPDB(['1p38', '1r39'], folder=TEMPDIR, copy=True)
        self.assertIsInstance(fns, list, 
            'failed to return a list of filenames')
        self.assertTrue(all([os.path.isfile(fn) for fn in fns]),
            'fetching multiple PDB files failed')
        self.filenames.extend(fns)
        
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
        
    def tearDown(self):
        """Remove downloaded files from disk."""
        
        for fn in self.filenames:
            if os.path.isfile(fn):
                os.remove(fn)
        
class TestParsePDB(unittest.TestCase):
    
    def setUp(self):
        """Set PDB file data and parse the PDB file."""
        
        self.pdb = PDB_FILES['multi_model_truncated']
        self.ag = parsePDB(self.pdb['path'])
         
    def testReturnType(self):
        """Test the outcome of a simple parsing scenario."""
        
        self.assertIsInstance(self.ag, prody.AtomGroup,
            'parsePDB failed to return an AtomGroup instance')
    
    def testNumOfAtoms(self):
        """Test the number of parsed atoms."""
        
        self.assertEqual(self.ag.getNumOfAtoms(), self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms')
    
    def testNumOfCoordsets(self):
        """Test the number of parsed coordinate sets."""
        
        self.assertEqual(self.ag.getNumOfCoordsets(), self.pdb['models'],
            'parsePDB failed to parse correct number of coordinate sets '
            '(models)')

    def testAtomGroupName(self):
        """Test the name of the parsed :class:`~prody.atomic.AtomGroup` 
        instance."""
        
        self.assertEqual(self.ag.getName(), 
             os.path.splitext(os.path.split(self.pdb['path'])[1])[0],
            'failed to set AtomGroup name based on filename')

    def testPDBArgument(self):
        """Test outcome of invalid *pdb* arguments."""
        
        self.assertRaises(IOError, parsePDB, self.pdb['path'] + '.gz')
        self.assertRaises(TypeError, parsePDB, None)

    def testModelArgument(self):
        """Test outcome of valid and invalid *model* arguments."""
        
        path = self.pdb['path']
        self.assertRaises(TypeError, parsePDB, path, model='0')
        self.assertRaises(ValueError, parsePDB, path, model=-1)
        self.assertRaises(proteins.PDBParserError, parsePDB, path, 
                          model=self.pdb['models']+1)
        self.assertIsNone(parsePDB(path, model=0),
            'parsePDB failed to parse no coordinate sets')
        self.assertEqual(parsePDB(path, model=1).getNumOfCoordsets(), 1,
            'parsePDB failed to parse the first coordinate set')
        self.assertEqual(parsePDB(path, model=self.pdb['models'])
            .getNumOfCoordsets(), 1,
            'parsePDB failed to parse the last coordinate set')

    def testNameArgument(self):
        """Test outcome of *name* argument."""
        
        path = self.pdb['path']
        name = 'small protein'    
        self.assertEqual(parsePDB(path, name=name).getName(), 
             name, 'parsePDB failed to set user given name')

        name = 1999
        self.assertEqual(parsePDB(path, name=name).getName(), 
             str(name), 'parsePDB failed to set user given non-string name')
            
    def testChainArgument(self):
        """Test outcome of valid and invalid *chain* arguments."""
        
        path = self.pdb['path']
        self.assertRaises(TypeError, parsePDB, path, chain=['A'])
        self.assertRaises(ValueError, parsePDB, path, chain='')
        self.assertIsNone(parsePDB(path, chain='$'))
        self.assertEqual(parsePDB(path, chain='A')
            .getNumOfAtoms(), self.pdb['atoms'],
            'parsePDB failed to parse correct number of atoms when chain is '
            'specified')

    def testSubsetArgument(self):
        """Test outcome of valid and invalid *subset* arguments."""

        path = self.pdb['path']
        self.assertRaises(TypeError, parsePDB, path, subset=['A'])
        self.assertRaises(ValueError, parsePDB, path, subset='')
        self.assertEqual(parsePDB(path, subset='ca').getNumOfAtoms(), 10)
        self.assertEqual(parsePDB(path, subset='bb').getNumOfAtoms(), 49)

    def testAgArgument(self):
        pass
    
    def testBiomolArgument(self):
        pass

    def testSecondaryArgument(self):
        pass

class TestParsePDBHeaderOnly(unittest.TestCase):
    
    def setUp(self):
        self.header = parsePDB('data/proteins_nmr_2k39_models_1to3.pdb', 
                                     header=True, model=0)

    def testHeaderType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')
        
    def testHeaderContent(self):
        self.assertEqual(self.header.get('resolution'), 'NOT APPLICABLE',
            'mistakes in header dictionary content')
        self.assertEqual(self.header.get('classification'), 
            'SIGNALING PROTEIN',
            'mistakes in header dictionary content')
        self.assertEqual(self.header.get('experiment'), 'SOLUTION NMR',
            'mistakes in header dictionary content')
        self.assertEqual(self.header.get('deposition_date'), '25-APR-08',
            'mistakes in header dictionary content')
        self.assertEqual(self.header.get('identifier'), '2K39',
            'mistakes in header dictionary content')
        self.assertEqual(self.header.get('source'), 
            'MOL_ID: 1;  ORGANISM_SCIENTIFIC: XENOPUS LAEVIS;  '
            'ORGANISM_COMMON: AFRICAN CLAWED FROG;  '
            'EXPRESSION_SYSTEM: ESCHERICHIA COLI',
            'mistakes in header dictionary content')
        self.assertEqual(self.header.get('title'), 
            'RECOGNITION DYNAMICS UP TO MICROSECONDS REVEALED FROM RDC  '
            'DERIVED UBIQUITIN ENSEMBLE IN SOLUTION',
            'mistakes in header dictionary content')

    def tearDown(self):
        
        self.header = None


class TestParsePDBHeaderAndAllModels(unittest.TestCase):

    def setUp(self):
        self.atomgroup, self.header = \
            parsePDB('data/proteins_nmr_2k39_models_1to3.pdb', 
                           header=True)

    def testAtomGroupType(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')
        self.assertIsInstance(self.atomgroup, prody.AtomGroup,
            'atom group type is incorrect')
        
    def testAtomGroupContent(self):
        
        self.assertEqual(self.atomgroup.getNumOfAtoms(), 1231,
            'incorrect number of atoms')
        self.assertEqual(self.atomgroup.getNumOfCoordsets(), 3,
            'incorrect number of coordinate sets (models)')

    def tearDown(self):
        
        self.header = None
        self.atomgroup = None


class TestParsePDBAltloc(unittest.TestCase):
    
    def setUp(self):
        
        self.pdbfile = 'data/proteins_altloc_1ejg.pdb'
    
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
    
    def setUp(self):
        """Setup the testing framework."""

        self.pdbs = [PDB_FILES['dssp']]
    
    @unittest.skipIf(prody.which('dssp') is None, 'dssp is not found')
    def testDSSPBridgePartners(self):
        """Check if the DSSP bridge-partners were correctly parsed and 
        assigned."""

        for pdb in self.pdbs:
            prot_ag = parsePDB(pdb['path'], folder=TEMPDIR)
            dssp = execDSSP(pdb['path'], outputdir=TEMPDIR)
            parseDSSP(dssp, prot_ag, parseall=True)
    
            # Map a dssp_resnum to its Residue object.
            dssp_dict = {}

            for chain in prot_ag.select("protein").getHierView():
                for res in chain:
                    dssp_resnum = res.getAttribute("dssp_resnum")[0]
                    dssp_dict[dssp_resnum] = res

            for res in dssp_dict.itervalues():
                bp1 = res.getAttribute("dssp_bp1")[0]
                bp2 = res.getAttribute("dssp_bp2")[0]

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
