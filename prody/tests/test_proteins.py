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

"""This module contains unit tests for :mod:`~prody.proteins` module."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import os
import os.path
import unittest
import prody
from prody.tests import TEMPDIR

prody.changeVerbosity('none')

class TestFetchPDB(unittest.TestCase):
    
    def setUp(self):
        
        self.filenames = []
    
    def testFetchingSingleFile(self):
        fn = prody.fetchPDB('1p38', copy=True)
        self.assertTrue(os.path.isfile(fn), 
                        'fetching a single PDB file failed')
        self.filenames.append(fn)
        
    def testFetchingMultipleFiles(self): 
        
        fns = prody.fetchPDB(['1p38', '1r39'], copy=True)
        self.assertIsInstance(fns, list, 
            'failed to return a list of filenames')
        self.assertTrue(all([os.path.isfile(fn) for fn in fns]),
            'fetching multiple PDB files failed')
        self.filenames.extend(fns)
        
    def testDecompressing(self):
        
        fns = prody.fetchPDB(['1p38', '1r39'], compressed=False)
        self.assertTrue(all([os.path.isfile(fn) for fn in fns]),
            'fetching decompressed PDB files failed')
        self.assertTrue(all([os.path.splitext(fn)[1] != '.gz' for fn in fns]),
            'decompressing PDB files failed')
        self.filenames.extend(fns)

    def testInvalidPDBIdentifier(self):
        
        self.assertIsNone(prody.fetchPDB('XXXXX'),
            'failed to return None for an invalid PDB identifier')
            
        self.assertFalse(all(prody.fetchPDB(['XXXXX', '654654', '-/-*/+', ''])),
            'failed to return None for invalid PDB identifiers')
        

    def tearDown(self):
        
        for fn in self.filenames:
            if os.path.isfile(fn):
                os.remove(fn)
        

class TestParsePDBHeaderOnly(unittest.TestCase): 
    
    def setUp(self):
        self.header = prody.parsePDB('data/proteins_nmr_2k39_models_1to3.pdb', 
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
            prody.parsePDB('data/proteins_nmr_2k39_models_1to3.pdb', 
                           header=True)

    def testReturnTypes(self):
        self.assertIsInstance(self.header, dict,
            'header type is incorrect')
        self.assertIsInstance(self.atomgroup, prody.AtomGroup,
            'atom group type is incorrect')
        
    def testAtomGroupContent(self):
        self.assertEqual(len(self.atomgroup), 1231,
            'len() function reports incorrect number of atoms')
        self.assertEqual(self.atomgroup.getNumOfAtoms(), 1231,
            'incorrect number of atoms')
        self.assertEqual(self.atomgroup.getNumOfCoordsets(), 3,
            'incorrect number of coordinate sets (models)')

    def tearDown(self):
        
        self.header = None
        self.atomgroup = None


class TestParsePDBSubset(unittest.TestCase):
    pass

class TestParsePDBChain(unittest.TestCase):
    pass

class TestParsePDBAltloc(unittest.TestCase):
    
    def setUp(self):
        
        self.pdbfile = 'data/proteins_altloc_1ejg.pdb'
    
    def testAltlocNone(self):
        
        self.assertEqual(len(prody.parsePDB(self.pdbfile)), 637,
            'failed to parse unspecified alternate locations correctly')
    
    def testAltlocA(self):
        
        self.assertEqual(len(prody.parsePDB(self.pdbfile, altloc='A')), 637,
            'failed to parse alternate locations A correctly')
        
    def testAltlocB(self):
        
        self.assertEqual(len(prody.parsePDB(self.pdbfile, altloc='B')), 634,
            'failed to parse alternate locations B correctly')

    def testAltlocC(self):
        
        self.assertEqual(len(prody.parsePDB(self.pdbfile, altloc='C')), 496,
            'failed to parse alternate locations C correctly')


class TestParsePDBName(unittest.TestCase):
    
    def setUp(self):
        
        self.name = 'small protein'
        self.pdbfile = 'data/proteins_small_3nir.pdb'
        self.filename = os.path.splitext(os.path.split(self.pdbfile)[1])[0]
    
    def testDefaultName(self):
        
        self.assertEqual(prody.parsePDB(self.pdbfile).getName(), self.filename,
            'failed to set AtomGroup name based on filename')

    def testUserGivenName(self):
        
        self.assertEqual(prody.parsePDB(self.pdbfile, name=self.name).getName(), 
                         self.name,
            'failed to set user given name as the AtomGroup name')


class TestParsePDBBiomolecule(unittest.TestCase):
    pass

class TestParsePDBSecondary(unittest.TestCase):
    pass

class TestParsePSF(unittest.TestCase):
    pass

class TestParsePSFandPDB(unittest.TestCase):
    pass
 
# Kian's parsing snippet, modify as you see fit.
class TestPDBParsing(unittest.TestCase):
    def setUp(self):
        """Setup the testing framework.

        """

        self.pdb_ids = ["1r19"]

        return

    def test_pdb_parsing(self):
        """Check if each of the PDB id's specified in self.pdb_ids is
        parsed "correctly".

        Note:
            Given the discrepencies in a number of PDB files, defining a
            "correctly" parsed PDB file is non-trivial.
        """

        from prody.proteins import parsePDB, execDSSP, parseDSSP

        for pdb_id in self.pdb_ids:
            prot_ag = parsePDB(pdb_id, folder=TEMPDIR)
            dssp = execDSSP(pdb_id, outputdir=TEMPDIR)
            parseDSSP(dssp, prot_ag, parseall=True)

        return

    def test_dssp_bridge_partners(self):
        """Check if the DSSP bridge-partners were correctly parsed and assigned.

        """

        from prody.proteins import parsePDB, execDSSP, parseDSSP

        for pdb_id in self.pdb_ids:
            prot_ag = parsePDB(pdb_id, folder=TEMPDIR)
            dssp = execDSSP(pdb_id, outputdir=TEMPDIR)
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
                msg_ = "BP1 (dssp_resnum: %d) of %s is missing" % (bp1, str(res))
                self.assertTrue(dssp_dict.has_key(bp1), msg=msg_)

            if bp2 != 0:
                msg_ = "BP2 (dssp_resnum: %d) of %s is missing" % (bp2, str(res))
                self.assertTrue(dssp_dict.has_key(bp2), msg=msg_)

        return


if __name__ == '__main__':
    unittest.main()
