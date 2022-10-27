"""This module contains unit tests for :mod:`prody.database.pfam` module."""

from numpy import tile, array, arange, ones
from numpy.testing import assert_allclose

from prody.tests import unittest
from prody.database.pfam import searchPfam
from prody.database.pfam import fetchPfamMSA
from prody.database.pfam import parsePfamPDBs

from prody.atomic.selection import Selection

import os
import shutil

from prody import LOGGER
LOGGER.verbosity = 'none'

class TestSearchPfam(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.workdir = 'pfam_search_tests'
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        os.chdir(self.workdir)

        self.queries = ['P19491', 'GRIA2_RAT', '6qkcB', '6qkcI', 
                        'VQVLLTTIGAFAAFGLMTIAISTDYWLYTRGLTHSGLWRICCLEGLK'\
                            'RGVCVKINHFAEYLLRVVRASSIFPILSAILLLLGGVCVAASR'\
                            'VYKSKRNIILGAGILFVAAGLSNIIGVIVYISANAGKNHYSYG'\
                            'WSFYFGGLSFILAEVIGVLAVNIYIERSR']

    def testUniprotAccMulti(self):
        """Test the outcome of a simple search scenario using a Uniprot Accession 
        for a multi-domain protein, AMPAR GluA2."""

        a = searchPfam(self.queries[0])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), 
                           ['PF00060', 'PF01094', 'PF10613'],
                           'searchPfam failed to return the right domain family IDs')

    def testUniprotIdMulti(self):
        """Test the outcome of a simple search scenario using a Uniprot ID
        for a multi-domain protein, AMPAR GluA2."""

        a = searchPfam(self.queries[1])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), 
                           ['PF00060', 'PF01094', 'PF10613'],
                           'searchPfam failed to return the right domain family IDs')
        
    def testPdbIdChMulti(self):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID for the same multi-domain protein from specifying chain B."""

        a = searchPfam(self.queries[2])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), 
                           ['PF00060', 'PF01094', 'PF10613'],
                           'searchPfam failed to return the right domain family IDs for AMPAR')
        
    def testPdbIdChSingle(self):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID to get the single domain protein TARP g8 from chain I."""

        a = searchPfam(self.queries[3])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), 
                           ['PF00822'],
                           'searchPfam failed to return the right domain family IDs for TARP')
        
    def testSeqSingle(self):
        """Test the outcome of a simple search scenario using the sequence 
        of the single domain protein TARP g8 from 6qkc chain I."""

        a = searchPfam(self.queries[4])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), 
                           ['PF00822'],
                           'searchPfam failed to return the right domain family IDs for TARP')


class TestFetchPfamMSA(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.query = 'PF00822'

        self.workdir = 'pfam_msa_tests'
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        os.chdir(self.workdir)

    def testDefault(self):
        """Test the outcome of fetching the domain MSA for claudins
        with default parameters."""

        b = fetchPfamMSA(self.query)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_full.sth')
        
        self.assertTrue(os.path.exists(b))


    def testSeed(self):
        """Test the outcome of fetching the domain MSA for claudins
        with the alignment type argument set to seed"""

        b = fetchPfamMSA(self.query, "seed")

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))

        
    def testFormat(self):
        """Test the outcome of fetching the domain MSA for claudins
        with keyword argument format set to fasta."""

        b = fetchPfamMSA(self.query, format="fasta")

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_full.fasta')
        
        self.assertTrue(os.path.exists(b))


    def testFolder(self):
        """Test the outcome of fetching the domain MSA for claudins
        with keyword folder set to a folder that is made especially."""

        folder = "new_folder"
        os.mkdir(folder)
        b = fetchPfamMSA(self.query, folder=folder)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'new_folder/PF00822_full.slx')
        
        self.assertTrue(os.path.exists(b))
    
    @classmethod
    def tearDownClass(self):
        os.chdir('..')
        shutil.rmtree(self.workdir)
            

class TestParsePfamPDBs(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        self.queries = ['PF20446', 'Q57ZF2_TRYB2', 'VAS1_BOVIN']

        self.workdir = 'pfam_pdb_tests'
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)
        os.chdir(self.workdir)

    def testPfamIdDefault(self):
        """Test the outcome of parsing PDBs for a tiny family 
        of ABC class ATPase N-terminal domains (5 members) 
        with the Pfam ID and default parameters."""

        b = parsePfamPDBs(self.queries[0])

        self.assertIsInstance(b, list,
            'fetchPfamMSA failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'fetchPfamMSA failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 5,
            'fetchPfamMSA failed to return a list of length 5')


    def testUniprotDefault(self):
        """Test the outcome of parsing PDBs for a tiny family 
        of ABC class ATPase N-terminal domains (5 members) 
        with the Uniprot long ID and default parameters."""

        b = parsePfamPDBs(self.queries[1])

        self.assertIsInstance(b, list,
            'fetchPfamMSA failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'fetchPfamMSA failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 5,
            'fetchPfamMSA failed to return a list of length 5')

        
    def testMultiDomainDefault(self):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1, 
        which has two domains but few relatives. Default parameters should 
        return Selection objects containing the first domain."""

        b = parsePfamPDBs(self.queries[2])

        self.assertIsInstance(b, list,
            'fetchPfamMSA failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'fetchPfamMSA failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 7,
            'fetchPfamMSA failed to return a list of length 7')
        
        self.assertEqual(b[0].getResnums()[0], 262,
            'fetchPfamMSA failed to return a first Selection with first resnum 262')        

    def testMultiDomainStart1(self):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1, 
        which has two domains but few relatives. Using start=1 should be like default and 
        return Selection objects containing the first domain."""

        b = parsePfamPDBs(self.queries[2], start=1)

        self.assertIsInstance(b, list,
            'fetchPfamMSA failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'fetchPfamMSA failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 7,
            'fetchPfamMSA failed to return a list of length 7')
        
        self.assertEqual(b[0].getResnums()[0], 262,
            'fetchPfamMSA failed to return a first Selection with first resnum 262')  
        
    def testMultiDomainStart2(self):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1, 
        which has two domains but few relatives. Setting start to 418 should 
        return Selection objects containing the second domain."""

        b = parsePfamPDBs(self.queries[2], start=418)

        self.assertIsInstance(b, list,
            'fetchPfamMSA failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'fetchPfamMSA failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 23,
            'fetchPfamMSA failed to return a list of length 23')
        
        self.assertEqual(b[0].getResnums()[0], 418,
            'fetchPfamMSA failed to return a first Selection with first resnum 418')  
        
    @classmethod
    def tearDownClass(self):
        os.chdir('..')
        shutil.rmtree(self.workdir)

