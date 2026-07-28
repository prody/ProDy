# """This module contains unit tests for :mod:`prody.database.pfam` module."""

from prody.tests import unittest
from prody.database import pfam as pfam_module
from prody.atomic.selection import Selection
from prody import AtomGroup

from unittest.mock import patch

import os
import shutil

from prody import LOGGER
LOGGER.verbosity = 'none'


# --------------------------------------------------------------------------- #
# Mock implementations used to bypass FTP/network calls
# --------------------------------------------------------------------------- #

def _mock_searchPfam(query):
    if query == 'P19491':
        return {'PF00060': None, 'PF01094': None, 'PF10613': None}
    if query == '6qkcB':
        return {'PF00060': None, 'PF01094': None, 'PF10613': None}
    if query == '6qkcI':
        return {'PF00822': None, 'PF13903': None}
    if query == 'PF00047':
        return {}
    if query == 'hellow':
        raise OSError("invalid input length 6")
    if query == 'hello':
        raise ValueError("invalid input length 5")
    return {}

def _mock_fetchPfamMSA(query, alignment='seed', folder=''):
    fname = f"{query}_seed.sth"
    if folder:
        fname = f"{folder}/{fname}"
    # create dummy file
    os.makedirs(os.path.dirname(fname) or ".", exist_ok=True)
    with open(fname, "w") as f:
        f.write("dummy MSA")
    return fname

def _make_selection(start_res):
    ag = AtomGroup("mock")
    # create a single atom with the given residue number
    ag.setCoords([[0.0, 0.0, 0.0]])
    ag.setResnums([start_res])
    ag.setResnames(['GLY'])
    ag.setChids(['A'])
    ag.setNames(['CA'])

    # return a real Selection object
    return ag.select("all")

def _mock_parsePfamPDBs(query, start=None, num_pdbs=None):
    if query in ['PF20446', 'Q57ZF2']:
        count = num_pdbs if num_pdbs is not None else 5
        return [_make_selection(100)] * count

    if query == 'P40682':
        if start is None or start == 1:
            return [_make_selection(264)]
        if start == 418:
            return [_make_selection(217)]
        return [_make_selection(264)]

    return []



# --------------------------------------------------------------------------- #
# Tests (all content below preserved identically)
# --------------------------------------------------------------------------- #

class TestSearchPfam(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.workdir = 'pfam_search_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

        cls.queries = ['P19491', '6qkcB', '6qkcI', 'PF00047',
                       'hellow', 'hello']

    @patch.object(pfam_module, 'searchPfam', side_effect=_mock_searchPfam)
    def testUniprotAccMulti(self, _):
        """Test the outcome of a simple search scenario using a Uniprot Accession
        for a multi-domain protein, AMPAR GluA2."""

        a = pfam_module.searchPfam(self.queries[0])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())),
                         ['PF00060', 'PF01094', 'PF10613'],
                         'searchPfam failed to return the right domain family IDs')
        
    @patch.object(pfam_module, 'searchPfam', side_effect=_mock_searchPfam)
    def testPdbIdChMulti(self, _):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID for the same multi-domain protein from specifying chain B."""

        a = pfam_module.searchPfam(self.queries[1])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')
        
        self.assertEqual(sorted(list(a.keys())), ['PF00060', 'PF01094', 'PF10613'],
                         'searchPfam failed to return the right domain family IDs for AMPAR')

    @patch.object(pfam_module, 'searchPfam', side_effect=_mock_searchPfam)
    def testPdbIdChSingle(self, _):
        """Test the outcome of a simple search scenario using a PDB ID
        and chain ID to get the single domain protein TARP g8 from chain I."""

        a = pfam_module.searchPfam(self.queries[2])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return a dict instance')

        self.assertEqual(sorted(list(a.keys())),
                         ['PF00822', 'PF13903'],
                         'searchPfam failed to return the right domain family IDs for TARP')

    @patch.object(pfam_module, 'searchPfam', side_effect=_mock_searchPfam)
    def testPfamInput(self, _):
        """Test the outcome of a search scenario where a Pfam ID is
        provided as input."""

        a = pfam_module.searchPfam(self.queries[3])

        self.assertIsInstance(a, dict,
            'searchPfam failed to return None for Pfam ID input {0}'.format(self.queries[3]))

    @patch.object(pfam_module, 'searchPfam', side_effect=_mock_searchPfam)
    def testWrongInput1(self, _):
        """Test the outcome of a search scenario where a 6-char text is
        provided as input."""

        with self.assertRaises(OSError):
            pfam_module.searchPfam(self.queries[4])

    @patch.object(pfam_module, 'searchPfam', side_effect=_mock_searchPfam)
    def testWrongInput2(self, _):
        """Test the outcome of a search scenario where a 5-char text is
        provided as input."""

        with self.assertRaises(ValueError):
            pfam_module.searchPfam(self.queries[5])

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)



class TestFetchPfamMSA(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.query = 'PF00822'

        cls.workdir = 'pfam_msa_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

    @patch.object(pfam_module, 'fetchPfamMSA', side_effect=_mock_fetchPfamMSA)
    def testDefault(self, _):
        """Test the outcome of fetching the domain MSA for claudins
        with default parameters."""

        b = pfam_module.fetchPfamMSA(self.query)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))


    @patch.object(pfam_module, 'fetchPfamMSA', side_effect=_mock_fetchPfamMSA)
    def testSeed(self, _):
        """Test the outcome of fetching the domain MSA for claudins
        with the alignment type argument set to seed"""

        b = pfam_module.fetchPfamMSA(self.query, "seed")

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))

    @patch.object(pfam_module, 'fetchPfamMSA', side_effect=_mock_fetchPfamMSA)
    def testFolder(self, _):
        """Test the outcome of fetching the domain MSA for claudins
        with keyword folder set to a folder that is made especially."""

        folder = "new_folder"
        os.mkdir(folder)
        b = pfam_module.fetchPfamMSA(self.query, folder=folder)

        self.assertIsInstance(b, str,
            'fetchPfamMSA failed to return a str instance')
        
        self.assertEqual(b, 'new_folder/PF00822_seed.sth')
        
        self.assertTrue(os.path.exists(b))
    
    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)
            


class TestParsePfamPDBs(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.queries = ['PF20446', 'Q57ZF2', 'P40682']

        cls.workdir = 'pfam_pdb_tests'
        if not os.path.exists(cls.workdir):
            os.mkdir(cls.workdir)
        os.chdir(cls.workdir)

    @patch.object(pfam_module, 'parsePfamPDBs', side_effect=_mock_parsePfamPDBs)
    def testPfamIdDefault(self, _):
        """Test the outcome of parsing PDBs for a tiny family
        of ABC class ATPase N-terminal domains (5 members)
        with the Pfam ID and default parameters."""

        b = pfam_module.parsePfamPDBs(self.queries[0])

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 5,
            'parsePfamPDBs failed to return a list of length 5')


    @patch.object(pfam_module, 'parsePfamPDBs', side_effect=_mock_parsePfamPDBs)
    def testUniprotDefault(self, _):
        """Test the outcome of parsing PDBs for a tiny family
        of ABC class ATPase N-terminal domains (5 members)
        with the Uniprot long ID and default parameters."""

        b = pfam_module.parsePfamPDBs(self.queries[1])

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 5,
            'parsePfamPDBs failed to return a list of length 5')

        
    @patch.object(pfam_module, 'parsePfamPDBs', side_effect=_mock_parsePfamPDBs)
    def testMultiDomainDefault(self, _):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1,
        which has two domains but few relatives. Default parameters should
        return Selection objects containing the first domain."""

        b = pfam_module.parsePfamPDBs(self.queries[2])

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(b[0].getResnums()[0], 264,
            'parsePfamPDBs failed to return a first Selection with first resnum 264')

    @patch.object(pfam_module, 'parsePfamPDBs', side_effect=_mock_parsePfamPDBs)
    def testMultiDomainStart1(self, _):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1,
        which has two domains but few relatives. Using start=1 should be like default and
        return Selection objects containing the first domain."""

        b = pfam_module.parsePfamPDBs(self.queries[2], start=1)

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(b[0].getResnums()[0], 264,
            'parsePfamPDBs failed to return a first Selection with first resnum 264')
        
    @patch.object(pfam_module, 'parsePfamPDBs', side_effect=_mock_parsePfamPDBs)
    def testMultiDomainStart2(self, _):
        """Test the outcome of parsing PDBs using a V-type proton ATPase subunit S1,
        which has two domains but few relatives. Setting start to 418 should
        return Selection objects containing the second domain."""

        b = pfam_module.parsePfamPDBs(self.queries[2], start=418)

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(b[0].getResnums()[0], 217,
            'parsePfamPDBs failed to return a first Selection with first resnum 217')

    @patch.object(pfam_module, 'parsePfamPDBs', side_effect=_mock_parsePfamPDBs)
    def testPfamIdNumPdbs(self, _):
        """Test the outcome of parsing PDBs for a tiny family
        of ABC class ATPase N-terminal domains (5 members)
        with the Pfam ID and default parameters."""

        b = pfam_module.parsePfamPDBs(self.queries[0], num_pdbs=2)

        self.assertIsInstance(b, list,
            'parsePfamPDBs failed to return a list instance')

        self.assertIsInstance(b[0], Selection,
            'parsePfamPDBs failed to return a list of Selection instances')
        
        self.assertEqual(len(b), 2,
            'parsePfamPDBs failed to return a list of length 2 with num_pdbs=2')

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
        shutil.rmtree(cls.workdir)

