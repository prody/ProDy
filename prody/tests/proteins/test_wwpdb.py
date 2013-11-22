"""This module contains unit tests for :mod:`~prody.proteins`."""

import os

from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'


class TestFTP(unittest.TestCase):

    def setUp(self):

        self.pdb = ['1ubi', '1aar', 'arg', 1234]
        self.fns = []
        self.len = [683, 1218, None, None]
        self.fetch = fetchPDBviaFTP
        self.protocol = 'FTP'


    @dec.slow
    def testCompressed(self):

        self.fns = self.fetch(*self.pdb, folder=TEMPDIR)

        for fn, pdb, n_atoms in zip(self.fns, self.pdb, self.len):
            if fn is None:
                continue
            self.assertTrue(os.path.isfile(fn),
                'failed to fetch PDB file via ' + self.protocol)
            atoms = parsePDB(fn)
            self.assertIsNotNone(atoms, 'PDB file ({}) empty, download via {}'
                                        .format(fn, self.protocol))
            self.assertEqual(n_atoms, len(atoms))


    @dec.slow
    def testDecompressed(self):

        self.fns = self.fetch(*self.pdb, folder=TEMPDIR, compressed=False)

        for fn, pdb, n_atoms in zip(self.fns, self.pdb, self.len):
            if fn is None:
                continue
            self.assertTrue(os.path.isfile(fn),
                'failed to fetch PDB file via ' + self.protocol)
            atoms = parsePDB(fn)
            self.assertIsNotNone(atoms, 'PDB file ({}) empty, download via {}'
                                        .format(fn, self.protocol))
            self.assertEqual(n_atoms, len(atoms))


    def tearDown(self):


        for fn in self.fns:
            if fn is None:
                continue
            try:
                #os.remove(fn)
                pass
            except:
                pass

class TestHTTP(TestFTP):

    def setUp(self):

        TestFTP.setUp(self)
        self.fetch = fetchPDBviaHTTP
        self.protocol = 'HTTP'

