__author__ = 'Ahmet Bakan, Anindita Dutta'

from prody.tests import TestCase
import os
from os.path import join

from numpy import array, log, zeros, char
from numpy.testing import assert_array_equal

from prody.tests.datafiles import *
from prody.tests import TEMPDIR
from prody import MSA, MSAFile, parseMSA, LOGGER, writeMSA
from prody.utilities import createStringIO, importDec
dec = importDec()

LOGGER.verbosity = None

FASTA = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
SELEX = parseMSA(pathDatafile('msa_Cys_knot.slx'))
STOCK = parseMSA(pathDatafile('msa_Cys_knot.sth'))
FASTA_LIST = list(MSAFile(pathDatafile('msa_Cys_knot.fasta')))
SELEX_LIST = list(MSAFile(pathDatafile('msa_Cys_knot.sth')))
STOCK_LIST = list(MSAFile(pathDatafile('msa_Cys_knot.sth')))

class TestMSAFile(TestCase):

    def testMSAFile(self):

        self.assertListEqual(FASTA_LIST, SELEX_LIST)
        self.assertListEqual(FASTA_LIST, STOCK_LIST)

    def testWriteFasta(self):

        fasta = createStringIO()
        with MSAFile(fasta, 'w', format='fasta') as out:
            for seq in MSAFile(pathDatafile('msa_Cys_knot.fasta')):
                out.write(seq)
        fasta.seek(0)
        fasta_list = list(MSAFile(fasta, format='fasta'))
        self.assertListEqual(FASTA_LIST, fasta_list)

    def testWriteSelex(self):

        selex = createStringIO()
        with MSAFile(selex, 'w', format='selex') as out:
            for seq in MSAFile(pathDatafile('msa_Cys_knot.slx')):
                out.write(seq)
        selex.seek(0)
        selex_list = list(MSAFile(selex, format='selex'))
        self.assertListEqual(SELEX_LIST, selex_list)

    def testWriteStockholm(self):

        stock = createStringIO()
        with MSAFile(stock, 'w', format='stockholm') as out:
            for seq in MSAFile(pathDatafile('msa_Cys_knot.sth'), split=False):
                out.write(seq)
        stock.seek(0)
        stock_list = list(MSAFile(stock, format='stockholm'))
        self.assertListEqual(STOCK_LIST, stock_list)

class TestParseMSA(TestCase):

    def testArray(self):

        assert_array_equal(FASTA._getArray(), SELEX._getArray())
        assert_array_equal(SELEX._getArray(), STOCK._getArray())

    def testIterator(self):

        self.assertListEqual(list(FASTA), list(SELEX))
        self.assertListEqual(list(FASTA), list(STOCK))

    def testMapping(self):

        self.assertDictEqual(FASTA._mapping, SELEX._mapping)
        self.assertDictEqual(FASTA._mapping, STOCK._mapping)

class TestWriteMSA(TestCase):

    def testSelex(self):
        filename = writeMSA(join(TEMPDIR, 'test.slx'), SELEX)
        selex = parseMSA(pathDatafile(filename))
        self.assertListEqual(list(SELEX), list(selex))
        if os.path.isfile(filename):
            os.remove(filename)

    def testFasta(self):
        filename = writeMSA(join(TEMPDIR, 'test.fasta.gz'), FASTA)
        fasta = list(MSAFile(pathDatafile(filename)))
        self.assertListEqual(list(FASTA), list(fasta))
        if os.path.isfile(filename):
            os.remove(filename)
