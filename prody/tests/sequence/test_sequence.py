from prody.tests import TestCase

from prody.tests.datafiles import *

from prody import LOGGER, Sequence

LOGGER.verbosity = None


class TestSequence(TestCase):


    def testStringConversion(self):

        seq = 'SOME-PROTEIN-SEQUENCE'
        self.assertEqual(str(Sequence(seq)), seq)
