"""This module contains unit tests for :mod:`~prody.ensemble`."""

from os.path import join
from prody.tests import TestCase

from numpy.testing import assert_equal, assert_allclose

from prody import DCDFile, writeDCD, parseDCD

from prody.tests import TEMPDIR
from prody.tests.ensemble import ALLATOMS, ENSEMBLE, RTOL, ATOL, DCD

class TestDCDFile(TestCase):

    def setUp(self):

        self.dcd = join(TEMPDIR, 'temp.dcd')

    def testWriteDCD(self):
        dcd = writeDCD(self.dcd, ALLATOMS)
        self.assertEqual(dcd, self.dcd, 'failed to write DCD file')

    def testParseDCD(self):
        e = parseDCD(writeDCD(self.dcd, ALLATOMS))
        assert_equal(e._getCoordsets(), DCD._getCoordsets(),
                     err_msg='failed to parse DCD file correctly')

    def testWrite(self):
        dcd = DCDFile(self.dcd, 'w')
        dcd.write(ENSEMBLE.getCoordsets())
        dcd.close()
        e = parseDCD(self.dcd)
        assert_allclose(e._getCoordsets(), ENSEMBLE._getCoordsets(),
                        rtol=RTOL, atol=ATOL,
                        err_msg='failed to parse DCD file correctly')

    def testWriteModeAppend(self):
        dcd = DCDFile(writeDCD(self.dcd, ENSEMBLE), 'a')
        dcd.write(ENSEMBLE.getCoordsets())
        dcd.close()
        e = parseDCD(self.dcd)
        n_csets = len(ENSEMBLE)
        coordsets = e._getCoordsets()
        assert_equal(coordsets, coordsets,
                     'failed to parse DCD file correctly')
        assert_allclose(coordsets[:n_csets], ENSEMBLE._getCoordsets(),
                        rtol=RTOL, atol=ATOL,
                        err_msg='failed to parse DCD file correctly')
