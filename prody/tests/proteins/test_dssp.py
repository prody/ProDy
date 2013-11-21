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

from numpy.testing import *

from prody import *
from prody import LOGGER
from prody.utilities import which
from prody.tests import TEMPDIR, unittest
from prody.tests.datafiles import *

LOGGER.verbosity = 'none'


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
            dssp = execDSSP(pathDatafile(pdb['file']), outputdir=TEMPDIR,
                            stderr=False)
            parseDSSP(dssp, prot_ag, parseall=True)

            # Map a dssp_resnum to its Residue object.
            dssp_dict = {}

            for chain in prot_ag.select("protein").getHierView():
                for res in chain:
                    dssp_resnum = res.getData("dssp_resnum")[0]
                    dssp_dict[dssp_resnum] = res

            for res in dssp_dict.values():
                bp1 = res.getData("dssp_bp1")[0]
                bp2 = res.getData("dssp_bp2")[0]

                if bp1 != 0:
                    msg_ = "BP1 (dssp_resnum: %d) of %s is missing" % \
                        (bp1, str(res))
                    self.assertIn(bp1, dssp_dict, msg=msg_)

                if bp2 != 0:
                    msg_ = "BP2 (dssp_resnum: %d) of %s is missing" % \
                        (bp2, str(res))
                    self.assertIn(bp2, dssp_dict, msg=msg_)
