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

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from unittest import TestCase

from numpy.testing import assert_equal

from prody.tests.test_datafiles import *

from prody import MSAFile, parseMSA

class TestMSAFile(TestCase):
    
    def testMSAFile(self):
        
        fasta = list(MSAFile(pathDatafile('msa_Cys_knot.fasta'))) 
        
        self.assertListEqual(fasta, 
                             list(MSAFile(pathDatafile('msa_Cys_knot.sth'))))
        
        self.assertListEqual(fasta,
                             list(MSAFile(pathDatafile('msa_Cys_knot.sth'))))


class TestParseMSA(TestCase):
    
    def testArray(self):
        
        #fasta = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
        selex = parseMSA(pathDatafile('msa_Cys_knot.slx'))
        stockholm = parseMSA(pathDatafile('msa_Cys_knot.sth'))
        
        #asser_equal(fasta._getArray(), selex._getArray())
        assert_equal(selex._getArray(), stockholm._getArray())

    def testIterator(self):
        
        #fasta = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
        selex = parseMSA(pathDatafile('msa_Cys_knot.slx'))
        stockholm = parseMSA(pathDatafile('msa_Cys_knot.sth'))
        
        #asser_equal(fasta._getArray(), selex._getArray())
        self.assertListEqual(list(selex), list(stockholm))
