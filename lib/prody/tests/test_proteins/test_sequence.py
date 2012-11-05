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

from numpy import array, log, zeros
from numpy.testing import assert_equal, assert_array_almost_equal, dec

from prody.tests.test_datafiles import *

from prody import MSAFile, parseMSA, LOGGER, calcShannonEntropy

LOGGER.verbosity = None

class TestMSAFile(TestCase):
    
    def testMSAFile(self):
        
        fasta = list(MSAFile(pathDatafile('msa_Cys_knot.fasta'))) 
        
        self.assertListEqual(fasta, 
                             list(MSAFile(pathDatafile('msa_Cys_knot.sth'))))
        
        self.assertListEqual(fasta,
                             list(MSAFile(pathDatafile('msa_Cys_knot.sth'))))


class TestParseMSA(TestCase):
    
    def testArray(self):
        
        fasta = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
        selex = parseMSA(pathDatafile('msa_Cys_knot.slx'))
        stockholm = parseMSA(pathDatafile('msa_Cys_knot.sth'))
        
        assert_equal(fasta._getArray(), selex._getArray())
        assert_equal(selex._getArray(), stockholm._getArray())

    def testIterator(self):
        
        fasta = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
        selex = parseMSA(pathDatafile('msa_Cys_knot.slx'))
        stockholm = parseMSA(pathDatafile('msa_Cys_knot.sth'))
        
        self.assertListEqual(list(fasta), list(selex))
        self.assertListEqual(list(fasta), list(stockholm))

    def testMapping(self):
        
        fasta = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
        selex = parseMSA(pathDatafile('msa_Cys_knot.slx'))
        stockholm = parseMSA(pathDatafile('msa_Cys_knot.sth'))
        
        self.assertDictEqual(fasta._mapping, selex._mapping)
        self.assertDictEqual(fasta._mapping, stockholm._mapping)


class TestCalcShannonEntropy(TestCase):


    def testSixSequences(self):
        
        msa = array([list('AAAAaaaaAAAAaaaa'), 
                     list('AAACaaacAAACaaac'),
                     list('AACDaacdAACDaacd'),
                     list('ACCEacceacceACCE'),
                     list('ACDFacdfacdfACDF'),
                     list('ACDGacdgacdgACDG')])

        expect = -log(1. / array([1, 2, 3, 6] * 4)) 
        result = calcShannonEntropy(msa)
        assert_array_almost_equal(expect, result)

    def testTwenty(self):
        
        msa = array([[char] for char in 'ACDEFGHIKLMNPQRSTVWY'])

        expect = -log(1. / 20)
        result = calcShannonEntropy(msa)
        assert_array_almost_equal(expect, result)

    
    def testSmallProbability(self):
        
        msa = zeros((1000000,1), '|S1')
        msa[0] = 'A'
        msa[1:] = 'C'
        expect = array([1., 999999.]) / 1000000
        expect = - (expect * log(expect)).sum() 
        result = calcShannonEntropy(msa)
        assert_array_almost_equal(expect, result)


    def testAmbiguous(self):
        
        msa = array([list('bjzxBJZX'),
                     list('bjzxBJZX'),])

        expect = -log(1. / array([2, 2, 2, 20] * 2)) 
        result = calcShannonEntropy(msa)
        assert_array_almost_equal(expect, result)

    def testGapDividend(self):
        
        msa = array([list('AAAA'), 
                     list('AAAC'),
                     list('AACD'),
                     list('ACCE'),
                     list('ACDF'),
                     list('ACDG'),
                     list('----')])

        expect = -log(1. / array([1, 2, 3, 6])) 
        result = calcShannonEntropy(msa, dividend=True)
        assert_array_almost_equal(expect, result)
        
"""
    def testSixSequences3(self):
        
        msa = array([list('AAAA'), 
                     list('AAAB'),
                     list('AABC'),
                     list('ABBD'),
                     list('ABCE'),
                     list('ABCF')])

        expect = -log(1. / array([1, 2, 3, 6])) 
        result = calcInfoEntropy(msa)
        assert_array_almost_equal(expect, result)
"""
