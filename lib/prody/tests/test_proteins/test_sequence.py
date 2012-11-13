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

from StringIO import StringIO

from numpy import array, log, zeros, char
from numpy.testing import assert_array_equal, assert_array_almost_equal, dec

from prody.tests.test_datafiles import *

from prody import MSAFile, parseMSA, LOGGER, calcShannonEntropy, calcMutualInfo
from prody import calcMSAOccupancy

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
        
        fasta = StringIO()
        with MSAFile(fasta, 'w', format='fasta') as out:
            for label, seq in MSAFile(pathDatafile('msa_Cys_knot.fasta'), 
                                      split=False):
                out.write(label, seq)
        fasta.seek(0)
        fasta_list = list(MSAFile(fasta))
        self.assertListEqual(FASTA_LIST, fasta_list)        
   
    def testWriteSelex(self):
        
        selex = StringIO()
        with MSAFile(selex, 'w', format='selex') as out:
            for label, seq in MSAFile(pathDatafile('msa_Cys_knot.slx'), 
                                      split=False):
                out.write(label, seq)
        selex.seek(0)
        selex_list = list(MSAFile(selex))
        self.assertListEqual(SELEX_LIST, selex_list) 
        
    def testWriteStockholm(self):
        
        stock = StringIO()
        with MSAFile(stock, 'w', format='stock') as out:
            for label, seq in MSAFile(pathDatafile('msa_Cys_knot.sth'), 
                                      split=False):
                out.write(label, seq)
        stock.seek(0)
        stock_list = list(MSAFile(stock))
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
        result = calcShannonEntropy(msa, omitgaps=True)
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


class TestCalcMutualInfo(TestCase):
  

    def testSixSequences(self):
        
        msa = array([list('ACCA'), 
                     list('ACDA'),
                     list('ACEC'),
                     list('ACGC')])

        expect = array([[0., 0., 0., 0.],
                        [0., 0., 0., 0.],
                        [0., 0., 0., log(2.)],
                        [0., 0., log(2.), 0.],]) 
        result = calcMutualInfo(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


    def testTwenty(self):
        
        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, s] for s in seq])

        expect = log(20.)
        expect = array([[0., expect],
                        [expect, 0.]])
        result = calcMutualInfo(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')
    

    def testTwentyReversed(self):
        
        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, seq[-i-1]] for i, s in enumerate(seq)])

        expect = log(20.)
        expect = array([[0., expect],
                        [expect, 0.]])
        result = calcMutualInfo(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


    def testAmbiguity(self):
        
        msa = array([list('OX'),
                     list('XO')])

        expect = array([[0., log(2.)],
                        [log(2.), 0.]]) 
        result = calcMutualInfo(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


    def testNoAmbiguity(self):
        
        msa = array([list('OX'),
                     list('XO')])

        expect = array([[0., log(2.)],
                        [log(2.), 0.]]) 
        result = calcMutualInfo(msa, ambiquity=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, ambiquity=False, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


    def testAmbiguity2(self):
        
        msa = array([list('AB'),
                     list('BZ')])
        expect = (2 * .25 * log(.25 / .5 / .25) + 
                  4 * .125 * log(.125 / .25 / .25))
        expect = array([[0., expect],
                        [expect, 0.]]) 
        result = calcMutualInfo(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')
        

    def testAmbiguity3(self):
        
        msa = array([list('XX')])

        expect = zeros((2, 2)) 
        result = calcMutualInfo(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


    def testAmbiguity4(self):
        
        msa = array([list('Bb'),
                     list('jJ'),
                     list('Zz'),])

        expect = log((1./12) / (1./6) / (1./6))
        expect = array([[0., expect],
                        [expect, 0.]]) 
        result = calcMutualInfo(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


    def testAmbiguity5(self):

        expect = array([[0., 0.],
                        [0., 0.]]) 

        for seq in ['bx', 'Xb', 'jX', 'Xj', 'xz', 'ZX',
                    'bj', 'jb', 'bz', 'zb', 'jz', 'zj']:
            msa = array([list(seq)])
            result = calcMutualInfo(msa, debug=False)
            assert_array_almost_equal(expect, result, err_msg=seq + ' failed')


    def testAmbiguity6(self):

        expect = zeros((2, 2)) 

        for seq in ['bb', 'jj', 'zz']:
            msa = array([list(seq)])
            result = calcMutualInfo(msa, debug=False)
            assert_array_almost_equal(expect, result, err_msg=seq + ' failed')

    def testAmbiguity7(self):
        
        msa = array([list('bx'),
                     list('xb')])
        expect = (72 * 0.0125 * log(0.0125/0.0250/0.275) + 
                  4 * 0.0250 * log(0.0250/0.275/0.275))
        expect = array([[0., expect],
                        [expect, 0.]]) 
        result = calcMutualInfo(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testInf(self):
        
        msa = zeros((500, 10), '|S1')
        msa.fill('.')
        msa[95,8] = 's'
        msa[95,9] = 'i'
        expect = zeros((10,10))
        expect[8,9] = expect[9,8] = 0.002 * log(500.) + 0.998 * log(1. / 0.998)
        result = calcMutualInfo(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = calcMutualInfo(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


class TestCalcMSAOccupancy(TestCase):
    
    def testResidueOccupancy(self):
        
        assert_array_equal(calcMSAOccupancy(FASTA, 'residue'),
                           char.isalpha(FASTA._msa).sum(0))
    
    def testSequenceOccupancy(self):
        
        assert_array_equal(calcMSAOccupancy(FASTA, 'sequence'),
                           char.isalpha(FASTA._msa).sum(1))
        
