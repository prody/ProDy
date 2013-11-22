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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from prody.tests import TestCase

from numpy import array, log, zeros, char, ones, fromfile
from numpy.testing import assert_array_equal, assert_array_almost_equal

from prody.tests.datafiles import *

from prody import LOGGER, calcShannonEntropy, buildMutinfoMatrix, parseMSA
from prody import calcMSAOccupancy, buildSeqidMatrix, uniqueSequences
from prody import buildOMESMatrix, buildSCAMatrix, calcMeff
from prody import buildDirectInfoMatrix

LOGGER.verbosity = None

FASTA = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
FASTA_ALPHA = char.isalpha(FASTA._msa)
FASTA_UPPER = char.upper(FASTA._msa)

FASTA_NUMBER, FASTA_LENGTH = FASTA_ALPHA.shape
FASTA_EYE = zeros((FASTA_NUMBER, FASTA_NUMBER))
for i in range(FASTA_NUMBER):
    FASTA_EYE[i, i] = 1
    for j in range(i + 1, FASTA_NUMBER):
        score = 0.0
        ncols = 0
        for k in range(FASTA_LENGTH):
            if FASTA_ALPHA[i, k] or FASTA_ALPHA[j, k]:
                if FASTA_UPPER[i, k] == FASTA_UPPER[j, k]:
                    score += 1
                ncols += 1
        FASTA_EYE[i, j] = FASTA_EYE[j, i] = score / ncols


class TestCalcShannonEntropy(TestCase):

    def testSixSequences(self):

        msa = array([list('AAAAaaaaAAAAaaaa'),
                     list('AAACaaacAAACaaac'),
                     list('AACDaacdAACDaacd'),
                     list('ACCEacceacceACCE'),
                     list('ACDFacdfacdfACDF'),
                     list('ACDGacdgacdgACDG')], dtype='|S1')

        expect = -log(1. / array([1, 2, 3, 6] * 4))
        result = calcShannonEntropy(msa)
        assert_array_almost_equal(expect, result)

    def testTwenty(self):

        msa = array([[char] for char in 'ACDEFGHIKLMNPQRSTVWY'], dtype='|S1')

        expect = -log(1. / 20)
        result = calcShannonEntropy(msa)
        assert_array_almost_equal(expect, result)

    def testSmallProbability(self):

        msa = zeros((1000000, 1), '|S1')
        msa[0] = 'A'
        msa[1:] = 'C'
        expect = array([1., 999999.]) / 1000000
        expect = - (expect * log(expect)).sum()
        result = calcShannonEntropy(msa)
        assert_array_almost_equal(expect, result)

    def testAmbiguous(self):

        msa = array([list('bjzxBJZX'),
                     list('bjzxBJZX'), ], dtype='|S1')

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
                     list('----')], dtype='|S1')

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
                     list('ACGC')], dtype='|S1')

        expect = array([[0., 0., 0., 0.],
                        [0., 0., 0., 0.],
                        [0., 0., 0., log(2.)],
                        [0., 0., log(2.), 0.], ])
        result = buildMutinfoMatrix(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testTwenty(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, s] for s in seq], dtype='|S1')

        expect = log(20.)
        expect = array([[0., expect],
                        [expect, 0.]])
        result = buildMutinfoMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testTwentyReversed(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, seq[-i-1]] for i, s in enumerate(seq)], dtype='|S1')

        expect = log(20.)
        expect = array([[0., expect],
                        [expect, 0.]])
        result = buildMutinfoMatrix(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity(self):

        msa = array([list('OX'),
                     list('XO')], dtype='|S1')

        expect = array([[0., log(2.)],
                        [log(2.), 0.]])
        result = buildMutinfoMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testNoAmbiguity(self):

        msa = array([list('OX'),
                     list('XO')], dtype='|S1')

        expect = array([[0., log(2.)],
                        [log(2.), 0.]])
        result = buildMutinfoMatrix(msa, ambiquity=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, ambiquity=False, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity2(self):

        msa = array([list('AB'),
                     list('BZ')], dtype='|S1')
        expect = (2 * .25 * log(.25 / .5 / .25) +
                  4 * .125 * log(.125 / .25 / .25))
        expect = array([[0., expect],
                        [expect, 0.]])
        result = buildMutinfoMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity3(self):

        msa = array([list('XX')], dtype='|S1')

        expect = zeros((2, 2))
        result = buildMutinfoMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity4(self):

        msa = array([list('Bb'),
                     list('jJ'),
                     list('Zz'), ], dtype='|S1')

        expect = log((1./12) / (1./6) / (1./6))
        expect = array([[0., expect],
                        [expect, 0.]])
        result = buildMutinfoMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity5(self):

        expect = array([[0., 0.],
                        [0., 0.]])

        for seq in ['bx', 'Xb', 'jX', 'Xj', 'xz', 'ZX',
                    'bj', 'jb', 'bz', 'zb', 'jz', 'zj']:
            msa = array([list(seq)], dtype='|S1')
            result = buildMutinfoMatrix(msa, debug=False)
            assert_array_almost_equal(expect, result, err_msg=seq + ' failed')

    def testAmbiguity6(self):

        expect = zeros((2, 2))

        for seq in ['bb', 'jj', 'zz']:
            msa = array([list(seq)], dtype='|S1')
            result = buildMutinfoMatrix(msa, debug=False)
            assert_array_almost_equal(expect, result, err_msg=seq + ' failed')

    def testAmbiguity7(self):

        msa = array([list('bx'),
                     list('xb')], dtype='|S1')
        expect = (72 * 0.0125 * log(0.0125/0.0250/0.275) +
                  4 * 0.0250 * log(0.0250/0.275/0.275))
        expect = array([[0., expect],
                        [expect, 0.]])
        result = buildMutinfoMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testInf(self):

        msa = zeros((500, 10), '|S1')
        msa.fill('.')
        msa[95, 8] = 's'
        msa[95, 9] = 'i'
        expect = zeros((10, 10))
        expect[8, 9] = expect[9, 8] = 0.002 * log(500.) + .998 * log(1. / .998)
        result = buildMutinfoMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildMutinfoMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testNorm(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, seq[-i-1]] for i, s in enumerate(seq)], dtype='|S1')

        expect = 1.
        expect = array([[0., expect],
                        [expect, 0.]])
        result = buildMutinfoMatrix(msa, norm=True)
        assert_array_almost_equal(expect, result, err_msg='norm failed')

    def testNorm2(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, 'O' if i % 2 else 'U'] for i, s in enumerate(seq)],
                    dtype='|S1')

        expect = log(1./20. / (1./20. * 1./2.)) / (-log(1./20.))
        expect = array([[0., expect],
                        [expect, 0.]])
        result = buildMutinfoMatrix(msa, norm=True)
        assert_array_almost_equal(expect, result, err_msg='norm failed')


class TestCalcMSAOccupancy(TestCase):

    def testResidueCount(self):

        assert_array_equal(calcMSAOccupancy(FASTA, 'residue', count=1),
                           FASTA_ALPHA.sum(0))

    def testSequenceCount(self):

        assert_array_equal(calcMSAOccupancy(FASTA, 'sequence', count=1),
                           FASTA_ALPHA.sum(1))

    def testResidueOccupancy(self):

        assert_array_equal(calcMSAOccupancy(FASTA, 'residue'),
                           FASTA_ALPHA.sum(0) / (FASTA.numSequences() * 1.0))

    def testSequenceOccupancy(self):

        assert_array_equal(calcMSAOccupancy(FASTA, 'sequence'),
                           FASTA_ALPHA.sum(1) / (FASTA.numResidues() * 1.0))


class TestIdentity(TestCase):

    def testIdentityMatrix(self):

        assert_array_almost_equal(FASTA_EYE, buildSeqidMatrix(FASTA))

    def testIdentityMatrixNonTurbo(self):

        assert_array_almost_equal(FASTA_EYE,
                                  buildSeqidMatrix(FASTA, turbo=False))


class TestUnique(TestCase):

    def testUnique(self):

        seqid = 0.98
        unique = ones(FASTA_NUMBER, bool)
        for i in range(FASTA_NUMBER):
            if not unique[i]:
                continue
            for j in range(i+1, FASTA_NUMBER):
                if FASTA_EYE[i, j] >= seqid:
                    unique[j] = False

        assert_array_equal(unique, uniqueSequences(FASTA, seqid))

    def testUnique2(self):

        seqid = 0.5
        unique = ones(FASTA_NUMBER, bool)
        for i in range(FASTA_NUMBER):
            if not unique[i]:
                continue
            for j in range(i+1, FASTA_NUMBER):
                if FASTA_EYE[i, j] >= seqid:
                    unique[j] = False

        assert_array_equal(unique, uniqueSequences(FASTA, seqid))


class TestCalcOMES(TestCase):

    def testZero(self):

        msa = array([list('ACCA'),
                     list('ACDA'),
                     list('ACCC'),
                     list('ACDC')], dtype='|S1')

        expect = array([[0., 0., 0., 0.],
                        [0., 0., 0., 0.],
                        [0., 0., 0., 0.],
                        [0., 0., 0., 0.], ])
        result = buildOMESMatrix(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testFourSequences(self):

        msa = array([list('ACCA'),
                     list('ACDA'),
                     list('ACDC'),
                     list('ACDC')], dtype='|S1')

        expect = array([[0., 0., 0., 0.],
                        [0., 0., 0., 0.],
                        [0., 0., 0., 4./3],
                        [0., 0., 4./3, 0.], ])
        result = buildOMESMatrix(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testTwenty(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, s] for s in seq], dtype='|S1')

        expect = array([[0., 380.],
                        [380., 0.]])
        result = buildOMESMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testTwentyReversed(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, seq[-i-1]] for i, s in enumerate(seq)], dtype='|S1')

        expect = array([[0., 380.],
                        [380., 0.]])
        result = buildOMESMatrix(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity(self):

        msa = array([list('OX'),
                     list('XO')], dtype='|S1')

        expect = array([[0., 2.],
                        [2., 0.]])
        result = buildOMESMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testNoAmbiguity(self):

        msa = array([list('OX'),
                     list('XO')], dtype='|S1')

        expect = array([[0., 2.],
                        [2., 0.]])
        result = buildOMESMatrix(msa, ambiquity=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, ambiquity=False, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity2(self):

        msa = array([list('AB'),
                     list('BZ')], dtype='|S1')
        expect = array([[0., 2.],
                        [2., 0.]])
        result = buildOMESMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity3(self):

        msa = array([list('XX')], dtype='|S1')

        expect = zeros((2, 2))
        result = buildOMESMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity4(self):

        msa = array([list('Bb'),
                     list('jJ'),
                     list('Zz'), ], dtype='|S1')

        expect = array([[0., 6.],
                        [6., 0.]])
        result = buildOMESMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testAmbiguity5(self):

        expect = array([[0., 0.],
                        [0., 0.]])

        for seq in ['bx', 'Xb', 'jX', 'Xj', 'xz', 'ZX',
                    'bj', 'jb', 'bz', 'zb', 'jz', 'zj']:
            msa = array([list(seq)], dtype='|S1')
            result = buildOMESMatrix(msa, debug=False)
            assert_array_almost_equal(expect, result, err_msg=seq + ' failed')

    def testAmbiguity6(self):

        expect = zeros((2, 2))

        for seq in ['bb', 'jj', 'zz']:
            msa = array([list(seq)], dtype='|S1')
            result = buildOMESMatrix(msa, debug=False)
            assert_array_almost_equal(expect, result, err_msg=seq + ' failed')

    def testAmbiguity7(self):

        msa = array([list('bx'),
                     list('xb')], dtype='|S1')
        expect = array([[0., 162./121],
                        [162./121, 0.]])
        result = buildOMESMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testInf(self):

        msa = zeros((500, 10), '|S1')
        msa.fill('.')
        msa[95, 8] = 's'
        msa[95, 9] = 'i'
        expect = zeros((10, 10))
        expect[8, 9] = expect[9, 8] = 500.
        result = buildOMESMatrix(msa, debug=False)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildOMESMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


class TestCalcSCA(TestCase):

    def testZero(self):

        msa = array([list('ACCD'),
                     list('ACDD'),
                     list('ACCC'),
                     list('ACDC')], dtype='|S1')

        expect = array([log(0.975/.025)*.5, log(0.95/.05)*.5])
        weight = ((expect ** 2).sum())**.5
        expect = expect / weight * array([log(0.975/.025), log(0.95/.05)])
        expect = (expect ** 2).mean() - (expect.mean()) ** 2
        expect = array([[0., 0., 0., 0.],
                        [0., 0., 0., 0.],
                        [0., 0., expect, 0.],
                        [0., 0., 0., expect], ])
        result = buildSCAMatrix(msa)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildSCAMatrix(msa, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')

    def testMATLAB(self):

        sca = fromfile(pathDatafile('msa_Cys_knot_sca.dat'))
        expect = sca.reshape((10, 10))
        fasta = FASTA[:, :10]
        result = buildSCAMatrix(fasta, turbo=True)
        assert_array_almost_equal(expect, result, err_msg='turbo failed')
        result = buildSCAMatrix(fasta, turbo=False)
        assert_array_almost_equal(expect, result, err_msg='w/out turbo failed')


class TestCalcMeff(TestCase):

    def testZero1(self):

        msa = array([list('ACCD')] * 100, dtype='|S1')
        expect = 1.
        result = calcMeff(msa)
        assert_array_almost_equal(expect, result)

        result = calcMeff(msa, seqid=0.9)
        assert_array_almost_equal(expect, result)

        result = calcMeff(msa, seqid=0.4)
        assert_array_almost_equal(expect, result)

        expect = (1., zeros((100)) + 1./100)
        result = calcMeff(msa, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        result = calcMeff(msa, seqid=0.9, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        result = calcMeff(msa, seqid=0.4, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

    def testZero2(self):

        msa = array([list('AACC')] * 50 + [list('AADD')] * 50, dtype='|S1')
        expect = 2.
        result = calcMeff(msa)
        assert_array_almost_equal(expect, result)

        result = calcMeff(msa, seqid=0.9)
        assert_array_almost_equal(expect, result)

        expect = 1.
        result = calcMeff(msa, seqid=0.4)
        assert_array_almost_equal(expect, result)

        expect = (2., zeros((100)) + 1./50)
        result = calcMeff(msa, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        result = calcMeff(msa, seqid=0.9, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        expect = (1., zeros((100)) + 1./100)
        result = calcMeff(msa, seqid=0.4, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

    def testTwenty(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, s] for s in seq], dtype='|S1')
        expect = 20.
        result = calcMeff(msa)
        assert_array_almost_equal(expect, result)

        result = calcMeff(msa, seqid=0.9)
        assert_array_almost_equal(expect, result)

        result = calcMeff(msa, seqid=0.4)
        assert_array_almost_equal(expect, result)

        expect = (20., ones(20))
        result = calcMeff(msa, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        result = calcMeff(msa, seqid=0.9, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        result = calcMeff(msa, seqid=0.4, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

    def testTwentyReversed(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, seq[-i-1]] for i, s in enumerate(seq)], dtype='|S1')
        expect = 20.
        result = calcMeff(msa)
        assert_array_almost_equal(expect, result)

        result = calcMeff(msa, seqid=0.9)
        assert_array_almost_equal(expect, result)

        result = calcMeff(msa, seqid=0.4)
        assert_array_almost_equal(expect, result)

        expect = (20., ones(20))
        result = calcMeff(msa, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        result = calcMeff(msa, seqid=0.9, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        result = calcMeff(msa, seqid=0.4, weight=True)
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

    def testMATLAB(self):

        expect = 1.8416666666666664e+01
        result = calcMeff(FASTA, refine=True)
        assert_array_almost_equal(expect, result)
        result = calcMeff(FASTA, refine=True, weight=True)
        expect = (expect,
                  array([1./3, 1./3, 1./4, 1./2, 1., 1., 1., 1., 1./3, 1./3,
                         1./3, 1./2, 1./2, 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                         1., 1./2, 1./2], dtype='float'))
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        expect = 23
        result = calcMeff(FASTA, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result)
        result = calcMeff(FASTA, seqid=0.9, refine=True, weight=True)
        expect = (expect,
                  array([0.5, 0.5, 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.5,
                         0.5, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                         1.], dtype='float'))
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')

        expect = 4.66689144189144
        result = calcMeff(FASTA, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result)
        result = calcMeff(FASTA, seqid=0.4, refine=True, weight=True)
        expect = (expect,
                  array([0.071428571, 0.071428571, 0.066666667, 0.083333333,
                         0.076923077, 0.066666667, 0.076923077, 0.090909091,
                         0.090909091, 0.071428571, 0.076923077, 0.076923077,
                         0.1, 0.2, 0.166666667, 0.25, 1., 0.25, 0.25, 0.25,
                         0.25, 0.333333333, 0.5, 0.071428571,
                         0.125, ], dtype='float'))
        assert_array_almost_equal(expect[0], result[0],
                                  err_msg='weight failed')
        assert_array_almost_equal(expect[1], result[1],
                                  err_msg='weight failed')


class TestDirectInfo(TestCase):

    def testZero(self):

        msa = array([list('ACCY')] * 100, dtype='|S1')
        expect = array([[0., 0.66325166608, 0.66325166608, 0.66222154839],
                        [0.66325166608, 0., 0.66325166608, 0.66222154839],
                        [0.66325166608, 0.66325166608, 0., 0.66222154839],
                        [0.66222154839, 0.66222154839, 0.66222154839, 0.], ])
        result = buildDirectInfoMatrix(msa)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.3
        expect = array([[0.,  0.53477644,  0.53477644,  0.58513139],
                        [0.53477644,  0.,  0.53477644,  0.58513139],
                        [0.53477644,  0.53477644,  0.,  0.58513139],
                        [0.58513139,  0.58513139,  0.58513139,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.8
        expect = array([[0.,  0.41397564,  0.41397564,  0.41295744],
                        [0.41397564,  0.,  0.41397564,  0.41295744],
                        [0.41397564,  0.41397564,  0.,  0.41295744],
                        [0.41295744,  0.41295744,  0.41295744,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testZero1(self):

        msa = array([list('ACCD')] * 100, dtype='|S1')
        expect = array([[0.,  0.13010138,  0.13010138,  0.13008827],
                        [0.13010138,  0.,  0.13010138,  0.13008827],
                        [0.13010138,  0.13010138,  0.,  0.13008827],
                        [0.13008827,  0.13008827,  0.13008827,  0.]])
        result = buildDirectInfoMatrix(msa)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.3
        expect = array([[0.,  0.23225802,  0.23225802,  0.23225606],
                        [0.23225802,  0.,  0.23225802,  0.23225606],
                        [0.23225802,  0.23225802,  0.,  0.23225606],
                        [0.23225606,  0.23225606,  0.23225606,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.8
        expect = array([[0.,  0.03910082,  0.03910082,  0.03910178],
                        [0.03910082,  0.,  0.03910082,  0.03910178],
                        [0.03910082,  0.03910082,  0.,  0.03910178],
                        [0.03910178,  0.03910178,  0.03910178,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testZero2(self):

        msa = array([list('AAYY')] * 50 + [list('YYDD')] * 50, dtype='|S1')
        expect = array([[0., 1.0248086877, 1.0001784999, 1.0001784999],
                        [1.0248086877, 0., 1.0001784999, 1.0001784999],
                        [1.0001784999, 1.0001784999, 0., 1.0248086877],
                        [1.0001784999, 1.0001784999, 1.0248086877, 0.], ])
        result = buildDirectInfoMatrix(msa)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.3
        expect = array([[0.,  1.30637709,  1.00041168,  1.00041168],
                       [1.30637709,  0.,  1.00041168,  1.00041168],
                       [1.00041168,  1.00041168,  0.,  1.30637709],
                       [1.00041168,  1.00041168,  1.30637709,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.8
        expect = array([[0.,  0.55686383,  0.55605426,  0.55605426],
                       [0.55686383,  0.,  0.55605426,  0.55605426],
                       [0.55605426,  0.55605426,  0.,  0.55686383],
                       [0.55605426,  0.55605426,  0.55686383,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testZero3(self):

        msa = array([list('AACC')] * 50 + [list('CCDD')] * 50, dtype='|S1')
        expect = array([[0.,  0.23179074,  0.23178758,  0.23178758],
                        [0.23179074,  0.,  0.23178758,  0.23178758],
                        [0.23178758,  0.23178758,  0.,  0.23178758],
                        [0.23178758,  0.23178758,  0.23178758,  0.]])
        result = buildDirectInfoMatrix(msa)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.3
        expect = array([[0.,  0.54966352,  0.54965557,  0.54965557],
                       [0.54966352,  0.,  0.54965557,  0.54965557],
                       [0.54965557,  0.54965557,  0.,  0.54965565],
                       [0.54965557,  0.54965557,  0.54965565,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.8
        expect = array([[0.,  0.04349673,  0.0434971,  0.0434971],
                        [0.04349673,  0.,  0.0434971,  0.0434971],
                        [0.0434971,  0.0434971,  0.,  0.04349711],
                        [0.0434971,  0.0434971,  0.04349711,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testZero4(self):

        msa = array([list('AAYY')] * 50 + [list('AADD')] * 50, dtype='|S1')
        expect = array([[0.,  0.68512031,  0.54656421,  0.54656421],
                        [0.68512031,  0.,  0.54656421,  0.54656421],
                        [0.54656421,  0.54656421,  0.,  1.05788193],
                        [0.54656421,  0.54656421,  1.05788193,  0.]])
        result = buildDirectInfoMatrix(msa)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.3
        expect = array([[0.,  0.65929972,  0.52603587,  0.52603587],
                       [0.65929972,  0.,  0.52603587,  0.52603587],
                       [0.52603587,  0.52603587,  0.,  1.99332504],
                       [0.52603587,  0.52603587,  1.99332504,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.8
        expect = array([[0.,  0.48628027,  0.27348214,  0.27348214],
                       [0.48628027,  0.,  0.27348214,  0.27348214],
                       [0.27348214,  0.27348214,  0.,  0.65235465],
                       [0.27348214,  0.27348214,  0.65235465,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testTwenty(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, s] for s in seq], dtype='|S1')
        expect = array([[0., 3.0302471958885744],
                        [3.0302471958885744, 0.]])
        result = buildDirectInfoMatrix(msa)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.3
        expect = array([[0.,  3.0757238],
                        [3.0757238,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.8
        expect = array([[0.,  1.89488914],
                        [1.89488914,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testTwentyReversed(self):

        seq = 'ACDEFGHIKLMNPQRSTVWY'
        msa = array([[s, seq[-i-1]] for i, s in enumerate(seq)], dtype='|S1')
        expect = array([[0., 3.030670764986982],
                        [3.030670764986982, 0.]])
        result = buildDirectInfoMatrix(msa)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.3
        expect = array([[0.,  3.12974624],
                        [3.12974624,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.3)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.3, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        # pseudo_weight=0.8
        expect = array([[0.,  1.89488914],
                        [1.89488914,  0.]])
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(msa, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.9, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.9, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
        result = buildDirectInfoMatrix(msa, seqid=0.4, pseudo_weight=0.8)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(
            msa, seqid=0.4, pseudo_weight=0.8, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testMATLAB8(self):

        di = fromfile(pathDatafile('msa_Cys_knot_di.dat'))
        expect = di.reshape((8, 8))
        fasta = FASTA[:, :8]
        result = buildDirectInfoMatrix(fasta)
        assert_array_almost_equal(
            expect, result, err_msg='w/out refine failed')
        result = buildDirectInfoMatrix(fasta, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')

    def testMATLAB10(self):

        di = fromfile(pathDatafile('msa_Cys_knot_di.dat'))
        expect = di.reshape((8, 8))
        fasta = FASTA[:, :10]
        result = buildDirectInfoMatrix(fasta, refine=True)
        assert_array_almost_equal(expect, result, err_msg='refine failed')
