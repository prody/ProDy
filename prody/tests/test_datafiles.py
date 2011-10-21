#!/usr/bin/python
# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import sys
import os.path
import unittest
import inspect
import tempfile

import prody


try:
    import nose
except ImportError:
    NONOSE = True
    NONOSE_MSG = 'nose could not be imported'
else:
    NONOSE = False
    NONOSE_MSG = ''

TEMPDIR = tempfile.gettempdir()
DATA_FILES = {
    'multi_model_truncated': {
        'pdb': '2k39',
        'file': 'pdb2k39_truncated.pdb',
        'atoms': 167,
        'models': 3
    },
    'dssp': {
        'pdb': '1r19',
        'file': 'pdb1r19_dssp.pdb',
        'atoms': 8216,
        'models': 1
    },
    'oneatom': {
        'pdb': '1ejg',
        'file': 'pdb1ejg_oneatom.pdb',
        'atoms': 1,
        'models': 1
    },
    '1ubi': {
        'pdb': '1ubi',
        'file': 'pdb1ubi_ca.pdb',
        'n_atoms': 76,
        'models': 1
    },
    'anm1ubi_hessian': {
        'file': 'anm1ubi_hessian.coo',
    },
    'anm1ubi_evalues': {
        'file': 'anm1ubi_evalues.dat',
    },
    'anm1ubi_vectors': {
        'file': 'anm1ubi_vectors.dat'
    },
    'gnm1ubi_kirchhoff': {
        'file': 'gnm1ubi_kirchhoff.coo',
    },
    'gnm1ubi_evalues': {
        'file': 'gnm1ubi_evalues.dat',
    },
    'gnm1ubi_vectors': {
        'file': 'gnm1ubi_vectors.dat'
    },
}


PARSERS = {
    '.pdb': prody.parsePDB,
    '.coo': prody.parseSparseMatrix,
    '.dat': prody.parseArray,
    '.gz': lambda fn, **kwargs: PARSERS[os.path.splitext(fn)[1]](fn, **kwargs)      
}
TESTS_PATH = os.path.abspath(os.path.split(inspect.getfile(
                                                   inspect.currentframe()))[0])
__all__ = ['parseDatafile', 'getDatafilePath', 
           'DATA_FILES', 'TEMPDIR', 'NONOSE', 'NONOSE_MSG']


def getDatafilePath(filename):

    assert isinstance(filename, str), 'filename must be a string'
    fn = os.path.join(TESTS_PATH, 'data', filename)
    assert os.path.isfile(fn), 'No such file: "{0:s}"'.format(fn)
    return fn

def parseDatafile(filename, **kwargs):
    """*filename* must be present in :file:`prody/tests/data` folder."""
    
    if filename in DATA_FILES:
        filename = DATA_FILES[filename]['file']
    fn = getDatafilePath(filename)
    return PARSERS[os.path.splitext(fn)[1]](fn, **kwargs)

for name, value in DATA_FILES.iteritems():
    value['path'] = getDatafilePath(value['file'])


class TestDatafilesMeta(type):
    
    def __init__(cls, name, bases, dict):
        
        for name, value in DATA_FILES.iteritems():
            fn = value['file']
            def testFunction(self, filename=fn, **kwargs):
                
                self.assertTrue(os.path.isfile(getDatafilePath(filename)))
            
            testFunction.__name__ = 'testDatafile_{0:s}'.format(name)
            testFunction.__doc__ = 'Test presence of "{0:s}"'.format(fn)
            setattr(cls, testFunction.__name__, testFunction)


class TestDatafiles(unittest.TestCase):

    """Test presence of data files."""
    
    __metaclass__ = TestDatafilesMeta
    
