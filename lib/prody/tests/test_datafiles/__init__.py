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

"""This module defines unit test data files and functions to access them.
Data files used in tests are truncated PDB files, e.g. most of atoms and/or
models and/or header sections are removed for having a compact installation
package that contains test modules and files as well."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'


from os.path import join, isfile, split, splitext
from unittest import TestCase

from numpy import array

from prody import parsePDB, parseDCD, parseSparseMatrix, parseArray
from prody.tests import TEMPDIR, TESTDIR


DATA_FILES = {
    'multi_model_truncated': {
        'pdb': '2k39',
        'file': 'pdb2k39_truncated.pdb',
        'atoms': 167,
        'models': 3,
        'rmsd_all_aligned': array([0.000, 1.380, 1.745]),
        'rmsd_bb_aligned': array([0.000, 0.367, 0.395]),
        'rmsd_ca_aligned': array([0.000, 0.393, 0.411]),
        'rmsd_noh_aligned': array([0.000, 0.985, 1.285]),
        'rmsd_all': array([0.000, 1.539, 1.964]),
        'rmsd_bb': array([0.000, 0.654, 0.678]),
        'rmsd_ca': array([0.000, 0.680, 0.690]),
        'rmsd_noh': array([0.000, 1.191, 1.543]),
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
    '1ejg': {
        'pdb': '1ejg',
        'file': 'pdb1ejg.pdb',
        'atoms': 637,
        'models': 1
    },
    '1ubi': {
        'pdb': '1ubi',
        'file': 'pdb1ubi.pdb',
        'n_atoms': 683,
        'models': 1
    },
    '1ubi_ca': {
        'pdb': '1ubi',
        'file': 'pdb1ubi_ca.pdb',
        'n_atoms': 76,
        'models': 1
    },
    'dcd': {
        'file': 'dcd2k39_truncated.dcd',
        'atoms': 167,
        'models': 3
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
    'sti': {
        'file': 'xmlSTI.xml'
    },
    '3mht': {
        'file': 'pdb3mht.pdb'
    },
    'Fasta': {
        'file': 'msa_Cys_knot.fasta'
    },
    'Selex': {
        'file': 'msa_Cys_knot.slx'
    },
    'Stockholm': {
        'file': 'msa_Cys_knot.sth'
    }
}


PARSERS = {
    '.dcd': parseDCD, '.pdb': parsePDB,
    '.coo': parseSparseMatrix, '.dat': parseArray,
    '.gz': lambda fn, **kwargs: PARSERS[splitext(fn)[1]](fn, **kwargs)      
}


__all__ = ['parseDatafile', 'pathDatafile', 'DATA_FILES']


def pathDatafile(filename):

    try:
        filename = DATA_FILES[filename]['file']
    except KeyError:
        pass
    assert isinstance(filename, str), 'filename must be a string'
    fn = join(TESTDIR, 'test_datafiles', filename)
    assert isfile(fn), 'No such file: "{0:s}"'.format(fn)
    return fn


def parseDatafile(filename, **kwargs):
    """*filename* must be present in :file:`prody/tests/test_datafiles`."""
    
    if filename in DATA_FILES:
        filename = DATA_FILES[filename]['file']
    fn = pathDatafile(filename)
    return PARSERS[splitext(fn)[1]](fn, **kwargs)


for name, value in DATA_FILES.iteritems():
    value['path'] = pathDatafile(value['file'])


class TestDatafiles(TestCase):

    """Test presence of data files."""
    
    pass

for name, value in DATA_FILES.items():
    fn = value['file']
    def func(self, filename=fn, **kwargs):
        
        self.assertTrue(isfile(pathDatafile(filename)))
    
    func.__name__ = 'testDatafile_{0:s}'.format(name)
    func.__doc__ = 'Test presence of "{0:s}"'.format(fn)
    setattr(TestDatafiles, func.__name__, func)
del func
