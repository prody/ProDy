"""This module defines unit test data files and functions to access them.
Data files used in tests are truncated PDB files, e.g. most of atoms and/or
models and/or header sections are removed for having a compact installation
package that contains test modules and files as well."""


from os.path import join, isfile, split, splitext
from prody.tests import TestCase

from numpy import array
import numpy as np

from prody import parsePDB, parseDCD, parseMMCIF
from prody import parseSparseMatrix, parseArray, loadModel
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
    '2nwl': {
        'pdb': '2nwl',
        'file': 'pdb2nwl-opm.pdb',
        'n_atoms': 5227,
        'models': 1
    },
    '2k39_ca': {
        'pdb': '2k39',
        'file': 'pdb2k39_ca.pdb',
        'n_atoms': 76,
        'models': 116
    },
    '1ubi_ca': {
        'pdb': '1ubi',
        'file': 'pdb1ubi_ca.pdb',
        'n_atoms': 76,
        'models': 1
    },
    '2gb1_truncated': {
        'pdb': '2gb1',
        'file': 'pdb2gb1_truncated.pdb',
        'n_atoms': 26,
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
    'rtb2gb1_hessian': {
        'file': 'rtb2gb1_hessian.coo'
    },
    'rtb2gb1_project': {
        'file': 'rtb2gb1_project.coo'
    },
    'pca2k39_cov': {
        'file': 'pca2k39_cov.coo',
    },
    'pca2k39_evalues': {
        'file': 'pca2k39_evalues.dat',
    },
    'pca2k39_vectors': {
        'file': 'pca2k39_vectors.dat'
    },
    'commute1ubi': {
        'file': 'commute1ubi.dat'
    },
    'hit1ubi': {
        'file': 'hit1ubi.dat'
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
    },
    'SCA': {
        'file': 'msa_Cys_knot_sca.dat'
    },
    'DI': {
        'file': 'msa_Cys_knot_di.dat'
    },
    'RTER': {
        'file': 'pdbRTER.pdb'
    },
    'five_digits': {
        'file': 'pdb1tw7_step3_charmm2namd.pdb'
    },
    'hex': {
        'file': 'pdb1tw7_step3_charmm2namd_doubled_hex.pdb'
    },
    'h36': {
        'file': 'pdb1tw7_step3_charmm2namd_doubled_h36.pdb'
    },
    'multi_model_cif': {
        'pdb': '6yfy',
        'file': 'mmcif_6yfy.cif',
        'atoms': 1460,
        'chainA_atoms': 187,
        'ca_atoms': 36,
        'bb_atoms': 144,
        'models': 26,
    },
    'long_chid_cif': {
        'pdb': '6zu5',
        'file': 'mmcif_6zu5.cif',
        'atoms': 165175,
        'chain_SX0_atoms': 1089,
    },
    '3hsy': {
        'pdb': '3hsy',
        'file': 'pdb3hsy.pdb',
        'atoms': 6058,
        'ca_atoms': 730,
        'models': 1,
        'biomols': 1
    },
    '3o21': {
        'pdb': '3o21',
        'file': 'pdb3o21.pdb',
        'atoms': 12793,
        'ca_atoms': 1489,
        'models': 1,
        'biomols': 2
    },
    '3p3w': {
        'pdb': '3p3w',
        'file': 'pdb3p3w.pdb',
        'atoms': 11484,
        'ca_atoms': 1482,
        'models': 1,
        'biomols': 2
    },
}


PARSERS = {
    '.dcd': parseDCD, '.pdb': parsePDB, '.cif': parseMMCIF,
    '.coo': parseSparseMatrix, '.dat': parseArray,
    '.txt': np.loadtxt,
    '.gz': lambda fn, **kwargs: PARSERS[splitext(fn)[1]](fn, **kwargs)
}


__all__ = ['parseDatafile', 'pathDatafile', 'DATA_FILES']


def pathDatafile(filename):

    try:
        filename = DATA_FILES[filename]['file']
    except KeyError:
        pass
    assert isinstance(filename, str), 'filename must be a string'
    fn = join(TESTDIR, 'datafiles', filename)
    assert isfile(fn), 'No such file: "{0:s}"'.format(fn)
    return fn


def parseDatafile(filename, **kwargs):
    """*filename* must be present in :file:`prody/tests/datafiles`."""

    if filename in DATA_FILES:
        filename = DATA_FILES[filename]['file']
    fn = pathDatafile(filename)
    return PARSERS[splitext(fn)[1]](fn, **kwargs)


for name, value in DATA_FILES.items():
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
