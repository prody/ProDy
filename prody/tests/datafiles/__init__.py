"""This module defines unit test data files and functions to access them.
Data files used in tests are truncated PDB files, e.g. most of atoms and/or
models and/or header sections are removed for having a compact installation
package that contains test modules and files as well."""


from os.path import join, isfile, split, splitext
from prody.tests import TestCase

from numpy import array
import numpy as np

from prody import parsePDB, parseDCD, parseMMCIF, parseMMTF
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
        'mmtf': '1ubi.mmtf',        
        'n_atoms': 683,
        'models': 1
    },
    '1ubi_addH': {
        'pdb': '1ubi',
        'file': 'pdb_addH_1ubi.pdb',   
        'n_atoms': 1474,
        'models': 1,
        'n_h': 791
    },
    '1ubi_mmtf': {
        'pdb': '1UBI',
        'file': '1ubi.mmtf',
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
    '2k39_mmtf': {
        'pdb': '2K39',
        'file': '2k39.mmtf',
        'n_atoms': 1231,
        'models': 116
    },    
    '3enl_pdb': {
        'pdb': '3enl',
        'file': 'pdb3enl.pdb',
        'n_atoms': 7294,
        'models': 1
    },
    '3enl_addH': {
        'pdb': '3enl',
        'file': 'addH_pdb3enl.pdb',   
        'n_atoms': 7641,
        'models': 1,
        'n_h': 3999
    }, 
    '3enl_mmtf': {
        'pdb': '3ENL',
        'file': 'mmtf3enl.mmtf',
        'n_atoms': 7294,
        'models': 1
    },   
    '1pwc_pdb': {
        'pdb': '1pwc',
        'file': '1pwc.pdb',
        'n_atoms': 3129,
        'models': 1
    },    
    '1pwc_mmtf': {
        'pdb': '1PWC',
        'file': '1pwc.mmtf',
        'n_atoms': 3129,
        'models': 1
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
    'lda2k39_evalues': {
        'file': 'lda2k39_evalues.dat',
    },
    'lda2k39_vectors': {
        'file': 'lda2k39_vectors.dat'
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
    'MSA': {
        'file': 'msa_3hsyA_3o21A.fasta'
    },
    'MSA_MEW': {
        'file': 'msa_3hsyA_3o21A_new.fasta'
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
    'hex_ter': {
        'file': 'pdb4v8r_hex.pdb'
    },
    'h36_ter': {
        'file': 'pdb4v8r_h36.pdb'
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
    'biomols_cif': {
        'pdb': '3o21',
        'file': 'mmcif_3o21.cif',
        'atoms': 12793,
        'bm0_atoms': 6281,
        'num_chains_united': 4,
        'bm_chains_united': [2, 2],
        'bm_chains_alone': [9, 10],
        'chainA_atoms_united': 3170,
        'chainA_atoms_alone': 3025,
        'ca_atoms': 1489,
        'bb_atoms': 5956,
        'models': 1,
        'unobs_B_start': 'G-------------------------------NQNTTEK-'
    },
    'pymol_cif': {
        'pdb': '3o21',
        'file': 'mmcif_3o21_pymol.cif',
        'atoms': 12793,
        'last_coords': array([ 50.327, -11.971, -19.976])
    },
    'big_biomols_cif': {
        'pdb': '7cth',
        'file': 'mmcif_7cth.cif',
        'atoms': 16647,
        'bm0_atoms': 16647 * 60,
        'num_chains': 14,
        'bm0_chains': 840,
    },
    'boltz_wrapped_line_cif':{
        'file': 'NRAS_BRAF_GDP_model_0.cif',
        'atoms': 7333,
        'num_chains': 3
    },
    'long_chid_cif': {
        'pdb': '6zu5',
        'file': 'mmcif_6zu5.cif',
        'atoms': 165175,
        'segment_SX0_atoms': 1089,
    },
    'chimerax_cif': {
        'pdb': '1ake',
        'file': 'mmcif_1ake_chimerax.cif',
        'biomols': 2,
        'bm0_atoms': 1954
    },
    '6zu5_sel': {
        'pdb': '6zu5',
        'file': '6zu5_sel_SE0_SF0_10-20.pdb',
        'atoms': 22,
        'segment_SF0_atoms': 11,
        'chid_order': ['Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z',
                       'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'],
        'sorted_order': ['B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B',
                         'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z']
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
        'chainA_atoms': 3170,
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
    '6flr': {
        'pdb': '6flr',
        'file': 'pdb6flr.pdb',
        'atoms_single': 6073,
        'ca_atoms_single': 741,
        'atoms_altloc': 6086,
        'ca_atoms_altloc': 743,
        'num_altlocs': 2,
        'biomols': 1,
        'anisousA': array([[1.3327, 0.8826, 0.7048, 0.0444, 0.3365, 0.1618]]),
        'anisousB': array([[1.3306, 0.8825, 0.7058, 0.0438, 0.338 , 0.1608]])
    },
    '6flr_sel': {
        'pdb': '6flr',
        'file': 'pdb6flr_B234.pdb',
        'atoms': 1,
        'ca_atoms': 1,
        'models': 2,
        'biomols': 1,
        'anisousA': array([[1.3327, 0.8826, 0.7048, 0.0444, 0.3365, 0.1618]]),
        'anisousB': array([[1.3306, 0.8825, 0.7058, 0.0438, 0.338 , 0.1608]])
    },
    'cif_6flr': {
        'pdb': '6flr',
        'file': 'cif_6flr.cif',
        'atoms_single': 6073,
        'ca_atoms_single': 741,
        'atoms_altloc': 6086,
        'ca_atoms_altloc': 743,
        'num_altlocs': 2,
        'biomols': 1,
        'anisousA': array([[1.3327, 0.8826, 0.7048, 0.0444, 0.3365, 0.1618]]),
        'anisousB': array([[1.3306, 0.8825, 0.7058, 0.0438, 0.338 , 0.1608]])
    },
    'gromacs': {
        'pdb': '6fpj',
        'file': 'pdb6fpj_Bb_fixed_solv_ions.pdb',
        'atoms': 83473,
        'ca_atoms': 758,
        'models': 1,
        'biomols': 1,
        'protein': 12194
    },
    'cath': {
        'file': 'cath.xml',
    },    
    '2k39_insty': {
        'file': '2k39_insty.pdb'
    },    
    '2k39_insty_first': {
        'file': '2k39_insty_first.pdb'
    },    
    '2k39_insty_dcd': {
        'file': '2k39_insty.dcd'
    },
    '2k39_hbs': {
        'file': '2k39_hbs.npy'
    },    
    '2k39_sbs': {
        'file': '2k39_sbs.npy'
    },    
    '2k39_rib': {
        'file': '2k39_rib.npy'
    },    
    '2k39_PiStack': {
        'file': '2k39_PiStack.npy'
    },
    '2k39_PiCat': {
        'file': '2k39_PiCat.npy'
    },
    '2k39_hph': {
        'file': '2k39_hph.npy'
    },
    '2k39_disu': {
        'file': '2k39_disu.npy'
    },
    '3o21_disu': {
        'file': '3o21_disu.npy'
    },
    '2k39_all': {
        'file': '2k39_all.npy'
    },
    '2k39_hph2': {
        'file': '2k39_hph2.npy'
    },    
    '2k39_all2': {
        'file': '2k39_all2.npy'
    },
    'pdb4ake_fixed': {
        'file': '4akeA_alg_fixed.pdb'
    },
    'mrc1ake': {
        'file': '1ake.mrc'
    },
    'pdb1ake': {
        'file': 'pdb1ake.pdb'
    },
    'pdb7pbl': {
        'file': 'pdb7pbl.pdb',
        'atoms': 31396,
        'nucleoside': 252
    },
    'pqrUnknown': {
        'file': 'pqr_snippet1.pqr',
        'atoms': 5,
        'models': 1
    },
    'pqrTranscomp': {
        'file': 'pqr_snippet2_transcomp.pqr',
        'atoms': 5,
        'models': 1
    },
    'pqrFpocket': {
        'file': 'pqr_snippet3_fpocket.pqr',
        'atoms': 5,
        'models': 1
    },
    'pqrPymol': {
        'file': 'pqr_snippet4_pymol.pqr',
        'atoms': 5,
        'models': 1
    },
    'anmd': {
        'pdb': '1ubi',
        'file': '1ubi_anmd_mode1_ens.pdb',
        'n_atoms': 683,
        'models': 5
    },
    'probes': {
        'file': 'probes_in.pdb',
        'n_atoms': 4,
        'long_resname': 'ACET',
        'short_resname': 'ACE'
    }
}


PARSERS = {
    '.dcd': parseDCD, '.pdb': parsePDB, '.cif': parseMMCIF,
    '.mmtf': parseMMTF,
    '.coo': parseSparseMatrix, '.dat': parseArray,
    '.txt': np.loadtxt,
    '.npy': lambda fn, **kwargs: np.load(fn, allow_pickle=True),
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
