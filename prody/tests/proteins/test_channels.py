"""This module contains unit tests for signature cavity/channel detection in
:mod:`~prody.proteins.channels`."""

import os
import shutil
import tempfile
from unittest.mock import MagicMock, patch

from prody import LOGGER, parsePDB
from prody.proteins.channels import calcSignatureCavities, calcSignatureChannels
from prody.tests import unittest
from prody.tests.datafiles import pathDatafile

LOGGER.verbosity = 'none'


def _seedLocalStructure(unaligned_dir, fixture_file, pdb_id):
    """Copy a local datafile fixture into *unaligned_dir* under the filename
    ProDy's local-PDB resolution expects for *pdb_id*, so parsePDB(pdb_id, ...)
    resolves it without touching the network."""

    shutil.copy(pathDatafile(fixture_file),
                os.path.join(unaligned_dir, 'pdb' + pdb_id + '.pdb'))


class TestCalcSignatureCavitiesInput(unittest.TestCase):

    def testRejectsTooShortStructureString(self):
        """A structures entry that is too short to contain a 4-character PDB
        code plus a chain identifier should be rejected before any PDB
        parsing is attempted."""

        with self.assertRaises(ValueError) as ctx:
            calcSignatureCavities(structures=['54'], ref_structure='5a46A')

        self.assertIn('PDB code', str(ctx.exception))
        self.assertIn('chain', str(ctx.exception))

    def testRejectsTooShortRefStructureString(self):
        """A ref_structure that is too short to contain a 4-character PDB
        code plus a chain identifier should be rejected the same way as a
        malformed structures entry, before any PDB parsing is attempted."""

        with self.assertRaises(ValueError) as ctx:
            calcSignatureCavities(structures=['5a46A'], ref_structure='54')

        self.assertIn('PDB code', str(ctx.exception))
        self.assertIn('chain', str(ctx.exception))


class TestCalcSignatureCavitiesAlignment(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='sigcav_')
        unaligned_dir = os.path.join(self.out_dir, 'unaligned')
        os.makedirs(unaligned_dir, exist_ok=True)
        # Reference (1ubiA) and one homolog (1ubjA, same coordinates under a
        # different fake code) -- both resolved locally, no network.
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubi')
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubj')

    def tearDown(self):
        shutil.rmtree(self.out_dir, ignore_errors=True)

    def testUsesCealignNotExternalTMalign(self):
        """Superposing a homolog onto the reference must not depend on an
        external TMalign binary -- it should succeed using ProDy's own
        CE-align, even though TMalign isn't installed on this machine."""

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        self.assertIsNotNone(result)

    def testSignaturePdbIsThresholdNamedAndIncludesReference(self):
        """calcSignatureCavities collapses thresholding into the same call
        (no separate getSignatureCavities step), defaults to a 20% threshold,
        and writes a combined PDB naming that threshold, containing the
        reference structure's own atoms alongside any consensus voxels."""

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        self.assertEqual(os.path.basename(result['signature_pdb']),
                          'signature_cavity_pt20.pdb')
        self.assertTrue(os.path.exists(result['signature_pdb']))

        written = parsePDB(result['signature_pdb'])
        self.assertIsNotNone(written)
        self.assertIsNotNone(written.select('protein'),
            'combined signature PDB must include the reference structure, '
            'not just voxel pseudo-atoms')

    def testResidueMappingNearSignatureVoxels(self):
        """When at least one voxel passes the (default 20%) threshold, the
        reference residues within distA=4.5 A of it must be reported, and
        mapped onto each homologue via msa_mappings, with a residue-summary
        file written alongside the signature PDB."""

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        written = parsePDB(result['signature_pdb'])
        fil = written.select('resname FIL')
        n_voxels = 0 if fil is None else len(fil)

        if n_voxels == 0:
            self.skipTest('No signature voxels found for this fixture pair; '
                           'residue-mapping is exercised by the real-data run.')

        self.assertIn('signature_cavity_residues', result)
        residues = result['signature_cavity_residues']
        self.assertIn('ref', residues)
        self.assertIn('1ubjA', residues)

        self.assertIn('signature_cavity_residues_path', result)
        self.assertTrue(os.path.exists(result['signature_cavity_residues_path']))

    def testEntropyColoringAndMsaOutput(self):
        """The reference's B-factor column must carry per-residue Shannon
        entropy (not left at the original crystallographic B-factors), and
        the ref-anchored alignment used to compute it must be written out as
        its own FASTA file alongside the numeric resnum mapping."""

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        self.assertIn('msa_fasta_path', result)
        self.assertTrue(os.path.exists(result['msa_fasta_path']))
        with open(result['msa_fasta_path']) as f:
            fasta = f.read()
        self.assertIn('>1ubiA', fasta)
        self.assertIn('>1ubjA', fasta)

        self.assertTrue(os.path.exists(os.path.join(self.out_dir, 'aligned_resnums.msa')))

        ref_written = parsePDB(os.path.join(self.out_dir, 'aligned', 'ref.pdb'))
        betas = ref_written.select('protein and name CA').getBetas()
        # 1ubjA is an identical copy of the reference -- every aligned column
        # has zero variation, so entropy must be uniformly ~0, not left at
        # whatever crystallographic B-factors pdb1ubi.pdb originally had.
        self.assertTrue((betas < 1e-6).all())


class TestCalcSignatureChannels(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='sigchan_')
        unaligned_dir = os.path.join(self.out_dir, 'unaligned')
        os.makedirs(unaligned_dir, exist_ok=True)
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubi')
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubj')

    def tearDown(self):
        shutil.rmtree(self.out_dir, ignore_errors=True)

    def testDispatchesToChannelFamilyNotCavityFamily(self):
        """calcSignatureChannels must reuse the same CE-align pipeline as
        calcSignatureCavities but dispatch to the channel-detection family,
        producing a distinctly-named signature_channel_pt20.pdb with
        'channel_records' (not 'cavity_records') in the result."""

        result = calcSignatureChannels(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        self.assertIsNotNone(result)
        self.assertIn('channel_records', result)
        self.assertNotIn('cavity_records', result)

        self.assertEqual(os.path.basename(result['signature_pdb']),
                          'signature_channel_pt20.pdb')
        self.assertTrue(os.path.exists(result['signature_pdb']))

        written = parsePDB(result['signature_pdb'])
        self.assertIsNotNone(written)
        self.assertIsNotNone(written.select('protein'),
            'combined signature PDB must include the reference structure, '
            'not just voxel pseudo-atoms')

    def testEntropyColoringAndMsaOutput(self):
        """calcSignatureChannels must get entropy coloring and MSA output for
        free from the shared pipeline, same as calcSignatureCavities."""

        result = calcSignatureChannels(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        self.assertIn('msa_fasta_path', result)
        self.assertTrue(os.path.exists(result['msa_fasta_path']))
        with open(result['msa_fasta_path']) as f:
            fasta = f.read()
        self.assertIn('>1ubiA', fasta)
        self.assertIn('>1ubjA', fasta)

        ref_written = parsePDB(os.path.join(self.out_dir, 'aligned', 'ref.pdb'))
        betas = ref_written.select('protein and name CA').getBetas()
        self.assertTrue((betas < 1e-6).all())

    @patch('prody.database.searchDali')
    def testDaliFilterKwargsOverridesDefaultCutoffs(self, mock_search_dali):
        """dali_filter_kwargs must be accepted by calcSignatureChannels too,
        merged over the default cutoffs and passed through to
        DaliRecord.filter, same as calcSignatureCavities."""

        mock_rec = MagicMock()
        mock_rec.isSuccess = True
        mock_rec.filter.return_value = ['1ubjA']
        mock_search_dali.return_value = mock_rec

        calcSignatureChannels(ref_structure='1ubiA',
                               output_path=self.out_dir,
                               dali_filter_kwargs={'cutoff_rmsd': 1.0})

        mock_rec.filter.assert_called_once_with(
            cutoff_len=0.5, cutoff_rmsd=1.0, cutoff_Z=8, stringency=True)

    def testMsaFastaAcceptedForEntropySource(self):
        """calcSignatureChannels must accept msa_fasta and use it as the
        entropy source, same as calcSignatureCavities."""

        ag = parsePDB(pathDatafile('pdb1ubi.pdb'))
        ref_seq = ag.select('protein and name CA').getSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>1ubiA\n{0}\n'.format(ref_seq))
            f.write('>unrelated_seq\n{0}\n'.format(ref_seq))

        result = calcSignatureChannels(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir, msa_fasta=msa_path)

        self.assertIsNotNone(result)
        self.assertNotIn('msa_fasta_path', result)

        ref_written = parsePDB(os.path.join(self.out_dir, 'aligned', 'ref.pdb'))
        betas = ref_written.select('protein and name CA').getBetas()
        self.assertTrue((betas < 1e-6).all())


class _FakeDaliRecord(object):
    """Minimal stand-in for DaliRecord: succeeds after a fixed number of
    fetch() calls, and returns a fixed filter() result."""

    def __init__(self, success_after=0, filter_result=None):
        self._success_after = success_after
        self.fetch_calls = 0
        self.filter_result = filter_result if filter_result is not None else []
        self.filter_calls = []

    @property
    def isSuccess(self):
        return self.fetch_calls >= self._success_after

    def fetch(self):
        self.fetch_calls += 1

    def filter(self, **kwargs):
        self.filter_calls.append(kwargs)
        return self.filter_result


class TestCalcSignatureCavitiesDaliDiscovery(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='sigcav_dali_')
        unaligned_dir = os.path.join(self.out_dir, 'unaligned')
        os.makedirs(unaligned_dir, exist_ok=True)
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubi')
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubj')

    def tearDown(self):
        shutil.rmtree(self.out_dir, ignore_errors=True)

    @patch('prody.database.searchDali')
    def testDaliDiscoveryUsesFilteredHitsAsStructures(self, mock_search_dali):
        """Omitting structures must trigger a DALI search against the
        reference, filtered with the agreed cutoffs, and feed the surviving
        hits into the same pipeline a user-supplied list would use."""

        mock_rec = MagicMock()
        mock_rec.isSuccess = True
        mock_rec.filter.return_value = ['1ubjA']
        mock_search_dali.return_value = mock_rec

        result = calcSignatureCavities(ref_structure='1ubiA',
                                        output_path=self.out_dir)

        mock_search_dali.assert_called_once_with('1ubi', 'A')
        mock_rec.filter.assert_called_once_with(
            cutoff_len=0.5, cutoff_rmsd=2.0, cutoff_Z=8, stringency=True)
        self.assertIn('1ubjA', result['msa_mappings'])

    @patch('prody.database.searchDali')
    def testDaliFilterKwargsOverridesDefaultCutoffs(self, mock_search_dali):
        """dali_filter_kwargs must be merged over the default cutoffs and
        passed through to DaliRecord.filter, letting a caller loosen or
        tighten the DALI search without editing the library."""

        mock_rec = MagicMock()
        mock_rec.isSuccess = True
        mock_rec.filter.return_value = ['1ubjA']
        mock_search_dali.return_value = mock_rec

        calcSignatureCavities(ref_structure='1ubiA',
                               output_path=self.out_dir,
                               dali_filter_kwargs={'cutoff_rmsd': 1.0})

        mock_rec.filter.assert_called_once_with(
            cutoff_len=0.5, cutoff_rmsd=1.0, cutoff_Z=8, stringency=True)

    @patch('prody.database.searchDali')
    def testDaliFilterKwargsHasNoEffectWithExplicitStructures(self, mock_search_dali):
        """dali_filter_kwargs must be ignored -- and DALI never invoked at
        all -- when an explicit structures list is given."""

        calcSignatureCavities(structures=['1ubjA'], ref_structure='1ubiA',
                               output_path=self.out_dir,
                               dali_filter_kwargs={'cutoff_rmsd': 1.0})

        mock_search_dali.assert_not_called()

    def testDaliRetryLoopFetchesUntilSuccess(self):
        """The isSuccess/fetch retry loop must keep calling fetch() until the
        record reports success before ever calling filter()."""

        fake_rec = _FakeDaliRecord(success_after=2, filter_result=['1ubjA'])
        with patch('prody.database.searchDali', return_value=fake_rec):
            calcSignatureCavities(ref_structure='1ubiA', output_path=self.out_dir)

        self.assertEqual(fake_rec.fetch_calls, 2)
        self.assertEqual(len(fake_rec.filter_calls), 1)

    @patch('prody.proteins.channels._warn')
    def testDaliZeroHitsReturnsNoneAndWarnsWithFilterSettings(self, mock_warn):
        """When filter() returns no hits, calcSignatureCavities must return
        None rather than proceed with a degenerate ensemble, and the warning
        must name the actual cutoffs used so the caller can act on it."""

        fake_rec = _FakeDaliRecord(success_after=0, filter_result=[])
        with patch('prody.database.searchDali', return_value=fake_rec):
            result = calcSignatureCavities(ref_structure='1ubiA',
                                            output_path=self.out_dir)

        self.assertIsNone(result)
        self.assertTrue(mock_warn.called)
        warned_msg = mock_warn.call_args[0][0]
        self.assertIn('cutoff_len=0.5', warned_msg)
        self.assertIn('cutoff_rmsd=2.0', warned_msg)
        self.assertIn('cutoff_Z=8', warned_msg)


class TestCalcSignatureCavitiesUserSuppliedMsa(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='sigcav_msa_')
        unaligned_dir = os.path.join(self.out_dir, 'unaligned')
        os.makedirs(unaligned_dir, exist_ok=True)
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubi')
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubj')

    def tearDown(self):
        shutil.rmtree(self.out_dir, ignore_errors=True)

    def _refSequence(self):
        ag = parsePDB(pathDatafile('pdb1ubi.pdb'))
        ca = ag.select('protein and name CA')
        return ca.getSequence()

    def testExactLabelMatchUsesSuppliedMsaForEntropy(self):
        """When the user's FASTA has a row exactly labeled ref_structure,
        that row anchors entropy directly -- no best-effort warning is
        needed, and the auto pseudo-MSA build/write is skipped entirely."""

        ref_seq = self._refSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>1ubiA\n{0}\n'.format(ref_seq))
            f.write('>unrelated_seq\n{0}\n'.format(ref_seq))

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir, msa_fasta=msa_path)

        self.assertIsNotNone(result)
        self.assertNotIn('msa_fasta_path', result)
        self.assertFalse(os.path.exists(
            os.path.join(self.out_dir, 'signature_cavity_msa.fasta')))

        ref_written = parsePDB(os.path.join(self.out_dir, 'aligned', 'ref.pdb'))
        betas = ref_written.select('protein and name CA').getBetas()
        self.assertTrue((betas < 1e-6).all())

    @patch('prody.proteins.channels._warn')
    def testBestEffortMatchWarnsWhenNoExactLabel(self, mock_warn):
        """When no row is labeled exactly ref_structure, the best-matching
        row by sequence identity against the reference is used instead,
        and a warning is emitted since the match was inferred."""

        ref_seq = self._refSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>some_other_label\n{0}\n'.format(ref_seq))

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir, msa_fasta=msa_path)

        self.assertIsNotNone(result)
        self.assertTrue(mock_warn.called)

        ref_written = parsePDB(os.path.join(self.out_dir, 'aligned', 'ref.pdb'))
        betas = ref_written.select('protein and name CA').getBetas()
        self.assertTrue((betas < 1e-6).all())

    def testBestEffortMatchExceedingMismatchThresholdRaises(self):
        """When even the best-matching row has >30% mismatch against the
        reference, calcSignatureCavities must raise ValueError rather than
        silently colouring from an untrustworthy alignment."""

        ref_seq = self._refSequence()
        unrelated_seq = 'A' * len(ref_seq)
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>some_other_label\n{0}\n'.format(unrelated_seq))

        with self.assertRaises(ValueError) as ctx:
            calcSignatureCavities(
                structures=['1ubjA'], ref_structure='1ubiA',
                output_path=self.out_dir, msa_fasta=msa_path)

        self.assertIn('mismatch', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
