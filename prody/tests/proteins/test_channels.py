"""This module contains unit tests for signature cavity/channel detection in
:mod:`~prody.proteins.channels`."""

import json
import os
import shutil
import tempfile
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

from prody import LOGGER, parsePDB, writePDB
from prody.proteins import channels
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


class TestCalcSignatureSelection(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='sigsel_')
        unaligned_dir = os.path.join(self.out_dir, 'unaligned')
        os.makedirs(unaligned_dir, exist_ok=True)
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubi')
        # 1ubj is 1ubi rigidly moved, waters included, so a kept hetero
        # group only matches the reference if the superposition carries it.
        moved = parsePDB(pathDatafile('pdb1ubi.pdb'))
        rot = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        moved.setCoords(moved.getCoords().dot(rot.T) + np.array([5., -3., 2.]))
        writePDB(os.path.join(unaligned_dir, 'pdb1ubj.pdb'), moved)

    def tearDown(self):
        shutil.rmtree(self.out_dir, ignore_errors=True)

    def _aligned(self, name):
        return parsePDB(os.path.join(self.out_dir, 'aligned', name))

    def testDefaultSelectionExcludesHeteroAtoms(self):
        """Without selection, only protein atoms reach the aligned and
        signature files, as before the parameter existed."""

        result = calcSignatureChannels(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        for ag in (self._aligned('1ubjA.pdb'), self._aligned('ref.pdb'),
                   parsePDB(result['signature_pdb'])):
            self.assertIsNone(ag.select('water'))
            self.assertIsNotNone(ag.select('protein'))

    def testSelectionKeepsHeteroAtomsInAlignedFrame(self):
        """A selection that keeps waters writes them to the aligned and
        signature files, superposed with the protein, while the reference
        entropy betas stay on protein atoms only."""

        result = calcSignatureChannels(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir, selection='protein or water')

        original_waters = parsePDB(pathDatafile('pdb1ubi.pdb')).select('water')
        ref_written = self._aligned('ref.pdb')
        mobile_written = self._aligned('1ubjA.pdb')

        for ag in (ref_written, mobile_written, parsePDB(result['signature_pdb'])):
            self.assertEqual(len(ag.select('water')), len(original_waters))

        self.assertLess(np.abs(
            mobile_written.select('water').getCoords()
            - ref_written.select('water').getCoords()).max(), 0.05)

        self.assertTrue((ref_written.select('protein').getBetas() < 1e-6).all())
        np.testing.assert_allclose(ref_written.select('water').getBetas(),
                                    original_waters.getBetas())

    def testSelectionMatchingNothingInReferenceRaises(self):
        """A selection that picks no reference atoms must fail loudly."""

        with self.assertRaises(ValueError) as ctx:
            calcSignatureChannels(
                structures=['1ubjA'], ref_structure='1ubiA',
                output_path=self.out_dir, selection='resname HEM')

        self.assertIn('selection', str(ctx.exception))


class TestCalcSignatureEntropyFile(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='sigent_')
        unaligned_dir = os.path.join(self.out_dir, 'unaligned')
        os.makedirs(unaligned_dir, exist_ok=True)
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubi')
        _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', '1ubj')
        self.ref_ca = parsePDB(pathDatafile('pdb1ubi.pdb')).select(
            'protein and name CA')

    def tearDown(self):
        shutil.rmtree(self.out_dir, ignore_errors=True)

    def _checkFile(self, result, kind):
        path = os.path.join(self.out_dir, 'signature_{0}_entropy.csv'.format(kind))
        self.assertEqual(os.path.realpath(result['entropy_path']),
                          os.path.realpath(path))
        on_disk = pd.read_csv(path)
        self.assertEqual(list(on_disk.columns), ['resnum', 'resname', 'entropy'])
        self.assertEqual(len(on_disk), self.ref_ca.numAtoms())
        self.assertEqual(list(on_disk['resnum']), list(self.ref_ca.getResnums()))
        self.assertEqual(list(on_disk['resname']), list(self.ref_ca.getResnames()))
        pd.testing.assert_frame_equal(result['entropy'], on_disk)
        return on_disk

    def testEnsembleMsaEntropyIsWrittenForEveryReferenceResidue(self):
        result = calcSignatureChannels(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        on_disk = self._checkFile(result, 'channel')
        self.assertTrue((on_disk['entropy'] < 1e-6).all())

    def testUncoveredResiduesAreNaN(self):
        """Residues beyond the end of the user's MSA row get NaN entropy."""

        covered = 70
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>1ubiA\n{0}\n'.format(self.ref_ca.getSequence()[:covered]))
            f.write('>other\n{0}\n'.format(self.ref_ca.getSequence()[:covered]))

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir, msa_fasta=msa_path)

        on_disk = self._checkFile(result, 'cavity')
        self.assertTrue(on_disk['entropy'][:covered].notna().all())
        self.assertTrue(on_disk['entropy'][covered:].isna().all())


class TestCalcSignatureRunSummary(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='sigsum_')
        unaligned_dir = os.path.join(self.out_dir, 'unaligned')
        os.makedirs(unaligned_dir, exist_ok=True)
        for pdb_id in ('1ubi', '1ubj', '1ubk'):
            _seedLocalStructure(unaligned_dir, 'pdb1ubi.pdb', pdb_id)

    def tearDown(self):
        shutil.rmtree(self.out_dir, ignore_errors=True)

    def _run(self, detect_channels, n_found=1, **kwargs):
        """Run with detection stubbed: *detect_channels* is called in place of
        calcChannels and may raise. *n_found* rows come back from the
        parameter and residue extraction."""

        with patch.object(channels, 'calcChannels', detect_channels), \
                patch.object(channels, 'getChannelParameters',
                             return_value=([10.0] * n_found, [1.5] * n_found,
                                           [100.0] * n_found)), \
                patch.object(channels, 'getChannelResidueNames',
                             return_value=['0: A1'] * n_found):
            result = calcSignatureChannels(
                structures=['1ubjA', '1ubkA'], ref_structure='1ubiA',
                output_path=self.out_dir, **kwargs)
        with open(result['run_summary_path']) as f:
            return result, json.load(f)

    def testSummaryRecordsReferenceParametersMsaAndStructures(self):
        result, summary = self._run(MagicMock(return_value=(['ch'], {})),
                                    selection='protein', voxel_res=0.75)

        self.assertEqual(
            os.path.realpath(result['run_summary_path']),
            os.path.realpath(os.path.join(
                self.out_dir, 'signature_channel_run_summary.json')))
        self.assertEqual(summary['kind'], 'channel')

        ref_seq = parsePDB(pathDatafile('pdb1ubi.pdb')).select(
            'protein and name CA').getSequence()
        self.assertEqual(summary['reference'], {
            'label': '1ubiA', 'sequence': ref_seq,
            'msa_anchor_row': None, 'msa_anchor_source': None})

        self.assertEqual(summary['parameters'], {
            'selection': 'protein', 'voxel_res': 0.75, 'threshold': 20,
            'distA': 4.5, 'dali_filter_kwargs': None, 'detection_kwargs': {}})
        self.assertEqual(summary['msa']['source'], 'ensemble')
        self.assertEqual(os.path.realpath(summary['msa']['path']),
                          os.path.realpath(result['msa_fasta_path']))

        self.assertEqual([s['label'] for s in summary['structures']],
                          ['1ubjA', '1ubkA'])
        for entry in summary['structures']:
            self.assertEqual(
                sorted(entry), ['aligned', 'coverage_pct', 'detection',
                                 'identity_pct', 'label', 'message',
                                 'n_objects', 'rmsd_A'])

    def testDetectionStatusOk(self):
        _, summary = self._run(MagicMock(return_value=(['ch'], {})))

        for entry in summary['structures']:
            self.assertTrue(entry['aligned'])
            self.assertEqual(entry['detection'], 'ok')
            self.assertEqual(entry['n_objects'], 1)
            self.assertIsNotNone(entry['rmsd_A'])
            self.assertIsNotNone(entry['coverage_pct'])
            self.assertIsNotNone(entry['identity_pct'])

    def testDetectionStatusNoneFoundWhenNoRowsComeBack(self):
        _, summary = self._run(MagicMock(return_value=([], {})), n_found=0)

        for entry in summary['structures']:
            self.assertTrue(entry['aligned'])
            self.assertEqual(entry['detection'], 'none_found')
            self.assertEqual(entry['n_objects'], 0)

    def testDetectionStatusErrorCarriesExceptionMessage(self):
        _, summary = self._run(MagicMock(side_effect=RuntimeError('boom')))

        for entry in summary['structures']:
            self.assertTrue(entry['aligned'])
            self.assertEqual(entry['detection'], 'error')
            self.assertEqual(entry['message'], 'boom')
            self.assertEqual(entry['n_objects'], 0)

    def testStructureThatFailsAlignmentIsSkippedWithReason(self):
        real_map = channels.mapChainOntoChain
        calls = []

        def failSecond(mobile, ref, **kwargs):
            calls.append(mobile)
            return None if len(calls) == 2 else real_map(mobile, ref, **kwargs)

        with patch.object(channels, 'mapChainOntoChain', failSecond):
            _, summary = self._run(MagicMock(return_value=(['ch'], {})))

        good, bad = summary['structures']
        self.assertTrue(good['aligned'])
        self.assertEqual(good['detection'], 'ok')
        self.assertFalse(bad['aligned'])
        self.assertEqual(bad['detection'], 'skipped')
        self.assertEqual(bad['n_objects'], 0)
        self.assertIn('CE-align could not map 1ubkA', bad['message'])
        self.assertIsNone(bad['rmsd_A'])

    def testUserMsaAnchorAndKwargsAreRecorded(self):
        ref_seq = parsePDB(pathDatafile('pdb1ubi.pdb')).select(
            'protein and name CA').getSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>1ubiA\n{0}\n>other\n{0}\n'.format(ref_seq))

        _, summary = self._run(MagicMock(return_value=(['ch'], {})),
                               msa_fasta=msa_path, dali_filter_kwargs={'cutoff_Z': 5},
                               probe_in=1.7)

        self.assertEqual(summary['msa'], {'source': 'user', 'path': msa_path})
        self.assertEqual(summary['reference']['msa_anchor_row'], '1ubiA')
        self.assertEqual(summary['reference']['msa_anchor_source'], 'label')
        self.assertEqual(summary['parameters']['detection_kwargs'],
                          {'probe_in': 1.7})
        self.assertEqual(summary['parameters']['dali_filter_kwargs'],
                          {'cutoff_Z': 5})

    def testNoImageFilesAreWritten(self):
        calcSignatureChannels(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir)

        written = [name for _, _, names in os.walk(self.out_dir)
                   for name in names if name.lower().endswith('.png')]
        self.assertEqual(written, [])


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

    @patch('prody.proteins.channels._warn')
    def testUniProtAccessionMatchAnchorsSuppliedMsa(self, mock_warn):
        """A Pfam-style row labeled by UniProt accession (not by
        ref_structure) must still anchor entropy directly, resolved via the
        reference PDB's own DBREF record (1UBI chain A -> P62988) -- this is
        the real-world Pfam case, where the label is never a pdbid+chain
        string, and no exact-label branch can ever fire. No best-effort
        warning should be emitted, since the match isn't a guess."""

        ref_seq = self._refSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>P62988.1/1-76\n{0}\n'.format(ref_seq))
            f.write('>UBIQ_UNRELATED/1-76\n{0}\n'.format(ref_seq))

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir, msa_fasta=msa_path)

        self.assertIsNotNone(result)
        self.assertFalse(mock_warn.called)

        ref_written = parsePDB(os.path.join(self.out_dir, 'aligned', 'ref.pdb'))
        betas = ref_written.select('protein and name CA').getBetas()
        self.assertTrue((betas < 1e-6).all())

    def testMsaRefLabelExplicitOverrideAnchorsSuppliedMsa(self):
        """msa_ref_label lets the caller name the anchor row directly,
        bypassing exact-label/accession/best-effort search entirely --
        useful when the caller already resolved the right row themselves
        (e.g. via a prior searchPfam call)."""

        ref_seq = self._refSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>totally_unrelated_label/1-76\n{0}\n'.format(ref_seq))
            f.write('>filler/1-76\n{0}\n'.format(ref_seq))

        result = calcSignatureCavities(
            structures=['1ubjA'], ref_structure='1ubiA',
            output_path=self.out_dir, msa_fasta=msa_path,
            msa_ref_label='totally_unrelated_label/1-76')

        self.assertIsNotNone(result)

        ref_written = parsePDB(os.path.join(self.out_dir, 'aligned', 'ref.pdb'))
        betas = ref_written.select('protein and name CA').getBetas()
        self.assertTrue((betas < 1e-6).all())

    def testMsaRefLabelNotFoundRaises(self):
        """An msa_ref_label that names no row in msa_fasta must raise
        ValueError rather than silently falling through to another
        matching strategy."""

        ref_seq = self._refSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            f.write('>some_label/1-76\n{0}\n'.format(ref_seq))
            f.write('>filler/1-76\n{0}\n'.format(ref_seq))

        with self.assertRaises(ValueError) as ctx:
            calcSignatureCavities(
                structures=['1ubjA'], ref_structure='1ubiA',
                output_path=self.out_dir, msa_fasta=msa_path,
                msa_ref_label='does_not_exist')

        self.assertIn('msa_ref_label', str(ctx.exception))

    def testOversizedUnmatchedMsaRaisesInsteadOfScanningEveryRow(self):
        """When no row is labeled ref_structure, no msa_ref_label is given,
        and no UniProt accession match is found, calcSignatureCavities must
        refuse to pairwise-align the reference against every row once the
        MSA is too large to search exhaustively (500+ rows) -- this is the
        scaling failure a real 'full' Pfam family alignment would otherwise
        hit, so it must fail fast with a clear message instead of hanging."""

        ref_seq = self._refSequence()
        msa_path = os.path.join(self.out_dir, 'user.fasta')
        with open(msa_path, 'w') as f:
            for i in range(501):
                f.write('>UNRELATED{0}_SPECIES/1-76\n{1}\n'.format(i, ref_seq))

        with self.assertRaises(ValueError) as ctx:
            calcSignatureCavities(
                structures=['1ubjA'], ref_structure='1ubiA',
                output_path=self.out_dir, msa_fasta=msa_path)

        self.assertIn('too many to search', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
