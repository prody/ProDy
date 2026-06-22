"""This module contains unit tests for :mod:`~prody.proteins.pdbclusters`."""

import os

from numpy.testing import *

from prody.utilities import importDec
dec = importDec()

from prody import fetchPDBClusters, loadPDBClusters, listPDBCluster
from prody import LOGGER, SETTINGS, getPackagePath
from prody.proteins import pdbclusters
from prody.tests import unittest, TEMPDIR

LOGGER.verbosity = 'none'

SQID = 40

# Cluster files are stored under getPackagePath(); when the package path is not
# set/writable (e.g. a fresh CI environment) prody prompts for it via input(),
# which fails under captured output ("reading from stdin while output is
# captured!").  Point the package path at the writable test TEMPDIR so the tests
# never block on a prompt; the original value is restored afterwards.
_ORIGINAL_PACKAGE_PATH = None


def setUpModule():
    global _ORIGINAL_PACKAGE_PATH
    _ORIGINAL_PACKAGE_PATH = SETTINGS.get('package_path', None)
    SETTINGS['package_path'] = TEMPDIR


def tearDownModule():
    SETTINGS['package_path'] = _ORIGINAL_PACKAGE_PATH


class TestFetchPDBClusters(unittest.TestCase):
    """Test downloading and querying the RCSB sequence clusters.  RCSB serves
    these as :file:`clusters-by-entity-{sqid}.txt` (the legacy
    :file:`bc-{sqid}.out` files were deprecated)."""

    @dec.slow
    def testFetchSavesUsableFile(self):

        fetchPDBClusters(SQID)
        filename = os.path.join(getPackagePath(), 'pdbclusters',
                                'bc-{0}.out.gz'.format(SQID))
        self.assertTrue(os.path.isfile(filename),
                        'fetchPDBClusters did not save a cluster file')
        self.assertGreater(os.path.getsize(filename), 100000,
                           'saved cluster file is implausibly small '
                           '(likely an error page, not cluster data)')

    @dec.slow
    def testListPDBClusterEntityFormat(self):

        fetchPDBClusters(SQID)
        # 101M is myoglobin; its single polymer entity (1) belongs to a large
        # cluster at low sequence identity.
        cluster = listPDBCluster('101M', 1, sqid=SQID)
        self.assertIsNotNone(cluster,
                             'listPDBCluster found no cluster for 101M entity 1')
        self.assertIn(('101M', '1'), cluster,
                      'query entity is missing from its own cluster')
        self.assertGreater(len(cluster), 1,
                           'cluster unexpectedly contains a single member')
        # Every member must be an (identifier, entity) pair, even for computed
        # structure models whose identifiers contain underscores (AF_..., MA_...).
        self.assertTrue(all(len(member) == 2 for member in cluster),
                        'cluster members are not (identifier, entity) pairs')
        # New data clusters by polymer entity, so the trailing field is an
        # entity number (e.g. "1"), not a chain letter (e.g. "A").
        self.assertTrue(all(ent.isdigit() for _, ent in cluster),
                        'cluster members are not in IDENTIFIER_ENTITY format')


class TestLoadPDBClustersFailure(unittest.TestCase):
    """A failed download must raise a clear error rather than crashing later in
    :func:`os.path.getmtime` with a :exc:`FileNotFoundError`."""

    def testFailedDownloadRaisesIOError(self):

        sqid = 50
        filename = os.path.join(getPackagePath(), 'pdbclusters',
                                'bc-{0}.out.gz'.format(sqid))
        orig_openURL = pdbclusters.openURL
        orig_clusters = pdbclusters.PDB_CLUSTERS[sqid]
        existed = os.path.isfile(filename)
        if existed:
            os.rename(filename, filename + '.bak')

        def _fail(*args, **kwargs):
            raise IOError('forced download failure')

        pdbclusters.openURL = _fail
        pdbclusters.PDB_CLUSTERS[sqid] = None
        try:
            self.assertRaises(IOError, loadPDBClusters, sqid)
        finally:
            pdbclusters.openURL = orig_openURL
            pdbclusters.PDB_CLUSTERS[sqid] = orig_clusters
            if existed:
                os.rename(filename + '.bak', filename)


if __name__ == '__main__':
    unittest.main()
