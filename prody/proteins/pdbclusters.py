# -*- coding: utf-8 -*-
"""This module defines functions for handling PDB sequence clusters."""

import os.path

from numpy import array, abs

from prody import LOGGER, SETTINGS, getPackagePath
from prody.utilities import openFile, openURL, pystr, isListLike

from numbers import Integral

__all__ = ['fetchPDBClusters', 'loadPDBClusters', 'listPDBCluster']

PDB_CLUSTERS = {30: None, 40: None, 50: None, 70: None,
                90: None, 95: None, 100: None}
PDB_CLUSTERS_UPDATE_WARNING = True
PDB_CLUSTERS_SQIDS = array(list(PDB_CLUSTERS))
PDB_CLUSTERS_SQIDS.sort()
PDB_CLUSTERS_SQID_STR = ', '.join([str(key) for key in PDB_CLUSTERS_SQIDS])

def loadPDBClusters(sqid=None):
    """Load previously fetched PDB sequence clusters from disk to memory."""

    PDB_CLUSTERS_PATH = os.path.join(getPackagePath(), 'pdbclusters')
    if sqid is None:
        sqid_list = list(PDB_CLUSTERS)
        LOGGER.info('Loading all PDB sequence clusters.')
    else:
        assert isinstance(sqid, Integral), 'sqid must be an integer'
        if sqid not in PDB_CLUSTERS:
            raise ValueError('PDB cluster data is not available for sequence '
                             'identity {0}%, try one of {1}'
                             .format(sqid, PDB_CLUSTERS_SQID_STR))
        LOGGER.info('Loading PDB sequence clusters for sequence identity '
                    '{0}.'.format(sqid))
        sqid_list = [sqid]
    global PDB_CLUSTERS_UPDATE_WARNING
    for sqid in sqid_list:
        filename = os.path.join(PDB_CLUSTERS_PATH,
                                'bc-{0}.out.gz'.format(sqid))
        if not os.path.isfile(filename):
            fetchPDBClusters(sqid)
        if not os.path.isfile(filename):
            raise IOError('PDB sequence clusters for {0}% sequence identity '
                          'could not be downloaded.'.format(sqid))

        if PDB_CLUSTERS_UPDATE_WARNING:
            import time
            diff = (time.time() - os.path.getmtime(filename)) / 604800.
            if diff > 1.:
                LOGGER.warning('PDB sequence clusters are {0:.1f} week(s) old,'
                               ' call `fetchPDBClusters` to receive updates.'
                               .format(diff))
                PDB_CLUSTERS_UPDATE_WARNING = False
        inp = openFile(filename)
        clusters_str = pystr(inp.read())

        clusters = []
        for cluster_str in clusters_str.split('\n'):
            cluster_str = cluster_str.strip()
            if len(cluster_str):
                # split on the last underscore only: identifiers may themselves
                # contain underscores (e.g. computed structure models such as
                # ``AF_AFP12345F1``), so each member is an (identifier, entity)
                # pair where *entity* is the trailing entity number.
                cluster = [tuple(item.rsplit('_', 1))
                           for item in cluster_str.split()]
                clusters.append(cluster)

        PDB_CLUSTERS[sqid] = clusters
        inp.close()

    if sqid is None:
        return PDB_CLUSTERS
    else:
        return clusters


def listPDBCluster(pdb, entity=1, sqid=95):
    """Returns the PDB sequence cluster that contains polymer *entity* in
    structure *pdb* for sequence identity level *sqid*.  The cluster is
    returned as a list of tuples, e.g. ``[('1XXX', '1'), ]``.  Note that RCSB
    now clusters by polymer *entity* (not by individual chain), so *entity* is
    an entity identifier such as ``1`` and the same PDB identifier may appear
    more than once across clusters if the structure has multiple entities.

    Before this function is used, :func:`fetchPDBClusters` needs to be called.
    This function will load the PDB sequence clusters for *sqid* automatically
    using :func:`loadPDBClusters`."""

    assert isinstance(pdb, str) and len(pdb) == 4, \
        'pdb must be 4 char long string'
    entity = str(entity)
    try:
        sqid = int(sqid)
    except TypeError:
        raise TypeError('sqid must be an integer')
    if not (30 <= sqid <= 100):
        raise ValueError('sqid must be between 30 and 100')
    sqid = PDB_CLUSTERS_SQIDS[abs(PDB_CLUSTERS_SQIDS-sqid).argmin()]
    clusters = PDB_CLUSTERS[sqid]
    if clusters is None:
        loadPDBClusters(sqid)
        clusters = PDB_CLUSTERS[sqid]
    pdb_entity = (pdb.upper(), entity)

    for cluster in clusters:
        if pdb_entity in cluster:
            return cluster
    return

def fetchPDBClusters(sqid=None):
    """Retrieve PDB sequence clusters.  PDB sequence clusters are results of
    the regular clustering of polymer entities in the PDB.  They are available
    at: https://cdn.rcsb.org/resources/sequence/clusters/ as
    :file:`clusters-by-entity-{sqid}.txt` (the legacy :file:`bc-{sqid}.out`
    files have been deprecated by RCSB).

    This function will download the cluster data and save it after compressing
    in your home directory in :file:`.prody/pdbclusters`.  Cluster data can be
    loaded using :func:`loadPDBClusters` function and be accessed using
    :func:`listPDBCluster`."""

    if sqid is not None:
        if isListLike(sqid):
            for s in sqid:
                if s not in PDB_CLUSTERS:
                    raise ValueError('sqid must be one or more of ' + PDB_CLUSTERS_SQID_STR)
            keys = list(sqid)
        else:
            if sqid not in PDB_CLUSTERS:
                raise ValueError('sqid must be one or more of ' + PDB_CLUSTERS_SQID_STR)
            keys = [sqid]
    else:
        keys = list(PDB_CLUSTERS)

    PDB_CLUSTERS_PATH = os.path.join(getPackagePath(), 'pdbclusters')
    if not os.path.isdir(PDB_CLUSTERS_PATH):
        os.mkdir(PDB_CLUSTERS_PATH)
    LOGGER.progress('Downloading sequence clusters', len(keys),
                    '_prody_fetchPDBClusters')
    count = 0
    for i, x in enumerate(keys):
        # RCSB deprecated the legacy bc-*.out files; the current data is served
        # as clusters-by-entity-*.txt (clustered by polymer entity, not chain).
        urlname = 'clusters-by-entity-{0}.txt'.format(x)
        url = ('https://cdn.rcsb.org/resources/sequence/clusters/' + urlname)
        try:
            inp = openURL(url)
        except IOError:
            LOGGER.warning('Clusters at {0}% sequence identity level could '
                           'not be downloaded.'.format(x))
            continue
        else:
            data = inp.read()
            inp.close()
            # Guard against empty bodies or HTML error pages (e.g. a 404 that
            # is returned with a 200 status) being saved as valid cluster data
            # and reported as a successful download.
            if not data or pystr(data[:200]).lstrip().startswith('<'):
                LOGGER.warning('Clusters at {0}% sequence identity level could '
                               'not be downloaded.'.format(x))
                continue
            out = openFile('bc-{0}.out.gz'.format(x), 'w',
                           folder=PDB_CLUSTERS_PATH)
            out.write(data)
            out.close()
            count += 1
        LOGGER.update(i, label='_prody_fetchPDBClusters')
    LOGGER.finish()
    if len(keys) == count:
        LOGGER.info('All selected PDB clusters were downloaded successfully.')
    elif count == 0:
        LOGGER.warn('PDB clusters could not be downloaded.')
