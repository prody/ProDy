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

"""This module defines functions for handling PDB sequence clusters."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

from numpy import array, abs

from prody import LOGGER, SETTINGS, getPackagePath
from prody.utilities import openFile

__all__ = ['fetchPDBClusters', 'loadPDBClusters', 'listPDBCluster', 
           'getPDBCluster',]

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
        assert isinstance(sqid, int), 'sqid must be an integer' 
        if sqid not in PDB_CLUSTERS:
            raise ValueError('PDB cluster data is not available for sequence '
                             'identity {0:d}%, try one of {1:s}'
                             .format(sqid, PDB_CLUSTERS_SQID_STR))
        LOGGER.info('Loading PDB sequence clusters for sequence identity '
                    '{0:d}.'.format(sqid))
        sqid_list = [sqid]
    global PDB_CLUSTERS_UPDATE_WARNING
    for sqid in sqid_list:
        filename = os.path.join(PDB_CLUSTERS_PATH, 
                                'bc-{0:d}.out.gz'.format(sqid))
        if not os.path.isfile(filename):
            fetchPDBClusters(sqid)
            
        if PDB_CLUSTERS_UPDATE_WARNING:
            import time
            diff = (time.time() - os.path.getmtime(filename)) / 604800.
            if diff > 1.:
                LOGGER.warning('PDB sequence clusters are {0:.1f} week(s) old,'
                               ' call `fetchPDBClusters` to receive updates.'
                               .format(diff))
                PDB_CLUSTERS_UPDATE_WARNING = False
        inp = openFile(filename)
        PDB_CLUSTERS[sqid] = inp.read()
        inp.close()

def getPDBCluster(pdb, ch, sqid=95):
    """Deprecated for removal in v1.4, use :func:`clustPDB` instead."""
    
    from prody import deprecate
    deprecate('getPDBCluster', 'listPDBCluster')
    
    return listPDB(pdb, ch, sqid)
    
def listPDBCluster(pdb, ch, sqid=95):
    """Return the PDB sequence cluster that contains chain *ch* in structure 
    *pdb* for sequence identity level *sqid*.  PDB sequence cluster will be 
    returned in as a list of tuples, e.g. ``[('1XXX', 'A'), ]``.  Note that 
    PDB clusters individual chains, so the same PDB identifier may appear 
    twice in the same cluster if the corresponding chain is present in the 
    structure twice.    
    
    Before this function is used, :func:`fetchPDBClusters` needs to be called. 
    This function will load the PDB sequence clusters for *sqid* automatically 
    using :func:`loadPDBClusters`."""

    assert isinstance(pdb, str) and len(pdb) == 4, \
        'pdb must be 4 char long string'
    assert isinstance(ch, str) and len(ch) == 1, \
        'ch must be a one char long string'
    try:
        sqid = int(sqid)
    except TypeError:
        raise TypeError('sqid must be an integer')
    if not (30 <= sqid <= 100):
        raise ValueError('sqid must be between 30 and 100')
    sqid = PDB_CLUSTERS_SQIDS[abs(PDB_CLUSTERS_SQIDS-sqid).argmin()]
    PDB_CLUSTERS_PATH = os.path.join(getPackagePath(), 'pdbclusters')
    clusters = PDB_CLUSTERS[sqid]
    if clusters is None: 
        loadPDBClusters(sqid)
        clusters = PDB_CLUSTERS[sqid]
    pdb_ch = pdb.upper() + '_' + ch.upper()
    index = clusters.index(pdb_ch)
    maxlen = clusters.index('\n') 
    end = clusters.find('\n', index)
    start = clusters.rfind('\n', index-maxlen, end)+1
    cluster = clusters[start:end]
    return [tuple(item.split('_')) for item in cluster.split()] 

def fetchPDBClusters(sqid=None):
    """Retrieve PDB sequence clusters.  PDB sequence clusters are results of 
    the weekly clustering of protein chains in the PDB generated by blastclust. 
    They are available at FTP site: ftp://resources.rcsb.org/sequence/clusters/
    
    This function will download about 10 Mb of data and save it after 
    compressing in your home directory in :file:`.prody/pdbclusters`.
    Compressed files will be less than 4 Mb in size.  Cluster data can 
    be loaded using :func:`loadPDBClusters` function and be accessed 
    using :func:`listPDBCluster`."""
    
    if sqid is not None:
        if sqid not in PDB_CLUSTERS:
            raise ValueError('sqid must be one of ' + PDB_CLUSTERS_SQID_STR)
        keys = [sqid]
    else:
        keys = list(PDB_CLUSTERS)
    
    import urllib2
    PDB_CLUSTERS_PATH = os.path.join(getPackagePath(), 'pdbclusters')
    if not os.path.isdir(PDB_CLUSTERS_PATH):
        os.mkdir(PDB_CLUSTERS_PATH)
    LOGGER.progress('Downloading sequence clusters', len(PDB_CLUSTERS),
                    '_prody_fetchPDBClusters')
    count = 0
    for i, x in enumerate(keys):
        filename = 'bc-{0:d}.out'.format(x)
        url = ('ftp://resources.rcsb.org/sequence/clusters/' + filename)
        try:
            inp = urllib2.urlopen(url)
        except urllib2.HTTPError:
            LOGGER.warning('Clusters at {0:d}% sequence identity level could '
                           'not be downloaded.')
            continue
        else:
            out = openFile(filename+'.gz', 'w', folder=PDB_CLUSTERS_PATH) 
            out.write(inp.read())
            inp.close()
            out.close()
            count += 1
        LOGGER.update(i, '_prody_fetchPDBClusters')
    LOGGER.clear()
    if len(PDB_CLUSTERS) == count:
        LOGGER.info('All PDB clusters were downloaded successfully.')
    elif count == 0:
        LOGGER.warn('PDB clusters could not be downloaded.')
