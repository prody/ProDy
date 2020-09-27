# -*- coding: utf-8 -*-
"""This module defines functions for fetching and parsing files in the PDB  
for the Chemical Compound Dictionary (CCD_). 

.. _CCD: https://www.wwpdb.org/data/ccd
"""

import os.path

from prody import LOGGER, SETTINGS, getPackagePath
from prody.utilities import openFile, openURL, pystr, isListLike
from prody.proteins import parseSTAR

__all__ = ['fetchCCDviaHTTP', 'parseCCD']

url = 'http://ligand-expo.rcsb.org/dictionaries/Components-pub.cif'
#ftp = 'ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz'

def fetchCCDviaHTTP(**kwargs):
    """Retrieve the whole Chemical Component Dictionary (CCD) resource.
    """
    CCD_PATH = os.path.join(getPackagePath(), 'ccd')

    success = 0
    failure = 0
    try:
        handle = openURL(url)
    except Exception as err:
        LOGGER.warn('download failed ({1}).'.format(str(err)))
        failure += 1
    else:
        data = handle.read()
        if len(data):
            if not os.path.isdir(CCD_PATH):
                os.mkdir(CCD_PATH)
            filename = CCD_PATH + '/Components-pub.cif'

            with open(filename, 'w+b') as outfile:
                outfile.write(data)

            success += 1
        else:
            failure += 1
    LOGGER.debug('CCD download via HTTP completed ({0} downloaded, '
                 '{1} failed).'.format(success, failure))

def parseCCD(*ids, **kwargs):
    """"Parse data from the Chemical Component Dictionary (CCD) resource

    :arg ids: one CCD identifier or a list of them.
        If **None** is provided then all of them are returned.
    :type ids: str, tuple, list, :class:`~numpy.ndarray`, **None**

    Returns :class:`.StarDataBlock` object or list of them.
    """
    CCD_PATH = os.path.join(getPackagePath(), 'ccd')
    
    n_ids = len(ids)
    if n_ids == 1:
        if isListLike(ids[0]):
            ids = ids[0]
            n_ids = len(ids)

    if n_ids == 1:
        ids = list(ids)

    filename = CCD_PATH + '/Components-pub.cif'
    if not os.path.isfile(filename):
        fetchCCDviaHTTP(**kwargs)

    data = parseSTAR(filename, shlex=True)
    ret = []
    for id in ids:
        try:
            ret.append(data[id])
        except ValueError:
            LOGGER.warn('id {0} not found in CCD data '
                        'so appending None'.format(id))
            ret.append(None)

    if n_ids == 1:
        return ret[0]

    return ret