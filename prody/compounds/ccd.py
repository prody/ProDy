# -*- coding: utf-8 -*-
"""This module defines functions for fetching and parsing files in the PDB  
for the Chemical Compound Dictionary (CCD_). 

.. _CCD: https://www.wwpdb.org/data/ccd
"""

from prody import LOGGER, PY3K
from prody.utilities import openURL, isListLike
from prody.proteins.starfile import parseSTARLines, StarDict

__all__ = ['parseCCD']

def parseCCD(ids):
    """Retrieve the whole Chemical Component Dictionary (CCD) resource.
    """
    if isListLike(ids):
        n_ids = len(ids)
    else:
        ids = [ids]
        n_ids = 1

    ret = []
    for id in ids:
        id_url = 'http://ligand-expo.rcsb.org/reports/{0}/{1}/{1}.cif'.format(id[0],
                                                                              id)
        try:
            handle = openURL(id_url)
        except Exception as err:
            LOGGER.warn('download failed ({1}).'.format(str(err)))
        else:
            data = handle.read()
            if len(data):
                if PY3K:
                    data = data.decode()

                parsingDict, prog = parseSTARLines(data.split('\n'), shlex=True)
                        
                star_dict = StarDict(parsingDict, prog, id)
                ret.append(star_dict[id])
            else:
                ret.append(None)
                LOGGER.warn('Could not parse CCD data for {0}'.format(id))

    if n_ids == 1:
        return ret[0]

    return ret
