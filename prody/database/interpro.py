# -*- coding: utf-8 -*-
"""This module defines functions for interfacing Interpro database."""

__author__ = 'James Krieger'

import json
from os.path import isfile
from prody import LOGGER, PY3K

__all__ = ['searchInterpro']

prefix = 'https://www.ebi.ac.uk/interpro/wwwapi/entry/'

def searchInterpro(query, **kwargs):
    """Returns Interpro search results in a list of dictionaries.  
    
    Matching family accessions as keys will map to various properties,
    including start and end residue positions.

    :arg query: UniProt ID or PDB identifier with or without a
        chain identifier, e.g. ``'1mkp'`` or ``'1mkpA'``.  
        UniProt ID of the specified chain, or the first
        protein chain will be used for searching the Pfam database
    :type query: str

    :arg timeout: timeout for blocking connection attempt in seconds, default
        is 60
    :type timeout: int
    """
    import requests
    
    LOGGER.timeit('_interpro')
    timeout = int(kwargs.get('timeout', 60))

    if len(query) == 4:
        url = prefix + "all/structure/pdb/" + query
        
    elif len(query) == 5:
        accession = None
        
        from prody import parsePDBHeader
        try:
            polymers = parsePDBHeader(query[:4], 'polymers')
        except Exception as err:
            raise ValueError('failed to parse header for {0} ({1})'
                             .format(query[:4], str(err)))

        chid = query[4:].upper()
        
        for poly in polymers:
            if chid and poly.chid != chid:
                continue
            for dbref in poly.dbrefs:
                if dbref.database != 'UniProt':
                    continue
                accession = dbref.accession
                LOGGER.info('UniProt accession {0} for {1} chain '
                            '{2} will be used.'
                            .format(accession, query[:4], poly.chid))
                break
            if accession is not None:
                break
            
        if accession is None:
            raise ValueError('A UniProt accession for PDB {0} could not be '
                             'parsed.'.format(repr(query)))
        else:
            url = prefix + "all/protein/uniprot/" + accession
        
    else:
        url = prefix + "all/protein/uniprot/" + query

    LOGGER.debug('Retrieving Interpro search results: ' + url)
    result = None
    sleep = 2
    while LOGGER.timing('_interpro') < timeout:
        try:
            result = requests.get(url, verify=False).content
        except Exception:
            pass
        else:
            if result not in ['PEND','RUN']:
                break
        
        sleep = 20 if int(sleep * 1.5) >= 20 else int(sleep * 1.5)
        LOGGER.sleep(int(sleep), '. Trying to reconnect...')

    if not result:
        raise IOError('Interpro search timed out or failed to parse results, '
                      ' check URL: ' + url)
    else:
        LOGGER.report('Interpro search completed in %.2fs.', '_interpro')

    if PY3K:
        result = result.decode()
        
    try:
        result = json.loads(result)
    except Exception as err:
        raise ValueError('failed to parse results as json, check URL: ' + url)

    return result["results"]
