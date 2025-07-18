# -*- coding: utf-8 -*-
"""This module defines functions for interfacing Pfam database."""

__author__ = 'Anindita Dutta, Ahmet Bakan, Cihan Kaya, James Krieger'

import re
from numbers import Integral

import numpy as np
from os.path import join
from io import BytesIO

from prody import LOGGER, PY3K
from prody.utilities import makePath, gunzip
from prody.utilities import relpath
from prody.proteins import parsePDB

import json

__all__ = ['searchPfam', 'fetchPfamMSA', 'parsePfamPDBs']

FASTA = 'fasta'
SELEX = 'selex'
STOCKHOLM = 'stockholm'

DOWNLOAD_FORMATS = set(['seed', 'full', 'uniprot', 
                        #'ncbi', 'metagenomics',
                        #'rp15', 'rp35', 'rp55', 'rp75'
                        ])
FORMAT_OPTIONS = ({'format': set([FASTA, SELEX, STOCKHOLM]),
                  'order': set(['tree', 'alphabetical']),
                  'inserts': set(['lower', 'upper']),
                  'gaps': set(['mixed', 'dots', 'dashes', 'none'])})

prefix = 'https://www.ebi.ac.uk/interpro/wwwapi/entry/'

def searchPfam(query, **kwargs):
    """Returns Pfam search results in a dictionary.  Matching Pfam accession
    as keys will map to evalue, alignment start and end residue positions.

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

    seq = ''.join(query.split())

    LOGGER.timeit('_pfam')
    timeout = int(kwargs.get('timeout', 60))

    if len(seq) <= 5:
        accession = None
        from prody import parsePDBHeader
        try:
            polymers = parsePDBHeader(seq[:4], 'polymers')
        except Exception as err:
            raise ValueError('failed to parse header for {0} ({1})'
                                .format(seq[:4], str(err)))
        else:
            chid = seq[4:].upper()

        for poly in polymers:
            if chid and poly.chid != chid:
                continue
            for dbref in poly.dbrefs:
                if dbref.database != 'UniProt':
                    continue
                accession = dbref.accession
                LOGGER.info('UniProt accession {0} for {1} chain '
                            '{2} will be used.'
                            .format(accession, seq[:4], poly.chid))
                break
            if accession is not None:
                break
        if accession is None:
            raise ValueError('A UniProt accession for PDB {0} could not be '
                                'parsed.'.format(repr(seq)))
        else:
            url = prefix + "all/protein/uniprot/" + accession

    elif seq.startswith('PF'):
        url = prefix + "pfam/" + seq
    else:
        url = prefix + "all/protein/uniprot/" + seq

    LOGGER.debug('Retrieving Pfam search results: ' + url)
    xml = None
    sleep = 2
    while LOGGER.timing('_pfam') < timeout:
        try:
            xml = requests.get(url, verify=False).content
        except Exception:
            pass
        else:
            if xml not in ['PEND','RUN']:
                break
        
        sleep = 20 if int(sleep * 1.5) >= 20 else int(sleep * 1.5)
        LOGGER.sleep(int(sleep), '. Trying to reconnect...')

    if not xml:
        raise IOError('Pfam search timed out or failed to parse results '
                      'XML, check URL: ' + url)
    else:
        LOGGER.report('Pfam search completed in %.2fs.', '_pfam')

    if PY3K:
        xml = xml.decode()
    else:
        xml = xml.encode()

    if xml.find('There was a system error on your last request.') > 0 or xml.find("Error") > 0:
        LOGGER.warn('No Pfam matches found for: ' + seq)
        return None
    elif xml.find('No valid UniProt accession or ID') > 0:
        raise ValueError('No valid UniProt accession or ID for: ' + seq)

    try:
        root = json.loads(xml)
    except Exception as err:
        raise ValueError('failed to parse results XML, check URL: ' + url)

    matches = dict()

    if seq.startswith('PF'):
        metadata = root['metadata']
        matches.setdefault(str(seq), dict(metadata.items()))
        return matches

    for entry in root["results"]:
        try:
            metadata = entry["metadata"]
            accession = metadata["accession"]
            if isinstance(metadata["member_databases"], dict):
                accessions2 = [list(value.keys())
                               for value in metadata["member_databases"].values()]
            else:
                accessions2 = []
        except KeyError:
            raise ValueError('failed to parse accessions from results, check URL: ' + url)

        pfamAccessions = []
        if re.search('PF[0-9]{5}$', accession):
            pfamAccessions.append(accession)
        else:
            for accession in np.array(accessions2).flatten():
                if (re.search('PF[0-9]{5}$', accession) and
                    entry["proteins"][0]["entry_protein_locations"] is not None):
                    pfamAccessions.append(accession)

        if len(pfamAccessions) == 0:
            continue

        for accession in pfamAccessions:
            match = matches.setdefault(str(accession), dict(metadata.items()))
            
            other_data = entry["proteins"]
            locations = match.setdefault("locations", [])
            for item1 in other_data:
                for key, value in item1.items():
                    if key == "entry_protein_locations":
                        for item2 in value:
                            new_dict = {}
                            for item3 in value[0]["fragments"]:
                                new_dict["start"] = item3["start"]
                                new_dict["end"] = item3["end"]
                                new_dict["score"] = item2["score"]
                                locations.append(new_dict)    

    query = 'Query ' + repr(query)

    if matches:
        LOGGER.info(query + ' matched {0} Pfam families.'.format(len(matches)))
    else:
        LOGGER.info(query + ' did not match any Pfam families.')
    return matches


def fetchPfamMSA(acc, alignment='seed', compressed=False, **kwargs):
    """Returns a path to the downloaded Pfam MSA file.

    :arg acc: Pfam ID or Accession Code
    :type acc: str

    :arg alignment: alignment type, one of ``'full'``, ``'seed'`` (default),
         ``'ncbi'``, ``'metagenomics'``, ``'rp15'``, ``'rp35'``, ``'rp55'``,
         ``'rp75'`` or ``'uniprot'`` where rp stands for representative 
         proteomes. InterPro Pfam seems to only have seed alignments
         easily accessible in most cases

    :arg compressed: gzip the downloaded MSA file, default is **False**

    :arg timeout: timeout for blocking connection attempt in seconds, default
        is 60

    :arg outname: out filename, default is input ``'acc_alignment.format'``
    :type outname: str

    :arg folder: output folder, default is ``'.'``
    :type folder: str
    """
    
    import requests

    if not re.search('(?<=PF)[0-9]{5}$', acc):
        raise ValueError('{0} is not a valid Pfam ID or Accession Code'
                         .format(repr(acc)))

    if alignment not in DOWNLOAD_FORMATS:
        raise ValueError('alignment must be one of full, seed, or uniprot')

    url = (prefix + "/pfam/" + acc + 
            '/?annotation=alignment:' + alignment + '&download')
    extension = '.sth'

    LOGGER.timeit('_pfam')
    timeout = kwargs.get('timeout', 60)
    response = None
    sleep = 2
    while LOGGER.timing('_pfam') < timeout:
        try:
            response = requests.get(url, verify=False).content
        except Exception:
            pass
        else:
            break
        
        sleep = 20 if int(sleep * 1.5) >= 20 else int(sleep * 1.5)
        LOGGER.sleep(int(sleep), '. Trying to reconnect...')

    outname = kwargs.get('outname', None)
    if not outname:
        outname = acc
    folder = str(kwargs.get('folder', '.'))
    filepath = join(makePath(folder), outname + '_' + alignment + extension)
    if compressed:
        filepath = filepath + '.gz'
        f_out = open(filepath, 'wb')
        f_out.write(response)
        f_out.close()
    else:
        gunzip(response, filepath)

    filepath = relpath(filepath)
    LOGGER.info('Pfam MSA for {0} is written as {1}.'
                .format(acc, filepath))

    return filepath

def parsePfamPDBs(query, data=None, **kwargs):
    """Returns a list of :class:`.AtomGroup` objects containing sections of chains 
    that correspond to a particular PFAM domain family. These are defined by 
    alignment start and end residue numbers.

    :arg query: Pfam ID, UniProt ID or PDB ID
        If a PDB ID is provided the corresponding UniProt ID is used.
        If this returns multiple matches then start or end must also be provided.
        This query is also used for label refinement of the Pfam domain MSA.
    :type query: str

    :arg data: If given the data list from the Pfam mapping table will 
        be output through this argument.
    :type data: list

    :keyword start: Residue number for defining the start of the domain.
        The PFAM domain that starts closest to this will be selected. 
        Default is **1**
    :type start: int

    :keyword end: Residue number for defining the end of the domain.
        The PFAM domain that ends closest to this will be selected. 
    :type end: int

    :keyword num_pdbs: Number of pdb chains to parse before stopping
        Default behaviour is to pass all of them
    :type num_pdbs: int
    """

    if data is None: data = []

    only_parse = kwargs.pop('only_parse', False)
    
    start = kwargs.pop('start', 1)
    end = kwargs.pop('end', None)
    num_pdbs = kwargs.pop('num_pdbs', None)

    if len(query) > 4 and query.startswith('PF'):
        pfam_acc = query
    else:
        if not isinstance(start, Integral) and not isinstance(end, Integral):
            raise ValueError('Please provide an integer for start or end '
                             'when using a UniProt ID or PDB ID.')

        pfam_matches = searchPfam(query, **kwargs)
        keys = list(pfam_matches.keys())

        if isinstance(start, Integral):
            try:
                start_diff = []
                for i, key in enumerate(pfam_matches):
                    start_diff.append(int(pfam_matches[key]['locations'][0]['start']) - start)
                start_diff = np.array(start_diff)
                pfam_acc = keys[np.nonzero(abs(start_diff) == min(abs(start_diff)))[0][0]]
            except KeyError:
                start_diff = []
                for i, key in enumerate(pfam_matches):
                    start_diff.append(int(pfam_matches[key]['locations']['ali_start']) - start)
                start_diff = np.array(start_diff)
                pfam_acc = keys[np.nonzero(abs(start_diff) == min(abs(start_diff)))[0][0]]

        if isinstance(end, Integral):
            end_diff = []
            for i, key in enumerate(pfam_matches):
                end_diff.append(int(pfam_matches[key]['locations'][0]['end']) - end)
            end_diff = np.array(end_diff)
            pfam_acc = keys[np.nonzero(abs(end_diff) == min(abs(end_diff)))[0][0]]

    from ftplib import FTP

    data_stream = BytesIO()
    ftp_host = 'ftp.ebi.ac.uk'
    ftp = FTP(ftp_host)
    ftp.login()
    ftp.cwd('pub/databases/Pfam/mappings')
    ftp.retrbinary('RETR pdb_pfam_mapping.txt', data_stream.write)
    ftp.quit()
    rawdata = data_stream.getvalue()
    data_stream.close()

    if PY3K:
        rawdata = rawdata.decode()

    fields = rawdata.split('\n')[1].split('\t')
    
    data_dicts = []
    for line in rawdata.split('\n'):
        if line.find(pfam_acc) != -1:
            data_dicts.append({})
            for j, entry in enumerate(line.strip().split('\t')):
                data_dicts[-1][fields[j]] = entry.strip(';')

    pdb_ids = [data_dict['PDB'] for data_dict in data_dicts][:num_pdbs]
    if len(pdb_ids) == 1:
        pdb_ids = pdb_ids[0]
    chains = [data_dict['CHAIN'] for data_dict in data_dicts][:num_pdbs]
    if len(chains) == 1:
        chains = chains[0]

    header = kwargs.pop('header', False)
    model = kwargs.get('model', None)
    results = parsePDB(pdb_ids, chain=chains, header=True,
                       unite_chains=True, **kwargs)

    ags, headers = results
    ags, headers = list(ags), list(headers)

    if model == 0:
        LOGGER.info('only header is requested and returned')
        return results

    if header:
        results = (ags, headers)
    else:
        results = ags

    if only_parse:
        return [result for result in results if result is not None]

    LOGGER.progress('Extracting Pfam domains...', len(ags))
    for i, ag in enumerate(ags):
        if ag is None:
            continue

        LOGGER.update(i)
        data_dict = data_dicts[i]
        ags[i] = ag.select('resnum {0} to {1}'.format(
                            int(data_dict['PDB_START']),
                            int(data_dict['PDB_END'])))
    LOGGER.finish()

    if isinstance(data, list):
        data.extend(data_dicts)
    else:
        LOGGER.warn('data should be a list in order to get output')
    
    return [result for result in results if result is not None]

