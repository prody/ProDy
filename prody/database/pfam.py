# -*- coding: utf-8 -*-
"""This module defines functions for interfacing Pfam database."""

__author__ = 'Anindita Dutta, Ahmet Bakan, Cihan Kaya, James Krieger'

import re
from numbers import Integral

import numpy as np
from os.path import join, isfile
from io import BytesIO

from prody import LOGGER, PY3K
from prody.utilities import makePath, openURL, gunzip, openFile, dictElement
from prody.utilities import relpath
from prody.proteins import parsePDB

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

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

old_prefix = 'https://pfam.xfam.org/'
prefix = 'https://pfam-legacy.xfam.org/'
new_prefix = 'https://www.ebi.ac.uk/interpro/wwwapi/entry/'

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

    import xml.etree.cElementTree as ET
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
            url = new_prefix + "all/protein/uniprot/" + accession

    else:
        url = new_prefix + "all/protein/uniprot/" + seq

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

    if xml.find('There was a system error on your last request.') > 0:
        LOGGER.warn('No Pfam matches found for: ' + seq)
        return None
    elif xml.find('No valid UniProt accession or ID') > 0:
        try:
            url = prefix + 'protein/' + accession + '?output=xml'
            LOGGER.debug('Retrieving Pfam search results: ' + url)
            xml = openURL(url, timeout=timeout).read()
        except:
            raise ValueError('No valid UniProt accession or ID for: ' + seq)
        
        if xml.find('No valid UniProt accession or ID') > 0:
            try:
                ag = parsePDB(seq, subset='ca')
                ag_seq = ag.getSequence()
                return searchPfam(ag_seq)
            except:
                try:
                    url = 'https://uniprot.org/uniprot/' + accession + '.xml'
                    xml = openURL(url, timeout=timeout).read()
                    if len(xml) > 0:
                        root = ET.XML(xml)
                        accession = root[0][0].text

                        url = prefix + 'protein/' + accession + '?output=xml'
                        LOGGER.debug('Retrieving Pfam search results: ' + url)
                        xml = openURL(url, timeout=timeout).read()
                    else:
                        raise ValueError('No valid UniProt accession or ID for: ' + seq)
                except:
                    raise ValueError('No valid UniProt accession or ID for: ' + seq)

    try:
        root = json.loads(xml)
        #return root
    except Exception as err:
        raise ValueError('failed to parse results XML, check URL: ' + url)

    matches = dict()
    for entry in root["results"]:
        try:
            metadata = entry["metadata"]
            accession = metadata["accession"]
        except KeyError:
            raise ValueError('failed to parse accessions from results, check URL: ' + url)

        if not re.search('PF[0-9]{5}$', accession):
            continue

        match = matches.setdefault(accession, dict(metadata.items()))
        
        other_data = entry["proteins"]
        locations = match.setdefault("locations", [])
        for item1 in other_data:
            for key, value in item1.items():
                if key == "entry_protein_locations":
                    locations.append(value)    

    query = 'Query ' + repr(query)

    if matches:
        LOGGER.info(query + ' matched {0} Pfam families.'.format(len(matches)))
    else:
        LOGGER.info(query + ' did not match any Pfam families.')
    return matches


def fetchPfamMSA(acc, alignment='full', compressed=False, **kwargs):
    """Returns a path to the downloaded Pfam MSA file.

    :arg acc: Pfam ID or Accession Code
    :type acc: str

    :arg alignment: alignment type, one of ``'full'`` (default), ``'seed'``,
         ``'ncbi'``, ``'metagenomics'``, ``'rp15'``, ``'rp35'``, ``'rp55'``,
         ``'rp75'`` or ``'uniprot'`` where rp stands for representative 
         proteomes

    :arg compressed: gzip the downloaded MSA file, default is **False**

    *Alignment Options*

    :arg format: a Pfam supported MSA file format, one of ``'selex'``,
        (default), ``'stockholm'`` or ``'fasta'``

    :arg order: ordering of sequences, ``'tree'`` (default) or
        ``'alphabetical'``

    :arg inserts: letter case for inserts, ``'upper'`` (default) or ``'lower'``

    :arg gaps: gap character, one of ``'dashes'`` (default), ``'dots'``,
        ``'mixed'`` or **None** for unaligned

    *Other Options*

    :arg timeout: timeout for blocking connection attempt in seconds, default
        is 60

    :arg outname: out filename, default is input ``'acc_alignment.format'``

    :arg folder: output folder, default is ``'.'``"""
    
    import requests

    # url = prefix + 'family/acc?id=' + acc
    # handle = openURL(url, timeout=int(kwargs.get('timeout', 60)))
    orig_acc = acc
    # acc = handle.readline().strip()
    # if PY3K:
    #     acc = acc.decode()
    url_flag = False

    if not re.search('(?<=PF)[0-9]{5}$', acc):
        raise ValueError('{0} is not a valid Pfam ID or Accession Code'
                         .format(repr(orig_acc)))

    if alignment not in DOWNLOAD_FORMATS:
        raise ValueError('alignment must be one of full, seed,'
                         #' ncbi or'
                         #' metagenomics'
                         ' or uniprot')
    # if alignment == 'ncbi' or alignment == 'metagenomics' or alignment == 'uniprot':
    #     #url = (prefix + 'family/' + acc + '/alignment/' +
    #     #       alignment + '/gzipped')
    #     url = (new_prefix + acc + 
    #            '/?annotation=alignment:' + alignment + '&download')
    #     url_flag = True
    #     extension = '.sth'
    # else:
    if not kwargs:
        #url = (prefix + 'family/' + acc + '/alignment/' +
        #       alignment + '/gzipped')
        url = (new_prefix + "/pfam/" + acc + 
                '/?annotation=alignment:' + alignment + '&download')
        url_flag = True
        extension = '.sth'
    else:
        raise ValueError('kwargs are not supported for Interpro Pfam')
    #     align_format = kwargs.get('format', 'selex').lower()

    #     if align_format not in FORMAT_OPTIONS['format']:
    #         raise ValueError('alignment format must be of type selex'
    #                             ' stockholm or fasta. MSF not supported')

    #     if align_format == SELEX:
    #         align_format, extension = 'pfam', '.slx'
    #     elif align_format == FASTA:
    #         extension = '.fasta'
    #     else:
    #         extension = '.sth'

    #     gaps = str(kwargs.get('gaps', 'dashes')).lower()
    #     if gaps not in FORMAT_OPTIONS['gaps']:
    #         raise ValueError('gaps must be of type mixed, dots, dashes, '
    #                             'or None')

    #     inserts = kwargs.get('inserts', 'upper').lower()
    #     if(inserts not in FORMAT_OPTIONS['inserts']):
    #         raise ValueError('inserts must be of type lower or upper')

    #     order = kwargs.get('order', 'tree').lower()
    #     if order not in FORMAT_OPTIONS['order']:
    #         raise ValueError('order must be of type tree or alphabetical')

    #     url = (prefix + 'family/' + acc + '/alignment/'
    #             + alignment + '/format?format=' + align_format +
    #             '&alnType=' + alignment + '&order=' + order[0] +
    #             '&case=' + inserts[0] + '&gaps=' + gaps + '&download=1')

    LOGGER.timeit('_pfam')
    timeout = kwargs.get('timeout', 60)
    response = None
    sleep = 2
    try_error = 3
    while LOGGER.timing('_pfam') < timeout:
        try:
            response = requests.get(url, verify=False).content
        except Exception:
            pass
        else:
            break
        
        sleep = 20 if int(sleep * 1.5) >= 20 else int(sleep * 1.5)
        LOGGER.sleep(int(sleep), '. Trying to reconnect...')

    # response = openURL(url, timeout=int(kwargs.get('timeout', 60)))
    outname = kwargs.get('outname', None)
    if not outname:
        outname = orig_acc
    folder = str(kwargs.get('folder', '.'))
    filepath = join(makePath(folder), outname + '_' + alignment + extension)
    if compressed:
        filepath = filepath + '.gz'
        if url_flag:
            f_out = open(filepath, 'wb')
        else:
            f_out = openFile(filepath, 'wb')
        # f_out.write(response.read())
        f_out.write(response)
        f_out.close()
    else:
        if url_flag:
            gunzip(response, filepath)
        else:
            with open(filepath, 'wb') as f_out:
                # f_out.write(response.read())
                f_out.write(response)

    filepath = relpath(filepath)
    LOGGER.info('Pfam MSA for {0} is written as {1}.'
                .format(orig_acc, filepath))

    return filepath

def parsePfamPDBs(query, data=[], **kwargs):
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
    """
    
    start = kwargs.pop('start', 1)
    end = kwargs.pop('end', None)

    if len(query) > 4 and query.startswith('PF'):
        pfam_acc = query
    else:
        pfam_matches = searchPfam(query, **kwargs)
        keys = list(pfam_matches.keys())

        if isinstance(start, Integral):
            try:
                start_diff = []
                for i, key in enumerate(pfam_matches):
                    start_diff.append(int(pfam_matches[key]['locations'][0]['start']) - start)
                start_diff = np.array(start_diff)
                pfam_acc = keys[np.where(abs(start_diff) == min(abs(start_diff)))[0][0]]
            except KeyError:
                start_diff = []
                for i, key in enumerate(pfam_matches):
                    start_diff.append(int(pfam_matches[key]['locations']['ali_start']) - start)
                start_diff = np.array(start_diff)
                pfam_acc = keys[np.where(abs(start_diff) == min(abs(start_diff)))[0][0]]

        elif isinstance(end, Integral):
            end_diff = []
            for i, key in enumerate(pfam_matches):
                end_diff.append(int(pfam_matches[key]['locations'][0]['end']) - end)
            end_diff = np.array(end_diff)
            pfam_acc = keys[np.where(abs(end_diff) == min(abs(end_diff)))[0][0]]

        else:
            raise ValueError('Please provide an integer for start or end '
                             'when using a UniProt ID or PDB ID.')

    from ftplib import FTP
    from .uniprot import queryUniprot

    data_stream = BytesIO()
    ftp_host = 'ftp.ebi.ac.uk'
    ftp = FTP(ftp_host)
    ftp.login()
    ftp.cwd('pub/databases/Pfam/current_release')
    ftp.retrbinary('RETR pdbmap.gz', data_stream.write)
    ftp.quit()
    zip_data = data_stream.getvalue()
    data_stream.close()

    rawdata = gunzip(zip_data)
    if PY3K:
        rawdata = rawdata.decode()

    fields = ['PDB_ID', 'chain', 'nothing', 'PFAM_Name', 'PFAM_ACC', 
              'UniprotAcc', 'UniprotResnumRange']
    
    data_dicts = []
    for line in rawdata.split('\n'):
        if line.find(pfam_acc) != -1:
            data_dicts.append({})
            for j, entry in enumerate(line.strip().split('\t')):
                data_dicts[-1][fields[j]] = entry.strip(';')

    pdb_ids = [data_dict['PDB_ID'] for data_dict in data_dicts]
    chains = [data_dict['chain'] for data_dict in data_dicts]

    header = kwargs.pop('header', False)
    model = kwargs.get('model', None)
    results = parsePDB(pdb_ids, chain=chains, header=True, **kwargs)

    ags, headers = results
    ags, headers = list(ags), list(headers)

    if model == 0:
        LOGGER.info('only header is requested and returned')
        return results

    if header:
        results = (ags, headers)
    else:
        results = ags

    LOGGER.progress('Extracting Pfam domains...', len(ags))
    comma_splitter = re.compile(r'\s*,\s*').split
    no_info = []
    for i, ag in enumerate(ags):
        LOGGER.update(i)
        data_dict = data_dicts[i]
        pfamRange = data_dict['UniprotResnumRange'].split('-')
        uniprotAcc = data_dict['UniprotAcc']
        try:
            uniData = queryUniprot(uniprotAcc)
        except:
            LOGGER.warn('No Uniprot record found for {0}'.format(data_dict['PDB_ID']))
            continue

        resrange = None
        found = False
        for key, value in uniData.items():
            if not key.startswith('dbReference'):
                continue
            try:
                pdbid = value['PDB']
            except:
                continue
            if pdbid != data_dict['PDB_ID']:
                continue
            pdbchains = value['chains']

            # example chain strings: "A=27-139, B=140-150" or "A/B=27-150"
            pdbchains = comma_splitter(pdbchains)
            for chain in pdbchains:
                chids, resrange = chain.split('=')
                chids = [chid.strip() for chid in chids.split('/')]
                if data_dict['chain'] in chids:
                    resrange = resrange.split('-')
                    found = True
                    break
            if found:
                break

        if found:
            header = headers[i]
            chain_accessions = [dbref.accession 
                                for dbref in header[data_dict['chain']].dbrefs]
            try:
                if len(chain_accessions) > 0:
                    right_part = np.where(np.array(chain_accessions) == 
                                        data_dict['UniprotAcc'])[0][0]
                else:
                    raise ValueError('There is no accession for a chain in the Header')
            except:
                LOGGER.warn('Could not map domains in {0}'
                            .format(data_dict['PDB_ID'] 
                            + data_dict['chain']))
                no_info.append(i)
                continue

            right_dbref = header[data_dict['chain']].dbrefs[right_part]
            chainStart = ag.select('chain {0}'.format(data_dict['chain'])
                                  ).getResnums()[0]
            missing = chainStart - right_dbref.first[0]
            partStart = ag.getResindices()[np.where(ag.getResnums() == 
                                           right_dbref.first[0] + missing)][0]
            pfStart, pfEnd = int(pfamRange[0]), int(pfamRange[1])
            uniStart, uniEnd = int(resrange[0]), int(resrange[1])

            resiStart = pfStart - uniStart + partStart - missing
            resiEnd = pfEnd - uniStart + partStart - missing
            ags[i] = ag.select('resindex {0} to {1}'.format(
                            resiStart, resiEnd)) 
        else:
            no_info.append(i)
    LOGGER.finish()

    for i in reversed(no_info):
        ags.pop(i)
        if header:
            headers.pop(i)

    if isinstance(data, list):
        data.extend(data_dicts)
    else:
        LOGGER.warn('data should be a list in order to get output')
    
    return results

