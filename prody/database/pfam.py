# -*- coding: utf-8 -*-
"""This module defines functions for interfacing Pfam database."""

__author__ = 'Anindita Dutta, Ahmet Bakan, Cihan Kaya'

import re
from numbers import Integral

import numpy as np
import os
from os.path import join, isfile
from io import BytesIO
import zlib

from prody import LOGGER, PY3K
from prody.utilities import makePath, openURL, gunzip, openFile, dictElement
from prody.utilities import relpath
from prody.proteins import parsePDB
from prody.sequence import parseMSA, refineMSA, MSA
from .uniprot import queryUniprot

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

__all__ = ['searchPfam', 'fetchPfamMSA', 'searchUniprotID', 'parsePfamPDBs']

FASTA = 'fasta'
SELEX = 'selex'
STOCKHOLM = 'stockholm'

DOWNLOAD_FORMATS = set(['seed', 'full', 'ncbi', 'metagenomics',
                        'rp15', 'rp35', 'rp55', 'rp75', 'uniprot'])
FORMAT_OPTIONS = ({'format': set([FASTA, SELEX, STOCKHOLM]),
                  'order': set(['tree', 'alphabetical']),
                  'inserts': set(['lower', 'upper']),
                  'gaps': set(['mixed', 'dots', 'dashes', 'none'])})

MINSEQLEN = 16

prefix = 'https://pfam.xfam.org/'

def searchPfam(query, **kwargs):
    """Returns Pfam search results in a dictionary.  Matching Pfam accession
    as keys will map to evalue, alignment start and end residue positions.

    :arg query: UniProt ID, PDB identifier, a protein sequence, or a sequence
        file. Sequence queries must not contain without gaps and must be at
        least 16 characters long
    :type query: str

    :arg timeout: timeout for blocking connection attempt in seconds, default
        is 60
    :type timeout: int

    *query* can also be a PDB identifier, e.g. ``'1mkp'`` or ``'1mkpA'`` with
    chain identifier.  UniProt ID of the specified chain, or the first
    protein chain will be used for searching the Pfam database."""

    if isfile(query):
        from prody.sequence import MSAFile
        try:
            seq = next(MSAFile(query))
        except:
            with openFile(query) as inp:
                seq = ''.join(inp.read().split())
        else:
            seq = seq[0][1]
        if not seq.isalpha():
            raise ValueError('could not parse a sequence without gaps from ' +
                             query)
    else:
        seq = ''.join(query.split())

    import xml.etree.cElementTree as ET
    LOGGER.timeit('_pfam')
    timeout = int(kwargs.get('timeout', 60))
    if len(seq) >= MINSEQLEN:
        if not seq.isalpha():
            raise ValueError(repr(seq) + ' is not a valid sequence')
        fseq = '>Seq\n' + seq
        parameters = { 'hmmdb' : 'pfam', 'seq': fseq }
        enc_params = urllib.urlencode(parameters).encode('utf-8')
        request = urllib2.Request('https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan', enc_params)

        results_url = urllib2.urlopen(request).geturl()

        #res_params = { 'output' : 'xml' }
        res_params = { 'format' : 'tsv' }
        enc_res_params = urllib.urlencode(res_params)
        #modified_res_url = results_url + '?' + enc_res_params
        modified_res_url = results_url.replace('results','download') + '?' + enc_res_params

        result_request = urllib2.Request(modified_res_url) 
        # url = ( urllib2.urlopen(request).geturl() + '?output=xml') 
        LOGGER.debug('Submitted Pfam search for sequence "{0}...".'
                     .format(seq[:MINSEQLEN]))

        #xml = urllib2.urlopen(result_request).read()
        tsv = urllib2.urlopen(result_request).read()
        # openURL(url, timeout=timeout).read()
        
        # try:
        #     root = ET.XML(xml)
        # except Exception as err:
        #     raise ValueError('failed to parse results XML, check URL: ' + modified_res_url)

        matches = {}
        #for child in root[0]:
            #if child.tag == 'hits':
                # accession = child.get('acc')
                # pfam_id = accession.split('.')[0]
                # matches[pfam_id]={}
                # matches[pfam_id]['accession']=accession
                # matches[pfam_id]['class']='Domain'
                # matches[pfam_id]['id']=child.get('name')
                # matches[pfam_id]['locations']={}
                # matches[pfam_id]['locations']['ali_end']=child[0].get('alisqto')
                # matches[pfam_id]['locations']['ali_start']=child[0].get('alisqfrom')
                # matches[pfam_id]['locations']['bitscore']=child[0].get('bitscore')
                # matches[pfam_id]['locations']['end']=child[0].get('alisqto')
                # matches[pfam_id]['locations']['evalue']=child.get('evalue')
                # matches[pfam_id]['locations']['evidence']='hmmer v3.0'
                # matches[pfam_id]['locations']['hmm_end']=child[0].get('alihmmto')
                # matches[pfam_id]['locations']['hmm_start']=child[0].get('alihmmfrom')
                # matches[pfam_id]['locations']['significant']=child[0].get('significant')    
                # matches[pfam_id]['locations']['start']=child[0].get('alisqfrom')
                # matches[pfam_id]['type']='Pfam-A'
        # return matches

        lines = tsv.split('\n')
        keys = lines[0].split('\t')
        root = {}
        for i, line in enumerate(lines[1:-1]):
            root[i] = {}
            for j, key in enumerate(keys):
                root[i][key] = line.split('\t')[j]

        for child in root.values():
            accession = child['Family Accession']
            pfam_id = accession.split('.')[0]
            matches[pfam_id]={}
            matches[pfam_id]['accession'] = accession
            matches[pfam_id]['class'] = 'Domain'
            matches[pfam_id]['id'] = child['Family id']
            matches[pfam_id]['locations'] = {}
            matches[pfam_id]['locations']['ali_end'] = child['Ali. End']
            matches[pfam_id]['locations']['ali_start'] = child['Ali. Start']
            matches[pfam_id]['locations']['bitscore'] = child['Bit Score']
            matches[pfam_id]['locations']['end'] = child['Env. End']
            matches[pfam_id]['locations']['cond_evalue'] = child['Cond. E-value']
            matches[pfam_id]['locations']['ind_evalue'] = child['Ind. E-value']
            matches[pfam_id]['locations']['evidence'] = 'hmmer v3.0'
            matches[pfam_id]['locations']['hmm_end'] = child['Model End']
            matches[pfam_id]['locations']['hmm_start'] = child['Model Start']
            #matches[pfam_id]['locations']['significant'] = child['significant']   
            matches[pfam_id]['locations']['start'] = child['Env. Start']
            matches[pfam_id]['type'] = 'Pfam-A'
        return matches

    else:
        if len(seq) <= 5:
            idcode = None
            from prody import parsePDBHeader
            try:
                polymers = parsePDBHeader(seq[:4], 'polymers')
            except Exception as err:
                LOGGER.warn('failed to parse header for {0} ({1})'
                            .format(seq[:4], str(err)))
            else:
                chid = seq[4:].upper()
 
            for poly in polymers:
                if chid and poly.chid != chid:
                    continue
                for dbref in poly.dbrefs:
                    if dbref.database != 'UniProt':
                        continue
                    idcode = dbref.idcode
                    LOGGER.info('UniProt ID code {0} for {1} chain '
                                '{2} will be used.'
                                .format(idcode, seq[:4], poly.chid))
                    break
                if idcode is not None:
                    break
            if idcode is None:
                LOGGER.warn('A UniProt ID code for PDB {0} could not be '
                            'parsed.'.format(repr(seq)))
                url = prefix + 'protein/' + seq + '?output=xml'
            else:
                url = prefix + 'protein/' + idcode + '?output=xml'

        else:
            url = prefix + 'protein/' + seq + '?output=xml'

    LOGGER.debug('Retrieving Pfam search results: ' + url)
    xml = None
    while LOGGER.timing('_pfam') < timeout:
        try:
            xml = openURL(url, timeout=timeout).read()
        except Exception:
            pass
        else:
            if xml not in ['PEND','RUN']:
                break

    if not xml:
        raise IOError('Pfam search timed out or failed to parse results '
                      'XML, check URL: ' + url)
    else:
        LOGGER.report('Pfam search completed in %.2fs.', '_pfam')

    if xml.find(b'There was a system error on your last request.') > 0:
        LOGGER.warn('No Pfam matches found for: ' + seq)
        return None

    try:
        root = ET.XML(xml)
    except Exception as err:
        raise ValueError('failed to parse results XML, check URL: ' + url)

    if len(seq) >= MINSEQLEN:
        try:
            xml_matches = root[0][0][0][0]
        except IndexError:
            raise ValueError('failed to parse results XML, check URL: ' + url)
    else:
        key = '{' + prefix + '}'
        results = dictElement(root[0], key)
        try:
            xml_matches = results['matches']
        except KeyError:
            raise ValueError('failed to parse results XML, check URL: ' + url)

    matches = dict()
    for child in xml_matches:

        try:
            accession = child.attrib['accession'][:7]
        except KeyError:
            raise ValueError('failed to parse results XML, check URL: ' + url)

        if not re.search('^P(F|B)[0-9]{5}$', accession):
            raise ValueError('{0} does not match pfam accession'
                             ' format'.format(accession))

        match = matches.setdefault(accession, dict(child.items()))
        locations = match.setdefault('locations', [])
        for loc in child:
            locations.append(dict(loc.items()))

    if len(seq) < MINSEQLEN:
        query = 'Query ' + repr(query)
    else:
        query = 'Query sequence'

    if matches:
        LOGGER.info(query + ' matched {0} Pfam families.'.format(len(matches)))
    else:
        LOGGER.info(query + ' did not match any Pfam families.')
    return matches

def searchUniprotID(query, search_b=False, skip_a=False, **kwargs):
    """Returns Pfam search results in a dictionary.  Matching Pfam accession
    as keys will map to evalue, alignment start and end residue positions.

    :arg query: UniProt ID, PDB identifier, protein sequence, or a 
        sequence file. Sequence queries must not contain gaps and 
        must be at least 16 characters long
    :type query: str

    :arg search_b: search Pfam-B families when **True**
    :type search_b: bool

    :arg skip_a: do not search Pfam-A families when **True**
    :type skip_a: bool

    :arg ga: use gathering threshold when **True**
    :type ga: bool

    :arg evalue: user specified e-value cutoff, must be smaller than 10.0
    :type evalue: float

    :arg timeout: timeout for blocking connection attempt in seconds, default
        is 60
    :type timeout: int

    *query* can also be a PDB identifier, e.g. ``'1mkp'`` or ``'1mkpA'`` with
    chain identifier.  UniProt ID of the specified chain, or the first
    protein chain will be used for searching the Pfam database."""

    query = str(query)
    seq = ''.join(query.split())

    import xml.etree.cElementTree as ET
    LOGGER.timeit('_pfam')
    timeout = int(kwargs.get('timeout', 60))
    url = prefix + 'protein/' + seq + '?output=xml'

    LOGGER.debug('Retrieving Pfam search results: ' + url)
    xml = None
    while LOGGER.timing('_pfam') < timeout:
        try:
            xml = openURL(url, timeout=timeout).read()
        except Exception:
            pass
        else:
            if xml:
                break

    if not xml:
        raise IOError('Pfam search timed out or failed to parse results '
                      'XML, check URL: ' + url)
    else:
        LOGGER.report('Pfam search completed in %.2fs.', '_pfam')

    if xml.find(b'There was a system error on your last request.') > 0:
        LOGGER.warn('No Pfam matches found for: ' + seq)
        return None

    try:
        root = ET.XML(xml)
    except Exception as err:
        raise ValueError('failed to parse results XML, check URL: ' + url)

    result = root[0].get('id')
    return result

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

    url = prefix + 'family/acc?id=' + acc
    handle = openURL(url)
    orig_acc = acc
    acc = handle.readline().strip()
    if PY3K:
        acc = acc.decode()
    url_flag = False

    if not re.search('(?<=PF)[0-9]{5}$', acc):
        raise ValueError('{0} is not a valid Pfam ID or Accession Code'
                         .format(repr(orig_acc)))

    if alignment not in DOWNLOAD_FORMATS:
        raise ValueError('alignment must be one of full, seed, ncbi or'
                         ' metagenomics')
    if alignment == 'ncbi' or alignment == 'metagenomics' or alignment == 'uniprot':
        url = (prefix + 'family/' + acc + '/alignment/' +
               alignment + '/gzipped')
        url_flag = True
        extension = '.sth'
    else:
        if not kwargs:
            url = (prefix + 'family/' + acc + '/alignment/' +
                   alignment + '/gzipped')
            url_flag = True
            extension = '.sth'
        else:
            align_format = kwargs.get('format', 'selex').lower()

            if align_format not in FORMAT_OPTIONS['format']:
                raise ValueError('alignment format must be of type selex'
                                 ' stockholm or fasta. MSF not supported')

            if align_format == SELEX:
                align_format, extension = 'pfam', '.slx'
            elif align_format == FASTA:
                extension = '.fasta'
            else:
                extension = '.sth'

            gaps = str(kwargs.get('gaps', 'dashes')).lower()
            if gaps not in FORMAT_OPTIONS['gaps']:
                raise ValueError('gaps must be of type mixed, dots, dashes, '
                                 'or None')

            inserts = kwargs.get('inserts', 'upper').lower()
            if(inserts not in FORMAT_OPTIONS['inserts']):
                raise ValueError('inserts must be of type lower or upper')

            order = kwargs.get('order', 'tree').lower()
            if order not in FORMAT_OPTIONS['order']:
                raise ValueError('order must be of type tree or alphabetical')

            url = (prefix + 'family/' + acc + '/alignment/'
                   + alignment + '/format?format=' + align_format +
                   '&alnType=' + alignment + '&order=' + order[0] +
                   '&case=' + inserts[0] + '&gaps=' + gaps + '&download=1')

    response = openURL(url, timeout=int(kwargs.get('timeout', 60)))
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
        f_out.write(response.read())
        f_out.close()
    else:
        if url_flag:
            gunzip(response.read(), filepath)
        else:
            with open(filepath, 'wb') as f_out:
                f_out.write(response.read())

    filepath = relpath(filepath)
    LOGGER.info('Pfam MSA for {0} is written as {1}.'
                .format(orig_acc, filepath))

    return filepath

def parsePfamPDBs(query, data=[], **kwargs):
    """Returns a list of AtomGroups containing sections of chains that 
    correspond to a particular PFAM domain family. These are defined by 
    alignment start and end residue numbers.

    :arg query: UniProt ID or PDB ID
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
        pfam_matches = searchPfam(query)
        keys = list(pfam_matches.keys())

        if isinstance(start, Integral):
            start_diff = []
            for i, key in enumerate(pfam_matches):
                start_diff.append(int(pfam_matches[key]['locations'][0]['start']) - start)
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
    results = parsePDB(*pdb_ids, chain=chains, header=True, **kwargs)

    ags, headers = results
    ags, headers = list(ags), list(headers)

    if model == 0:
        LOGGER.info('only header is requested and returned')
        return results

    if header:
        results = (ags, headers)
    else:
#        ags = results
#        ags = list(ags)
        results = ags

    LOGGER.progress('Extracting Pfam domains...', len(ags))
    comma_splitter = re.compile(r'\s*,\s*').split
    no_info = []
    for i, ag in enumerate(ags):
        LOGGER.update(i)
        data_dict = data_dicts[i]
        pfamRange = data_dict['UniprotResnumRange'].split('-')
        uniprotAcc = data_dict['UniprotAcc']
        uniData = queryUniprot(uniprotAcc)
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

        multi = False
        if found:
            header = headers[i]
            chain_accessions = [dbref.accession 
                                for dbref in header[data_dict['chain']].dbrefs]
            if len(chain_accessions) > 1:
                multi = True
            else:
                multi = False
            right_part = np.where(np.array(chain_accessions) == 
                                  data_dict['UniprotAcc'])[0][0]
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

