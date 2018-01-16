# -*- coding: utf-8 -*-
"""This module defines functions for interfacing Pfam database."""

__author__ = 'Anindita Dutta, Ahmet Bakan, Cihan Kaya'

import re

from os.path import join, isfile
from prody import LOGGER, PY3K
from prody.utilities import makePath, openURL, gunzip, openFile, dictElement
from prody.utilities import relpath

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

__all__ = ['searchPfam', 'fetchPfamMSA', 'searchUniprotID']

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


def searchPfam(query, **kwargs):
    """Returns Pfam search results in a dictionary.  Matching Pfam accession
    as keys will map to evalue, alignment start and end residue positions.

    :arg query: UniProt ID, PDB identifier, protein sequence, or a sequence
        file, sequence queries must not contain without gaps and must be at
        least 16 characters long
    :type query: str

    :arg timeout: timeout for blocking connection attempt in seconds, default
        is 60
    :type timeout: int

    *query* can also be a PDB identifier, e.g. ``'1mkp'`` or ``'1mkpA'`` with
    chain identifier.  UniProt ID of the specified chain, or the first
    protein chain will be used for searching the Pfam database."""

    prefix = '{https://pfam.xfam.org/}'
    query = str(query)
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

        results_url = urllib2.urlopen(request).getheader('location')

        res_params = { 'output' : 'xml' }
        enc_res_params = urllib.urlencode(res_params)
        modified_res_url = results_url + '?' + enc_res_params

        result_request = urllib2.Request(modified_res_url) 
        # url = ( urllib2.urlopen(request).geturl() + '?output=xml') 
        LOGGER.debug('Submitted Pfam search for sequence "{0}...".'
                     .format(seq[:MINSEQLEN]))

        xml = urllib2.urlopen(result_request).read()
        # openURL(url, timeout=timeout).read()
        
        try:
            root = ET.XML(xml)
        except Exception as err:
            raise ValueError('failed to parse results XML, check URL: ' + url)
        matches = {}
        for child in root[0]:
            if child.tag == 'hits':
                accession = child.get('acc')
                pfam_id = accession.split('.')[0]
                matches[pfam_id]={}
                matches[pfam_id]['accession']=accession
                matches[pfam_id]['class']='Domain'
                matches[pfam_id]['id']=child.get('name')
                matches[pfam_id]['locations']={}
                matches[pfam_id]['locations']['ali_end']=child[0].get('alisqto')
                matches[pfam_id]['locations']['ali_start']=child[0].get('alisqfrom')
                matches[pfam_id]['locations']['bitscore']=child[0].get('bitscore')
                matches[pfam_id]['locations']['end']=child[0].get('alisqto')
                matches[pfam_id]['locations']['evalue']=child.get('evalue')
                matches[pfam_id]['locations']['evidence']='hmmer v3.0'
                matches[pfam_id]['locations']['hmm_end']=child[0].get('alihmmto')
                matches[pfam_id]['locations']['hmm_start']=child[0].get('alihmmfrom')
                matches[pfam_id]['locations']['significant']=child[0].get('significant')    
                matches[pfam_id]['locations']['start']=child[0].get('alisqfrom')
                matches[pfam_id]['type']='Pfam-A'
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
                url = 'https://pfam.xfam.org/protein/' + seq + '?output=xml'
            else:
                url = ('https://pfam.xfam.org/protein/' +
                       idcode + '?output=xml')

        else:
            url = 'https://pfam.xfam.org/protein/' + seq + '?output=xml'

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
	key = '{' + root.items()[1][1].split()[0] + '}'
	print key
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

    :arg query: UniProt ID, PDB identifier, protein sequence, or a sequence
        file, sequence queries must not contain without gaps and must be at
        least 16 characters long
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

    prefix = '{http://pfam.xfam.org/}'
    query = str(query)
    seq = ''.join(query.split())

    import xml.etree.cElementTree as ET
    LOGGER.timeit('_pfam')
    timeout = int(kwargs.get('timeout', 60))
    url = 'http://pfam.xfam.org/protein/' + seq + '?output=xml'

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

    url = 'http://pfam.xfam.org/family/acc?id=' + acc
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
        url = ('http://pfam.xfam.org/family/' + acc + '/alignment/' +
               alignment + '/gzipped')
        url_flag = True
        extension = '.sth'
    else:
        if not kwargs:
            url = ('http://pfam.xfam.org/family/' + acc + '/alignment/' +
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

            url = ('http://pfam.xfam.org/family/' + acc + '/alignment/'
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


if __name__ == '__main__':

    from prody import *
    #filepath = fetchPfamMSA('PF00497',alignment='seed',compressed=False,
    #                         format='stockholm',gaps='dashes')
    #print filepath
    #results = list(iterSequences(filepath))

    #filepath1 = fetchPfamMSA('PF00007',alignment='seed',compressed=True,
    #                         timeout=5)
    #filepath2 = fetchPfamMSA('PF00007',alignment='seed',compressed=True,
    #                         format='selex')
    #filepath3 = fetchPfamMSA('PF00007',alignment='seed',compressed=False,
    #                         format='fasta', outname='mymsa')
    #results_sth = list(MSAFile(filepath1))
    #results_slx = list(MSAFile(filepath2))
    #results_fasta = list(MSAFile(filepath3))
    #from Bio import AlignIO
    #alignment = AlignIO.read(filepath3,'fasta')
    #results_obj = list(MSAFile(alignment))
    #import numpy
    #numpy.testing.assert_equal(results_fasta,results_obj)

    #matches1 = searchPfam('P12821')
    #matches1 = searchPfam('P08581')
    #matches2 = searchPfam('test.seq')

    """matches2 = searchPfam('PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGSSEPCALCSLHSIGKIGGAQNRSYSKLLC\
GLLAERLRISPDRVYINYYDMNAANVGWNNSTFA', evalue=2, skipAs=True)
    """
    #matches3 = searchPfam('NSIQIGGLFPRGADQEYSAFRVGMVQFSTSEFRLTPHIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFC\
#GTLHVSFITPSFPTDGTHPFVIQMRPDLKGALLSLIEYYQWDKFAYLYDSDRGLSTLQAVLDSAAEKKWQVTAINVGNINNDKKDETYRSLFQDLELKKERRVILDCERDKVNDIVDQVITIGKHVKGYHYIIANLGFTDGDLLKIQFGGAEVSGFQIVD\
#YDDSLVSKFIERWSTLEEKEYPGAHTATIKYTSALTYDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEI\
#ERALKQVQVEGLSGNIKFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVLTEDDTSGLEQKTVVVTTILESPYVMMK\
#ANHAALAGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKIWNGMVGELVYGKADIAIAPLTITLVREEVI\
#DFSKPFMSLGISIMIKKPQKSKPGVFSFLDPLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTN\
#EFGIFNSLWFSLGAFMQQGADISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDLSKQTEIAYGT\
#LDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGKYAYLLESTMNEYIEQRKPCDTMKVGGNLDS\
#KGYGIATPKGSSLGTPVNLAVLKLSEQGLLDKLKNKWWYDKGECGAKDSGSKEKTSALSLSNVAGVFYILVGGLGLAMLV\
#ALIEFCYKSRAEAKRMKGLVPRG', delay=10, evalue=2, searchBs=True)
