# -*- coding: utf-8 -*-
"""This module defines functions for blast searching the Protein Data Bank."""

import os.path
import numpy as np
from prody import LOGGER, PY3K
from prody.utilities import dictElement, openURL, which
from prody.sequence import parseMSA, MSA, Sequence

import platform, os, re, sys, time, urllib

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

from prody.atomic import Atomic
from prody.proteins.pdbfile import parsePDB

__all__ = ['PDBBlastRecord', 'blastPDB',]

def blastPDB(sequence, filename=None, **kwargs):
    """Returns a :class:`PDBBlastRecord` instance that contains results from
    blast searching *sequence* against the PDB using NCBI blastp.

    :arg sequence: an object with an associated sequence string 
         or a sequence string itself
    :type sequence: :class:`Atomic`, :class:`Sequence`, or str

    :arg filename: a *filename* to save the results in XML format
    :type filename: str

    *hitlist_size* (default is ``250``) and *expect* (default is ``1e-10``)
    search parameters can be adjusted by the user.  *sleep* keyword argument
    (default is ``2`` seconds) determines how long to wait to reconnect for
    results.  Sleep time is multiplied by 1.5 when results are not ready.  
    *timeout* (default is 120 s) determines when to give up waiting for the results.
    """

    if sequence == 'runexample':
        sequence = ('ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPI'
                    'SDHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMN'
                    'DAYDIVKMKKSNISPNFNFMGQLLDFERTL')

    elif isinstance(sequence, Atomic):
        sequence = sequence.calpha.getSequence()

    elif isinstance(sequence, Sequence):
        sequence = str(sequence)

    elif isinstance(sequence, str):
        if len(sequence) == 4 or len(sequence) == 5:
            sequence = parsePDB(sequence)
            sequence = ag.calpha.getSequence()
        sequence = ''.join(sequence.split())

    else:
        raise TypeError('sequence must be Atomic, Sequence, or str not {0}'
                        .format(type(sequence)))

    headers = {'User-agent': 'ProDy'}
    query = [('DATABASE', 'pdb'), ('ENTREZ_QUERY', '(none)'),
             ('PROGRAM', 'blastp'),]

    expect = float(kwargs.pop('expect', 10e-10))
    if expect <= 0:
        raise ValueError('expect must be a positive number')
    query.append(('EXPECT', expect))
    hitlist_size = int(kwargs.pop('hitlist_size', 250))
    if hitlist_size <= 0:
        raise ValueError('expect must be a positive integer')
    query.append(('HITLIST_SIZE', hitlist_size))
    query.append(('QUERY', sequence))
    query.append(('CMD', 'Put'))

    sleep = float(kwargs.pop('sleep', 2))
    timeout = float(kwargs.pop('timeout', 120))

    try:
        import urllib.parse
        urlencode = lambda data: bytes(urllib.parse.urlencode(data), 'utf-8')
    except ImportError:
        from urllib import urlencode

    url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

    data = urlencode(query)
    LOGGER.timeit('_prody_blast')
    LOGGER.info('Blast searching NCBI PDB database for "{0}..."'
                .format(sequence[:5]))
    handle = openURL(url, data=data, headers=headers)

    html = handle.read()
    index = html.find(b'RID =')
    if index == -1:
        raise Exception('NCBI did not return expected response.')
    else:
        last = html.find(b'\n', index)
        rid = html[index + len('RID ='):last].strip()

    index = html.find(b'RTOE =')
    if index == -1:
        rtoe = None # This is not used
    else:
        last = html.find(b'\n', index)
        rtoe = int(html[index + len('RTOE ='):last].strip())

    query = [('ALIGNMENTS', 500), ('DESCRIPTIONS', 500),
             ('FORMAT_TYPE', 'XML'), ('RID', rid), ('CMD', 'Get')]
    data = urlencode(query)

    while True:
        LOGGER.sleep(int(sleep), 'to reconnect NCBI for search results.')
        LOGGER.write('Connecting to NCBI for search results...')
        handle = openURL(url, data=data, headers=headers)
        results = handle.read()
        index = results.find(b'Status=')
        LOGGER.clear()
        if index < 0:
            break
        last = results.index(b'\n', index)
        status = results[index+len('Status='):last].strip()
        if status.upper() == 'READY':
            break
        sleep = int(sleep * 1.5)
        if LOGGER.timing('_prody_blast') > timeout:
            LOGGER.warn('Blast search time out.')
            return None
    LOGGER.clear()
    LOGGER.report('Blast search completed in %.1fs.', '_prody_blast')

    try:
        ext_xml = filename.lower().endswith('.xml')
    except AttributeError:
        pass
    else:
        if not ext_xml:
            filename += '.xml'
        out = open(filename, 'w')
        out.write(results)
        out.close()
        LOGGER.info('Results are saved as {0}.'.format(repr(filename)))

    return PDBBlastRecord(results, sequence)


class PDBBlastRecord(object):

    """A class to store results from blast searches."""


    def __init__(self, xml, sequence=None):
        """Instantiate a PDBBlastRecord object instance.

        :arg xml: blast search results in XML format or an XML file that
            contains the results
        :type xml: str

        :arg sequence: query sequence
        :type sequence: str
        """

        if sequence:
            try:
                sequence = ''.join(sequence.split())
                _ = sequence.isalpha()
            except AttributeError:
                raise TypeError('sequence must be a string')
            else:
                if not _:
                    raise ValueError('not a valid protein sequence')
        self._sequence = sequence

        import xml.etree.cElementTree as ET
        if len(xml) < 100:
            if os.path.isfile(xml):
                xml = ET.parse(xml)
                root = xml.getroot()
            else:
                raise ValueError('xml is not a filename and does not look like'
                                 ' a valid XML string')
        else:
            root = ET.XML(xml)

        root = dictElement(root, 'BlastOutput_')
        if root['db'] != 'pdb':
            raise ValueError('blast search database in xml must be "pdb"')
        if root['program'] != 'blastp':
            raise ValueError('blast search program in xml must be "blastp"')
        self._param = dictElement(root['param'][0], 'Parameters_')

        query_len = int(root['query-len'])
        if sequence and len(sequence) != query_len:
            raise ValueError('query-len and the length of the sequence do not '
                             'match, xml data may not be for given sequence')
        hits = []
        for iteration in root['iterations']:
            for hit in dictElement(iteration, 'Iteration_')['hits']:
                hit = dictElement(hit, 'Hit_')
                data = dictElement(hit['hsps'][0], 'Hsp_')
                for key in ['align-len', 'gaps', 'hit-frame', 'hit-from',
                            'hit-to', 'identity', 'positive', 'query-frame',
                            'query-from', 'query-to']:
                    data[key] = int(data[key])
                data['query-len'] = query_len
                for key in ['evalue', 'bit-score', 'score']:
                    data[key] = float(data[key])
                p_identity = 100.0 * data['identity'] / (data['query-to'] -
                                                    data['query-from'] + 1)
                data['percent_identity'] = p_identity
                p_overlap = (100.0 * (data['align-len'] - data['gaps']) /
                              query_len)
                data['percent_coverage'] = p_overlap
                
                for item in (hit['id'] + hit['def']).split('>gi'):
                    head, title = item.split(None, 1)
                    head = head.split('|')
                    pdb_id = head[-2].lower()
                    chain_id = head[-1][:1]
                    pdbch = dict(data)
                    pdbch['pdb_id'] = pdb_id
                    pdbch['chain_id'] = chain_id
                    pdbch['title'] = (head[-1][1:] + title).strip()
                    hits.append((p_identity, p_overlap, pdbch))
        hits.sort(key=lambda hit: hit[0], reverse=True)
        self._hits = hits

    def getSequence(self):
        """Returns the query sequence that was used in the search."""

        return self._sequence

    def getParameters(self):
        """Returns parameters used in blast search."""

        return dict(self._param)

    def getHits(self, percent_identity=0., percent_overlap=0., chain=False):
        """Returns a dictionary in which PDB identifiers are mapped to structure
        and alignment information.

        :arg percent_identity: PDB hits with percent sequence identity equal
            to or higher than this value will be returned, default is ``0.``
        :type percent_identity: float
        :arg percent_overlap: PDB hits with percent coverage of the query
          sequence equivalent or better will be returned, default is ``0.``
        :type percent_overlap: float
        :arg chain: if chain is **True**, individual chains in a PDB file
          will be considered as separate hits , default is **False**
        :type chain: bool"""

        assert isinstance(percent_identity, (float, int)), \
            'percent_identity must be a float or an integer'
        assert isinstance(percent_overlap, (float, int)), \
            'percent_overlap must be a float or an integer'
        assert isinstance(chain, bool), 'chain must be a boolean'

        hits = {}
        for p_identity, p_overlap, hit in self._hits:
            if p_identity < percent_identity:
                break
            if p_overlap < percent_overlap:
                continue
            if chain:
                key = (hit['pdb_id'], hit['chain_id'])
            else:
                key = hit['pdb_id']
            if not key in hits:
                hits[key] = hit
        return hits

    def getBest(self):
        """Returns a dictionary containing structure and alignment information
        for the hit with highest sequence identity."""

        return dict(self._hits[0][2])

    def writeSequences(self, filename, **kwargs):
        """
        Returns a plot that contains a dendrogram of the sequence similarities among
        the sequences in given hit list. 

        :arg hits: A dictionary that contains hits that are obtained from a blast record object. 
        :type hits: dict

        Arguments of getHits can be parsed as kwargs.
        """
        if not filename.lower().endswith('.fasta'):
            filename += '.fasta'
    
        with open(filename,'w') as f_out:
            for z in self.getHits(**kwargs):
                f_out.write(">" + str(z) + "\n")
                f_out.write(hits[z]['hseq'])
                f_out.write("\n")

        return filename
