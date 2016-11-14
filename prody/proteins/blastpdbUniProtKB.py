# -*- coding: utf-8 -*-
"""This module defines functions for blast searching Protein Data Bank."""

import os.path

from prody import LOGGER
from prody.utilities import dictElement, openURL

__all__ = ['PDBBlastRecord', 'blastPDB']

def blastPDBUniProtKB(sequence, filename=None, **kwargs):
    """Returns a :class:`PDBBlastRecord` instance that contains results from
    blast searching of ProteinDataBank database *sequence* using NCBI blastp.

    :arg sequence: single-letter code amino acid sequence of the protein
        without any gap characters, all white spaces will be removed
    :type sequence: str

    :arg filename: a *filename* to save the results in XML format
    :type filename: str

    *hitlist_size* (default is ``250``) and *expect* (default is ``1e-10``)
    search parameters can be adjusted by the user.  *sleep* keyword argument
    (default is ``2`` seconds) determines how long to wait to reconnect for
    results.  Sleep time is doubled when results are not ready.  *timeout*
    (default is 120s) determines when to give up waiting for the results.
    """

    if sequence == 'runexample':
        sequence = ('ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPI'
                    'SDHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMN'
                    'DAYDIVKMKKSNISPNFNFMGQLLDFERTL')
    else:
        try:
            sequence = ''.join(sequence.split())
            _ = sequence.isalpha()
        except AttributeError:
            raise TypeError('sequence must be a string')
        else:
            if not _:
                raise ValueError('not a valid protein sequence')
    headers = {'User-agent': 'ProDy'}

    query = [('DATABASE', 'swissprot'), ('ENTREZ_QUERY', '(none)'),
             ('PROGRAM', 'blastp'),]
    expect = float(kwargs.pop('expect', 10e-10))
    if expect <= 0:
        raise ValueError('expect must be a positive number')
    query.append(('EXPECT', expect))
    hitlist_size = int(kwargs.pop('hitlist_size', 250))
    if hitlist_size <= 0:
        raise ValueError('expect must be a positive integer')
    psiblast = 'true'
    query.append(('RUN_PSIBLAST', psiblast))
    query.append(('HITLIST_SIZE', hitlist_size))
    query.append(('QUERY', sequence))
    query.append(('CMD', 'Put'))

    sleep = float(kwargs.pop('sleep', 2))
    timeout = float(kwargs.pop('timeout', 120))

    if kwargs:
        LOGGER.warn('Keyword argument(s) {0} are not used.'
                    .format(', '.join([repr(key) for key in kwargs])))

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
    index = html.find(b'name="RID" type="hidden" value="')
    if index == -1:
        raise Exception('NCBI did not return expected response.')
    else:
        last = html.find(b'>',index)
        rid = html[index + len('name="RID" type="hidden" value="'):last-1].strip()

    index = html.find(b'name="RTOE" type="hidden" value="')
    if index == -1:
        rtoe = None # This is not used
    else:
        last = html.find(b'>', index)
        rtoe = html[index + len('name="RTOE" type="hidden" value="'):last-1].strip()

    query = [('ALIGNMENTS', 500), ('DESCRIPTIONS', 500),
             ('FORMAT_TYPE', 'XML'), ('RID', rid), ('CMD', 'Get')]
    data = urlencode(query)

    while True:
        LOGGER.sleep(int(sleep), 'to reconnect NCBI for search results.')
        LOGGER.write('Connecting NCBI for search results...')
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
    return SwissProtBlastRecord(results, sequence)

class SwissProtBlastRecord(object):

    """A class to store results from ProteinDataBank blast search."""


    def __init__(self, xml, sequence=None):
        """Instantiate a PDBlast object instance.

        :arg xml: blast search results in XML format or an XML file that
            contains the results
        :type xml: str

        :arg sequence: query sequence
        :type sequence: str"""

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
        if root['db'] != 'swissprot':
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
                data['percent_overlap'] = p_overlap
                for item in (hit['id'] + ' ' + hit['def']).split('>gi'):
                    #>gi|1633462|pdb|4AKE|A Chain A, Adenylate Kinase
                    #                        __________TITLE__________
                    head, title = item.split(None, 1)
                    head = head.split('|')
                    seqInfo = dict(data)
                    accession = head[-2].split(".")[0]
                    protName = head[-1].split("_")[0]
                    species = head[-1].split("_")[1]
                    seqInfo['accession'] = accession
                    seqInfo['protName'] = protName
                    seqInfo['species'] = species
                    hits.append((p_identity, p_overlap, seqInfo))
        hits.sort(key=lambda hit: hit[0], reverse=True)
        self._hits = hits

    def getSequence(self):
        """Returns the query sequence that was used in the search."""

        return self._sequence

    def getParameters(self):
        """Returns parameters used in blast search."""

        return dict(self._param)

    def getHits(self, percent_identity=90., percent_overlap=70.):
        """Returns a dictionary in which PDB identifiers are mapped to structure
        and alignment information.

        :arg percent_identity: PDB hits with percent sequence identity equal
            to or higher than this value will be returned, default is ``90.0``
        :type percent_identity: float
        :arg percent_overlap: PDB hits with percent coverage of the query
          sequence equivalent or better will be returned, default is ``70.0``
        :type percent_overlap: float
        :arg chain: if chain is **True**, individual chains in a PDB file
          will be considered as separate hits , default is **False**
        :type chain: bool"""

        assert isinstance(percent_identity, (float, int)), \
            'percent_identity must be a float or an integer'
        assert isinstance(percent_overlap, (float, int)), \
            'percent_overlap must be a float or an integer'

        hits = {}
        for p_identity, p_overlap, hit in self._hits:
            if p_identity < percent_identity:
                break
            if p_overlap < percent_overlap:
                continue
            if not key in hits:
                hits[key] = hit
        return hits

    def getBest(self):
        """Returns a dictionary containing structure and alignment information
        for the hit with highest sequence identity."""

        return dict(self._hits[0][2])



