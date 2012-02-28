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

"""This module defines functions for accessing Protein Data Bank."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path
import shutil
from glob import glob

import numpy as np

from prody.tools import makePath, gunzip, relpath, openFile
from prody import getPackagePath

__all__ = ['PDBBlastRecord', 'blastPDB', 'fetchPDB', 
           'getPDBLocalFolder', 'getPDBMirrorPath', 'getWWPDBFTPServer', 
           'setPDBLocalFolder', 'setPDBMirrorPath', 'setWWPDBFTPServer',
           'fetchPDBClusters', 'loadPDBClusters', 'getPDBCluster',]
           
pkg = __import__(__package__)
LOGGER = pkg.LOGGER
SETTINGS = pkg.SETTINGS
           
PDB_CLUSTERS = {30: None, 40: None, 50: None, 70: None, 
                90: None, 95: None, 100: None}
PDB_CLUSTERS_UPDATE_WARNING = True

_PDB_EXTENSIONS = set(['.pdb', '.PDB', '.gz', '.GZ', '.ent', '.ENT', 
                       '.pdb.gz', '.PDB.GZ', '.ent.gz', '.ENT.GZ',
                       '.xml', '.XML', '.xml.gz', '.XML.GZ',
                       '.cif', '.CIF', '.cif.gz', '.CIF.GZ',])
_PDB_FORMATS = set(['pdb', 'cif', 'xml'])

_WWPDB_RCSB = ('RCSB PDB (USA)', 'ftp.wwpdb.org', '/pub/pdb/')
_WWPDB_PDBe = ('PDBe (Europe)', 'ftp.ebi.ac.uk', '/pub/databases/rcsb/pdb/')
_WWPDB_PDBj = ('PDBj (Japan)', 'pdb.protein.osaka-u.ac.jp', '/pub/pdb/')

WWPDB_FTP_SERVERS = {
    'rcsb'   : _WWPDB_RCSB,
    'usa'    : _WWPDB_RCSB,
    'us'     : _WWPDB_RCSB,
    'pdbe'   : _WWPDB_PDBe,
    'euro'   : _WWPDB_PDBe,
    'europe' : _WWPDB_PDBe,
    'eu'     : _WWPDB_PDBe,
    'pdbj'   : _WWPDB_PDBj,
    'japan'  : _WWPDB_PDBj,
    'jp'     : _WWPDB_PDBj,
}




PROTSEQ_ALPHABET = set('ARNDCQEGHILKMFPSTWYVBJOUXZ-')

def checkSequence(sequence):
    """Check validity of a protein sequence.  If a valid sequence, return
    after standardizing it (make all upper case, remove spaces, etc.), 
    otherwise return ``False``."""
    
    if isinstance(sequence, str):
        sequence = ''.join(sequence.split()).upper()
        if PROTSEQ_ALPHABET.issuperset(set(sequence)):
            return sequence
    return False

def blastPDB(sequence, filename=None, **kwargs):
    """Blast search *sequence* in ProteinDataBank database using NCBI blastp.
        
    :arg sequence: Single-letter code amino acid sequence of the protein.
    :type sequence: str 
    
    :arg filename: Provide a *filename* to save the results in XML format. 
    :type filename: str
    
    Results are returned in a :class:`PDBBlastRecord` instance.
    
    User can adjust *hitlist_size* (default is 250) and *expect* (default is 
    1e-10) parameter values.
    
    User can also set the time interval for checking the results using
    *sleep* keyword argument.  Default value for *sleep* is 2 seconds.
    Finally, user can set a *timeout* value after which the function 
    will return ``None`` if the results are not ready.  Default for
    *timeout* is 30 seconds."""
    
    if kwargs.pop('runexample', False):
        sequence = 'ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPI'\
                   'SDHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMN'\
                   'DAYDIVKMKKSNISPNFNFMGQLLDFERTL'
    else:
        sequence = checkSequence(sequence)
        
    if not sequence:
        raise ValueError('not a valid protein sequence')

    query = [('DATABASE', 'pdb'), ('ENTREZ_QUERY', '(none)'),
             ('PROGRAM', 'blastp'),] 
    expect = kwargs.pop('expect', 10e-10)
    assert isinstance(expect, (float, int)), 'expect must be a float'
    assert expect > 0, 'expect must be a positive number'
    query.append(('EXPECT', expect))
    hitlist_size = kwargs.pop('hitlist_size', 250)
    assert isinstance(hitlist_size, int), 'hitlist_size must be an integer'
    assert hitlist_size > 0, 'expect must be a positive integer'
    query.append(('HITLIST_SIZE', hitlist_size))
    query.append(('QUERY', sequence))
    query.append(('CMD', 'Put'))
    
    sleep = float(kwargs.pop('sleep', 2))
    timeout = float(kwargs.pop('timeout', 20))
    
    if kwargs:
        LOGGER.warning("Keyword argument(s) '{0:s}' are not used."
                       .format("', '".join(kwargs.keys())))

    import urllib, urllib2
    
    url = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi'
    
    data = urllib.urlencode(query)
    LOGGER.timeit()
    LOGGER.info('Blast searching NCBI PDB database for "{0:s}..."'
                .format(sequence[:5]))
    request = urllib2.Request(url, data, {'User-agent': 'ProDy'})
    handle = urllib2.urlopen(request)
    
    html = handle.read()
    index = html.find('RID =')
    if index == -1:
        raise Exception('NCBI did not return expected response.')
    else:
        last = html.find('\n', index)
        rid = html[index + len('RID ='):last].strip()

    index = html.find('RTOE =')
    if index == -1:
        rtoe = None # This is not used
    else:
        last = html.find('\n', index)
        rtoe = int(html[index + len('RTOE ='):last].strip())

    query = [('ALIGNMENTS', 500), ('DESCRIPTIONS', 500), 
             ('FORMAT_TYPE', 'XML'), ('RID', rid), ('CMD', 'Get')]
    data = urllib.urlencode(query)
    
    while True:
        LOGGER.sleep(int(sleep), ' to connect NCBI for search results.')
        LOGGER.write('Connecting NCBI for search results...')
        request = urllib2.Request(url, data, {'User-agent': 'ProDy'})
        handle = urllib2.urlopen(request)
        results = handle.read()
        index = results.find('Status=')
        LOGGER.clear()
        if index < 0:
            break
        last = results.index('\n', index)
        status = results[index+len('Status='):last].strip()
        if status.upper() == 'READY':
            break
        sleep *= 2
        if LOGGER.timing() > timeout:
            LOGGER.warning('Blast search time out.')
            return None
    LOGGER.clear()
    LOGGER.timing('Blast search completed in %.1fs.')
    if filename is not None:
        filename = str(filename)
        if not filename.lower().endswith('.xml'):
                filename += '.xml'        
        out = open(filename, 'w')
        out.write(results)
        out.close()
        LOGGER.info('Results are saved as {0:s}.'.format(filename))
    return PDBBlastRecord(sequence, results)

def children2dict(element, prefix=None):
    dict_ = {}
    length = False
    if isinstance(prefix, str):
        length = len(prefix)
    for child in element:
        tag = child.tag
        if length and tag.startswith(prefix):
            tag = tag[length:]
        if len(child) == 0:
            dict_[tag] = child.text
        else:
            dict_[tag] = child
    return dict_

class PDBBlastRecord(object):

    """A class to store results from ProteinDataBank blast search."""
    
    __slots__ = ['_param', '_sequence', '_hits']

    def __init__(self, sequence, xml):
        """Instantiate a PDBlast object instance.
        
        :arg sequence: query sequence
        :type xml: str
        :arg xml: blast search results in XML format or an XML file that 
            contains the results
        :type xml: str
        
        """
        
        sequence = checkSequence(sequence)
        if not sequence:
            raise ValueError('not a valid protein sequence')
        self._sequence = sequence
        
        import xml.etree.cElementTree as ET
        assert isinstance(xml, str), 'xml must be a string'
        if len(xml) < 100:
            if os.path.isfile(xml):
                xml = ET.parse(xml)
                xml.getroot()
            else:
                raise ValueError('xml is not a filename and does not look like'
                                 ' a valid XML string')
        else:
            root = ET.XML(xml)
        
        root = children2dict(root, 'BlastOutput_')
        if root['db'] != 'pdb':
            raise ValueError('blast search database in xml must be "pdb"')
        if root['program'] != 'blastp':
            raise ValueError('blast search program in xml must be "blastp"')
        self._param = children2dict(root['param'][0], 'Parameters_')
        query_length = int(root['query-len'])
        if len(sequence) != query_length:
            raise ValueError('query-len and the length of the sequence do not '
                             'match, xml data may not be for the given' 
                             'sequence')
        hits = [] 
        for iteration in root['iterations']:
            for hit in children2dict(iteration, 'Iteration_')['hits']:
                hit = children2dict(hit, 'Hit_')
                data = children2dict(hit['hsps'][0], 'Hsp_')
                for key in ['align-len', 'gaps', 'hit-frame', 'hit-from',
                            'hit-to', 'identity', 'positive', 'query-frame',
                            'query-from', 'query-to']:
                    data[key] = int(data[key])
                for key in ['evalue', 'bit-score', 'score']:
                    data[key] = float(data[key])
                p_identity = 100.0 * data['identity'] / data['align-len']
                data['percent_identity'] = p_identity
                p_coverage = 100.0 * (data['align-len'] - data['gaps']) / \
                    query_length
                data['percent_coverage'] = p_coverage  
                for item in (hit['id'] + hit['def']).split('>gi'):
                    #>gi|1633462|pdb|4AKE|A Chain A, Adenylate Kinase
                    #                        __________TITLE__________
                    head, title = item.split(None, 1)
                    head = head.split('|')
                    pdb_id = head[-2].lower() 
                    chain_id = head[-1][0]
                    pdbch = dict(data)
                    pdbch['pdb_id'] = pdb_id
                    pdbch['chain_id'] = chain_id
                    pdbch['title'] = (head[-1][1:] + title).strip()
                    hits.append((p_identity, p_coverage, pdbch))
        hits.sort(reverse=True)
        self._hits = hits
        
        
    def getSequence(self):    
        """Return the query sequence that was used in the search."""
        
        return self._sequence
        
    def getParameters(self):
        """Return parameters used in blast search."""
        
        return self._param
        
    def getHits(self, percent_identity=90., percent_coverage=70., chain=False):
        """Return a dictionary that contains hits.
        
        Returns a dictionary whose keys are PDB identifiers.  Each dictionary
        entry is also a dictionary that contains information on the blast hit
        and alignment. 
        
        :arg percent_identity: PDB hits with percent sequence identity equal 
            to or higher than this value will be returned, default is ``90.0``.
        :type percent_identity: float
        :arg percent_coverage: PDB hits with percent coverage of the query 
          sequence equivalent or better will be returned, default is ``70.0``.
        :type percent_coverage: float
        :arg chain: if chain is set ``True``, individual chains in a PDB file
          will be considered as separate hits when they match the query
          sequence, default is ``False``.
        :type chain: bool
        
        """
        
        assert isinstance(percent_identity, (float, int)), \
            'percent_identity must be a float or an integer'
        assert isinstance(percent_coverage, (float, int)), \
            'percent_coverage must be a float or an integer'
        assert isinstance(chain, bool), 'chain must be a boolean'
        
        hits = {}
        for p_identity, p_coverage, hit in self._hits:
            if p_identity < percent_identity:
                break
            if p_coverage < percent_coverage:
                continue
            if chain:
                key = (hit['pdb_id'], hit['chain_id'])
            else:
                key = hit['pdb_id']
            if not key in hits: 
                hits[key] = hit
        return hits
    
    def getBest(self):
        """Returns a dictionary for the hit with highest sequence identity."""
        
        return dict(self._hits[0][2])
    



