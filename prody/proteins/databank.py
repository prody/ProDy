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

from prody.atomic import ATOMIC_DATA_FIELDS
from prody.atomic import AtomGroup
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

def getPDBLocalFolder():
    """Return the path to a local PDB folder and folder structure specifier. 
    If a local folder is not set, ``None`` will be returned."""

    folder = SETTINGS.get('pdb_local_folder')
    if folder is not None:
        if isinstance(folder, str) and os.path.isdir(folder):
            return folder, SETTINGS.get('pdb_local_divided', True)
        else:
            LOGGER.warning('PDB local folder "{0:s}" is not a accessible.'
                           .format(folder))

def setPDBLocalFolder(folder, divided=False):
    """Set a local PDB folder.  Setting a local PDB folder will make 
    :func:`fetchPDB` function to seek that folder for presence of requested
    PDB files.  Also, files downloaded from `wwPDB <http://www.wwpdb.org/>`_ 
    FTP servers will be saved in this folder.  This may help users to store 
    PDB files in a single place and have access to them in different working 
    directories.
    
    If *divided* is ``True``, the divided folder structure of wwPDB servers 
    will be assumed when reading from and writing to the local folder.  For 
    example, a structure with identifier **1XYZ** will be present as 
    :file:`pdblocalfolder/yz/pdb1xyz.pdb.gz`. 
    
    If *divided* is ``False``, a plain folder structure will be expected and 
    adopted when saving files.  For example, the same structure will be 
    present as :file:`pdblocalfolder/1xyz.pdb.gz`.
    
    Finally, in either case, lower case letters will be used and compressed
    files will be stored."""
    
    if not isinstance(folder, str):
        raise TypeError('folder must be a string')
    assert isinstance(divided, bool), 'divided must be a boolean'
    if os.path.isdir(folder):
        folder = os.path.abspath(folder)
        LOGGER.info('Local PDB folder is set: "{0:s}"'.format(folder))
        if divided:
            LOGGER.info('When using local PDB folder, wwPDB divided '
                        'folder structure will be assumed.')
        else:
            LOGGER.info('When using local PDB folder, a plain folder structure '
                        'will be assumed.')
        SETTINGS['pdb_local_folder'] = folder
        SETTINGS['pdb_local_divided'] = divided
        SETTINGS.save()
    else:
        raise IOError('No such directory: "{0:s}"'.format(folder))

def getPDBMirrorPath():
    """Return the path to a local PDB mirror, or ``None`` if a mirror path is 
    not set."""

    path = SETTINGS.get('pdb_mirror_path')
    if isinstance(path, str):
        if os.path.isdir(path):
            return path
        else:
            LOGGER.warning('PDB mirror path "{0:s}" is not a accessible.'
                           .format(path))

def setPDBMirrorPath(path):
    """Set the path to a local PDB mirror."""
    
    if not isinstance(path, str):
        raise TypeError('path must be a string')
    if os.path.isdir(path):
        path = os.path.abspath(path)
        LOGGER.info('Local PDB mirror path is set: "{0:s}"'.format(path))
        SETTINGS['pdb_mirror_path'] = path
        SETTINGS.save()
    else:
        raise IOError('No such directory: "{0:s}"'.format(path))

def setWWPDBFTPServer(key):
    """Set the `wwPDB <http://www.wwpdb.org/>`_ FTP server used for downloading
    PDB structures when needed.  Use one of the following keywords for setting 
    a different server.
    
    +---------------------------+-----------------------------+
    | wwPDB FTP server          | *Key* (case insensitive)    |
    +===========================+=============================+
    | RCSB PDB (USA) (default)  | RCSB, USA, US               |
    +---------------------------+-----------------------------+
    | PDBe (Europe)             | PDBe, Europe, Euro, EU      |
    +---------------------------+-----------------------------+
    | PDBj (Japan)              | PDBj, Japan, Jp             |
    +---------------------------+-----------------------------+
    """
    
    server = WWPDB_FTP_SERVERS.get(key.lower())
    if server is not None:
        SETTINGS['wwpdb_ftp'] = server
        SETTINGS.save()
    else:
        LOGGER.warning('{0:s} is not a valid key.'.format(key))

def getWWPDBFTPServer():
    """Return a tuple containing name, host, and path of the currently 
    set `wwPDB <http://www.wwpdb.org/>`_ FTP server."""
    
    server = SETTINGS.get('wwpdb_ftp', None)
    if server is None:
        LOGGER.warning('A wwPDB FTP server is not set, default FTP server '
                       'RCSB PDB is used. Use `setWWPDBFTPServer` function '
                       'to set a server close to your location.')
        return _WWPDB_RCSB
    else:
        if server[2].endswith('data/structures/divided/pdb/'):
            return (server[0], server[1], 
                    server[2][:-len('data/structures/divided/pdb/')])
        else:
            return server

def fetchPDB(pdb, folder='.', compressed=True, copy=False, **kwargs):
    """Retrieve PDB, PDBML, or mmCIF file(s) for specified *pdb* identifier(s).  
    *pdb* may be a string or a list.  The function will return a filename or a 
    list of filenames depending on input (see :ref:`fetchpdb` for examples).  

    If *compressed* is ``False``, all files will be decompressed.  If *copy* is 
    ``True``, all files from local PDB mirror will copied to the user specified 
    *folder*.  *format* keyword argument can be used to retrieve `PDBML 
    <http://pdbml.pdb.org/>`_ and `mmCIF <http://mmcif.pdb.org/>`_ files:  
    ``format="cif"`` will fetch an mmCIF file (e.g. :file:`1XXX.cif.gz`), 
    similarly ``format="xml"`` will fetch a PDBML file.  If PDBML header file 
    is desired, ``format="xml", noatom=True`` will do the job (e.g. 
    :file:`1XXX-noatom.xml.gz`)
    
    The order of file search operations are as follows:  First, files are 
    sought in *folder*.  Second, local PDB mirror will be sought, if one is 
    set by the user (see :func:`setPDBMirrorPath`).  Then, local PDB folder
    will be sought, if one is  set by the user (see :func:`setPDBLocalFolder`).
    Finally, if files are not found locally, they will be downloaded one of 
    wwPDB FTP servers (use :func:`setWWPDBFTPServer` to specify one close to 
    you)."""
    
    if isinstance(pdb, str):
        identifiers = [pdb]
    elif isinstance(pdb, list):
        identifiers = pdb
    else:
        raise TypeError('pdb may be a string or a list of strings')
        
    assert isinstance(folder, str), 'folder must be a string'
    assert isinstance(compressed, bool), 'compressed must be a boolean'
    assert isinstance(copy, bool), 'copy must be a boolean'
    format = kwargs.pop('format', 'pdb')
    assert isinstance(format, str), 'format must be a string'
    format = format.lower()
    assert format in _PDB_FORMATS, '"{0:s}" is not valid format'.format(format)
    noatom = kwargs.pop('noatom', False) 
    assert isinstance(noatom, bool), 'noatom must be a boolean'
    if kwargs:
        raise TypeError('"{0:s}" is not a valid keyword argument for this' 
                        'function'.format(kwargs.iterkeys().next()))
    if folder != '.':
        folder = makePath(folder)
    if not os.access(folder, os.W_OK):
        raise IOError('permission to write in {0:s} is denied, please '
                      'specify another folder'.format(folder))
    
    filenames = []
    exists = 0
    success = 0
    failure = 0
    download = False
    if format == 'pdb':
        divided = 'data/structures/divided/pdb'
        pdbext = '.ent.gz'
        extensions = ['.ent', '.pdb'] # '.pdb' should be the last item
        prefix = 'pdb'
    elif format == 'xml':
        if noatom:
            divided = 'data/structures/divided/XML-noatom'
            pdbext = '-noatom.xml.gz'
            extensions = ['-noatom.xml']
        else:
            divided = 'data/structures/divided/XML'
            pdbext = '.xml.gz'
            extensions = ['.xml']
        prefix = ''
    else:
        divided = 'data/structures/divided/mmCIF'
        pdbext = '.cif.gz'
        extensions = ['.cif'] # '.pdb' should be the last item
        prefix = ''
    
    pdbfnmap = {}
    for extension in extensions:
        for pdbfn in glob(os.path.join(folder, '*' + extension + '*')): 
            if os.path.splitext(pdbfn)[1] in _PDB_EXTENSIONS:
                pdbfnmap[os.path.split(pdbfn)[1].split('.')[0].lower()] = pdbfn
        for pdbfn in glob(os.path.join(folder, '*' + extension.upper() + '*')):
            if os.path.splitext(pdbfn)[1] in _PDB_EXTENSIONS:
                pdbfnmap[os.path.split(pdbfn)[1].split('.')[0].lower()] = pdbfn
                
    for i, pdbid in enumerate(identifiers):
        # Check validity of identifiers
        if not isinstance(pdbid, str):
            LOGGER.debug('{0:s} is not a valid identifier.'.format(pdbid))
            filenames.append(None)
            failure += 1 
            continue
        pdbid = pdbid.strip().lower()
        if not (len(pdbid) == 4 and pdbid.isalnum()):
            LOGGER.debug('{0:s} is not a valid identifier.'.format(pdbid))
            filenames.append(None)
            failure += 1 
            continue
        # Check if file exists in working directory
        identifiers[i] = pdbid
        if noatom:
            fn = pdbfnmap.get(pdbid + '-noatom', None)
        else:
            fn = pdbfnmap.get(pdbid, None) or pdbfnmap.get('pdb'+pdbid, None)
        if fn:
            fn = relpath(fn)
            if not compressed:
                temp, ext = os.path.splitext(fn) 
                if ext == '.gz':
                    fn = gunzip(fn, temp)
            filenames.append(fn)
            LOGGER.debug('{0:s} ({1:s}) is found in the working directory.'
                         .format(pdbid, fn))
            exists += 1
            continue
        # Check the PDB mirror
        mirror_path = getPDBMirrorPath()
        if mirror_path is not None and os.path.isdir(mirror_path):
            fn = os.path.join(mirror_path, divided, pdbid[1:3], 
                              prefix + pdbid + pdbext)
            if os.path.isfile(fn):
                if copy or not compressed:
                    if compressed:
                        filename = os.path.join(folder, pdbid + extension + 
                                                        '.gz')
                        shutil.copy(fn, filename)
                    else:
                        filename = os.path.join(folder, pdbid + extension)
                        gunzip(fn, filename)
                    filenames.append(filename)
                    LOGGER.debug('{0:s} copied from local mirror ({1:s})'
                                 .format(pdbid, filename))
                    success += 1
                else:
                    filenames.append(fn)
                    
                    LOGGER.debug('{0:s} ({1:s}...{2:s}) is found in the local '
                                'mirror.'.format(pdbid, 
                                fn[:fn[1:].index(os.path.sep)+2], fn[-15:]))
                    exists += 1
                continue
        # Check the PDB mirror
        local_folder = getPDBLocalFolder()
        if format and local_folder:
            local_folder, is_divided = local_folder
            if is_divided:
                fn = os.path.join(local_folder, pdbid[1:3], 
                                  'pdb' + pdbid + '.pdb.gz')
            else:
                fn = os.path.join(local_folder, pdbid + '.pdb.gz')
                
            if os.path.isfile(fn):
                if copy or not compressed:
                    if compressed:
                        filename = os.path.join(folder, pdbid + extension + 
                                                        '.gz')
                        shutil.copy(fn, filename)
                    else:
                        filename = os.path.join(folder, pdbid + extension)
                        gunzip(fn, filename)
                    filenames.append(filename)
                    LOGGER.debug('{0:s} copied from local PDB folder ({1:s})'
                                 .format(pdbid, filename))
                    success += 1
                else:
                    filenames.append(fn)
                    
                    LOGGER.debug('{0:s} ({1:s}...{2:s}) is found in the PDB '
                                'local folder.'.format(pdbid, 
                                fn[:fn[1:].index(os.path.sep)+2], fn[-15:]))
                    exists += 1
                continue

        filenames.append(pdbid)
        download = True
    if download:
        from ftplib import FTP
        ftp_name, ftp_host, ftp_path = getWWPDBFTPServer()
        LOGGER.debug('Connecting wwPDB FTP server {0:s}.'.format(ftp_name))
        if format == 'pdb' and not copy and local_folder:
            folder = local_folder
            compressed = True
            if is_divided:
                getfn = lambda folder, pdbid, ext: \
                    os.path.join(makePath(os.path.join(local_folder, 
                                            pdbid[1:3])), 'pdb' + pdbid + ext)
            else:
                getfn = lambda folder, pdbid, ext: os.path.join(folder,
                                                                pdbid + ext)
                
        else: 
            getfn = lambda folder, pdbid, ext: os.path.join(folder, 
                                                            pdbid + ext)
        try:
            ftp = FTP(ftp_host)
        except Exception as error:
            raise type(error)('FTP connection problem, potential reason: '
                              'no internet connectivity')
        else:
            #ftp_path = os.path.join(ftp_path, divided)
            ftp.login('')
            for i, pdbid in enumerate(identifiers):
                if pdbid != filenames[i]:
                    continue
                filename = getfn(folder, pdbid, extension)
                if compressed:
                    filename += '.gz'

                pdbfile = open(filename, 'w+b')
                fn = prefix + pdbid + pdbext
                try:
                    ftp.cwd(ftp_path)
                    ftp.cwd(divided)
                    ftp.cwd(pdbid[1:3])
                    ftp.retrbinary('RETR ' + fn, pdbfile.write)
                except Exception as error:
                    pdbfile.close()
                    os.remove(filename)
                    if fn in ftp.nlst():
                        LOGGER.debug('{0:s} download failed ({1:s}). It '
                                     'is possible that you don\'t have '
                                     'rights to download .gz files in the '
                                     'current network.'.format(pdbid, 
                                     str(error)))
                    else:
                        LOGGER.debug('{0:s} download failed. {1:s} does not '
                                     'exist on {2:s}.'
                                     .format(fn, pdbid, ftp_host))
                    failure += 1
                    filenames[i] = None 
                else:
                    pdbfile.close()
                    if not compressed:
                        gunzip(filename)
                    filename = relpath(filename)
                    LOGGER.debug('{0:s} downloaded ({1:s})'
                                 .format(pdbid, filename))
                    success += 1
                    filenames[i] = filename
            ftp.quit()
    if len(identifiers) == 1:
        return filenames[0]    
    else:
        LOGGER.info('PDB download completed ({2:d} found, '
                    '{0:d} downloaded, {1:d} failed).'
                    .format(success, failure, exists))
        return filenames

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
    


def loadPDBClusters(sqid=None):
    """Load previously fetched PDB sequence clusters from disk to memory."""

    PDB_CLUSTERS_PATH = os.path.join(getPackagePath(), 'pdbclusters')
    if sqid is None:
        sqid_list = PDB_CLUSTERS.keys()
        LOGGER.info('Loading all PDB sequence clusters.')
    else:
        assert isinstance(sqid, int), 'sqid must be an integer' 
        if sqid not in PDB_CLUSTERS:
            keys = PDB_CLUSTERS.keys()
            keys.sort()
            raise ValueError('PDB cluster data is not available for sequence '
                             'identity {0:d}%, try one of {1:s}'
                             .format(sqid, ', '.join([str(x) for x in keys])))
        LOGGER.info('Loading PDB sequence clusters for sequence identity '
                    '{0:d}.'.format(sqid))
        sqid_list = [sqid]
    global PDB_CLUSTERS_UPDATE_WARNING
    for sqid in sqid_list:
        filename = os.path.join(PDB_CLUSTERS_PATH, 
                                'bc-{0:d}.out.gz'.format(sqid))
        if not os.path.isfile(filename):
            raise IOError('Local copy of PDB sequence clusters is not found, '
                          'call `fetchPDBClusters` function.')
        if PDB_CLUSTERS_UPDATE_WARNING:
            import time
            diff = (time.time() - os.path.getmtime(filename)) / 604800.
            if diff > 1.:
                LOGGER.warning('PDB sequence clusters are {0:.1f} week(s) old,'
                               ' call `fetchPDBClusters` to receive updates.'
                               .format(diff))
                PDB_CLUSTERS_UPDATE_WARNING = False
        inp = openFile(filename)
        PDB_CLUSTERS[sqid] = inp.read()
        inp.close()

def getPDBCluster(pdb, ch, sqid=95):
    """Return the PDB sequence cluster for chain *ch* in structure *pdb*
    that chains sharing sequence identity *sqid* or more.  PDB sequence cluster
    will be returned in the form of a list of tuples, e.g. 
    ``[('1XXX', 'A'), ('2YYY', 'A')]``.  Note that PDB clusters chains, so
    the same PDB identifier may appear twice in the same cluster if the 
    corresponding chain is present in the structure twice.    
    
    Before this function is used, :func:`fetchPDBClusters` needs to be called. 
    This function will load the PDB sequence clusters for *sqid* automatically 
    using :func:`loadPDBClusters`."""

    assert isinstance(pdb, str) and len(pdb) == 4, \
        'pdb must be 4 char long string'
    assert isinstance(ch, str) and len(ch) == 1, \
        'ch must be a one char long string'
    assert isinstance(sqid, int), 'sqid must be an integer'
    PDB_CLUSTERS_PATH = os.path.join(getPackagePath(), 'pdbclusters')
    if sqid not in PDB_CLUSTERS:
        keys = PDB_CLUSTERS.keys()
        keys.sort()
        raise ValueError('PDB cluster data is not available for sequence '
                         'identity {0:d}%, try one of {1:s}'
                         .format(sqid, ', '.join([str(x) for x in keys])))
    clusters = PDB_CLUSTERS[sqid]
    if clusters is None: 
        loadPDBClusters(sqid)
        clusters = PDB_CLUSTERS[sqid]
    pdb_ch = pdb.upper() + '_' + ch.upper()
    index = clusters.index(pdb_ch)
    maxlen = clusters.index('\n') 
    end = clusters.find('\n', index)
    start = clusters.rfind('\n', index-maxlen, end)+1
    cluster = clusters[start:end]
    return [tuple(item.split('_')) for item in cluster.split()] 

def fetchPDBClusters():
    """Downloads PDB sequence clusters.  PDB sequence clusters are results of 
    the weekly clustering of protein chains in the PDB generated by blastclust. 
    They are available at FTP site: ftp://resources.rcsb.org/sequence/clusters/
    
    This function will download about 10 Mb of data and save it after 
    compressing in your home directory in :file:`.prody/pdbclusters`.
    Compressed files will be less than 4 Mb in size.  Cluster data can 
    be loaded using :func:`loadPDBClusters` function and be accessed 
    using :func:`getPDBCluster`."""
    
    import urllib2
    PDB_CLUSTERS_PATH = os.path.join(getPackagePath(), 'pdbclusters')
    if not os.path.isdir(PDB_CLUSTERS_PATH):
        os.mkdir(PDB_CLUSTERS_PATH)
    LOGGER.progress('Downloading sequence clusters', len(PDB_CLUSTERS))
    count = 0
    for i, x in enumerate(PDB_CLUSTERS.keys()):
        filename = 'bc-{0:d}.out'.format(x)
        url = ('ftp://resources.rcsb.org/sequence/clusters/' + filename)
        try:
            inp = urllib2.urlopen(url)
        except urllib2.HTTPError:
            LOGGER.warning('Clusters at {0:d}% sequence identity level could '
                           'not be downloaded.')
            continue
        else:
            out = openFile(filename+'.gz', 'w', folder=PDB_CLUSTERS_PATH) 
            out.write(inp.read())
            inp.close()
            out.close()
            count += 1
        LOGGER.update(i)
    LOGGER.clear()
    if len(PDB_CLUSTERS) == count:
        LOGGER.info('All PDB clusters were downloaded successfully.')
    elif count == 0:
        LOGGER.warning('PDB clusters could not be downloaded.')
