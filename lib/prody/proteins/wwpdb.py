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

"""This module defines functions for accessing wwPDB servers."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from os import getcwd
from glob import glob
from os.path import sep as pathsep
from os.path import isdir, isfile, join, split, splitext

from prody import LOGGER, SETTINGS
from prody.utilities import makePath, gunzip, relpath, copyFile, openURL


__all__ = ['wwPDBServer', 'getWWPDBFTPServer', 'setWWPDBFTPServer',
           'fetchPDBviaFTP', 'fetchPDBviaHTTP']

           
_PDB_EXTENSIONS = set(['.pdb', '.PDB', '.gz', '.GZ', '.ent', '.ENT', 
                       '.pdb.gz', '.PDB.GZ', '.ent.gz', '.ENT.GZ'])
_XML_EXTENSIONS = set(['.xml', '.XML', '.xml.gz', '.XML.GZ'])
_CIF_EXTENSIONS = set(['.cif', '.CIF', '.cif.gz', '.CIF.GZ'])

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

_URL_US = lambda pdb: ('http://www.rcsb.org/pdb/files/%s.pdb.gz' % 
                       pdb.upper())
_URL_EU = lambda pdb: ('http://www.ebi.ac.uk/pdbe-srv/view/files/%s.ent.gz' % 
                       pdb.lower())
_URL_JP = lambda pdb: ('http://www.pdbj.org/pdb_all/pdb%s.ent.gz' % 
                       pdb.lower())
WWPDB_HTTP_URL = {
    'rcsb'   : _URL_US,
    'usa'    : _URL_US,
    'us'     : _URL_US,
    'pdbe'   : _URL_EU,
    'euro'   : _URL_EU,
    'europe' : _URL_EU,
    'eu'     : _URL_EU,
    'pdbj'   : _URL_JP,
    'japan'  : _URL_JP,
    'jp'     : _URL_JP,
}

def wwPDBServer(*key):
    """Set/get `wwPDB`_ FTP/HTTP server location used for downloading PDB 
    structures.  Use one of the following keywords for setting a server:
    
    +---------------------------+-----------------------------+
    | wwPDB FTP server          | *Key* (case insensitive)    |
    +===========================+=============================+
    | RCSB PDB (USA) (default)  | RCSB, USA, US               |
    +---------------------------+-----------------------------+
    | PDBe (Europe)             | PDBe, Europe, Euro, EU      |
    +---------------------------+-----------------------------+
    | PDBj (Japan)              | PDBj, Japan, Jp             |
    +---------------------------+-----------------------------+
    
    .. _wwPDB: http://www.wwpdb.org/"""
    
    if not key:
        server = SETTINGS.get('wwpdb', None)
        if server is None:
            LOGGER.warn('A wwPDB server is not set, default server RCSB PDB '
                        'is used. You may use `wwPDBServer` function to set '
                        'a physically close server location.')
        return server
    elif len(key) == 1:
        try:
            key = key[0].lower()
        except AttributeError:
            raise TypeError('key must be a string')
        if key in WWPDB_FTP_SERVERS:
            SETTINGS['wwpdb'] = key
            SETTINGS.save()
        else:
            raise ValueError('{0:s} is not a valid key'.format(repr(key)))
    else:
        raise TypeError('one key argument is expected, {0:d} given'
                        .format(len(key)))


def setWWPDBFTPServer(key):
    """Deprecated for removal in v1.4, use :func:`wwPDBServer` instead."""
    
    return wwPDBServer(key)


def getWWPDBFTPServer():
    """Deprecated for removal in v1.4, use :func:`wwPDBServer` instead."""
    
    return wwPDBServer()


def checkIdentifiers(*pdb):
    """Check whether items in *identifiers* look like PDB identifiers.  Replace
    those that do not look like one with **None** or raise a :exc:`TypeError` 
    when a non string is encountered."""
    
    identifiers = []
    for pid in pdb:
        try:        
            pid = pid.strip().lower()
        except AttributeError:
            LOGGER.warn('{0:s} is not a valid identifier.'.format(repr(pid)))
            identifiers.append(None)
        else:        
            if not (len(pid) == 4 and pid.isalnum()):
                LOGGER.warn('{0:s} is not a valid identifier.'
                            .format(repr(pid)))
                identifiers.append(None)
            else:
                identifiers.append(pid)
    return identifiers


def fetchPDBviaFTP(*pdb, **kwargs):
    """Retrieve PDB (default), PDBML, or mmCIF file(s) for specified *pdb* 
    identifier(s) and return path(s).  Downloaded files will be stored in 
    local PDB folder, if one is set using :meth:`.pathPDBFolder`, and copied
    into *folder*, if specified by the user.  If no destination folder is 
    specified, files will be saved in the current working directory.  If 
    *compressed* is **False**, decompressed files will be copied into 
    *folder*.  *format* keyword argument can be used to retrieve
    `PDBML <http://pdbml.pdb.org/>`_ and `mmCIF <http://mmcif.pdb.org/>`_ 
    files: ``format='cif'`` will fetch an mmCIF file, and ``format='xml'`` 
    will fetch a PDBML file.  If PDBML header file is desired, ``noatom=True`` 
    argument will do the job."""

    if kwargs.get('check', True):
        identifiers = checkIdentifiers(*pdb)
    else:
        identifiers = list(pdb)

    output_folder = kwargs.pop('folder', None)
    compressed = bool(kwargs.pop('compressed', True))
    format = str(kwargs.pop('format', 'pdb')).lower()
    if not format in _PDB_FORMATS:
        raise ValueError(repr(format) + ' is not valid format')
    noatom = bool(kwargs.pop('noatom', False)) 
    
    if format == 'pdb':
        ftp_divided = 'data/structures/divided/pdb'
        ftp_pdbext = '.ent.gz'
        ftp_prefix = 'pdb'
        extension = '.pdb'
    elif format == 'xml':
        if noatom:
            ftp_divided = 'data/structures/divided/XML-noatom'
            ftp_pdbext = '-noatom.xml.gz'
            extension = '-noatom.xml'
        else:
            ftp_divided = 'data/structures/divided/XML'
            ftp_pdbext = '.xml.gz'
            extension = '.xml'
        ftp_prefix = ''
    elif format == 'cif':
        ftp_divided = 'data/structures/divided/mmCIF'
        ftp_pdbext = '.cif.gz'
        ftp_prefix = ''
        extension = '.cif'
    else:
        raise ValueError('{0:s} is not a recognized format'
                         .format(repr(format)))

    local_folder = pathPDBFolder()

    if format == 'pdb' and local_folder:
        local_folder, is_divided = local_folder
        if is_divided:
            getPath = lambda pdb: join(makePath(join(local_folder, pdb[1:3])), 
                                       'pdb' + pdb + '.pdb.gz')
        else:
            getPath = lambda pdb: join(local_folder, pdb + '.pdb.gz')
        if output_folder is None:
            second = lambda filename, pdb: filename
        else:      
            if compressed:
                second = lambda filename, pdb: (copyFile(filename,
                            join(output_folder, pdb + extension + '.gz')))
            else:
                second = lambda filename, pdb: gunzip(filename,
                            join(output_folder, pdb + extension))
            
    else:
        if output_folder is None:
            output_folder = getcwd()
        if compressed:
            getPath = lambda pdb: join(output_folder, pdb + extension + '.gz')
            second = lambda filename, pdb: filename
        else:
            getPath = lambda pdb: join(output_folder, pdb + extension)
            second = lambda filename, pdb: gunzip(getPath(pdb), getPath(pdb))
        
    
    ftp_name, ftp_host, ftp_path = WWPDB_FTP_SERVERS[wwPDBServer() or 'us']
    LOGGER.debug('Connecting wwPDB FTP server {0:s}.'.format(ftp_name))

    from ftplib import FTP
    try:
        ftp = FTP(ftp_host)
    except Exception as error:
        raise type(error)('FTP connection problem, potential reason: '
                          'no internet connectivity')
    else:
        success = 0
        failure = 0
        filenames = []
        ftp.login('')
        for pdb in identifiers:
            if pdb is None:
                filenames.append(None)
            
            data = []
            ftp_fn = ftp_prefix + pdb + ftp_pdbext
            try:
                ftp.cwd(ftp_path)
                ftp.cwd(ftp_divided)
                ftp.cwd(pdb[1:3])
                ftp.retrbinary('RETR ' + ftp_fn, data.append)
            except Exception as error:
                if ftp_fn in ftp.nlst():
                    LOGGER.warn('{0:s} download failed ({1:s}). It is '
                                'possible that you do not have rights to '
                                'download .gz files in the current network.'
                                .format(pdb, str(error)))
                else:
                    LOGGER.warn('{0:s} download failed. {1:s} does not exist '
                                'on {2:s}.'.format(ftp_fn, pdb, ftp_host))
                failure += 1
                filenames.append(None)
            else:
                if len(data):
                    filename = getPath(pdb)
        
                    with open(filename, 'w+b') as pdbfile:
                        write = pdbfile.write
                        [write(block) for block in data]

                    filename = relpath(second(filename, pdb))
                    LOGGER.debug('{0:s} downloaded ({1:s})'
                                 .format(pdb, filename))
                    success += 1
                    filenames.append(filename)
                else:
                    LOGGER.warn('{0:s} download failed, reason unknown.'
                                .format(pdb))
                    failure += 1
                    filenames.append(None)

        ftp.quit()

    if len(identifiers) == 1:
        return filenames[0]    
    else:
        if kwargs.get('report', True):
            LOGGER.debug('PDB download via FTP completed ({0:d} downloaded, '
                         '{1:d} failed).'.format(success, failure))
        return filenames


def fetchPDBviaHTTP(*pdb, **kwargs):
    """Retrieve PDB file(s) for specified *pdb* identifier(s) and return 
    path(s).  Downloaded files will be stored in local PDB folder, if one 
    is set using :meth:`.pathPDBFolder`, and copied into *folder*, if 
    specified by the user.  If no destination folder is specified, files 
    will be saved in the current working directory.  If *compressed* is 
    **False**, decompressed files will be copied into *folder*."""

    if kwargs.get('check', True):
        identifiers = checkIdentifiers(*pdb)
    else:
        identifiers = list(pdb)

    output_folder = kwargs.pop('folder', None)
    compressed = bool(kwargs.pop('compressed', True))

    extension = '.pdb'
    local_folder = pathPDBFolder()
    if local_folder:
        local_folder, is_divided = local_folder
        if is_divided:
            getPath = lambda pdb: join(makePath(join(local_folder, pdb[1:3])), 
                                       'pdb' + pdb + '.pdb.gz')
        else:
            getPath = lambda pdb: join(local_folder, pdb + '.pdb.gz')
        if output_folder is None:
            second = lambda filename, pdb: filename
        else:      
            if compressed:
                second = lambda filename, pdb: (copyFile(filename,
                            join(output_folder, pdb + extension + '.gz')))
            else:
                second = lambda filename, pdb: gunzip(filename,
                            join(output_folder, pdb + extension))
            
    else:
        if output_folder is None:
            output_folder = getcwd()
        if compressed:
            getPath = lambda pdb: join(output_folder, pdb + extension + '.gz')
            second = lambda filename, pdb: filename
        else:
            getPath = lambda pdb: join(output_folder, pdb + extension)
            second = lambda filename, pdb: gunzip(getPath(pdb), getPath(pdb))
        
    
    getURL = WWPDB_HTTP_URL[wwPDBServer() or 'us']

    success = 0
    failure = 0
    filenames = []
    for pdb in identifiers:
        if pdb is None:
            filenames.append(None)
        try:
            handle = openURL(getURL(pdb))
        except Exception as err:
            LOGGER.warn('{0:s} download failed ({0:s}).'.format(pdb, str(err)))
            failure += 1
            filenames.append(None)
        else:
            data = handle.read() 
            if len(data):
                filename = getPath(pdb)
    
                with open(filename, 'w+b') as pdbfile:
                    write = pdbfile.write
                    [write(block) for block in data]

                filename = relpath(second(filename, pdb))
                LOGGER.debug('{0:s} downloaded ({1:s})'
                             .format(pdb, filename))
                success += 1
                filenames.append(filename)
            else:
                LOGGER.warn('{0:s} download failed, reason unknown.'
                            .format(pdb))
                failure += 1
                filenames.append(None)

    if len(identifiers) == 1:
        return filenames[0]    
    else:
        if kwargs.get('report', True):
            LOGGER.debug('PDB download via HTTP completed ({0:d} downloaded, '
                         '{1:d} failed).'.format(success, failure))
        return filenames

if __name__ == '__main__':
    
    pdbids = ['1mkp', '1zz2', 'nano']
    
    for gzip in [False, True]:
        fetchPDBviaFTP(*pdbids, compressed=gzip, folder='.')
        fetchPDBviaFTP(*pdbids, compressed=gzip, folder='.', format='cif')
        fetchPDBviaFTP(*pdbids, compressed=gzip, folder='.', format='xml')
        fetchPDBviaFTP(*pdbids, compressed=gzip, folder='.', format='xml', 
                       noatom=1)
        fetchPDBviaHTTP(*pdbids, compressed=gzip, folder='.')
    from glob import glob
    from os import remove
    for pdb in pdbids:
        fns = glob(pdb + '.*')
        print pdb, '>', ', '.join(fns)
        for fn in fns:
            remove(fn)
