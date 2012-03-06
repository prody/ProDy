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

"""This module defines functions for accessing wwPDB FTP servers."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path
import shutil
from glob import glob

from prody.tools import makePath, gunzip, relpath
from localpdb import getPDBLocalFolder, getPDBMirrorPath

__all__ = ['getWWPDBFTPServer', 'setWWPDBFTPServer', 'fetchPDB',]

pkg = __import__(__package__)
LOGGER = pkg.LOGGER
SETTINGS = pkg.SETTINGS
           
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
    assert format in _PDB_FORMATS, '{0:s} is not valid format'.format(
                                                                repr(format))
    noatom = kwargs.pop('noatom', False) 
    assert isinstance(noatom, bool), 'noatom must be a boolean'
    if kwargs:
        raise TypeError('{0:s} is not a valid keyword argument for this' 
                        'function'.format(repr(kwargs.iterkeys().next())))
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
