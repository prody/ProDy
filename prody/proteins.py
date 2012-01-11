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

"""This module defines classes and functions to fetch, parse, and write 
structural data files, execute structural analysis programs, and to access 
and search structural databases, e.g. `ProteinDataBank <http://wwpdb.org>`_.

Protein structural data
===============================================================================

Access databases 
-------------------------------------------------------------------------------

Following functions are provided for access to protein structural data:

=========================  ====================================================
Function                   Description
=========================  ====================================================
:func:`blastPDB`           blast search NCBI PDB database
:func:`fetchPDB`           retrieve PDB/PDBML/mmCIF files from wwPDB
:func:`fetchLigandData`    retrieve ligand from Ligand-Expo 
:func:`fetchPDBClusters`   retrieve PDB sequence cluster data from wwPDB
:func:`getPDBCluster`      access PDB sequence clusters
:func:`setPDBLocalFolder`  set a local folder for storing PDB files
:func:`setPDBMirrorPath`   set a local PDB mirror path
:func:`setWWPDBFTPServer`  set a wwPDB FTP server for downloads 
:func:`getPDBLocalFolder`  get preset local PDB folder
:func:`getPDBMirrorPath`   get preset local PDB mirror path
:func:`getWWPDBFTPServer`  get preset wwPDB FTP server
=========================  ====================================================


Parse/write files
-------------------------------------------------------------------------------

Following ProDy functions are for parsing and writing atomic data:

======================  =======================================================
Function                Description
======================  =======================================================
:func:`parsePDB`        parse atomic data from files in :file:`.pdb` format
:func:`parsePSF`        parse atomic data from files in :file:`.psf` format
:func:`parsePQR`        parse atomic data from files in :file:`.pqr` format
:func:`parsePDBHeader`  parse header data from :file:`.pdb` files 
:func:`parsePDBStream`  parse atomic data from a stream in :file:`.pdb` format
:func:`parseDSSP`       parse structural data from :program:`dssp` output
:func:`parseSTRIDE`     parse structural data from :program:`stride` output
:func:`writePDB`        write atomic data to a file in :file:`.pdb` format
:func:`writePQR`        write atomic data to a file in :file:`.pqr` format
:func:`writePDBStream`  write atomic data in :file:`.pdb` format to a stream
======================  =======================================================

.. seealso::
   
   Atom data (coordinates, atom names, residue names, etc.) parsed from 
   PDB/PSF/PQR files are stored in :class:`~prody.atomic.AtomGroup` instances.  
   See :mod:`~prody.atomic` module documentation for more details. 


Store data
-------------------------------------------------------------------------------

Following classes are for storing meta, structural, and/or search data: 

=======================  ======================================================
Function                 Description
=======================  ======================================================
:class:`Chemical`        store PDB chemical (heterogen) component data
:class:`Polymer`         store PDB polymer (macromolecule) component data
:class:`PDBBlastRecord`  store and evaluate NCBI PDB blast search results 
=======================  ======================================================

Execute programs
-------------------------------------------------------------------------------

Following functions can be used to execute structural analysis programs from 
within Python:

========================  =====================================================
Function                  Description
========================  =====================================================
:func:`execDSSP`          execute :program:`dssp`
:func:`execSTRIDE`        execute :program:`stride`
:func:`performDSSP`       execute :program:`dssp` and parse results
:func:`performSTRIDE`     execute :program:`stride` and parse results
========================  =====================================================


Edit structures
-------------------------------------------------------------------------------

Following functions allow editing structures using structural data from PDB 
header records:

=========================  ====================================================
Function                   Description
=========================  ====================================================
:func:`assignSecstr`       add secondary structure data from header to atoms
:func:`buildBiomolecules`  build biomolecule data based on header records
=========================  ====================================================


Visualization
-------------------------------------------------------------------------------

:func:`showProtein` function can be used to take a quick look at protein 
structures. 

:mod:`prody.proteins`
===============================================================================

.. doctest::
    :hide:
        
    >>> from prody import *
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    

.. plot::
   :nofigs: 
   :context: 
    
   from prody import *
   import matplotlib.pyplot as plt
   import numpy as np
    
   plt.close('all')    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path
import os
import shutil
import sys
from glob import glob
from collections import defaultdict

plt = None
import numpy as np

from tools import *

BioBlast = None

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS
from prody.atomic import *


PDB_CLUSTERS = {30: None, 40: None, 50: None, 70: None, 
                90: None, 95: None, 100: None}
PDB_CLUSTERS_UPDATE_WARNING = True


__all__ = ['Chemical', 'Polymer', 'PDBBlastRecord',
           'assignSecondaryStructure', 'assignSecstr',
           'applyBiomolecularTransformations', 'buildBiomolecules',
           'blastPDB', 'fetchPDB', 
           'getPDBLocalFolder', 'getPDBMirrorPath', 'getWWPDBFTPServer', 
           'setPDBLocalFolder', 'setPDBMirrorPath', 'setWWPDBFTPServer',
           'parsePDBStream', 'parsePDB', 'parsePSF', 'parsePQR',
           'writePDBStream', 'writePDB', 'writePQR', 
           'parsePDBHeader',
           'fetchLigandData',
           'execDSSP', 'parseDSSP', 'performDSSP',
           'execSTRIDE', 'parseSTRIDE', 'performSTRIDE',
           'fetchPDBClusters', 'loadPDBClusters', 'getPDBCluster',
           'showProtein']

class PDBParseError(Exception):    
    pass

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
    If a local folder is not set, ``None`` will be returned.
    
    .. versionadded:: 0.9"""

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
    
    .. versionadded:: 0.9
    
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
    not set.
    
    .. versionadded:: 0.6.1"""

    path = SETTINGS.get('pdb_mirror_path')
    if isinstance(path, str):
        if os.path.isdir(path):
            return path
        else:
            LOGGER.warning('PDB mirror path "{0:s}" is not a accessible.'
                           .format(path))

def setPDBMirrorPath(path):
    """Set the path to a local PDB mirror.  
    
    .. versionadded:: 0.6.1"""
    
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
    PDB structures when needed.
    
    .. versionadded:: 0.6.1
    
    Use one of the following keywords for setting a different server.
    
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
    set `wwPDB <http://www.wwpdb.org/>`_ FTP server.
    
    .. versionadded:: 0.6.1"""
    
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

    .. versionchanged:: 0.8
       *compressed* and *copy* argument is introduced.  
    
    .. versionchanged:: 0.8.2
       When *compressed* is false, compressed files found in *folder* or 
       local PDB mirror are decompressed.
    
    .. versionadded:: 0.9
       File discovery is improved to handle a local PDB folder. See 
       :func:`setPDBLocalFolder` method for details.  
       *format* and *noatom* keyword arguments are added.

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

_parsePQRdoc = """
    :arg title: Title of the AtomGroup instance.  When ``None`` is passed,
        AtomGroup is named after the PDB filename.  
    :type title: str
    
    :arg ag: :class:`~prody.atomic.AtomGroup` instance for storing data parsed 
        from PDB file.  Number of atoms in *ag* and number of atoms parsed from
        the PDB file must be the same.  Atoms in *ag* and the PDB file must be 
        in the same order.  Non-coordinate data stored in *ag* will be 
        overwritten with those parsed from the file. 
    :type ag: :class:`~prody.atomic.AtomGroup`

    :arg chain: Chain identifiers for parsing specific chains, e.g. 
        ``chain='A'``, ``chain='B'``, ``chain='DE'``. By default all chains
        are parsed.        
    :type chain: str

    :arg subset: A predefined keyword to parse subset of atoms.
        Valid keywords are ``"calpha"`` (``"ca"``) or ``"backbone"`` 
        (``"bb"``), or ``None`` (read all atoms), e.g. ``subset='bb'``
    :type subset: str
"""

_parsePDBdoc = _parsePQRdoc + """
    :arg model: model index or None (read all models), 
        e.g. ``model=10``
    :type model: int, list

    :arg header: If ``True`` PDB header content will be parsed and returned.
    :type header: bool

    :arg altloc: If a location indicator is passed, such as ``'A'`` or ``'B'``, 
         only indicated alternate locations will be parsed as the single 
         coordinate set of the AtomGroup.  If ``True`` all alternate locations 
         will be parsed and each will be appended as a distinct coordinate set.
         Default is ``"A"``.
    :type altloc: str

    :arg biomol: If ``True``, return biomolecule obtained by transforming the
        coordinates using information from header section.
    :type biomol: False

    :arg secondary: If ``True``, parse the secondary structure information
        from header section and assign data to atoms.
    :type secondary: False
        

    If ``model=0`` and ``header=True``, return header 
    dictionary only.
    
    .. versionchanged:: 0.6
       Default behavior for parsing alternate locations have changed. 
       Alternate locations indicated by ``A`` are parsed.
    
    .. versionchanged:: 0.7.1
       *name* is now a keyword argument.  *biomol* and *secondary* keyword 
       arguments makes the parser return the biological molecule and/or
       assign secondary structure information.
       
    .. versionchanged:: 0.8.1
       User can pass an :class:`~prody.atomic.AtomGroup` instance as *ag*.

    .. versionchanged:: 0.8.2
       *chain* and *subset* arguments are appended to atom group name, e.g. for 
       ``('1mkp', chain='A', subset='calpha')`` name will be ``"1mkp_A_ca"``.

    .. versionchanged:: 0.9
       *name* keyword argument is renamed as *title* argument..


    Note that this function does not evaluate ``CONECT`` records.
    
    """
    
_PDBSubsets = {'ca': 'ca', 'calpha': 'ca', 'bb': 'bb', 'backbone': 'bb'}

def parsePDB(pdb, **kwargs):
    """Return an :class:`~prody.atomic.AtomGroup` and/or dictionary containing 
    header data parsed from a PDB file. 
    
    This function extends :func:`parsePDBStream`.
    
    |example| See :ref:`parsepdb` for a detailed example.
    
    :arg pdb: A valid PDB identifier or filename.  
        If needed, PDB files are downloaded using :func:`fetchPDB()` function.
    """
    
    title = kwargs.get('title', kwargs.get('name'))
    if not os.path.isfile(pdb):
        if len(pdb) == 4 and pdb.isalnum():
            if title is None:
                title = pdb.lower()
                kwargs['title'] = title
            filename = fetchPDB(pdb)
            if filename is None:
                raise IOError('PDB file for {0:s} could not be downloaded.'
                              .format(pdb))
            pdb = filename
        else:
            raise IOError('{0:s} is not a valid filename or a valid PDB '
                          'identifier.'.format(pdb))
    if title is None:
        fn, ext = os.path.splitext(os.path.split(pdb)[1])
        if ext == '.gz':
            fn, ext = os.path.splitext(fn)
        title = fn.lower()
        kwargs['title'] = title
    pdb = openFile(pdb)
    result = parsePDBStream(pdb, **kwargs)
    pdb.close()
    return result

parsePDB.__doc__ += _parsePDBdoc
    
def parsePDBStream(stream, **kwargs):
    """Return an :class:`~prody.atomic.AtomGroup` and/or 
    dictionary containing header data parsed from a stream of PDB lines. 
    
    :arg stream: Anything that implements the method readlines() 
        (e.g. :class:`file`, buffer, stdin).
    """
    
    model = kwargs.get('model')
    header = kwargs.get('header', False)
    assert isinstance(header, bool), 'header must be a boolean'
    chain = kwargs.get('chain')
    subset = kwargs.get('subset')
    altloc = kwargs.get('altloc', 'A')
    if model is not None:
        if isinstance(model, int):
            if model < 0:
                raise ValueError('model must be greater than 0')
        else:
            raise TypeError('model must be an integer, {0:s} is invalid'
                            .format(str(model)))
    title_suffix = ''
    if subset is not None: 
        if not isinstance(subset, str):
            raise TypeError('subset must be a string')
        elif subset.lower() not in _PDBSubsets:
            raise ValueError('"{0:s}" is not a valid subset'.format(subset))
        title_suffix = '_' + _PDBSubsets[subset]
    if chain is not None:
        if not isinstance(chain, str):
            raise TypeError('chain must be a string')
        elif len(chain) == 0:
            raise ValueError('chain must not be an empty string')
        title_suffix = '_' + chain + title_suffix
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, prody.AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    else:
        ag = prody.AtomGroup(str(kwargs.get('title', kwargs.get('name', 
                                                'Unknown'))) + title_suffix)
        n_csets = 0
    
    biomol = kwargs.get('biomol', False)
    secondary = kwargs.get('secondary', False)
    split = 0
    hd = None
    if model != 0:
        LOGGER.timeit()
        lines = stream.readlines()
        if header or biomol or secondary:
            hd, split = _getHeaderDict(lines)
        _parsePDBLines(ag, lines, split, model, chain, subset, altloc)
        if ag.numAtoms() > 0:
            LOGGER.timing('{0:d} atoms and {1:d} coordinate set(s) were '
                          'parsed in %.2fs.'.format(ag.numAtoms(), 
                           ag.numCoordsets() - n_csets))
        else:
            ag = None
            LOGGER.warning('Atomic data could not be parsed, please '
                           'check the input file.')
    elif header:
        hd, split = _getHeaderDict(stream)

    if secondary:
        try:
            ag = assignSecstr(hd, ag)
        except:
            raise PDBParseError('secondary structure assignments could not '
                                 'be made, check input file')
    if biomol:
        try:
            ag = buildBiomolecules(hd, ag)
        except:
            raise PDBParseError('biomolecule could not be generated, check'
                                 'input file')
        else:
            if isinstance(ag, list):
                LOGGER.info('Biomolecular transformations were applied, {0:d} '
                            'biomolecule(s) are returned.'.format(len(ag)))
            else:
                LOGGER.info('Biomolecular transformations were applied to the '
                            'coordinate data.')
    if model != 0:
        if header:
            return ag, hd
        else:
            return ag
    else:
        return hd

parsePDBStream.__doc__ += _parsePDBdoc


def parsePQR(filename, **kwargs):
    """Return an :class:`~prody.atomic.AtomGroup` containing data parsed 
    from PDB lines. 
    
    :arg filename: a PQR filename
    :type filename: str
    """
    
    title = kwargs.get('title', kwargs.get('name'))
    model = 1
    header = False
    chain = kwargs.get('chain')
    subset = kwargs.get('subset')
    altloc = kwargs.get('altloc', 'A')
    if not os.path.isfile(filename):
        raise IOError('No such file: "{0:s}"'.format(filename))
    if title is None:
        fn, ext = os.path.splitext(os.path.split(filename)[1])
        if ext == '.gz':
            fn, ext = os.path.splitext(fn)
        title = fn.lower()
    title_suffix = ''
    if subset is not None:
        if not isinstance(subset, str):
            raise TypeError('subset must be a string')
        elif subset.lower() not in _PDBSubsets:
            raise ValueError('"{0:s}" is not a valid subset'.format(subset))
        title_suffix = '_' + _PDBSubsets[subset]
    if chain is not None:
        if not isinstance(chain, str):
            raise TypeError('chain must be a string')
        elif len(chain) == 0:
            raise ValueError('chain must not be an empty string')
        title_suffix = '_' + chain + title_suffix
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, prody.AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    else:
        ag = prody.AtomGroup(title + title_suffix)
        n_csets = 0
        
    pqr = openFile(filename)
    lines = pqr.readlines()
    pqr.close()
    LOGGER.timeit()
    ag = _parsePDBLines(ag, lines, split=0, model=1, chain=chain, 
                        subset=subset, altloc_torf=False, format='pqr')
    if ag.numAtoms() > 0:
        LOGGER.timing('{0:d} atoms and {1:d} coordinate sets were '
                      'parsed in %.2fs.'.format(ag.numAtoms(), 
                         ag.numCoordsets() - n_csets))
        return ag
    else:
        return None

parsePQR.__doc__ += _parsePQRdoc

def _parsePDBLines(atomgroup, lines, split, model, chain, subset, 
                   altloc_torf, format='pdb'):
    """Return an AtomGroup. See also :func:`parsePDBStream()`.
    
    :arg lines: PDB/PQR lines 
    :arg split: starting index for coordinate data lines"""
    
    format = format.upper()
    if format == 'PDB':
        isPDB = True
    else:
        isPDB = False
    
    if subset is not None:
        subset = subset.lower()
        if subset in ('calpha', 'ca'):
            subset = set(('CA',))
        elif subset in ('backbone', 'bb'):
            subset = set(prody.getBackboneAtomNames())
        only_subset = True
        protein_resnames = set(prody.getKeywordResnames('protein'))
    else:
        only_subset = False
    if chain is None:
        only_chains = False
    else:
        only_chains = True
    onlycoords = False
    n_atoms = atomgroup.numAtoms()
    if n_atoms > 0:
        asize = n_atoms
    else:
        # most PDB files contain less than 99999 atoms
        asize = min(len(lines) - split, 99999)
    alength = asize
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['chain'].dtype)
    hetero = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['hetero'].dtype)
    altlocs = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['altloc'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['serial'].dtype)
    if isPDB:
        segnames = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['segment'].dtype)
        elements = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['element'].dtype)
        bfactors = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['beta'].dtype)
        occupancies = np.zeros(asize, 
                               dtype=ATOMIC_DATA_FIELDS['occupancy'].dtype)
        secondary = None
        anisou = None
        siguij = None
    else:
        charges = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['charge'].dtype)
        radii = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['radius'].dtype)
        type_ = 'PDB'
        
    asize = 2000 # increase array length by this much when needed 
        
    start = split
    stop = len(lines)
    nmodel = 0
    # if a specific model is requested, skip lines until that one
    if isPDB and model is not None and model != 1:
        for i in range(split, len(lines)):
            if lines[i][:5] == 'MODEL':
                nmodel += 1
                if model == nmodel:
                    start = i+1 
                    stop = len(lines)
                    break
        if nmodel != model:
            raise PDBParseError('model {0:d} is not found'.format(model))
    if isinstance(altloc_torf, str): 
        if altloc_torf.strip() != 'A':
            LOGGER.info('Parsing alternate locations {0:s}.'
                        .format(altloc_torf))
            which_altlocs = ' ' + ''.join(altloc_torf.split())
        else:
            which_altlocs = ' A'
        altloc_torf = False
    else:
        which_altlocs = ' A'
        altloc_torf = True
        
    acount = 0
    altloc = defaultdict(list)
    i = start
    END = False
    while i < stop:
        line = lines[i]
        startswith = line[0:6]
        
        if startswith == 'ATOM  ' or startswith == 'HETATM':
            if only_subset:
                atomname = line[12:16].strip()
                resname = line[17:21].strip()
                if not (atomname in subset and resname in protein_resnames): 
                    i += 1
                    continue
            else:
                atomname = line[12:16]
                resname = line[17:21]
                
            chid = line[21]
            if only_chains:
                if not chid in chain:
                    i += 1
                    continue
            alt = line[16]
            if alt not in which_altlocs:
                altloc[alt].append((line, i))
                i += 1
                continue
            try:
                coordinates[acount, 0] = line[30:38]
                coordinates[acount, 1] = line[38:46]
                coordinates[acount, 2] = line[46:54]
            except:
                if acount >= n_atoms > 0:
                    if nmodel ==0:
                        raise ValueError(format + 'file and AtomGroup ag must '
                                         'have same number of atoms')
                    LOGGER.warning('Discarding model {0:d}, which contains '
                                   'more atoms than first model does.'
                                   .format(nmodel+1))
                    acount = 0
                    nmodel += 1
                    coordinates = np.zeros((n_atoms, 3), dtype=float)
                    while lines[i][:6] != 'ENDMDL':
                        i += 1
                else:
                    raise PDBParseError('invalid or missing coordinate(s) at '
                                         'line {0:d}.'.format(i+1))
            if onlycoords:
                acount += 1
                i += 1
                continue
            
            serials[acount] = line[6:11]
            altlocs[acount] = alt 
            atomnames[acount] = atomname
            resnames[acount] = resname
            chainids[acount] = chid
            resnums[acount] = line[22:26]#.split()[0])
            icodes[acount] = line[26]
            if isPDB:
                try:
                    occupancies[acount] = line[54:60]
                except:
                    LOGGER.warning('failed to parse occupancy at line {0:d}'
                                   .format(i))
                try:
                    bfactors[acount] = line[60:66]
                except:
                    LOGGER.warning('failed to parse beta-factor at line {0:d}'
                                   .format(i))
                hetero[acount] = startswith[0] == 'H'
                segnames[acount] = line[72:76]
                elements[acount] = line[76:78]
            else:
                try:
                    charges[acount] = line[54:62]
                except:
                    LOGGER.warning('failed to parse charge at line {0:d}'
                                   .format(i))
                try:
                    radii[acount] = line[62:69]
                except:
                    LOGGER.warning('failed to parse radius at line {0:d}'
                                   .format(i))
            acount += 1
            if n_atoms == 0 and acount >= alength:
                # if arrays are short extend them with zeros
                alength += asize
                coordinates = np.concatenate(
                    (coordinates, np.zeros((asize, 3), float)))
                atomnames = np.concatenate((atomnames,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['name'].dtype)))
                resnames = np.concatenate((resnames,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['resname'].dtype)))
                resnums = np.concatenate((resnums,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['resnum'].dtype)))
                chainids = np.concatenate((chainids,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['chain'].dtype)))
                hetero = np.concatenate((hetero,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['hetero'].dtype)))
                altlocs = np.concatenate((altlocs,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['altloc'].dtype)))
                icodes = np.concatenate((icodes,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['icode'].dtype)))
                serials = np.concatenate((serials,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['serial'].dtype)))
                if isPDB:
                    bfactors = np.concatenate((bfactors,
                        np.zeros(asize, ATOMIC_DATA_FIELDS['beta'].dtype)))
                    occupancies = np.concatenate((occupancies,
                        np.zeros(asize, ATOMIC_DATA_FIELDS['occupancy'].dtype)))
                    segnames = np.concatenate((segnames,
                        np.zeros(asize, ATOMIC_DATA_FIELDS['segment'].dtype)))
                    elements = np.concatenate((elements,
                        np.zeros(asize, ATOMIC_DATA_FIELDS['element'].dtype)))
                    if anisou is not None:
                        anisou = np.concatenate((anisou, np.zeros((asize, 6), 
                            ATOMIC_DATA_FIELDS['anisou'].dtype)))
                    if siguij is not None:
                        siguij = np.concatenate((siguij, np.zeros((asize, 6), 
                            ATOMIC_DATA_FIELDS['siguij'].dtype)))
                else:
                    charges = np.concatenate((charges,
                        np.zeros(asize, ATOMIC_DATA_FIELDS['charge'].dtype)))
                    radii = np.concatenate((radii,
                        np.zeros(asize, ATOMIC_DATA_FIELDS['radius'].dtype)))
        #elif startswith == 'END   ' or startswith == 'CONECT':
        #    i += 1
        #    break
        elif startswith == 'ENDMDL' or startswith[:3] == 'END':
            if acount == 0:
                # If there is no atom record between ENDMDL & END skip to next
                i += 1
                continue
            if model is not None:
                i += 1
                break
            diff = stop - i - 1
            if diff < acount:
                END = True
            if onlycoords:
                if acount < n_atoms:
                    LOGGER.warning('Discarding model {0:d}, which contains '
                                   '{1:d} fewer atoms than the first model '
                                   'does.'.format(nmodel+1, n_atoms-acount))
                else:
                    coordsets[nmodel] = coordinates
                    nmodel += 1
                acount = 0
                if not END:
                    coordinates = coordsets[nmodel]
            else:
                if acount != n_atoms > 0:
                    raise ValueError('PDB file and AtomGroup ag must have '
                                    'same number of atoms')
                # this is where to decide if more coordsets should be expected
                if END: 
                    atomgroup.setCoords(coordinates[:acount])
                else:
                    coordsets = np.zeros((diff/acount+1, acount, 3))
                    coordsets[0] = coordinates[:acount]
                    onlycoords = True
                if not only_subset:
                    atomnames = np.char.strip(atomnames[:acount])
                    resnames = np.char.strip(resnames[:acount])
                atomgroup.setNames(atomnames[:acount])
                atomgroup.setResnames(resnames[:acount])
                atomgroup.setResnums(resnums[:acount])
                atomgroup.setChids(chainids[:acount])
                atomgroup.setHeteros(hetero[:acount])
                atomgroup.setAltlocs(altlocs[:acount])
                atomgroup.setIcodes(np.char.strip(icodes[:acount]))
                atomgroup.setSerials(serials[:acount])
                if isPDB:
                    atomgroup.setBetas(bfactors[:acount])
                    atomgroup.setOccupancies(occupancies[:acount])
                    atomgroup.setSegnames(np.char.strip(segnames[:acount]))
                    atomgroup.setElements(np.char.strip(elements[:acount]))
                    if anisou is not None:
                        atomgroup.setAnisous(anisou[:acount] / 10000)
                    if siguij is not None:
                        atomgroup.setAnistds(siguij[:acount] / 10000)
                else:
                    atomgroup.setCharges(charges[:acount])
                    atomgroup.setRadii(radii[:acount])
                    
                nmodel += 1
                n_atoms = acount 
                acount = 0
                coordinates = np.zeros((n_atoms, 3), dtype=float)
                if altloc and altloc_torf:
                    _evalAltlocs(atomgroup, altloc, chainids, resnums, 
                                 resnames, atomnames)
                    altloc = defaultdict(list)
                if END:
                    break
        elif isPDB and startswith == 'ANISOU':
            if anisou is None:
                anisou = True
                anisou = np.zeros((alength, 6), 
                    dtype=ATOMIC_DATA_FIELDS['anisou'].dtype)
            try:
                index = acount - 1
                anisou[index, 0] = line[28:35]
                anisou[index, 1] = line[35:42]
                anisou[index, 2] = line[43:49]
                anisou[index, 3] = line[49:56]
                anisou[index, 4] = line[56:63]
                anisou[index, 5] = line[63:70]
            except:
                LOGGER.warning('failed to parse anisotropic temperature '
                    'factors at line {0:d}'.format(i))
        elif isPDB and startswith =='SIGUIJ':
            if siguij is None:
                siguij = np.zeros((alength, 6), 
                    dtype=ATOMIC_DATA_FIELDS['siguij'].dtype)
            try:
                index = acount - 1
                siguij[index, 0] = line[28:35]
                siguij[index, 1] = line[35:42]
                siguij[index, 2] = line[43:49]
                siguij[index, 3] = line[49:56]
                siguij[index, 4] = line[56:63]
                siguij[index, 5] = line[63:70]
            except:
                LOGGER.warning('failed to parse standard deviations of '
                    'anisotropic temperature factors at line {0:d}'.format(i))
        elif startswith =='SIGATM':
            pass
        i += 1
    if onlycoords:
        if acount == atomgroup.numAtoms():
            coordsets[nmodel] = coordinates
            nmodel += 1
        if nmodel == coordsets.shape[0]:
            atomgroup.setCoords(coordsets)
        else:
            atomgroup.setCoords(coordsets[:nmodel])
    elif not END:
        # this means last line wast an ATOM line, so atomgroup is not decorated
        atomgroup.setCoords(coordinates[:acount])
        if not only_subset:
            atomnames = np.char.strip(atomnames[:acount])
            resnames = np.char.strip(resnames[:acount])
        atomgroup.setNames(atomnames[:acount])
        atomgroup.setResnames(resnames[:acount])
        atomgroup.setResnums(resnums[:acount])
        atomgroup.setChids(chainids[:acount])
        atomgroup.setHeteros(hetero[:acount])
        atomgroup.setAltlocs(altlocs[:acount])
        atomgroup.setIcodes(np.char.strip(icodes[:acount]))
        atomgroup.setSerials(serials[:acount])
        if isPDB:
            if anisou is not None:
                atomgroup.setAnisous(anisou[:acount] / 10000)
            if siguij is not None:
                atomgroup.setAnistds(siguij[:acount] / 10000)
            atomgroup.setSegnames(np.char.strip(segnames[:acount]))
            atomgroup.setElements(np.char.strip(elements[:acount]))
            atomgroup.setBetas(bfactors[:acount])
            atomgroup.setOccupancies(occupancies[:acount])
        else:
            atomgroup.setCharges(charges[:acount])
            atomgroup.setRadii(radii[:acount])
            

    if altloc and altloc_torf:
        _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames)
                
    return atomgroup

def _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames):
    altloc_keys = altloc.keys()
    altloc_keys.sort()
    indices = {}
    for key in altloc_keys:
        xyz = atomgroup.getCoords()
        success = 0
        lines = altloc[key]
        for line, i in lines:
            aan = line[12:16].strip()
            arn = line[17:21].strip()
            ach = line[21]
            ari = int(line[22:26].split()[0])
            rn, ids, ans = indices.get((ach, ari), (None, None, None))
            if ids is None:
                ids = indices.get(ach, None)
                if ids is None:
                    ids = (chainids == ach).nonzero()[0]
                    indices[ach] = ids
                ids = ids[resnums[ids] == ari]
                if len(ids) == 0:
                    LOGGER.warning("failed to parse alternate location '{0:s}'" 
                                   " at line {1:d}, residue does not exist as "
                                   "altloc 'A'".format(key, i+1))
                    continue
                rn = resnames[ids[0]]
                ans = atomnames[ids]
                indices[(ach, ari)] = (rn, ids, ans)
            if rn != arn:
                LOGGER.warning("failed to parse alternate location '{0:s}' at "
                               "line {1:d}, residue names do not match " 
                               "(expected '{2:s}', parsed '{3:s}')"
                               .format(key, i+1, rn, arn))
                continue
            index = ids[(ans == aan).nonzero()[0]]
            if len(index) != 1:
                LOGGER.warning("failed to parse alternate location '{0:s}' at "
                               "line {1:d}, could not identify matching atom "
                               "('{2:s}' not found in the residue)"
                               .format(key, i+1, aan))
                continue
            try:
                xyz[index[0], 0] = float(line[30:38])
                xyz[index[0], 1] = float(line[38:46])
                xyz[index[0], 2] = float(line[46:54])
            except:
                LOGGER.warning('failed to parse alternate location {0:s} at '
                               'line {1:d}, could not read coordinates'
                               .format(key, i+1))
                continue
            success += 1
        LOGGER.info('{0:d} out of {1:d} alternate location {2:s} lines were '
                    'parsed successfully.'.format(success, len(lines), key))
        if success > 0:
            LOGGER.info('Alternate location {0:s} is appended as a coordinate '
                        'set to the atom group.'
                        .format(key, atomgroup.getTitle()))
            atomgroup.addCoordset(xyz)

class Chemical(object):
    
    """A data structure for storing information on chemical components 
    (or heterogens) in PDB structures.
    
    .. versionadded:: 0.9
    
    A :class:`Chemical` instance has the following attributes:
        
    =========== ===== =========================================================
    Attribute   Type  Description (RECORD TYPE)
    =========== ===== =========================================================
    identifier  str   chemical component identifier (or residue name) (HET)
    name        str   chemical name (HETNAM)
    chain       str   chain identifier (HET)
    number      int   residue (or sequence) number (HET)
    icode       str   insertion code (HET)
    n_atoms     int   number of atoms present in the structure (HET)
    description str   description of the chemical component (HET)
    synonyms    list  synonyms (HETSYN)
    formula     str   chemical formula (FORMUL)
    pdbentry    str   PDB entry that chemical data is extracted from
    =========== ===== =========================================================
    
    Chemical class instances can be obtained as follows:
        
    >>> chemical = parsePDBHeader('1zz2', 'chemicals')[0]
    >>> chemical
    <Chemical: B11 (1ZZ2_A_362)>
    >>> print(chemical.name)
    N-[3-(4-FLUOROPHENOXY)PHENYL]-4-[(2-HYDROXYBENZYL) AMINO]PIPERIDINE-1-SULFONAMIDE
    >>> chemical.n_atoms
    33
    >>> len(chemical)
    33
    
    """
    
    __slots__ = ['identifier', 'name', 'chain', 'resnum', 'icode', 
                 'n_atoms', 'description', 'synonyms', 'formula', 'pdbentry']
    
    def __init__(self, identifier):
        
        #: residue name (or chemical component identifier)
        self.identifier = identifier
        #: chemical name
        self.name = None
        #: chain identifier
        self.chain = None
        #: residue (or sequence) number
        self.resnum = None
        #: insertion code
        self.icode = None
        #: number of atoms present in the structure
        self.n_atoms = None
        #: description of the chemical component
        self.description = None
        #: list of synonyms
        self.synonyms = None
        #: chemical formula
        self.formula = None
        #: PDB entry that chemical data is extracted from
        self.pdbentry = None
        
    def __str__(self):
        return self.identifier
    
    def __repr__(self):
        return '<Chemical: {0:s} ({1:s}_{2:s}_{3:d})>'.format(
                    self.identifier, self.pdbentry, self.chain, self.resnum)

    def __len__(self):
        return self.n_atoms

_PDB_DBREF = { 
    'GB': 'GenBank',
    'PDB': 'ProteinDataBank',
    'UNP': 'UniProt',
    'NORINE': 'Norine',
    'UNIMES': 'UNIMES'
}

class Polymer(object):
    
    """A data structure for storing information on polymer components 
    (protein or nucleic) of PDB structures.
    
    .. versionadded:: 0.9
    
    A :class:`Polymer` instance has the following attributes:
        
    ============= ====== ======================================================
    Attribute     Type   Description (RECORD TYPE)
    ============= ====== ======================================================
    identifier    str    chain identifier
    name          str    name of the polymer (macromolecule) (COMPND)
    fragment      str    specifies a domain or region of the molecule (COMPND)
    synonyms      list   list of synonyms for the polymer (COMPND)
    ec            list   list of associated Enzyme Commission numbers (COMPND)
    engineered    bool   indicates that the polymer was produced using 
                         recombinant technology or by purely chemical synthesis
                         (COMPND)
    mutation      bool   indicates presence of a mutation (COMPND)
    comments      str    additional comments
    sequence      str    polymer chain sequence (SEQRES)
    sqfirst       tuple  (resnum, icode) of the *first* residue in structure
    sqlast        tuple  (resnum, icode) of the *last* residue in structure
    dbabbr        str    reference sequence database abbreviation (DBREF[1|2])
    dbname        str    reference sequence database name (DBREF[1|2])
    dbidentifier  str    sequence database identification code (DBREF[1|2])
    dbaccession   str    sequence database accession code (DBREF[1|2])
    dbfirst       tuple  (resnum, icode) of the *first* residue in database
    dblast        tuple  (resnum, icode) of the *last* residue in database
    different     list   | differences from database sequence (SEQADV)
                         | when different residues are present, they will be 
                           each will be represented as: ``(residueName, 
                           residueNumberInsertionCode, dbResidueName, 
                           dbResidueNumber, comment)``
    modified      list   | modified residues (SEQMOD)
                         | when modified residues are present, each will be 
                           represented as: ``(residueName, 
                           residueNumberInsertionCode, standardName, comment)``
    pdbentry      str    PDB entry that polymer data is extracted from
    ============= ====== ======================================================
    
    Polymer class instances can be obtained as follows:
    
    >>> polymer = parsePDBHeader('2k39', 'polymers')[0]
    >>> polymer
    <Polymer: UBIQUITIN (2K39_A)>
    >>> print(polymer.pdbentry)
    2K39
    >>> print(polymer.identifier)
    A
    >>> print(polymer.name)
    UBIQUITIN
    >>> print(polymer.sequence)
    MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
    >>> len(polymer.sequence)
    76
    >>> len(polymer)
    76
    >>> print(polymer.dbname)
    UniProt
    >>> print(polymer.dbaccession)
    P62972
    >>> print(polymer.dbidentifier)
    UBIQ_XENLA
    
    """
    
    __slots__ = ['identifier', 'name', 'fragment', 'synonyms', 'ec', 
                 'engineered', 'mutation', 'comments', 'sequence', 'pdbentry', 
                 'dbabbr', 'dbname', 'dbidentifier', 'dbaccession', 
                 'modified', 'different',
                 'sqfirst', 'sqlast', 'dbfirst', 'dblast']
    
    def __init__(self, identifier):
        
        #: chain identifier
        self.identifier = identifier
        #: name of the polymer (macromolecule)
        self.name = ''
        #: specifies a domain or region of the molecule
        self.fragment = None
        #: list of synonyms for the molecule
        self.synonyms = None
        #: list of associated Enzyme Commission numbers
        self.ec = None
        self.engineered = None
        """indicates that the molecule was produced using recombinant 
        technology or by purely chemical synthesis"""
        #: indicates presence of a mutation
        self.mutation = None
        #: additional comments
        self.comments = None
        #: polymer chain sequence
        self.sequence = ''
        #: reference sequence database abbreviation 
        self.dbabbr = None
        #: reference sequence database name 
        self.dbname = None
        #: sequence database identification code
        self.dbidentifier = None
        #: sequence database accession code 
        self.dbaccession = None
        #: ``(resnum, icode)`` of the *first* residue in structure
        self.sqfirst = None
        #: ``(resnum, icode)`` of the *last* residue in structure
        self.sqlast = None
        #: ``(resnum, icode)`` of the *first* residue in database
        self.dbfirst = None
        #: ``(resnum, icode)`` of the *last* residue in database
        self.dblast = None
        #: modified residues
        self.modified = None
        #: differences from database reference sequence
        self.different = None
        #: PDB entry that polymer data is extracted from        
        self.pdbentry = None
        
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return '<Polymer: {0:s} ({1:s}_{2:s})>'.format(self.name, 
                                                self.pdbentry, self.identifier)

    def __len__(self): 
        return len(self.sequence)
    
_START_COORDINATE_SECTION = set(['ATOM  ', 'MODEL ', 'HETATM'])

def cleanString(string, nows=False):
    """*nows* is no white space."""
    
    if nows:
        return ''.join(string.strip().split())
    else:
        return ' '.join(string.strip().split())

def parsePDBHeader(pdb, *keys):
    """Return header data dictionary for *pdb*.  This function is equivalent to 
    ``parsePDB(pdb, header=True, model=0, meta=False)``, likewise *pdb* may be 
    an identifier or a filename.
    
    .. versionadded:: 0.9
    
    List of header records that are parsed. 
    
    ============ ================= ============================================
    Record type  Dictionary key(s)  Description 
    ============ ================= ============================================
    HEADER       | classification  | molecule classification 
                 | deposition_date | deposition date
                 | identifier      | PDB identifier
    TITLE        title             title for the experiment or analysis
    SPLIT        split             list of PDB entries that make up the whole
                                   structure when combined with this one
    COMPND       polymers          see :class:`Polymer`
    EXPDTA       experiment        information about the experiment
    NUMMDL       n_models          number of models
    MDLTYP       model_type        additional structural annotation
    AUTHOR       authors           list of contributors
    JRNL         reference         reference information dictionary:
                                     * *authors*: list of authors
                                     * *title*: title of the article
                                     * *editors*: list of editors
                                     * *issn*:
                                     * *reference*: journal, vol, issue, etc.
                                     * *publisher*: publisher information
                                     * *pmid*: pubmed identifier
                                     * *doi*: digital object identifier 
    DBREF[1|2]   polymers          see :class:`Polymer`
    SEQADV       polymers          see :class:`Polymer`
    SEQRES       polymers          see :class:`Polymer`
    MODRES       polymers          see :class:`Polymer`
    HELIX        polymers          see :class:`Polymer`
    SHEET        polymers          see :class:`Polymer`
    HET          chemicals         see :class:`Chemical`
    HETNAM       chemicals         see :class:`Chemical`
    HETSYN       chemicals         see :class:`Chemical`
    FORMUL       chemicals         see :class:`Chemical`
    REMARK 2     resolution        resolution of structures, when applicable
    REMARK 4     version           PDB file version
    REMARK 350   biomoltrans       biomolecular transformation lines 
                                   (unprocessed)
    ============ ================= ============================================
    
    Header records that are not parsed are: OBSLTE, CAVEAT, SOURCE, KEYWDS, 
    REVDAT, SPRSDE, SSBOND, LINK, CISPEP, CRYST1, ORIGX1, ORIGX2, ORIGX3, 
    MTRIX1, MTRIX2, MTRIX3, and REMARK X not mentioned above.
    
    
    Usage examples:
        
    >>> 
    
    """
    
    if not os.path.isfile(pdb):
        if len(pdb) == 4 and pdb.isalnum():
            filename = fetchPDB(pdb)
            if filename is None:
                raise IOError('PDB file for {0:s} could not be downloaded.'
                              .format(pdb))
            pdb = filename
        else:
            raise IOError('{0:s} is not a valid filename or a valid PDB '
                          'identifier.'.format(pdb))
    pdb = openFile(pdb)
    header, _ = _getHeaderDict(pdb, *keys)
    pdb.close()
    return header

def _getHeaderDict(stream, *keys):
    """Return header data in a dictionary.  *stream* may be a list of PDB lines
    or a stream."""
    
    lines = defaultdict(list)
    for loc, line in enumerate(stream):
        startswith = line[0:6]
        if startswith in _START_COORDINATE_SECTION:
            break
        lines[startswith].append((loc, line))
    for i, line in lines['REMARK']:
        lines[line[:10]].append((i, line))
    
    pdbid = _PDB_HEADER_MAP['identifier'](lines)
    lines['pdbid'] = pdbid
    if keys:
        keys = list(keys)
        for k, key in enumerate(keys):
            if key in _PDB_HEADER_MAP:
                value = _PDB_HEADER_MAP[key](lines)
                keys[k] = value
            else:
                raise KeyError('"{0:s}" is not a valid header data identifier'
                               .format(key))
            if key in ('chemicals', 'polymers'):
                for component in value:
                    component.pdbentry = pdbid
        if len(keys) == 1:
            return keys[0], loc
        else:
            return tuple(keys), loc
    else:
        header = {}
        keys = _PDB_HEADER_MAP.iterkeys()
        for key, func in _PDB_HEADER_MAP.iteritems():
            value = func(lines)
            if value is not None:
                header[key] = value        
        for chem in header.get('chemicals', []):
            chem.pdbentry = pdbid
        for poly in header.get('polymers', []):
            poly.pdbentry = pdbid
        return header, loc

def _getBiomoltrans(lines): 

    applyToChains = (' ')
    biomolecule = defaultdict(list)
    currentBiomolecule = '1'
    for i, line in lines['REMARK 350']:
        
        if line[13:18] == 'BIOMT':
            biomt = biomolecule[currentBiomolecule]
            if len(biomt) == 0:
                biomt.append(applyToChains)
            biomt.append(line[23:])
        elif line[11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
            applyToChains = line[41:].replace(' ', 
                                              '').strip().split(',')
        elif line[11:23] == 'BIOMOLECULE:': 
            currentBiomolecule = line.split()[-1]
    return dict(biomolecule)
        
        
def _getResolution(lines): 

    for i, line in lines['REMARK   2']:
        if 'RESOLUTION' in line:
            try:
                return float(line[23:30])
            except:
                return None

def _getHelix(lines):
    
    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    helix = {}
    for i, line in lines['HELIX ']:
        try:
            chid = line[19]
            #        helix class,      serial number,   identifier
            value = (int(line[38:40]), int(line[7:10]), line[11:14].strip())
        except:
            continue
        
        initICode = line[25]
        initResnum = int(line[21:25])
        if initICode != ' ':
            for icode in alphas[alphas.index(initICode):]:
                helix[(chid, initResnum, icode)] = value
            initResnum += 1    
        endICode = line[37]
        endResnum = int(line[33:37])
        if endICode != ' ':
            for icode in alphas[:alphas.index(endICode)+1]:
                helix[(chid, endResnum, icode)] = value
            endResnum -= 1    
        for resnum in range(initResnum, endResnum+1):
            helix[(chid, resnum, '')] = value
    return helix

def _getSheet(lines):
    
    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    sheet = {}
    for i, line in lines['SHEET ']:
        try:
            chid = line[21]
            value = (int(line[38:40]), int(line[7:10]), 
                     line[11:14].strip())
        except:
            continue

        initICode = line[26]
        initResnum = int(line[22:26])
        if initICode != ' ':
            for icode in alphas[alphas.index(initICode):]:
                sheet[(chid, initResnum, icode)] = value
            initResnum += 1    
        endICode = line[37]
        endResnum = int(line[33:37])
        if endICode != ' ':
            for icode in alphas[:alphas.index(endICode)+1]:
                sheet[(chid, endResnum, icode)] = value
            endResnum -= 1    
        for resnum in range(initResnum, endResnum+1):
            sheet[(chid, resnum, '')] = value
    return sheet


def _getReference(lines):
    """Return a reference of the PDB entry."""

    ref = {}
    title = ''
    authors = []
    editors = []
    reference = ''
    publisher = ''
    for i, line in lines['JRNL  ']:
        try:
            what = line.split(None, 2)[1]
        except:
            continue
        if what == 'AUTH':
            authors.extend(line[19:].strip().split(','))
        elif what == 'TITL':
            title += line[19:]
        elif what == 'EDIT':
            editors.extend(line[19:].strip().split(','))
        elif what == 'REF':
            reference += line[19:]
        elif what == 'PUBL':
            publisher += line[19:]
        elif what == 'REFN':
            ref['issn'] = line[19:].strip()
        elif what == 'PMID':
            ref['pmid'] = line[19:].strip()
        elif what == 'DOI':
            ref['doi'] = line[19:].strip()
    ref['authors'] = authors
    ref['title'] = cleanString(title)
    ref['editors'] = editors
    ref['reference'] = cleanString(reference)
    ref['publisher'] = cleanString(publisher)

    return ref


def _getPolymers(lines):
    """Return list of polymers (macromolecules)."""
    
    pdbid = lines['pdbid']
    polymers = dict()
    for i, line in lines['SEQRES']:
        ch = line[11]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        poly.sequence += ''.join(prody.compare.getSequence(line[19:].split()))
    for i, line in lines['DBREF ']:
        ch = line[12]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        poly.dbabbr = line[26:32].strip()
        poly.dbname = _PDB_DBREF.get(poly.dbabbr, 'Unknown')
        if poly.dbabbr == 'PDB':
            poly.dbaccession = poly.dbidentifier = pdbid
            if not (pdbid == line[33:37] == line[42:46] == line[7:11]):
                LOGGER.warning('wrong accession code and identifier for chain '
                               '{2:s} ({0:s}:{1:d})'.format(pdbid, i, ch))
        else:
            poly.dbaccession = line[33:41].strip()
            poly.dbidentifier = line[42:54].strip()
            if pdbid == poly.dbaccession or pdbid == poly.dbidentifier:
                LOGGER.warning('wrong database abbreviation for chain '
                               '{2:s} ({0:s}:{1:d})'.format(pdbid, i, ch))
        try:
            poly.sqfirst = (int(line[14:18]), line[18])
        except:
            LOGGER.warning('failed to parse first residue number for sequence '
                           '({0:s}:{1:d})'.format(pdbid, i))
        try:
            poly.sqlast = (int(line[20:24]), line[24])
        except:
            LOGGER.warning('failed to parse last residue number for sequence '
                           '({0:s}:{1:d})'.format(pdbid, i))
        try:
            poly.dbfirst = (int(line[56:60]), line[60])
        except:
            LOGGER.warning('failed to parse first residue number for database '
                           '({0:s}:{1:d})'.format(pdbid, i))
        try:
            poly.dblast = (int(line[62:67]), line[67])
        except:
            LOGGER.warning('failed to parse last residue number for database '
                           '({0:s}:{1:d})'.format(pdbid, i))
    for i, line in lines['DBREF1']:
        ch = line[12]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        poly.dbabbr = line[26:32].strip()
        poly.dbname = _PDB_DBREF.get(poly.dbabbr, 'Unknown')
        poly.dbidentifier = line[47:67].strip()
        try:
            poly.sqfirst = (int(line[14:18]), line[18])
        except:
            LOGGER.warning('failed to parse first residue number for sequence '
                           '({0:s}:{1:d})'.format(pdbid, i))
        try:
            poly.sqlast = (int(line[20:24]), line[24])
        except:
            LOGGER.warning('failed to parse last residue number for sequence '
                           '({0:s}:{1:d})'.format(pdbid, i))
    for i, line in lines['DBREF2']:
        ch = line[12]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        poly.dbaccession = line[18:40].strip()
        try:
            poly.dbfirst = (int(line[45:55]), '')
        except:
            LOGGER.warning('failed to parse first residue number for database '
                           '({0:s}:{1:d})'.format(pdbid, i))
        try:
            poly.dblast = (int(line[57:67]), '')
        except:
            LOGGER.warning('failed to parse last residue number for database '
                           '({0:s}:{1:d})'.format(pdbid, i))

    for i, line in lines['MODRES']:
        ch = line[16]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        if poly.modified is None:
            poly.modified = []
        poly.modified.append((line[12:15].strip(), line[18:22].strip() + 
                   line[22].strip(), line[24:27].strip(), line[29:70].strip()))
    for i, line in lines['SEQADV']:
        ch = line[16]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        if poly.different is None:
            poly.different = []
        dbabbr = line[24:28].strip()
        if poly.dbabbr != dbabbr:
            LOGGER.warning("reference database mismatch in SEQADV, expected "
                           "'{0:s}' parsed '{1:s}' ({2:s}:{3:d})"
                           .format(poly.dbabbr, dbabbr, pdbid, i))
            continue
        dbaccession = line[29:38].strip() 
        if poly.dbaccession != dbaccession:
            LOGGER.warning("database idcode mismatch in SEQADV, expected "
                           "'{0:s}' parsed '{1:s}' ({2:s}:{3:d})"
                           .format(poly.dbaccession, dbaccession, pdbid, i))
            continue
        poly.different.append((line[12:15].strip(), line[18:22].strip() + 
            line[22].strip(), line[39:42].strip(), line[43:48].strip(), 
            line[49:70].strip()))
    
    string = ''
    for i, line in lines['COMPND']:
        string += line[10:]
    string = cleanString(string)
    dict_ = {}
    for molecule in string.split('MOL_ID'):
        dict_.clear()
        for token in molecule.strip().split(';'):
            if not token:
                continue
            items = token.split(':', 1)
            if len(items) == 2:
                key, value = items 
            dict_[key.strip()] = value.strip()
        chains = dict_.pop('CHAIN', '').strip()
        if not chains:
            continue
        for ch in chains.split(','):
            ch = ch.strip()
            poly = polymers.get(ch, Polymer(ch))
            polymers[ch] = poly
            poly.name = dict_.get('MOLECULE', '').strip()
            poly.fragment = dict_.get('FRAGMENT', '').strip()
            poly.comments = dict_.get('OTHER_DETAILS', '').strip()
            poly.synonyms = [s.strip() 
                             for s in dict_.get('SYNONYMS', '').split(',')]
            poly.ec = [s.strip() 
                       for s in dict_.get('EC', '').split(',')]
            poly.engineered = dict_.get('ENGINEERED', '').strip() == 'YES'
            poly.mutation = dict_.get('MUTATION', '').strip() == 'YES'
        
    return polymers.values()    

def _getChemicals(lines):
    """Return list of chemical components (heterogens)."""
    
    chemicals = defaultdict(list)
    chem_names = defaultdict(str)
    chem_synonyms = defaultdict(str)
    chem_formulas = defaultdict(str)
    for i, line in lines['HET   ']:
        chem = Chemical(line[7:10].strip())
        chem.chain = line[12].strip()
        chem.resnum = int(line[13:17])
        chem.icode = line[17].strip()
        chem.n_atoms = int(line[20:25])
        chem.description = line[30:70].strip()
        chemicals[chem.identifier].append(chem)
    for i, line in lines['HETNAM']:
        chem = line[11:14].strip()
        chem_names[chem] += line[15:70].rstrip()
    for i, line in lines['HETSYN']:
        chem = line[11:14].strip()
        chem_synonyms[chem] += line[15:70].rstrip()
    for i, line in lines['FORMUL']:
        chem = line[12:15].strip()
        chem_formulas[chem] += line[18:70].rstrip()

    for chem, name in chem_names.iteritems():
        name = cleanString(name)
        for chem in chemicals[chem]:
            chem.name = name
    for chem, formula in chem_formulas.iteritems():
        formula = cleanString(formula)
        for chem in chemicals[chem]:
            chem.formula = formula
    for chem, synonyms in chem_synonyms.iteritems():
        synonyms = cleanString(synonyms)
        synonyms = synonyms.split(';')
        for chem in chemicals[chem]:
            chem.synonyms = synonyms
    
    alist = []
    for chem in chemicals.itervalues():
        for chem in chem:
            alist.append(chem)
    return alist 


def _getVersion(lines):
    
    for i, line in lines['REMARK   4']:
        if 'COMPLIES' in line:
            try:
                # Return a string, because floating makes 3.20, 3.2 or
                # may arise problems if wwPDB uses a version number like 3.30.1
                return line.split('V.')[1].split(',')[0].strip()
            except:
                return None

def _getNumModels(lines):

    # "NUMMDL", Integer, 11 - 14: Number of models.
    line = lines['NUMMDL']
    if line:
        i, line = line[0]
        try:
            header['n_models'] = int(line[10:14])
        except:
            return None

# Make sure that lambda functions defined below won't raise exceptions
_PDB_HEADER_MAP = {
    'helix': _getHelix,
    'sheet': _getSheet,
    'chemicals': _getChemicals,
    'polymers': _getPolymers,
    'reference': _getReference,
    'resolution': _getResolution,
    'biomoltrans': _getBiomoltrans,
    'version': _getVersion,
    'deposition_date': lambda lines: lines['HEADER'][0][1][50:59].strip() 
                                    if lines['HEADER'] else None,
    'classification': lambda lines: lines['HEADER'][0][1][10:50].strip() 
                                        if lines['HEADER'] else None,
    'identifier': lambda lines: lines['HEADER'][0][1][62:66].strip() 
                                    if lines['HEADER'] else None,
    'title': lambda lines: cleanString(
                ''.join([line[1][10:].rstrip() for line in lines['TITLE ']])
                ) if lines['TITLE '] else None,
    'experiment': lambda lines: cleanString(
                ''.join([line[1][10:].rstrip() for line in lines['EXPDTA']])
                ) if lines['EXPDTA'] else None,
    'authors': lambda lines: cleanString(
                ''.join([line[1][10:].rstrip() for line in lines['AUTHOR']]),
                True).split(',') if lines['AUTHOR'] else None,
    'split': lambda lines: (' '.join([line[1][11:].rstrip() 
                                      for line in lines['SPLIT ']])).split()
                             if lines['SPLIT '] else None,
    'model_type': lambda lines: cleanString(
                   ''.join([line[1][10:].rstrip() for line in lines['MDLTYP']])
                   ) if lines['MDLTYP'] else None,
    'n_models': _getNumModels,
}


def parsePSF(filename, title=None, ag=None):
    """Return an :class:`~prody.atomic.AtomGroup` instance storing data 
    parsed from X-PLOR format PSF file *filename*.  If *title* is not given, 
    *filename* will be set as the title of the :class:`AtomGroup` instance.  
    An :class:`AtomGroup` instance may be provided as *ag* argument.  When 
    provided, *ag* must have the same number of atoms in the same order as 
    the file.  Data from PSF file will be added to *ag*.  This may overwrite 
    present data if it overlaps with PSF file content.
    
    .. versionadded:: 0.8.1
    
    Note that this function does not evaluate bonds, angles, dihedrals, and
    impropers sections.
    
    """
    
    if ag is not None:
        if not isinstance(ag, prody.AtomGroup):
            raise TypeError('ag must be an AtomGroup instance') 
    
    psf = openFile(filename)
    line = psf.readline()
    i_line = 1
    while line:
        line = line.strip()
        if line.strip().endswith('!NATOM'):
            n_atoms = int(line.split('!')[0])
            break
        line = psf.readline()
        i_line += 1
    if title is None:
        title = os.path.splitext(os.path.split(filename)[1])[0]
    else:
        title = str(title)
    if ag is None:
        ag = prody.AtomGroup(title)
    else:
        if n_atoms != ag.numAtoms():
            raise ValueError('ag and PSF file must have same number of atoms')
        
    serials = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['serial'].dtype)
    segnames = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['segment'].dtype)
    resnums = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['resnum'].dtype)
    resnames = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['resname'].dtype)
    atomnames = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['name'].dtype)
    atomtypes = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['type'].dtype)
    charges = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['charge'].dtype)
    masses = np.zeros(n_atoms, ATOMIC_DATA_FIELDS['mass'].dtype)
    
    lines = psf.readlines(71 * (n_atoms + 5))
    index = True
    split = False
    if len(lines) < n_atoms:
        raise IOError('number of lines in PSF is less than the number of '
                      'atoms')
    
    for i, line in enumerate(lines):
        if i == n_atoms:
            break
        i_line += 1
        if index:
            try:            
                serials[i] = int(line[:8])
                segnames[i] = line[9:13].strip()
                resnums[i] = int(line[14:19])
                resnames[i] = line[19:23].strip()
                atomnames[i] = line[24:28].strip()
                atomtypes[i] = line[29:35].strip()
                charges[i] = float(line[35:44])
                masses[i] = float(line[50:60])
            except:
                LOGGER.warning('line {0:d} in {1:s} is not formatted as '
                               'expected, trying alternate method for the '
                               'rest of the file.'.format(i_line, filename))
                split = True
                index = False
        if split:
            try:
                items = line.strip().split()
                serials[i] = int(items[0])
                segnames[i] = items[1]
                resnums[i] = int(items[2])
                resnames[i] = items[3]
                atomnames[i] = items[4]
                atomtypes[i] = items[5]
                charges[i] = float(items[6])
                masses[i] = float(items[7])
            except:
                IOError('line {0:d} in {1:s} could not be parsed. Please '
                        'report this error.'.format(i_line, filename))
    psf.close()
    ag.setSerials(serials)
    ag.setSegnames(segnames)
    ag.setResnums(resnums)
    ag.setResnames(resnames)
    ag.setNames(atomnames)
    ag.setTypes(atomtypes)
    ag.setCharges(charges)
    ag.setMasses(masses)
    return ag

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
        """Return parameters used in blast search.
        
        .. versionadded: 0.8.3"""
        
        return self._param
        
    def getHits(self, percent_identity=90., percent_coverage=70., chain=False):
        """Return a dictionary that contains hits.
        
        Returns a dictionary whose keys are PDB identifiers.  Each dictionary
        entry is also a dictionary that contains information on the blast hit
        and alignment. 
        
        .. versionadded:: 0.8.2
           *chain* argument is added to allow user retrieve individual
           chains as hits.

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
    
    This method uses :meth:`qdblast` function in :mod:`NCBIWWW` module of 
    Biopython.  It will use *blastp* program and search *pdb* database.
    Results are parsed using :meth:`NCBIXML.parse` and passed to a
    :class:`PDBBlastRecord`
    
    User can adjust *hitlist_size* (default is 250) and *expect* (default is 
    1e-10) parameter values.
    
    User can also set the time interval for checking the results using
    *sleep* keyword argument.  Default value for *sleep* is 2 seconds.
    Finally, user can set a *timeout* value after which the function 
    will return ``None`` if the results are not ready.  Default for
    *timeout* is 30 seconds. 
    
    .. versionchanged:: 0.8.3
       *sequence* argument must be a string.  Biopython objects are not 
       accepted.  *sleep* and *timeout* arguments are added.

    """
    
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

_writePDBdoc = """
    :arg atoms: Atomic data container.
    :type atoms: :class:`~prody.atomic.Atomic` 
    
    :arg model: Model index or list of model indices.
    :type model: int, list
        
    If *models* is ``None``, all coordinate sets will be written. Model 
    indices start from 1.
    
    *atoms* instance must at least contain coordinates and atom names data.
    
    """

def writePDBStream(stream, atoms, model=None):
    """Write *atoms* in PDB format to a *stream*.
    
    :arg stream: anything that implements the method write() 
        (e.g. file, buffer, stdout)
    
    """
    
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms does not have a valid type')
    if isinstance(atoms, prody.Atom):
        atoms = prody.Selection(atoms.getAtomGroup(), [atoms.getIndex()], 
                                atoms.getACSI(), 
                                'index ' + str(atoms.getIndex()))

    if model is None:
        model = np.arange(atoms.numCoordsets(), dtype=np.int)
    elif isinstance(model, int):
        model = np.array([model], np.int) -1
    elif isinstance(model, list):
        model = np.array(model, np.int) -1
    else:
        raise TypeError('model must be an integer or a list of integers')
    if model.min() < 0 or model.max() >= atoms.numCoordsets():
        raise ValueError('model index or indices is not valid')
        
    n_atoms = atoms.numAtoms()
    atomnames = atoms.getNames()
    if atomnames is None:
        raise RuntimeError('atom names are not set')
    for i, an in enumerate(atomnames):
        lenan = len(an)
        if lenan < 4:
            atomnames[i] = ' ' + an
        elif lenan > 4:
            atomnames[i] = an[:4]
    resnames = atoms._getResnames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms
    resnums = atoms._getResnums()
    if resnums is None:
        resnums = np.ones(n_atoms, int)
    chainids = atoms._getChids()
    if chainids is None: 
        chainids = np.zeros(n_atoms, '|S1')
    occupancies = atoms._getOccupancies()
    if occupancies is None:
        occupancies = np.zeros(n_atoms, float)
    bfactors = atoms._getBetas()
    if bfactors is None:
        bfactors = np.zeros(n_atoms, float)
    icodes = atoms._getIcodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, '|S1')
    hetero = ['ATOM'] * n_atoms 
    heteroflags = atoms._getHeteros()
    if heteroflags is not None:
        hetero = np.array(hetero, '|S6')
        hetero[heteroflags] = 'HETATM'
    elements = atoms._getElements()
    if elements is None:
        elements = np.zeros(n_atoms, '|S1')
    altlocs = atoms._getAltlocs()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, '|S1')
    segments = atoms._getSegnames()
    if segments is None:
        segments = np.zeros(n_atoms, '|S6')
    
    acsi = atoms.getACSI()
    multi = False
    if len(model) > 1:
        multi = True
    format = ('{0:6s}{1:5d} {2:4s}{3:1s}' +
              '{4:4s}{5:1s}{6:4d}{7:1s}   ' + 
              '{8:8.3f}{9:8.3f}{10:8.3f}' +
              '{11:6.2f}{12:6.2f}      ' +
              '{13:4s}{14:2s}\n').format
    write = stream.write
    for m in model:
        if multi:
            stream.write('MODEL{0:9d}\n'.format(m+1))
        atoms.setACSI(m)
        coords = atoms._getCoords()
        for i, xyz in enumerate(coords):
            write(format(hetero[i], i+1, atomnames[i], altlocs[i], 
                         resnames[i], chainids[i], int(resnums[i]), 
                         icodes[i], 
                         xyz[0], xyz[1], xyz[2], 
                         occupancies[i], bfactors[i],  
                         segments[i], elements[i].rjust(2)))
        if multi:
            write('ENDMDL\n')
            altlocs = np.zeros(n_atoms, '|S1')
    atoms.setACSI(acsi)

writePDBStream.__doc__ += _writePDBdoc

def writePDB(filename, atoms, model=None):
    """Write *atoms* in PDB format to a file with name *filename*.  Returns 
    *filename* upon success.  If *filename* ends with :file:`.gz`, a compressed
    file will be written."""
    
    out = openFile(filename, 'w')
    writePDBStream(out, atoms, model)
    out.close()
    return filename

writePDB.__doc__ += _writePDBdoc

def writePQR(filename, atoms):
    """Write *atoms* in PQR format to a file with name *filename*.  Only 
    current coordinate set is written.  Returns *filename* upon success.  If 
    *filename* ends with :file:`.gz`, a compressed file will be written.."""
    
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms does not have a valid type')
    if isinstance(atoms, prody.Atom):
        atoms = prody.Selection(atoms.getAtomGroup(), [atoms.getIndex()], 
                                atoms.getACSI(), 
                                'index ' + str(atoms.getIndex()))
    stream = openFile(filename, 'w')
    n_atoms = atoms.numAtoms()
    atomnames = atoms.getNames()
    if atomnames is None:
        raise RuntimeError('atom names are not set')
    for i, an in enumerate(atomnames):
        lenan = len(an)
        if lenan < 4:
            atomnames[i] = ' ' + an
        elif lenan > 4:
            atomnames[i] = an[:4]
    resnames = atoms._getResnames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms
    resnums = atoms._getResnums()
    if resnums is None:
        resnums = np.ones(n_atoms, int)
    chainids = atoms._getChids()
    if chainids is None: 
        chainids = np.zeros(n_atoms, '|S1')
    charges = atoms._getCharges()
    if charges is None:
        charges = np.zeros(n_atoms, float)
    radii = atoms._getRadii()
    if radii is None:
        radii = np.zeros(n_atoms, float)
    icodes = atoms._getIcodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, '|S1')
    hetero = ['ATOM'] * n_atoms 
    heteroflags = atoms._getHeteros()
    if heteroflags is not None:
        hetero = np.array(hetero, '|S6')
        hetero[heteroflags] = 'HETATM'
    altlocs = atoms._getAltlocs()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, '|S1')
    
    format = ('{0:6s}{1:5d} {2:4s}{3:1s}' +
              '{4:4s}{5:1s}{6:4d}{7:1s}   ' + 
              '{8:8.3f}{9:8.3f}{10:8.3f}' +
              '{11:8.4f}{12:7.4f}\n').format
    coords = atoms._getCoords()
    write = stream.write
    for i, xyz in enumerate(coords):
        write(format(hetero[i], i+1, atomnames[i], altlocs[i], 
                     resnames[i], chainids[i], int(resnums[i]), 
                     icodes[i], xyz[0], xyz[1], xyz[2], charges[i], radii[i]))
    write('TER\nEND')
    stream.close()
    return filename

mapHelix = {
     1: 'H', # 4-turn helix (alpha helix)
     2: '',  # other helix, Right-handed omega
     3: 'I', # 5-turn helix (pi helix)
     4: '',  # other helix, Right-handed gamma
     5: 'G', # 3-turn helix (3-10 helix)
     6: '',  # Left-handed alpha
     7: '',  # Left-handed omega
     8: '',  # Left-handed gamma
     9: '',  # 2 - 7 ribbon/helix
    10: '',  # Polyproline
}

def assignSecondaryStructure(header, atoms, coil=False):
    """Deprecated, use :func:`assignSecstr`."""
    
    prody.deprecate('assignSecondaryStructure', 'assignSecstr')
    return assignSecstr(header, atoms, coil)

def assignSecstr(header, atoms, coil=False):
    """Assign secondary structure to *atoms* from *header* dictionary.

    *header* must be a dictionary parsed using the :func:`parsePDB`.
    *atoms* may be an instance of :class:`~prody.atomic.AtomGroup`, 
    :class:`~prody.atomic.Selection`, :class:`~prody.atomic.Chain` or 
    :class:`~prody.atomic.Residue`. 

    The Dictionary of Protein Secondary Structure, in short DSSP, type 
    single letter codes assignments are used:     
    
      * **G** = 3-turn helix (310 helix). Min length 3 residues.
      * **H** = 4-turn helix (alpha helix). Min length 4 residues.
      * **I** = 5-turn helix (pi helix). Min length 5 residues.
      * **T** = hydrogen bonded turn (3, 4 or 5 turn)
      * **E** = extended strand in parallel and/or anti-parallel 
        beta-sheet conformation. Min length 2 residues.
      * **B** = residue in isolated beta-bridge (single pair beta-sheet 
        hydrogen bond formation)
      * **S** = bend (the only non-hydrogen-bond based assignment).
      * **C** = residues not in one of above conformations.
    
    
    See http://en.wikipedia.org/wiki/Protein_secondary_structure#The_DSSP_code
    for more details. 
    
    Following PDB helix classes are omitted: 
    
      * Right-handed omega (2, class number)
      * Right-handed gamma (4)
      * Left-handed alpha (6)
      * Left-handed omega (7)
      * Left-handed gamma (8)
      * 2 - 7 ribbon/helix (9)
      * Polyproline (10)
      
    .. versionchanged:: 0.7
       Secondary structures are assigned to all atoms in a residue.  Amino acid
       residues without any secondary structure assignments in the header 
       section will be assigned coil (C) conformation.  This can be prevented
       by passing ``coil=False`` argument.
       
    .. versionchanged:: 0.8
       Default value for *coil* argument is set to ``False``.
       
    .. versionchanged:: 0.8.2
       Raises :class:`ValueError` when *header* does not contain 
       secondary structure data.
    
    """
    
    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')
    helix = header.get('helix', {})
    sheet = header.get('sheet', {})
    if len(helix) == 0 and len(sheet) == 0:
        raise ValueError('header does not contain secondary structure data')

    ssa = atoms.getSecstrs()
    if ssa is None:
        if isinstance(atoms, prody.AtomGroup):
            ag = atoms
        else:
            ag = atoms.getAtomGroup()
        ag.setSecstrs(np.zeros(ag.numAtoms(), 
                            ATOMIC_DATA_FIELDS['secondary'].dtype))
    atoms.select('protein').setSecstrs('C')
    hierview = atoms.getHierView()
    count = 0
    for key, value in helix.iteritems():
        res = hierview.getResidue(*key)
        if res is None:
            continue
        res.setSecstrs(mapHelix[value[0]])
        count += 1
    for key, res in sheet.iteritems():
        res = hierview.getResidue(*key)
        if res is None:
            continue
        res.setSecstrs('E')
        count += 1
    LOGGER.info('Secondary structures were assigned to {0:d} residues.'
                .format(count))
    return atoms        

def fetchLigandData(cci, save=False, folder='.'):
    """Fetch ligand data from `Ligand Expo <http://ligand-expo.rcsb.org/>`_.

    .. versionadded:: 0.8.1
    
    .. versionchanged:: 0.8.2
       URL of the XML file is returned in the dictionary with key ``url``.
    
    *cci* may be 3-letter chemical component identifier or a valid XML 
    filename.  If ``save=True`` is passed, XML file will be saved in the 
    specified *folder*. 
    
    This function is compatible with PDBx/PDBML v 4.0.
    
    Returns a dictionary that contains ligand data.  Ligand atom data with 
    *model* and *ideal* coordinate sets are also stored in this dictionary.
    Note that this dictionary will contain data that is present in the XML
    file and all Ligand Expo XML files do not contain every possible data
    field.  So, it may be better if you use :meth:`dict.get` instead of
    indexing the dictionary, e.g. to retrieve formula weight (or relative
    molar mass) of the chemical component use ``data.get('formula_weight')``
    instead of ``data['formula_weight']`` to avoid exceptions when this
    data field is not found in the XML file.
    
    Following example downloads data for ligand STI (a.k.a. Gleevec and
    Imatinib) and calculates RMSD between model (X-ray structure 1IEP) and 
    ideal (energy minimized) coordinate sets::
    
      ligand_data = fetchLigandData('STI')
      print ligand_data['model_coordinates_db_code'] 
      # 1IEP
      ligand_model = ligand_data['model'] 
      ligand_ideal = ligand_data['ideal']
      transformation = superpose(ligand_ideal.noh, ligand_model.noh)
      print( calcRMSD(ligand_ideal.noh, ligand_model.noh).round(2) )
      # 2.27
    """

    
    if not isinstance(cci, str):
        raise TypeError('cci must be a string')
    if os.path.isfile(cci):
        inp = openFile(cci)
        xml = inp.read()
        inp.close()
    elif len(cci) > 4 or not cci.isalnum(): 
        raise ValueError('cci must be 3-letters long and alphanumeric or '
                         'a valid filename')
    else:
        #'http://www.pdb.org/pdb/files/ligand/{0:s}.xml'
        url = ('http://ligand-expo.rcsb.org/reports/{0[0]:s}/{0:s}/{0:s}.xml'
               .format(cci.upper()))
        import urllib2
        try:
            inp = urllib2.urlopen(url)
        except urllib2.HTTPError:
            raise IOError('XML file for ligand {0:s} is not found'.format(cci))
        else:
            xml = inp.read()
            inp.close()
        if save:
            out = openFile(cci+'.xml', mode='w', folder=folder)
            out.write(xml)
            out.close()

    import xml.etree.cElementTree as ET

    root = ET.XML(xml)
    if root.get('{http://www.w3.org/2001/XMLSchema-instance}schemaLocation') \
        != 'http://pdbml.pdb.org/schema/pdbx-v40.xsd pdbx-v40.xsd':
        LOGGER.warning('XML does not seem to be in PDBx/PDBML v 4.0 format, '
                       'resulting dictionary may not have all possible data')
    ns = root.tag[:root.tag.rfind('}')+1]
    len_ns = len(ns)
    dict_ = {'url': url}

    for child in list(root.find(ns + 'chem_compCategory')[0]):
        tag = child.tag[len_ns:]
        if tag.startswith('pdbx_'):
            tag = tag[5:]
        dict_[tag] = child.text
    dict_['formula_weight'] = float(dict_.get('formula_weight')) 

    identifiers_and_descriptors = []
    results = root.find(ns + 'pdbx_chem_comp_identifierCategory')
    if results:
        identifiers_and_descriptors.extend(results)
    results = root.find(ns + 'pdbx_chem_comp_descriptorCategory')
    if results:
        identifiers_and_descriptors.extend(results)
    for child in identifiers_and_descriptors:
        program = child.get('program').replace(' ', '_')
        type_ = child.get('type').replace(' ', '_')
        dict_[program + '_' + type_] = child[0].text
        dict_[program + '_version'] = child.get('program_version')

    dict_['audits'] = [(audit.get('action_type'), audit.get('date'))
        for audit in list(root.find(ns + 'pdbx_chem_comp_auditCategory'))]

    atoms = list(root.find(ns + 'chem_comp_atomCategory'))
    n_atoms = len(atoms)
    ideal_coords = np.zeros((n_atoms, 3))
    model_coords = np.zeros((n_atoms, 3))
    
    atomnames = np.zeros(n_atoms, dtype=ATOMIC_DATA_FIELDS['name'].dtype) 
    elements = np.zeros(n_atoms, dtype=ATOMIC_DATA_FIELDS['element'].dtype)
    resnames = np.zeros(n_atoms, dtype=ATOMIC_DATA_FIELDS['resname'].dtype)
    charges = np.zeros(n_atoms, dtype=ATOMIC_DATA_FIELDS['charge'].dtype)
    
    resnums = np.ones(n_atoms, dtype=ATOMIC_DATA_FIELDS['charge'].dtype)
    
    alternate_atomnames = np.zeros(n_atoms, 
                                        dtype=ATOMIC_DATA_FIELDS['name'].dtype)
    leaving_atom_flags = np.zeros(n_atoms, np.bool)
    aromatic_flags = np.zeros(n_atoms, np.bool)
    stereo_configs = np.zeros(n_atoms, np.bool)
    ordinals = np.zeros(n_atoms, int)
    
    for i, atom in enumerate(atoms):
        data = dict([(child.tag[len_ns:], child.text) for child in list(atom)])

        atomnames[i] = data.get('pdbx_component_atom_id', 'X')
        elements[i] = data.get('type_symbol', 'X')
        resnames[i] = data.get('pdbx_component_comp_id', 'UNK')
        charges[i] = float(data.get('charge', 0))
        
        alternate_atomnames[i] = data.get('alt_atom_id', 'X')
        leaving_atom_flags[i] = data.get('pdbx_leaving_atom_flag') == 'Y'
        aromatic_flags[i] = data.get('pdbx_atomatic_flag') == 'Y'
        stereo_configs[i] = data.get('pdbx_stereo_config') == 'Y'
        ordinals[i] = int(data.get('pdbx_ordinal',0))
        
        model_coords[i, 0] = float(data.get('model_Cartn_x', 0))
        model_coords[i, 1] = float(data.get('model_Cartn_y', 0))
        model_coords[i, 2] = float(data.get('model_Cartn_z', 0))
        ideal_coords[i, 0] = float(data.get('pdbx_model_Cartn_x_ideal', 0))
        ideal_coords[i, 1] = float(data.get('pdbx_model_Cartn_y_ideal', 0))
        ideal_coords[i, 2] = float(data.get('pdbx_model_Cartn_z_ideal', 0))

    model = AtomGroup(cci + ' model')
    model.setCoords(model_coords)
    model.setNames(atomnames)
    model.setResnames(resnames)
    model.setResnums(resnums)
    model.setElements(elements)
    model.setCharges(charges)
    model.setData('leaving_atom_flags', leaving_atom_flags)
    model.setData('aromatic_flags', aromatic_flags)
    model.setData('stereo_configs', stereo_configs)
    model.setData('ordinals', ordinals)
    model.setData('alternate_atomnames', alternate_atomnames)
    dict_['model'] = model
    ideal = model.copy()
    ideal.setTitle(cci + ' ideal')
    ideal.setCoords(ideal_coords)
    dict_['ideal'] = ideal

    return dict_      

def applyBiomolecularTransformations(header, atoms, biomol=None):
    """Deprecated, use :func:`buildBiomolecules`."""
    
    prody.deprecate('applyBiomolecularTransformations', 'buildBiomolecules')
    return buildBiomolecules(header, atoms, biomol)

def buildBiomolecules(header, atoms, biomol=None):
    """Return *atoms* after applying biomolecular transformations from *header*
    dictionary.
    
    .. versionchanged:: 0.7.1
       Biomolecular transformations are applied to all coordinate sets in the
       molecule.
    
    .. versionchanged:: 0.8.2
       Raises :class:`ValueError` when *header* does not contain biomolecular 
       transformations.

    Some PDB files contain transformations for more than 1 biomolecules.  A 
    specific set of transformations can be choosen using *biomol* argument.
    Transformation sets are identified by numbers, e.g. ``"1"``, ``"2"``, ...
    
    If multiple biomolecular transformations are provided in the *header*
    dictionary, biomolecules will be returned as 
    :class:`~prody.atomic.AtomGroup` instances in a :func:`list`.  

    If the resulting biomolecule has more than 26 chains, the molecular 
    assembly will be split into multiple :class:`~prody.atomic.AtomGroup`
    instances each containing at most 26 chains.  These 
    :class:`~prody.atomic.AtomGroup` instances will be returned in a tuple.
   
    Note that atoms in biomolecules are ordered according to chain identifiers.
   
    """

    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must be an Atomic instance')
    biomt = header.get('biomoltrans')
    if not isinstance(biomt, dict) or len(biomt) == 0:
        raise ValueError("header doesn't contain biomolecular transformations")
    
    if not isinstance(atoms, prody.AtomGroup):
        atoms = atoms.copy()
    biomols = []
    if biomol is None: 
        keys = biomt.keys()
    else:
        biomol = str(biomol)
        if biomol in biomt:
            keys = [biomol]
        else:
            LOGGER.warning('Transformations for biomolecule {0:s} was not '
                           'found in the header dictionary.'.format(biomol))
            return None

    c_max = 26
    keys.sort()
    for i in keys: 
        chids = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'*20)
        ags = []
        mt = biomt[i]
        # mt is a list, first item is list of chain identifiers
        # following items are lines corresponding to transformation
        # mt must have 3n + 1 lines
        if (len(mt) - 1) % 3 != 0:
            LOGGER.warning('Biomolecular transformations {0:s} were not '
                           'applied'.format(i))
            continue
        chids_used = []
        for times in range((len(mt) - 1)/ 3):
            rotation = np.zeros((3,3))
            translation = np.zeros(3)
            line = np.fromstring(mt[times*3+1], sep=' ')
            rotation[0,:] = line[:3]
            translation[0] = line[3]
            line = np.fromstring(mt[times*3+2], sep=' ')
            rotation[1,:] = line[:3]
            translation[1] = line[3]
            line = np.fromstring(mt[times*3+3], sep=' ')
            rotation[2,:] = line[:3]
            translation[2] = line[3]
            t = prody.Transformation(rotation, translation)
            
            for chid in mt[0]:
                chids_used.append(chid)
                newag = atoms.copy('chain ' + chid)
                if newag is None:
                    continue
                
                for acsi in range(newag.numCoordsets()):
                    newag.setACSI(acsi)
                    newag = t.apply(newag)
                newag.setACSI(0)
                ags.append(newag)
        if ags:
            # Handles the case when there is more atom groups than the number
            # of chain identifiers
            if len(chids_used) != len(set(chids_used)):
                for newag in ags:
                    newag.select('all').setChids(chids.pop(0))
            if len(ags) <= c_max:
                ags_ = [ags]
            else:
                ags_ = []
                for j in range(len(ags)/c_max + 1):
                    ags_.append(ags[j*c_max: (j+1)*c_max])
            parts = []
            for k, ags in enumerate(ags_):
                if not ags:
                    continue
                newag = ags.pop(0)
                while ags:
                    newag += ags.pop(0)
                if len(ags_) > 1:
                    newag.setTitle('{0:s} biomolecule {1:s} part {2:d}'
                                   .format(atoms.getTitle(), i, k+1))
                else:
                    newag.setTitle('{0:s} biomolecule {1:s}'
                                   .format(atoms.getTitle(), i))
                parts.append(newag)
            if len(parts) == 1:
                biomols.append(newag)
            else:
                biomols.append(tuple(parts))
    if biomols:
        if len(biomols) == 1:
            return biomols[0]
        else:
            return biomols
    else:
        return None

def execDSSP(pdb, outputname=None, outputdir=None, stderr=True):
    """Execute DSSP for given *pdb*.  *pdb* can be a PDB identifier or a PDB 
    file path.  If *pdb* is a compressed file, it will be decompressed using
    Python :mod:`gzip` library.  When no *outputname* is given, output name 
    will be :file:`pdb.dssp`.  :file:`.dssp` extension will be appended 
    automatically to *outputname*.  If :file:`outputdir` is given, DSSP 
    output and uncompressed PDB file will be written into this folder.
    Upon successful execution of :command:`dssp pdb > out` command, output
    filename is returned.  On Linux platforms, when *stderr* is false, 
    standard error messages are suppressed, i.e.
    ``dssp pdb > outputname 2> /dev/null``.
    
    For more information on DSSP see http://swift.cmbi.ru.nl/gv/dssp/.
    If you benefited from DSSP, please consider citing [WK83]_.
    
    .. versionadded:: 0.8
    
    .. versionchanged:: 0.9.2
       *stderr* keyword argument is added. 
    """
    
    dssp = which('dssp')
    if dssp is None:
        raise EnvironmentError('command not found: dssp executable is not '
                               'found in one of system paths')
    assert outputname is None or isinstance(outputname, str),\
        'outputname must be a string'
    assert outputdir is None or isinstance(outputdir, str),\
        'outputdir must be a string'
    if not os.path.isfile(pdb):
        pdb = fetchPDB(pdb, compressed=False)
    if pdb is None:
        raise ValueError('pdb is not a valid PDB identifier or filename')
    if os.path.splitext(pdb)[1] == '.gz':
        if outputdir is None:
            pdb = gunzip(pdb, os.path.splitext(pdb)[0])
        else:
            pdb = gunzip(pdb, os.path.join(outputdir, 
                                os.path.split(os.path.splitext(pdb)[0])[1]))
    if outputdir is None:
        outputdir = '.'
    if outputname is None:
        out = os.path.join(outputdir,
                        os.path.splitext(os.path.split(pdb)[1])[0] + '.dssp')
    else:
        out = os.path.join(outputdir, outputname + '.dssp')
        
    if not stderr and PLATFORM != 'Windows':
        status = os.system('{0:s} {1:s} > {2:s} 2> /dev/null'.format(
                            dssp, pdb, out))
    else:
        status = os.system('{0:s} {1:s} > {2:s}'.format(dssp, pdb, out))

    if status == 0:
        return out
    
def parseDSSP(dssp, ag, parseall=False):
    """Parse DSSP data from file *dssp* into :class:`~prody.atomic.AtomGroup`
    instance *ag*.  DSSP output file must be in the new format used from July 
    1995 and onwards.  When *dssp* file is parsed, following attributes are 
    added to *ag*:
        
    * *dssp_resnum*: DSSP's sequential residue number, starting at the first 
      residue actually in the data set and including chain breaks; this number 
      is used to refer to residues throughout.
    
    * *dssp_acc*: number of water molecules in contact with this residue \*10. 
      or residue water exposed surface in Angstrom^2.
      
    * *dssp_kappa*: virtual bond angle (bend angle) defined by the three Cα 
      atoms of residues I-2,I,I+2.  Used to define bend (structure code 'S').
        
    * *dssp_alpha*: virtual torsion angle (dihedral angle) defined by the four 
      Cα atoms of residues I-1,I,I+1,I+2.Used to define chirality (structure 
      code '+' or '-').
        
    * *dssp_phi* and *dssp_psi*: IUPAC peptide backbone torsion angles 

    The following attributes are parsed when ``parseall=True`` is passed: 

    * *dssp_bp1*, *dssp_bp2*, and *dssp_sheet_label*: residue number of first 
      and second bridge partner followed by one letter sheet label
      
    * *dssp_tco*: cosine of angle between C=O of residue I and C=O of residue 
      I-1.  For α-helices, TCO is near +1, for β-sheets TCO is near -1.  Not 
      used for structure definition.
      
    * *dssp_NH_O_1_index*, *dssp_NH_O_1_energy*, etc.: hydrogen bonds; e.g. 
      -3,-1.4 means: if this residue is residue i then N-H of I is h-bonded to 
      C=O of I-3 with an electrostatic H-bond energy of -1.4 kcal/mol.  There 
      are two columns for each type of H-bond, to allow for bifurcated H-bonds.
        
    See http://swift.cmbi.ru.nl/gv/dssp/DSSP_3.html for details.
    
    .. versionadded:: 0.8"""
    
    if not os.path.isfile(dssp):
        raise IOError('{0:s} is not a valid file path'.format(dssp))
    if not isinstance(ag, prody.AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')
        
    dssp = open(dssp)
    
    n_atoms = ag.numAtoms()
    NUMBER = np.zeros(n_atoms, int)
    SHEETLABEL = np.zeros(n_atoms, '|S1')
    ACC = np.zeros(n_atoms, float)
    KAPPA = np.zeros(n_atoms, float)
    ALPHA = np.zeros(n_atoms, float)
    PHI = np.zeros(n_atoms, float)
    PSI = np.zeros(n_atoms, float)

    if parseall:
        BP1 = np.zeros(n_atoms, int)
        BP2 = np.zeros(n_atoms, int)
        NH_O_1 = np.zeros(n_atoms, int)
        NH_O_1_nrg = np.zeros(n_atoms, float)
        O_HN_1 = np.zeros(n_atoms, int)
        O_HN_1_nrg = np.zeros(n_atoms, float)
        NH_O_2 = np.zeros(n_atoms, int)
        NH_O_2_nrg = np.zeros(n_atoms, float)
        O_HN_2 = np.zeros(n_atoms, int)
        O_HN_2_nrg = np.zeros(n_atoms, float)
        TCO = np.zeros(n_atoms, float)

    ag.setSecstrs(np.zeros(n_atoms, 
                              dtype=ATOMIC_DATA_FIELDS['secondary'].dtype))
    for line in dssp:
        if line.startswith('  #  RESIDUE'):
            break
    for line in dssp:
        if line[13] == '!':
            continue
        res = ag[(line[11], int(line[5:10]), line[10].strip())]
        if res is None:
            continue
        indices = res.getIndices()
        res.setSecstrs(line[16].strip())
        NUMBER[indices] = int(line[:5])
        SHEETLABEL[indices] = line[33].strip()
        ACC[indices] = int(line[35:38])
        KAPPA[indices] = float(line[91:97])
        ALPHA[indices] = float(line[97:103])
        PHI[indices] = float(line[103:109])
        PSI[indices] = float(line[109:115])

        if parseall:
            BP1[indices] = int(line[25:29])
            BP2[indices] = int(line[29:33])
            NH_O_1[indices] = int(line[38:45])
            NH_O_1_nrg[indices] = float(line[46:50]) 
            O_HN_1[indices] = int(line[50:56])
            O_HN_1_nrg[indices] = float(line[57:61])
            NH_O_2[indices] = int(line[61:67])
            NH_O_2_nrg[indices] = float(line[68:72])
            O_HN_2[indices] = int(line[72:78])
            O_HN_2_nrg[indices] = float(line[79:83])
            TCO[indices] = float(line[85:91])
    
    ag.setData('dssp_resnum', NUMBER)
    ag.setData('dssp_sheet_label', SHEETLABEL)
    ag.setData('dssp_acc', ACC)
    ag.setData('dssp_kappa', KAPPA)
    ag.setData('dssp_alpha', ALPHA)
    ag.setData('dssp_phi', PHI)
    ag.setData('dssp_psi', PSI)

    if parseall:
        ag.setData('dssp_bp1', BP1)
        ag.setData('dssp_bp2', BP2)
        ag.setData('dssp_NH_O_1_index', NH_O_1)
        ag.setData('dssp_NH_O_1_energy', NH_O_1_nrg)
        ag.setData('dssp_O_NH_1_index', O_HN_1)
        ag.setData('dssp_O_NH_1_energy', O_HN_1_nrg)    
        ag.setData('dssp_NH_O_2_index', NH_O_2)
        ag.setData('dssp_NH_O_2_energy', NH_O_2_nrg)
        ag.setData('dssp_O_NH_2_index', O_HN_2)
        ag.setData('dssp_O_NH_2_energy', O_HN_2_nrg)
        ag.setData('dssp_tco', TCO)
    return ag

def performDSSP(pdb, parseall=False, stderr=True):
    """Perform DSSP calculations and parse results.  DSSP data is returned 
    in an :class:`~prody.atomic.AtomGroup` instance.  See also :func:`execDSSP` 
    and :func:`parseDSSP`.
    
    .. versionadded:: 0.8
    
    .. versionchanged:: 0.9.2
       Added *stderr* argument, see :func:`execDSSP` for details."""
    
    pdb = fetchPDB(pdb, compressed=False)
    return parseDSSP(execDSSP(pdb, stderr=stderr), parsePDB(pdb), parseall)
    
def execSTRIDE(pdb, outputname=None, outputdir=None):
    """Execute STRIDE program for given *pdb*.  *pdb* can be an identifier or 
    a PDB file path.  If *pdb* is a compressed file, it will be decompressed 
    using Python :mod:`gzip` library.  When no *outputname* is given, output 
    name will be :file:`pdb.stride`.  :file:`.stride` extension will be 
    appended automatically to *outputname*.  If :file:`outputdir` is given, 
    STRIDE output and uncompressed PDB file will be written into this folder.
    Upon successful execution of :command:`stride pdb > out` command, output
    filename is returned. 
    
    For more information on STRIDE see http://webclu.bio.wzw.tum.de/stride/.
    If you benefited from STRIDE, please consider citing [DF95]_.
    
    .. versionadded:: 0.8"""
    
    stride = which('stride')
    if stride is None:
        raise EnvironmentError('command not found: stride executable is not '
                               'found in one of system paths')
    assert outputname is None or isinstance(outputname, str),\
        'outputname must be a string'
    assert outputdir is None or isinstance(outputdir, str),\
        'outputdir must be a string'
    if not os.path.isfile(pdb):
        pdb = fetchPDB(pdb, compressed=False)
    if pdb is None:
        raise ValueError('pdb is not a valid PDB identifier or filename')
    if os.path.splitext(pdb)[1] == '.gz':
        if outputdir is None:
            pdb = gunzip(pdb, os.path.splitext(pdb)[0])
        else:
            pdb = gunzip(pdb, os.path.join(outputdir, 
                                os.path.split(os.path.splitext(pdb)[0])[1]))
    if outputdir is None:
        outputdir = '.'
    if outputname is None:
        out = os.path.join(outputdir,
                        os.path.splitext(os.path.split(pdb)[1])[0] + '.stride')
    else:
        out = os.path.join(outputdir, outputname + '.stride')
        
    status = os.system('{0:s} {1:s} > {2:s}'.format(stride, pdb, out))
    if status == 0:
        return out
    
def parseSTRIDE(stride, ag):
    """Parse STRIDE output from file *stride* into 
    :class:`~prody.atomic.AtomGroup` instance *ag*.  STRIDE output file must 
    be in the new format used from July 1995 and onwards.  When *stride* file 
    is parsed, following attributes are added to *ag*:
        
    * *stride_resnum*: STRIDE's sequential residue number, starting at the 
      first residue actually in the data set.
    
    * *stride_phi*, *stride_psi*: peptide backbone torsion angles phi and psi
    
    * *stride_area*: residue solvent accessible area
    
    .. versionadded:: 0.8"""
    
    if not os.path.isfile(stride):
        raise IOError('{0:s} is not a valid file path'.format(stride))
    if not isinstance(ag, prody.AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')
        
    stride = open(stride)
    
    n_atoms = ag.numAtoms()
    NUMBER = np.zeros(n_atoms, int)
    AREA = np.zeros(n_atoms, float)
    PHI = np.zeros(n_atoms, float)
    PSI = np.zeros(n_atoms, float)

    ag.setSecstrs(np.zeros(n_atoms),
                     dtype=ATOMIC_DATA_FIELDS['secondary'].dtype)
    for line in stride:
        if not line.startswith('ASG '):
            continue
        res = ag[(line[9], int(line[10:15]), line[15].strip())]
        if res is None:
            continue
        indices = res.getIndices()
        res.setSecstrs(line[24].strip())
        NUMBER[indices] = int(line[16:20])
        PHI[indices] = float(line[42:49])
        PSI[indices] = float(line[52:59])
        AREA[indices] = float(line[64:69])
    ag.setData('stride_resnum', NUMBER)
    ag.setData('stride_phi', PHI)
    ag.setData('stride_psi', PSI)
    ag.setData('stride_area', AREA)
    return ag

def performSTRIDE(pdb):
    """Perform STRIDE calculations and parse results.  STRIDE data is 
    returned in an :class:`~prody.atomic.AtomGroup` instance.  See also 
    :func:`execSTRIDE` and :func:`parseSTRIDE`.
    
    .. versionadded:: 0.8"""
    
    pdb = fetchPDB(pdb, compressed=False)
    return parseSTRIDE(execSTRIDE(pdb), parsePDB(pdb))

def fetchPDBClusters():
    """Downloads PDB sequence clusters.  PDB sequence clusters are results of 
    the weekly clustering of protein chains in the PDB generated by blastclust. 
    They are available at FTP site: ftp://resources.rcsb.org/sequence/clusters/
    
    This function will download about 10 Mb of data and save it after 
    compressing in your home directory in :file:`.prody/pdbclusters`.
    Compressed files will be less than 4 Mb in size.  Cluster data can 
    be loaded using :func:`loadPDBClusters` function and be accessed 
    using :func:`getPDBCluster`.
    
    .. versionadded:: 0.8.2"""
    
    import urllib2
    PDB_CLUSTERS_PATH = os.path.join(prody.getPackagePath(), 'pdbclusters')
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

def loadPDBClusters(sqid=None):
    """Load previously fetched PDB sequence clusters from disk to memory.
    
    .. versionadded:: 0.8.2"""

    PDB_CLUSTERS_PATH = os.path.join(prody.getPackagePath(), 'pdbclusters')
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
    using :func:`loadPDBClusters`.
    
    .. versionadded:: 0.8.2"""

    assert isinstance(pdb, str) and len(pdb) == 4, \
        'pdb must be 4 char long string'
    assert isinstance(ch, str) and len(ch) == 1, \
        'ch must be a one char long string'
    assert isinstance(sqid, int), 'sqid must be an integer'
    PDB_CLUSTERS_PATH = os.path.join(prody.getPackagePath(), 'pdbclusters')
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

def showProtein(*atoms, **kwargs):
    """Show protein representation using :meth:`~mpl_toolkits.mplot3d.Axes3D`.
    This function is designed for generating a quick view of the contents of a 
    :class:`~prody.atomic.AtomGroup` or :class:`~prody.atomic.Selection`.
           
    .. versionadded:: 0.9
    
    Protein atoms matching ``"calpha"`` selection are displayed using solid 
    lines by picking a random and unique color per chain.  Line with can 
    be adjusted using *lw* argument, e.g. ``lw=12``. Default width is 4.  
    Chain colors can be overwritten using chain identifier as in ``A='green'``.  
      
    Water molecule oxygen atoms are represented by red colored circles.  Color 
    can be changed using *water* keyword argument, e.g. ``water='aqua'``.
    Water marker and size can be changed using *wmarker* and *wsize* keywords, 
    defaults values are ``wmarker='.', wsize=6``.
    
    Hetero atoms matching ``"hetero and noh"`` selection are represented by 
    circles and unique colors are picked at random on a per residue basis.  
    Colors can be customized using residue name as in ``NAH='purple'``.  Note 
    that this will color all distinct residues with the same name in the same 
    color.  Hetero atom marker and size can be changed using *hmarker* and 
    *hsize* keywords, default values are ``hmarker='o', hsize=6``. 

    ProDy will set the size of axis so the representation is not distorted when
    the figure window is close to a square.  Colors are picked at random,
    except for water oxygens which will always be colored red.
    

    .. plot::
       :context:
       :include-source:
       
       p38 = parsePDB('1p38')
       p38inh = parsePDB('1zz2')
       matchAlign(p38inh, p38)
       plt.figure(figsize=(5,4))
       showProtein(p38, p38inh)
       plt.legend(prop={'size': 10})
        
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    """
    
    alist = atoms
    for atoms in alist:    
        if not isinstance(atoms, prody.Atomic):
            raise TypeError('atoms must be an Atomic instance')
    if not plt: prody.importPyPlot()
    cf = plt.gcf()
    show = None
    for child in cf.get_children():
        if isinstance(child, Axes3D):
            show = child
            break 
    if show is None:
        show = Axes3D(cf)
    from matplotlib import colors
    cnames = dict(colors.cnames)
    wcolor = kwargs.get('water', 'red').lower()
    avoid = np.array(colors.hex2color(cnames.pop(wcolor, cnames.pop('red'))))
    for cn, val in cnames.items():
        clr = np.array(colors.hex2color(val))
        if clr.sum() > 2.4:
            cnames.pop(cn)
        elif np.abs(avoid - clr).sum() <= 0.6:
            cnames.pop(cn)
    cnames = cnames.keys()
    import random
    random.shuffle(cnames)
    min_ = list()
    max_ = list()
    for atoms in alist:
        if isinstance(atoms, prody.AtomGroup):
            title = atoms.getTitle()
        else:
            title = atoms.getAtomGroup().getTitle()
        calpha = atoms.select('calpha')
        if calpha:
            for ch in prody.HierView(calpha, chain=True):
                xyz = ch._getCoords()
                chid = ch.getIdentifier()
                show.plot(xyz[:,0], xyz[:,1], xyz[:,2], 
                          label=title + '_' + chid,
                          color=kwargs.get(chid, cnames.pop()).lower(),
                          lw=kwargs.get('lw', 4))
        water = atoms.select('water and noh')
        if water: 
            xyz = atoms.select('water')._getCoords()
            show.plot(xyz[:,0], xyz[:,1], xyz[:,2], label=title + '_water',
                      color=wcolor, 
                      ls='None', marker=kwargs.get('wmarker', '.'), 
                      ms=kwargs.get('wsize', 6))
        hetero = atoms.select('not protein and not nucleic and not water')
        if hetero:
            for res in prody.HierView(hetero).iterResidues():
                xyz = res._getCoords()
                resname = res.getName()
                resnum = str(res.getNumber())
                chid = res.getChid()
                show.plot(xyz[:,0], xyz[:,1], xyz[:,2], ls='None',
                          color=kwargs.get(resname, cnames.pop()).lower(), 
                          label=title + '_' + chid + '_' + resname + resnum, 
                          marker=kwargs.get('hmarker', 'o'), 
                          ms=kwargs.get('hsize', 6))
        xyz = atoms._getCoords()
        min_.append(xyz.min(0))
        max_.append(xyz.max(0))

    show.set_xlabel('x')
    show.set_ylabel('y')
    show.set_zlabel('z')
    min_ = np.array(min_).min(0)
    max_ = np.array(max_).max(0)
    center = (max_ + min_) / 2
    half = (max_ - min_).max() / 2
    show.set_xlim3d(center[0]-half, center[0]+half)
    show.set_ylim3d(center[1]-half, center[1]+half)
    show.set_zlim3d(center[2]-half, center[2]+half)
    if kwargs.get('legend', False):
        show.legend(prop={'size': 10})
    return show
