# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

"""This module defines classes and functions to fetch, parse, 
and write PDB files, and also to blast search 
`ProteinDataBank <http://wwpdb.org>`_.

Classes
-------

  * :class:`PDBBlastRecord`

Functions
---------

  * :func:`applyBiomolecularTransformations`
  * :func:`assignSecondaryStructure`
  * :func:`blastPDB`
  * :func:`fetchPDB`
  * :func:`getPDBMirrorPath`
  * :func:`getWWPDBFTPServer`
  * :func:`setPDBMirrorPath`
  * :func:`setWWPDBFTPServer` 
  * :func:`parsePDB`
  * :func:`parsePDBStream`
  * :func:`parsePSF`
  * :func:`writePDB`
  * :func:`writePDBStream`
  
  * :func:`fetchLigandData`
  * :func:`fetchPDBClusters`
  * :func:`loadPDBClusters`
  * :func:`getPDBCluster`

  * :func:`execDSSP`
  * :func:`parseDSSP`
  * :func:`performDSSP`
  
  * :func:`execSTRIDE`
  * :func:`parseSTRIDE`
  * :func:`performSTRIDE`

  
   
   
.. doctest::
    :hide:
        
    >>> from prody import *
    >>> import numpy as np
    >>> allmodels = parsePDB('2k39', subset='ca')
    >>> model10 = parsePDB('2k39', subset='ca', model=10)
    >>> np.all(allmodels.getCoordsets(9) == model10.getCoordinates())
    True
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import gzip
import os.path
import time
import os
import shutil
from glob import glob
from collections import defaultdict

ET = None
urllib2 = None

import numpy as np

BioBlast = None


import prody
LOGGER = prody.LOGGER
from prody.atomic import *


PDB_CLUSTERS = {30: None, 40: None, 50: None, 70: None, 
                90: None, 95: None, 100: None}
PDB_CLUSTERS_UPDATE_WARNING = True


__all__ = ['PDBBlastRecord', 
           'assignSecondaryStructure',
           'applyBiomolecularTransformations',
           'blastPDB', 'fetchPDB', 
           'getPDBMirrorPath', 'getWWPDBFTPServer', 
           'setPDBMirrorPath', 'setWWPDBFTPServer',
           'parsePDBStream', 'parsePDB', 'parsePSF',
           'writePDBStream', 'writePDB',
           'fetchLigandData',
           'execDSSP', 'parseDSSP', 'performDSSP',
           'execSTRIDE', 'parseSTRIDE', 'performSTRIDE',
           'fetchPDBClusters', 'loadPDBClusters', 'getPDBCluster',
           ]

class PDBParserError(Exception):    
    pass


def _makePath(path):
    """Make all directories that does not exist in a given path."""
    if os.path.isabs(path):
        path = os.path.relpath(path)
    if not os.path.isdir(path):
        dirs = path.split(os.sep)
        for i in range(len(dirs)):
            dirname = os.sep.join(dirs[:i+1])
            try:
                if not os.path.isdir(dirname): 
                    os.mkdir(dirname)
            except OSError:
                raise OSError('{0:s} could not be created, please '
                              'specify another path'.format(path))
                return os.getcwd()
    return os.path.join(os.getcwd(), path)


_pdb_extensions = set(['.pdb', '.PDB', '.gz', '.GZ', '.ent', '.ent.gz'])
_WWPDB_RCSB = ('RCSB PDB (USA)', 'ftp.wwpdb.org', 
    '/pub/pdb/data/structures/divided/pdb/')
_WWPDB_PDBe = ('PDBe (Europe)', 'ftp.ebi.ac.uk', 
    '/pub/databases/rcsb/pdb/data/structures/divided/pdb/')
_WWPDB_PDBj = ('PDBj (Japan)', 'pdb.protein.osaka-u.ac.jp', 
    '/pub/pdb/data/structures/divided/pdb/')
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

def getPDBMirrorPath():
    """Set the path to a local PDB mirror.
    
    .. versionadded:: 0.6.1
    
    """

    path = prody._ProDySettings.get('pdb_mirror_path')
    if isinstance(path, str):
        if os.path.isdir(path):
            return path
        else:
            LOGGER.warning('PDB mirror path "{0:s}" is not a accessible.'
                           .format(path))


def setPDBMirrorPath(path):
    """Return the path to a local PDB mirror.
    
    .. versionadded:: 0.6.1
    
    Returns ``None`` if a PDB local mirror path is not set.
    
    """
    
    path = str(path)
    if os.path.isdir(path):
        prody._ProDySettings['pdb_mirror_path'] = path
        prody._saveProDySettings()
    else:
        LOGGER.warning('{0:s} is not a valid path.'.format(path))


def setWWPDBFTPServer(key):
    """Set the PDB FTP server used for downloading PDB structures when needed.
    
    .. versionadded:: 0.6.1
    
    Use one of the following keywords for setting a different server.
    
    +---------------------------+-----------------------------+
    | WWPDB FTP server          | *Key* (case insensitive)    |
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
        prody._ProDySettings['wwpdb_ftp'] = server
        prody._saveProDySettings()
    else:
        LOGGER.warning('{0:s} is not a valid key.'.format(key))

def getWWPDBFTPServer():
    """Return a tuple containing name, host, and path of the currently 
    set WWPDB FTP server.
    
    .. versionadded:: 0.6.1
    
    """
    
    server = prody._ProDySettings.get('wwpdb_ftp')
    if server is None:
        LOGGER.warning('A WWPDB FTP server is not set by the user. '
                       'Default FTP server RCSB PDB is returned. Use '
                       'setWWPDBFTPServer function for choosing a server '
                       'physically close to your location.')
        return _WWPDB_RCSB
    else:
        return server

def fetchPDB(pdb, folder='.', compressed=True, copy=False):
    """Return the path(s) to PDB file(s) for specified identifier(s).

    *pdb* may be a list of PDB identifiers or an identifier string.

    .. versionchanged:: 0.8
       *compressed* and *copy* argument is introduced.  
    
    .. versionchanged:: 0.8.2
       When *compressed* is false, compressed files found in *folder* or 
       local PDB mirror are decompressed.

    If *folder* already contains a PDB file matching given identifier, a 
    file will not be downloaded and the path to the existing file
    will be returned.
    
    If a file matching the given PDB identifier is not found in *folder*,
    the PDB file will be sought in the local PDB mirror, if a local
    mirror is set by the user. When PDB is found in local repository, 
    the path to the file will be returned. 
    
    Finally, if PDB file is not found in *folder* or local mirror,
    it will be downloaded from the user-selected WWPDB FTP server.
    User can set one of the WWPDB FTP servers using :func:`setWWPDBFTPServer`. 
    Downloaded files will be saved in *folder*. FTP servers provides gunzipped 
    PDB files.
    
    If *compressed* argument is set to ``False``, downloaded files will be 
    decompressed. When a compressed file is found in the *folder*, it will
    also be decompressed.
    
    For PDB files found in a local mirror of PDB, setting *copy* ``True`` will
    copy them from the mirror to the user specified *folder*.  
        
    """
    
    if isinstance(pdb, str):
        identifiers = [pdb]
    elif isinstance(pdb, list):
        identifiers = pdb
    else:
        raise TypeError('pdb may be a string or a list of strings')
        
        
    if folder != '.':
        folder = _makePath(folder)
    if not os.access(folder, os.W_OK):
        raise IOError('permission to write in {0:s} is denied, please '
                      'specify another folder'.format(folder))
    
    filenames = []
    exists = 0
    success = 0
    failure = 0
    download = False

    pdbfnmap = {}
    for pdbfn in glob(os.path.join(folder, '*.pdb*')): 
        if os.path.splitext(pdbfn)[1] in _pdb_extensions:
            pdbfnmap[os.path.split(pdbfn)[1].split('.')[0].lower()] = pdbfn
    for pdbfn in glob(os.path.join(folder, '*.PDB*')):
        if os.path.splitext(pdbfn)[1] in _pdb_extensions:
            pdbfnmap[os.path.split(pdbfn)[1].split('.')[0].lower()] = pdbfn
                
    mirror_path = getPDBMirrorPath()
    for i, pdbid in enumerate(identifiers):
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
        identifiers[i] = pdbid
        fn = pdbfnmap.get(pdbid, None)
        if fn:
            fn = os.path.relpath(fn)
            if not compressed:
                temp, ext = os.path.splitext(fn) 
                if ext == '.gz':
                    fn = gunzip(fn, temp)
            filenames.append(fn)
            LOGGER.debug('{0:s} ({1:s}) is found in the working directory.'
                         .format(pdbid, fn))
            exists += 1
            continue
        if mirror_path is not None and os.path.isdir(mirror_path):
            fn = os.path.join(mirror_path, 'data/structures/divided/pdb',
                    pdbid[1:3], 'pdb' + pdbid + '.ent.gz')
            if os.path.isfile(fn):
                if copy or not compressed:
                    if compressed:
                        filename = os.path.join(folder, pdbid + '.pdb.gz')
                        shutil.copy(fn, filename)
                    else:
                        filename = os.path.join(folder, pdbid + '.pdb')
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
        filenames.append(pdbid)
        download = True
    if download:
        from ftplib import FTP
        ftp_name, ftp_host, ftp_path = getWWPDBFTPServer()
        LOGGER.debug('Connecting WWPDB FTP server {0:s}.'.format(ftp_name))
        try:
            ftp = FTP(ftp_host)
        except Exception as error:
            raise type(error)('FTP connection problem, potential reason: '
                              'no internet connectivity')
        else:
            ftp.login('')
            for i, pdbid in enumerate(identifiers):
                if pdbid != filenames[i]:
                    continue
                if compressed:
                    filename = os.path.join(folder, pdbid + '.pdb.gz')
                else:
                    filename = os.path.join(folder, pdbid + '.pdb')
                pdbfile = open(filename, 'w+b')
                try:
                    ftp.cwd(os.path.join(ftp_path, pdbid[1:3]))
                    ftp.retrbinary('RETR pdb{0:s}.ent.gz'.format(pdbid), 
                                   pdbfile.write)
                except Exception as error:
                    pdbfile.close()
                    os.remove(filename)
                    if 'pdb{0:s}.ent.gz'.format(pdbid) in ftp.nlst():
                        LOGGER.debug('{0:s} download failed ({1:s}). It '
                                     'is possible that you don\'t have '
                                     'rights to download .gz files in the '
                                     'current network.'.format(pdbid, 
                                     str(error)))
                    else:
                        LOGGER.debug('{0:s} download failed. pdb{0:s}.ent.'
                                     'gz does not exist on ftp.wwpdb.org.'
                                     .format(pdbid))
                    failure += 1
                    filenames[i] = None 
                else:
                    pdbfile.close()
                    if not compressed:
                        gunzip(filename)
                    filename = os.path.relpath(filename)
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

def gunzip(filename, outname=None):
    """Decompresses *filename* and saves as *outname*. 
    
    When *outname* is ``None``, *filename* is used as the output name.
    
    Returns output filename upon successful completion.
    
    """

    if not isinstance(filename, str):
        raise TypeError('filename must be a string')
    if not os.path.isfile(filename):
        raise ValueError('{0:s} does not exist'.format(filename))
    if outname is None: 
        outname = filename
    inp = gzip.open(filename, 'rb')
    data = inp.read()
    inp.close()
    out = open(outname, 'w')
    out.write(data)
    out.close()
    return outname

_parsePDBdoc = """
    :arg model: model index or None (read all models), 
        e.g. ``model=10``
    :type model: int, list

    :arg header: If ``True`` PDB header content will be parsed and returned.
    :type header: bool

    :arg chain: Chain identifiers for parsing specific chains, e.g. 
        ``chain='A'``, ``chain='B'``, ``chain='DE'``. By default all chains
        are parsed.        
    :type chain: str

    :arg subset: A predefined keyword to parse subset of atoms.
        Valid keywords are ``"calpha"`` (``"ca"``) or ``"backbone"`` 
        (``"bb"``), or ``None`` (read all atoms), e.g. ``subset='bb'``
    :type subset: str

    :arg altloc: If a location indicator is passed, such as ``'A'`` or ``'B'``, 
         only indicated alternate locations will be parsed as the single 
         coordinate set of the AtomGroup. If ``True`` all alternate locations 
         will be parsed and each will be appended as a distinct coordinate set.
         Default is ``"A"``.
    :type altloc: str

    :arg name: Name of the AtomGroup instance. When ``None`` is passed,
        AtomGroup is named after the PDB filename.  
    :type name: str
    
    :arg biomol: If ``True``, return biomolecule obtained by transforming the
        coordinates using information from header section.
    :type biomol: False

    :arg secondary: If ``True``, parse the secondary structure information
        from header section and assign data to atoms.
    :type secondary: False
    
    :arg ag: :class:`~prody.atomic.AtomGroup` instance for storing data parsed 
        from PDB file. Number of atoms in *ag* and number of atoms parsed from
        the PDB file must be the same. Atoms in *ag* and the PDB file must be 
        in the same order. Non-coordinate data stored in *ag* will be 
        overwritten with those parsed from the file. 
    :type ag: :class:`~prody.atomic.AtomGroup`
    

    If ``model=0`` and ``header=True``, return header 
    dictionary only.
    
    .. versionchanged:: 0.6
       Default behavior for parsing alternate locations have changed. 
       Alternate locations indicated by ``A`` are parsed.
    
    .. versionchanged:: 0.7.1
       *name* is now a keyword argument. *biomol* and *secondary* keyword 
       arguments makes the parser return the biological molecule and/or
       assign secondary structure information.
       
    .. versionchanged:: 0.8.1
       User can pass an :class:`~prody.atomic.AtomGroup` instance as *ag*.

    """
    
_PDBSubsets = ['ca', 'calpha', 'bb', 'backbone']

def parsePDB(pdb, model=None, header=False, chain=None, subset=None, 
             altloc='A', **kwargs):
    """Return an :class:`~prody.atomic.AtomGroup` and/or 
    dictionary containing header data parsed from a stream of PDB lines. 
    
    This function extends :func:`parsePDBStream`.
    
    |example| See :ref:`parsepdb` for a detailed example.
    
    :arg pdb: A valid PDB identifier or filename.  
        If needed, PDB files are downloaded using :func:`fetchPDB()` function.  
        
    """
    
    name = kwargs.get('name', None)
    if not os.path.isfile(pdb):
        if len(pdb) == 4 and pdb.isalnum():
            if name is None:
                name = pdb.lower()
                kwargs['name'] = name
            filename = fetchPDB(pdb)
            if filename is None:
                raise IOError('PDB file for {0:s} could not be downloaded.'
                              .format(pdb))
            pdb = filename
        else:
            raise IOError('{0:s} is not a valid filename or a valid PDB '
                          'identifier.'.format(pdb))
    if name is None:
        fn, ext = os.path.splitext(os.path.split(pdb)[1])
        if ext == '.gz':
            fn, ext = os.path.splitext(fn)
        name = fn.lower()
        kwargs['name'] = name
    if pdb.endswith('.gz'):
        pdb = gzip.open(pdb)
    else:
        pdb = open(pdb)
    result = parsePDBStream(pdb, model, header, chain, subset, 
                            altloc, **kwargs)
    pdb.close()
    return result

parsePDB.__doc__ += _parsePDBdoc
    
def parsePDBStream(stream, model=None, header=False, chain=None, subset=None, 
                   altloc='A', **kwargs):
    """Return an :class:`~prody.atomic.AtomGroup` and/or 
    dictionary containing header data parsed from a stream of PDB lines. 
    
    :arg stream: Anything that implements the method readlines() 
        (e.g. :class:`file`, buffer, stdin).

    """
    
    if model is not None:
        if isinstance(model, int):
            if model < 0:
                raise ValueError('model must be greater than 0')
        else:
            raise TypeError('model must be an integer, {0:s} is invalid'
                            .format(str(model)))
    if subset is not None: 
        if not isinstance(subset, str):
            raise TypeError('subset must be a string')
        elif subset.lower() not in _PDBSubsets:
            raise ValueError('"{0:s}" is not a valid subset'.format(subset))
    if chain is not None:
        if not isinstance(chain, str):
            raise TypeError('chain must be a string')
        elif len(chain) == 0:
            raise ValueError('chain must not be an empty string')
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, prody.AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.getNumOfCoordsets()
    else:
        ag = prody.AtomGroup(kwargs.get('name', 'unknown'))
        n_csets = 0
    biomol = kwargs.get('biomol', False)
    secondary = kwargs.get('secondary', False)
    split = 0
    hd = None
    lines = stream.readlines()
    if header or biomol or secondary:
        hd, split = _getHeaderDict(lines)
    if model != 0:
        start = time.time()
        _parsePDBLines(ag, lines, split, model, chain, subset, altloc)
        if ag.getNumOfAtoms() > 0:
            LOGGER.info('{0:d} atoms and {1:d} coordinate sets were '
                        'parsed in {2:.2f}s.'.format(ag.getNumOfAtoms(), 
                         ag.getNumOfCoordsets() - n_csets, time.time()-start))
        else:
            ag = None
            LOGGER.warning('Atomic data could not be parsed, please '
                           'check the input file.')
    if secondary:
        try:
            ag = assignSecondaryStructure(hd, ag)
        except:
            raise PDBParserError('secondary structure assignments could not '
                                 'be made, check input file')
    if biomol:
        try:
            ag = applyBiomolecularTransformations(hd, ag)
        except:
            raise PDBParserError('biomolecule could not be generated, check'
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

def _parsePDBLines(atomgroup, lines, split, model, chain, subset, altloc_torf):
    """Return an AtomGroup. See also :func:`parsePDBStream()`.
    
    :arg lines: PDB lines 
    :arg split: starting index for coordinate data lines"""
    
    if subset is not None:
        subset = subset.lower()
        if subset in ('calpha', 'ca'):
            subset = set(('CA',))
        elif subset in ('backbone', 'bb'):
            subset = set(prody.getBackboneAtomNames())
        only_subset = True
        protein_resnames = set(prody.getProteinResidueNames())
    else:
        only_subset = False
    if chain is None:
        only_chains = False
    else:
        only_chains = True
    onlycoords = False
    n_atoms = atomgroup.getNumOfAtoms()
    if n_atoms > 0:
        asize = n_atoms
    else:
        # most PDB files contain less than 99999 atoms
        asize = min(len(lines) - split, 99999)
    alength = asize
    coordinates = np.zeros((asize, 3), dtype=np.float64)
    atomnames = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['chain'].dtype)
    bfactors = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['beta'].dtype)
    occupancies = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['occupancy'].dtype)
    hetero = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['hetero'].dtype)
    altlocs = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['altloc'].dtype)
    segnames = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['segment'].dtype)
    elements = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['element'].dtype)
    secondary = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['secondary'].dtype)
    anisou = np.zeros((asize, 6), dtype=ATOMIC_DATA_FIELDS['anisou'].dtype)
    siguij = np.zeros((asize, 6), dtype=ATOMIC_DATA_FIELDS['siguij'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_DATA_FIELDS['serial'].dtype)
    asize = 2000 # increase array length by this much when needed 
        
    start = split
    stop = len(lines)
    nmodel = 0
    if model is not None and model != 1:
        for i in range(split, len(lines)):
            if lines[i][:5] == 'MODEL':
                nmodel += 1
                if model == nmodel:
                    start = i+1 
                    stop = len(lines)
                    break
        if nmodel != model:
            raise PDBParserError('model {0:d} is not found'.format(model))
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
    is_anisou = False
    is_siguij = False
    is_scndry = False
    altloc = defaultdict(list)
    i = start
    while i < stop:
        line = lines[i]
        startswith = line[0:6]
        
        if startswith == 'ATOM  ' or startswith == 'HETATM':
            atomname = line[12:16].strip()
            resname = line[17:21].strip()
            if only_subset:
                if not (atomname in subset and resname in protein_resnames): 
                    i += 1
                    continue
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
                coordinates[acount, 0] = float(line[30:38])
                coordinates[acount, 1] = float(line[38:46])
                coordinates[acount, 2] = float(line[46:54])
            except:
                if acount >= n_atoms > 0:
                    if nmodel ==0:
                        raise ValueError('PDB file and AtomGroup ag must have '
                                         'same number of atoms')
                    LOGGER.warning('Discarding model {0:d}, which contains '
                                   'more atoms than first model does.'
                                   .format(nmodel+1))
                    acount = 0
                    nmodel += 1
                    coordinates = np.zeros((n_atoms, 3), dtype=np.float64)
                    while lines[i][:6] != 'ENDMDL':
                        i += 1
                else:
                    raise PDBParserError('invalid or missing coordinate(s) at '
                                         'line {0:d}.'.format(i+1))
            if onlycoords:
                acount += 1
                i += 1
                continue
            
            serials[acount] = int(line[6:11])
            altlocs[acount] = alt 
            atomnames[acount] = atomname
            
            # resname = line[17:21], but some are 4 chars long
            resnames[acount] = resname
            chainids[acount] = chid
            resnums[acount] = int(line[22:26].split()[0])
            icodes[acount] = line[26].strip()
            try:
                occupancies[acount] = float(line[54:60])
            except:
                LOGGER.warning('failed to parse occupancy at line {0:d}'
                               .format(i))
            try:
                bfactors[acount] = float(line[60:66])
            except:
                LOGGER.warning('failed to parse beta-factor at line {0:d}'
                               .format(i))
            
            if startswith[0] == 'H':
                hetero[acount] = True

            segnames[acount] = line[72:76].strip()
            elements[acount] = line[76:78].strip()
            acount += 1
            if n_atoms == 0 and acount >= alength:
                alength += asize
                coordinates = np.concatenate(
                    (coordinates, np.zeros((asize, 3), np.float64)))
                atomnames = np.concatenate((atomnames,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['name'].dtype)))
                resnames = np.concatenate((resnames,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['resname'].dtype)))
                resnums = np.concatenate((resnums,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['resnum'].dtype)))
                chainids = np.concatenate((chainids,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['chain'].dtype)))
                bfactors = np.concatenate((bfactors,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['beta'].dtype)))
                occupancies = np.concatenate((occupancies,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['occupancy'].dtype)))
                hetero = np.concatenate((hetero,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['hetero'].dtype)))
                altlocs = np.concatenate((altlocs,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['altloc'].dtype)))
                segnames = np.concatenate((segnames,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['segment'].dtype)))
                elements = np.concatenate((elements,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['element'].dtype)))
                secondary = np.concatenate((secondary,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['secondary'].dtype)))
                anisou = np.concatenate((anisou,
                    np.zeros((asize, 6), ATOMIC_DATA_FIELDS['anisou'].dtype)))
                siguij = np.concatenate((siguij,
                    np.zeros((asize, 6), ATOMIC_DATA_FIELDS['siguij'].dtype)))
                icodes = np.concatenate((icodes,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['icode'].dtype)))
                serials = np.concatenate((serials,
                    np.zeros(asize, ATOMIC_DATA_FIELDS['serial'].dtype)))
        #elif startswith == 'END   ' or startswith == 'CONECT':
        #    i += 1
        #    break
        elif startswith == 'ENDMDL' or startswith[:3] == 'END':
            if acount == 0:
                # If there is no atom record between ENDMDL and END skip to next
                i += 1
                continue
            elif model is not None:
                i += 1
                break
            elif onlycoords:
                if acount < n_atoms:
                    LOGGER.warning('Discarding model {0:d}, which contains '
                                   '{1:d} fewer atoms than the first model '
                                   'does.'.format(nmodel+1, n_atoms-acount))
                else:
                    atomgroup.addCoordset(coordinates)
                nmodel += 1
                acount = 0
                coordinates = np.zeros((n_atoms, 3), dtype=np.float64)
            else:
                if acount != n_atoms > 0:
                    raise ValueError('PDB file and AtomGroup ag must have '
                                    'same number of atoms')
                atomgroup.addCoordset(coordinates[:acount])
                atomgroup.setAtomNames(atomnames[:acount])
                atomgroup.setResidueNames(resnames[:acount])
                atomgroup.setResidueNumbers(resnums[:acount])
                atomgroup.setChainIdentifiers(chainids[:acount])
                atomgroup.setTempFactors(bfactors[:acount])
                atomgroup.setOccupancies(occupancies[:acount])
                atomgroup.setHeteroFlags(hetero[:acount])
                atomgroup.setAltLocIndicators(altlocs[:acount])
                atomgroup.setSegmentNames(segnames[:acount])
                atomgroup.setElementSymbols(elements[:acount])
                atomgroup.setInsertionCodes(icodes[:acount])
                atomgroup.setSerialNumbers(serials[:acount])
                if is_scndry:
                    atomgroup.setSecondaryStrs(secondary[:acount])
                if is_anisou:
                    atomgroup.setAnisoTempFactors(anisou[:acount] / 10000)
                if is_siguij:
                    atomgroup.setAnisoStdDevs(siguij[:acount] / 10000)
                
                onlycoords = True
                nmodel += 1
                n_atoms = acount 
                acount = 0
                coordinates = np.zeros((n_atoms, 3), dtype=np.float64)
                if altloc and altloc_torf:
                    _evalAltlocs(atomgroup, altloc, chainids, resnums, 
                                 resnames, atomnames)
                    altloc = defaultdict(list)
        elif startswith == 'ANISOU':
            is_anisou = True
            try:
                index = acount - 1
                anisou[index, 0] = float(line[28:35])
                anisou[index, 1] = float(line[35:42])
                anisou[index, 2] = float(line[43:49])
                anisou[index, 3] = float(line[49:56])
                anisou[index, 4] = float(line[56:63])
                anisou[index, 5] = float(line[63:70])
            except:
                LOGGER.warning('failed to parse anisotropic temperature '
                    'factors at line {0:d}'.format(i))
        elif startswith =='SIGUIJ':
            is_siguij = True
            try:
                index = acount - 1
                siguij[index, 0] = float(line[28:35])
                siguij[index, 1] = float(line[35:42])
                siguij[index, 2] = float(line[43:49])
                siguij[index, 3] = float(line[49:56])
                siguij[index, 4] = float(line[56:63])
                siguij[index, 5] = float(line[63:70])
            except:
                LOGGER.warning('failed to parse standard deviations of '
                    'anisotropic temperature factors at line {0:d}'.format(i))
        elif startswith =='SIGATM':
            pass
        i += 1
    if onlycoords:
        if acount == atomgroup.getNumOfAtoms():
            atomgroup.addCoordset(coordinates)
    else:            
        atomgroup.setCoordinates(coordinates[:acount])
        atomgroup.setAtomNames(atomnames[:acount])
        atomgroup.setResidueNames(resnames[:acount])
        atomgroup.setResidueNumbers(resnums[:acount])
        atomgroup.setChainIdentifiers(chainids[:acount])
        atomgroup.setTempFactors(bfactors[:acount])
        atomgroup.setOccupancies(occupancies[:acount])
        atomgroup.setHeteroFlags(hetero[:acount])
        atomgroup.setAltLocIndicators(altlocs[:acount])
        atomgroup.setSegmentNames(segnames[:acount])
        atomgroup.setElementSymbols(elements[:acount])
        atomgroup.setInsertionCodes(icodes[:acount])
        atomgroup.setSerialNumbers(serials[:acount])
        if is_scndry:
            atomgroup.setSecondaryStrs(secondary[:acount])
        if is_anisou:
            atomgroup.setAnisoTempFactors(anisou[:acount] / 10000)
        if is_siguij:
            atomgroup.setAnisoStdDevs(siguij[:acount] / 10000)

    if altloc and altloc_torf:
        _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames)
                
    return atomgroup
    
def _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames):
    altloc_keys = altloc.keys()
    altloc_keys.sort()
    indices = {}
    for key in altloc_keys:
        xyz = atomgroup.getCoordinates()
        success = 0
        lines = altloc[key]
        for line, i in lines:
            #-->
            #try:
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
                    LOGGER.warning('failed to parse alternate location {0:s} '
                                   'at line {1:d}, residue does not exist as '
                                   'altloc A'.format(key, i+1))
                    continue
                rn = resnames[ids[0]]
                ans = atomnames[ids]
                indices[(ach, ari)] = (rn, ids, ans)
            if rn != arn:
                LOGGER.warning('failed to parse alternate location {0:s} at '
                               'line {1:d}, residue names do not match '
                               '(expected {2:s}, parsed {3:s})'
                               .format(key, i+1, rn, arn))
                continue
            index = ids[(ans == aan).nonzero()[0]]
            if len(index) != 1:
                LOGGER.warning('failed to parse alternate location {0:s} at '
                               'line {1:d}, could not identify matching atom '
                               '({2:s} not found in the residue)'
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
            #except Exception as exception:
            #    print i, line
            #    print exception
            #-->
        LOGGER.info('{0:d} out of {1:d} alternate location {2:s} lines were '
                    'parsed successfully.'.format(success, len(lines), key))
        if success > 0:
            LOGGER.info('Alternate location {0:s} is appended as a coordinate '
                        'set to the atom group.'
                        .format(key, atomgroup.getName()))
            atomgroup.addCoordset(xyz)
    
def _getHeaderDict(lines):
    """Return header data in a dictionary."""
    header = {}
    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    helix = {}
    sheet = {}
    biomolecule = defaultdict(list)
    sequences = defaultdict(str)
    applyToChains = (' ')
    i = 0
    for i, line in enumerate(lines):        
        line = lines[i]
        startswith = line[0:6]
        if startswith == 'ATOM  ' or startswith == 'HETATM' or \
           startswith == 'MODEL ':
            if biomolecule:
                header['biomolecular_transformations'] = biomolecule
            if sequences:
                header['sequences'] = dict(sequences)
            if helix:
                header['helix'] = helix
            if sheet:
                header['sheet'] = sheet
            return header, i
        
        if startswith == 'HELIX ':
            try:
                chid = line[19]
                # helix class, serial number, identifier
                value = (int(line[38:40]), int(line[7:10]), 
                         line[11:14].strip())
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
        elif startswith == 'SHEET ':
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
            
        elif startswith == 'HEADER':
            header['deposition_date'] = line[50:59].strip()
            header['classification'] = line[10:50].strip()
            header['identifier'] = line[62:66]
        elif startswith == 'TITLE ':
            temp = header.get('title', '')
            temp += ' ' + line[10:]
            header['title'] = temp.strip()
        elif startswith == 'EXPDTA':
            temp = header.get('experiment', '')
            temp += ' ' + line[10:]
            header['experiment'] = temp.strip()
        elif startswith == 'AUTHOR':
            temp = header.get('authors', [])
            temp += line[10:].strip().split(',')
            while '' in temp:
                temp.pop(temp.index(''))
            header['authors'] = temp
        elif startswith == 'SOURCE':
            temp = header.get('source', '')
            temp += ' ' + line[10:]
            header['source'] = temp.strip()
        elif startswith == 'SEQRES':
            sequences[line[11]] += ''.join(
                                prody.compare.getSequence(line[19:].split()))
        elif startswith == 'REMARK':
            nmbr = line[7:10]
            if nmbr == '  2':
                if 'RESOLUTION' in line:
                    header['resolution'] = line[23:].strip()[:-1]    
            elif nmbr == '350':
                if line[13:18] == 'BIOMT':
                    biomt = biomolecule[currentBiomolecule]
                    if len(biomt) == 0:
                        biomt.append(applyToChains)
                    biomt.append(line[23:])
                elif line[11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                    applyToChains = line[41:].replace(' ', '').strip().split(',')
                elif line[11:23] == 'BIOMOLECULE:': 
                    currentBiomolecule = line.split()[-1]
        elif startswith == 'COMPND':
            temp = header.get('compounds', '')
            temp += ' ' + line[10:]
            header['compounds'] = temp.strip()
        elif startswith == 'JRNL  ':
            temp =  header.get('reference', {})
            items = line.split()
            if items[1] == 'AUTH':
                tmp = temp.get('authors', [])
                tmp += line[19:].strip().split(',')
                while '' in tmp:
                    tmp.pop(tmp.index(''))
                temp['authors'] = tmp
            elif items[1] == 'TITL':
                tmp = temp.get('title', '')
                tmp += ' ' + line[19:]
                temp['title'] = tmp.strip()
            elif items[1] == 'REF':
                tmp = temp.get('reference', '')
                tmp += ' ' + line[19:]
                temp['reference'] = ' '.join(tmp.split())
            elif items[1] == 'PMID':
                temp['pubmed'] = line[19:].strip()
            elif items[1] == 'DOI':
                temp['doi'] = line[19:].strip()
            header['reference'] = temp
    if biomolecule:
        header['biomolecular_transformations'] = biomolecule
    if sequences:
        header['sequences'] = dict(sequences)
    if helix:
        header['helix'] = helix
    if sheet:
        header['sheet'] = sheet
    return header, i

def parsePSF(filename, name=None, ag=None):
    """Return an :class:`~prody.atomic.AtomGroup` instance storing data 
    parsed from X-PLOR format PSF file *filename*. If *name* is not given, 
    *filename* will be set as the name of the :class:`AtomGroup` instance. 
    An :class:`AtomGroup` instance may be provided as *ag* argument. When 
    provided, *ag* must have the same number of atoms in the same order as 
    the file. Data from PSF file will be added to *ag*. This may overwrite 
    present data if it overlaps with PSF file content.
    
    .. versionadded:: 0.8.1
    
    """
    
    if ag is not None:
        if not isinstance(ag, prody.AtomGroup):
            raise TypeError('ag must be an AtomGroup instance') 
    
    psf = open(filename)
    line = psf.readline()
    i_line = 1
    while line:
        line = line.strip()
        if line.strip().endswith('!NATOM'):
            n_atoms = int(line.split('!')[0])
            break
        line = psf.readline()
        i_line += 1
    if name is None:
        name = os.path.splitext(os.path.split(filename)[1])[0]
    else:
        name = str(name)
    if ag is None:
        ag = prody.AtomGroup(name)
    else:
        if n_atoms != ag.getNumOfAtoms():
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
    ag.setSerialNumbers(serials)
    ag.setSegmentNames(segnames)
    ag.setResidueNumbers(resnums)
    ag.setResidueNames(resnames)
    ag.setAtomNames(atomnames)
    ag.setAtomTypes(atomtypes)
    ag.setCharges(charges)
    ag.setMasses(masses)
    return ag

class PDBBlastRecord(object):

    """A class to store results from ProteinDataBank blast search."""
    
    __slots__ = ['_record', '_sequence']

    def __init__(self, sequence, blast_record):
        """Instantiate a PDBlast object instance.
        
        :arg sequence: string of one-letter amino acid code
        :arg blast_record: record from parsing results using NCBIXML.parse
        
        """
        self._sequence = sequence
        self._record = blast_record
        
    def getSequence(self):    
        """Return the query sequence that was used in the search."""
        return self._sequence
        
    def printSummary(self):
        """Prints a summary of the results to the screen."""
        record = self._record
        LOGGER.info('Blast results summary:')
        LOGGER.info('  {0:s}: {1}'.format('database', record.database))
        LOGGER.info('  {0:s}: {1}'.format('expect', record.expect))
        LOGGER.info('  {0:s}: {1}'.format('matrix', record.matrix))
        LOGGER.info('  {0:s}: {1}'.format('query', record.query))
        LOGGER.info('  {0:s}: {1}'.format('query length', record.query_length))
        LOGGER.info('  {0:s}: {1}'.format('blast version', record.version))

    def getHits(self, percent_identity=90., percent_coverage=70.):
        """Return a dictionary that contains hits.
        
        Returns a dictionary whose keys are PDB identifiers. Each dictionary
        entry is also a dictionary that contains information on the blast hit
        and alignment. 

        :arg percent_identity: PDB hits with percent sequence identity equal 
            to or higher than this value will be returned.
        :type percent_identity: float, default is 90.0
        :arg percent_coverage: PDB hits with percent coverage of the query 
          sequence equivalent or better will be returned.
        :type percent_coverage: float, default is 70.0
        
        """
        
        query_length = self._record.query_length
        
        pdb_hits = {}
        for alignment in self._record.alignments:
            #Stores information about one hit in the alignments section.
            #######self.sqid = [ 0.0 ]
            #A list because getPercentSequenceIdentity returns a list
            hsp = alignment.hsps[0]
            p_identity = (100.0 * hsp.identities / query_length)
            p_coverage = (100.0 * hsp.align_length / query_length)
            if p_identity > percent_identity and p_coverage > percent_coverage:
                items = alignment.title.split('>gi')
                for item in items:
                    #>gi|1633462|pdb|4AKE|A Chain A, Adenylate Kinase
                    #                        __________TITLE__________
                    title = item[item.find(' ')+1:].encode('utf-8').strip()
                    pdb_id = item.split()[0].split('|')[3].encode('utf-8')
                    pdb_id = pdb_id.lower()
                    # if pdb_id is extracted, but a better alignment exists
                    if (pdb_id not in pdb_hits or 
                        p_identity > pdb_hits[pdb_id]['percent_identity']):
                        hit = {}
                        hit['pdb_id'] = pdb_id
                        hit['chain_id'] = title[6]
                        hit['percent_identity'] = p_identity
                        hit['percent_coverage'] = p_coverage
                        hit['bits'] = hsp.bits
                        hit['score'] = hsp.score
                        hit['expect'] = hsp.expect
                        hit['align_length'] = hsp.align_length
                        hit['identities'] = hsp.identities
                        hit['positives'] = hsp.positives
                        hit['gaps'] = hsp.gaps
                        hit['pdb_title'] = title
                        hit['query'] = hsp.query
                        hit['query_start'] = hsp.query_start
                        hit['query_end'] = hsp.query_end
                        hit['sbjct'] = hsp.sbjct
                        hit['sbjct_start'] = hsp.sbjct_start
                        hit['sbjct_end'] = hsp.sbjct_end
                        hit['match'] = hsp.match
                        pdb_hits[pdb_id] = hit
        LOGGER.info('{0:d} hits were identified.'
                         .format(len(pdb_hits)))
        return pdb_hits
    
    def getBest(self):
        """Returns a dictionary for the hit with highest sequence identity."""
        
        from heapq import heappop, heappush
        
        query_length = self._record.query_length
        
        pdb_hits = {}
        hit_list = []
        for alignment in self._record.alignments:
            hsp = alignment.hsps[0]
            p_identity = (100.0 * hsp.identities / query_length)
            p_coverage = (100.0 * hsp.align_length / query_length)
            items = alignment.title.split('>gi')
            for item in items:
                title = item[item.find(' ')+1:].encode('utf-8').strip()
                pdb_id = item.split()[0].split('|')[3].encode('utf-8')
                pdb_id = pdb_id.lower()
                if (pdb_id not in pdb_hits or 
                    p_identity > pdb_hits[pdb_id]['percent_identity']):
                    hit = {}
                    hit['pdb_id'] = pdb_id
                    hit['chain_id'] = title[6]
                    hit['percent_identity'] = p_identity
                    hit['percent_coverage'] = p_coverage
                    hit['bits'] = hsp.bits
                    hit['score'] = hsp.score
                    hit['expect'] = hsp.expect
                    hit['align_length'] = hsp.align_length
                    hit['identities'] = hsp.identities
                    hit['positives'] = hsp.positives
                    hit['gaps'] = hsp.gaps
                    hit['pdb_title'] = title
                    hit['query'] = hsp.query
                    hit['query_start'] = hsp.query_start
                    hit['query_end'] = hsp.query_end
                    hit['sbjct'] = hsp.sbjct
                    hit['sbjct_start'] = hsp.sbjct_start
                    hit['sbjct_end'] = hsp.sbjct_end
                    hit['match'] = hsp.match
                    pdb_hits[pdb_id] = hit
                    heappush(hit_list, (-p_identity, hit))
        return heappop(hit_list)[1]

def blastPDB(sequence, filename=None, **kwargs):
    """Blast search ProteinDataBank for a given sequence and return results.
        
    :arg sequence: Single-letter code amino acid sequence of the protein.
    :type sequence: str, :class:`Bio.SeqRecord.SeqRecord`, or 
        :class:`Bio.Seq.Seq` 
    
    :arg filename: Provide a *filename* to save the results in XML format. 
    :type filename: str
    
    This method uses :meth:`qdblast` function in :mod:`NCBIWWW` module of 
    Biopython. It will use *blastp* program and search *pdb* database.
    Results are parsed using :meth:`NCBIXML.parse` and passed to a
    :class:`PDBBlastRecord`
    
    User can adjust *hitlist_size* (default is 250) and *expect* (default is 
    1e-10) parameter values. 

    """
    if BioBlast is None: prody.importBioBlast()
    if not isinstance(sequence, str):
        try:
            sequence = sequence.data
        except AttributeError:    
            try:
                sequence = str(sequence.seq)
            except AttributeError:    
                pass
        if not isinstance(sequence, str):
            raise TypeError('sequence must be a string or a Biopython object')
    if not sequence:
        raise ValueError('sequence cannot be an empty string')

    if 'hitlist_size' not in kwargs:
        kwargs['hitlist_size'] = 250
    if 'expect' not in kwargs:
        kwargs['expect'] = 1e-10

    sequence = ''.join(sequence.split())
    LOGGER.info('Blasting ProteinDataBank for "{0:s}...{1:s}"'
                .format(sequence[:10], sequence[-10:]))
    start = time.time()
    results = BioBlast.qblast('blastp', 'pdb', sequence, format_type='XML', 
                              **kwargs)

    LOGGER.info('Blast search completed in {0:.2f}s.'
                .format(time.time()-start))
                
    if filename is not None:
        filename = str(filename)
        if not filename.lower().endswith('.xml'):
                filename += '.xml'        
        out = open(filename, 'w')
        for line in self._results:
            out.write(line)
        out.close()
        LOGGER.info('Results are saved as {0:s}.'.format(filename))    
        results.reset()
    record = BioBlast.parse(results).next()
    return PDBBlastRecord(sequence, record)

_writePDBdoc = """
    :arg atoms: Atomic data container.
    :type atoms: :class:`~prody.atomic.Atomic` 
    
    :arg model: Model index or list of model indices.
    :type model: int, list
        
    If *models* is ``None``, all coordinate sets will be written. Model 
    indices start from 1.
    
    *atoms* instance must at least contain coordinates and atom names data.
    
    """

def writePDBStream(stream, atoms, model=None, sort=False):
    """Write *atoms* in PDB format to a *stream*.
    
    :arg stream: anything that implements the method write() 
        (e.g. file, buffer, stdout)
    
    """
    
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms does not have a valid type')
    if isinstance(atoms, prody.Atom):
        atoms = prody.Selection(atoms.getAtomGroup(), [atoms.getIndex()], 
                                atoms.getActiveCoordsetIndex(), 
                                'index ' + str(atoms.getIndex()))

    if model is None:
        model = np.arange(atoms.getNumOfCoordsets())
    elif isinstance(model, int):
        model = np.array([model], np.int64) -1
    elif isinstance(model, list):
        model = np.array(model, np.int64) -1
    else:
        raise TypeError('model must be an integer or a list of integers')
    if model.min() < 0 or model.max() >= atoms.getNumOfCoordsets():
        raise ValueError('model index or indices is not valid')
        
    n_atoms = atoms.getNumOfAtoms()
    atomnames = atoms.getAtomNames()
    if atomnames is None:
        raise RuntimeError('atom names are not set')
    for i, an in enumerate(atomnames):
        lenan = len(an)
        if lenan < 4:
            atomnames[i] = ' ' + an
        elif lenan > 4:
            atomnames[i] = an[:4]
    resnames = atoms._getResidueNames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms
    resnums = atoms._getResidueNumbers()
    if resnums is None:
        resnums = np.ones(n_atoms, np.int64)
    chainids = atoms._getChainIdentifiers()
    if chainids is None: 
        chainids = np.zeros(n_atoms, '|S1')
    occupancies = atoms._getOccupancies()
    if occupancies is None:
        occupancies = np.zeros(n_atoms, np.float64)
    bfactors = atoms._getTempFactors()
    if bfactors is None:
        bfactors = np.zeros(n_atoms, np.float64)
    icodes = atoms._getInsertionCodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, '|S1')
    hetero = ['ATOM'] * n_atoms 
    heteroflags = atoms._getHeteroFlags()
    if heteroflags is not None:
        hetero = np.array(hetero, '|S6')
        hetero[heteroflags] = 'HETATM'
    elements = atoms._getElementSymbols()
    if elements is None:
        elements = np.zeros(n_atoms, '|S1')
    altlocs = atoms._getAltLocIndicators()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, '|S1')
    segments = atoms._getSegmentNames()
    if segments is None:
        segments = np.zeros(n_atoms, '|S6')
    
    acsi = atoms.getActiveCoordsetIndex()
    multi = False
    if len(model) > 1:
        multi = True
    line = ('{0:6s}{1:5d} {2:4s}{3:1s}' +
            '{4:4s}{5:1s}{6:4d}{7:1s}   ' + 
            '{8:8.3f}{9:8.3f}{10:8.3f}' +
            '{11:6.2f}{12:6.2f}      ' +
            '{13:4s}{14:2s}\n')
    for m in model:
        if multi:
            stream.write('MODEL{0:9d}\n'.format(m+1))
        atoms.setActiveCoordsetIndex(m)
        coords = atoms._getCoordinates()
        for i, xyz in enumerate(coords):
            stream.write(line.format(hetero[i], i+1, atomnames[i], altlocs[i], 
                                     resnames[i], chainids[i], int(resnums[i]), 
                                     icodes[i], 
                                     xyz[0], xyz[1], xyz[2], 
                                     occupancies[i], bfactors[i],  
                                     segments[i], elements[i].rjust(2)))
        if multi:
            stream.write('ENDMDL\n')
            altlocs = np.zeros(n_atoms, '|S1')
    atoms.setActiveCoordsetIndex(acsi)

writePDBStream.__doc__ += _writePDBdoc

def writePDB(filename, atoms, model=None, sort=None):
    """Write *atoms* in PDB format to a file with name *filename*.
    
    Returns *filename* if file is succesfully written. 
   
    """
    
    out = open(filename, 'w')
    writePDBStream(out, atoms, model, sort)
    out.close()
    return filename

writePDB.__doc__ += _writePDBdoc

mapHelix = {
1: 'H', # 4-turn helix (alpha helix)
2: '', # other helix, Right-handed omega
3: 'I', # 5-turn helix (pi helix)
4: '', # other helix, Right-handed gamma
5: 'G', # 3-turn helix (3-10 helix)
6: '', # Left-handed alpha
7: '', # Left-handed omega
8: '', # Left-handed gamma
9: '', # 2 - 7 ribbon/helix
10: '' # Polyproline
}

def assignSecondaryStructure(header, atoms, coil=False):
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
       Secondary structures are assigned to all atoms in a residue. Amino acid
       residues without any secondary structure assignments in the header 
       section will be assigned coil (C) conformation. This can be prevented
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

    ssa = atoms.getSecondaryStrs()
    if ssa is None:
        if isinstance(atoms, prody.AtomGroup):
            ag = atoms
        else:
            ag = atoms.getAtomGroup()
        ag.setSecondaryStrs(np.zeros(ag.getNumOfAtoms(), 
                            ATOMIC_DATA_FIELDS['secondary'].dtype))
    atoms.select('protein').setSecondaryStrs('C')
    hierview = atoms.getHierView()
    count = 0
    for key, value in helix.iteritems():
        res = hierview.getResidue(*key)
        if res is None:
            continue
        res.setSecondaryStrs(mapHelix[value[0]])
        count += 1
    for key, res in sheet.iteritems():
        res = hierview.getResidue(*key)
        if res is None:
            continue
        res.setSecondaryStrs('E')
        count += 1
    LOGGER.info('Secondary structures were assigned to {0:d} residues.'
                .format(count))
    return atoms
            

def fetchLigandData(cci, save=False, folder='.'):
    """Fetch ligand data from Ligand Expo (http://ligand-expo.rcsb.org/).

    .. versionadded:: 0.8.1
    
    .. versionchanged:: 0.8.2
       URL of the XML file is returned in the dictionary with key ``url``.
    
    *cci* may be 3-letter chemical component identifier or a valid XML 
    filename. If ``aave=True`` is passed, XML file will be saved in the 
    specified *folder*. 
    
    This function is compatible with PDBx/PDBML v 4.0.
    
    Returns a dictionary that contains ligand data. Ligand atom data with 
    *model* and *ideal* coordinate sets are also stored in this dictionary.
    Note that this dictionary will contain data that is present in the XML
    file and all Ligand Expo XML files do not contain every possible data
    field. So, it may be better if you use :meth:`dict.get` instead of
    indexing the dictionary, e.g. to retrieve formula weight (or relative
    molar mass) of the chemical component use ``data.get('formula_weight')``
    instead of ``data['formula_weight']`` to avoid exceptions when this
    data field is not found in the XML file.
    
    Following example downloads data for ligand STI (a.k.a. Gleevec and
    Imatinib) and calculates RMSD between model (X-ray) and ideal (energy 
    minimized) coordinate sets.
    
    >>> ligand_data = fetchLigandData('STI')
    >>> ligand_data['model_coordinates_db_code'] 
    '1IEP'
    
    Model (X-ray) coordinates are from structure **1IEP**.
    
    >>> ligand_model = ligand_data['model'] 
    >>> ligand_ideal = ligand_data['ideal']
    >>> transformation = superpose(ligand_ideal.noh, ligand_model.noh)
    >>> print( calcRMSD(ligand_ideal.noh, ligand_model.noh).round(2) )
    2.27
    """

    global urllib2, ET
    if not isinstance(cci, str):
        raise TypeError('cci must be a string')
    elif os.path.isfile(cci):
        inp = open(cci)
        xml = inp.read()
        inp.close()
    elif len(cci) > 4 or not cci.isalnum(): 
        raise ValueError('cci must be 3-letters long and alphanumeric or '
                         'a valid filename')
    else:
        #'http://www.pdb.org/pdb/files/ligand/{0:s}.xml'
        url = ('http://ligand-expo.rcsb.org/reports/{0[0]:s}/{0:s}/{0:s}.xml'
               .format(cci.upper()))
        if urllib2 is None:
            import urllib2
            prody.proteins.urllib2 = urllib2
        try:
            inp = urllib2.urlopen(url)
        except urllib2.HTTPError:
            raise IOError('XML file for ligand {0:s} is not found'.format(cci))
        else:
            xml = inp.read()
            inp.close()
        if save:
            out = open(os.path.join(folder, cci+'.xml'), 'w')
            out.write(xml)
            out.close()

    if ET is None:
        import xml.etree.cElementTree as ET
        prody.proteins.ET = ET

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
    ordinals = np.zeros(n_atoms, np.int64)
    
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
    model.setCoordinates(model_coords)
    model.setAtomNames(atomnames)
    model.setResidueNames(resnames)
    model.setResidueNames(resnums)
    model.setElementSymbols(elements)
    model.setCharges(charges)
    model.setAttribute('leaving_atom_flags', leaving_atom_flags)
    model.setAttribute('aromatic_flags', aromatic_flags)
    model.setAttribute('stereo_configs', stereo_configs)
    model.setAttribute('ordinals', ordinals)
    model.setAttribute('alternate_atomnames', alternate_atomnames)
    dict_['model'] = model
    ideal = model.copy()
    ideal.setName(cci + ' ideal')
    ideal.setCoordinates(ideal_coords)
    dict_['ideal'] = ideal

    return dict_      

def applyBiomolecularTransformations(header, atoms, biomol=None):
    """Return *atoms* after applying biomolecular transformations from *header*
    dictionary.
    
    .. versionchanged:: 0.7.1
       Biomolecular transformations are applied to all coordinate sets in the
       molecule.
    
    .. versionchanged:: 0.8.2
       Raises :class:`ValueError` when *header* does not contain 
       biomolecular transformations.

    Some PDB files contain transformations for more than 1 biomolecules. A 
    specific set of transformations can be choosen using *biomol* argument.
    Transformation sets are identified by integer numbers, e.g. 1, 2, ...
    
    If multiple biomolecular transformations are provided in the *header*
    dictionary, biomolecules will be returned as 
    :class:`~prody.atomic.AtomGroup` instances in a :func:`list`.  

    If the resulting biomolecule has more than 26 chains, the molecular 
    assembly will be split into multiple :class:`~prody.atomic.AtomGroup`
    instances each containing at most 26 chains. These 
    :class:`~prody.atomic.AtomGroup` instances will be returned in a tuple.
   
    Note that atoms in biomolecules are ordered according to chain identifiers.
   
    """

    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must be an Atomic instance')
    biomt = header.get('biomolecular_transformations')
    if not isinstance(biomt, dict) or len(biomt) == 0:
        raise ValueError('header does not contain biomolecular' 
                         'transformations')
    
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
                
                for acsi in range(newag.getNumOfCoordsets()):
                    newag.setActiveCoordsetIndex(acsi)
                    newag = t.apply(newag)
                newag.setActiveCoordsetIndex(0)
                ags.append(newag)
        if ags:
            # Handles the case when there is more atom groups than the number
            # of chain identifiers
            if len(chids_used) != len(set(chids_used)):
                for newag in ags:
                    newag.select('all').setChainIdentifiers(chids.pop(0))
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
                    newag.setName('{0:s} biomolecule {1:s} part {2:d}'
                                  .format(atoms.getName(), i, k+1))
                else:
                    newag.setName('{0:s} biomolecule {1:s}'
                                  .format(atoms.getName(), i))
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

def execDSSP(pdb, outputname=None, outputdir=None):
    """Execute DSSP for given *pdb*. *pdb* can be a PDB identifier or a PDB 
    file path. If *pdb* is a compressed file, it will be decompressed using
    Python :mod:`gzip` library. When no *outputname* is given, output name 
    will be :file:`pdbidentifier.dssp`. :file:`.dssp` extension will be 
    appended automatically to *outputname*. If :file:`outputdir` is given,
    DSSP output and uncompressed PDB file will be written into this folder.
    Upon successful execution of :command:`dssp pdb > out` command, output
    filename is returned. 
    
    For more information on DSSP see http://swift.cmbi.ru.nl/gv/dssp/.
    If you benefited from DSSP, please consider citing [WK83]_.
    
    .. versionadded:: 0.8"""
    
    dssp = prody.which('dssp')
    if dssp is None:
        raise EnvironmentError('command not found: dssp executable is not '
                               'found in one of system paths')
    
    if not os.path.isfile(pdb):
        pdb = fetchPDB(pdb, compressed=False)
    if pdb is None:
        raise ValueError('pdb is not a valid PDB identifier or filename')
    if os.path.splitext(pdb)[1] == '.gz':
        if outputdir:
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
        
    status = os.system('{0:s} {1:s} > {2:s}'.format(dssp, pdb, out))
    if status == 0:
        return out
    
def parseDSSP(dssp, ag, parseall=False):
    """Parse DSSP data from file *dssp* into :class:`~prody.atomic.AtomGroup`
    instance *ag*. DSSP output file must be in the new format used from July 
    1995 and onwards. When *dssp* file is parsed, following attributes are 
    added to *ag*:
        
    * *dssp_resnum*: DSSP's sequential residue number, starting at the first 
      residue actually in the data set and including chain breaks; this number 
      is used to refer to residues throughout.
    
    * *dssp_acc*: number of water molecules in contact with this residue \*10. 
      or residue water exposed surface in Angstrom^2.
      
    * *dssp_kappa*: virtual bond angle (bend angle) defined by the three Cα 
      atoms of residues I-2,I,I+2. Used to define bend (structure code 'S').
        
    * *dssp_alpha*: virtual torsion angle (dihedral angle) defined by the four 
      Cα atoms of residues I-1,I,I+1,I+2.Used to define chirality (structure 
      code '+' or '-').
        
    * *dssp_phi* and *dssp_psi*: IUPAC peptide backbone torsion angles 

    The following attributes are parsed when ``parseall=True`` is passed: 

    * *dssp_bp1*, *dssp_bp2*, and *dssp_sheet_label*: residue number of first 
      and second bridge partner followed by one letter sheet label
      
    * *dssp_tco*: cosine of angle between C=O of residue I and C=O of residue 
      I-1. For α-helices, TCO is near +1, for β-sheets TCO is near -1. Not 
      used for structure definition.
      
    * *dssp_NH_O_1_index*, *dssp_NH_O_1_energy*, etc.: hydrogen bonds; e.g. 
      -3,-1.4 means: if this residue is residue i then N-H of I is h-bonded to 
      C=O of I-3 with an electrostatic H-bond energy of -1.4 kcal/mol. There 
      are two columns for each type of H-bond, to allow for bifurcated H-bonds.
        
    See http://swift.cmbi.ru.nl/gv/dssp/DSSP_3.html for details.
    
    .. versionadded:: 0.8"""
    
    if not os.path.isfile(dssp):
        raise IOError('{0:s} is not a valid file path'.format(dssp))
    if not isinstance(ag, prody.AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')
        
    dssp = open(dssp)
    
    n_atoms = ag.getNumOfAtoms()
    NUMBER = np.zeros(n_atoms, np.int64)
    SHEETLABEL = np.zeros(n_atoms, '|S1')
    ACC = np.zeros(n_atoms, np.float64)
    KAPPA = np.zeros(n_atoms, np.float64)
    ALPHA = np.zeros(n_atoms, np.float64)
    PHI = np.zeros(n_atoms, np.float64)
    PSI = np.zeros(n_atoms, np.float64)

    if parseall:
        BP1 = np.zeros(n_atoms, np.int64)
        BP2 = np.zeros(n_atoms, np.int64)
        NH_O_1 = np.zeros(n_atoms, np.int64)
        NH_O_1_nrg = np.zeros(n_atoms, np.float64)
        O_HN_1 = np.zeros(n_atoms, np.int64)
        O_HN_1_nrg = np.zeros(n_atoms, np.float64)
        NH_O_2 = np.zeros(n_atoms, np.int64)
        NH_O_2_nrg = np.zeros(n_atoms, np.float64)
        O_HN_2 = np.zeros(n_atoms, np.int64)
        O_HN_2_nrg = np.zeros(n_atoms, np.float64)
        TCO = np.zeros(n_atoms, np.float64)

    ag.setSecondaryStrs(np.zeros(n_atoms))
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
        res.setSecondaryStrs(line[16].strip())
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
            NH_O_1[indices] = int(line[40:45])
            NH_O_1_nrg[indices] = float(line[47:51]) 
            O_HN_1[indices] = int(line[51:56])
            O_HN_1_nrg[indices] = float(line[57:61])
            NH_O_2[indices] = int(line[62:67])
            NH_O_2_nrg[indices] = float(line[68:72])
            O_HN_2[indices] = int(line[73:78])
            O_HN_2_nrg[indices] = float(line[79:83])
            TCO[indices] = float(line[85:91])
    
    ag.setAttribute('dssp_resnum', NUMBER)
    ag.setAttribute('dssp_sheet_label', SHEETLABEL)
    ag.setAttribute('dssp_acc', ACC)
    ag.setAttribute('dssp_kappa', KAPPA)
    ag.setAttribute('dssp_alpha', ALPHA)
    ag.setAttribute('dssp_phi', PHI)
    ag.setAttribute('dssp_psi', PSI)

    if parseall:
        ag.setAttribute('dssp_bp1', BP1)
        ag.setAttribute('dssp_bp2', BP2)
        ag.setAttribute('dssp_NH_O_1_index', NH_O_1)
        ag.setAttribute('dssp_NH_O_1_energy', NH_O_1_nrg)
        ag.setAttribute('dssp_O_NH_1_index', O_HN_1)
        ag.setAttribute('dssp_O_NH_1_energy', O_HN_1_nrg)    
        ag.setAttribute('dssp_NH_O_2_index', NH_O_2)
        ag.setAttribute('dssp_NH_O_2_energy', NH_O_2_nrg)
        ag.setAttribute('dssp_O_NH_2_index', O_HN_2)
        ag.setAttribute('dssp_O_NH_2_energy', O_HN_2_nrg)
        ag.setAttribute('dssp_tco', TCO)
    return ag

def performDSSP(pdb, parseall=False):
    """Perform DSSP calculations and parse results. DSSP data is returned 
    in an :class:`~prody.atomic.AtomGroup` instance. See also :func:`execDSSP` 
    and :func:`parseDSSP`.
    
    .. versionadded:: 0.8"""
    
    pdb = fetchPDB(pdb, compressed=False)
    return parseDSSP(execDSSP(pdb), parsePDB(pdb), parseall)
    
def execSTRIDE(pdb, outputname=None, outputdir=None):
    """Execute STRIDE program for given *pdb*. *pdb* can be an identifier 
    or a PDB file path. If *pdb* is a compressed file, it will be decompressed 
    using Python :mod:`gzip` library. When no *outputname* is given, output 
    name will be :file:`pdbidentifier.stride`. :file:`.stride` extension will 
    be appended automatically to *outputname*. If :file:`outputdir` is given, 
    STRIDE output and uncompressed PDB file will be written into this folder.
    Upon successful execution of :command:`stride pdb > out` command, output
    filename is returned. 
    
    For more information on STRIDE see http://webclu.bio.wzw.tum.de/stride/.
    If you benefited from STRIDE, please consider citing [DF95]_.
    
    .. versionadded:: 0.8"""
    
    stride = prody.which('stride')
    if stride is None:
        raise EnvironmentError('command not found: stride executable is not '
                               'found in one of system paths')
    
    if not os.path.isfile(pdb):
        pdb = fetchPDB(pdb, compressed=False)
    if pdb is None:
        raise ValueError('pdb is not a valid PDB identifier or filename')
    if os.path.splitext(pdb)[1] == '.gz':
        if outputdir:
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
    :class:`~prody.atomic.AtomGroup` instance *ag*. STRIDE output file must 
    be in the new format used from July 1995 and onwards. When *stride* file 
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
    
    n_atoms = ag.getNumOfAtoms()
    NUMBER = np.zeros(n_atoms, np.int64)
    AREA = np.zeros(n_atoms, np.float64)
    PHI = np.zeros(n_atoms, np.float64)
    PSI = np.zeros(n_atoms, np.float64)

    ag.setSecondaryStrs(np.zeros(n_atoms))
    for line in stride:
        if not line.startswith('ASG '):
            continue
        res = ag[(line[9], int(line[10:15]), line[15].strip())]
        if res is None:
            continue
        indices = res.getIndices()
        res.setSecondaryStrs(line[24].strip())
        NUMBER[indices] = int(line[16:20])
        PHI[indices] = float(line[42:49])
        PSI[indices] = float(line[52:59])
        AREA[indices] = float(line[64:69])
    ag.setAttribute('stride_resnum', NUMBER)
    ag.setAttribute('stride_phi', PHI)
    ag.setAttribute('stride_psi', PSI)
    ag.setAttribute('stride_area', AREA)
    return ag

def performSTRIDE(pdb):
    """Perform STRIDE calculations and parse results. STRIDE data is 
    returned in an :class:`~prody.atomic.AtomGroup` instance. See also 
    :func:`execSTRIDE` and :func:`parseSTRIDE`.
    
    .. versionadded:: 0.8"""
    
    pdb = fetchPDB(pdb, compressed=False)
    return parseSTRIDE(execSTRIDE(pdb), parsePDB(pdb))




def fetchPDBClusters():
    """Downloads PDB sequence clusters. PDB sequence clusters are results of 
    the weekly clustering of protein chains in the PDB generated by blastclust. 
    They are available at FTP site: ftp://resources.rcsb.org/sequence/clusters/
    
    This function will download about 10 Mb of data and save it after 
    compressing in your home directory in :file:`.prody/pdbclusters`.
    Compressed files will be less than 4 Mb in size. Cluster data can 
    be loaded using :func:`loadPDBClusters` function and be accessed 
    using :func:`getPDBCluster`.
    
    .. versionadded:: 0.8.2"""
    
    global urllib2
    if urllib2 is None:
        import urllib2
        prody.proteins.urllib2 = urllib2
    PDB_CLUSTERS_PATH = os.path.join(prody.getPackagePath(), 'pdbclusters')
    if not os.path.isdir(PDB_CLUSTERS_PATH):
        os.mkdir(PDB_CLUSTERS_PATH)
    progress = prody.ProDyProgress(len(PDB_CLUSTERS))
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
            out = gzip.open(os.path.join(PDB_CLUSTERS_PATH, 
                                         filename+'.gz'), 'w')
            out.write(inp.read())
            inp.close()
            out.close()
        progress.report(i)
    progress.clear()

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
            diff = (time.time() - os.path.getmtime(filename)) / 604800.
            if diff > 1.:
                LOGGER.warning('PDB sequence clusters are {0:.1f} week old, '
                               'call `fetchPDBClusters` to receive updates.'
                               .format(diff))
                PDB_CLUSTERS_UPDATE_WARNING = False
        inp = gzip.open(filename)
        
        clusters = {}
        for line in inp.readlines():
            cluster = []
            for item in line.split():
                item = tuple(item.split('_'))
                cluster.append(item) 
                clusters[item] = cluster
        inp.close()

        PDB_CLUSTERS[sqid] = clusters

def getPDBCluster(pdb, ch, sqid=95):
    """Return the PDB sequence cluster for chain *ch* in structure *pdb*
    that chains sharing sequence identity *sqid* or more. PDB sequence cluster
    will be returned in the form of a list of tuples, e.g. 
    ``[('1XXX', 'A'), ('2YYY', 'A')]``. Note that PDB clusters chains, so
    the same PDB identifier may appear twice in the same cluster if the 
    corresponding chain is present in the structure twice.    
    
    Before this function is used, :func:`fetchPDBClusters` needs to be called. 
    This function will load the PDB sequence clusters for *sqid* automatically 
    using :func:`loadPDBClusters`.
    
    .. versionadded:: 0.8.2"""

    assert isinstance(pdb, str), 'pdb must be a string'
    assert isinstance(ch, str), 'pdb must be a string'
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
    return list(clusters[(pdb.upper(), ch.upper())])
