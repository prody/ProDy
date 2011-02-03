# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan
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
and write PDB files, and also to blast search `ProteinDataBank <http://wwpdb.org>`_.

Classes
-------

  * :class:`PDBFetcher`
  * :class:`RCSB_PDBFetcher`
  * :class:`PDBlastRecord`

Functions
---------

  * :func:`applyBiomolecularTransformations`
  * :func:`assignSecondaryStructure`
  * :func:`fetchPDB` 
  * :func:`blastPDB`
  * :func:`parsePDB`
  * :func:`parsePDBStream`
  * :func:`writePDB`
  * :func:`writePDBStream`
   
   
.. doctest::
    :hide:
        
    >>> import numpy as np
    >>> allmodels = parsePDB('2k39', subset='ca')
    >>> model10 = parsePDB('2k39', subset='ca', model=10)
    >>> np.all(allmodels.getCoordsets(9) == model10.getCoordinates())
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import gzip
import os.path
import time
import os
from glob import glob
from collections import defaultdict

import numpy as np

import prody
from prody import ProDyLogger as LOGGER
from prody.atomic import *


__all__ = ['PDBlastRecord', 'PDBFetcher', 'RCSB_PDBFetcher', 
           'assignSecondaryStructure',
           'applyBiomolecularTransformations',
           'blastPDB', 
           'fetchPDB',
           'parsePDBStream', 'parsePDB', 
           'writePDBStream', 'writePDB',
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


class PDBFetcher(object):
    """Base class for PDB fetcher classes."""
    
    
    @staticmethod
    def fetch(pdb, folder='.'):
        """Fetch and save PDB file(s) into *folder*."""
        pass


class RCSB_PDBFetcher(PDBFetcher):
    """A class to fetch PDB files from FTP server of RCSB."""
    
    @staticmethod
    def fetch(pdb, folder='.'):
        """Fetch pdb file(s) from RCSB PDB - ftp.wwpdb.org.
        
        Downloads PDB files by establishing an FTP connection to *ftp.wwpdb.org*.

        Downloaded files will be saved in *folder*. FTP server provides 
        gunzipped PDB files. If *folder* already contains a PDB file 
        matching an identifier, a file will not be downloaded, but will be 
        contained in the returned list.
        
        *pdb* may be a list of PDB identifiers or an identifier string.
            
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
            pdbfnmap[os.path.split(pdbfn)[1].split('.')[0].lower()] = pdbfn
        for pdbfn in glob(os.path.join(folder, '*.PDB*')): 
            pdbfnmap[os.path.split(pdbfn)[1].split('.')[0].lower()] = pdbfn

        for i, pdbid in enumerate(identifiers):
            pdbid = pdbid.strip().lower()
            if not (isinstance(pdbid, str) and len(pdbid) == 4 and pdbid.isalnum()):
                LOGGER.debug('{0:s} is not a valid identifier.'.format(pdbid))
                filenames.append(None)
                failure += 1 
                continue
            identifiers[i] = pdbid
            try:
                fn = pdbfnmap[pdbid]
            except KeyError:
                pass
            else:
                fn = os.path.relpath(fn)
                filenames.append(fn)
                LOGGER.debug('{0:s} ({1:s}) is found in the target directory.'
                             .format(pdbid, fn))
                exists += 1
                continue
            filenames.append(pdbid)
            download = True
        if download:
            from ftplib import FTP
            try:
                ftp = FTP('ftp.wwpdb.org')
            except Exception as error:
                raise type(error)('FTP connection problem, potential reason: '
                                  'no internet connectivity')
            else:
                ftp.login('')
                for i, pdbid in enumerate(identifiers):
                    if pdbid != filenames[i]:
                        continue
                    filename = os.path.join(folder, pdbid + '.pdb.gz')
                    pdbfile = open(filename, 'w+b')
                    try:
                        ftp.cwd('/pub/pdb/data/structures/divided/pdb/{0:s}'
                                .format(pdbid[1:3]))
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
 
DEFAULT_PDBFetcher = RCSB_PDBFetcher
 
def fetchPDB(pdb, folder='.', fetcher=None):
    """Fetch PDB files using default PDB fetcher class.
    
    If *fetcher* is ``None``, default fetcher will be used. Default fetcher 
    (:data:`~prody.proteins.DEFAULT_PDBFetcher`) is :class:`RCSB_PDBFetcher`.
    
    """
    if fetcher is None:
        global DEFAULT_PDBFetcher
        return DEFAULT_PDBFetcher.fetch(pdb, folder)
    else:
        return fetcher.fetch(pdb, folder)

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

    If ``model=0`` and ``header=True``, return header 
    dictionary only."""
    
_PDBSubsets = ['ca', 'calpha', 'bb', 'backbone']

def parsePDB(pdb, model=None, header=False, chain=None, subset=None, 
             altloc='A', name=None):
    """Return an :class:`~prody.atomic.AtomGroup` and/or 
    dictionary containing header data parsed from a stream of PDB lines. 
    
    This function extends :func:`parsePDBStream`.
    
    |example| See :ref:`parsepdb` for a detailed example.
    
    :arg pdb: A valid PDB identifier or filename.  
        If needed, PDB files are downloaded using :func:`fetchPDB()` function.  
        
    """
    if not os.path.isfile(pdb):
        if len(pdb) == 4 and pdb.isalnum():
            download = RCSB_PDBFetcher.fetch(pdb)
            if not download:
                raise PDBParserError('PDB file for {0:s} could not be '
                                   'downloaded.'.format(pdb))
            pdb = download
        else:
            raise PDBParserError('{0:s} is not a valid filename or a valid '
                               'PDB identifier.'.format(pdb))
    name, temp = os.path.splitext(os.path.split(pdb)[1])
    if temp == '.gz':
        name, temp = os.path.splitext(name)
        pdb = gzip.open(pdb)
    else:
        pdb = open(pdb)
    name = name.lower()
    result = parsePDBStream(pdb, model, header, chain, subset, altloc, name)
    return result

parsePDB.__doc__ += _parsePDBdoc
    
def parsePDBStream(stream, model=None, header=False, chain=None, subset=None, 
                   altloc='A', name=None):
    """Return an :class:`~prody.atomic.AtomGroup` and/or 
    dictionary containing header data parsed from a stream of PDB lines. 
    
    :arg stream: Anything that implements the method readlines() 
        (e.g. :class:`file`, buffer, stdin).

    """
    if model is not None:
        if isinstance(model, int):
            if model < 0:
                raise PDBParserError('model must be greater than 0')
        else:
            raise TypeError('model must be an integer, {0:s} is invalid'
                            .format(str(model)))
    if subset is not None: 
        if not isinstance(subset, str):
            raise TypeError('subset must be a string')
        elif subset not in _PDBSubsets:
            raise ValueError('"{0:s}" is not a valid subset'.format(subset))

    lines = stream.readlines()

    split = 0
    hd = None
    ag = None
    if header:
        hd, split = _getHeaderDict(lines)
    
    if model != 0:
        start = time.time()
        ag = _getAtomGroup(lines, split, model, chain, subset, altloc)
        if name is None:
            ag.setName('unknown')
        else:
            ag.setName(name)
        LOGGER.info('{0:d} atoms and {1:d} coordinate sets were '
                    'parsed in {2:.2f}s.'.format(ag._n_atoms, ag._n_coordsets, 
                                                 time.time()-start))
    if ag is not None and hd is not None:
        return ag, hd
    elif ag is not None:
        return ag
    else:
        return hd

parsePDBStream.__doc__ += _parsePDBdoc

def _getAtomGroup(lines, split, model, chain, subset, altloc_torf):
    """Return an AtomGroup. See also :func:`parsePDBStream()`.
    
    :arg lines: Lines from a PDB file.
    :arg split: Starting index for lines related to coordinate data.
    """
    asize = len(lines) - split
    alength = 5000
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
    asize = 2000
    
    onlycoords = False
        
    if subset is not None:
        subset = subset.lower()
        if subset in ('calpha', 'ca'):
            subset = ('CA',)
        elif subset in ('backbone', 'bb'):
            subset = ('CA', 'C', 'N', 'O')
        only_subset = True
    else:    
        only_subset = False
    if isinstance(chain, str):
        only_chains = True
    else:
        only_chains = False
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
        LOGGER.info('Parsing alternate locations {0:s}.'.format(altloc_torf))
        which_altlocs = ' ' + ''.join(altloc_torf.split())
        altloc_torf = False
    else:
        which_altlocs = ' A'
        
    acount = 0
    is_anisou = False
    is_siguij = False
    is_scndry = False
    altloc = defaultdict(list)
    n_atoms = 0
    i = start
    while i < stop:
        line = lines[i]
        startswith = line[0:6]
        
        if startswith == 'ATOM  ' or startswith == 'HETATM':
            atomname = line[12:16].strip()
            if only_subset:
                if not atomname in subset: 
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
            except Exception:
                if acount >= n_atoms:
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
            
            altlocs[acount] = alt 
            atomnames[acount] = atomname
            
            # resname = line[17:21], but some are 4 chars long
            resnames[acount] = line[17:21].strip()
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
            if acount >= alength:
                alength += asize
                coordinates = np.concatenate(
                    (coordinates, np.zeros((asize, 3), np.float64)))
                atomnames = np.concatenate(
                    (atomnames, np.zeros(asize, ATOMIC_DATA_FIELDS['name'].dtype)))
                resnames = np.concatenate( 
                    (resnames, np.zeros(asize, ATOMIC_DATA_FIELDS['resname'].dtype)))
                resnums = np.concatenate( 
                    (resnums, np.zeros(asize, ATOMIC_DATA_FIELDS['resnum'].dtype)))
                chainids = np.concatenate( 
                    (chainids, np.zeros(asize, ATOMIC_DATA_FIELDS['chain'].dtype)))
                bfactors = np.concatenate( 
                    (bfactors, np.zeros(asize, ATOMIC_DATA_FIELDS['beta'].dtype)))
                occupancies = np.concatenate( 
                    (occupancies, np.zeros(asize, ATOMIC_DATA_FIELDS['occupancy'].dtype)))
                hetero = np.concatenate( 
                    (hetero, np.zeros(asize, ATOMIC_DATA_FIELDS['hetero'].dtype)))
                altlocs = np.concatenate( 
                    (altlocs, np.zeros(asize, ATOMIC_DATA_FIELDS['altloc'].dtype)))
                segnames = np.concatenate( 
                    (segnames, np.zeros(asize, ATOMIC_DATA_FIELDS['segment'].dtype)))
                elements = np.concatenate(
                    (elements, np.zeros(asize, ATOMIC_DATA_FIELDS['element'].dtype)))
                secondary = np.concatenate(
                    (secondary, np.zeros(asize, ATOMIC_DATA_FIELDS['secondary'].dtype)))
                anisou = np.concatenate(
                    (anisou, np.zeros((asize, 6), ATOMIC_DATA_FIELDS['anisou'].dtype)))
                siguij = np.concatenate(
                    (siguij, np.zeros((asize, 6), ATOMIC_DATA_FIELDS['siguij'].dtype)))
                icodes = np.concatenate(
                    (icodes, np.zeros(asize, ATOMIC_DATA_FIELDS['icode'].dtype)))
        elif startswith == 'END   ' or startswith == 'CONECT':
            i += 1
            break
        elif startswith == 'ENDMDL':
            if model is not None:
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
                atomgroup = prody.AtomGroup('')
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
                    _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames)
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
        atomgroup = prody.AtomGroup('')
        atomgroup.setCoordinates(coordinates[:acount])
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
                    LOGGER.warning('failed to parse alternate location {0:s} at line {1:d}, residue does not exist as altloc A'.format(key, i+1))
                    continue
                rn = resnames[ids[0]]
                ans = atomnames[ids]
                indices[(ach, ari)] = (rn, ids, ans)
            if rn != arn:
                LOGGER.warning('failed to parse alternate location {0:s} at line {1:d}, residue names do not match (expected {2:s}, parsed {3:s})'.format(key, i+1, rn, arn))
                continue
            index = ids[(ans == aan).nonzero()[0]]
            if len(index) != 1:
                LOGGER.warning('failed to parse alternate location {0:s} at line {1:d}, could not identify matching atom ({2:s} not found in the residue)'.format(key, i+1, aan))
                continue
            try:
                xyz[index[0], 0] = float(line[30:38])
                xyz[index[0], 1] = float(line[38:46])
                xyz[index[0], 2] = float(line[46:54])
            except:
                LOGGER.warning('failed to parse alternate location {0:s} at line {1:d}, could not read coordinates'.format(key, i+1))
                continue
            success += 1
            #except Exception as exception:
            #    print i, line
            #    print exception
            #-->
        LOGGER.info('{0:d} out of {1:d} alternate location {2:s} lines were parsed successfully.'.format(success, len(lines), key))
        if success > 0:
            LOGGER.info('Alternate location {0:s} is appended as a coordinate set to the atom group.'.format(key, atomgroup.getName()))
            atomgroup.addCoordset(xyz)
    
def _getHeaderDict(lines):
    """Return header data in a dictionary."""
    header = {}
    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    helix = {}
    sheet = {}
    biomolecule = defaultdict(list)
    applyToChains = (' ')
    for i in xrange(len(lines)):        
        line = lines[i]
        startswith = line[0:6]
        if startswith == 'ATOM  ' or startswith == 'HETATM' or \
           startswith == 'MODEL ':
            if biomolecule:
                header['biomolecular_transformations'] = biomolecule
            if helix:
                header['helix'] = helix
            if sheet:
                header['sheet'] = sheet
            return header, i
        
        if startswith == 'HELIX ':
            chid = line[19]
            # helix class, serial number, identifier
            value = (int(line[38:40]), int(line[7:10]), line[11:14].strip())
            
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
            chid = line[21]
            # helix class, serial number, identifier
            value = (int(line[38:40]), int(line[7:10]), line[11:14].strip())
            
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
    if helix:
        header['helix'] = helix
    if sheet:
        header['sheet'] = sheet
    return header, i

class PDBlastRecord(object):

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
        """Return the sequence that was used in the search."""
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

        :arg percent_identity: PDB hits with percent identity equivalent or 
          better than this value will be returned.
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
        
    Sequence must be a protein sequence in one letter amino acid code.
    
    This method uses :meth:`qdblast` function in :mod:`NCBIWWW` module of 
    Biopython. It will use *blastp* program and search *pdb* database.
    Results are parsed using :meth:`NCBIXML.parse` and passed to a
    :class:`PDBlastRecord`
    
    User can adjust *hitlist_size* (default is 250) and *expect* (default is 
    1e-10) parameter values. 


    """
    if prody.NCBIWWW is None: prody.importBioBlast()
    if not isinstance(sequence, str):
        raise TypeError('sequence must be a string')
    elif not sequence:
        raise ValueError('sequence cannot be an empty string')

    if 'hitlist_size' not in kwargs:
        kwargs['hitlist_size'] = 250
    if 'expect' not in kwargs:
        kwargs['expect'] = 1e-10

    sequence = ''.join(sequence.split())
    LOGGER.info('Blasting ProteinDataBank for "{0:s}...{1:s}"'
                .format(sequence[:10], sequence[-10:]))
    start = time.time()
    results = prody.NCBIWWW.qblast('blastp', 'pdb', sequence, format_type='XML', **kwargs)

    LOGGER.info('Blast search completed in {0:.2f}s.'
                .format(time.time()-start))
                
    
    if isinstance(filename, str):
        if not filename.lower().endswith('.xml'):
                filename += '.xml'
        else: 
            raise TypeError('filename must be a string')
        
        if compressed: 
            filename += '.gz'
            out = gzip.open(filename, 'w')
        else: 
            out = open(filename, 'w')
        for line in self._results:
            out.write(line)
        out.close()
        LOGGER.info('Results are saved as {0:s}.'.format(filename))    
    
    results.reset()
    record = prody.NCBIXML.parse(results).next()
    return PDBlastRecord(sequence, record)

def writePDBStream(stream, atoms, model=None, sort=False):
    """Write *atoms* in PDB format to a *stream*.
    
    :arg stream: anything that implements the method write() (e.g. file, buffer, stdout)
    
    :arg atoms: Atomic data container.
    :type atoms: :class:`~prody.atomic.Atomic` 
    
    :arg model: Model index or list of model indices.
    :type model: int, list
        
    :arg sort: if True, atoms will be sorted Chain, ResidueNumber, AtomName (not working yet).
    
    If *models* = ``None``, all coordinate sets will be outputted. Model indices
    start from 1.
    
    *atoms* instance must at least contain coordinates and atom names data.
    
    """
    if not isinstance(atoms, (prody.AtomGroup, prody.AtomSubset, prody.AtomMap)):
        raise TypeError('atoms does not have a valid type')
    if isinstance(atoms, prody.Atom):
        atoms = prody.Selection(atoms.getAtomGroup(), [atoms.getIndex()], 
                                atoms.getActiveCoordsetIndex(), 'index ' + str(atoms.getIndex()))

    if model is None:
        model = range(atoms.getNumOfCoordsets())
    elif isinstance(model, int):
        model = np.array([model]) -1
    elif isinstance(models, int):
        model = np.array(model) -1
    else:
        raise TypeError('models must be an integer or a list of integers')


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
    resnames = atoms.getResidueNames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms
    resnums = atoms.getResidueNumbers()
    if resnums is None:
        resnums = np.ones(n_atoms, np.int64)
    chainids = atoms.getChainIdentifiers()
    if chainids is None: 
        chainids = np.zeros(n_atoms, '|S1')
    occupancies = atoms.getOccupancies()
    if occupancies is None:
        occupancies = np.zeros(n_atoms, np.float64)
    bfactors = atoms.getTempFactors()
    if bfactors is None:
        bfactors = np.zeros(n_atoms, np.float64)
    icodes = atoms.getInsertionCodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, '|S1')
    hetero = ['ATOM'] * n_atoms 
    heteroflags = atoms.getHeteroFlags()
    if heteroflags is not None:
        hetero = np.array(hetero, '|S6')
        hetero[heteroflags] = 'HETATM'
    elements = atoms.getElementSymbols()
    if elements is None:
        elements = np.zeros(n_atoms, '|S1')
    altlocs = atoms.getAltLocIndicators()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, '|S1')
    segments = atoms.getSegmentNames()
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
        coords = atoms.getCoordinates()
        for i, xyz in enumerate(coords):
            stream.write(line.format(hetero[i], i+1, atomnames[i], altlocs[i], 
                                     resnames[i], chainids[i], int(resnums[i]), icodes[i], 
                                     xyz[0], xyz[1], xyz[2], 
                                     occupancies[i], bfactors[i],  
                                     segments[i], elements[i].rjust(2)))
        if multi:
            stream.write('ENDMDL\n')
            altlocs = np.zeros(n_atoms, '|S1')
    atoms.setActiveCoordsetIndex(acsi)
    
def writePDB(filename, atoms, model=None, sort=None):
    """Write *atoms* in PDB format to a file with name *filename*.
    
    Returns *filename* is file is succesfully written. 
    See :func:`writePDBStream` for more information.
    """
    out = open(filename, 'w')
    writePDBStream(out, atoms, model, sort)
    out.close()
    return filename

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

def assignSecondaryStructure(header, atoms):
    """Assign secondary structure to alpha carbons in *atoms* from *header* 
    dictionary.

    *header* must be a dictionary parsed using the :func:`parsePDB`.
    *atoms* may be an instance of :class:`~prody.atomic.AtomGroup`, 
    :class:`~prody.atomic.Selection`, :class:`~prody.atomic.Chain` or :class:`~prody.atomic.Residue`. 

    The Dictionary of Protein Secondary Structure, in short DSSP, type 
    single letter codes assignments are used:     
    
      * **G** = 3-turn helix (310 helix). Min length 3 residues.
      * **H** = 4-turn helix (alpha helix). Min length 4 residues.
      * **I** = 5-turn helix (pi helix). Min length 5 residues.
      * **T** = hydrogen bonded turn (3, 4 or 5 turn)
      * **E** = extended strand in parallel and/or anti-parallel 
        beta-sheet conformation. Min length 2 residues.
      * **B** = residue in isolated beta-bridge (single pair beta-sheet hydrogen bond formation)
      * **S** = bend (the only non-hydrogen-bond based assignment).
    
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
    
    """
    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')
    helix = header.get('helix', {})
    sheet = header.get('sheet', {})
    if len(helix) == 0 and len(sheet) == 0:
        LOGGER.warning('header does not contain secondary structure data')
        return None
    ssa = atoms.getSecondaryStructureAssignments()
    if ssa is None:
        if isinstance(atoms, prody.AtomGroup):
            ag = atoms
        else:
            ag = atoms.getAtomGroup()
        ag.setSecondaryStructureAssignments(np.zeros(ag.getNumOfAtoms(), ATOMIC_DATA_FIELDS['secondary'].dtype))
            
    hierview = atoms.getHierView()
    count = 0
    for key, value in helix.iteritems():
        res = hierview.getResidue(*key)
        if res is None:
            continue
        atom = res.getAtom('CA')
        if atom is not None:
            count += 1
            atom.setSecondaryStructureAssignment(mapHelix[value[0]])
    for key, res in sheet.iteritems():
        res = hierview.getResidue(*key)
        if res is None:
            continue
        atom = res.getAtom('CA')
        if atom is not None:
            count += 1
            atom.setSecondaryStructureAssignment('E')
    LOGGER.info('Secondary structures were assigned to {0:d} residues.'.format(count))
    return atoms
            
            
            
def applyBiomolecularTransformations(header, atoms, biomol=None):
    """Returns *atoms* after applying biomolecular transformations from *header* 
    dictionary.
    
    Some PDB files contain transformations for more than 1 biomolecules. A 
    specific set of transformations can be choosen using *biomol* argument.
    Transformation sets are identified by integer numbers, e.g. 1, 2, ...
    
    """
    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')
    biomt = header.get('biomolecular_transformations', {})
    if len(biomt) == 0:
        LOGGER.warning('header does not contain biomolecular transformations')
        return None
    chids = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') 
    if isinstance(atoms, prody.AtomGroup):
        ag = atoms
    else:
        ag = atoms.getAtomGroup()
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

    keys.sort()
    for i in keys: 
        ags = []
        mt = biomt[i]
        if (len(mt) - 1) % 3 != 0:
            LOGGER.warning('Biomolecular {0:s} transformations were not applied'.format(i))
            continue
        for times in range((len(mt) - 1)/ 3):
            for chid in mt[0]:
                newag = ag.copy(atoms.select('chain ' + chid))
                newag.select('all').setChainIdentifiers(chids.pop(0))
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
                newag = t.apply(newag)
                ags.append(newag)
        if ags:
            newag = ags.pop(0)
            while ags:
                newag += ags.pop(0)
            newag.setName('{0:s} biomolecule {1:s}'.format(ag.getName(), i))
            biomols.append(newag)
    if biomols:
        if len(biomols) == 1:
            return biomols[0]
        else:
            return biomols
    else:
        return None
