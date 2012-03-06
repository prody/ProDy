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

"""This module defines functions for parsing and writing PDB files."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from collections import defaultdict
import os.path


import numpy as np

from prody.atomic import Atomic, Atom, AtomGroup
from prody.atomic import getBackboneAtomNames, getKeywordResnames 
from prody.atomic import ATOMIC_FIELDS
from prody.tools import openFile

from wwpdbftp import fetchPDB
from header import getHeaderDict, buildBiomolecules, assignSecstr

__all__ = ['parsePDBStream', 'parsePDB', 'parsePQR',
           'writePDBStream', 'writePDB',]

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

class PDBParseError(Exception):    
    pass
    
_parsePQRdoc = """
    :arg title: Title of the AtomGroup instance.  When ``None`` is passed,
        AtomGroup is named after the PDB filename.  
    :type title: str
    
    :arg ag: :class:`~.AtomGroup` instance for storing data parsed 
        from PDB file.  Number of atoms in *ag* and number of atoms parsed from
        the PDB file must be the same.  Atoms in *ag* and the PDB file must be 
        in the same order.  Non-coordinate data stored in *ag* will be 
        overwritten with those parsed from the file. 
    :type ag: :class:`~.AtomGroup`

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
        
    If ``model=0`` and ``header=True``, return header dictionary only.
    
    Note that this function does not evaluate ``CONECT`` records.
    
    """
    
_PDBSubsets = {'ca': 'ca', 'calpha': 'ca', 'bb': 'bb', 'backbone': 'bb'}

def parsePDB(pdb, **kwargs):
    """Return an :class:`~.AtomGroup` and/or dictionary containing 
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
    """Return an :class:`~.AtomGroup` and/or dictionary containing header data 
    parsed from a stream of PDB lines. 
    
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
            raise ValueError('{0:s} is not a valid subset'
                             .format(repr(subset)))
        title_suffix = '_' + _PDBSubsets[subset]
    if chain is not None:
        if not isinstance(chain, str):
            raise TypeError('chain must be a string')
        elif len(chain) == 0:
            raise ValueError('chain must not be an empty string')
        title_suffix = '_' + chain + title_suffix
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    else:
        ag = AtomGroup(str(kwargs.get('title', kwargs.get('name', 
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
            hd, split = getHeaderDict(lines)
        _parsePDBLines(ag, lines, split, model, chain, subset, altloc)
        if ag.numAtoms() > 0:
            LOGGER.timing('{0:d} atoms and {1:d} coordinate set(s) were '
                          'parsed in %.2fs.'.format(ag.numAtoms(), 
                           ag.numCoordsets() - n_csets))
        else:
            ag = None
            LOGGER.warn('Atomic data could not be parsed, please '
                        'check the input file.')
    elif header:
        hd, split = getHeaderDict(stream)

    if secondary:
        ag = assignSecstr(hd, ag)
    if biomol:
        ag = buildBiomolecules(hd, ag)

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
    """Return an :class:`~.AtomGroup` containing data parsed 
    from PDB lines. 
    
    :arg filename: a PQR filename
    :type filename: str"""
    
    title = kwargs.get('title', kwargs.get('name'))
    model = 1
    header = False
    chain = kwargs.get('chain')
    subset = kwargs.get('subset')
    altloc = kwargs.get('altloc', 'A')
    if not os.path.isfile(filename):
        raise IOError('No such file: {0:s}'.format(repr(filename)))
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
            raise ValueError('{0:s} is not a valid subset'
                             .format(repr(subset)))
        title_suffix = '_' + _PDBSubsets[subset]
    if chain is not None:
        if not isinstance(chain, str):
            raise TypeError('chain must be a string')
        elif len(chain) == 0:
            raise ValueError('chain must not be an empty string')
        title_suffix = '_' + chain + title_suffix
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    else:
        ag = AtomGroup(title + title_suffix)
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
                   altloc_torf, format='PDB'):
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
            subset = set(getBackboneAtomNames())
        only_subset = True
        protein_resnames = set(getKeywordResnames('protein'))
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
    addcoords = False
    if atomgroup.numCoordsets() > 0:
        addcoords = True
    alength = asize
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    hetero = np.zeros(asize, dtype=ATOMIC_FIELDS['hetero'].dtype)
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    if isPDB:
        segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
        elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)
        bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
        occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)
        secondary = None
        anisou = None
        siguij = None
    else:
        charges = np.zeros(asize, dtype=ATOMIC_FIELDS['charge'].dtype)
        radii = np.zeros(asize, dtype=ATOMIC_FIELDS['radius'].dtype)
        
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
                    LOGGER.warn('Discarding model {0:d}, which contains more '
                            'atoms than first model does.'.format(nmodel+1))
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
                    LOGGER.warn('failed to parse occupancy at line {0:d}'
                                .format(i))
                try:
                    bfactors[acount] = line[60:66]
                except:
                    LOGGER.warn('failed to parse beta-factor at line {0:d}'
                                .format(i))
                hetero[acount] = startswith[0] == 'H'
                segnames[acount] = line[72:76]
                elements[acount] = line[76:78]
            else:
                try:
                    charges[acount] = line[54:62]
                except:
                    LOGGER.warn('failed to parse charge at line {0:d}'
                                .format(i))
                try:
                    radii[acount] = line[62:69]
                except:
                    LOGGER.warn('failed to parse radius at line {0:d}'
                                .format(i))
            acount += 1
            if n_atoms == 0 and acount >= alength:
                # if arrays are short extend them with zeros
                alength += asize
                coordinates = np.concatenate(
                    (coordinates, np.zeros((asize, 3), float)))
                atomnames = np.concatenate((atomnames,
                    np.zeros(asize, ATOMIC_FIELDS['name'].dtype)))
                resnames = np.concatenate((resnames,
                    np.zeros(asize, ATOMIC_FIELDS['resname'].dtype)))
                resnums = np.concatenate((resnums,
                    np.zeros(asize, ATOMIC_FIELDS['resnum'].dtype)))
                chainids = np.concatenate((chainids,
                    np.zeros(asize, ATOMIC_FIELDS['chain'].dtype)))
                hetero = np.concatenate((hetero,
                    np.zeros(asize, ATOMIC_FIELDS['hetero'].dtype)))
                altlocs = np.concatenate((altlocs,
                    np.zeros(asize, ATOMIC_FIELDS['altloc'].dtype)))
                icodes = np.concatenate((icodes,
                    np.zeros(asize, ATOMIC_FIELDS['icode'].dtype)))
                serials = np.concatenate((serials,
                    np.zeros(asize, ATOMIC_FIELDS['serial'].dtype)))
                if isPDB:
                    bfactors = np.concatenate((bfactors,
                        np.zeros(asize, ATOMIC_FIELDS['beta'].dtype)))
                    occupancies = np.concatenate((occupancies,
                        np.zeros(asize, ATOMIC_FIELDS['occupancy'].dtype)))
                    segnames = np.concatenate((segnames,
                        np.zeros(asize, ATOMIC_FIELDS['segment'].dtype)))
                    elements = np.concatenate((elements,
                        np.zeros(asize, ATOMIC_FIELDS['element'].dtype)))
                    if anisou is not None:
                        anisou = np.concatenate((anisou, np.zeros((asize, 6), 
                            ATOMIC_FIELDS['anisou'].dtype)))
                    if siguij is not None:
                        siguij = np.concatenate((siguij, np.zeros((asize, 6), 
                            ATOMIC_FIELDS['siguij'].dtype)))
                else:
                    charges = np.concatenate((charges,
                        np.zeros(asize, ATOMIC_FIELDS['charge'].dtype)))
                    radii = np.concatenate((radii,
                        np.zeros(asize, ATOMIC_FIELDS['radius'].dtype)))
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
                    LOGGER.warn('Discarding model {0:d}, which contains '
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
                    if addcoords:
                        atomgroup.addCoordset(coordinates[:acount])
                    else:
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
                    dtype=ATOMIC_FIELDS['anisou'].dtype)
            try:
                index = acount - 1
                anisou[index, 0] = line[28:35]
                anisou[index, 1] = line[35:42]
                anisou[index, 2] = line[43:49]
                anisou[index, 3] = line[49:56]
                anisou[index, 4] = line[56:63]
                anisou[index, 5] = line[63:70]
            except:
                LOGGER.warn('failed to parse anisotropic temperature '
                    'factors at line {0:d}'.format(i))
        elif isPDB and startswith =='SIGUIJ':
            if siguij is None:
                siguij = np.zeros((alength, 6), 
                    dtype=ATOMIC_FIELDS['siguij'].dtype)
            try:
                index = acount - 1
                siguij[index, 0] = line[28:35]
                siguij[index, 1] = line[35:42]
                siguij[index, 2] = line[43:49]
                siguij[index, 3] = line[49:56]
                siguij[index, 4] = line[56:63]
                siguij[index, 5] = line[63:70]
            except:
                LOGGER.warn('failed to parse standard deviations of '
                    'anisotropic temperature factors at line {0:d}'.format(i))
        elif startswith =='SIGATM':
            pass
        i += 1
    if onlycoords:
        if acount == atomgroup.numAtoms():
            coordsets[nmodel] = coordinates
            nmodel += 1
        if nmodel == coordsets.shape[0]:
            if addcoords:
                atomgroup.addCoordset(coordsets)
            else:
                atomgroup.setCoords(coordsets)
        else:
            if addcoords:
                atomgroup.addCoordset(coordsets[:nmodel])
            else:
                atomgroup.setCoords(coordsets[:nmodel])
    elif not END:
        # this means last line wast an ATOM line, so atomgroup is not decorated
        if addcoords:
            atomgroup.addCoordset(coordinates[:acount])
        else:
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
                    LOGGER.warn("failed to parse altloc {0:s} at line {1:d}, "
                                "residue does not exist as altloc 'A'".format(
                                repr(key), i+1))
                    continue
                rn = resnames[ids[0]]
                ans = atomnames[ids]
                indices[(ach, ari)] = (rn, ids, ans)
            if rn != arn:
                LOGGER.warn("failed to parse altloc {0:s} at line {1:d}, "
                            "residue names do not match (expected {2:s}, "
                            "parsed {3:s})".format(repr(key), i+1, repr(rn), 
                                                   repr(arn)))
                continue
            index = ids[(ans == aan).nonzero()[0]]
            if len(index) != 1:
                LOGGER.warn("failed to parse altloc {0:s} at line {1:d}, could"
                            " not identify matching atom ({2:s} not found in "
                            "the residue)".format(repr(key), i+1, repr(aan)))
                continue
            try:
                xyz[index[0], 0] = float(line[30:38])
                xyz[index[0], 1] = float(line[38:46])
                xyz[index[0], 2] = float(line[46:54])
            except:
                LOGGER.warn('failed to parse altloc {0:s} at line {1:d}, could'
                            ' not read coordinates'.format(repr(key), i+1))
                continue
            success += 1
        LOGGER.info('{0:d} out of {1:d} altloc {2:s} lines were parsed.'
                    .format(success, len(lines), repr(key)))
        if success > 0:
            LOGGER.info('Altloc {0:s} is appended as a coordinate set to the '
                        'atom group.'.format(repr(key), atomgroup.getTitle()))
            atomgroup.addCoordset(xyz)

_writePDBdoc = """
    :arg atoms: Atomic data container.
    :type atoms: :class:`~.Atomic` 
    
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
    
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms does not have a valid type')
    if isinstance(atoms, Atom):
        atoms = atoms.select('all')

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
    
    acsi = atoms.getACSIndex()
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
        atoms.setACSIndex(m)
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
    atoms.setACSIndex(acsi)

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
            

