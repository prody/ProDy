# -*- coding: utf-8 -*-
"""This module defines functions for parsing and writing `PDB files`_.

.. _PDB files: http://www.wwpdb.org/documentation/format32/v3.2.html"""

from collections import defaultdict
import os.path
import time
from numbers import Integral

import numpy as np

from prody.atomic import AtomGroup, Atom, Selection
from prody.atomic import flags
from prody.atomic import ATOMIC_FIELDS
from prody.utilities import openFile, isListLike
from prody.utilities.misctools import decToHybrid36, hybrid36ToDec, packmolRenumChains
from prody import LOGGER, SETTINGS

from .header import getHeaderDict, buildBiomolecules, assignSecstr, isHelix, isSheet
from .localpdb import fetchPDB
from .ciffile import parseMMCIF
from .emdfile import parseEMD

__all__ = ['parsePDBStream', 'parsePDB', 'parseChainsList', 'parsePQR',
           'writePDBStream', 'writePDB', 'writeChainsList', 'writePQR',
           'writePQRStream']

MAX_N_ATOM = 99999 
MAX_N_RES = 9999

class PDBParseError(Exception):
    pass

_parsePQRdoc = """
    :arg title: title of the :class:`.AtomGroup` instance, default is the
        PDB filename or PDB identifier
    :type title: str

    :arg ag: :class:`.AtomGroup` instance for storing data parsed from PDB
        file, number of atoms in *ag* and number of atoms parsed from the PDB
        file must be the same and atoms in *ag* and those in PDB file must
        be in the same order.  Non-coordinate data stored in *ag* will be
        overwritten with those parsed from the file.
    :type ag: :class:`.AtomGroup`

    :arg chain: chain identifiers for parsing specific chains, e.g.
        ``chain='A'``, ``chain='B'``, ``chain='DE'``, by default all
        chains are parsed
    :type chain: str

    :arg subset: a predefined keyword to parse subset of atoms, valid keywords
        are ``'calpha'`` (``'ca'``), ``'backbone'`` (``'bb'``), or **None**
        (read all atoms), e.g. ``subset='bb'``
    :type subset: str
"""

_parsePDBdoc = _parsePQRdoc + """
    :arg model: model index or None (read all models), e.g. ``model=10``
    :type model: int, list

    :arg header: if **True** PDB header content will be parsed and returned
    :type header: bool

    :arg altloc: if a location indicator is passed, such as ``'A'`` or ``'B'``,
         only indicated alternate locations will be parsed as the single
         coordinate set of the AtomGroup,  if *altloc* is set **True** all
         alternate locations will be parsed and each will be appended as a
         distinct coordinate set, default is ``"A"``
    :type altloc: str

    :arg biomol: if **True**, biomolecules are obtained by transforming the
        coordinates using information from header section will be returned.
        This option uses :func:`.buildBiomolecules` and as noted there, 
        atoms in biomolecules are ordered according to the original chain 
        IDs. Chains may have the same chain ID, in which case they are given 
        different segment names.
        Default is **False**
    :type biomol: bool

    :arg secondary: if **True**, secondary structure information from header
        section will be assigned to atoms.
        Default is **False**
    :type secondary: bool

    If ``model=0`` and ``header=True``, return header dictionary only.

    """

_PDBSubsets = {'ca': 'ca', 'calpha': 'ca', 'bb': 'bb', 'backbone': 'bb'}

def parsePDB(*pdb, **kwargs):
    """Returns an :class:`.AtomGroup` and/or dictionary containing header data
    parsed from a PDB file.

    This function extends :func:`.parsePDBStream`.

    See :ref:`parsepdb` for a detailed usage example.

    :arg pdb: one PDB identifier or filename, or a list of them.
        If needed, PDB files are downloaded using :func:`.fetchPDB()` function.
    
    You can also provide arguments that you would like passed on to fetchPDB().

    :arg extend_biomol: whether to extend the list of results with a list
        rather than appending, which can create a mixed list, 
        especially when biomol=True.
        Default value is False to reproduce previous behaviour.
        This value is ignored when result is not a list (header=True or model=0).
    :type extend_biomol: bool 
    """
    extend_biomol = kwargs.pop('extend_biomol', False)

    n_pdb = len(pdb)
    if n_pdb == 0:
        raise ValueError('Please provide a PDB ID or filename')

    if n_pdb == 1:
        if isListLike(pdb[0]) or isinstance(pdb[0], dict):
            pdb = pdb[0]
            n_pdb = len(pdb)

    if n_pdb == 1:
        return _parsePDB(pdb[0], **kwargs)
    else:
        results = []
        lstkwargs = {}
        for key in kwargs:
            argval = kwargs.get(key)
            if np.isscalar(argval):
                argval = [argval]*n_pdb
            lstkwargs[key] = argval

        start = time.time()
        LOGGER.progress('Retrieving {0} PDB structures...'
                    .format(n_pdb), n_pdb, '_prody_parsePDB')
        for i, p in enumerate(pdb):
            kwargs = {}
            for key in lstkwargs:
                kwargs[key] = lstkwargs[key][i]
            c = kwargs.get('chain','')
            LOGGER.update(i, 'Retrieving {0}...'.format(p+c), 
                          label='_prody_parsePDB')
            result = _parsePDB(p, **kwargs)
            if not isinstance(result, tuple):
                if isinstance(result, dict):
                    result = (None, result)
                else:
                    result = (result, None)
            results.append(result)

        results = list(zip(*results))
        LOGGER.finish()

        for i in reversed(range(len(results))):
            if all(j is None for j in results[i]):
                results.pop(i)
        if len(results) == 1:
            results = results[0]
        results = list(results)

        model = kwargs.get('model')
        header = kwargs.get('header', False)
        if model != 0 and header:
            numPdbs = len(results[0])
        else:
            numPdbs = len(results)

            if extend_biomol:
                results_old = results
                results = []
                for entry in results_old:
                    if isinstance(entry, AtomGroup):
                        results.append(entry)
                    else:
                        results.extend(entry)

        LOGGER.info('{0} PDBs were parsed in {1:.2f}s.'
                     .format(numPdbs, time.time()-start))

        return results

def _getPDBid(pdb):
    l = len(pdb)
    if l == 4:
        pdbid, chain = pdb, ''
    elif l == 5:
        pdbid = pdb[:4]; chain = pdb[-1]
    elif ':' in pdb:
        i = pdb.find(':')
        pdbid = pdb[:i]; chain = pdb[i+1:]
    else:
        raise IOError('{0} is not a valid filename or a valid PDB '
                      'identifier.'.format(pdb))
    if not pdbid.isalnum():
        raise IOError('{0} is not a valid filename or a valid PDB '
                      'identifier.'.format(pdb))
    if chain != '' and not chain.isalnum():
        raise IOError('{0} is not a valid chain identifier.'.format(chain))
    return pdbid, chain

def _parsePDB(pdb, **kwargs):
    title = kwargs.get('title', None)
    chain = ''
    if not os.path.isfile(pdb):
        pdb, chain = _getPDBid(pdb)
        if title is None:
            title = pdb
            kwargs['title'] = title
        filename = fetchPDB(pdb, **kwargs)
        if filename is None:
            try:
                LOGGER.warn("Trying to parse mmCIF file instead")
                return parseMMCIF(pdb+chain, **kwargs)
            except:
                try:
                    LOGGER.warn("Trying to parse EMD file instead")
                    return parseEMD(pdb+chain, **kwargs)
                except:
                    raise IOError('PDB file for {0} could not be downloaded.'
                                .format(pdb))
        pdb = filename
    if title is None:
        title, ext = os.path.splitext(os.path.split(pdb)[1])
        if ext == '.gz':
            title, ext = os.path.splitext(title)
        if len(title) == 7 and title.startswith('pdb'):
            title = title[3:]
        kwargs['title'] = title

    if pdb.endswith('.pdb') or pdb.endswith('.pdb.gz'):
        stream = openFile(pdb, 'rt')
        if chain != '':
            kwargs['chain'] = chain
        result = parsePDBStream(stream, **kwargs)
        stream.close()
        return result
    else:
        try:
            LOGGER.warn("Trying to parse as mmCIF file instead")
            return parseMMCIF(pdb, **kwargs)
        except KeyError:
            try:
                LOGGER.warn("Trying to parse as EMD file instead")
                return parseEMD(pdb, **kwargs)
            except ValueError:
                LOGGER.warn("Could not parse anything so returning None")
                return None
            except:
                raise IOError('PDB file for {0} could not be downloaded.'
                              .format(pdb))

parsePDB.__doc__ += _parsePDBdoc

def parsePDBStream(stream, **kwargs):
    """Returns an :class:`.AtomGroup` and/or dictionary containing header data
    parsed from a stream of PDB lines.

    :arg stream: Anything that implements the method ``readlines``
        (e.g. :class:`file`, buffer, stdin)""" 

    model = kwargs.get('model')
    header = kwargs.get('header', False)
    assert isinstance(header, bool), 'header must be a boolean'
    chain = kwargs.get('chain')
    subset = kwargs.get('subset')
    altloc = kwargs.get('altloc', 'A')
    packmol = kwargs.get('packmol', False)

    auto_bonds = SETTINGS.get('auto_bonds')
    get_bonds = kwargs.get('bonds', auto_bonds)

    if model is not None:
        if isinstance(model, Integral):
            if model < 0:
                raise ValueError('model must be greater than 0')
        else:
            raise TypeError('model must be an integer, {0} is invalid'
                            .format(str(model)))
    title_suffix = ''
    if subset:
        try:
            subset = _PDBSubsets[subset.lower()]
        except AttributeError:
            raise TypeError('subset must be a string')
        except KeyError:
            raise ValueError('{0} is not a valid subset'
                             .format(repr(subset)))
        title_suffix = '_' + subset
    if chain is not None:
        if not isinstance(chain, str):
            raise TypeError('chain must be a string')
        elif len(chain) == 0:
            raise ValueError('chain must not be an empty string')
        title_suffix = chain + title_suffix
    ag = kwargs.pop('ag', None)
    if ag is not None:
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    elif model != 0:
        ag = AtomGroup(str(kwargs.get('title', 'Unknown')) + title_suffix)
        n_csets = 0

    biomol = kwargs.get('biomol', False)
    auto_secondary = SETTINGS.get('auto_secondary')
    secondary = kwargs.get('secondary', auto_secondary)
    split = 0
    hd = None
    if model != 0:
        LOGGER.timeit()
        try:
            lines = stream.readlines()
        except AttributeError as err:
            try:
                lines = stream.read().split('\n')
            except AttributeError:
                raise err
        if not len(lines):
            raise ValueError('empty PDB file or stream')
        if header or biomol or secondary:
            hd, split = getHeaderDict(lines)
        bonds = [] if get_bonds else None
        _parsePDBLines(ag, lines, split, model, chain, subset, altloc, bonds=bonds)
        if bonds:
            try:
                ag.setBonds(bonds)
            except ValueError:
                LOGGER.warn('Bonds read from CONECT records do not apply to subset so were not added')
        if ag.numAtoms() > 0:
            LOGGER.report('{0} atoms and {1} coordinate set(s) were '
                          'parsed in %.2fs.'.format(ag.numAtoms(),
                          ag.numCoordsets() - n_csets))
        else:
            ag = None
            LOGGER.warn('Atomic data could not be parsed, please '
                        'check the input file.')
    elif header:
        hd, split = getHeaderDict(stream)

    if ag is not None and isinstance(hd, dict):
        if secondary:
            if auto_secondary:
                try:
                    ag = assignSecstr(hd, ag)
                except ValueError:
                    pass
            else:
                ag = assignSecstr(hd, ag)
        if biomol:
            ag = buildBiomolecules(hd, ag)

            if isinstance(ag, list):
                LOGGER.info('Biomolecular transformations were applied, {0} '
                            'biomolecule(s) are returned.'.format(len(ag)))
            else:
                LOGGER.info('Biomolecular transformations were applied to the '
                            'coordinate data.')

    if packmol:
        ag = packmolRenumChains(ag)

    if model != 0:
        if header:
            return ag, hd
        else:
            return ag
    else:
        return hd

parsePDBStream.__doc__ += _parsePDBdoc


def parsePQR(filename, **kwargs):
    """Returns an :class:`.AtomGroup` containing data parsed from PDB lines.

    :arg filename: a PQR filename
    :type filename: str"""

    title = kwargs.get('title', kwargs.get('name'))
    chain = kwargs.get('chain')
    subset = kwargs.get('subset')
    if not os.path.isfile(filename):
        raise IOError('No such file: {0}'.format(repr(filename)))
    if title is None:
        fn, ext = os.path.splitext(os.path.split(filename)[1])
        if ext == '.gz':
            fn, ext = os.path.splitext(fn)
        title = fn.lower()
    title_suffix = ''
    if subset:
        try:
            subset = _PDBSubsets[subset.lower()]
        except AttributeError:
            raise TypeError('subset must be a string')
        except KeyError:
            raise ValueError('{0} is not a valid subset'
                             .format(repr(subset)))
        title_suffix = '_' + subset
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

    pqr = openFile(filename, 'rt')
    lines = pqr.readlines()
    pqr.close()
    LOGGER.timeit()
    ag = _parsePDBLines(ag, lines, split=0, model=1, chain=chain,
                        subset=subset, altloc_torf=False, format='pqr')
    if ag.numAtoms() > 0:
        LOGGER.report('{0} atoms and {1} coordinate sets were '
                      'parsed in %.2fs.'.format(ag.numAtoms(),
                      ag.numCoordsets() - n_csets))
        return ag
    else:
        return None

parsePQR.__doc__ += _parsePQRdoc

def _parsePDBLines(atomgroup, lines, split, model, chain, subset,
                   altloc_torf, format='PDB', bonds=None):
    """Returns an AtomGroup. See also :func:`.parsePDBStream()`.

    :arg lines: PDB/PQR lines
    :arg split: starting index for coordinate data lines"""

    format = format.upper()
    isPDB = format == 'PDB'

    num_ters = 0

    if subset:
        if subset == 'ca':
            subset = set(('CA',))
        elif subset in 'bb':
            subset = flags.BACKBONE
        only_subset = True
        protein_resnames = flags.AMINOACIDS
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
        asize = len(lines) - split
    addcoords = False
    if atomgroup.numCoordsets() > 0:
        addcoords = True
    alength = asize
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    hetero = np.zeros(asize, dtype=bool)
    termini = np.zeros(asize, dtype=bool)
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    charges = np.zeros(asize, dtype=ATOMIC_FIELDS['charge'].dtype)
    if isPDB:
        segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
        elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)
        bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
        occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)
        anisou = None
        siguij = None
    else:
        radii = np.zeros(asize, dtype=ATOMIC_FIELDS['radius'].dtype)

    asize = 2000 # increase array length by this much when needed

    start = split
    stop = len(lines)
    nmodel = 0
    bonds_serial = []
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
            raise PDBParseError('model {0} is not found'.format(model))

    if isinstance(altloc_torf, str):
        if altloc_torf == 'all':
            which_altlocs = 'all'
        elif altloc_torf.strip() != 'A':
            LOGGER.info('Parsing alternate locations {0}.'
                        .format(altloc_torf))
            which_altlocs = ' ' + ''.join(altloc_torf.split())
        else:
            which_altlocs = ' A'
        altloc_torf = False
    else:
        which_altlocs = ' A'
        altloc_torf = True

    acount = 0
    coordsets = None
    altloc = defaultdict(list)
    i = start
    END = False
    warned_5_digit = False
    dec = True
    while i < stop:
        line = lines[i]
        if not isPDB:
            fields = line.split()
            if len(fields) == 10:
                fields.insert(4, '')
            elif len(fields) != 11:
                LOGGER.warn('wrong number of fields for PQR format at line %d'%i)
                i += 1
                continue

        if isPDB:
            startswith = line[0:6].strip()
        else:
            startswith = fields[0]

        if startswith == 'ATOM' or startswith == 'HETATM':
            if isPDB:
                atomname = line[12:16].strip()
                resname = line[17:20].strip()
            else:
                atomname= fields[2]
                resname = fields[3]

            if only_subset:
                if not (atomname in subset and resname in protein_resnames):
                    i += 1
                    continue

            if isPDB:
                chid = line[20:22].strip()
            else:
                chid = fields[4]

            if only_chains:
                if not chid in chain:
                    i += 1
                    continue

            if isPDB:
                alt = line[16]
                if alt not in which_altlocs and which_altlocs != 'all':
                    altloc[alt].append((line, i))
                    i += 1
                    continue
            else:
                alt = ' '
            try:
                if isPDB:
                    coordinates[acount, 0] = line[30:38]
                    coordinates[acount, 1] = line[38:46]
                    coordinates[acount, 2] = line[46:54]
                else:
                    coordinates[acount, 0] = fields[6]
                    coordinates[acount, 1] = fields[7]
                    coordinates[acount, 2] = fields[8]
            except:
                if acount >= n_atoms > 0:
                    if nmodel == 0:
                        raise ValueError(format + 'file and AtomGroup ag must '
                                         'have same number of atoms')
                    LOGGER.warn('Discarding model {0}, which contains {1} more '
                                'atoms than first model does.'
                                .format(nmodel+1,acount-n_atoms+1))
                    acount = 0
                    nmodel += 1
                    coordinates = np.zeros((n_atoms, 3), dtype=float)
                    if isPDB:
                        while lines[i][:6] != 'ENDMDL':
                            i += 1
                else:
                    raise PDBParseError('invalid or missing coordinate(s) at '
                                         'line {0}'.format(i+1))
            if onlycoords:
                acount += 1
                i += 1
                continue

            serial_str = line[6:11] if isPDB else fields[1]
            try:
                serials[acount] = int(serial_str)
            except ValueError:
                try:
                    isnumeric = np.alltrue([x.isdigit() for x in serial_str])
                    if not isnumeric and serial_str == serial_str.upper():
                        serials[acount] = hybrid36ToDec(serial_str)
                    else:
                        # lower case is found in hexadecimal PDB files
                        serials[acount] = int(serial_str, 16)
                except ValueError:
                    if acount > 0:
                        LOGGER.warn('failed to parse serial number in line {0}. Assigning it by incrementing.'
                                    .format(i))
                        serials[acount] = serials[acount-1]+1
                    else:
                        LOGGER.warn('failed to parse serial number in line {0}. Assigning it as 1.'
                                    .format(i))
                        serials[acount] = 1

            altlocs[acount] = alt
            atomnames[acount] = atomname
            resnames[acount] = resname
            chainids[acount] = chid
            if isPDB:
                resnum_str = line[22:26]
                icode = line[26]
                if dec:
                    try:
                        resnum = int(resnum_str)
                    except ValueError:
                        dec = False

                    if icode.isdigit():
                        if not warned_5_digit:
                            LOGGER.warn('parsed 5 digit residue number including numeric insertion code')
                            warned_5_digit = True
                        resnum = int(str(resnum) + icode)
                    else:
                        icodes[acount] = icode

                if dec and acount > 2 and resnums[acount-2] > resnum and resnums[acount-2] >= MAX_N_RES:
                    dec = False

                if not dec:
                    resnum = resnum_str
                    try:
                        isnumeric = np.alltrue([x.isdigit() or x==' ' for x in resnum_str])
                        if not isnumeric and resnum_str == resnum_str.upper():
                            resnum = hybrid36ToDec(resnum_str, resnum=True)
                        else:
                            # lower case is found in hexadecimal PDB files
                            resnum = int(resnum_str, 16)

                    except ValueError:
                        if acount > 0:
                            LOGGER.warn('failed to parse residue number in line {0}. Assigning it by incrementing.'
                                        .format(i))
                            resnum = resnums[acount-1]+1
                        else:
                            LOGGER.warn('failed to parse residue number in line {0}. Assigning it as 1.'
                                        .format(i))
                            resnum = 1

                resnums[acount] = resnum
            else:
                resnum = fields[5]
                if resnum[-1].isalpha():
                    icode = resnum[-1]
                else:
                    icode = ' '
                try:
                    resnums[acount] = int(resnum)
                except ValueError:
                    try:
                        resnums[acount] = int(resnum, 16)
                    except ValueError:
                        LOGGER.warn('failed to parse residue number in line {0}. Assigning it by incrementing.'
                                    .format(i))
                        resnums[acount] = resnums[acount-1]+1
                icodes[acount] = icode

            if isPDB:
                try:
                    occupancies[acount] = line[54:60]
                except:
                    LOGGER.warn('failed to parse occupancy at line {0}'
                                .format(i))
                try:
                    bfactors[acount] = line[60:66]
                except:
                    LOGGER.warn('failed to parse beta-factor at line {0}'
                                .format(i))
                hetero[acount] = startswith[0] == 'H'
                segnames[acount] = line[72:76]
                elements[acount] = line[76:78]
                try:
                    charges[acount] = int(line[79] + line[78])
                except:
                    charges[acount] = 0
            else:
                try:
                    charges[acount] = fields[9]
                except:
                    LOGGER.warn('failed to parse charge at line {0}'
                                .format(i))
                try:
                    radii[acount] = fields[10]
                except:
                    LOGGER.warn('failed to parse radius at line {0}'
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
                hetero = np.concatenate((hetero, np.zeros(asize, bool)))
                termini = np.concatenate((termini, np.zeros(asize, bool)))
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
        elif startswith == 'CONECT':
            if bonds is not None:
                atom_serial = line[6:11]
                bonded1_serial = line[11:16]
                bonds_serial.append([int(atom_serial), int(bonded1_serial)])

                bonded2_serial = line[16:21]
                if len(bonded2_serial.strip()):
                    bonds_serial.append([int(atom_serial), int(bonded2_serial)])

                bonded3_serial = line[21:26]
                if len(bonded3_serial.strip()):
                    bonds_serial.append([int(atom_serial), int(bonded3_serial)])

                bonded4_serial = line[26:31]  # fixed typo
                if len(bonded4_serial.strip()):
                    bonds_serial.append([int(atom_serial), int(bonded4_serial)])

        elif not onlycoords and (startswith == 'TER   ' or
            startswith.strip() == 'TER'):
            termini[acount - 1] = True
            num_ters += 1
        elif startswith == 'ENDMDL' or startswith[:3] == 'END':
            if acount == 0:
                # If there is no atom record between ENDMDL & END skip to next
                i += 1
                continue
            if model is not None:
                if bonds is None:
                    i += 1
                    break
                else:
                    i += 1
                    for j in range(i, stop):
                        if lines[j].startswith('CONECT'):
                            i = j
                            break
                    continue
            diff = stop - i - 1
            END = diff < acount
            if coordsets is not None:
                END = END or nmodel >= coordsets.shape[0]
            if onlycoords:
                if acount < n_atoms:
                    LOGGER.warn('Discarding model {0}, which contains '
                                '{1} fewer atoms than the first model '
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
                    coordinates.resize((acount, 3), refcheck=False)
                    if addcoords:
                        atomgroup.addCoordset(coordinates)
                    else:
                        atomgroup._setCoords(coordinates)
                else:
                    coordsets = np.zeros((int(diff//acount+1), acount, 3))
                    coordsets[0] = coordinates[:acount]
                    onlycoords = True
                atomnames.resize(acount, refcheck=False)
                resnames.resize(acount, refcheck=False)
                resnums.resize(acount, refcheck=False)
                chainids.resize(acount, refcheck=False)
                hetero.resize(acount, refcheck=False)
                termini.resize(acount, refcheck=False)
                altlocs.resize(acount, refcheck=False)
                icodes.resize(acount, refcheck=False)
                serials.resize(acount, refcheck=False)
                if not only_subset:
                    atomnames = np.char.strip(atomnames)
                    resnames = np.char.strip(resnames)
                atomgroup.setNames(atomnames)
                atomgroup.setResnames(resnames)
                atomgroup.setResnums(resnums)
                atomgroup.setChids(chainids)
                atomgroup.setFlags('hetatm', hetero)
                atomgroup.setFlags('pdbter', termini)
                atomgroup.setAltlocs(altlocs)
                atomgroup.setIcodes(np.char.strip(icodes))
                atomgroup.setSerials(serials)
                if isPDB:
                    bfactors.resize(acount, refcheck=False)
                    occupancies.resize(acount, refcheck=False)
                    segnames.resize(acount, refcheck=False)
                    elements.resize(acount, refcheck=False)
                    atomgroup.setBetas(bfactors)
                    atomgroup.setOccupancies(occupancies)
                    atomgroup.setSegnames(np.char.strip(segnames))
                    atomgroup.setElements(np.char.strip(elements))
                    from prody.utilities.misctools import getMasses
                    atomgroup.setMasses(getMasses(np.char.strip(elements)))
                    if anisou is not None:
                        anisou.resize((acount, 6), refcheck=False)
                        atomgroup.setAnisous(anisou / 10000)
                    if siguij is not None:
                        siguij.resize((acount, 6), refcheck=False)
                        atomgroup.setAnistds(siguij / 10000)
                else:
                    radii.resize(acount, refcheck=False)
                    atomgroup.setRadii(radii)
                charges.resize(acount, refcheck=False)
                atomgroup.setCharges(charges)

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

            alt = line[16]
            if alt not in which_altlocs and which_altlocs != 'all':
                altloc[alt].append((line, i))
                i += 1
                continue
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
                    'factors at line {0}'.format(i))
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
                    'anisotropic temperature factors at line {0}'.format(i))
        elif startswith =='SIGATM':
            pass
        i += 1
    if onlycoords:
        if acount == atomgroup.numAtoms():
            coordsets[nmodel] = coordinates
            nmodel += 1
        del coordinates
        coordsets.resize((nmodel, atomgroup.numAtoms(), 3), refcheck=False)
        if addcoords:
            atomgroup.addCoordset(coordsets)
        else:
            atomgroup._setCoords(coordsets)
    elif not END:
        # this means last line was an ATOM line, so atomgroup is not decorated
        coordinates.resize((acount, 3), refcheck=False)
        if addcoords:
            atomgroup.addCoordset(coordinates)
        else:
            atomgroup._setCoords(coordinates)
        atomnames.resize(acount, refcheck=False)
        resnames.resize(acount, refcheck=False)
        resnums.resize(acount, refcheck=False)
        chainids.resize(acount, refcheck=False)
        hetero.resize(acount, refcheck=False)
        termini.resize(acount, refcheck=False)
        altlocs.resize(acount, refcheck=False)
        icodes.resize(acount, refcheck=False)
        serials.resize(acount, refcheck=False)
        if not only_subset:
            atomnames = np.char.strip(atomnames)
            resnames = np.char.strip(resnames)
        atomgroup.setNames(atomnames)
        atomgroup.setResnames(resnames)
        atomgroup.setResnums(resnums)
        atomgroup.setChids(chainids)
        atomgroup.setFlags('hetatm', hetero)
        atomgroup.setFlags('pdbter', termini)
        atomgroup.setAltlocs(altlocs)
        atomgroup.setIcodes(np.char.strip(icodes))
        atomgroup.setSerials(serials)
        if isPDB:
            if anisou is not None:
                anisou.resize((acount, 6), refcheck=False)
                atomgroup.setAnisous(anisou / 10000)
            if siguij is not None:
                siguij.resize((acount, 6), refcheck=False)
                atomgroup.setAnistds(siguij / 10000)
            bfactors.resize(acount, refcheck=False)
            occupancies.resize(acount, refcheck=False)
            segnames.resize(acount, refcheck=False)
            elements.resize(acount, refcheck=False)
            atomgroup.setSegnames(np.char.strip(segnames))
            atomgroup.setElements(np.char.strip(elements))
            from prody.utilities.misctools import getMasses
            atomgroup.setMasses(getMasses(np.char.strip(elements)))
            atomgroup.setBetas(bfactors)
            atomgroup.setOccupancies(occupancies)
        else:
            radii.resize(acount, refcheck=False)
            atomgroup.setRadii(radii)
        charges.resize(acount, refcheck=False)
        atomgroup.setCharges(charges)

    # rematch the bond-atom serial numbers to atom ids
    serial_to_id = {int(serial): aidx for aidx, serial in enumerate(serials)}
    for bond in bonds_serial:
        try:
            bonds.append([serial_to_id[bond[0]], serial_to_id[bond[1]]])
        except KeyError:
            LOGGER.warn("Bond connecting atom serial numbers {0}"
                        " and {1} contains atoms not included in the model"
                        .format(bond[0], bond[1])
                        )

    if altloc and altloc_torf:
        _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames)

    return atomgroup

def _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames):
    altloc_keys = list(altloc)
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
                    LOGGER.warn("failed to parse altloc {0} at line {1}, "
                                "residue not present for altloc 'A'".format(
                                repr(key), i+1))
                    continue
                rn = resnames[ids[0]]
                ans = atomnames[ids]
                indices[(ach, ari)] = (rn, ids, ans)
            if rn != arn:
                LOGGER.warn("failed to parse altloc {0} at line {1}, "
                            "residue name mismatch (expected {2}, "
                            "parsed {3})".format(repr(key), i+1, repr(rn),
                                                   repr(arn)))
                continue
            index = ids[(ans == aan).nonzero()[0]]
            if len(index) != 1:
                LOGGER.warn("failed to parse altloc {0} at line {1}, atom"
                            " {2} not found in the residue"
                            .format(repr(key), i+1, repr(aan)))
                continue
            try:
                xyz[index[0], 0] = float(line[30:38])
                xyz[index[0], 1] = float(line[38:46])
                xyz[index[0], 2] = float(line[46:54])
            except:
                LOGGER.warn('failed to parse altloc {0} at line {1}, could'
                            ' not read coordinates'.format(repr(key), i+1))
                continue
            success += 1
        LOGGER.info('{0} out of {1} altloc {2} lines were parsed.'
                    .format(success, len(lines), repr(key)))
        if success > 0:
            LOGGER.info('Altloc {0} is appended as a coordinate set to '
                        'atomgroup {1}.'.format(repr(key), atomgroup.getTitle()))
            atomgroup.addCoordset(xyz, label='altloc ' + key)


#HELIXLINE = ('HELIX  %3d %3s %-3s %1s %4d%1s %-3s %1s %4d%1s%2d'
#             '                               %5d\n')

HELIXLINE = ('HELIX  {serNum:3d} {helixID:>3s} '
             '{initResName:<3s} {initChainID:1s} {initSeqNum:4d}{initICode:1s} '
             '{endResName:<3s} {endChainID:1s} {endSeqNum:4d}{endICode:1s}'
             '{helixClass:2d}                               {length:5d}\n')

SHEETLINE = ('SHEET  {strand:3d} {sheetID:>3s}{numStrands:2d} '
             '{initResName:3s} {initChainID:1s}{initSeqNum:4d}{initICode:1s} '
             '{endResName:3s} {endChainID:1s}{endSeqNum:4d}{endICode:1s}{sense:2d} \n')

PDBLINE_LT100K = ('%-6s%5d %-4s%1s%-3s%2s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '%4s%2s%2s\n')

# Residue number
PDBLINE_GE10K = ('%-6s%5d %-4s%1s%-3s%2s%4x%1s   '
                 '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                 '%4s%2s%2s\n')

# Serial number
PDBLINE_GE100K = ('%-6s%5x %-4s%1s%-3s%2s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '%4s%2s%2s\n')

# Both
PDBLINE_GE100K_GE10K = ('%-6s%5x %-4s%1s%-3s%2s%4x%1s   '
                        '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                        '%4s%2s%2s\n')

# All cases
PDBLINE_H36 = ('%-6s%5s %-4s%1s%-3s%2s%4s%1s   '
               '%8.3f%8.3f%8.3f%6.2f%6.2f      '
               '%4s%2s%2s\n')

ANISOULINE_LT100K = ('%-6s%5d %-4s%1s%-3s%2s%4d%1s '
                     '%7d%7d%7d%7d%7d%7d  '
                     '%4s%2s%2s\n')

# Residue number
ANISOULINE_GE10K = ('%-6s%5d %-4s%1s%-3s%2s%4x%1s '
                    '%7d%7d%7d%7d%7d%7d  '
                    '%4s%2s%2s\n')

# Serial number
ANISOULINE_GE100K = ('%-6s%5x %-4s%1s%-3s%2s%4d%1s '
                     '%7d%7d%7d%7d%7d%7d  '
                     '%4s%2s%2s\n')

# Both
ANISOULINE_GE100K_GE10K = ('%-6s%5x %-4s%1s%-3s%2s%4x%1s '
                           '%7d%7d%7d%7d%7d%7d  '
                           '%4s%2s%2s\n')

# All cases
ANISOULINE_H36 = ('%-6s%5s %-4s%1s%-3s%2s%4s%1s '
                  '%7d%7d%7d%7d%7d%7d  '
                  '%4s%2s%2s\n')

_writePDBdoc = """

    :arg atoms: an object with atom and coordinate data

    :arg csets: coordinate set indices, default is all coordinate sets

    :arg beta: a list or array of number to be outputted in beta column

    :arg occupancy: a list or array of number to be outputted in occupancy
        column

    :arg hybrid36: whether to use hybrid36 format for atoms with serial
        greater than 99999. Hexadecimal is used otherwise.
        Default is False
    :type hybrid36: bool 
    """

def writeChainsList(chains, filename):
    """
    Write a text file containing a list of chains that can be parsed.

    :arg chains: a list of :class:`.Chain` objects
    :type chains: list

    :arg filename: the name of the file to be written
    :type filename: str
    """

    fo = open(filename,'w')
    for chain in chains:
        fo.write(chain.getTitle() + ' ' + chain.getChid() + '\n')
    fo.close()
    return

def parseChainsList(filename):
    """
    Parse a set of PDBs and extract chains based on a list in a text file.

    :arg filename: the name of the file to be read
    :type filename: str

    Returns: lists containing an :class:'.AtomGroup' for each PDB, 
    the headers for those PDBs, and the requested :class:`.Chain` objects
    """

    fi = open(filename,'r')
    lines = fi.readlines()
    fi.close()

    pdb_ids = []
    ags = []
    headers = []
    chains = []
    num_lines = len(lines)
    LOGGER.progress('Starting', num_lines, '_prody_parseChainsList')
    for i, line in enumerate(lines):
        LOGGER.update(i, 'Parsing lines...', label='_prody_parseChainsList')
        pdb_id = line.split()[0].split('_')[0]
        if not pdb_id in pdb_ids:
            pdb_ids.append(pdb_id)

            ag, header = parsePDB(pdb_id, compressed=False, \
                                  subset=line.split()[0].split('_')[1], header=True)

            ags.append(ag)
            headers.append(header)

        chains.append(ag.getHierView()[line.strip().split()[1]])

    LOGGER.finish()
    LOGGER.info('{0} PDBs have been parsed and {1} chains have been extracted. \
                '.format(len(ags),len(chains)))

    return ags, headers, chains

def writePDBStream(stream, atoms, csets=None, **kwargs):
    """Write *atoms* in PDB format to a *stream*.

    :arg stream: anything that implements a :meth:`write` method (e.g. file,
        buffer, stdout)

    :arg renumber: whether to renumber atoms with serial indices
        Default is **True**
    :type renumber: bool
    """

    renumber = kwargs.get('renumber', True)

    remark = str(atoms)
    try:
        coordsets = atoms._getCoordsets(csets)
    except AttributeError:
        try:
            coordsets = atoms._getCoords()
        except AttributeError:
            raise TypeError('atoms must be an object with coordinate sets')
        if coordsets is not None:
            coordsets = [coordsets]
    else:
        if coordsets.ndim == 2:
            coordsets = [coordsets]
    if coordsets is None:
        raise ValueError('atoms does not have any coordinate sets')

    try:
        acsi = atoms.getACSIndex()
    except AttributeError:
        try:
            atoms = atoms.getAtoms()
        except AttributeError:
            raise TypeError('atoms must be an Atomic instance or an object '
                            'with `getAtoms` method')
        else:
            if atoms is None:
                raise ValueError('atoms is not associated with an Atomic '
                                 'instance')
            try:
                acsi = atoms.getACSIndex()
            except AttributeError:
                raise TypeError('atoms does not have a valid type')

    try:
        atoms.getIndex()
    except AttributeError:
        pass
    else:
        atoms = atoms.select('all')

    n_atoms = atoms.numAtoms()

    hybrid36 = kwargs.get('hybrid36', False)

    occupancy = kwargs.get('occupancy')
    if occupancy is None:
        occupancies = atoms._getOccupancies()
        if occupancies is None:
            occupancies = np.zeros(n_atoms, float)
    else:
        occupancies = np.array(occupancy)
        if len(occupancies) != n_atoms:
            raise ValueError('len(occupancy) must be equal to number of atoms')

    beta = kwargs.get('beta')
    if beta is None:
        bfactors = atoms._getBetas()
        if bfactors is None:
            bfactors = np.zeros(n_atoms, float)
    else:
        bfactors = np.array(beta)
        if len(bfactors) != n_atoms:
            raise ValueError('len(beta) must be equal to number of atoms')

    atomnames = atoms.getNames()
    if atomnames is None:
        raise ValueError('atom names are not set')
    for i, an in enumerate(atomnames):
        if len(an) < 4:
            atomnames[i] = ' ' + an

    s_or_u = np.array(['a']).dtype.char

    altlocs = atoms._getAltlocs()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, s_or_u + '1')

    resnames = atoms._getResnames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms

    chainids = atoms._getChids()
    if chainids is None:
        chainids = np.zeros(n_atoms, s_or_u + '1')

    resnums = atoms._getResnums()
    if resnums is None:
        resnums = np.ones(n_atoms, int)

    serials = atoms._getSerials()
    if serials is None or renumber:
        serials = np.arange(n_atoms, dtype=int) + 1

    icodes = atoms._getIcodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, s_or_u + '1')

    hetero = ['ATOM'] * n_atoms
    heteroflags = atoms._getFlags('hetatm')
    if heteroflags is None:
        heteroflags = atoms._getFlags('hetero')
    if heteroflags is not None:
        hetero = np.array(hetero, s_or_u + '6')
        hetero[heteroflags] = 'HETATM'

    elements = atoms._getElements()
    if elements is None:
        elements = np.zeros(n_atoms, s_or_u + '1')
    else:
        elements = np.char.rjust(elements, 2)

    segments = atoms._getSegnames()
    if segments is None:
        segments = np.zeros(n_atoms, s_or_u + '6')

    charges = atoms._getCharges()
    charges2 = np.empty(n_atoms, s_or_u + '2')
    if charges is not None:
        for i, charge in enumerate(charges):
            charges2[i] = str(abs(int(charge)))

            if np.sign(charge) == -1:
                charges2[i] += '-'
            else:
                charges2[i] += '+'

            if charges2[i] == '0+':
                charges2[i] = '  '

    anisous = atoms._getAnisous()
    if anisous is not None:
        anisous = np.array(anisous * 10000, dtype=int)

    # write remarks
    stream.write('REMARK {0}\n'.format(remark))

    # write secondary structures (if any)
    secondary = kwargs.get('secondary', True)
    secstrs = atoms._getSecstrs()
    if secstrs is not None and secondary:
        secindices = atoms._getSecindices()
        secclasses = atoms._getSecclasses()
        secids = atoms._getSecids()

        # write helices
        for i in range(1,max(secindices)+1):
            torf = np.logical_and(isHelix(secstrs), secindices==i)
            if torf.any():
                helix_resnums = resnums[torf]
                helix_chainids = chainids[torf]
                helix_resnames = resnames[torf]
                helix_secclasses = secclasses[torf]
                helix_secids = secids[torf]
                helix_icodes = icodes[torf]
                L = helix_resnums[-1] - helix_resnums[0] + 1

                stream.write(HELIXLINE.format(serNum=i, helixID=helix_secids[0], 
                            initResName=helix_resnames[0], initChainID=helix_chainids[0], 
                            initSeqNum=helix_resnums[0], initICode=helix_icodes[0],
                            endResName=helix_resnames[-1], endChainID=helix_chainids[-1], 
                            endSeqNum=helix_resnums[-1], endICode=helix_icodes[-1],
                            helixClass=helix_secclasses[0], length=L))

        # write strands
        torf_all_sheets = isSheet(secstrs)
        sheet_secids = secids[torf_all_sheets]

        for sheet_id in np.unique(sheet_secids):
            torf_strands_in_sheet = np.logical_and(torf_all_sheets, secids==sheet_id)
            strand_indices = secindices[torf_strands_in_sheet]
            numStrands = len(np.unique(strand_indices))

            for i in np.unique(strand_indices):
                torf_strand = np.logical_and(torf_strands_in_sheet, secindices==i)
                strand_resnums = resnums[torf_strand]
                strand_chainids = chainids[torf_strand]
                strand_resnames = resnames[torf_strand]
                strand_secclasses = secclasses[torf_strand]
                strand_icodes = icodes[torf_strand]

                stream.write(SHEETLINE.format(strand=i, sheetID=sheet_id, numStrands=numStrands,
                            initResName=strand_resnames[0], initChainID=strand_chainids[0], 
                            initSeqNum=strand_resnums[0], initICode=strand_icodes[0],
                            endResName=strand_resnames[-1], endChainID=strand_chainids[-1], 
                            endSeqNum=strand_resnums[-1], endICode=strand_icodes[-1],
                            sense=strand_secclasses[0]))
        pass

    # write atoms
    multi = len(coordsets) > 1
    write = stream.write
    for m, coords in enumerate(coordsets):
        if multi:
            write('MODEL{0:9d}\n'.format(m+1))

        if not hybrid36:
            # We need to check whether serial and residue numbers become hexadecimal
            reached_max_n_atom = False
            reached_max_n_res = False

            pdbline = PDBLINE_LT100K
            anisouline = ANISOULINE_LT100K
        else:
            warned_hybrid36 = False

        warned_5_digit = False

        for i, xyz in enumerate(coords):
            if hybrid36:
                pdbline = PDBLINE_H36
                anisouline = ANISOULINE_H36

                if not warned_hybrid36:
                    LOGGER.warn('hybrid36 format is being used')
                    warned_hybrid36 = True

            else:
                if not (reached_max_n_atom or reached_max_n_res) and (i == MAX_N_ATOM or serials[i] > MAX_N_ATOM):
                    reached_max_n_atom = True
                    pdbline = PDBLINE_GE100K
                    anisouline = ANISOULINE_GE100K
                    LOGGER.warn('Indices are exceeding 99999 and hexadecimal format is being used for indices')

                elif not (reached_max_n_atom or reached_max_n_res) and resnums[i] > MAX_N_RES:
                    reached_max_n_res = True
                    pdbline = PDBLINE_GE10K
                    anisouline = ANISOULINE_GE10K
                    LOGGER.warn('Resnums are exceeding 9999 and hexadecimal format is being used for resnums')

                elif reached_max_n_atom and not reached_max_n_res and resnums[i] > MAX_N_RES:
                    reached_max_n_res = True
                    pdbline = PDBLINE_GE100K_GE10K
                    anisouline = ANISOULINE_GE100K_GE10K
                    LOGGER.warn('Resnums are exceeding 9999 and hexadecimal format is being used for indices and resnums')

                elif reached_max_n_res and not reached_max_n_atom and (i == MAX_N_ATOM or serials[i] > MAX_N_ATOM):
                    reached_max_n_atom = True
                    pdbline = PDBLINE_GE100K_GE10K
                    anisouline = ANISOULINE_GE100K_GE10K
                    LOGGER.warn('Indices are exceeding 99999 and hexadecimal format is being used for indices and resnums')

            if hybrid36:
                serial = decToHybrid36(serials[i])
                resnum = decToHybrid36(resnums[i], resnum=True)
            else:
                serial = serials[i]
                resnum = resnums[i]

            if pdbline == PDBLINE_LT100K or hybrid36:
                if len(str(resnum)) == 5:
                    if icodes[i] == '':
                        icodes[i] = str(resnum)[4]

                        if not warned_5_digit:
                            LOGGER.warn('Storing 5-digit resnums using insertion codes')
                            warned_5_digit = True
                    else:
                        LOGGER.warn('Truncating 5-digit resnum as insertion code is busy.')

                    resnum = int(str(resnum)[:4])

                elif len(str(resnum)) > 5:
                    if not warned_5_digit:
                        LOGGER.warn('Truncating {0}-digit resnum as too long to be '
                                    'supported by insertion code.'.format(len(str(resnum))))
                        warned_5_digit = True

                    resnum = int(str(resnum)[:4])
            else:
                final_resnum = '%4x' % int(resnum)

                if len(str(final_resnum)) == 5:
                    if icodes[i] == '':
                        icodes[i] = str(final_resnum)[4]

                        if not warned_5_digit:
                            LOGGER.warn('Storing 5-digit hex resnums using insertion codes')
                            warned_5_digit = True
                    else:
                        LOGGER.warn('Truncating 5-digit hex resnum as insertion code is busy.')

                    resnum = int(str(final_resnum)[:4], 16)

                elif len(str(final_resnum)) > 5:
                    if not warned_5_digit:
                        LOGGER.warn('Truncating {0}-digit hex resnum ({1}) as too long to be '
                                    'supported by insertion code.'.format(len(str(final_resnum)), 
                                                                          final_resnum))
                        warned_5_digit = True

                    resnum = int(str(final_resnum)[:4], 16)


            write(pdbline % (hetero[i], serial,
                             atomnames[i], altlocs[i],
                             resnames[i], chainids[i], resnum,
                             icodes[i],
                             xyz[0], xyz[1], xyz[2],
                             occupancies[i], bfactors[i],
                             segments[i], elements[i], charges2[i]))

            if anisous is not None:
                anisou = anisous[i]

                write(anisouline % ("ANISOU", serial,
                                    atomnames[i], altlocs[i],
                                    resnames[i], chainids[i], resnum,
                                    icodes[i],
                                    anisou[0], anisou[1], anisou[2],
                                    anisou[3], anisou[4], anisou[5],
                                    segments[i], elements[i], charges2[i]))

            if atoms.getFlags('pdbter') is not None and atoms.getFlags('pdbter')[i]:
                write('TER\n')

        if multi:
            write('ENDMDL\n')
            altlocs = np.zeros(n_atoms, s_or_u + '1')

writePDBStream.__doc__ += _writePDBdoc

def writePDB(filename, atoms, csets=None, autoext=True, **kwargs):
    """Write *atoms* in PDB format to a file with name *filename* and return
    *filename*.  If *filename* ends with :file:`.gz`, a compressed file will
    be written.

    :arg renumber: whether to renumber atoms with serial indices
        Default is **True**
    :type renumber: bool
    """

    if not (filename.endswith('.pdb') or filename.endswith('.pdb.gz') or
            filename.endswith('.ent') or filename.endswith('.ent.gz')):
        filename += '.pdb'
    out = openFile(filename, 'wt')
    writePDBStream(out, atoms, csets, **kwargs)
    out.close()
    return filename

writePDB.__doc__ += _writePDBdoc + """
    :arg autoext: when not present, append extension :file:`.pdb` to *filename*
"""

def writePQRStream(stream, atoms, **kwargs):
    if isinstance(atoms, Atom):
        atoms = Selection(atoms.getAtomGroup(), [atoms.getIndex()],
                          atoms.getACSIndex(),
                          'index ' + str(atoms.getIndex()))
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

    s_or_u = np.array(['a']).dtype.char


    write = stream.write

    calphas = atoms.ca
    ssa = calphas.getSecstrs()
    helix = []
    sheet = []
    if ssa is not None:
        ss_prev = ssa[0]
        ss_start = 0
        ss_end = 1
        for i, ss in enumerate(ssa):
            if ss != ss_prev:
                # store prev secstr and prepare for next
                ss_end = i-1
                init = calphas[ss_start]
                end = calphas[ss_end]
                length = ss_end - ss_start + 1

                entries = [init.getSecindex(), init.getSecid(),
                           init.getResname(), init.getChid(), 
                           init.getResnum(), init.getIcode(),
                           end.getResname(), end.getChid(),
                           end.getResnum(), end.getIcode(),
                           init.getSecclass()]

                if ssa[ss_end] == 'H':
                    helix.append(["HELIX "] + entries +
                                 ['', length])

                elif ssa[ss_end] == 'E':
                    sheet.append(["SHEET "] + entries)

                ss_start = i
                ss_prev = ss

    format_helix = ('{0:6s} {1:3d} {2:3s} ' +
                    '{3:3s} {4:1s} {5:4d}{6:1s} ' +
                    '{7:3s} {8:1s} {9:4d}{10:1s} ' +
                    '{11:2d} {12:30s} {13:5d}\n').format
    for line in helix:
        write(format_helix(*line))


    sorted_sheet = sorted(sheet, key=lambda item: (item[2], item[1]))
    sheet_prev = 'A'
    num_strands_list = []
    for i, item1 in enumerate(sorted_sheet):
        if item1[2] != sheet_prev:
            num_strands = sorted_sheet[i-1][1]
            num_strands_list.append(num_strands)

            sheet_prev = item1[2]

            for item2 in sorted_sheet[i-num_strands:i]:
                item2.append(num_strands)

    num_strands = item1[1]
    for item2 in sorted_sheet[i-num_strands+1:]:
        item2.append(num_strands)

    format_sheet = ('{0:6s} {1:3d} {2:3s}{12:2d} ' +
                    '{3:3s} {4:1s}{5:4d}{6:1s}' +
                    '{7:3s} {8:1s}{9:4d}{10:1s}' +
                    '{11:2d}\n').format

    for i, line in enumerate(sorted_sheet):
        write(format_sheet(*line))

    resnames = atoms._getResnames()
    if resnames is None:
        resnames = ['UNK'] * n_atoms
    resnums = atoms._getResnums()
    if resnums is None:
        resnums = np.ones(n_atoms, int)
    chainids = atoms._getChids()
    if chainids is None:
        chainids = np.zeros(n_atoms, s_or_u + '1')
    charges = atoms._getCharges()
    if charges is None:
        charges = np.zeros(n_atoms, float)
    radii = atoms._getRadii()
    if radii is None:
        radii = np.zeros(n_atoms, float)
    icodes = atoms._getIcodes()
    if icodes is None:
        icodes = np.zeros(n_atoms, s_or_u + '1')
    hetero = ['ATOM'] * n_atoms
    heteroflags = atoms._getFlags('hetatm')
    if heteroflags is None:
        heteroflags = atoms._getFlags('hetero')
    if heteroflags is not None:
        hetero = np.array(hetero, s_or_u + '6')
        hetero[heteroflags] = 'HETATM'
    altlocs = atoms._getAltlocs()
    if altlocs is None:
        altlocs = np.zeros(n_atoms, s_or_u + '1')

    format = ('{0:6s} {1:5d} {2:4s} {3:1s}' +
              '{4:4s} {5:1s} {6:4d} {7:1s}   ' +
              '{8:8.3f} {9:8.3f} {10:8.3f}' +
              '{11:8.4f} {12:7.4f}\n').format
    coords = atoms._getCoords()

    for i, xyz in enumerate(coords):
        write(format(hetero[i], i+1, atomnames[i], altlocs[i],
                     resnames[i], chainids[i], int(resnums[i]),
                     icodes[i], xyz[0], xyz[1], xyz[2], charges[i], radii[i]))

def writePQR(filename, atoms, **kwargs):
    """Write *atoms* in PQR format to a file with name *filename*.  Only
    current coordinate set is written.  Returns *filename* upon success.  If
    *filename* ends with :file:`.gz`, a compressed file will be written."""

    stream = openFile(filename, 'w')
    writePQRStream(stream, atoms, **kwargs)
    stream.close()
    return filename
