# -*- coding: utf-8 -*-
"""This module defines functions for parsing and writing `PDB files`_.

.. _PDB files: http://www.wwpdb.org/documentation/format32/v3.2.html"""

from collections import defaultdict
import os.path


import numpy as np

from prody.atomic import AtomGroup
from prody.atomic import flags
from prody.atomic import ATOMIC_FIELDS
from prody.utilities import openFile
from prody import LOGGER, SETTINGS

from .header import getHeaderDict, buildBiomolecules, assignSecstr, isHelix, isSheet
from .localpdb import fetchPDB

__all__ = ['parsePDBStream', 'parsePDB', 'parsePQR',
           'writePDBStream', 'writePDB',]

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

    :arg biomol: if **True**, biomolecule obtained by transforming the
        coordinates using information from header section will be returned
    :type biomol: False

    :arg secondary: if **True**, secondary structure information from header
        section will be assigned atoms
    :type secondary: False

    If ``model=0`` and ``header=True``, return header dictionary only.

    Note that this function does not evaluate ``CONECT`` records.

    """

_PDBSubsets = {'ca': 'ca', 'calpha': 'ca', 'bb': 'bb', 'backbone': 'bb'}

def parsePDB(pdb, **kwargs):
    """Returns an :class:`.AtomGroup` and/or dictionary containing header data
    parsed from a PDB file.

    This function extends :func:`.parsePDBStream`.

    See :ref:`parsepdb` for a detailed usage example.

    :arg pdb: a PDB identifier or a filename
        If needed, PDB files are downloaded using :func:`.fetchPDB()` function.
        You can also provide arguments that you would like passed on to fetchPDB().
    """
    title = kwargs.get('title', None)
    if not os.path.isfile(pdb):
        if len(pdb) == 4 and pdb.isalnum():
            if title is None:
                title = pdb
                kwargs['title'] = title
            filename = fetchPDB(pdb, report=True, **kwargs)
            if filename is None:
                raise IOError('PDB file for {0} could not be downloaded.'
                              .format(pdb))
            pdb = filename
        else:
            raise IOError('{0} is not a valid filename or a valid PDB '
                          'identifier.'.format(pdb))
    if title is None:
        title, ext = os.path.splitext(os.path.split(pdb)[1])
        if ext == '.gz':
            title, ext = os.path.splitext(title)
        if len(title) == 7 and title.startswith('pdb'):
            title = title[3:]
        kwargs['title'] = title
    pdb = openFile(pdb, 'rt')
    result = parsePDBStream(pdb, **kwargs)
    pdb.close()
    return result

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
    if model is not None:
        if isinstance(model, int):
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
        title_suffix = '_' + chain + title_suffix
    ag = None
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    elif model != 0:
        ag = AtomGroup(str(kwargs.get('title', 'Unknown')) + title_suffix)
        n_csets = 0

    biomol = kwargs.get('biomol', False)
    auto_secondary = None
    secondary = kwargs.get('secondary')
    if not secondary:
        auto_secondary = SETTINGS.get('auto_secondary')
        secondary = auto_secondary
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
        _parsePDBLines(ag, lines, split, model, chain, subset, altloc)
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
    model = 1
    header = False
    chain = kwargs.get('chain')
    subset = kwargs.get('subset')
    altloc = kwargs.get('altloc', 'A')
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
                   altloc_torf, format='PDB'):
    """Returns an AtomGroup. See also :func:`.parsePDBStream()`.

    :arg lines: PDB/PQR lines
    :arg split: starting index for coordinate data lines"""

    format = format.upper()
    if format == 'PDB':
        isPDB = True
    else:
        isPDB = False

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
    hetero = np.zeros(asize, dtype=bool)
    termini = np.zeros(asize, dtype=bool)
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    if isPDB:
        segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
        elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)
        bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
        occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)
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
            raise PDBParseError('model {0} is not found'.format(model))
    if isinstance(altloc_torf, str):
        if altloc_torf.strip() != 'A':
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
                    LOGGER.warn('Discarding model {0}, which contains more '
                            'atoms than first model does.'.format(nmodel+1))
                    acount = 0
                    nmodel += 1
                    coordinates = np.zeros((n_atoms, 3), dtype=float)
                    while lines[i][:6] != 'ENDMDL':
                        i += 1
                else:
                    raise PDBParseError('invalid or missing coordinate(s) at '
                                         'line {0}.'.format(i+1))
            if onlycoords:
                acount += 1
                i += 1
                continue

            try:
                serials[acount] = line[6:11]
            except ValueError:
                try:
                    serials[acount] = int(line[6:11], 16)
                except ValueError:
                    LOGGER.warn('Failed to parse serial number in line {0}.'
                                .format(i))
                    serials[acount] = serials[acount-1]+1
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
            else:
                try:
                    charges[acount] = line[54:62]
                except:
                    LOGGER.warn('failed to parse charge at line {0}'
                                .format(i))
                try:
                    radii[acount] = line[62:69]
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
        #elif startswith == 'END   ' or startswith == 'CONECT':
        #    i += 1
        #    break
        elif not onlycoords and (startswith == 'TER   ' or
            startswith.strip() == 'TER'):
            termini[acount - 1] = True
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
                    charges.resize(acount, refcheck=False)
                    radii.resize(acount, refcheck=False)
                    atomgroup.setCharges(charges)
                    atomgroup.setRadii(radii)

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
            charges.resize(acount, refcheck=False)
            radii.resize(acount, refcheck=False)
            atomgroup.setCharges(charges)
            atomgroup.setRadii(radii)

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
            LOGGER.info('Altloc {0} is appended as a coordinate set to the '
                        'atom group.'.format(repr(key), atomgroup.getTitle()))
            atomgroup.addCoordset(xyz, label='altloc ' + key)

PDBLINE = ('{0:6s}{1:5d} {2:4s}{3:1s}'
           '{4:4s}{5:1s}{6:4d}{7:1s}   '
           '{8:8.3f}{9:8.3f}{10:8.3f}'
           '{11:6.2f}{12:6.2f}      '
           '{13:4s}{14:2s}\n')

HELIXLINE = ('%-6s %3d %-3s %-3s %1s %4d%1s %-3s %1s %4d%1s%2d'
             '                              %5d\n')

PDBLINE_LT100K = ('%-6s%5d %-4s%1s%-4s%1s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '%4s%2s\n')

PDBLINE_GE100K = ('%-6s%5x %-4s%1s%-4s%1s%4d%1s   '
                  '%8.3f%8.3f%8.3f%6.2f%6.2f      '
                  '%4s%2s\n')


_writePDBdoc = """

    :arg atoms: an object with atom and coordinate data

    :arg csets: coordinate set indices, default is all coordinate sets

    :arg beta: a list or array of number to be outputted in beta column

    :arg occupancy: a list or array of number to be outputted in occupancy
        column
    """

def writePDBStream(stream, atoms, csets=None, **kwargs):
    """Write *atoms* in PDB format to a *stream*.

    :arg stream: anything that implements a :meth:`write` method (e.g. file,
        buffer, stdout)"""

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

    # write remarks
    stream.write('REMARK {0}\n'.format(remark))

    # write secondary structures (if any)
    secstrs = atoms._getSecstrs()
    if secstrs is not None:
        secindices = atoms._getSecindices()
        secclasses = atoms._getSecclasses()
        secids = atoms._getSecids()

        # write helices
        for i in range(1,max(secindices)+1):
            torf = np.logical_or(isHelix(secstrs), secindices==i)
            helix_resnums = resnums[torf]
            helix_chainids = chainids[torf]
            helix_resname = resnames[torf]

        pass

    # write atoms
    multi = len(coordsets) > 1
    write = stream.write
    for m, coords in enumerate(coordsets):
        pdbline = PDBLINE_LT100K
        if multi:
            write('MODEL{0:9d}\n'.format(m+1))
        for i, xyz in enumerate(coords):
            if i == 99999:
                pdbline = PDBLINE_GE100K
            write(pdbline % (hetero[i], i+1,
                         atomnames[i], altlocs[i],
                         resnames[i], chainids[i], resnums[i],
                         icodes[i],
                         xyz[0], xyz[1], xyz[2],
                         occupancies[i], bfactors[i],
                         segments[i], elements[i]))
        if multi:
            write('ENDMDL\n')
            altlocs = np.zeros(n_atoms, s_or_u + '1')

writePDBStream.__doc__ += _writePDBdoc

def writePDB(filename, atoms, csets=None, autoext=True, **kwargs):
    """Write *atoms* in PDB format to a file with name *filename* and return
    *filename*.  If *filename* ends with :file:`.gz`, a compressed file will
    be written."""

    if not ('.pdb' in filename or '.pdb.gz' in filename or
             '.ent' in filename or '.ent.gz' in filename):
        filename += '.pdb'
    out = openFile(filename, 'wt')
    writePDBStream(out, atoms, csets, **kwargs)
    out.close()
    return filename

writePDB.__doc__ += _writePDBdoc + """
    :arg autoext: when not present, append extension :file:`.pdb` to *filename*
"""
