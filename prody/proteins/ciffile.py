# -*- coding: utf-8 -*-
"""This module defines functions for parsing `mmCIF files`_.

.. _mmCIF files: http://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html"""


from collections import OrderedDict
import os.path
import numpy as np

from prody.atomic import AtomGroup
from prody.atomic import flags
from prody.atomic import ATOMIC_FIELDS
from prody.utilities import openFile
from prody import LOGGER, SETTINGS

from .localpdb import fetchPDB
from .starfile import parseSTARLines, StarDict, parseSTARSection
from .cifheader import getCIFHeaderDict
from .header import buildBiomolecules, assignSecstr, isHelix, isSheet

__all__ = ['parseMMCIFStream', 'parseMMCIF']


class MMCIFParseError(Exception):
    pass


_parseMMCIFdoc = """
    :arg title: title of the :class:`.AtomGroup` instance, default is the
        PDB filename or PDB identifier
    :type title: str

    :arg chain: chain identifiers for parsing specific chains, e.g.
        ``chain='A'``, ``chain='B'``, ``chain='DE'``, by default all
        chains are parsed
    :type chain: str

    :arg subset: a predefined keyword to parse subset of atoms, valid keywords
        are ``'calpha'`` (``'ca'``), ``'backbone'`` (``'bb'``), or **None**
        (read all atoms), e.g. ``subset='bb'``
    :type subset: str

    :arg model: model index or None (read all models), e.g. ``model=10``
    :type model: int, list

    :arg altloc: if a location indicator is passed, such as ``'A'`` or ``'B'``,
         only indicated alternate locations will be parsed as the single
         coordinate set of the AtomGroup, if *altloc* is set **True** all
         alternate locations will be parsed and each will be appended as a
         distinct coordinate set, default is ``"A"``
    :type altloc: str
    """

_PDBSubsets = {'ca': 'ca', 'calpha': 'ca', 'bb': 'bb', 'backbone': 'bb'}


def parseMMCIF(pdb, **kwargs):
    """Returns an :class:`.AtomGroup` and/or a :class:`.StarDict` containing header data
    parsed from an mmCIF file. If not found, the mmCIF file will be downloaded
    from the PDB. It will be downloaded in uncompressed format regardless of
    the compressed keyword.

    This function extends :func:`.parseMMCIFStream`.

    :arg pdb: a PDB identifier or a filename
        If needed, mmCIF files are downloaded using :func:`.fetchPDB()` function.
    :type pdb: str

    :arg chain: comma separated string or list-like of chain IDs
    :type chain: str, tuple, list, :class:`~numpy.ndarray`
    """
    chain = kwargs.pop('chain', None)
    title = kwargs.get('title', None)
    auto_bonds = SETTINGS.get('auto_bonds')
    get_bonds = kwargs.get('bonds', auto_bonds)
    if get_bonds:
        LOGGER.warn('Parsing struct_conn information from mmCIF is current unsupported and no bond information is added to the results')
    if not os.path.isfile(pdb):
        if len(pdb) == 5 and pdb.isalnum():
            if chain is None:
                chain = pdb[-1]
                pdb = pdb[:4]
            else:
                raise ValueError('Please provide chain as a keyword argument '
                                 'or part of the PDB ID, not both')
        else:
            chain = chain

        if len(pdb) == 4 and pdb.isalnum():
            if title is None:
                title = pdb
                kwargs['title'] = title

            if os.path.isfile(pdb + '.cif'):
                filename = pdb + '.cif'
            elif os.path.isfile(pdb + '.cif.gz'):
                filename = pdb + '.cif.gz'
            else:
                filename = fetchPDB(pdb, report=True,
                                    format='cif', compressed=False)
                if filename is None:
                    raise IOError('mmCIF file for {0} could not be downloaded.'
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
    cif = openFile(pdb, 'rt')
    result = parseMMCIFStream(cif, chain=chain, **kwargs)
    cif.close()
    return result


def parseMMCIFStream(stream, **kwargs):
    """Returns an :class:`.AtomGroup` and/or a class:`.StarDict` 
    containing header data parsed from a stream of CIF lines.

    :arg stream: Anything that implements the method ``readlines``
        (e.g. :class:`file`, buffer, stdin)
        
    """

    model = kwargs.get('model')
    subset = kwargs.get('subset')
    chain = kwargs.get('chain')
    altloc = kwargs.get('altloc', 'A')
    header = kwargs.get('header', False)
    assert isinstance(header, bool), 'header must be a boolean'

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
    auto_secondary = SETTINGS.get('auto_secondary')
    secondary = kwargs.get('secondary', auto_secondary)
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
            hd = getCIFHeaderDict(lines)

        _parseMMCIFLines(ag, lines, model, chain, subset, altloc)

        if ag.numAtoms() > 0:
            LOGGER.report('{0} atoms and {1} coordinate set(s) were '
                          'parsed in %.2fs.'.format(ag.numAtoms(),
                                                    ag.numCoordsets() - n_csets))
        else:
            ag = None
            LOGGER.warn('Atomic data could not be parsed, please '
                        'check the input file.')
    elif header:
        hd = getCIFHeaderDict(stream)

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


parseMMCIFStream.__doc__ += _parseMMCIFdoc


def _parseMMCIFLines(atomgroup, lines, model, chain, subset,
                     altloc_torf):
    """Returns an AtomGroup. See also :func:`.parsePDBStream()`.

    :arg lines: mmCIF lines
    """

    if subset is not None:
        if subset == 'ca':
            subset = set(('CA',))
        elif subset in 'bb':
            subset = flags.BACKBONE
        protein_resnames = flags.AMINOACIDS

    asize = 0
    i = 0
    models = []
    nModels = 0
    fields = OrderedDict()
    fieldCounter = -1
    foundAtomBlock = False
    doneAtomBlock = False
    start = 0
    stop = 0
    while not doneAtomBlock:
        line = lines[i]
        if line[:11] == '_atom_site.':
            fieldCounter += 1
            fields[line.split('.')[1].strip()] = fieldCounter

        if line.startswith('ATOM ') or line.startswith('HETATM'):
            if not foundAtomBlock:
                foundAtomBlock = True
                start = i
            models.append(line.split()[fields['pdbx_PDB_model_num']])
            if len(models) == 1 or (models[asize] != models[asize-1]):
                nModels += 1
            asize += 1
        else:
            if foundAtomBlock:
                doneAtomBlock = True
                stop = i
        i += 1

    new_start = start
    new_stop = stop
    if model is not None:
        startp1 = start + 1
        for i in range(start, stop):
            if int(models[i-startp1]) != model and int(models[i+1-startp1]) == model:
                new_start = i
            if model != nModels and (int(models[i-startp1]) == model and int(models[i+1-startp1]) != model):
                new_stop = i
                break
        if not str(model) in models:
            raise MMCIFParseError('model {0} is not found'.format(model))

    start = new_start
    stop = new_stop

    addcoords = False
    if atomgroup.numCoordsets() > 0:
        addcoords = True

    if isinstance(altloc_torf, str):
        if altloc_torf == 'all':
            which_altlocs = 'all'
        elif altloc_torf.strip() != 'A':
            LOGGER.info('Parsing alternate locations {0}.'
                        .format(altloc_torf))
            which_altlocs = '.' + ''.join(altloc_torf.split())
        else:
            which_altlocs = '.A'
        altloc_torf = False
    else:
        which_altlocs = '.A'
        altloc_torf = True

    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
    hetero = np.zeros(asize, dtype=bool)
    termini = np.zeros(asize, dtype=bool)
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)
    bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
    occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)

    n_atoms = atomgroup.numAtoms()
    if n_atoms > 0:
        asize = n_atoms

    acount = 0
    for line in lines[start:stop]:
        startswith = line.split()[fields['group_PDB']]

        atomname = line.split()[fields['auth_atom_id']]
        if atomname.startswith('"') and atomname.endswith('"'):
            atomname = atomname[1:-1]
        resname = line.split()[fields['auth_comp_id']]

        if subset is not None:
            if not (atomname in subset and resname in protein_resnames):
                continue

        chID = line.split()[fields['auth_asym_id']]
        if chain is not None:
            if isinstance(chain, str):
                chain = chain.split(',')
            if not chID in chain:
                continue

        segID = line.split()[fields['label_asym_id']]

        alt = line.split()[fields['label_alt_id']]
        if alt not in which_altlocs and which_altlocs != 'all':
            continue

        if alt == '.':
            alt = ' '

        coordinates[acount] = [line.split()[fields['Cartn_x']],
                               line.split()[fields['Cartn_y']],
                               line.split()[fields['Cartn_z']]]
        atomnames[acount] = atomname
        resnames[acount] = resname
        resnums[acount] = line.split()[fields['auth_seq_id']]
        chainids[acount] = chID
        segnames[acount] = segID
        hetero[acount] = startswith == 'HETATM' # True or False

        if chainids[acount] != chainids[acount-1]: 
            termini[acount-1] = True

        altlocs[acount] = alt
        
        try:
            icodes[acount] = line.split()[fields['pdbx_PDB_ins_code']]
        except KeyError:
            icodes[acount] = ''

        if icodes[acount] == '?': 
            icodes[acount] = ''

        serials[acount] = line.split()[fields['id']]
        elements[acount] = line.split()[fields['type_symbol']]
        bfactors[acount] = line.split()[fields['B_iso_or_equiv']]
        occupancies[acount] = line.split()[fields['occupancy']]

        acount += 1

    if model is None:
        modelSize = acount//nModels
    else:
        modelSize = acount

    if addcoords:
        atomgroup.addCoordset(coordinates[:modelSize])
    else:
        atomgroup._setCoords(coordinates[:modelSize])

    atomgroup.setNames(atomnames[:modelSize])
    atomgroup.setResnames(resnames[:modelSize])
    atomgroup.setResnums(resnums[:modelSize])
    atomgroup.setSegnames(segnames[:modelSize])
    atomgroup.setChids(chainids[:modelSize])
    atomgroup.setFlags('hetatm', hetero[:modelSize])
    atomgroup.setFlags('pdbter', termini[:modelSize])
    atomgroup.setAltlocs(altlocs[:modelSize])
    atomgroup.setIcodes(icodes[:modelSize])
    atomgroup.setSerials(serials[:modelSize])

    atomgroup.setElements(elements[:modelSize])
    from prody.utilities.misctools import getMasses
    atomgroup.setMasses(getMasses(elements[:modelSize]))
    atomgroup.setBetas(bfactors[:modelSize])
    atomgroup.setOccupancies(occupancies[:modelSize])

    anisou = None
    siguij = None
    try:
        data = parseSTARSection(lines, "_atom_site_anisotrop")
        x = data[0] # check if data has anything in it
    except IndexError:
        LOGGER.warn("No anisotropic B factors found")
    else:
        anisou = np.zeros((acount, 6),
                          dtype=ATOMIC_FIELDS['anisou'].dtype)
        
        if "_atom_site_anisotrop.U[1][1]_esd" in data[0].keys():
            siguij = np.zeros((acount, 6),
                              dtype=ATOMIC_FIELDS['siguij'].dtype)

        for entry in data:
            try:
                index = np.where(atomgroup.getSerials() == int(
                    entry["_atom_site_anisotrop.id"]))[0][0]
            except:
                continue
            
            anisou[index, 0] = entry['_atom_site_anisotrop.U[1][1]']
            anisou[index, 1] = entry['_atom_site_anisotrop.U[2][2]']
            anisou[index, 2] = entry['_atom_site_anisotrop.U[3][3]']
            anisou[index, 3] = entry['_atom_site_anisotrop.U[1][2]']
            anisou[index, 4] = entry['_atom_site_anisotrop.U[1][3]']
            anisou[index, 5] = entry['_atom_site_anisotrop.U[2][3]'] 

            if siguij is not None:
                try:
                    siguij[index, 0] = entry['_atom_site_anisotrop.U[1][1]_esd']
                    siguij[index, 1] = entry['_atom_site_anisotrop.U[2][2]_esd']
                    siguij[index, 2] = entry['_atom_site_anisotrop.U[3][3]_esd']
                    siguij[index, 3] = entry['_atom_site_anisotrop.U[1][2]_esd']
                    siguij[index, 4] = entry['_atom_site_anisotrop.U[1][3]_esd']
                    siguij[index, 5] = entry['_atom_site_anisotrop.U[2][3]_esd']
                except:
                    pass

        atomgroup.setAnisous(anisou) # no division needed anymore

        if not np.any(siguij):
            atomgroup.setAnistds(siguij)  # no division needed anymore

    if model is None:
        for n in range(1, nModels):
            atomgroup.addCoordset(coordinates[n*modelSize:(n+1)*modelSize])

    return atomgroup
