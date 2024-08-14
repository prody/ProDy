# -*- coding: utf-8 -*-
"""This module defines functions for parsing MMTF files and structure objects.

.. MMTF files: https://mmtf.rcsb.org"""

from numbers import Number
import os.path
from collections import OrderedDict

from prody.atomic import AtomGroup
from prody.atomic import flags
from prody.atomic import ATOMIC_FIELDS

from prody.utilities import openFile, isListLike, copy
from prody import LOGGER
from prody.utilities.misctools import getMasses
from .header import getHeaderDict, buildBiomolecules

import struct as st
import numpy as np

__all__ = ['parseMMTF', 'writeMMTF']

_parseMMTFdoc = """
    :arg chain: Chain identifier(s) to parse (e.g., 'A' for chain A). If not provided,
        all chains are parsed. If a PDB ID is used, chain can also be specified as part
        of the ID (e.g., '1abcA' to parse chain A).
    :type chain: str, optional

    :arg title: Title to assign to the resulting AtomGroup. If not provided, the
        title is extracted from the MMTF structure or set to the PDB ID.
    :type title: str, optional

    :arg subset: a predefined keyword to parse subset of atoms, valid keywords
        are ``'calpha'`` (``'ca'``), ``'backbone'`` (``'bb'``), or **None**
        (read all atoms), e.g. ``subset='bb'``
    :type subset: str

    :arg model: model index or None (read all models), e.g. ``model=10``
    :type model: int, list

    :arg altloc: if a location indicator is passed, such as ``'A'`` or ``'B'``,
         only indicated alternate locations will be parsed as the single
         coordinate set of the AtomGroup. If *altloc* is ``'all'`` then all
         alternate locations will be parsed and each will be appended as a
         distinct coordinate set, default is ``"A"``.  In the rare instance
         where all atoms have a location indicator specified and this does not
         match altloc, the first location indicator in the file is used.
    :type altloc: str
    """

def parseMMTF(mmtf_struc, **kwargs):
    """
    Parse an MMTF (Macromolecular Transmission Format) structure or fetch it from the PDB,
    and return an AtomGroup containing the parsed data.

    :param mmtf_struc: The MMTF structure to parse. It can be provided in one of the following ways:
        - A string representing a PDB ID (e.g., '1abc').
        - The filename of an MMTF file (ending with '.mmtf' or '.mmtf.gz').
        - An MMTF structure object (file parsed through mmtf-python).
    :type mmtf_struc: str or MMTF Structure object

    :param chain: Chain identifier(s) to parse (e.g., 'A' for chain A). If not provided,
        all chains are parsed. If a PDB ID is used, chain can also be specified as part
        of the ID (e.g., '1abcA' to parse chain A).
    :type chain: str, optional

    :return: An AtomGroup containing the parsed atomic data.
    :rtype: AtomGroup
    """
    try:
        from mmtf import fetch, parse
        import gzip, msgpack
        import mmtf
    except ImportError:
        raise ImportError("Install mmtf to read in mmtf structure objects (e.g. pip install mmtf-python)")

    chain = kwargs.pop('chain', None)
    title = kwargs.get('title', None)

    if type(mmtf_struc)==str:
        if not mmtf_struc.endswith(".mmtf") and not mmtf_struc.endswith(".mmtf.gz"): #entering just the pdb id
            if len(mmtf_struc) == 5 and mmtf_struc.isalnum():
                if chain is None:
                    chain = mmtf_struc[-1]
                    mmtf_struc = mmtf_struc[:4]
                else:
                    raise ValueError('Please provide chain as a keyword argument '
                                    'or part of the PDB ID, not both')
            else:
                chain = chain

            if len(mmtf_struc) == 4 and mmtf_struc.isalnum():
                if title is None:
                    title = mmtf_struc
                    kwargs['title'] = title
                structure = fetch(mmtf_struc)
                if structure is None:
                    raise IOError('mmtf file for {0} could not be downloaded.'
                                    .format(mmtf_struc))
                mmtf_struc = structure
            else:
                raise IOError('{0} is not a valid mmtf filename or a valid PDB '
                            'identifier.'.format(mmtf_struc))
            
        else: #entering the mmtf file name
            structure=parse(mmtf_struc)
            if structure is None:
                raise IOError('mmtf file for {0} could not be downloaded.'
                                .format(mmtf_struc))
            mmtf_struc = structure
            title = mmtf_struc.structure_id  
    elif type(mmtf_struc)==bytearray:
        data = gzip.decompress(mmtf_struc)
        unpack = msgpack.unpackb(data)
        dec = mmtf.MMTFDecoder()
        dec.decode_data(unpack)
        mmtf_struc = dec


    #if none of the above loops are entered, user should have passed a mmtf structure object
    if title is None:
        title = mmtf_struc.structure_id

    result = _parseMMTF(mmtf_struc, chain=chain, **kwargs)

    return result

parseMMTF.__doc__ += _parseMMTFdoc

def _has_biomol(hd):
    '''Return true if header actually has multiple transformations'''
    if 'biomoltrans' not in hd:
        return False
    for v in hd['biomoltrans'].values():
        if len(v) > 4:
            return True
    return False
    
def _parseMMTF(mmtf_struc, **kwargs):
    LOGGER.timeit()
    ag = AtomGroup()
    model = kwargs.get('model')
    subset = kwargs.get('subset')
    chain = kwargs.get('chain')
    header = kwargs.get('header', False)
    get_bonds = kwargs.get('bonds',False) 
    altloc_sel = kwargs.get('altloc', None)
    
    assert isinstance(header, bool), 'header must be a boolean'

    if model is not None:
        if isinstance(model, int):
            if model < 0:
                raise ValueError('model must be greater than 0')
        else:
            raise TypeError('model must be an integer, {0} is invalid'
                            .format(str(model)))
    title_suffix = ''

    biomol = kwargs.get('biomol', False)
    hd = set_header(mmtf_struc)
    ag = set_info(ag, mmtf_struc, get_bonds, altloc_sel)

    if ag.numAtoms() > 0:
            LOGGER.report('{0} atoms and {1} coordinate set(s) were '
                          'parsed in %.2fs.'.format(ag.numAtoms(),
                                                    ag.numCoordsets()))
    else:
        ag = None
        LOGGER.warn('Atomic data could not be parsed, please '
                    'check the input file.')
        
    
        
    if ag is not None and isinstance(hd, dict):
        if biomol and _has_biomol(hd):
            ag = buildBiomolecules(hd, ag)
            if isinstance(ag, list):
                LOGGER.info('Biomolecular transformations were applied, {0} '
                            'biomolecule(s) are returned.'.format(len(ag)))
            else:
                LOGGER.info('Biomolecular transformations were applied to the '
                            'coordinate data.')    

    if header:
        return ag, hd
    else:
        return ag


def set_header(data_api):

    #get the transform here and convert it to the format that prody wants
    chain_list = data_api.chain_name_list
    assembly = _bio_transform(data_api)
    header = {
    'r_work': data_api.r_work,
    'r_free': data_api.r_free,
    'resolution': data_api.resolution,
    'title': data_api.title,
    'deposition_date': data_api.deposition_date,
    'release_date': data_api.release_date,
    'experimental_methods': data_api.experimental_methods,
    'biomoltrans': assembly
    }

    return header

def _bio_transform(dec):
    ret = {}
    matrix_line_format_string = ' %9.6f %9.6f %9.6f      %9.5f            \n'
    for t in dec.bio_assembly:
        L = []
        for trans in t['transformList']:
            chis = trans['chainIndexList']
            chains = sorted(set([dec.chain_name_list[c] for c in chis]))
            m = trans['matrix']
            L += [chains,
                  matrix_line_format_string % tuple(m[:4]),
                  matrix_line_format_string % tuple(m[4:8]),
                  matrix_line_format_string % tuple(m[8:12])]
        ret[t['name']] = L
    return ret

def set_info(atomgroup, mmtf_data, get_bonds=False, altloc_sel=None):

    mmtfHETATMtypes = set([
        "D-SACCHARIDE",
        "D-SACCHARIDE 1,4 AND 1,4 LINKING",
        "D-SACCHARIDE 1,4 AND 1,6 LINKING",
        "L-SACCHARIDE",
        "L-SACCHARIDE 1,4 AND 1,4 LINKING",
        "L-SACCHARIDE 1,4 AND 1,6 LINKING",
        "NON-POLYMER",
        "OTHER",
        "PEPTIDE-LIKE",
        "SACCHARIDE"])

    asize = mmtf_data.num_atoms

    if len(mmtf_data.chains_per_model) > 1:
        # get number of atoms in first model
        asize = 0
        groupIndex = 0
        modelChainCount = mmtf_data.chains_per_model[0]
        for chainGroupCount in mmtf_data.groups_per_chain[:modelChainCount]:
            #traverse groups
            for _ in range(chainGroupCount):
                group = mmtf_data.group_list[mmtf_data.group_type_list[groupIndex]]        
                asize += len(group['atomNameList'])
                groupIndex += 1

    fields = OrderedDict()

    x = mmtf_data.x_coord_list
    y = mmtf_data.y_coord_list
    z = mmtf_data.z_coord_list

    if len(x) != mmtf_data.num_models*asize:
        LOGGER.warn('Multi-model MMTF files with different molecules not supported.  Keeping only first model')
        coords = np.array([x, y, z]).T[:asize].reshape(1, asize, 3)
    else:     
        coords = np.array([x, y, z]).T.reshape(mmtf_data.num_models, asize, 3)

    # Initialize arrays for atom properties
    atom_names = np.empty(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.empty(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.empty(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.empty(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    segnames = np.empty(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
    segnames[:] = ''
    elements = np.empty(asize, dtype=ATOMIC_FIELDS['element'].dtype)
    icodes = np.empty(asize, dtype=ATOMIC_FIELDS['icode'].dtype)

    hetero = np.zeros(asize, dtype=bool)
    termini = np.zeros(asize, dtype=bool)

    altlocs = np.array(mmtf_data.alt_loc_list[:asize], dtype=ATOMIC_FIELDS['altloc'].dtype)
    serials = np.array(mmtf_data.atom_id_list[:asize], dtype=ATOMIC_FIELDS['serial'].dtype)
    bfactors = np.array(mmtf_data.b_factor_list[:asize], dtype=ATOMIC_FIELDS['beta'].dtype)
    occupancies = np.array(mmtf_data.occupancy_list[:asize], dtype=ATOMIC_FIELDS['occupancy'].dtype)

    modelIndex = 0
    chainIndex = 0
    groupIndex = 0
    atomIndex = 0

    #traverse models
    for modelChainCount in mmtf_data.chains_per_model:
        chain_name_list = mmtf_data.chain_name_list

        #traverse chains
        for k in range(modelChainCount):
            chainGroupCount = mmtf_data.groups_per_chain[chainIndex]

            #traverse groups
            for _ in range(chainGroupCount):
                group = mmtf_data.group_list[mmtf_data.group_type_list[groupIndex]]
                groupAtomCount = len(group.get('atomNameList'))

                #traverse atoms
                for i in range(groupAtomCount):
                    atom_names[atomIndex] = group.get('atomNameList')[i]
                    elements[atomIndex] = group.get('elementList')[i]
                    resnames[atomIndex] = group.get('groupName')
                    resnums[atomIndex] = mmtf_data.group_id_list[groupIndex]
                    icodes[atomIndex] = mmtf_data.ins_code_list[groupIndex]
                    chainids[atomIndex] = chain_name_list[k]
                    if group.get('chemCompType') in mmtfHETATMtypes:
                        hetero[atomIndex] = True

                    atomIndex += 1

                groupIndex += 1
            chainIndex += 1
        modelIndex += 1
        break

    #detect termini based on chain changes
    termini[:-1] = chainids[1:] != chainids[:-1]
    termini[-1] = True  #the last atom is always a terminus

    mask = np.full(asize, True, dtype=bool)
    if altloc_sel != 'all':
        #mask out any unwanted alternative locations
        default_altloc = altloc_sel if altloc_sel != None else 'A'
        mask = (altlocs == '') | (altlocs == default_altloc)
        if np.all(mask == False) and altloc_sel == None and len(altlocs):
            #nothing selected, use first altloc; 6uwi
            mask = altlocs == altlocs[0]            
        
    atomgroup.setCoords(coords[:,mask])
    atomgroup.setNames(atom_names[mask])
    atomgroup.setResnums(resnums[mask])
    atomgroup.setResnames(resnames[mask])
    atomgroup.setChids(chainids[mask])
    atomgroup.setElements(elements[mask])
    atomgroup.setMasses(getMasses(elements[mask]))
    atomgroup.setBetas(bfactors[mask])
    atomgroup.setAltlocs(altlocs[mask])
    atomgroup.setOccupancies(occupancies[mask])
    atomgroup.setFlags('hetatm', hetero[mask])
    atomgroup.setFlags('pdbter', termini[mask])
    atomgroup.setFlags('selpdbter', termini[mask])
    atomgroup.setSerials(serials[mask])
    atomgroup.setIcodes(icodes[mask])
    atomgroup.setSegnames(segnames[mask])
    atomgroup.setTitle(mmtf_data.structure_id)
    
    if get_bonds and hasattr(mmtf_data,'bond_atom_list'):
        #have to remap any masked out atoms
        remaining = np.arange(asize)[mask]
        remap = np.full(asize,-1)
        remap[remaining] = np.arange(len(remaining))
        allbonds = np.array(mmtf_data.bond_atom_list).reshape(-1,2)
        nonpeptide = []
        #irgnore bonds between C and N in adjacent residues
        for a,b in allbonds:
            if a < asize and b < asize and mask[a] and mask[b] :
                if atom_names[a] != 'N' or atom_names[b] != 'C' or resnums[a]-resnums[b] != 1:
                    nonpeptide.append((remap[a],remap[b]))
        atomgroup.setBonds(nonpeptide)

    return atomgroup

def writeMMTF(filename, atoms, csets=None, autoext=True, **kwargs):
    """Write *atoms* in MMTF format to a file with name *filename* and return
    *filename*.  If *filename* ends with :file:`.gz`, a compressed file will
    be written.
    
    :arg atoms: an object with atom and coordinate data
    :type atoms: :class:`.Atomic`

    :arg csets: coordinate set indices, default is all coordinate sets

    :arg autoext: when not present, append extension :file:`.mmtf` to *filename*

    :keyword header: header to write too
    :type header: dict
    """
    try:
        from Bio.PDB.mmtf import MMTFIO
    except ImportError:
        raise ImportError('Biopython MMTFIO could not be imported. '
            'Reinstall ProDy or install Biopython and mmtf-python'
            'to solve the problem.')

    header = kwargs.get('header', None)

    if autoext and not filename.lower().endswith('.mmtf'):
        filename += '.mmtf'

    structure = atoms.toBioPythonStructure(header=header, csets=csets)
    io=MMTFIO()
    io.set_structure(structure)
    io.save(filename)
    return filename
