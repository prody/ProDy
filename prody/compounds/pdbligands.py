# -*- coding: utf-8 -*-
"""This module defines functions for fetching PDB ligand data."""

from os.path import isdir, isfile, join, split, splitext

import numpy as np

from prody import LOGGER, SETTINGS, getPackagePath, PY3K
from prody.atomic import AtomGroup, ATOMIC_FIELDS
from prody.utilities import openFile, makePath, openURL

__all__ = ['PDBLigandRecord', 'fetchPDBLigand', 'parsePDBLigand']


class PDBLigandRecord(object):
    """Class for handling the output of fetchPDBLigand"""
    def __init__(self, data):
        self._rawdata = data

    def getCanonicalSMILES(self):
        return self._rawdata['CACTVS_SMILES_CANONICAL']


def fetchPDBLigand(cci, filename=None):
    """Fetch PDB ligand data from PDB_ for chemical component *cci*.
    *cci* may be 3-letter chemical component identifier or a valid XML
    filename.  If *filename* is given, XML file will be saved with that name.

    If you query ligand data frequently, you may configure ProDy to save XML
    files in your computer.  Set ``ligand_xml_save`` option **True**, i.e.
    ``confProDy(ligand_xml_save=True)``.  Compressed XML files will be save
    to ProDy package folder, e.g. :file:`/home/user/.prody/pdbligands`.  Each
    file is around 5Kb when compressed.

    This function is compatible with PDBx/PDBML v 4.0.

    Ligand data is returned in a dictionary.  Ligand coordinate atom data with
    *model* and *ideal* coordinate sets are also stored in this dictionary.
    Note that this dictionary will contain data that is present in the XML
    file and all Ligand Expo XML files do not contain every possible data
    field.  So, it may be better if you use :meth:`dict.get` instead of
    indexing the dictionary, e.g. to retrieve formula weight (or relative
    molar mass) of the chemical component use ``data.get('formula_weight')``
    instead of ``data['formula_weight']`` to avoid exceptions when this data
    field is not found in the XML file.  URL and/or path of the XML file are
    returned in the dictionary with keys ``url`` and ``path``, respectively.

    Following example downloads data for ligand STI (a.k.a. Gleevec and
    Imatinib) and calculates RMSD between model (X-ray structure 1IEP) and
    ideal (energy minimized) coordinate sets:

    .. ipython:: python

       from prody import *
       ligand_data = fetchPDBLigand('STI')
       ligand_data['model_coordinates_db_code']
       ligand_model = ligand_data['model']
       ligand_ideal = ligand_data['ideal']
       transformation = superpose(ligand_ideal.noh, ligand_model.noh)
       calcRMSD(ligand_ideal.noh, ligand_model.noh)"""

    if not isinstance(cci, str):
        raise TypeError('cci must be a string')

    if isfile(cci):
        inp = openFile(cci)
        xml = inp.read()
        inp.close()
        url = None
        path = cci
        cci = splitext(splitext(split(cci)[1])[0])[0].upper()
    elif len(cci) > 4 or not cci.isalnum():
        raise ValueError('cci must be 3-letters long and alphanumeric or '
                         'a valid filename')
    else:
        xml = None
        cci = cci.upper()
        if SETTINGS.get('ligand_xml_save'):
            folder = join(getPackagePath(), 'pdbligands')
            if not isdir(folder):
                makePath(folder)
            xmlgz = path = join(folder, cci + '.xml.gz')
            if isfile(xmlgz):
                with openFile(xmlgz) as inp:
                    xml = inp.read()
        else:
            folder = None
            path = None

        url = ('http://files.rcsb.org/ligands/download/{0}'
               '.xml'.format(cci.upper()))
        if not xml:
            try:
                inp = openURL(url)
            except IOError:
                raise IOError('XML file for ligand {0} is not found online'
                              .format(cci))
            else:
                xml = inp.read()
                if PY3K:
                    xml = xml.decode()
                inp.close()

            if filename:
                out = openFile(filename, mode='w', folder=folder)
                out.write(xml)
                out.close()
            if SETTINGS.get('ligand_xml_save'):
                with openFile(xmlgz, 'w') as out:
                    out.write(xml)

    import xml.etree.cElementTree as ET

    root = ET.XML(xml)
    if (root.get('{http://www.w3.org/2001/XMLSchema-instance}'
                 'schemaLocation') !=
            'http://pdbml.pdb.org/schema/pdbx-v40.xsd pdbx-v40.xsd'):
        LOGGER.warn('XML is not in PDBx/PDBML v 4.0 format, resulting '
                    'dictionary may not contain all data fields')
    ns = root.tag[:root.tag.rfind('}')+1]
    len_ns = len(ns)
    dict_ = {'url': url, 'path': path}

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
                       for audit in
                       list(root.find(ns + 'pdbx_chem_comp_auditCategory'))]

    atoms = list(root.find(ns + 'chem_comp_atomCategory'))
    n_atoms = len(atoms)
    ideal_coords = np.zeros((n_atoms, 3))
    model_coords = np.zeros((n_atoms, 3))

    atomnames = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['name'].dtype)
    elements = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['element'].dtype)
    resnames = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['resname'].dtype)
    charges = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['charge'].dtype)

    resnums = np.ones(n_atoms, dtype=ATOMIC_FIELDS['charge'].dtype)

    alternate_atomnames = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['name'].dtype)
    leaving_atom_flags = np.zeros(n_atoms, bool)
    aromatic_flags = np.zeros(n_atoms, bool)
    stereo_configs = np.zeros(n_atoms, bool)
    ordinals = np.zeros(n_atoms, int)

    name2index = {}

    for i, atom in enumerate(atoms):
        data = dict([(child.tag[len_ns:], child.text) for child in list(atom)])

        name = data.get('pdbx_component_atom_id', 'X')
        name2index[name] = i
        atomnames[i] = name
        elements[i] = data.get('type_symbol', 'X')
        resnames[i] = data.get('pdbx_component_comp_id', 'UNK')
        charges[i] = float(data.get('charge', 0))

        alternate_atomnames[i] = data.get('alt_atom_id', 'X')
        leaving_atom_flags[i] = data.get('pdbx_leaving_atom_flag') == 'Y'
        aromatic_flags[i] = data.get('pdbx_atomatic_flag') == 'Y'
        stereo_configs[i] = data.get('pdbx_stereo_config') == 'Y'
        ordinals[i] = int(data.get('pdbx_ordinal', 0))

        model_coords[i, 0] = float(data.get('model_Cartn_x', 0))
        model_coords[i, 1] = float(data.get('model_Cartn_y', 0))
        model_coords[i, 2] = float(data.get('model_Cartn_z', 0))
        ideal_coords[i, 0] = float(data.get('pdbx_model_Cartn_x_ideal', 0))
        ideal_coords[i, 1] = float(data.get('pdbx_model_Cartn_y_ideal', 0))
        ideal_coords[i, 2] = float(data.get('pdbx_model_Cartn_z_ideal', 0))

    pdbid = dict_.get('model_coordinates_db_code')
    if pdbid:
        model = AtomGroup(cci + ' model ({0})'.format(pdbid))
    else:
        model = AtomGroup(cci + ' model')
    model.setCoords(model_coords)
    model.setNames(atomnames)
    model.setResnames(resnames)
    model.setResnums(resnums)
    model.setElements(elements)
    model.setCharges(charges)
    model.setFlags('leaving_atom_flags', leaving_atom_flags)
    model.setFlags('aromatic_flags', aromatic_flags)
    model.setFlags('stereo_configs', stereo_configs)
    model.setData('ordinals', ordinals)
    model.setData('alternate_atomnames', alternate_atomnames)
    dict_['model'] = model
    ideal = model.copy()
    ideal.setTitle(cci + ' ideal')
    ideal.setCoords(ideal_coords)
    dict_['ideal'] = ideal

    bonds = []
    warned = set()
    for bond in list(root.find(ns + 'chem_comp_bondCategory') or bonds):
        name_1 = bond.get('atom_id_1')
        name_2 = bond.get('atom_id_2')
        try:
            bonds.append((name2index[name_1], name2index[name_2]))
        except KeyError:
            if name_1 not in warned and name_1 not in name2index:
                warned.add(name_1)
                LOGGER.warn('{0} specified {1} in bond category is not '
                            'a valid atom name.'.format(repr(name_1), cci))
            if name_2 not in warned and name_2 not in name2index:
                warned.add(name_2)
                LOGGER.warn('{0} specified {1} in bond category is not '
                            'a valid atom name.'.format(repr(name_2), cci))
    if bonds:
        bonds = np.array(bonds, int)
        model.setBonds(bonds)
        ideal.setBonds(bonds)
    return dict_


def parsePDBLigand(cci, filename=None):
    """See :func:`.fetchPDBLigand`"""
    lig_dict = fetchPDBLigand(cci, filename)
    return PDBLigandRecord(lig_dict)

