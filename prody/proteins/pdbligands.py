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

"""This module defines functions for fetching PDB ligand data."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path

import numpy as np

from prody import LOGGER, SETTINGS, getPackagePath
from prody.atomic import AtomGroup, ATOMIC_FIELDS
from prody.utilities import openFile, makePath

__all__ = ['fetchPDBLigand']


def fetchPDBLigand(cci, filename=None):
    """Fetch PDB ligand data from `Ligand Expo <http://ligand-expo.rcsb.org/>`_
    for chemical component *cci*.  *cci* may be 3-letter chemical component 
    identifier or a valid XML filename.  If *filename* is given, XML file 
    will be saved with that name.
    
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
    instead of ``data['formula_weight']`` to avoid exceptions when this
    data field is not found in the XML file.  URL of the XML file is returned 
    in the dictionary with key ``url``.
    
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
      # 2.27"""

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
        xml = None
        cci = cci.upper()
        if SETTINGS.get('ligand_xml_save'):
            folder = os.path.join(getPackagePath(), 'pdbligands')
            if not os.path.isdir(folder):
                makePath(folder)
            xmlgz = os.path.join(folder, cci + '.xml.gz')
            if os.path.isfile(xmlgz):
                LOGGER.info('Parsing XML file for {0:s} from local folder.'
                             .format(cci))
                with openFile(xmlgz) as inp:
                    xml = inp.read()
        url = ('http://ligand-expo.rcsb.org/reports/{0[0]:s}/{0:s}/{0:s}'
               '.xml'.format(cci.upper()))
        if not xml:
            #'http://www.pdb.org/pdb/files/ligand/{0:s}.xml'
            import urllib2
            try:
                inp = urllib2.urlopen(url)
            except urllib2.HTTPError:
                raise IOError('XML file for ligand {0:s} is not found online'
                              .format(cci))
            else:
                xml = inp.read()
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
    
    atomnames = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['name'].dtype) 
    elements = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['element'].dtype)
    resnames = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['resname'].dtype)
    charges = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['charge'].dtype)
    
    resnums = np.ones(n_atoms, dtype=ATOMIC_FIELDS['charge'].dtype)
    
    alternate_atomnames = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['name'].dtype)
    leaving_atom_flags = np.zeros(n_atoms, np.bool)
    aromatic_flags = np.zeros(n_atoms, np.bool)
    stereo_configs = np.zeros(n_atoms, np.bool)
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
        ordinals[i] = int(data.get('pdbx_ordinal',0))
        
        model_coords[i, 0] = float(data.get('model_Cartn_x', 0))
        model_coords[i, 1] = float(data.get('model_Cartn_y', 0))
        model_coords[i, 2] = float(data.get('model_Cartn_z', 0))
        ideal_coords[i, 0] = float(data.get('pdbx_model_Cartn_x_ideal', 0))
        ideal_coords[i, 1] = float(data.get('pdbx_model_Cartn_y_ideal', 0))
        ideal_coords[i, 2] = float(data.get('pdbx_model_Cartn_z_ideal', 0))

    pdbid = dict_.get('model_coordinates_db_code')
    if pdbid:
        model = AtomGroup(cci + ' model ({0:s})'.format(pdbid))
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
    for bond in list(root.find(ns + 'chem_comp_bondCategory')):
        name_1 = bond.get('atom_id_1')
        name_2 = bond.get('atom_id_2')
        bonds.append((name2index[name_1], name2index[name_2]))
    bonds = np.array(bonds, int)
    model.setBonds(bonds)
    ideal.setBonds(bonds)
    return dict_      
