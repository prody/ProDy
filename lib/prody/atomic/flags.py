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

"""This module defines atom flags that are used in :ref:`selections`.
You can read this page in interactive sessions using ``help(flags)``.

.. _flags:

Atom flags 
===============================================================================

Flag labels can be used in atom selections:
    
>>> from prody import * 
>>> p = parsePDB('1ubi')
>>> p.select('protein')
<Selection: 'protein' from 1ubi (602 atoms)>

Flag labels can be combined with dot operator as follows to make selections:
    
>>> p.protein
<Selection: 'protein' from 1ubi (602 atoms)>
>>> p.protein.acidic # selects acidic residues
<Selection: '(acidic) and (protein)' from 1ubi (94 atoms)>

Flag labels can be prefixed with ``'is'`` to check whether all atoms in 
an :class:`.Atomic` instance are flagged the same way: 
    
>>> p.protein.ishetero
False
>>> p.water.ishetero
True

Flag labels can also be used to make quick atom counts:
    
>>> p.numAtoms()
683
>>> p.numAtoms('protein')
602
>>> p.numAtoms('water')
81

"""

from re import compile as recompile
from time import time
from collections import defaultdict

from numpy import array, ones, zeros

from prody import SETTINGS
from prody.utilities import joinLinks, joinTerms, tabulate, wrapText

__all__ = ['flagDefinition', 'listNonstdAAProps', 'getNonstdProperties', 
           'addNonstdAminoacid', 'delNonstdAminoacid',]


TIMESTAMP_KEY = 'flags_timestamp'
DEFINITIONS_KEY = 'flags_definitions'
NONSTANDARD_KEY = 'flags_nonstandard'

ALIASES = {}
PLANTERS = {}
EDITORS = {}
FIELDS = defaultdict(set) # flags that fill be nulled when a field changes
FIELDSDEFAULT = ['name', 'resname', 'resnum']
TIMESTAMP = 0

PDBLIGSUM = ('.. _{0:s}: '
             'http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId={0:s}\n')

STANDARDAA = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
              'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
              'TYR', 'VAL']

NONSTANDARD = {
    'ASX': set(['acyclic', 'surface', 'polar', 'medium',]),
    'GLX': set(['acyclic', 'surface', 'large', 'polar']),
    'CSO': set(['acyclic', 'neutral', 'surface', 'medium', 'polar']),
    'HIP': set(['cyclic', 'basic', 'surface', 'large', 'polar']),
    'HSD': set(['cyclic', 'basic', 'surface', 'large', 'polar']),
    'HSE': set(['cyclic', 'basic', 'surface', 'large', 'polar']),
    'HSP': set(['cyclic', 'acidic', 'surface', 'large', 'polar']),
    'SEC': set(['acyclic', 'neutral', 'buried', 'polar', 'medium']),
    'SEP': set(['acyclic', 'surface', 'acidic', 'large', 'polar']),
    'TPO': set(['acyclic', 'surface', 'acidic', 'large', 'polar']),
    'PTR': set(['cyclic', 'aromatic', 'surface', 'acidic', 'large', 'polar']),
    'XLE': set(['aliphatic', 'acyclic', 'buried', 'hydrophobic', 'large']),
    'XAA': set(),
}

DEFAULTS = {

    'stdaa': set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                  'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                  'THR', 'TRP', 'TYR', 'VAL']),
              
    'nonstdaa': set(NONSTANDARD),

    # GLY is excluded from this list
    'aliphatic': set(['ALA', 'ILE', 'LEU', 'VAL', 'PRO', 'GLY']),
    'aromatic': set(['HIS', 'PHE', 'TRP', 'TYR']),
    
    'acidic': set(['ASP', 'GLU']),
    'basic': set(['LYS', 'ARG', 'HIS']),
    'neutral': set(['ALA', 'ASN', 'CYS', 'GLN', 'GLY', 'ILE', 'LEU', 'MET', 
                    'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']),
    
    'buried': set(['ALA', 'CYS', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'VAL']),
    'surface': set(['ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'GLY', 'HIS', 'LYS', 
                    'PRO', 'SER', 'THR', 'TYR']),
    
    'cyclic': set(['HIS', 'PHE', 'PRO', 'TRP', 'TYR']),
    'acyclic': set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                    'ILE', 'LEU', 'LYS', 'MET', 'SER', 'THR', 'VAL']),
    
    'hydrophobic': set(['ALA', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 
                        'VAL']),
    'polar': set(['ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
                  'LYS', 'SER', 'THR', 'TYR']),
                    
    'small': set(['ALA', 'GLY', 'SER']),
    'medium': set(['ASN', 'ASP', 'CYS', 'PRO', 'THR', 'VAL']),
    'large': set(['ARG', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 
                  'PHE', 'TRP', 'TYR']),


    'hydrogen': '[0-9]?H.*', 
    'carbon': 'C.*', 
    'nitrogen': 'N.*',
    'oxygen': 'O.*', 
    'sulfur': 'S.*',
    
    'nucleobase': set(['GUN', 'ADE', 'CYT', 'THY', 'URA']),
    'nucleotide': set(['DA', 'DC', 'DG', 'DT', 'DU', 'A', 'C', 'G', 'T', 'U']),
    'nucleoside': set(['AMP', 'ADP', 'ATP', 'CDP', 'CTP', 'GMP', 'GDP', 'GTP', 
                       'TMP', 'TTP', 'UMP', 'UDP', 'UTP']),
                   
    'at': set(['ADE', 'A', 'THY', 'T']),
    'cg': set(['CYT', 'C', 'GUN', 'G']),
    'purine': set(['ADE', 'A', 'GUN', 'G']),
    'pyrimidine': set(['CYT', 'C', 'THY', 'T', 'URA', 'U']),
    
    'water': set(['HOH', 'DOD', 'WAT', 'TIP3', 'H2O', 'OH2', 'TIP', 'TIP2', 
                  'TIP4']),
    
    'ion': set(['AL', 'BA', 'CA', 'CD', 'CL', 'CO', 'CS', 'CU', 'CU1', 'CUA', 
                'HG', 'IN', 'IOD', 'K', 'MG', 'MN3', 'NA', 'PB', 'PT', 'RB', 
                'TB', 'TL', 'WO4', 'YB', 'ZN']),
    'ion_other': set(['CAL', 'CES', 'CLA', 'POT', 'SOD', 'ZN2']),
    
    
    'lipid': set(['GPE', 'LPP', 'OLA', 'SDS', 'STE']),
    'lipid_other': set(['DLPE', 'DMPC', 'LPPC', 'OLEO', 'PALM', 'PCGL', 'POPC', 
                        'POPE', 'STEA']),
              
    'sugar': set(['GLC', 'GLO', 'BGC']),
    'sugar_other': set(['AGLC']),
    
    'heme': set(['1FH', '2FH', 'DDH', 'DHE', 'HAS', 'HDD', 'HDE', 'HDM', 
                 'HEA', 'HEB', 'HEC', 'HEM', 'HEO', 'HES', 'HEV', 'NTE', 
                 'SRM', 'VER']),
    'heme_other': set(['HEMO', 'HEMR'])

}

DEFAULTS['backbone'] = DEFAULTS['bb'] = set(['CA', 'C', 'O', 'N'])
DEFAULTS['backbonefull'] = DEFAULTS['bbfull'] =  set(['CA', 'C', 'O', 'N', 'H', 
                                                      'H1', 'H2', 'H3', 'OXT'])

CATEGORIES = {
    'aromaticity': set(['aromatic', 'aliphatic']), 
    'cyclic': set(['acyclic', 'cyclic']),
    'charge': set(['acidic', 'basic', 'neutral']),
    'depth': set(['buried', 'surface']),
    'hydrophobicity': set(['hydrophobic', 'polar']),
    'size': set(['small', 'medium', 'large']),
}

CATEGORIZED = {}
for catgroup in CATEGORIES.values():
    for category in catgroup:
        CATEGORIZED[category] = list(DEFAULTS[category])
INDEPENDENT = list(CATEGORIZED) 
CATEGORIZED['charged'] = list(CATEGORIZED['acidic'])
CATEGORIZED['charged'].extend(CATEGORIZED['basic'])

for resi, cats in NONSTANDARD.items():
    for cat in cats:
        CATEGORIZED[cat].append(resi)

def addPlanter(func, *labels, **kwargs):
    
    aliases = kwargs.get('aliases', True)
    editor = kwargs.get('editor', False)
    for label in labels:
        PLANTERS[label] = func
        if aliases:
            ALIASES[label] = labels
        else:
            ALIASES[label] = [label]
        if editor:
            EDITORS[label] = editor
    for field in kwargs.get('fields', FIELDSDEFAULT):
        FIELDS[field].update(labels)


def updateNonstandard(nonstd):
    
    SETTINGS[NONSTANDARD_KEY] = nonstd
    SETTINGS[TIMESTAMP_KEY] = int(time())
    SETTINGS.save()
    updateDefinitions()


def changeDefinitions(**kwargs):
    
    defs = SETTINGS.get(DEFINITIONS_KEY, {})
    defs.update(kwargs)
    SETTINGS[DEFINITIONS_KEY] = defs
    SETTINGS[TIMESTAMP_KEY] = int(time())
    SETTINGS.save()
    updateDefinitions()

    
def resetDefinitions(flag):
    
    if flag == 'all':
        SETTINGS.pop(DEFINITIONS_KEY, None)
        SETTINGS.pop(NONSTANDARD_KEY, None)
        SETTINGS[TIMESTAMP_KEY] = int(time())
        SETTINGS.save()
        updateDefinitions()
    elif flag == 'nonstdaa': 
        SETTINGS.pop(NONSTANDARD_KEY, None)
        SETTINGS[TIMESTAMP_KEY] = int(time())
        SETTINGS.save()
        updateDefinitions()
    else:        
        try:
            SETTINGS.pop(DEFINITIONS_KEY, {}).pop(flag)
        except KeyError:
            pass
        else:
            SETTINGS[TIMESTAMP_KEY] = int(time())
            SETTINGS.save()
            updateDefinitions()


def changeResnames(flag, resnames):
    """Change the list of residue names associated with a flag.  *flag* must 
    be a string, and *resnames* may be a list, tuple, or set of strings.  The 
    existing list of residue names will be overwritten with the given residue
    names."""
    
    if flag not in EDITORS:
        raise ValueError('{0:s} is not an editable flag'.format(repr(flag)))
    
    resnames = [str(rn) for rn in resnames]
    changeDefinitions(**{flag: resnames})



def changeNameRegex(flag, regex):
    """Set regular expression used for flagging elements based on atom names.
    See :ref:`element-flags` for the default list of flags."""
    
    if flag not in EDITORS:
        raise ValueError('{0:s} is not an editable flag'.format(repr(flag)))

    try:
        recompile(regex)
    except Exception as err:
        raise ValueError('{0:s} is not a valid regular expression, {1:s}'
                         .format(repr(regex), str(err)))
    else:
        changeDefinitions(**{flag: regex})

def changeBackbone(flag, names):
    """Set protein :term:`backbone` or :term:`backbonefull` atom names.  
    *names* must be a list of atom names."""
    
    if flag not in EDITORS:
        raise ValueError('{0:s} is not an editable flag'.format(repr(flag)))
    
    names = [str(nm) for nm in names]
    
    if flag.endswith('full'):    
        changeDefinitions(bbfull=names, backbonefull=names)
    else:
        changeDefinitions(bb=names, backbone=names)
    

__doc__ += """

Protein
-------------------------------------------------------------------------------

.. glossary::

   protein
   aminoacid
      indicates the twenty standard amino acids (:term:`stdaa`) and some 
      non-standard amino acids (:term:`nonstdaa`) described below.  Residue 
      must also have an atom named ``'CA'`` in addition to having a qualifying
      residue name.


   stdaa
      {stdaa:s}


   nonstdaa
      indicates one of the following residues:

      ==========  ===================================================
      `ASX`_ (B)  asparagine or aspartic acid
      `GLX`_ (Z)  glutamine or glutamic acid
      `CSO`_ (C)  S-hydroxycysteine
      `HIP`_ (H)  ND1-phosphohistidine
       HSD   (H)  prototropic tautomer of histidine, H on ND1 (CHARMM)
       HSE   (H)  prototropic tautomer of histidine, H on NE2 (CHARMM)       
       HSP   (H)  protonated histidine 
      `SEC`_ (U)  selenocysteine
      `SEP`_ (S)  phosphoserine
      `TPO`_ (T)  phosphothreonine
      `PTR`_ (Y)  O-phosphotyrosine
       XLE   (J)  leucine or isoleucine 
       XAA   (X)  unspecified or unknown 
      ==========  ===================================================

      You can modify the list of non-standard amino acids using 
      :func:`addNonstdAminoacid`, :func:`delNonstdAminoacid`, and 
      :func:`listNonstdAAProps`.


   calpha
   ca
      Cα atoms of :term:`protein` residues, same as selection 
      ``'name CA and protein'``


   backbone
   bb
      non-hydrogen backbone atoms of :term:`protein` residues, same as
      selection ``'name CA C O N and protein'``


   backbonefull
   bbfull
      backbone atoms of :term:`protein` residues, same as selection  
      ``'name CA C O N H H1 H2 H3 OXT and protein'``


   sidechain
   sc
      side-chain atoms of :term:`protein` residues, same as selection 
      ``'protein and not backbonefull'``


""".format(

stdaa = wrapText('indicates the standard amino acid residues: ' + 
                 joinLinks(DEFAULTS['stdaa'], last=', and ', sort=True), 
                 subsequent_indent=' '*6),
)



cats = list(CATEGORIZED)
cats.sort()

for cat in cats:
    res = CATEGORIZED[cat]
    res.sort() 
    __doc__ += """
   {cat:s}
       {res:s}

""".format(cat=cat, 
res=wrapText('residues ' + ', '.join(res), subsequent_indent=' '*6))

for resi in DEFAULTS['stdaa']:
    __doc__ += PDBLIGSUM.format(resi)
for resi in DEFAULTS['nonstdaa']:
    __doc__ += PDBLIGSUM.format(resi)


__doc__ += """

Nucleic
-------------------------------------------------------------------------------

.. glossary::

   nucleic
      indicates :term:`nucleobase`, :term:`nucleotide`, and some
      :term:`nucleoside` derivatives that are described below, so it is same
      as ``'nucleobase or nucleotide or nucleoside'``.

   nucleobase
      indicates `ADE`_ (adenine), `GUN`_ (guanine), `CYT`_ 
      (cytosine), `THY`_ (thymine), and `URA`_ (uracil).

   nucleotide
      indicates residues with the following names:

      =======  ==================================
      `DA`_    2'-deoxyadenosine-5'-monophosphate
      `DC`_    2'-deoxycytidine-5'-monophosphate
      `DG`_    2'-deoxyguanosine-5'-monophosphate
      `DT`_    2'-deoxythymidine-5'-monophosphate
      `DU`_    2'-deoxyuridine-5'-monophosphate
      `A`_     adenosine-5'-monophosphate
      `C`_     cytidine-5'-monophosphate
      `G`_     guanosine-5'-monophosphate
      `T`_     2'-deoxythymidine-5'-monophosphate
      `U`_     uridine-5'-monophosphate
      =======  ==================================

   nucleoside
      indicates following nucleoside derivatives that are recognized by *PDB*:

      =======  ================================
      `AMP`_   adenosine monophosphate             
      `ADP`_   adenosine-5'-diphosphate
      `ATP`_   adenosine-5'-triphosphate
      `CDP`_   cytidine-5'-diphosphate
      `CTP`_   cytidine-5'-triphosphate
      `GMP`_   guanosine
      `GDP`_   guanosine-5'-diphosphate
      `GTP`_   guanosine-5'-triphosphate
      `TMP`_   thymidine-5'-phosphate
      `TTP`_   thymidine-5'-triphosphate
      `UMP`_   2'-deoxyuridine 5'-monophosphate
      `UDP`_   uridine 5'-diphosphate
      `UTP`_   uridine 5'-triphosphate
      =======  ================================


   at
      same as selection ``'resname ADE A THY T'``

   cg
      same as selection ``'resname CYT C GUN G'``

   purine
      same as selection ``'resname ADE A GUN G'``

   pyrimidine
      same as selection ``'resname CYT C THY T URA U'``

"""

for resi in DEFAULTS['nucleobase']:
    __doc__ += PDBLIGSUM.format(resi)

for resi in DEFAULTS['nucleotide']:
    __doc__ += PDBLIGSUM.format(resi)

for resi in DEFAULTS['nucleoside']:
    __doc__ += PDBLIGSUM.format(resi)

__doc__ += """

Heteros
-------------------------------------------------------------------------------

.. glossary::
   
   hetero
      indicates anything other than a :term:`protein` or a 
      :term:`nucleic` residue, i.e. ``'not (protein or nucleic)'``.


   hetatm
      is available when atomic data is parsed from a PDB or similar 
      format file and indicates atoms that are marked ``'HETATM'`` in the file. 


   water
      indices `HOH`_ and `DOD`_ recognized by *PDB* and also WAT, TIP3, H2O, 
      OH2, TIP, TIP2, and TIP4 recognized by molecular dynamics (MD) force 
      fields.  
    
      .. _HOH: http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId=HOH
    
      .. _DOD: http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId=DOD
    
      Previously used water types HH0, OHH, and SOL conflict with other 
      compounds in the *PDB*, so are removed from the definition of this flag.


   ion
      indicates the following ions most of which are recognized by the *PDB* 
      and others by MD force fields.

      =======  ==================  =====  ========  ==========
      \                            *PDB*  *Source*  *Conflict*
      `AL`_    aluminum            Yes
      `BA`_    barium              Yes
      `CA`_    calcium             Yes
      `CD`_    cadmium             Yes
      `CL`_    chloride            Yes
      `CO`_    cobalt (ii)         Yes
      `CS`_    cesium              Yes
      `CU`_    copper (ii)         Yes
      `CU1`_   copper (i)          Yes
      `CUA`_   dinuclear copper    Yes
      `HG`_    mercury (ii)        Yes
      `IN`_    indium (iii)        Yes
      `IOD`_   iodide              Yes
      `K`_     potassium           Yes
      `MG`_    magnesium           Yes
      `MN3`_   manganese (iii)     Yes
      `NA`_    sodium              Yes
      `PB`_    lead (ii)           Yes
      `PT`_    platinum (ii)       Yes
      `RB`_    rubidium            Yes
      `TB`_    terbium (iii)       Yes
      `TL`_    thallium (i)        Yes
      `WO4`_   thungstate (vi)     Yes
      `YB`_    ytterbium (iii)     Yes
      `ZN`_    zinc                Yes
      CAL      calcium             No     CHARMM    Yes
      CES      cesium              No     CHARMM    Yes
      CLA      chloride            No     CHARMM    Yes
      POT      potassium           No     CHARMM    Yes
      SOD      sodium              No     CHARMM    Yes
      ZN2      zinc                No     CHARMM    No
      =======  ==================  =====  ========  ==========

      Ion identifiers that are obsoleted by *PDB* (MO3, MO4, MO5, MO6, NAW, 
      OC7, and ZN1) are removed from this definition.


   lipid
      {lipid:s}


   sugar
      {sugar:s}


   heme
      {heme:s}


""".format(
 
lipid=wrapText('indicates ' + 
               joinLinks(DEFAULTS['lipid'], last=', and ', sort=True) + 
               ' from *PDB*, and also ' + ', '.join(DEFAULTS['lipid_other']) +
               ' from *CHARMM* force field.', subsequent_indent=' '*6),
           
sugar=wrapText('indicates ' + 
               joinLinks(DEFAULTS['sugar'], last=', and ', sort=True) + 
               ' from *PDB*, and also AGLC from *CHARMM*.', 
               subsequent_indent=' '*6),

heme=wrapText('indicates ' + 
              joinLinks(DEFAULTS['heme'], last=', and ', sort=True) +
              ' from *PDB*, and also HEMO and HEMR from CHARMM.',
              subsequent_indent=' '*6)
)


for resi in DEFAULTS['ion']:
    __doc__ += PDBLIGSUM.format(resi)

for resi in DEFAULTS['lipid']:
    __doc__ += PDBLIGSUM.format(resi)

for resi in DEFAULTS['sugar']:
    __doc__ += PDBLIGSUM.format(resi)

for resi in DEFAULTS['heme']:
    __doc__ += PDBLIGSUM.format(resi)


__doc__ += """

Elements
-------------------------------------------------------------------------------

Following elements found in proteins are recognized by applying regular 
expressions to atom names:

.. glossary::
   
   carbon
      carbon atoms, same as ``'name "C.*" and not ion'``

   nitrogen
      nitrogen atoms, same as ``'name "N.*" and not ion'``

   oxygen
      oxygen atoms, same as ``'name "O.*" and not ion'``

   sulfur
      sulfur atoms, same as ``'name "S.*" and not ion'``

   hydrogen
      hydrogen atoms, same as ``'name "[1-9]?H.*" and not ion'``

   noh
      non hydrogen atoms, same as ``'not hydrogen``


``'not ion'`` is appended to above definitions to avoid conflicts with 
:term:`ion` atoms.

Structure
-------------------------------------------------------------------------------

Following secondary structure flags are defined but before they can be used,
secondary structure assignments must be made. 

.. glossary::

   extended
      extended conformation, same as ``'secondary E'``

   helix
      α-helix conformation, same as ``'secondary H'``

   helix310
      3_10-helix conformation, same as ``'secondary G'``

   helixpi
      π-helix conformation, same as ``'secondary I'``

   turn
      hydrogen bonded turn conformation, same as ``'secondary T'``

   bridge
      isolated beta-bridge conformation, same as ``'secondary B'``

   bend
      bend conformation, same as ``'secondary S'``

   coil
      not in one of above conformations, same as ``'secondary C'``

Others
-------------------------------------------------------------------------------

.. glossary::
    
   dummy
      indicates dummy atoms in an :class:`.AtomMap`
   
   mapped
      indicates mapped atoms in an :class:`.AtomMap`

Functions
===============================================================================

Following functions can be used to customize flag definitions:
    
    * :func:`.flagDefinition`
    * :func:`.addNonstdAminoacid`
    * :func:`.delNonstdAminoacid`
    * :func:`.listNonstdAAProps`

"""

# join split lists
for key in ['ion', 'lipid', 'sugar', 'heme']:
    DEFAULTS[key].update(DEFAULTS.pop(key + '_other'))

DEFINITIONS = None
AMINOACIDS = None
BACKBONE = None

def updateDefinitions():
    """Update definitions and set some global variables.  This function must be
    called at the end of the module."""

    global DEFINITIONS, AMINOACIDS, BACKBONE, TIMESTAMP
    DEFINITIONS = {}
    user = SETTINGS.get('flag_definitions', {})
    
    # nucleics
    nucleic = set()
    for key in ['nucleobase', 'nucleoside', 'nucleotide']:
        aset = set(user.get(key, DEFAULTS[key]))
        nucleic.update(aset)
        DEFINITIONS[key] = aset
    DEFINITIONS['nucleic'] = nucleic
    
    # heteros
    for key in ['water', 'lipid', 'ion', 'sugar', 'heme', 
                 'at', 'cg', 'purine', 'pyrimidine',]:
        DEFINITIONS[key] = set(user.get(key, DEFAULTS[key]))
        
    DEFINITIONS['backbone'] = DEFINITIONS['bb'] = set(user.get(key, 
                                                           DEFAULTS['bb']))
    DEFINITIONS['backbonefull'] = DEFINITIONS['bbfull'] = set(user.get(key, 
                                                           DEFAULTS['bbfull']))

    # element regex
    for key in ['hydrogen', 'carbon', 'nitrogen', 'oxygen', 'sulfur']:
        DEFINITIONS[key] = recompile(user.get(key, DEFAULTS[key]))

    try:
        nonstd = SETTINGS[NONSTANDARD_KEY]
        
    except KeyError:
        nonstd = NONSTANDARD
        DEFINITIONS.update(CATEGORIZED)
    else:

        for cat in CATEGORIES:
            for key in CATEGORIES[cat]:
                DEFINITIONS[key] = set(DEFAULTS[key])

        DEFINITIONS['charged'] = set(DEFINITIONS['acidic'])
        DEFINITIONS['charged'].update(DEFINITIONS['basic'])

        for resi, props in nonstd.items():
            for prop in props: 
                DEFINITIONS[prop].add(resi)

    DEFINITIONS['stdaa'] = DEFAULTS['stdaa']
    DEFINITIONS['nonstdaa'] = set(nonstd)
    AMINOACIDS = set(DEFINITIONS['stdaa'])
    AMINOACIDS.update(DEFINITIONS['nonstdaa'])
    DEFINITIONS['protein'] = DEFINITIONS['aminoacid'] = AMINOACIDS
    
    BACKBONE = DEFINITIONS['bb']

    global TIMESTAMP
    TIMESTAMP = SETTINGS.get(TIMESTAMP_KEY, 0)


#==============================================================================
# PLANTERS
#==============================================================================

# flag planters are functions that take an AtomGroup instance
# and set flags for a given label using the atomic data stored in the AtomGroup
# a planter can also set the indices of subset atoms for calculated flag
# planters return the flags after calculation is over

# protein
#==============================================================================

def setCalpha(ag, label):
    
    flags = ag._getNames() == 'CA'
    indices = flags.nonzero()[0]
    if len(indices):
        torf = array([rn in AMINOACIDS for rn in ag._getResnames()[indices]])
        flags[indices] = torf
        ag._setFlags(label, flags)
        ag._setSubset(label, indices[torf])
    else:
        ag._setFlags(label, flags)
        ag._setSubset(label, array([]))
    return flags

addPlanter(setCalpha, 'ca', 'calpha')
    
def setProtein(ag, label):
    
    calpha = ag._getSubset('calpha')
    if len(calpha):
        torf = zeros(ag.numResidues(), bool)
        resindices = ag._getResindices()
        torf[resindices[calpha]] = True
        flags = torf[resindices]
    else:
        flags = zeros(ag.numAtoms(), bool)
    ag._setFlags('protein', flags)
    return flags

addPlanter(setProtein, 'protein', 'aminoacid')

# subsets
#==============================================================================

def setBackbone(ag, label):
    
    protein = ag._getSubset('protein')
    if len(protein):
        bb = DEFINITIONS[label] 
        names = ag._getNames()[protein]
        flags = zeros(ag.numAtoms(), bool)
        flags[protein] = [nm in bb for nm in names]
    else:
        flags = zeros(ag.numAtoms(), bool)
    ag._setFlags(label, flags)
    return flags

addPlanter(setBackbone, 'bb', 'backbone', editor=changeBackbone)
addPlanter(setBackbone, 'bbfull', 'backbonefull', editor=changeBackbone)

def setSidechain(ag, label):
    
    protein = ag._getSubset('protein')
    if len(protein):
        flags = zeros(ag.numAtoms(), bool)
        flags[protein] = True
        flags[ag._getSubset('bbfull')] = False
    else:
        flags = zeros(ag.numAtoms(), bool)
    ag._setFlags(label, flags)
    return flags

addPlanter(setSidechain, 'sc', 'sidechain')


def setCategories(ag, label):
    
    calpha = ag._getSubset('ca')
    if len(calpha):
        resnames = DEFINITIONS[label]
        residx = ag._getResindices()
        torf = zeros(ag.numResidues(), bool)
        torf[residx[calpha]] = [rn in resnames 
                                for rn in ag._getResnames()[calpha]]
        flags = torf[residx]
    else:
        flags = zeros(ag.numAtoms(), bool)
    ag._setFlags(label, flags)
    return flags
    
addPlanter(setCategories, 'stdaa', 'nonstdaa', *list(CATEGORIZED.keys()), 
           aliases=False)


# hetero, nucleic, water, etc.
#==============================================================================

def setResiflag(ag, label):    
    
    resnames = DEFINITIONS[label]
    flags = array([rn in resnames for rn in ag._getResnames()], bool)
    ag._setFlags(label, flags)
    return flags

addPlanter(setResiflag, 'nucleobase', 'nucleoside', 'nucleotide',
           'water', 'ion', 'lipid', 'sugar', 'heme', 'at', 'cg', 'purine', 
           'pyrimidine', aliases=False, editor=changeResnames)
addPlanter(setResiflag, 'nucleic') 


def setHetero(ag, label):
    
    flags = ones(ag.numAtoms(), bool)
    flags[ag._getSubset('protein')] = False
    flags[ag._getSubset('nucleic')] = False
    ag._setFlags('hetero', flags)
    return flags

addPlanter(setHetero, 'hetero',)


# element
#==============================================================================

def setElement(ag, label):
    
    match = DEFINITIONS[label].match 
    flags = array([match(nm) is not None for nm in ag._getNames()])
    flags[ag._getSubset('ion')] = False
    ag._setFlags(label, flags)
    return flags

addPlanter(setElement, 'hydrogen', 'carbon', 'nitrogen', 'oxygen', 'sulfur', 
           aliases=False, editor=changeNameRegex)

def setNoh(ag, label):
    
    flags = ones(ag.numAtoms(), bool)
    flags[ag._getSubset('hydrogen')] = False
    ag._setFlags(label, flags)
    return flags

addPlanter(setNoh, 'noh')


# secondary
#==============================================================================

SECONDARY = {'extended': 'E', 'helix': 'H', 'helix310': 'G', 'helixpi': 'I',
             'turn': 'T', 'bridge': 'B', 'bend': 'S', 'coil': 'C'}

def setSecondary(ag, label):
    
    secstrs = ag._getSecstrs()
    if secstrs is None:
        raise ValueError('no secondary structure assignments')
    flags = secstrs == SECONDARY[label]
    ag._setFlags(label, flags)
    return flags

addPlanter(setSecondary, *list(SECONDARY.keys()), aliases=False, 
           fields=['secondary']) 


# secondary
#==============================================================================

def setDummy(ag, label):
    
    if label == 'dummy':
        flags = zeros(len(ag), bool)
    else:
        flags = ones(len(ag), bool)
    ag._setFlags(label, flags)
    return flags

addPlanter(setDummy, 'dummy', 'mapped', aliases=False, fields=[]) 


#==============================================================================
# User functions
#==============================================================================

def flagDefinition(*arg, **kwarg):
    """Learn, change, or reset :ref:`flags` definitions.
    
    **Learn a definition**
    
    Calling this function with no arguments will return list of flag names
    whose definitions you can learn:
    
    >>> flagDefinition() # doctest: +ELLIPSIS
    ['acidic', 'acyclic', 'aliphatic', ..., 'sulfur', 'surface', 'water']
    
    Passing a flag name will return its definition:
    
    >>> flagDefinition('backbone')
    ['C', 'CA', 'N', 'O']
    >>> flagDefinition('hydrogen')
    '[0-9]?H.*'

    **Change a definition**    

    Calling the function with ``editable=True`` argument will return flag
    names those definitions that can be edited: 
        
    >>> flagDefinition(editable=True) # doctest: +ELLIPSIS
    ['at', 'backbone', 'backbonefull', ..., 'sugar', 'sulfur', 'water']

    Pass an editable flag name with its new definition:
    
    >>> flagDefinition(nitrogen='N.*')
    >>> flagDefinition(backbone=['CA', 'C', 'O', 'N'])
    >>> flagDefinition(nucleobase=['ADE', 'CYT', 'GUN', 'THY', 'URA'])
    
    Note that the type of the new definition must be the same as the type
    of the old definition.  Flags with editable definitions are: 
    {editable}
    
    **Reset definitions**
    
    Pass *reset* keyword as follows to restore all default definitions of 
    editable flags and also non-standard amino acids.
    
    >>> flagDefinition(reset='all')
    
    Or, pass a specific editable flag label to restore its definition:
        
    >>> flagDefinition(reset='nitrogen')"""

    if arg and kwarg:
        raise ValueError('only a single argument or a single keyword '
                           'argument is accepted at a time')
        
    elif arg:
        if len(arg) > 1:
            raise ValueError('specify only one flag label')
        arg = arg[0]
        try:
            definition = DEFINITIONS[arg]
        except KeyError:
            raise ValueError('{0:s} is not a valid flag label'
                               .format(repr(arg)))
        else:
            try:
                return definition.pattern
            except AttributeError:  
                definition = list(definition)
                definition.sort()
                return definition

    elif kwarg:
        if len(kwarg) > 1:
            raise ValueError('specify only one keyword argument')
        key, val = list(kwarg.items())[0]
        if key == 'editable' and val: 
            alist = list(EDITORS)
            alist.sort()
            return alist
        if key == 'reset' and val: 
            resetDefinitions(val)
            return  
        try:
            editor = EDITORS[key]
        except KeyError:
            raise ValueError('{0:s} is not an editable flag or a valid '
                               'keyword'.format(repr(key)))
        else:
            type_ = type(flagDefinition(key))
            if type(val) != type_:
                raise TypeError('expected {0:s} as type of new definition, '
                                  'found {1:s}'.format(repr(type_.__name__), 
                                  repr(type(val).__name__)))
            editor(key, val)
    else:
        alist = list(DEFINITIONS)
        alist.sort()
        return alist

flagDefinition.__doc__ = flagDefinition.__doc__.format(
    editable=wrapText(joinTerms(flagDefinition(editable=True), last=', and '), 
                      subsequent_indent=' '*4))


def addNonstdAminoacid(resname, *properties):
    """Add non-standard amino acid *resname* with *properties* selected from:
     
      * {props}
    
    >>> addNonstdAminoacid('PTR', 'acidic', 'aromatic', 'cyclic', 'large', 
    ... 'polar', 'surface')
    
    Default set of non-standard amino acids can be restored as follows:
    
    >>> flagDefinition(reset='nonstdaa')"""
    
    resname = str(resname)
    if len(resname) > 4:
        LOGGER.warn('Residue name {0:s} is unusually long.'
                    .format(repr(resname)))
    propset = set(properties)
    for cat, val in CATEGORIES.items():
        intersection = val.intersection(propset)
        if intersection:
            if len(intersection) > 1:
                raise ValueError('amino acid properties {0:s} cannot be '
                                   'present together'
                      .format(', '.join([repr(prp) for prp in intersection])))
            for prop in intersection:
                propset.remove(prop)
    if propset:
        raise ValueError('amino acid property {0:s} is not valid'
                           .format(repr(propset.pop())))
        
    nonstd = SETTINGS.get(NONSTANDARD_KEY, NONSTANDARD)
    nonstd[resname] = set(properties)
    updateNonstandard(nonstd)    
  
INDEPENDENT.sort()
addNonstdAminoacid.__doc__ = addNonstdAminoacid.__doc__.format(
    props='\n      * '.join([
        '*{0:s}*: {1:s}'.format(cat, joinTerms(terms, last=', or ', sort=True))
        for cat, terms in CATEGORIES.items()
    ])
)
    
    
def delNonstdAminoacid(resname):
    """Delete non-standard amino acid *resname*.
    
    >>> delNonstdAminoacid('PTR')
    >>> flagDefinition('nonstdaa') # doctest: +ELLIPSIS
    ['ASX', 'CSO', 'GLX', ..., 'TPO', 'XAA', 'XLE']
    
    Default set of non-standard amino acids can be restored as follows:
    
    >>> flagDefinition(reset='nonstdaa')"""
    
    
    nonstd = SETTINGS.get(NONSTANDARD_KEY, NONSTANDARD)
    try:
        nonstd.pop(resname)
    except KeyError:
        raise ValueError('{0:s} is not a non-standard residue name'
                           .format(repr(resname)))
    else:
        updateNonstandard(nonstd)


def getNonstdProperties(resname):
    """Deprecated for removal in v1.4, use :func:`listNonstdAAProps` instead.
    """
    
    from prody import deprecate
    deprecate('getNonstdProperties', 'listNonstdAAProps')
    
    return listNonstdAAProps(resname) 


def listNonstdAAProps(resname):
    """Return properties of non-standard amino acid *resname*.
    
    >>> listNonstdAAProps('PTR')
    ['acidic', 'aromatic', 'cyclic', 'large', 'polar', 'surface']"""
    
    try:
        alist = list(SETTINGS.get(NONSTANDARD_KEY, NONSTANDARD)[resname])
    except KeyError:
        raise ValueError('{0:s} is not a non-standard residue name'
                           .format(repr(resname)))
    else:
        alist.sort()
        return alist
        

updateDefinitions()
