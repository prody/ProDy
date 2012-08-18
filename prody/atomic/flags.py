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
from textwrap import wrap as textwrap
from collections import defaultdict

from numpy import array, ones, zeros

from prody import SETTINGS
from prody.utilities import tabulate, wrap

__all__ = ['flagDefinition', 'addNonstdAA', 'delNonstdAA', 'getNonstdAAs',]

ALIASES = {}
PLANTERS = {}
FIELDS = defaultdict(set)
FIELDSDEFAULT = ['name', 'resname', 'resnum']
TIMESTAMP = None

def timestamp():
    global TIMESTAMP
    TIMESTAMP = time()

PDBLIGSUM = ('.. _{0:s}: '
             'http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId={0:s}\n')

def addPlanter(func, *labels, **kwargs):
    
    aliases = kwargs.get('aliases', True)
    for label in labels:
        PLANTERS[label] = func
        if aliases:
            ALIASES[label] = labels
        else:
            ALIASES[label] = [label]
    for field in kwargs.get('fields', FIELDSDEFAULT):
        FIELDS[field].update(labels)

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
CATEGORIZED['charged'] = list(CATEGORIZED['acidic'])
CATEGORIZED['charged'].extend(CATEGORIZED['basic'])

for resi, cats in NONSTANDARD.iteritems():
    for cat in cats:
        CATEGORIZED[cat].append(resi)
    
__doc__ += """

Protein atoms
===============================================================================

.. glossary::

   protein
   aminoacid
      These flags indicate the twenty standard amino acids (:term:`stdaa`) 
      and some non-standard amino acids (:term:`nonstdaa`) described below.  
      Residue must also have an atom named ``'CA'`` in addition to having
      a qualifying residue name.


   stdaa
      {stdaa:s}


   nonstdaa
      This flag indicates one of the following residues:

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
      :func:`addNonstdAA`, :func:`delNonstdAA`, and :func:`getNonstdAAs`.


   ca
   calpha
      Cα atoms of :term:`protein` residues, same as selection 
      ``'name CA and protein'``


   bb
   backbone        
      non-hydrogen backbone atoms of :term:`protein` residues, same as
      selection ``'name CA C O N and protein'``


   bbfull
   backbonefull
      backbone atoms of :term:`protein` residues, same as selection  
      ``'name CA C O N H H1 H2 H3 OXT and protein'``


   sc
   sidechain
      side-chain atoms of :term:`protein` residues, same as selection 
      ``'not backbone'``


""".format(

stdaa = wrap('This flag indicates the standard amino acid residues: `' + 
             '`_, `'.join(DEFAULTS['stdaa']) + '`_', subsequent_indent=' '*6),
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
res=wrap('residues ' + ', '.join(res), subsequent_indent=' '*6))

for resi in DEFAULTS['stdaa']:
    __doc__ += PDBLIGSUM.format(resi)
for resi in DEFAULTS['nonstdaa']:
    __doc__ += PDBLIGSUM.format(resi)


__doc__ += """

Nucleic atoms
===============================================================================

.. glossary::

   nucleic
      This flag indicates :term:`nucleobase`, :term:`nucleotide`, and some 
      :term:`nucleoside` derivatives that are described below, so it is same
      as ``'nucleobase or nucleotide or nucleoside'``.

   nucleobase
      This flag indicates `ADE`_ (adenine), `GUN`_ (guanine), `CYT`_ 
      (cytosine), `THY`_ (thymine), and `URA`_ (uracil).

   nucleotide
      This flag indicates residues with the following names:

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
      This flag indicates following nucleoside derivatives that are recognized 
      by *PDB*:

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

Hetero atoms
===============================================================================

.. glossary::
   
   hetero
      This flag indicates anything other than a :term:`protein` or a 
      :term:`nucleic` residue, i.e. ``'not (protein or nucleic)'``.


   hetatm
      This flag is available when atomic data is parsed from a PDB or similar 
      format file and indicates atoms that are marked ``'HETATM'`` in the file. 


   water
      This flag definition includes `HOH`_ and `DOD`_ recognized by *PDB* and 
      also WAT, TIP3, H2O, OH2, TIP, TIP2, and TIP4 recognized by some molecular 
      dynamics (MD) force fields.  
    
      .. _HOH: http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId=HOH
    
      .. _DOD: http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId=DOD
    
      Previously used water types HH0, OHH, and SOL conflict with other 
      compounds in the *PDB*, so are removed from the definition of this flag.


   ion
      This flag definition includes the following ions most of which are 
      recognized by the *PDB* and others by some MD force fields.

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
      {lipid:s}


   heme
      {heme:s}


""".format(
 
lipid=wrap('This flag indicates `' + '`_, `'.join(DEFAULTS['lipid']) + 
           '`_ from *PDB* and ' + ', '.join(DEFAULTS['lipid_other']) +
           ' from *CHARMM* force field.', subsequent_indent=' '*6),
           
sugar = wrap('**sugar** flag indicates `' + '`_, `'.join(DEFAULTS['sugar']) + 
             '`_ from *PDB* and AGLC from *CHARMM*.', subsequent_indent=' '*6),

heme = wrap('This flag indicates `' + '`_, `'.join(DEFAULTS['heme']) + 
            '`_ from *PDB*, as well as HEMO and HEMR from CHARMM.',
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
===============================================================================

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


``'not ion'`` is appended to above definitions to avoid conflicts with ion
atoms.

Structure
===============================================================================

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

"""

# join split lists
for key in ['ion', 'lipid', 'sugar', 'heme']:
    DEFAULTS[key].update(DEFAULTS.pop(key + '_other'))

EDITABLE = set(['nucleobase', 'nucleoside', 'nucleotide',
                'at', 'cg', 'purine', 'pyrimidine', 
                'water', 'heme', 'lipid', 'ion', 'sugar',
                'backbone', 'backbonefull', 'bb', 'bbfull',
                'hydrogen', 'carbon', 'nitrogen', 'oxygen', 'sulfur'])

DEFINITIONS = None
AMINOACIDS = None
BACKBONE = None

def setDefinitions():
    

    global DEFINITIONS, AMINOACIDS, BACKBONE
    DEFINITIONS = {}
    defs = SETTINGS.get('flag_definitions', DEFAULTS)
    
    # nucleics
    nucleic = set()
    for key in ['nucleobase', 'nucleoside', 'nucleotide']:
        aset = set(defs[key])
        nucleic.update(aset)
        DEFINITIONS[key] = aset
    DEFINITIONS['nucleic'] = nucleic
    
    # heteros
    for key in ['water', 'lipid', 'ion', 'sugar', 'heme', 
                 'at', 'cg', 'purine', 'pyrimidine',]:
        DEFINITIONS[key] = set(defs[key])
        
    DEFINITIONS['backbone'] = DEFINITIONS['bb'] = set(defs['bb'])
    DEFINITIONS['backbonefull'] = DEFINITIONS['bbfull'] = set(defs['bbfull'])

    # element regex
    for key in ['hydrogen', 'carbon', 'nitrogen', 'oxygen', 'sulfur']:
        DEFINITIONS[key] = recompile(defs[key])
        

    try:
        nonstd = SETTINGS['nonstandard']
        
    except KeyError:
        nonstd = NONSTANDARD
        DEFINITIONS.update(CATEGORIZED)
    else:

        for key in CATEGORIES:
            DEFINITIONS[key] = set(DEFAULTS[key])

        DEFINITIONS['charged'] = set(DEFINITIONS['acidic'])
        DEFINITIONS['charged'].update(DEFINITIONS['basic'])

        for resi, props in nonstd.iteritems():
            for prop in props: 
                DEFINITIONS[prop].add(resi)

    DEFINITIONS['stdaa'] = DEFAULTS['stdaa']
    DEFINITIONS['nonstdaa'] = set(nonstd.keys())
    AMINOACIDS = set()
    AMINOACIDS.update(DEFAULTS['stdaa'])
    AMINOACIDS.update(nonstd.keys())
    DEFINITIONS['protein'] = DEFINITIONS['aminoacid'] = AMINOACIDS
    
    BACKBONE = DEFINITIONS['bb']


#==============================================================================
# PLANTERS
#==============================================================================

# protein
#==============================================================================

def setCalpha(ag, label):
    
    flags = ag._getNames() == 'CA'
    indices = flags.nonzero()[0]
    aliases = ALIASES['calpha']
    if len(indices):
        torf = array([rn in AMINOACIDS for rn in ag._getResnames()[indices]])
        flags[indices] = torf
        ag._setFlags(flags, *aliases)
        ag._setSubset(indices[torf], *aliases)
    else:
        ag._setFlags(flags, *aliases)
        ag._setSubset(array([]), *aliases)
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
    ag._setFlags(flags, *ALIASES['protein'])
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
    ag._setFlags(flags, *ALIASES[label])
    return flags

addPlanter(setBackbone, 'bb', 'backbone')
addPlanter(setBackbone, 'bbfull', 'backbonefull')

def setSidechain(ag, label):
    
    protein = ag._getSubset('protein')
    if len(protein):
        flags = zeros(ag.numAtoms(), bool)
        flags[protein] = True
        flags[ag._getSubset('bbfull')] = False
    else:
        flags = zeros(ag.numAtoms(), bool)
    ag._setFlags(flags, *ALIASES[label])
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
    ag._setFlags(flags, *ALIASES[label])
    return flags
    
addPlanter(setCategories, 'stdaa', 'nonstdaa', *CATEGORIZED.keys(), 
           aliases=False)


# hetero, nucleic, water, etc.
#==============================================================================

def setResiflag(ag, label):    
    
    resnames = DEFINITIONS[label]
    flags = array([rn in resnames for rn in ag._getResnames()], bool)
    ag._setFlags(flags, *ALIASES[label])
    return flags

addPlanter(setResiflag, 'nucleic', 'nucleobase', 'nucleoside', 'nucleotide',
           'water', 'ion', 'lipid', 'sugar', 'heme', 'at', 'cg', 'purine', 
           'pyrimidine', aliases=False)


def setHetero(ag, label):
    
    flags = ones(ag.numAtoms(), bool)
    flags[ag._getSubset('protein')] = False
    flags[ag._getSubset('nucleic')] = False
    ag._setFlags(flags, 'hetero')
    return flags

addPlanter(setHetero, 'hetero',)


# element
#==============================================================================

def setElement(ag, label):
    
    match = DEFINITIONS[label].match 
    flags = array([match(nm) is not None for nm in ag._getNames()])
    flags[ag._getSubset('ion')] = False
    ag._setFlags(flags, *ALIASES[label])
    return flags

addPlanter(setElement, 'hydrogen', 'carbon', 'nitrogen', 'oxygen', 'sulfur', 
           aliases=False)

def setNoh(ag, label):
    
    flags = ones(ag.numAtoms(), bool)
    flags[ag._getSubset('hydrogen')] = False
    ag._setFlags(flags, *ALIASES[label])
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
    ag._setFlags(flags, *ALIASES[label])
    return flags

addPlanter(setSecondary, *SECONDARY.keys(), aliases=False, 
           fields=['secondary']) 


#==============================================================================
# User functions
#==============================================================================

def flagDefinition(*arg, **kwarg):
    """Learn, change, or reset :ref:`flags` definitions.
    
    **Editable definitions**

    
    Call this function with no arguments:
    
    >>> flagDefinition()
    
    
    **Learn a definition**
    
    
    Pass an editable flag name:
    
    >>> flagDefinition('bb')
    ['C', 'CA', 'N', 'O']
    >>> flagDefinition('hydrogen')
    '[0-9]?H.*'

    **Change a definition**    
    

    Pass an editable flag name with its new definition:
    
    >>>
    
    Note that the type of the new definition must be the same as the type
    of the old definition.
    
    **Reset definitions**
    
    Pass *reset* keyword as follows to restore default definitions of editable
    flags.
    
    >>> flagDefinition(reset=True)
       
    """

    if arg and kwarg:
        raise ValueError('')
    elif arg:
        if len(arg) > 1:
            raise ValueError('specify only one flag label')
        arg = arg[0]
        try:
            definition = DEFINITIONS[arg]
        except KeyError:
            raise ValueError('{0:s} is not a flag label that can be '
                               'modified by the user')
        else:
            try:
                return definition.pattern
            except AttributeError:  
                definition = list(definition)
                definition.sort()
                return definition
    elif kwarg:
        pass
    else:
        return list(DEFINITIONS)

def addNonstdAA(resname, *props):
    """Add a non-standard amino acid. *resname* cannot be longer than four 
    characters."""
    
    pass
    
def delNonstdAA(resname, *props):
    """Add a non-standard amino acid. *resname* cannot be longer than four 
    characters."""
    
    pass

def getNonstdAAs(resname, *props):
    """Add a non-standard amino acid. *resname* cannot be longer than four 
    characters."""
    
    pass

def getKeywordResnames(keyword):
    """Return residue names associated with a keyword.
    
    >>> getKeywordResnames('acidic')
    ['ASP', 'GLU', 'PTR', 'SEP', 'TPO']"""
    
    assert isinstance(keyword, str), 'keyword must be a string instance'
    try:
        resnames = KEYWORD_RESNAMES[keyword]
        resnames.sort()
        return resnames  
    except KeyError:
        if keyword in KEYWORD_RESNAMES_READONLY:
            LOGGER.warn('{0:s} is defined as {1:s}'.format(repr(keyword), 
                                    repr(KEYWORD_RESNAMES_READONLY[keyword])))
        else:
            LOGGER.warn("{0:s} is not a keyword".format(repr(keyword)))

def setKeywordResnames(keyword, resnames):
    """Change the list of residue names associated with a keyword.  *keyword* 
    must be a string, and *resnames* may be a list, tuple, or set of strings. 
    The existing list of residue names will be overwritten with the given 
    residue names.  Note that changes in keyword definitions are not saved 
    permanently.
    
    >>> setKeywordResnames('acidic', ['ASP', 'GLU', 'PTR', 'SEP', 'TPO'])"""
    
    if not isinstance(keyword, str):
        raise TypeError('keyword must be a string')
    if not isinstance(resnames, (list, tuple, set)):
        raise TypeError('resnames must be a list, set, or tuple')
    if not areAllStrings(resnames):
        raise TypeError('all items in resnames must be string instances')
    
    if keyword in KEYWORD_RESNAMES_READONLY:
        LOGGER.warn("{0:s} is defined as {1:s} and cannot be changed directly"
            .format(repr(keyword), repr(KEYWORD_RESNAMES_READONLY[keyword])))
        return
    if keyword in KEYWORD_RESNAMES:
        for rn in resnames:
            if not isinstance(rn, str):
                raise TypeError('all items in resnames must be strings')
        KEYWORD_RESNAMES[keyword] = list(set(resnames))
        _setReadonlyResidueNames()
    else:
        raise ValueError("{0:s} is not a valid keyword".format(repr(keyword)))


def getElementRegex(name):
    """Return regular expression used for flagging elements based on atom 
    names.  See :ref:`element-flags` for the default list of flags.
    
    >>> getElementRegex('nitrogen')
    'N.*'"""
    
    try:
        return ELEMENT_REGEX[name]   
    except KeyError:
        LOGGER.warn('{0:s} is not a valid element flag name'.format(name))


def setElementRegex(name, regex):
    """Set regular expression used for flagging elements based on atom names.
    See :ref:`element-flags` for the default list of flags.
    
    >>> setElementRegex('nitrogen', 'N.*')"""
    

    if not name in ELEMENT_REGEX:
        raise ValueError('{0:s} is not a valid element flag name'
                           .format(repr(name)))
    try:
        ELEMENT_COMPILED[name] = recompile(regex)
    except Exception as err:
        raise ValueError("{0:s} is not a valid regular expression, {1:s}"
                         .format(repr(regex), str(err)))
    else:
        ELEMENT_REGEX[name] = regex
        SETTINGS['element_name_regex'] = ELEMENT_REGEX 
        timestamp()

def getBackboneNames(full=False):
    """Return protein backbone atom names.  When *full* is **True** argument
    returns atom names for *backbonefull* keyword.  See :ref:`subset-flags` 
    for details.
    
    >>> getBackboneNames()
    ['C', 'CA', 'N', 'O']"""
    
    if full:
        bban = list(BACKBONEFULL)
    else:
        bban = list(BACKBONE)
    bban.sort()
    return bban 

def setBackboneNames(names, full=False):
    """Set protein backbone atom names.  *names* must be a list of atom names.
    Setting *full* argument **True** will change atom names for *backbonefull* 
    flag.  See :ref:`subset-flags` for details.
    
    >>> setBackboneNames('CA C O N H H1 H2 H3 OXT'.split(), full=True)"""
    
    try:
        bb = set(names)
    except TypeError:
        raise TypeError('names must be a list of strings')
    else:
        if len(bb) != len(names):
            raise ValueError('names contains non-unique items')

    if full:    
        global BACKBONEFULL
        BACKBONEFULL = SETTINGS['flags_backbonefull'] = bb
    else:
        global BACKBONE
        BACKBONE = SETTINGS['flags_backbone'] = bb
    timestamp()

setDefinitions()
