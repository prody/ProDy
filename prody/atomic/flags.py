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

import re as RE
from time import time
from textwrap import wrap as textwrap
from collections import defaultdict

from numpy import array, ones, zeros

from prody import SETTINGS
from prody.utilities import tabulate, wrap

__all__ = ['addNonstdAA', 'delNonstdAA', 'getNonstdAAs',
           'addFlagResname', 'delFlagResname', 'getFlagResnames']

ALIASES = {}
PLANTERS = {}
FIELDS = defaultdict(set)
FIELDSDEFAULT = ['name', 'resname', 'resnum']
TIMESTAMP = None

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

CATEGORIES = {
    'aromaticity': set(['aromatic', 'aliphatic']), 
    'cyclic': set(['acyclic', 'cyclic']),
    'charge': set(['acidic', 'basic', 'neutral']),
    'depth': set(['buried', 'surface']),
    'hydrophobicity': set(['hydrophobic', 'polar']),
    'size': set(['small', 'medium', 'large']),
}


STANDARDAA = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
              'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
              'TYR', 'VAL']

__doc__ += """

Protein
===============================================================================

**protein** and **aminoacid** flags indicate the twenty standard amino 
acids and some non-standard amino acids described below.  A residue is 
considered an amino acid if its name is one of the following abbreviations
and it has an atom named ``'CA'``.
  
Amino acids
-------------------------------------------------------------------------------

""" + wrap('**stdaa** flag indicates the standard amino acid residues: `' + 
           '`_, `'.join(STANDARDAA) + '`_') + '\n\n'

for resi in STANDARDAA:
    __doc__ += PDBLIGSUM.format(resi)

__doc__ += """

Non-standard
-------------------------------------------------------------------------------

**nonstdaa** flag indicates one of the following residues:

==========  ===================================================================
`ASX`_ (B)  asparagine or aspartic acid
`GLX`_ (Z)  glutamine or glutamic acid
`CSO`_ (C)  S-hydroxycysteine
`HIP`_ (H)  ND1-phosphohistidine
 HSD   (H)  prototropic tautomer of histidine, H on ND1 (CHARMM force field)
 HSE   (H)  prototropic tautomer of histidine, H on NE2 (CHARMM force field)       
 HSP   (H)  protonated histidine 
`SEC`_ (U)  selenocysteine
`SEP`_ (S)  phosphoserine
`TPO`_ (T)  phosphothreonine
`PTR`_ (Y)  O-phosphotyrosine
 XLE   (J)  leucine or isoleucine 
 XAA   (X)  unspecified or unknown 
==========  ===================================================================

Some of these residues do not correspond to amino acids in the *PDB*, but 
flagging them as amino acids when they have Cα atoms prevents conflicts.

You can modify the list of non-standard amino acids using :func:`addNonstdAA`,
:func:`delNonstdAA`, and :func:`getNonstdAAs` functions.

"""

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
    'PTR': set(['cyclic', 'surface', 'acidic', 'large', 'polar']),
    'XLE': set(['aliphatic', 'acyclic', 'buried', 'hydrophobic', 'large']),
    'XAA': set(),
}

for resi in NONSTANDARD:
    __doc__ += PDBLIGSUM.format(resi)

__doc__ += """

Subsets of atoms
-------------------------------------------------------------------------------

Following flags indicate well defined subsets of protein atoms: 

==================  ===========================================================
**calpha**          Cα atoms of protein residues, same as 
                    ``'name CA and protein'``
**ca**              same as above
**backbone**        non-hydrogen backbone atoms of protein residues, same as
                    ``'name CA C O N and protein'``
**bb**              same as above
**backbonefull**    backbone atoms of protein residues, same as
                    ``'name CA C O N H H1 H2 H3 OXT and protein'``
**bbfull**          same as above 
**sidechain**       side-chain atoms of protein residues, same as
                    ``'not backbone'``
**sc**              same as above
==================  ===========================================================
"""


CATEGORIZED = {

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
}

_ = SETTINGS.get('nonstandard')
if _ is not None:
    NONSTANDARD = _

for resi, props in NONSTANDARD.iteritems():
    for prop in props: 
        CATEGORIZED[prop].add(resi)

CATEGORIZED['charged'] = set(CATEGORIZED['acidic'])
CATEGORIZED['charged'].update(CATEGORIZED['basic'])

__doc__ += """

Side-chain properties
-------------------------------------------------------------------------------

Amino acids are categorized (or flagged) as follows based on the 
physicochemical properties of their side-chains:

"""

cats = CATEGORIZED.keys()
cats.sort()

col1 = []
col2 = []
for cat in cats:
    col1.append('**' + cat + '**')
    resnames = list(CATEGORIZED[cat])
    resnames.sort()
    col2.append(', '.join(resnames))
__doc__ += tabulate(col1, col2, header=False)

'''
maxlen = max(map(len, cats)) + 6
wrapat = 79 - maxlen
table = '=' * (maxlen - 2) + '  ' + '=' * wrapat + '\n'
__doc__ += table
for cat in cats:
    __doc__ += ('**' + cat + '**').ljust(maxlen)
    resnames = list(CATEGORIZED[cat])
    resnames.sort()
    for i, line in enumerate(textwrap(', '.join(resnames), wrapat)): 
        __doc__ += ('\n' + ' ' * maxlen if i else '') + line
    __doc__ += '\n'
__doc__ += table
'''
AMINOACIDS = set()
AMINOACIDS.update(STANDARDAA)
AMINOACIDS.update(NONSTANDARD.keys())

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

BACKBONE = set(['CA', 'C', 'O', 'N'])
BACKBONEFULL = set(['CA', 'C', 'O', 'N', 'H', 'H1', 'H2', 'H3', 'OXT'])

def setBackbone(ag, label):
    
    protein = ag._getSubset('protein')
    if len(protein):
        bb = BACKBONEFULL if label.endswith('full') else BACKBONE
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

CATEGORIZED['stdaa'] = STANDARDAA
CATEGORIZED['nonstdaa'] = NONSTANDARD

def setCategories(ag, label):
    
    calpha = ag._getSubset('ca')
    if len(calpha):
        resnames = CATEGORIZED[label]
        residx = ag._getResindices()
        torf = zeros(ag.numResidues(), bool)
        torf[residx[calpha]] = [rn in resnames 
                                for rn in ag._getResnames()[calpha]]
        flags = torf[residx]
    else:
        flags = zeros(ag.numAtoms(), bool)
    ag._setFlags(flags, *ALIASES[label])
    return flags
    
addPlanter(setCategories, *CATEGORIZED.keys(), aliases=False)


"""Among these list of backbone atom names can be changed using 
:func:`setBackboneAtomNames`  and regular expressions for element types
can be changed using :func:`setAtomNameRegex`.

Below functions can be used to learn and change the definitions of 
some selection keywords:

  * Learn keyword definitions:
    
    * :func:`getAtomNameRegex`
    * :func:`getBackboneAtomNames`
    * :func:`getKeywordResnames` 
    
  * Change keyword definitions:
    
    * :func:`setAtomNameRegex`
    * :func:`setBackboneAtomNames`
    * :func:`setKeywordResnames`
"""

BACKBONE_ATOM_NAMES = set(('CA', 'N', 'C', 'O'))
BACKBONE_FULL_ATOM_NAMES = set(('CA', 'N', 'C', 'O', 
                                'H', 'H1', 'H2', 'H3', 'OXT'))

RESIFLAGS = {}

__doc__ += """

Nucleic
===============================================================================

**nucleic** flag indicates nucleobases, nucleotides, and some nucleoside
derivatives described below.

Nucleobases
-------------------------------------------------------------------------------

**nucleobase** flag indicates following residues that are recognized by the
*PDB*:

=======  ======================================================================
`ADE`_   adenine    
`GUN`_   guanine    
`CYT`_   cytosine   
`THY`_   thymine    
`URA`_   uracil     
=======  ======================================================================

"""

NUCLEIC = set()
RESIFLAGS['nucleobase'] = set(['GUN', 'ADE', 'CYT', 'THY', 'URA'])
NUCLEIC.update(RESIFLAGS['nucleobase'])

for resi in RESIFLAGS['nucleobase']:
    __doc__ += PDBLIGSUM.format(resi)

__doc__ += """
Nucleotides
-------------------------------------------------------------------------------

**nucleotide** flag indicates residues with the following names:

=======  ======================================================================
`DA`_    2'-deoxyadenosine-5'-monophosphate
`DC`_    2'-deoxycytidine-5'-monophosphate
`DG`_    2'-deoxyguanosine-5'-monophosphate
`DT`_    2'-deoxythymidine-5'-monophosphate
`DU`_    2'-deoxyuridine-5'-monophosphate
`A`_     adenosine-5'-monophosphate
`C`_     cytidine-5'-monophosphate
`G`_     guanosine-5'-monophosphate
`T`_     2'-deoxythymidine-5'-monophosphate (superceded by DT)
`U`_     uridine-5'-monophosphate
=======  ======================================================================

"""

RESIFLAGS['nucleotide'] = set(['DA', 'DC', 'DG', 'DT', 'DU', 
                               'A', 'C', 'G', 'T', 'U'])
NUCLEIC.update(RESIFLAGS['nucleotide'])

for resi in RESIFLAGS['nucleotide']:
    __doc__ += PDBLIGSUM.format(resi)

__doc__ += """

Nucleosides
-------------------------------------------------------------------------------

**nucleoside** flag indicates following nucleoside derivatives that are
recognized by *PDB*:

=======  ======================================================================
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
=======  ======================================================================

"""

RESIFLAGS['nucleoside'] = set(('AMP ADP ATP CDP CTP GMP GDP GTP TMP TTP '
                               'UMP UDP UTP').split()) 
NUCLEIC.update(RESIFLAGS['nucleoside'])

for resi in RESIFLAGS['nucleoside']:
    __doc__ += PDBLIGSUM.format(resi)


__doc__ += """

Others
-------------------------------------------------------------------------------

==============  ===============================================================
**at**          same as ``'resname ADE A THY T'``
**cg**          same as ``'resname CYT C GUN G'``
**purine**      same as ``'resname ADE A GUN G'``
**pyrimidine**  same as ``'resname CYT C THY T URA U'``
==============  ===============================================================
"""

RESIFLAGS['at'] = set(('ADE A THY T').split()) 
RESIFLAGS['cg'] = set(('CYT C GUN G').split())
RESIFLAGS['purine'] = set(('ADE A GUN G').split())
RESIFLAGS['pyrimidine'] = set(('CYT C THY T URA U').split())
RESIFLAGS['nucleic'] = NUCLEIC

def setHetero(ag, label):
    
    flags = ones(ag.numAtoms(), bool)
    flags[ag._getSubset('protein')] = False
    flags[ag._getSubset('nucleic')] = False
    ag._setFlags(flags, 'hetero')
    return flags

addPlanter(setHetero, 'hetero',)

__doc__ += """
Heteros
===============================================================================

**hetero** flag indicates anything other than a protein or a nucleic residue,
i.e. ``'not (protein or nucleic)'``.

In addition, **hetatm** flag will be if atomic data is parsed from a PDB or a 
similar file.  This flag will indicate atoms that are marked ``'HETATM'`` in
the data file. 

Water
-------------------------------------------------------------------------------

**water** flag definition includes `HOH`_ and `DOD`_ recognized by *PDB* and 
also WAT, TIP3, H2O, OH2, TIP, TIP2, and TIP4 recognized by some molecular 
dynamics (MD) force fields.  

.. _HOH: http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId=HOH
.. _DOD: http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId=DOD

Previously used water types HH0, OHH, and SOL conflict with other compounds in
the *PDB*, so are removed from the definition of this flag.
              
"""

RESIFLAGS['water'] = set(['HOH', 'DOD', 'WAT', 'TIP3', 'H2O', 'OH2', 
                          'TIP', 'TIP2', 'TIP4'])

resnames = ('AL BA CA CD CL CO CS CU CU1 CUA HG IN IOD K MG MN3 NA PB PT RB '
            'TB TL WO4 YB ZN CAL CES CLA POT SOD ZN2').split()
RESIFLAGS['ion'] = set(resnames)
resnames.sort()

__doc__ += """
Ions
-------------------------------------------------------------------------------

**ion** flag definition includes the following ions most of which are 
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

Ion identifiers that are obsoleted by *PDB* (MO3, MO4, MO5, MO6, NAW, OC7, and 
ZN1) are removed from this definition.

"""

for resi in resnames:
    __doc__ += PDBLIGSUM.format(resi)


resnames = ('GPE LPP SDS STE OLA').split()
RESIFLAGS['lipid'] = set(resnames)
resnames.sort()

__doc__ += """

Lipid
-------------------------------------------------------------------------------

""" + wrap('**lipid** flag indicates `' + '`_, `'.join(resnames) + 
           '`_ from *PDB* and LPPC DLPE DMPC POPC POPE PALM '
           'OLEO STEA PCGL from *CHARMM* force field.') + '\n\n'
resnames.extend('LPPC DLPE DMPC POPC POPE PALM OLEO STEA PCGL'.split())
for resi in resnames:
    __doc__ += PDBLIGSUM.format(resi)

resnames = 'GLC GLO BGC'.split()
RESIFLAGS['sugar'] = set(resnames)
resnames.sort()

__doc__ += """

Sugar
-------------------------------------------------------------------------------

""" + wrap('**sugar** flag indicates `' + '`_, `'.join(resnames) + 
           '`_ from *PDB* and AGLC from *CHARMM*.') + '\n\n'

for resi in resnames:
    __doc__ += PDBLIGSUM.format(resi)

resnames = ('1FH 2FH DDH DHE HAS HDD HDE HDM HEA HEB HEC HEM HEO HES HEV NTE '
            'SRM VER').split()
RESIFLAGS['heme'] = set(resnames)
resnames.sort()

__doc__ += """

Heme
-------------------------------------------------------------------------------

""" + wrap('**heme** flag indicates `' + '`_, `'.join(resnames) + 
           '`_ from *PDB*.') + '\n\n'

for resi in resnames:
    __doc__ += PDBLIGSUM.format(resi)


def setResiflag(ag, label):    
    
    resnames = RESIFLAGS[label]
    flags = array([rn in resnames for rn in ag._getResnames()], bool)
    ag._setFlags(flags, *ALIASES[label])
    return flags

for label in RESIFLAGS.keys():
    addPlanter(setResiflag, label)

__doc__ += """

Elements
===============================================================================

Following elements are recognized by applying regular expressions to atom 
names:

============  =================================================================
**carbon**    carbon atoms, same as ``'name "C.*" and not ion'``
**hydrogen**  hydrogen atoms, same as ``'name "[1-9]?H.*" and not ion'``
**noh**       non hydrogen atoms, same as ``'not hydrogen``
**nitrogen**  nitrogen atoms, same as ``'name "N.*" and not ion'``
**oxygen**    oxygen atoms, same as ``'name "O.*" and not ion'``
**sulfur**    sulfur atoms, same as ``'name "S.*" and not ion'``
============  =================================================================

"""

NAME_REGEX = {
    'carbon': RE.compile('C.*'),
    'hydrogen': RE.compile('[0-9]?H.*'),
    'nitrogen': RE.compile('N.*'),
    'oxygen': RE.compile('O.*'),
    'sulfur': RE.compile('S.*'),
}

def setElement(ag, label):
    
    match = NAME_REGEX[label].match 
    flags = array([match(nm) is not None for nm in ag._getNames()])
    flags[ag._getSubset('ion')] = False
    ag._setFlags(flags, *ALIASES[label])
    return flags

addPlanter(setElement, *NAME_REGEX.keys(), aliases=False)

def setNoh(ag, label):
    
    flags = ones(ag.numAtoms(), bool)
    flags[ag._getSubset('hydrogen')] = False
    ag._setFlags(flags, *ALIASES[label])
    return flags

addPlanter(setNoh, 'noh')

__doc__ += """

Structure
===============================================================================

Following secondary structure flags are defined but before they can be used,
secondary structure assignments must be made. 

==============  ===============================================================
**extended**    extended conformation, same as ``'secondary E'``
**helix**       α-helix conformation, same as ``'secondary H'``
**helix310**    3_10-helix conformation, same as ``'secondary G'``
**helixpi**     π-helix conformation, same as ``'secondary I'``
**turn**        hydrogen bonded turn conformation, same as ``'secondary T'``
**bridge**      isolated beta-bridge conformation, same as ``'secondary B'``
**bend**        bend conformation, same as ``'secondary S'``
**coil**        not in one of above conformations, same as ``'secondary C'``
==============  ===============================================================

"""

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

# definitions

def setDefinitions():
    global TIMESTAMP
    TIMESTAMP = time()

setDefinitions()


FLAGDOCS = ("``'" + "'``, ``'".join(RESIFLAGS.keys()) + 
            "'`` are acceptable *flag* arguments.")

    
def addFlagResname(flag, resname):
    """Add *resname* to *flag* definition."""
    
    pass

def delFlagResname(flag, resname):
    """Remove *resname* from *flag* definition."""
    
    pass

def getFlagResnames(flag):
    """Return list of residue names associated with *flag*."""
    
    pass

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

def getAtomNameRegex(name):
    """Return regular expression used for selecting common elements.
    
    >>> getAtomNameRegex('nitrogen')
    'N.*'"""
    
    assert isinstance(name, str), 'name must be a string instance'
    try:
        return KEYWORD_NAME_REGEX[name].pattern   
    except KeyError:
        LOGGER.warn('{0:s} is not a valid element'.format(name))

def setAtomNameRegex(name, regex):
    """Set regular expression used for selecting common elements.  Note that 
    changes in keyword definitions are not saved permanently.
    
    >>> setAtomNameRegex('nitrogen', 'N.*')"""
    
    assert isinstance(name, str), 'name must be a string instance'
    if not name in KEYWORD_NAME_REGEX:
        raise ValueError("{0:s} is not a valid keyword".format(repr(name)))
    if not isinstance(regex, str):
        raise TypeError("regex must be a string instance")
    try:
        regex = RE.compile(regex)
    except:
        raise ValueError("{0:s} is not a valid regular expression"
                         .format(repr(regex)))
    else:
        KEYWORD_NAME_REGEX[name] = regex

def getBackboneAtomNames(full=False):
    """Return protein backbone atom names.  ``full=True`` argument returns 
    atom names for *backbonefull* keyword.
    
    >>> getBackboneAtomNames()
    ['C', 'CA', 'N', 'O']"""
    
    assert isinstance(full, bool), 'full must be a boolean instance'
    if full:
        bban = list(BACKBONE_FULL_ATOM_NAMES)
    else:
        bban = list(BACKBONE_ATOM_NAMES)
    bban.sort()
    return bban 

def setBackboneAtomNames(backbone_atom_names, full=False):
    """Set protein backbone atom names.  Atom names for *backbonefull* keyword 
    can be set by passing ``full=True`` argument.  Note that changes in keyword
    definitions are not saved permanently."""
    
    if not isinstance(backbone_atom_names, (list, tuple, set)):
        raise TypeError('backbone_atom_names must be a list, tuple, or set')
    if not areAllStrings(backbone_atom_names):
        raise TypeError('all items in backbone_atom_names must be string '
                        'instances')
    assert isinstance(full, bool), 'full must be a boolean instance'
    if full:    
        global BACKBONE_FULL_ATOM_NAMES
        BACKBONE_FULL_ATOM_NAMES = set(backbone_atom_names)
    else:
        global BACKBONE_ATOM_NAMES
        BACKBONE_ATOM_NAMES = set(backbone_atom_names)
    _buildKeywordMap()
