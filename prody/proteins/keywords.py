# ProDy: A Python Package for Protein Structural Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
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
""":mod:`keywords` module defines methods to get and set definitions
of keywords used when making atom selections.

Functions:

  * Get:
    
    * :func:`getAromaticResidueNames`
    * :func:`getNucleicResidueNames`
    * :func:`getBackboneAtomNames`
    * :func:`getProteinResidueNames`
    * :func:`getHydrogenRegex`
    * :func:`getBasicResidueNames`
    * :func:`getAliphaticResidueNames`
    * :func:`getCyclicResidueNames`
    * :func:`getSmallResidueNames`
    * :func:`getWaterResidueNames`
    * :func:`getMediumResidueNames`
    * :func:`getAcidicResidueNames`
    
  * Set:
    
    * :func:`setAromaticResidueNames`
    * :func:`setNucleicResidueNames`
    * :func:`setBackboneAtomNames`
    * :func:`setProteinResidueNames`
    * :func:`setHydrogenRegex`
    * :func:`setBasicResidueNames`
    * :func:`setAliphaticResidueNames`
    * :func:`setCyclicResidueNames`
    * :func:`setSmallResidueNames`
    * :func:`setWaterResidueNames`
    * :func:`setMediumResidueNames`
    * :func:`setAcidicResidueNames`


"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

__all__ = ['getAromaticResidueNames', 'setAromaticResidueNames',
           'getNucleicResidueNames', 'setNucleicResidueNames',
           'getBackboneAtomNames', 'setBackboneAtomNames',
           'getProteinResidueNames', 'setProteinResidueNames',
           'getHydrogenRegex', 'setHydrogenRegex',
           'getBasicResidueNames', 'setBasicResidueNames',
           'getAliphaticResidueNames', 'setAliphaticResidueNames',
           'getCyclicResidueNames', 'setCyclicResidueNames',
           'getSmallResidueNames', 'setSmallResidueNames',
           'getWaterResidueNames', 'setWaterResidueNames',
           'getMediumResidueNames', 'setMediumResidueNames',
           'getAcidicResidueNames', 'setAcidicResidueNames',
           ]

import numpy as np

BACKBONE_ATOM_NAMES = ('CA', 'N', 'C', 'O') 
PROTEIN_RESIDUE_NAMES = ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 
        'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
        'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD', 'HSE', 'HSP')
        
WATER_RESIDUE_NAMES = ('HOH', 'WAT', 'TIP3', 'H2O')#, 'HH0', 'OHH', 'OH2', 'SOL', 'TIP', 'TIP2', 'TIP4')
NUCLEIC_RESIDUE_NAMES = ('GUA', 'ADE', 'CYT', 'THY', 'URA', 
                           'DA', 'DC', 'DG', 'DT')

HYDROGEN_REGEX = '[0-9]?H.*'
#HYDROGEN_REGEX_COMPILED = re.compile(HYDROGEN_REGEX)

ACIDIC_RESIDUE_NAMES = ('ASP', 'GLU')
BASIC_RESIDUE_NAMES = ('LYS', 'ARG', 'HIS', 'HSP')
ALIPHATIC_RESIDUE_NAMES = ('ALA', 'GLY', 'ILE', 'LEU', 'VAL')
AROMATIC_RESIDUE_NAMES = ('HIS', 'PHE', 'TRP', 'TYR')
SMALL_RESIDUE_NAMES = ('ALA', 'GLY', 'SER')
MEDIUM_RESIDUE_NAMES = ('VAL', 'THR', 'ASP', 'ASN', 'PRO', 'CYS')
HYDROPHOBIC_RESIDUE_NAMES = ()
CYCLIC_RESIDUE_NAMES = ('HIS', 'PHE', 'PRO', 'TRP', 'TYR')


def getBackboneAtomNames():
    """Return protein :term:`backbone` atom names."""
    return list(BACKBONE_ATOM_NAMES)

def setBackboneAtomNames(backbone_atom_names):
    """Set protein :term:`backbone` atom names."""
    if not isinstance(backbone_atom_names, (list, tuple)):
        raise TypeError('backbone_atom_names must be a list or a tuple')
    BACKBONE_ATOM_NAMES = tuple(backbone_atom_names)

def getProteinResidueNames():
    """Return :term:`protein` residue names."""
    return list(PROTEIN_RESIDUE_NAMES)

def setProteinResidueNames(protein_residue_names):
    """Set :term:`protein` residue names."""
    if not isinstance(protein_residue_names, (list, tuple)):
        raise TypeError('protein_residue_names must be a list or a tuple')
    PROTEIN_RESIDUE_NAMES = tuple(protein_residue_names)

def getAcidicResidueNames():
    """Return :term:`acidic` residue names."""
    return list(ACIDIC_RESIDUE_NAMES)

def setAcidicResidueNames(acidic_residue_names):
    """Set :term:`acidic` residue names."""
    if not isinstance(acidic_residue_names, (list, tuple)):
        raise TypeError('acidic_residue_names must be a list or a tuple')
    ACIDIC_RESIDUE_NAMES = tuple(acidic_residue_names)
    
def getBasicResidueNames():
    """Return :term:`basic` residue names."""
    return list(BASIC_RESIDUE_NAMES)

def setBasicResidueNames(basic_residue_names):
    """Set :term:`basic` residue names."""
    if not isinstance(basic_residue_names, (list, tuple)):
        raise TypeError('basic_residue_names must be a list or a tuple')
    BASIC_RESIDUE_NAMES = tuple(basic_residue_names)

def getAliphaticResidueNames():
    """Return :term:`aliphatic` residue names."""
    return list(ALIPHATIC_RESIDUE_NAMES)

def setAliphaticResidueNames(aliphatic_residue_names):
    """Set :term:`aliphatic` residue names."""
    if not isinstance(aliphatic_residue_names, (list, tuple)):
        raise TypeError('aliphatic_residue_names must be a list or a tuple')
    ALIPHATIC_RESIDUE_NAMES = tuple(aliphatic_residue_names)

def getAromaticResidueNames():
    """Return :term:`aromatic` residue names."""
    return list(AROMATIC_RESIDUE_NAMES)

def setAromaticResidueNames(aromatic_residue_names):
    """Set :term:`aromatic` residue names."""
    if not isinstance(aromatic_residue_names, (list, tuple)):
        raise TypeError('aromatic_residue_names must be a list or a tuple')
    AROMATIC_RESIDUE_NAMES = tuple(aromatic_residue_names)

def getSmallResidueNames():
    """Return :term:`small` residue names."""
    return list(SMALL_RESIDUE_NAMES)

def setSmallResidueNames(small_residue_names):
    """Set :term:`small` residue names."""
    if not isinstance(small_residue_names, (list, tuple)):
        raise TypeError('small_residue_names must be a list or a tuple')
    SMALL_RESIDUE_NAMES = tuple(small_residue_names)

def getMediumResidueNames():
    """Return :term:`medium` residue names."""
    return list(MEDIUM_RESIDUE_NAMES)

def setMediumResidueNames(medium_residue_names):
    """Set :term:`medium` residue names."""
    if not isinstance(medium_residue_names, (list, tuple)):
        raise TypeError('medium_residue_names must be a list or a tuple')
    MEDIUM_RESIDUE_NAMES = tuple(medium_residue_names)

def getCyclicResidueNames():
    """Return :term:`cyclic` residue names."""
    return list(CYCLIC_RESIDUE_NAMES)

def setCyclicResidueNames(cyclic_residue_names):
    """Set :term:`cyclic` residue names."""
    if not isinstance(cyclic_residue_names, (list, tuple)):
        raise TypeError('cyclic_residue_names must be a list or a tuple')
    CYCLIC_RESIDUE_NAMES = tuple(cyclic_residue_names)

def getWaterResidueNames():
    """Return :term:`water` residue names."""
    return list(WATER_RESIDUE_NAMES)

def setWaterResidueNames(water_residue_names):
    """Set :term:`water` residue names."""
    if not isinstance(water_residue_names, (list, tuple)):
        raise TypeError('water_residue_names must be a list or a tuple')
    WATER_RESIDUE_NAMES = tuple(water_residue_names)

def getNucleicResidueNames():
    """Return :term:`nucleic` residue names."""
    return list(NUCLEIC_RESIDUE_NAMES)

def setNucleicResidueNames(nucleic_residue_names):
    """Set :term:`nucleic` residue names."""
    if not isinstance(water_residue_names, (list, tuple)):
        raise TypeError('nucleic_residue_names must be a list or a tuple')
    NUCLEIC_RESIDUE_NAMES = tuple(nucleic_residue_names)

def getHydrogenRegex():
    """Return regular expression to match :term:`hydrogen` atom names."""
    return list(HYDROGEN_REGEX)

def setHydrogenRegex(hydrogen_regex):
    """Set regular expression to match :term:`hydrogen` atom names."""
    if not isinstance(hydrogen_regex, (str)):
        raise TypeError('hydrogen_regex must be a string')
    HYDROGEN_REGEX = hydrogen_regex
    #HYDROGEN_REGEX_COMPILED = re.compile(hydrogen_regex)
