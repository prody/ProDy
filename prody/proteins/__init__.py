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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np

DTYPES = {'coordinates': np.float64, 

          'atomnames': '|S6', 
          'altlocs': '|S1',
          'anisou': np.float64, 
          'chainids': '|S1', 
          'elements': '|S2', 
          'hetero': np.bool, 
          'occupancies': np.float64, 
          'resnames': '|S6',
          'resnums': np.int64, 
          'secondary': '|S1',
          'segnames': '|S6', 
          'siguij': np.float64,
          'bfactors': np.float64,
          'icodes': '|S1',
          
          'atomtypes': '|S6',
          'charges': np.float64, 
          'masses': np.float64, 
          'radii': np.float64, 
}

class Field(object):
    __slots__ = ('_name', '_var', '_dtype',  '_doc', '_doc_pl', '_meth', '_meth_pl')
    def __init__(self, name, var, dtype, doc, meth, doc_pl=None, meth_pl=None):
        self._name = name
        self._var = var
        self._dtype = dtype
        self._doc = doc
        if doc_pl is None:
            self._doc_pl = doc + 's'
        else:
            self._doc_pl = doc_pl
        self._meth = meth
        if meth_pl is None:
            self._meth_pl = meth + 's'
        else:
            self._meth_pl = meth_pl
        
    def name(self):
        return self._name
    name = property(name, doc='Atomic data field name.')
    def var(self):
        return self._var
    var = property(var, doc='Atomic data variable name.')
    def dtype(self):
        return self._dtype
    dtype = property(dtype, doc='Atomic data type (NumPy types).')
    def doc(self):
        return self._doc
    doc = property(doc, doc='Atomic data field name to be used in documentation.')
    def doc_pl(self):
        return self._doc_pl
    doc_pl = property(doc_pl, doc='Atomic data field name in plural form to be used in documentation.')
    def meth(self):
        return self._meth
    meth = property(meth, doc='Atomic data field name to be used in get/set methods.')
    def meth_pl(self):
        return self._meth_pl
    meth_pl = property(meth_pl, doc='Atomic data field name in plural form to be used in get/set methods.')

ATOMIC_DATA_FIELDS = {
    'name':      Field('name',      'names',       '|S6',      'atom name',                      'AtomName'),
    'altloc':    Field('altloc',    'altlocs',     '|S1',      'alternate location indicator',   'AltLocIndicator'),
    'anisou':    Field('anisou',    'anisou',      np.float64, 'anisotropic temperature factor', 'AnisoTempFactor'),
    'chain':     Field('chain',     'chids',       '|S1',      'chain identifier',               'ChainIdentifier'),
    'element':   Field('element',   'elements',    '|S2',      'element symbol',                 'ElementSymbol'),
    'hetero':    Field('hetero',    'hetero',      np.bool,    'hetero flag',                    'HeteroFlag'),
    'occupancy': Field('occupancy', 'occupancies', np.float64, 'occupancy value',                'Occupancy',      meth_pl='Occupancies'),
    'resname':   Field('resname',   'resnames',    '|S6',      'residue name',                   'ResidueName'),
    'resnum':    Field('resnum',    'resnums',     np.int64,   'residue number',                 'ResidueNumber'),
    'secondary': Field('secondary', 'secondary',   '|S1',      'secondary structure assignment', 'SecondaryStr'),
    'segment':   Field('segment',   'segments',    '|S6',      'segment name',                   'SegmentName'),
    'siguij':    Field('siguij',    'siguij',      np.float64, 'standard deviations for the anisotropic temperature factor',                   
                                                                                                 'AnisoStdDev'),
    'beta':      Field('beta',      'bfactors',    np.float64, 'temperature (B) factor',         'TempFactor'),
    'icode':     Field('icode',     'icodes',      '|S1',      'insertion code',                 'InsertionCode'),
    'type':      Field('type',      'types',       '|S6',      'atom type',                      'AtomType'),
    'charge':    Field('charge',    'charges',     np.float64, 'atomic partial charge',          'Charge'),
    'mass':      Field('mass',      'masses',      np.float64, 'atomic mass',                    'Mass', 'atomic masses', 'Masses'),
    'radius':    Field('radius',    'radii',       np.float64, 'atomic radius',                  'Radius', 'atomic radii', 'Radii'),
}

DATA_FIELD_ALIASES = {}

__all__ = []

import prodb
from prodb import *  
__all__ += prodb.__all__

from . import keywords
from keywords import *
__all__ += keywords.__all__

from . import measure
from measure import *
__all__ += measure.__all__

import select

from . import atomgroup 
from atomgroup import *
__all__ += atomgroup.__all__

from . import compare
from compare import *
__all__ += compare.__all__
