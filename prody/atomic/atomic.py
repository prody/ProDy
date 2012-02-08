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

"""This module defines atomic data fields and base classes."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'


__all__ = ['Atomic', 'MultiCoordset']

isMacro = lambda none: None
isKeyword = lambda none: None
isReserved = lambda none: None
SELECT = None

READONLY = set(['numbonds', 'resindex', 'chindex', 'segindex'])

class Field(object):
    
    """Atomic data field."""
    
    __slots__ = ['_name', '_var', '_dtype',  '_doc', '_doc_pl', 
                 '_meth', '_meth_pl', '_ndim', '_none', '_selstr',
                 '_depr', '_depr_pl', '_synonym', '_readonly', '_call']
                 
    def __init__(self, name, dtype, **kwargs):
        self._name = name
        self._dtype = dtype
        self._var = kwargs.get('var', name + 's')
        self._doc = kwargs.get('doc', name)
        self._ndim = kwargs.get('ndim', 1)
        self._meth = kwargs.get('meth', name.capitalize())
        self._doc_pl = kwargs.get('doc_pl', self._doc + 's')
        self._meth_pl = kwargs.get('meth_pl', self._meth + 's')
        self._none = kwargs.get('none')
        self._selstr = kwargs.get('selstr')
        self._depr = kwargs.get('depr')
        if self._depr is None:
            self._depr_pl = None
        else:
            self._depr_pl = kwargs.get('depr_pl', self._depr + 's')
        self._synonym = kwargs.get('synonym')
        self._readonly = kwargs.get('readonly', False)
        self._call = kwargs.get('call', None)
        
    def name(self):
        return self._name
    
    name = property(name, doc='Data field name used in atom selections.')
    
    def var(self):
        return self._var
    
    var = property(var, doc='Internal variable name.')
    
    def dtype(self):
        return self._dtype
    
    dtype = property(dtype, doc='Data type (primitive Python types).')
    
    def doc(self):
        return self._doc
    
    doc = property(doc, doc='Data field name, as used in documentation.')
    
    def doc_pl(self):
        return self._doc_pl
    
    doc_pl = property(doc_pl, doc='Plural form for documentation.')
    
    def meth(self):
        return self._meth
    
    meth = property(meth, doc='Atomic get/set method name.')
    
    def meth_pl(self):        
        return self._meth_pl
    
    meth_pl = property(meth_pl, doc='get/set method name in plural form.')
    
    def ndim(self):
        return self._ndim
    
    ndim = property(ndim, doc='Expected number of data array dimensions.')
    
    def none(self):
        return self._none
    
    none = property(none, doc='When to set the value of the variable to None.')
    
    def selstr(self):
        return self._selstr
    
    selstr = property(selstr, doc='Selection string examples.')
    
    def synonym(self):
        return self._synonym
    
    synonym = property(synonym, doc='Synonym used in atom selections.')
        
    def readonly(self):
        return self._readonly
    
    readonly = property(readonly, 
                        doc='Read-only attribute without a set method.')    
    def call(self):
        return self._call
    
    call = property(call, doc='List of AtomGroup methods to call.')    
    
    def depr(self):
        return self._depr
    
    depr = property(depr, doc='Deprecated method name.')
    
    def depr_pl(self):
        return self._depr_pl
    
    depr_pl = property(depr_pl, doc='Deprecated method name in plural form.')
    
    def getDocstr(self, meth, plural=True, selex=True):
        """Return documentation string for the field."""
        
        assert meth in ('set', 'get', '_get'), "meth must be 'set' or 'get'"
        assert isinstance(plural, bool), 'plural must be a boolean'
        assert isinstance(selex, bool), 'selex must be a boolean'
        
        if meth == 'get':
            if plural:
                docstr = 'Return a copy of {0:s}.'.format(self.doc_pl)
            else:
                docstr = 'Return {0:s} of the atom.'.format(self.doc)
        elif meth == 'set':
            if plural:
                docstr = 'Set {0:s}.'.format(self.doc_pl)
            else:
                docstr = 'Set {0:s} of the atom.'.format(self.doc)
        else:
            selex = False
            if plural:
                docstr = 'Return {0:s} array.'.format(self.doc_pl) 
            
        selstr = self.selstr
        if selex and selstr:
            if plural:
                doc = self.doc_pl
            else:
                doc = self.doc
            if '(' in doc:
                doc = doc[:doc.index('(')]
            selex = "'``, ``'".join(selstr)
            selex = ("  {0:s} can be used in atom selections, e.g. "
                     "``'{1:s}'``.").format(doc.capitalize(), selex)
            if self.synonym is not None:
                selex = selex + ('  Note that *{0:s}* is a synonym for '
                    '*{1:s}*.').format(self.synonym, self.name)
            return docstr + selex
        else:
            return docstr


ATOMIC_DATA_FIELDS = {
    'name':      Field('name', '|S6', selstr=('name CA CB',), depr='AtomName'),
    'altloc':    Field('altloc', '|S1', doc='alternate location indicator', 
                       selstr=('altloc A B', 'altloc _'), 
                       depr='AltLocIndicator'),
    'anisou':    Field('anisou', float, doc='anisotropic temperature factor', 
                       ndim=2, depr='AnisoTempFactor'),
    'chain':     Field('chain', '|S1', var='chids', doc='chain identifier', 
                       meth='Chid', none='hv', synonym='chid', 
                       selstr=('chain A', 'chid A B C', 'chain _'), 
                       depr='ChainIdentifier'),
    'element':   Field('element', '|S2', doc='element symbol', 
                       selstr=('element C O N',), depr='ElementSymbol'),
    'hetero':    Field('hetero', bool, doc='hetero flag', 
                       selstr=('hetero', 'hetero and not water'), 
                       depr='HeteroFlag'),
    'occupancy': Field('occupancy', float, var='occupancies', 
                       doc='occupancy value', meth_pl='Occupancies',
                       selstr=('occupancy 1', 'occupancy > 0')),
    'resname':   Field('resname', '|S6', doc='residue name', 
                       selstr=('resname ALA GLY',), depr='ResidueName'),
    'resnum':    Field('resnum', int, doc='residue number', none='hv',
                       selstr=('resnum 1 2 3', 'resnum 120A 120B', 
                               'resnum 10 to 20', 'resnum 10:20:2', 
                               'resnum < 10'), synonym='resid',
                       depr='ResidueNumber'),
    'secondary': Field('secondary', '|S1', var='secondaries', 
                       doc='secondary structure assignment', 
                       meth='Secstr', synonym='secstr',
                       selstr=('secondary H E', 'secstr H E'),  
                       depr='SecondaryStr'),
    'segment':   Field('segment', '|S6', doc='segment name', meth='Segname',
                       selstr=('segment PROT', 'segname PROT'), 
                       synonym='segname', depr='SegmentName'),
    'siguij':    Field('siguij', float, doc='standard deviations for '
                       'anisotropic temperature factor', meth='Anistd', ndim=2, 
                       depr='AnisoStdDev'),
    'serial':    Field('serial', int, doc='serial number (from file)', 
                       doc_pl='serial numbers (from file)', none='sn2i', 
                       selstr=('serial 1 2 3', 'serial 1 to 10', 
                       'serial 1:10:2', 'serial < 10'), depr='SerialNumber'),
    'beta':      Field('beta', float, doc='β-value (temperature factor)', 
                       doc_pl='β-values (or temperature factors)', 
                       selstr=('beta 555.55', 'beta 0 to 500', 'beta 0:500', 
                       'beta < 500'), depr='TempFactor'),
    'icode':     Field('icode', '|S1', doc='insertion code', none='hv', 
                       selstr=('icode A', 'icode _'), depr='InsertionCode'),
    'type':      Field('type', '|S6', selstr=('type CT1 CT2 CT3',), 
                       depr='AtomType'),
    'charge':    Field('charge', float, doc='partial charge',  
                       selstr=('charge 1', 'abs(charge) == 1', 'charge < 0')),
    'mass':      Field('mass', float, var='masses', doc_pl='masses', 
                       meth_pl='Masses', selstr=('12 <= mass <= 13.5',)),
    'radius':    Field('radius', float, var='radii', doc='radius',  
                       doc_pl='radii', meth_pl='Radii', 
                       selstr=('radii < 1.5', 'radii ** 2 < 2.3')),
    'resindex':  Field('resindex', int, var='resindices', doc='residue index',  
                       doc_pl='residue indices', meth_pl='Resindices',
                       selstr=('resindex 0'), readonly=True, 
                       call=['getHierView']),
    'chindex':   Field('chindex', int, var='chindices', doc='chain index',  
                       doc_pl='chain indices', meth_pl='Chindices',
                       selstr=('chindex 0'), readonly=True, 
                       call=['getHierView']),
    'segindex':  Field('segindex', int, var='segindices', doc='segment index',  
                       doc_pl='segment indices', meth_pl='Segindices',
                       selstr=('segindex 0'), readonly=True, 
                       call=['getHierView']),
}


ATOMIC_ATTRIBUTES = {}
for field in ATOMIC_DATA_FIELDS.values():
    ATOMIC_ATTRIBUTES[field.var] = field


def wrapGetMethod(fn):
    def getMethod(self):
        return fn(self)
    return getMethod


def wrapSetMethod(fn):
    def setMethod(self, data):
        return fn(self, data)
    return setMethod


class Atomic(object):
    
    """Base class for all atomic classes.  This class can be used for type
    checking:
    
    >>> from prody import *
    >>> ag = parsePDB('1aar')
    >>> isinstance(ag, Atomic)
    True
    >>> prot = ag.protein
    >>> isinstance(prot, Atomic)
    True"""
    
    def __getattribute__(self, name):
        
        try:
            return object.__getattribute__(self, name)

        except AttributeError:
            selstr = name
            items = name.split('_')
            word = items[0]
            if (isKeyword(word) or items == 'not' or isMacro(word)):
                selstr = ' '.join(items)
                return SELECT.select(self, selstr)

        raise AttributeError("'{0:s}' object has no attribute '{1:s}' "
                             "and '{2:s}' is not a valid selection string"
                             .format(self.__class__.__name__, name, selstr))
    
    def __setattr__(self, name, value):
        
        if isReserved(name):
            raise AttributeError("'{0:s}' is a reserved name".format(name))
        if isMacro(name):
            raise AttributeError("'{0:s}' is a selection macro".format(name))

        self.__dict__[name] = value

    def select(self, selstr, **kwargs):
        """Return atoms matching the criteria in *selstr*.
        
        .. seealso:: :meth:`~select.Select.select()` for more usage 
           details."""
        
        return SELECT.select(self, selstr, **kwargs)


class MultiCoordset(object):
    
    """Base class for all atomic classes handling multiple coordinate sets.
    This class can be used for type checking:
    
    >>> from prody import *
    >>> ag = parsePDB('1aar')
    >>> isinstance(ag, MultiCoordset)
    True
    >>> prot = ag.protein
    >>> isinstance(prot, MultiCoordset)
    True"""

    def __init__(self, acsi):
        
        self._acsi = acsi
    
    def getACSIndex(self):
        """Return index of the active coordinate set."""
        
        return self._acsi
