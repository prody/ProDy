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

"""This module handles individual sequences."""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import re

from numpy import arange, char, fromstring, zeros

from prody import LOGGER, PY3K

try:
    range = xrange
except NameError:
    pass

SPLITLABEL = re.compile('/*-*').split

__all__ = ['Sequence']

def splitSeqLabel(label):
   """Return label, starting residue number, and ending residue number parsed
   from sequence label."""

   try:
       idcode, start, end = SPLITLABEL(label)
   except Exception:
       return label, None, None
   else:   
       try:
           return idcode, int(start), int(end)
       except Exception:
           return label, None, None


class Sequence(object):
    
    """Handle individual sequences of an `.MSA` object""" 
    
    __slots__ = ['_msa', '_seq', '_index', '_label']
    
    def __init__(self, *args):
        """Depending on input arguments, instances may point to an 
        :class:`.MSA` object or store its own data:
        
        *MSA Pointer*

        An :class:`.MSA` instance and an index:

        >>> from prody import *
        >>> msa = parseMSA('piwi_seed.sth')
        >>> Sequence(msa, 0)
        <Sequence: YQ53_CAEEL (piwi_seed[0]; length 404; 328 residues and 76 gaps)>
        >>> msa[0]
        <Sequence: YQ53_CAEEL (piwi_seed[0]; length 404; 328 residues and 76 gaps)>
        
        *Independent*
        
        Instantiation with sequence and label (optional) string:
        
        >>> Sequence('SOME-SEQUENCE-STRING', 'MySeq/1-18')
        <Sequence: MySeq (length 20; 18 residues and 2 gaps)>"""
        
        if len(args) == 2:
            one, two = args
            try:
                one.lower, two.lower
            except AttributeError:
                self._msa = one
                self._index = two
                self._seq = self._label = None
            else:
                self._seq = fromstring(one, '|S1')
                self._label = two
                self._msa = self._index = None
        elif len(args) == 1:
            self._seq = fromstring(args[0], '|S1')
            self._msa = self._index = None
            self._label = ''
        else:
            raise ValueError('msa and index, or seq [and label] must be'
                             'specified')
    
    @property
    def _array(self):
        """Sequence data array."""
        
        return self._seq if self._msa is None else self._msa._msa[self._index]
        
    def __str__(self):
        
        if PY3K:
            return self._array.tostring().decode()
        else:
            return self._array.tostring()
        
    def __len__(self):
        
        return len(self._array)
                      
    def __repr__(self):
        
        msa = ''
        if self._msa is not None:
            msa = '{0}[{1}]; '.format(self._msa.getTitle(), self._index)
        return ('<Sequence: {0} ({1}length {2}; {3} residues and '
                '{4} gaps)>').format(self.getLabel(), msa, len(self), 
                self.numResidues(), self.numGaps())
    
    def __eq__(self, other):
        
        try:
            this = self._array
            that = other._array
            return this.shape == that.shape and (this == that).all()
        except AttributeError:
            return False 
    
    def getMSA(self):
        """Return :class:`.MSA` instance or **None**."""
        
        return self._msa
    
    def getIndex(self):
        """Return sequence index or **None**."""
        
        return self._index
    
    # This function should be able to update MSA._mapping and MSA._labels
    #def setLabel(self, label):
    #    """Set the label to be associated with object"""
    #    
    #    self._label = str(label)
        
    def getLabel(self, full=False):
        """Return label of the sequence."""
        
        label = self._label
        if label is None: label = self._msa._labels[self._index] 
        return label if full else splitSeqLabel(label)[0]
    
    def numGaps(self):
        """Return number of gap characters."""
        
        array = self._array
        return len(array) - sum(char.isalpha(array))
            
    def numResidues(self):
        """Return the number of alphabet characters."""
        
        return sum(char.isalpha(self._array))
        
    def getResnums(self, gaps=False):
        """Return list of residue numbers associated with non-gapped *seq*.  
        When *gaps* is **True**, return a list containing the residue numbers 
        with gaps appearing as **None**.  Residue numbers are inferred from the
        full label.  When label does not contain residue number information, 
        indices a range of numbers starting from 1 is returned."""
        
        title, start, end = splitSeqLabel(self.getLabel(True))
        try:
            start, end = int(start), int(end)
        except:
            LOGGER.info('Cannot parse label start, end values, Setting '
                        'resnums 1 to {0:d}'.format(self.numResidues()))
            start, end = 1, self.numResidues()
        else:
            if (end - start + 1) != self.numResidues():
                LOGGER.info('Label start-end position does not match '
                            'length of ungapped sequence. Setting '
                            'resnums 1 to {0:d}'.format(self.numResidues()))
                start, end = 1, self.numResidues()
        
        resnums = iter(range(start, end + 1))
        if gaps:
            return [next(resnums) if torf else None 
                    for torf in char.isalpha(self._array())]
        else:
            return list(resnums)
        
    def copy(self):
        """Return a copy of the instance that owns its sequence data."""
        
        return Sequence(str(self), self.getLabel())
