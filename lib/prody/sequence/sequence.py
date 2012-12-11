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

""" This module handles individual sequences and label of an `.MSA` object."""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import re
SPLITLABEL = re.compile('/*-*').split
from prody import LOGGER
import numpy as np

__all__ = ['Sequence', 'splitSeqLabel']

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
    
    def __init__(self, **kwargs):
        """ Takes the following input arguements:
        arg msa: a :class:`.MSA` ojbect for reference.
        type msa: class:`.MSA` object
        
        arg index: can specify index for a particular sequence of the *msa*
        type index: int
        
        arg label: can specify a label. If *msa* is given than *seq* is set
                    based on *label*
        type label: str
        
        arg seq: user specified sequence
        type seq: str
        
        *Instantiation*
        >>> seq = Sequence(msa=msaobj, index=0)
        >>> seq = Sequence(seq='SOME-SEQUENCE-STRING', label='MySeq/1-18')"""
        
        msa = kwargs.get('msa')
        label = kwargs.get('label')
        seq = kwargs.get('seq')
        index = kwargs.get('index')
        self._msa = None
        self._index = None
        self._label = None
        self._seq = None
        if msa is not None:
            try:
                self._msasize = msa.numSequences()
            except AttributeError:
                raise TypeError('msa must be an MSA object')
            else:
                self._msa = msa
                if index is None:
                    raise ValueError('Index missing. {0} should be '
                                     'accompanied by index'.format(str(msa)))
                else:
                    index = int(index)
                    if index >= self._msasize or index < 0:
                        raise ValueError('index cannot be < than 0 or > than '
                                         'the number of sequences in the msa')
                    self._index = index
        else:
            if label is not None and seq is not None:
                self._label = str(label)
                self._seq = str(seq)
            else:
                raise ValueError('Either msa and index, or label and sequence'
                                 ' need to be specified.')
        
    def __str__(self):
        
        if self._seq is None:
            return self._msa._getArray()[self._index].tostring()
        else:
            return self._seq
     
    def __len__(self):
        
        if self._seq is None:
            return self._msa.numResidues()
        else:        
            return len(self._seq)
                      
    def __repr__(self):
        
        return '<Sequence: {0} (length {1}; {2} residues and {3} gaps)>'.format(
                self.getLabel(), len(self), self.numResidues(), self.numGaps())
    
    def setMSA(self, msa):
        """ Set a reference to an `.MSA` object"""
        
        try:
            self._msaSize = msa.numSequences()
        except AttributeError:
            raise TypeError('msa must be an MSA object')
        else:
            self._msa = msa
    
    def getMSA(self):
        """ Get the msa reference"""
        
        return self._msa
    
    def getIndex(self):
        """ Get the index associated with object"""
        
        return self._index
        
    def setLabel(self, label):
        """Set the label to be associated with object"""
        
        self._label = str(label)
        
    def getLabel(self):
        """Get the label associated with object"""
        
        if self._label is None:
            return self._msa._labels[self._index]
        else:
            return self._label
    
    def setSeq(self, seq):
        """ Set the *seq* to be associated with object"""
        
        self._seq = str(seq)
        
    def numGaps(self):
        """Return the no of gaps associated with *seq*. If *seq* is not set,
        returns *None*"""
        
        return len(self) - sum(np.char.isalpha(list(str(self))))
            
            
    def numResidues(self):
        """Return the length of nongapped *seq*. If *seq* is not set,
        returns *None*"""
        
        return sum(np.char.isalpha(list(str(self))))
        
    def getResnums(self, gaps=False):
        """Return a list of residue numbers associated with nongapped *seq* if
        *gaps* is set to False (default). Else return a list containing the
        residue numbers with gaps appearing as *None*. Resiues numbers are
        based on *label* start-end position if given, otherwise the residues
        are indices starting at 1 to length of ungapped *seq*"""
        
        title, start, end = splitSeqLabel(self.getLabel())
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
        resnums = list(xrange(start, end + 1))
        if gaps:
            resnums = np.zeros(len(self), dtype='object')
            cols = np.char.isalpha(list(str(self)))
            resnums[cols] = list(xrange(start, end + 1))
            resnums[np.where(cols == False)[0]] = None
            resnums = list(resnums)
        return resnums
        
    def copy(self):
        
        return Sequence(label=getLabel(), seq=str(self))
        
    
        
        
        
            