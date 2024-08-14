# -*- coding: utf-8 -*-
"""This module handles individual sequences."""

from numpy import char, fromstring

from prody import LOGGER, PY3K
from prody.atomic import Atomic

from prody.utilities import splitSeqLabel

__all__ = ['Sequence']

class Sequence(object):

    """Handle individual sequences of an :class:`.MSA` object"""

    __slots__ = ['_msa', '_seq', '_index', '_label']

    def __init__(self, *args):
        """Depending on input arguments, instances may point to an
        :class:`.MSA` object or store its own data."""

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
            return self._array.tobytes().decode()
        else:
            return self._array.tobytes()

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
        """Returns :class:`.MSA` instance or **None**."""

        return self._msa

    def getArray(self):
        return self._array

    def getIndex(self):
        """Returns sequence index or **None**."""

        return self._index

    # This function should be able to update MSA._mapping and MSA._labels
    #def setLabel(self, label):
    #    """Set the label to be associated with object"""
    #
    #    self._label = str(label)

    def getLabel(self, full=False):
        """Returns label of the sequence."""

        label = self._label
        if label is None:
            label = self._msa._labels[self._index]
        return (label if full else splitSeqLabel(label)[0]).strip()

    def numGaps(self):
        """Returns number of gap characters."""

        array = self._array
        return len(array) - sum(char.isalpha(array))

    def numResidues(self):
        """Returns the number of alphabet characters."""

        return sum(char.isalpha(self._array))

    def getResnums(self, gaps=False, report_match=False):
        """Returns list of residue numbers associated with non-gapped *seq*.
        When *gaps* is **True**, return a list containing the residue numbers
        with gaps appearing as **None**.  
        
        Residue numbers are inferred from the full label if possible. 
        When the label does not contain residue number information, 
        a range of numbers starting from 1 is returned."""

        title, start, end = splitSeqLabel(self.getLabel(True))
        match = False
        try:
            start, end = int(start), int(end)
        except:
            LOGGER.info('Cannot parse start and end values from sequence label {0}. Setting '
                        'resnums 1 to {1:d}'.format(title, self.numResidues()))
            start, end = 1, self.numResidues()
        else:
            if (end - start + 1) != self.numResidues():
                LOGGER.info('Label {0} start-end entry does not match '
                            'length of ungapped sequence. Setting '
                            'resnums 1 to {1:d}'.format(title, self.numResidues()))
                start, end = 1, self.numResidues()
            else:
                LOGGER.info('Label {0} start-end entry matches '
                            'length of ungapped sequence. Setting '
                            'resnums {1:d} to {2:d}'.format(title, start, end))             
                match = True

        resnums = iter(range(start, end + 1))
        if gaps:
            result = [next(resnums) if torf else None
                      for torf in char.isalpha(self._array)]
        else:
            result = list(resnums)

        if report_match:
            return match, result
        return result

    def copy(self):
        """Returns a copy of the instance that owns its sequence data."""

        return Sequence(str(self), self.getLabel())
