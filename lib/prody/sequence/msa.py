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

"""This module defines MSA analysis functions."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import all, zeros, dtype, array, char, fromstring, arange

from .msafile import splitSeqLabel, MSAFile 

from prody import LOGGER

__all__ = ['MSA', 'refineMSA']

class MSA(object):
    
    """Store and manipulate multiple sequence alignments.
    
    >>> from prody import *
    >>> fetchPfamMSA('piwi', alignment='seed') # DOCTEST: +SKIP
    'piwi_seed.sth'
    >>> msafile = 'piwi_seed.sth'
    >>> msa = parseMSA(msafile)
    >>> msa
    <MSA: piwi_seed (20 sequences, 404 residues)>

    *Querying*
    
    You can query whether a sequence in contained in the instance using
    the UniProt identifier of the sequence as follows:
        
    >>> 'YQ53_CAEEL' in msa
    True
    
    *Indexing and slicing*
    
    Retrieve a sequence at a given index:
    
    >>> msa[0] # doctest: +ELLIPSIS
    ('YQ53_CAEEL', 'DIL...YK', 650, 977)
    
    Retrieve a sequence by UniProt ID:
    
    >>> msa['YQ53_CAEEL'] # doctest: +ELLIPSIS
    ('YQ53_CAEEL', 'DIL...YK', 650, 977)
    
    Slice an MSA instance:
    
    >>> msa[:2]
    <MSA: piwi_seed' (2 sequences, 404 residues)>
    
    Slice using a list of UniProt IDs:
    
    >>> msa[:2] == msa[['YQ53_CAEEL', 'Q21691_CAEEL']]
    True
    
    Retrieve a character or a slice of a sequence:

    >>> msa[0,0]
    'D'
    >>> msa[0,0:10]
    'DILVGIAR.E'
    
    Slice MSA rows and columns:
    
    >>> msa[:10,20:40]
    <MSA: piwi_seed' (10 sequences, 20 residues)>
    
    *Refinement*
    
    Columns in an MSA that are gaps for a given sequence can be eliminated
    from the data as follows:
        
    >>> msa[:, 'YQ53_CAEEL'] # doctest: +ELLIPSIS
    <MSA: piwi_seed' (20 sequences, 328 residues)>
    
    This operation removed 76 columns, which is the number of gaps in sequence
    with label ``'YQ53_CAEEL'``.
    
    *Selective parsing*
    
    Filtering and slicing available to :class:`.MSAFile` class can be used to 
    parse an MSA selectively, which may be useful in low memory situations:
        
    >>> msa = parseMSA(msafile, filter=lambda lbl, seq: 'ARATH' in lbl, 
    ...                slice=list(range(10)) + list(range(394,404)))
    >>> msa
    <MSA: piwi_seed (3 sequences, 20 residues)>

    Compare this to result from parsing the complete file:
    >>> parseMSA(msafile)
    <MSA: piwi_seed (20 sequences, 404 residues)>"""
    
    def __init__(self, msa, title='Unknown', labels=None, **kwargs):
        """*msa* must be a 2D Numpy character array. *labels* is a list of
        sequence labels (or titles).  *mapping* should map label or part of 
        label to sequence index in *msa* array. If *mapping* is not given,
        one will be build from *labels*."""
        
        try:
            ndim, dtype_, shape = msa.ndim, msa.dtype, msa.shape
        except AttributeError:
            raise TypeError('msa is not a Numpy array')

        self._aligned = aligned = kwargs.get('aligned', True)
        if aligned:
            if ndim != 2:
                raise ValueError('msa.dim must be 2')
            if dtype_ != dtype('|S1'):
                raise ValueError('msa must be a character array')
        numseq = shape[0]

        if labels and len(labels) != numseq:
            raise ValueError('len(labels) must be equal to number of '
                             'sequences')
        self._labels = labels
        mapping = kwargs.get('mapping')
        if mapping is None:
            if labels is not None:
                # map labels to sequence index
                self._mapping = mapping = {
                splitSeqLabel(label)[0]: i for i, label in enumerate(labels)
                }
        elif mapping:
            try:
                mapping['isdict']
            except KeyError:
                pass
            except Exception:
                raise TypeError('mapping must be a dictionary')
        self._mapping = mapping                
        if labels is None:
            self._labels = [None] * numseq
            
        self._msa = msa
        self._title = str(title) or 'Unknown'
        self._split = bool(kwargs.get('split', True))

        
    def __str__(self):
        
        return 'MSA ' + self._title
        
    def __repr__(self):
        
        if self._aligned:
            return '<MSA: {0} ({1} sequences, {2} residues)>'.format(
                    self._title, self.numSequences(), self.numResidues())
        else:
            return '<MSA: {0} ({1} sequences, not aligned)>'.format(
                    self._title, self.numSequences())
    
    def __getitem__(self, index):
        
        try:
            length = len(index)
        except TypeError: # type(index) -> int, slice
            rows, cols = index, None
        else:
            try:
                _ = index.strip
            except AttributeError: # type(index) -> tuple, list
                try:
                    _ = index.sort
                except AttributeError: # type(index) -> tuple
                    if length == 1:
                        rows, cols = index[0], None
                    elif length == 2:
                        rows, cols = index
                    else:
                        raise IndexError('invalid index: ' + repr(index))
                else: # type(index) -> list
                    rows, cols = index, None
            else: # type(index) -> str
                rows, cols = index, None 

        try: # ('PROT_HUMAN', )
            rows = self._mapping.get(rows, rows)
        except (TypeError, KeyError):
            mapping = self._mapping
            try:
                rows = [mapping[key] for key in rows]
            except (KeyError, TypeError):
                pass

        if cols is None:
            msa = self._msa[rows]
        else:
            if not self._aligned:
                raise ValueError('msa is not aligned, '
                                 'column indexing is not possible')
            try:
                cols = self._mapping[cols]
            except (KeyError, TypeError):
                pass
            else:
                cols = char.isalpha(self._msa[cols])
                
            try:
                msa = self._msa[rows, cols]
            except Exception:
                raise IndexError('invalid index: ' + str(index))
            
        try:
            shape, ndim = msa.shape, msa.ndim
        except AttributeError:
            return msa
        else:
            if ndim < 2:
                seq = msa.tostring()
                try:
                    seq.format
                except AttributeError:
                    seq = seq.decode('utf')
                if cols is None:
                    label, start, end = splitSeqLabel(self._labels[rows])
                    return label, seq, start, end
                else:
                    return seq
            else:
                if msa.base is not None:
                    msa = msa.copy()
                try:
                    labels = self._labels[rows]
                except TypeError:
                    temp = self._labels
                    labels = [temp[i] for i in rows]
                return MSA(msa, title=self._title + '\'', labels=labels) 
               
    def __iter__(self):
        
        if self._split:
            for i, label in enumerate(self._labels):
                label, start, end = splitSeqLabel(label)
                seq = self._msa[i].tostring()
                try:
                    seq.format
                except AttributeError:
                    seq = seq.decode('utf-8')
                yield label, seq, start, end
        else:
            for i, label in enumerate(self._labels):
                yield label, self._msa[i].tostring()
            
    
    def __contains__(self, key):
        
        try:
            return key in self._mapping
        except Exception:
            pass
        return False
    
    def __eq__(self, other):

        try:
            other = other._getArray()
        except AttributeError:
            return False
        
        try:
            return all(other == self._msa)
        except Exception:
            pass
        return False
        
    def _getSplit(self):
        
        return self._split
    
    def _setSplit(self, split):
        
        self._split = bool(split)
        
    split = property(_getSplit, _setSplit, 
                     doc='Return split label when iterating or indexing.')

    def isAligned(self):
        """Return **True** if MSA is aligned."""
        
        return self._aligned
        
    def numSequences(self):
        """Return number of sequences."""
        
        return self._msa.shape[0]

    def numResidues(self):
        """Return number of residues (or columns in the MSA), if MSA is 
        aligned."""
        
        if self._aligned:
            return self._msa.shape[1]
    
    def getTitle(self):
        """Return title of the instance."""
        
        return self._title
    
    def setTitle(self, title):
        """Set title of the instance."""
        
        self._title = str(title)
        
    def getLabel(self, index, full=False):
        """Return label of the sequence at given *index*.  Residue numbers will
        be removed from the sequence label, unless *full* is **True**."""
        
        index = self._mapping.get(index, index)
        if full:
            return self._labels[index]
        else:
            return splitSeqLabel(self._labels[index])[0]
                
    def getResnums(self, index):
        """Return starting and ending residue numbers (:term:`resnum`) for the
        sequence at given *index*."""

        index = self._mapping.get(index, index)
        return splitSeqLabel(self._labels[index])[1:]
    
    def getArray(self):
        """Return a copy of the MSA character array."""
        
        return self._msa.copy()
    
    def _getArray(self):
        """Return MSA character array."""
        
        return self._msa
    
    def getIndex(self, label):
        """Return index of the sequence entry with *label*."""
        
        try:
            index = self._mapping.get(label)
        except TypeError:
            raise TypeError('label must be a string')
        else:
            return index


def refineMSA(msa, label=None, row_occ=None, col_occ=None, **kwargs):
    """Refine *msa* by removing sequences (rows) and residues (columns) that 
    contain gaps.
    
    :arg msa: multiple sequence alignment
    :type msa: :class:`.MSA`
    
    :arg label: sequence label used for indexing
    :type label: str
    
    :arg row_occ: row occupancy, sequences with less occupancy will be 
        removed after *label* refinement is applied
    :type row_occ: float
    
    :arg col_occ: column occupancy, residue positions with less occupancy
        will be removed after other refinements are applied
    :type col_occ: float
    
    For Pfam MSA data, *label* is UniProt entry name for the protein.  You may
    also use PDB structure and chain identifiers, e.g. ``'1p38'`` or 
    ``'1p38A'``, for *label* argument and UniProt entry names will be parsed 
    using :func:`.parsePDBHeader` function (see also :class:`.Polymer` and 
    :class:`DBRef`)."""

    # if msa is a char array, it will be refined but label won't work
    try:    
        ndim, dtype_ = msa.ndim, msa.dtype
    except AttributeError:
        try:
            arr = msa._getArray()
        except AttributeError:
            raise TypeError('msa must be a character array or an MSA instance') 
        ndim, dtype_ = arr.ndim, arr.dtype
    else:
        arr, msa = msa, None
    
    if dtype('|S1') != dtype_:
        raise ValueError('msa must be a character array or an MSA instance')
    if ndim != 2:
        raise ValueError('msa must be a 2D array or an MSA instance')
            

    title = []
    cols = None
    if label is not None:
        try:
            upper, lower = label.upper(), label.lower()
        except AttributeError:
            raise TypeError('label must be a string')

        if msa is None:
            raise TypeError('msa must be an MSA instance, '
                            'label cannot be used')
        
        index = msa.getIndex(label)
        if index is None: index = msa.getIndex(upper)            
        if index is None: index = msa.getIndex(lower)
        
        if index is None and (len(label) == 4 or len(label) == 5):
            from prody import parsePDBHeader
            try:
                polymers = parsePDBHeader(label[:4], 'polymers')
            except Exception as err:
                LOGGER.warn('failed to parse header for {0} ({1})'
                            .format(label[:4], str(err)))
            else:
                chid = label[4:].upper()
                for poly in polymers:
                    if chid and poly.chid != chid:
                        continue
                    for dbref in poly.dbrefs:
                        if index is None: 
                            index = msa.getIndex(dbref.idcode)
                            if index is not None:
                                LOGGER.info('{0} idcode {1} for {2}{3}'
                                            'is found in {3}'.format(
                                            dbref.database, dbref.idcode,
                                            label[:4], poly.chid, str(msa)))
                        if index is None: 
                            index = msa.getIndex(dbref.accession)
                            if index is not None:
                                LOGGER.info('{0} idcode {1} for {2}{3}'
                                            ' is found in {3}'.format(
                                            dbref.database, dbref.accession,
                                            label[:4], poly.chid, str(msa)))
                    
            
        if index is None:        
            raise ValueError('label is not in msa, or msa is not indexed')
        title.append('label=' + label)
        cols = char.isalpha(arr[index]).nonzero()[0]
        arr = arr.take(cols, 1)
        
    rows = None
    from .analysis import calcMSAOccupancy
    if row_occ is not None:
        try:
            row_occ = float(row_occ)
        except Exception as err:
            raise TypeError('row_occ must be a float ({0})'.format(str(err)))
        assert 0. <= row_occ <= 1., 'row_occ must be between 0 and 1'
        
        rows = (calcMSAOccupancy(arr, 'row') >= row_occ).nonzero()[0]
        arr = arr[rows]
        title.append('row_occ>=' + str(row_occ))

    if col_occ is not None:
        try:
            col_occ = float(col_occ)
        except Exception as err:
            raise TypeError('col_occ must be a float ({0})'.format(str(err)))
        assert 0. <= col_occ <= 1., 'col_occ must be between 0 and 1'
        
        cols = (calcMSAOccupancy(arr, 'col') >= col_occ).nonzero()[0]
        arr = arr.take(cols, 1)
        title.append('col_occ>=' + str(col_occ))
        
    if not title:
        raise ValueError('label, row_occ, col_occ all cannot be None')
    
    # depending on slicing of rows, arr may not have it's own memory
    if arr.base is not None:
        arr = arr.copy()
    
    if msa is None:
        return arr
    else:
        if rows is None:
            from copy import copy
            labels = copy(msa._labels)
            mapping = copy(msa._mapping)
        else:
            labels = msa._labels
            labels = [labels[i] for i in rows]
            mapping = None
        return MSA(arr, title=msa.getTitle() + ' refined ({0})'
                   .format(', '.join(title)), labels=labels, mapping=mapping)
