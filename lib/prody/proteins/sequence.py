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

"""This module defines functions and classes for parsing, manipulating, and
analyzing multiple sequence alignments."""

__author__ = 'Anindita Dutta, Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Anindita Dutta, Ahmet Bakan'

__all__ = ['MSAFile', 'MSA', 'parseMSA', 
           'calcShannonEntropy', 'calcMutualInfo']

FASTA = 'fasta'
SELEX = 'selex'
STOCKHOLM = 'stockholm'

WSJOIN = ' '.join
ESJOIN = ''.join

import re
from os.path import isfile, splitext, split, getsize

from numpy import all, zeros, dtype, array, char

from prody import LOGGER
from prody.utilities import openFile

SPLITLABEL = re.compile('/*-*').split 

def splitLabel(label):
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


class MSAFile(object):
    
    """Yield tuples containing sequence id, sequence, start index, end index, 
    from an MSA file or object.
    
    >> msafile = fetchPfamMSA('piwi', alignment='seed')
    >> msa = MSAFile(msafile, filter=lambda tpl: 'ARATH' in tpl[0])
    >> for label, seq, rns, rne in msa: 
    ...     print label
    ... AGO6_ARATH
    ... AGO4_ARATH
    ... AGO10_ARATH"""
    
    def __init__(self, msa, aligned=True, filter=None, slice=None):
        """*msa* may be an MSA file in fasta, Stockholm, or Selex format or a
        Biopython MSA object.  If *aligned* is **True**, unaligned sequences 
        will cause an :exc:`IOError` exception.  *filter* is used for filtering
        sequences and must return a boolean for a single tuple argument that 
        contains ``(label, seq, start, end)``."""
        
        self._msa = None
        self._bio = None
        self._aligned = bool(aligned)
        self._lenseq = None
        self._numseq = None
        self._format = None
        self.setFilter(filter)
        self.setSlice(slice)
            
        if isfile(str(msa)):    
            self._msa = msa
            with openFile(msa) as msa: 
                line = msa.readline()
                if line[0] == '>':
                    LOGGER.info('Parsing MSA in fasta format')
                    self._iter = self._iterFasta
                    self._format = 'fasta'
                elif line[0] == '#' and 'STOCKHOLM' in line:
                    LOGGER.info('Parsing MSA in Stockholm format')
                    self._iter = self._iterStockholm
                    self._format = STOCKHOLM
                else:
                    LOGGER.info('Parsing MSA in Selex format')
                    self._iter = self._iterStockholm
                    self._format = SELEX
        else:    
            try:
                bio = iter(msa)
            except TypeError:
                raise TypeError('msa must be a filename or Bio object')
            else:
                record = bio.next()
                try:
                    label, seq = record.id, str(record.seq)
                except AttributeError:
                    raise TypeError('msa must be a filename or Bio object')
                else:
                    self._bio = iter(msa)
                    self._iter = self._iterBio
                    self._format = 'Biopython'
       
    def __del__(self):
        
        if self._msa:
            self._msa.close()
            
    def __iter__(self):
        
        filter = self._filter
        if filter is None:
            for label, seq in self._iter():
                label, start, end = splitLabel(label)
                yield label, seq, start, end
        else:
            for label, seq in self._iter():
                label, start, end = splitLabel(label)
                result = label, seq, start, end
                if filter(result):
                    yield result    
                             
    def _getFormat(self):
        """Return format of the MSA file."""
        
        return self._format
    
    format = property(_getFormat, None, doc="Format of the MSA file.")
            
    def _iterBio(self):
        """Yield sequences from a Biopython MSA object."""            
        
        aligned = self._aligned
        lenseq = self._lenseq
        slice = self._slice
        numseq = 0
        for record in self._bio:
            label, seq = record.id, str(record.seq)
            if not lenseq:
                self._lenseq = lenseq = len(seq)
            if aligned and lenseq != len(seq):
                raise IOError('sequence for {0:s} does not have '
                              'expected length {1:d}'
                              .format(label, lenseq))
            if slice:
                seq = seq[slice] 
            numseq += 1
            yield label, seq
        self._numseq = numseq

    def _iterFasta(self):
        """Yield sequences from an MSA file in fasta format."""

        aligned = self._aligned
        lenseq = self._lenseq
        slice = self._slice
        if slice is None:
            slice is False
        temp = []
        numseq = 0
        with openFile(self._msa) as msa: 
            label = msa.readline()[1:]
            for line in msa:
                if line[0] == '>':
                    seq = ESJOIN(temp)
                    if not lenseq:
                        self._lenseq = lenseq = len(seq)
                    if aligned and lenseq != len(seq):
                        raise IOError('sequence for {0:s} does not have '
                                      'expected length {1:d}'
                                      .format(label, lenseq))
                    if slice:
                        seq = seq[slice] 
                    numseq += 1
                    yield label, seq
                    temp = []
                    label = line[1:].strip()
                else:
                    temp.append(line.strip())
            numseq += 1
            yield label, ESJOIN(temp)
        self._numseq = numseq
    
    def _iterStockholm(self):
        """Yield sequences from an MSA file in Stockholm format."""


        aligned = self._aligned
        lenseq = self._lenseq
        slice = self._slice
        if slice is None:
            slice is False
        numseq = 0

        with openFile(self._msa) as msa:
            for line in msa:
                if line[0] == '#' or line[0] == '/':
                    continue
                items = line.split()
                label = WSJOIN(items[:-1])
                seq = items[-1]
                if not lenseq:
                    self._lenseq = lenseq = len(seq)
                if aligned and lenseq != len(seq):
                    raise IOError('sequence for {0:s} does not have '
                                  'expected length {1:d}'
                                  .format(label, lenseq))
                if slice:
                    seq = seq[slice] 
                numseq += 1
                yield label, seq
        self._numseq = numseq
        
    def numSequences(self):
        """Return number of sequences."""
        
        if self._numseq is None:
            for i in self._iter():
                pass
        return self._numseq
        
    def numResidues(self):
        """Return number of residues (or columns in the MSA)."""
        
        if self._lenseq is None:
            self._iter().next()
        return self._lenseq
    
    def getFilter(self):
        """Return function used for filtering sequences."""
        
        return self._filter
    
    def setFilter(self, filter):
        """Set function used for filtering sequences."""
        
        if filter is None:
            self._filter = None
            return
        
        if not callable(filter):
            raise TypeError('filter must be callable')
        
        try: 
            result = filter(('TEST_TITLE', 'SEQUENCE-WITH-GAPS', 1, 15))
        except Exception as err:
            raise TypeError('filter function must be not raise exceptions, '
                            'e.g. ' + str(err))
        else:
            try:
                result = result or not result
            except Exception as err: 
                raise ValueError('filter function must return a boolean, '
                                 'e.g. ' + str(err))
        self._filter = filter
    
    def getSlice(self):
        """Return object used to slice sequences."""
        
        return self._filter
    
    def setSlice(self, slice):
        """Set object used to slice sequences."""

        if slice is None:
            self._slice = None
            return
        
        seq = 'SEQUENCE'
        try: 
            result = seq[slice] 
        except Exception:
            raise TypeError('slice cannot be used for slicing sequences')
            arr = array(list(seq))
            try:
                result = arr[slice]
            except Exception:
                raise TypeError('slice cannot be used for slicing sequences')
            else:
                pass
                
        else:
            self._slice = slice
    
class MSA(object):
    
    """Store and manipulate multiple sequence alignments.
    
    >>> from prody import *
    >>> msa = parseMSA('piwi', alignment='seed')
    >>> msa
    <MSA: piwi (20 sequences, 404 residues)>

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
    <MSA: piwi' (2 sequences, 404 residues)>
    
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
    <MSA: piwi' (10 sequences, 20 residues)>
    
    *Refinement*
    
    Columns in an MSA that are gaps for a given sequence can be eliminated
    from the data as follows:
        
    >>> msa[:, 'YQ53_CAEEL'] # doctest: +ELLIPSIS
    <MSA: piwi' (20 sequences, 328 residues)>
    
    This operation removed 76 columns, which is the number of gaps in sequence
    with label ``'YQ53_CAEEL'``."""
    
    def __init__(self, msa, **kwargs):
        """*msa* may be an :class:`MSAFile` instance or an MSA file in a 
        supported format."""
        
        try:
            ndim, dtype_, shape = msa.ndim, msa.dtype, msa.shape
        except AttributeError:
            try:
                numseq, lenseq = msa.numSequences(), msa.numResidues()
            except AttributeError:
                try:
                    msa = MSAFile(msa)
                except Exception as err:
                    raise TypeError('msa was not recognized ({0:s})'
                                    .format(str(err)))
                else:
                    numseq, lenseq = msa.numSequences(), msa.numResidues()
            
            self._msa = msaarray = zeros((numseq, lenseq), dtype='|S1')
            self._labels = []
            labels = self._labels.append
            self._mapping = mapping = {}
            
            for i, (label, seq, start, end) in enumerate(msa):
                labels((label, start, end))
                mapping[label] = i
                msaarray[i] = list(seq)
        else:
            if ndim != 2:
                raise ValueError('msa.dim must be 2')
            if dtype_ != dtype('|S1'):
                raise ValueError('msa must be a character array')
            numseq = shape[0]
            self._labels = labels = kwargs.get('labels')
            if labels and len(self._labels) != numseq:
                raise ValueError('len(labels) must be equal to number of '
                                 'sequences')
            
            self._mapping = mapping = kwargs.get('mapping')
            if mapping is None and labels is not None:
                # map labels to sequence index
                self._mapping = mapping = {
                    splitLabel(label)[0]: i for i, label in enumerate(labels)
                }
                
            if labels is None:
                self._labels = [None] * numseq
                
            
            self._msa = msa
        self._title = kwargs.get('title', 'Unknown')
        
        
    def __str__(self):
        
        return 'MSA ' + self._title
        
    def __repr__(self):
        
        return '<MSA: {0:s} ({1:d} sequences, {2:d} residues)>'.format(
                self._title, self.numSequences(), self.numResidues())
    
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
            if ndim == 0:
                return msa
            elif ndim == 1:
                if cols is None:
                    label, start, end = splitLabel(self._labels[rows])
                    return label, msa.tostring(), start, end
                else:
                    return msa.tostring()
            else:
                try:
                    labels = self._labels[rows]
                except TypeError:
                    temp = self._labels
                    labels = [temp[i] for i in rows]
                return MSA(msa, title=self._title + '\'', labels=labels) 
               
    def __iter__(self):
        
        for i, label in enumerate(self._labels):
            label, start, end = splitLabel(label)
            yield label, self._msa[i].tostring(), start, end
    
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
    
    def numSequences(self):
        """Return number of sequences."""
        
        return self._msa.shape[0]

    def numResidues(self):
        """Return number of residues (or columns in the MSA)."""
        
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
            return splitLabel(self._labels[index])[0]
                
    def getResnums(self, index):
        """Return starting and ending residue numbers (:term:`resnum`) for the
        sequence at given *index*."""

        index = self._mapping.get(index, index)
        return splitLabel(self._labels[index])[1:]
    
    def getArray(self):
        """Return a copy of the MSA character array."""
        
        return self._msa.copy()
    
    def _getArray(self):
        """Return MSA character array."""
        
        return self._msa

    
def parseMSA(msa, **kwargs):
    """Return an :class:`MSA` instance that stores multiple sequence alignment
    and sequence labels parsed from Stockholm, Selex, or Fasta format *msa* 
    file.  If *msa* is a Pfam id code or accession, MSA file will be downloaded
    using :func:`fetchPfamMSA` with default parameters.  Note that *msa* may be
    a compressed file. Uncompressed MSA files are parsed using C code at a 
    fraction of the time it would take to parse compressed files in Python."""
    
    msa = str(msa)
    if isfile(msa):
        # if MSA is a compressed file or filter/slice is passed, use 
        #   Python parsers
        ext = splitext(msa)[1] 
        if ext == '.gz' or 'filter' in kwargs or 'slice' in kwargs:
            return MSA(msa, **kwargs)
        else:
            filename = msa
    else:    
        from prody.database import fetchPfamMSA
        try:
            filename = fetchPfamMSA(msa, **kwargs)
        except IOError:
            raise ValueError('msa must be an MSA filename or a Pfam accession')
        else:
            if 'compressed' in kwargs:
                return MSA(filename, **kwargs)

    msafile = MSAFile(filename)
    title = splitext(split(msa)[1])[0]
    format = msafile.format
    lenseq = msafile.numResidues()
    numseq = getsize(filename) / (lenseq + 10)
    del msafile
    
    if format == FASTA:
        from .msatools import parseFasta
        msaarr, labels, mapping = parseFasta(filename, lenseq, numseq)
    elif format == SELEX or format == STOCKHOLM:
        from .msatools import parseSelex
        msaarr, labels, mapping = parseSelex(filename, lenseq, numseq)
    else:
        raise IOError('MSA file format is not recognized')
    return MSA(msa=msaarr, title=title, labels=labels, mapping=mapping)

def calcShannonEntropy(msa, ambiguity=True, omitgaps=False):
    """Return Shannon entropy array calculated for *msa*, which may be 
    an :class:`MSA` instance or a 2D Numpy character array.  Implementation 
    is case insensitive and handles ambiguous amino acids as follows:
    
      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      
    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.
    
    All non-alphabet characters are considered as gaps, and they are handled 
    in two ways:  
      
      * as a distinct character with its own probability, by default 
      * non-existent, the probability of observing amino acids in a given
        column is adjusted, when *omitgaps* is **True**."""
    
    try:
        msa = msa._getArray()
    except AttributeError:
        pass
    
    try:
        dtype_, ndim, shape = msa.dtype, msa.ndim, msa.shape
    except AttributeError:
        raise TypeError('msa must be an MSA instance or a 2D character array')
        
    if dtype_ != dtype('|S1') or ndim != 2:
        raise TypeError('msa must be an MSA instance or a 2D character array')
        
    entropy = zeros(shape[1], float)
    from .msatools import calcShannonEntropy
    calcShannonEntropy(msa, entropy, ambiguity=bool(ambiguity),
                       omitgaps=bool(omitgaps))
    return entropy
    
    
def calcMutualInfo(msa, ambiguity=True, turbo=True, **kwargs):
    """Return mutual information matrix calculated for *msa*, which may be an 
    :class:`MSA` instance or a 2D Numpy character array.  Implementation 
    is case insensitive and handles ambiguous amino acids as follows:
    
      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      * Joint probability of observing a pair of ambiguous amino acids is 
        allocated to all potential combinations, e.g. probability of **XX**
        is allocated to 400 combinations of standard amino acids, similarly
        probability of **XB** is allocated to 40 combinations of *D* and *N*
        with the standard amino acids.      
      
    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.  All non-alphabet characters 
    are considered as gaps.
    
    By default, the will try to run in the *turbo* mode, which uses memory
    as large as the MSA array itself but runs four to five times faster.  If
    memory allocation fails, the implementation will switch to slower and
    memory efficient mode."""
    
    try:
        msa = msa._getArray()
    except AttributeError:
        pass
    
    try:
        dtype_, ndim, shape = msa.dtype, msa.ndim, msa.shape
    except AttributeError:
        raise TypeError('msa must be an MSA instance or a 2D character array')
        
    if dtype_ != dtype('|S1') or ndim != 2:
        raise TypeError('msa must be an MSA instance or a 2D character array')
        
    mutinfo = zeros((shape[1], shape[1]), float)
    from .msatools import calcMutualInfo
    turbo = bool(turbo)
    LOGGER.timeit('_mutinfo')
    turbo = calcMutualInfo(msa, mutinfo, ambiguity=bool(ambiguity), 
                           turbo=turbo, debug=kwargs.get('debug', False))
    if turbo:
        LOGGER.report('Mutual information matrix was calculated in turbo mode '
                      'in %.2fs.', '_mutinfo')
    else:
        LOGGER.report('Mutual information matrix was calculated '
                      'in %.2fs.', '_mutinfo')   
    return mutinfo
