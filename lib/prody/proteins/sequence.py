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
           'calcShannonEntropy', 'calcMutualInfo',
           'calcMSAOccupancy']

FASTA = 'fasta'
SELEX = 'selex'
STOCKHOLM = 'stockholm'

WSJOIN = ' '.join
ESJOIN = ''.join

NUMLINES = 1000

import re
from os.path import isfile, splitext, split, getsize

from numpy import all, zeros, dtype, array, char, fromstring

from prody import LOGGER, PY2K
from prody.utilities import openFile

if PY2K: range = xrange

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
    
    """Handle MSA files in FASTA, SELEX and Stockholm formats. 
    
    >>> from prody import * 
    >>> msafile = fetchPfamMSA('piwi', alignment='seed')

    *Reading a file*

    Iterating over a file will yield sequence id, sequence, residue start and 
    end indices:

    >>> msa = MSAFile(msafile)
    >>> for seq in msa: # doctest: +ELLIPSIS 
    ...     print seq
    ('YQ53_CAEEL', 'DILVGIAR.EKKP...NLAKRGRNNYK', 650, 977)
    ('Q21691_CAEEL', 'TIVFGIIA.EKRP...NLAKRGHNNYK', 673, 1001)
    ('AGO6_ARATH', 'FILCILPERKTSD...LAAAQVAQFTK', 541, 851)
    (...)
    ('O02095_CAEEL', 'QLLFFVVK..SRY...RYSQRGAMVLA', 574, 878)
    ('Q19645_CAEEL', 'PFVLFISD..DVP...ELAKRGTGLYK', 674, 996)
    ('O62275_CAEEL', 'TFVFIITD.DSIT...EYAKRGRNLWN', 594, 924)
    
    *Filtering sequences*
    
    Any function that takes label and sequence arguments and returns a boolean
    value can be used for filtering the sequences.  A sequence will be yielded 
    if the function returns **True**.  In the following example, sequences from
    organism *ARATH* are filtered:
    
    >>> msa = MSAFile(msafile, filter=lambda lbl, seq: 'ARATH' in lbl)
    >>> for seq in msa: # doctest: +ELLIPSIS 
    ...     print seq
    ('AGO6_ARATH', 'FIL...FTK', 541, 851)
    ('AGO4_ARATH', 'FIL...FMK', 577, 885)
    ('AGO10_ARATH', 'LLL...YLE', 625, 946)

    *Slicing sequences*
    
    A list of integers can be used to slice sequences as follows.
    
    >>> msa = MSAFile(msafile, slice=list(range(10)) + list(range(394,404)))
    >>> for seq in msa: # doctest: +ELLIPSIS 
    ...     print seq
    ('YQ53_CAEEL', 'DILVGIAR.ELAKRGRNNYK', 650, 977)
    ('Q21691_CAEEL', 'TIVFGIIA.ELAKRGHNNYK', 673, 1001)
    ('AGO6_ARATH', 'FILCILPERKAAAQVAQFTK', 541, 851)
    (...)
    ('O02095_CAEEL', 'QLLFFVVK..YSQRGAMVLA', 574, 878)
    ('Q19645_CAEEL', 'PFVLFISD..LAKRGTGLYK', 674, 996)
    ('O62275_CAEEL', 'TFVFIITD.DYAKRGRNLWN', 594, 924)"""
    
    def __init__(self, msa, mode='r', format=None, aligned=True, **kwargs):
        """*msa* may be a filename or a stream.  Multiple sequence alignment
        can be read or written in FASTA (:file:`.fasta`), Stockholm 
        (:file:`.sth`), or SELEX (:file:`.slx`) *format*.  For recognized 
        extensions, specifying *format* is not needed.  If *aligned* is 
        **True**, unaligned sequences in the file or stream will cause an 
        :exc:`IOError` exception.  *filter*, a function that returns a 
        boolean, can be used for filtering sequences.  *slice* can be used 
        to slice sequences, and is applied after filtering."""
        
        if mode not in 'rwa':
            raise ValueError("mode string must be one of 'r', 'w', or 'a', "
                             "not {0}".format(repr(mode)))
        
        self._bio = self._filename = None
        self._lenseq = None
        self._format = None
        self._closed = False
        self._readline = self._readlines = None

        if mode == 'r':
            filename = str(msa)
            if isfile(filename):
                self._filename = filename
                title, ext = splitext(split(msa)[1])
                if ext.lower() == '.gz':
                    title, ext = splitext(split(msa)[1])[0]
                self._title = title
                self._stream =  openFile(msa)
                self._readline = readline = self._stream.readline
                self._readlines = self._stream.readlines
            else:
                try:
                    self._readline = msa.readline
                except AttributeError:
                    pass
                else:
                    self._stream = msa
                    try:
                        self._readlines = msa.readlines
                    except AttributeError:
                        pass
                    
            if self._readline is not None: 
                self._firstline = line = readline()
                if line[0] == '>':
                    LOGGER.info('Opened MSA file in FASTA format')
                    self._iterator = self._iterFasta
                    self._format = FASTA
                elif line[0] == '#' and 'STOCKHOLM' in line:
                    LOGGER.info('Opened MSA file in Stockholm format')
                    self._iterator = self._iterStockholm
                    self._format = STOCKHOLM
                else:
                    LOGGER.info('Opened MSA file in SELEX format')
                    self._iterator = self._iterStockholm
                    self._format = SELEX

            else:        
                try:
                    bio = iter(msa)
                except TypeError:
                    raise TypeError('msa must be a filename, a stream, or '
                                    'Bio object')
                else:
                    record = bio.next()
                    try:
                        label, seq = record.id, str(record.seq)
                    except AttributeError:
                        raise TypeError('msa must be a filename, a stream,'
                                        ' or Bio object')
                    else:
                        self._bio = iter(msa)
                        self._iter = self._iterBio
                        self._format = 'Biopython'

        
            self._split = bool(kwargs.get('split', True))
            self.setFilter(kwargs.get('filter', None))
            self.setSlice(kwargs.get('slice', None))
            self._aligned = bool(aligned)
            self._iter = self._iterator()
        else:    
            pass
       
    def __del__(self):
        
        self.close()
            
    def __iter__(self):
        
        filter = self._filter
        slicer = self._slicer
        split = self._split
        if filter is None:
            for label, seq in self._iter:
                if split:
                    label, start, end = splitLabel(label)
                    yield label, slicer(seq), start, end
                else:
                    yield label, slicer(seq)
        else:
            for label, seq in self._iter:
                if filter(label, seq):
                    if split:
                        label, start, end = splitLabel(label)
                        yield label, slicer(seq), start, end    
                    else:
                        yield label, slicer(seq)
                             
    def __str__(self):
        
        return 'MSAFile ' + self._title 
    
    def __repr__(self):
        
        if self._closed:
            return '<MSAFile: {0:s} ({1:s}; closed)>'.format(self._title, 
                self._format) 
        else:
            return '<MSAFile: {0:s} ({1:s})>'.format(self._title, self._format) 
    
    def __enter__(self):

        return self

    def __exit__(self, type, value, tb):

        self.close()
    
    def _readlines(self, size=None):
        """Read multiple lines, in case stream does not have this method."""
        
        if self._closed:
            raise ValueError('I/O operation on closed file')
        import sys
        size = size or sys.maxint
        lines = []
        append = lines.append
        readline = stream.readline
        for i in range(size):
            append(readline)
        return lines
        
    def close(self):
        """Close the file.  This method will not affect a stream."""
        
        if self._filename:
            try:
                self._stream.close()
            except Exception:
                pass
            else:
                self._closed = True
    
    def _isClosed(self):
        
        return self._closed
    
    closed = property(_isClosed, None, doc="True for closed file.")
    
    def _getFormat(self):
        """Return format of the MSA file."""
        
        return self._format
    
    format = property(_getFormat, None, doc="Format of the MSA file.")
            
    
    def reset(self):
        """Return to the beginning of the file."""
        
        self._readline()
        self._iterator = self._iter()
    
    def _iterBio(self):
        """Yield sequences from a Biopython MSA object."""            
        
        aligned = self._aligned
        lenseq = self._lenseq
        numseq = 0
        
        for record in self._bio:
            label, seq = record.id, str(record.seq)
            if not lenseq:
                self._lenseq = lenseq = len(seq)
            if aligned and lenseq != len(seq):
                raise IOError('sequence for {0:s} does not have '
                              'expected length {1:d}'
                              .format(label, lenseq))
            numseq += 1
            yield label, seq

    def _iterFasta(self):
        """Yield sequences from a file or stream in FASTA format."""

        aligned = self._aligned
        lenseq = self._lenseq
        temp = []
        
        readlines = self._readlines
        line = self._firstline
        label = line[1:]
        lines = readlines(NUMLINES)
        while lines:
            for line in lines:
                if line[0] == '>':
                    seq = ESJOIN(temp)
                    if not lenseq:
                        self._lenseq = lenseq = len(seq)
                    if aligned and lenseq != len(seq):
                        raise IOError('sequence for {0:s} does not have '
                                      'expected length {1:d}'
                                      .format(label, lenseq))
                    yield label, seq
                    temp = []
                    label = line[1:].strip()
                else:
                    temp.append(line.strip())
            lines = readlines(NUMLINES)
        yield label, ESJOIN(temp)
            
    
    def _iterStockholm(self):
        """Yield sequences from an MSA file in Stockholm/SELEX format."""

        aligned = self._aligned
        lenseq = self._lenseq
        readlines = self._readlines

        lines = [self._firstline]
        lines.extend(readlines(NUMLINES))
        while lines:
            for line in lines:
                ch = line[0] 
                if ch == '#' or ch == '/':
                    continue
                items = line.split()
                if len(items) == 2:
                    label = items[0]
                    seq = items[1]
                else:
                    label = WSJOIN(items[:-1])
                    seq = items[-1]
                if not lenseq:
                    self._lenseq = lenseq = len(seq)
                if aligned and lenseq != len(seq):
                    raise IOError('sequence for {0:s} does not have '
                                  'expected length {1:d}'
                                  .format(label, lenseq))
                yield label, seq
            lines = readlines(NUMLINES)
    
    def getTitle(self):
        """Return title of the instance."""
        
        return self._title
    
    def setTitle(self, title):
        """Set title of the instance."""
        
        self._title = str(title)
    
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
            result = filter('TEST_TITLE', 'SEQUENCE-WITH-GAPS')
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
        
        return self._slice
    
    def setSlice(self, slice):
        """Set object used to slice sequences."""

        if slice is None:
            self._slice = None 
            self._slicer = lambda seq: seq
        else:
            seq = 'SEQUENCE' * 1000
            try: 
                result = seq[slice]
            except Exception:
                arr = fromstring(seq, '|S1')
                try:
                    result = arr[slice]
                except Exception:
                    raise TypeError('invalid slice: ' + repr(slice))
                else:
                    self._slice = slice
                    self._slicer = lambda seq, slc=slice: fromstring(seq,
                                                        '|S1')[slc].tostring()
            else:
                self._slice = slice
                self._slicer = lambda seq, slc=slice: seq[slc]


class MSA(object):
    
    """Store and manipulate multiple sequence alignments.
    
    >>> from prody import *
    >>> msafile = fetchPfamMSA('piwi', alignment='seed')
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
    
    Filtering and slicing available to :class:`MSAFile` class can be used to 
    parse an MSA selectively, which may be useful in low memory situations:
        
    >>> msa = MSA(msafile, filter=lambda lbl, seq: 'ARATH' in lbl, 
    ...           slice=list(range(10)) + list(range(394,404)))
    >>> msa
    <MSA: piwi_seed (3 sequences, 20 residues)>

    Compare this to result from parsing the complete file:
    >>> MSA(msafile)
    <MSA: piwi_seed (20 sequences, 404 residues)>"""
    
    def __init__(self, msa, **kwargs):
        """*msa* may be an :class:`MSAFile` instance or an MSA file in a 
        supported format."""
        
        try:
            ndim, dtype_, shape = msa.ndim, msa.dtype, msa.shape
        except AttributeError:
            try:
                numseq, lenseq = msa.numSequences, msa.numResidues
            except AttributeError:
                kwargs['split'] = False
                try:
                    msa = MSAFile(msa, **kwargs)
                except Exception as err:
                    raise TypeError('msa was not recognized ({0:s})'
                                    .format(str(err)))
            
            self._msa = []
            sappend = self._msa.append
            self._labels = []
            lappend = self._labels.append
            self._mapping = mapping = {}
            
            for i, (label, seq) in enumerate(msa):
                lappend(label)
                sappend(fromstring(seq, '|S1'))
                mapping[splitLabel(label)[0]] = i
            self._msa = array(self._msa, '|S1')
            self._title = kwargs.get('title', msa.getTitle())
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
    using :func:`.fetchPfamMSA` with default parameters.  Note that *msa* may be
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
    LOGGER.timeit('_parsemsa')
    msafile = MSAFile(filename)
    title = splitext(split(msa)[1])[0]
    format = msafile.format
    lenseq = len(next(iter(msafile))[1])
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
    msa = MSA(msa=msaarr, title=title, labels=labels, mapping=mapping)
    LOGGER.report('MSA of {1:d} residue long {0:d} sequence(s) was parsed in '
                  '%.2fs.'.format(*msaarr.shape), '_parsemsa') 
    return msa 

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
        
    from .msatools import calcShannonEntropy
    return calcShannonEntropy(msa, ambiguity=bool(ambiguity), 
                              omitgaps=bool(omitgaps))
    
    
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
        
    from .msatools import calcMutualInfo
    turbo = bool(turbo)
    LOGGER.timeit('_mutinfo')
    mutinfo = calcMutualInfo(msa, ambiguity=bool(ambiguity), 
                           turbo=turbo, debug=kwargs.get('debug', False))
    LOGGER.report('Mutual information matrix was calculated in %.2fs.', 
                  '_mutinfo')   
    return mutinfo


def calcMSAOccupancy(msa, occ='res'):
    """Return occupancy array calculated for residues (default) or sequences
    (``occ='seq'``) of *msa*, which may be an :class:`MSA` instance or a 2D 
    Numpy character array.  Implementation is case insensitive."""
    
    from .msatools import calcMSAOccupancy
    
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

    try:
        occ = occ.startswith('res')
    except AttributeError:
        raise TypeError('occ must be a string')
    return calcMSAOccupancy(msa, occ)
    
