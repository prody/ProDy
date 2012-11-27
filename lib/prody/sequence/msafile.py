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

__all__ = ['MSAFile', 'splitSeqLabel', 'parseMSA', 'writeMSA']

FASTA = 'FASTA'
SELEX = 'SELEX'
STOCKHOLM = 'Stockholm'
MSAFORMATS = {
    FASTA.lower(): FASTA,
    SELEX.lower(): SELEX,
    STOCKHOLM.lower(): STOCKHOLM,
}
EXTENSIONS = {
    FASTA: '.fasta',
    SELEX: '.slx',
    STOCKHOLM: '.sth',
    FASTA.lower(): '.fasta',
    SELEX.lower(): '.slx',
    STOCKHOLM.lower(): '.sth',
    '.sth': STOCKHOLM,
    '.slx': SELEX,
    '.fasta': FASTA
}

WSJOIN = ' '.join
ESJOIN = ''.join

NUMLINES = 1000
LEN_FASTA_LINE = 60
LEN_SELEX_LABEL = 31

import re
from os.path import isfile, splitext, split, getsize

from numpy import all, zeros, dtype, array, char, fromstring

from prody import LOGGER, PY2K
from prody.utilities import openFile

if PY2K: range = xrange

SPLITLABEL = re.compile('/*-*').split 


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


class MSAFile(object):
    
    """Handle MSA files in FASTA, SELEX and Stockholm formats. 
    
    >>> from prody import * 
    >>> msafile = fetchPfamMSA('piwi', alignment='seed')

    *Reading a file*

    Iterating over a file will yield sequence id, sequence, residue start and 
    end indices:

    >>> msa = MSAFile(msafile)
    >>> for seq in msa: # doctest: +ELLIPSIS 
    ...     print(seq)
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
    ...     print(seq)
    ('AGO6_ARATH', 'FIL...FTK', 541, 851)
    ('AGO4_ARATH', 'FIL...FMK', 577, 885)
    ('AGO10_ARATH', 'LLL...YLE', 625, 946)

    *Slicing sequences*
    
    A list of integers can be used to slice sequences as follows.
    
    >>> msa = MSAFile(msafile, slice=list(range(10)) + list(range(394,404)))
    >>> for seq in msa: # doctest: +ELLIPSIS 
    ...     print(seq)
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
        boolean, can be used for filtering sequences, see :meth:`setFilter`
        for details.  *slice* can be used to slice sequences, and is applied 
        after filtering, see :meth:`setSlice` for details."""
        
        if mode not in 'rwa':
            raise ValueError("mode string must be one of 'r', 'w', or 'a', "
                             "not {0}".format(repr(mode)))
        
        self._bio = self._filename = None
        self._lenseq = None
        self._format = None
        self._closed = False
        self._readline = self._readlines = None
        self._aligned = bool(aligned)

        if mode == 'r':
            filename = str(msa)
            if isfile(filename):
                self._filename = filename
                title, ext = splitext(split(msa)[1])
                if ext.lower() == '.gz':
                    title, ext = splitext(split(title)[1])
                self._title = title
                self._stream =  openFile(msa)
                self._readline = readline = self._stream.readline
                self._readlines = self._stream.readlines
            else:
                try:
                    self._readline = readline = msa.readline
                except AttributeError:
                    pass
                else:
                    self._stream = msa
                    try:
                        self._readlines = msa.readlines
                    except AttributeError:
                        pass
                    self._title = 'input stream'
                    
            if readline is not None: 
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
            self.setFilter(kwargs.get('filter', None),
                           kwargs.get('filter_full', False))
            self.setSlice(kwargs.get('slice', None))
            self._iter = self._iterator()
        else:
            if format is not None:
                try:
                    format = format.lower()
                except AttributeError:
                    raise ValueError('format must be a string')

                if format.startswith('f'):
                    self._format = format = FASTA
                elif format.startswith('st'):
                    self._format = format = STOCKHOLM
                else:
                    self._format = format = SELEX
            try:
                self._write = write = msa.write
            except AttributeError:
                try:
                    title, ext = splitext(split(msa)[1])
                    if ext:
                        if ext.lower() == '.gz':
                            title, ext = splitext(split(title)[1])
                    self._title = title
                    if ext:
                        if not format:
                            if ext.startswith('.f'):
                                self._format = format = FASTA
                            elif ext.startswith('.st'):
                                self._format = format = STOCKHOLM
                            else:
                                self._format = format = SELEX
                except Exception:
                    raise TypeError('msa must be a file name or a stream')
                else:
                    if not ext:
                        msa +=  EXTENSIONS[self._format]
                
                self._filename = msa
                self._stream = stream = openFile(msa, mode)
                self._write = write = stream.write
                
            else:
                try:
                    closed = msa.closed
                except AttributeError:
                    pass
                else:
                    if closed:
                        raise ValueError('msa stream must not be closed')
                self._title = 'output stream'
            
            if mode == 'w' and format == STOCKHOLM:
                write('# STOCKHOLM 1.0\n')
            if format.startswith('S'):
                self._selex_line = '{0:' + str(LEN_SELEX_LABEL) + 's} {1:s}\n'
        self._mode = mode
       
    def __del__(self):
        
        self.close()
            
    def __iter__(self):
        
        if self._mode != 'r':
            raise IOError('File not open for reading')
        
        filter = self._filter
        slicer = self._slicer
        split = self._split
        if filter is None:
            for label, seq in self._iter:
                if split:
                    label, start, end = splitSeqLabel(label)
                    yield label, slicer(seq), start, end
                else:
                    yield label, slicer(seq)
        else:
            if self._filter_full:
                for label, seq in self._iter:
                    if filter(label, seq):
                        if split:
                            label, start, end = splitSeqLabel(label)
                            yield label, slicer(seq), start, end    
                        else:
                            yield label, slicer(seq)
            else:                             
                for label, seq in self._iter:
                    lbl, start, end = splitSeqLabel(label)
                    if filter(lbl, seq):
                        if split:
                            yield lbl, slicer(seq), start, end    
                        else:
                            yield label, slicer(seq)
                
        
    def __str__(self):
        
        return 'MSAFile ' + self._title 
    
    def __repr__(self):
        
        if self._closed:
            return '<MSAFile: {0:s} ({1:s}; closed)>'.format(self._title, 
                self._format) 
        else:
            return '<MSAFile: {0:s} ({1:s}; mode {2:s})>'.format(
                    self._title, self._format, repr(self._mode)) 
    
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
            if not self._stream.closed:
                if (self._mode != 'r' and self._format == STOCKHOLM):
                    self._write('//\n')
                try:
                    self._stream.close()
                except Exception:
                    pass
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
    
    def isAligned(self):
        """Return **True** if MSA is aligned."""
        
        return self._aligned
    
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
        
    def getFilename(self):
        """Return filename, or **None** if instance is handling a stream."""
        
        return self._filename
    
    def getFormat(self):
        """Return file format."""
        
        return self._format
    
    def getFilter(self):
        """Return function used for filtering sequences."""
        
        return self._filter
    
    def setFilter(self, filter, filter_full=False):
        """Set function used for filtering sequences.  *filter* will be applied
        to split sequence label, by default.  If *filter_full* is **True**,
        filter will be applied to the full label.  """
        
        self._filter_full = bool(filter_full)
        
        if filter is None:
            self._filter = None
            return
        
        if not callable(filter):
            raise TypeError('filter must be callable')
        
        try: 
            result = filter('TEST_TITLE', 'SEQUENCE-WITH-GAPS')
        except Exception as err:
            raise TypeError('filter function must not raise exceptions, '
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
        """Set object used to *slice* sequences, which may be a :func:`slice`
        or a :func:`list` of numbers."""

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

    def write(self, label, sequence):
        """Write *label* and *sequence* into the MSA file."""
        
        if self._closed:
            raise ValueError('I/O operation on closed file')
        write = self._write
        
        try:
            _, _ = label.upper, sequence.upper
        except AttributeError:
            raise TypeError('label and sequence must be sequences')
        
        if self._lenseq is None:
            lenseq = self._lenseq = len(sequence)
        else:
            lenseq = self._lenseq
            if self._aligned and lenseq != self._lenseq:
                raise ValueError('writing an aligned MSA file, '
                                 'len(sequence) must be ' + str(lenseq))
        if self._format == FASTA:
            write('>')
            write(label)
            write('\n')
            beg = 0 
            end = LEN_FASTA_LINE
            lenseq = len(sequence)
            while beg < lenseq:
                write(sequence[beg:end])
                write('\n')
                beg += LEN_FASTA_LINE
                end += LEN_FASTA_LINE
        else:
            write(self._selex_line.format(label, sequence))


def parseMSA(filename, **kwargs):
    """Return an :class:`.MSA` instance that stores multiple sequence alignment
    and sequence labels parsed from Stockholm, SELEX, or FASTA format 
    *filename* file, which may be a compressed file. Uncompressed MSA files 
    are parsed using C code at a fraction of the time it would take to parse 
    compressed files in Python."""
    
    from .msa import MSA

    try:
        fileok = isfile(filename)
    except TypeError:
        raise TypeError('filename must be a string')
    else:
        if not fileok:
            raise IOError('[Errno 2] No such file or directory: ' + 
                          repr(filename))
    
    # if MSA is a compressed file or filter/slice is passed, use 
    #   Python parsers

    LOGGER.timeit('_parsemsa')
    
    title, ext = splitext(filename)
    aligned = kwargs.get('aligned', True)
    if (ext.lower() == '.gz' or 'filter' in kwargs or 'slice' in kwargs or
        not aligned):
        if ext.lower() == '.gz':
            title = spliext(title)
        msa = MSAFile(filename, split=False, **kwargs)
        seqlist = []
        sappend = seqlist.append
        labels = []
        lappend = labels.append
        mapping = {}
        maxlen = 0
        for i, (label, seq) in enumerate(msa):
            lappend(label)
            if aligned:
                sappend(fromstring(seq, '|S1'))
            else:
                if len(seq) > maxlen:
                    maxlen = len(seq)
                sappend(seq)
            mapping[splitSeqLabel(label)[0]] = i
        if aligned:
            msaarr = array(seqlist, '|S1')
        else:
            msaarr = array(seqlist, '|S' + str(maxlen))
    else:
        msafile = MSAFile(filename)
        format = msafile.format
        lenseq = len(next(iter(msafile))[1])
        numseq = getsize(filename) / (lenseq + 10)
        del msafile
        
        if format == FASTA:
            from .msaio import parseFasta
            msaarr, labels, mapping = parseFasta(filename, lenseq, numseq)
        elif format == SELEX or format == STOCKHOLM:
            from .msaio import parseSelex
            msaarr, labels, mapping = parseSelex(filename, lenseq, numseq)
        else:
            raise IOError('MSA file format is not recognized')
    msa = MSA(msa=msaarr, title=title, labels=labels, mapping=mapping)

    if aligned:
        LOGGER.report('{0:d} sequence(s) with {1:d} residues were parsed in '
                      '%.2fs.'.format(*msaarr.shape), '_parsemsa') 
    else:
        LOGGER.report('{0:d} sequence(s) were parsed in %.2fs.'
                      .format(*msaarr.shape), '_parsemsa') 
    return msa


def writeMSA(filename, msa, **kwargs):
    """Return *filename* containing *msa*, a :class:`.MSA` or :class:`.MSAFile`
    instance, in the specified *format*, which can be *SELEX*, *Stockholm*, or 
    *FASTA*.  If *compressed* is **True** or *filename* ends with :file:`.gz`, 
    a compressed file will be written.  :class:`.MSA` instances will be written
    using C function into uncompressed files."""
    
    fntemp, ext = splitext(filename)
    ext = ext.lower() 
    compressed = kwargs.get('compressed', ext == '.gz')
    if compressed and ext != '.gz':
        filename += '.gz'
    format = kwargs.get('format', None)
    if format:
        try:
            format = MSAFORMATS[format.lower()]
        except KeyError:
            raise ValueError('format {0:s} is not recognized'
                             .format(repr(format)))
    else:
        if ext == '.gz':
            ext = splitext(fntemp)[1].lower()
        
        try:
            format = EXTENSIONS[ext]
        except KeyError:
            raise ValueError('format is not specified, and file extension '
                             '{0:s} is not recognized'.format(repr(ext)))
        
    fast = False
    try:
        seqarr, _, _ = msa._getArray(), msa.numResidues(), msa.numSequences()
    except AttributeError:
        try:
            _ = msa.getFormat(), msa.getFilename(), msa.getFilter()
        except AttributeError:
            raise ValueError('msa must be an MSA or MSAFile instance, not '
                             .format(type(msa).__name__))
        else:
            seqiter = msa
            
    else:
        seqiter = msa
        fast = True
        
    if not fast or compressed:
        split, seqiter._split = seqiter._split, False
        with MSAFile(filename, 'w', format=format) as out:
            write = out.write
            [write(label, seq) for label, seq in seqiter]
        seqiter._split = split
    else:
        from prody.utilities import backupFile
        backupFile(filename)
        if format == FASTA:
            from .msaio import writeFasta
            writeFasta(filename, msa._labels, seqarr, 
                       kwargs.get('line_length', LEN_FASTA_LINE))
        else:
            from .msaio import writeSelex
            writeSelex(filename, msa._labels, seqarr, 
                   stockholm=format != SELEX,
                   label_length=kwargs.get('label_length', LEN_SELEX_LABEL))
    return filename
