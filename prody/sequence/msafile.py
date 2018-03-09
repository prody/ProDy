# -*- coding: utf-8 -*-
"""This module defines functions and classes for parsing, manipulating, and
analyzing multiple sequence alignments."""

__author__ = 'Anindita Dutta, Ahmet Bakan'

from os.path import isfile, splitext, split, getsize

from numpy import array, fromstring, empty

from .sequence import splitSeqLabel, Sequence

from prody import LOGGER
from prody.utilities import openFile

__all__ = ['MSAFile', 'splitSeqLabel', 'parseMSA', 'writeMSA']

FASTA = 'FASTA'
SELEX = 'SELEX'
STOCKHOLM = 'Stockholm'
CLUSTAL = 'CLUSTAL'
PIR = 'PIR'
MSAFORMATS = {
    FASTA.lower(): FASTA,
    SELEX.lower(): SELEX,
    STOCKHOLM.lower(): STOCKHOLM,
    CLUSTAL.lower(): CLUSTAL,
    PIR.lower(): PIR,
}
MSAEXTMAP = {
    FASTA: '.fasta',
    SELEX: '.slx',
    STOCKHOLM: '.sth',
    CLUSTAL: '.aln',
    PIR: '.ali',
    FASTA.lower(): '.fasta',
    SELEX.lower(): '.slx',
    STOCKHOLM.lower(): '.sth',
    CLUSTAL.lower(): '.aln', 
    PIR.lower(): '.ali',
    '.sth': STOCKHOLM,
    '.slx': SELEX,
    '.fasta': FASTA,
    '.aln': CLUSTAL,
    '.ali': PIR,
}

WSJOIN = ' '.join
ESJOIN = ''.join

NUMLINES = 1000
LEN_FASTA_LINE = 60
LEN_SELEX_LABEL = 31

try:
    range = xrange
except NameError:
    pass


class MSAFile(object):

    """Handle MSA files in FASTA, SELEX, CLUSTAL and Stockholm formats."""

    def __init__(self, msa, mode='r', format=None, aligned=True, **kwargs):
        """*msa* may be a filename or a stream.  Multiple sequence alignments
        can be read from or written in FASTA (:file:`.fasta`), Stockholm
        (:file:`.sth`), CLUSTAL (:file:`.aln`), or SELEX (:file:`.slx`) *format*.  
        For specified extensions, *format* argument is not needed. If *aligned* is
        **True**, unaligned sequences in the file or stream will cause an
        :exc:`IOError` exception.  *filter*, a function that returns a
        boolean, can be used for filtering sequences, see :meth:`setFilter`
        for details.  *slice* can be used to slice sequences, and is applied
        after filtering, see :meth:`setSlice` for details."""

        if mode[0] not in 'rwa':
            raise ValueError("mode string must be one of 'r', 'w', or 'a', "
                             "not {0}".format(repr(mode)))

        if 'b' in mode:
            mode = mode.replace('b', '')
        if 't' not in mode:
            mode += 't'

        self._format = None
        if format is not None:
            try:
                self._format = format = MSAFORMATS[format.lower()]
            except AttributeError:
                raise TypeError('format argument must be a string')
            except KeyError:
                raise ValueError('format argument is not recognized')

        self._filename = filename = None
        if mode.startswith('r'):
            try:
                torf = isfile(msa)
            except:
                pass
            else:
                if torf:
                    self._filename = filename = msa
        else:
            try:
                msa.lower, msa.strip
            except AttributeError:
                pass
            else:
                self._filename = filename = msa

        if filename is not None:
            self._filename = filename
            title, ext = splitext(split(filename)[1])
            if ext.lower() == '.gz':
                title, ext = splitext(split(title)[1])
            if format is None:
                try:
                    self._format = format = MSAEXTMAP[ext.lower()]
                except KeyError:
                    raise TypeError('format is not specified and could not be '
                                    'determined from file extension')
            self._title = title
            self._stream = openFile(msa, mode)

        else:
            if self._format is None:
                raise ValueError('format must be specified when msa is a '
                                 'stream')
            self._stream = msa
            self._title = 'stream'
            try:
                closed = self._stream.closed
            except AttributeError:
                closed = self._stream.myfileobj.closed
            if closed:
                raise ValueError('msa stream must not be closed')


        self._lenseq = None
        self._closed = False
        self._readline = None
        self._aligned = bool(aligned)

        if mode.startswith('r'):
            self._readline = self._stream.readline
            try:
                self._readlines = self._stream.readlines
            except AttributeError:
                pass

            self.setFilter(kwargs.get('filter', None),
                           kwargs.get('filter_full', False))
            self.setSlice(kwargs.get('slice', None))
            self._iter = self._itermap[format](self)
        else:
            try:
                self._write = write = self._stream.write
            except AttributeError:
                raise TypeError('msa must be a filename or a stream with '
                                'write method')
            if mode.startswith('w') and format == STOCKHOLM:
                write('# STOCKHOLM 1.0\n')
            if format.startswith('S'):
                self._selex_line = '{0:' + str(LEN_SELEX_LABEL) + 's} {1}\n'

        self._mode = mode

    def __del__(self):

        self.close()

    def __iter__(self):

        if not self._mode.startswith('r'):
            raise IOError('File not open for reading')

        filter = self._filter
        slicer = self._slicer
        if filter is None:
            for seq, label in self._iter:
                yield Sequence(slicer(seq), label)
        else:
            if self._filter_full:
                for seq, label in self._iter:
                    if filter(label, seq):
                        yield Sequence(slicer(seq), label)
            else:
                for seq, label in self._iter:
                    if filter(splitSeqLabel(label)[0], seq):
                        yield Sequence(slicer(seq), label)

    def __str__(self):

        return 'MSAFile ' + self._title

    def __repr__(self):

        if self._closed:
            return '<MSAFile: {0} ({1}; closed)>'.format(self._title,
                                                         self._format)
        else:
            return ('<MSAFile: {0} ({1}; mode {2})>'
                    ).format(self._title, self._format, repr(self._mode))

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
        readline = self._stream.readline
        for i in range(size):
            line = readline()
            if line:
                append(line)
            else:
                break
        return lines

    def close(self):
        """Close the file.  This method will not affect a stream."""

        if self._filename is None:
            self._closed = True
            return

        if not self._mode.startswith('r') and self._format == STOCKHOLM:
            try:
                self._write('//\n')
            except ValueError:
                LOGGER.info('Failed to write terminal slash characters to '
                            'closed file.')

        try:
            self._stream.close()
        except Exception:
            pass
        self._closed = True

    def _isClosed(self):

        return self._closed

    closed = property(_isClosed, None, doc="True for closed file.")

    def _getFormat(self):
        """Returns format of the MSA file."""

        return self._format

    format = property(_getFormat, None, doc="Format of the MSA file.")

    def reset(self):
        """Returns to the beginning of the file."""

        self._readline()
        self._iterator = self._iter()

    def isAligned(self):
        """Returns **True** if MSA is aligned."""

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
                raise IOError('sequence for {0} does not have '
                              'expected length {1}'
                              .format(label, lenseq))
            numseq += 1
            yield seq, label

    def _iterFasta(self):
        """Yield sequences from a file or stream in FASTA format."""

        aligned = self._aligned
        lenseq = self._lenseq
        temp = []

        label = ''
        lines = []
        while not label.startswith('>'):
            if not lines:
                lines = self._readlines(NUMLINES)
            label = lines.pop(0)
        label = label[1:].strip()

        while lines:
            for line in lines:
                if line.startswith('>'):
                    seq = ESJOIN(temp)
                    if not lenseq:
                        self._lenseq = lenseq = len(seq)
                    if aligned and lenseq != len(seq):
                        raise IOError('sequence for {0} does not have '
                                      'expected length {1}'
                                      .format(label, lenseq))
                    yield seq, label
                    temp = []
                    label = line[1:].strip()
                else:
                    temp.append(line.strip())
            lines = self._readlines(NUMLINES)
        yield ESJOIN(temp), label

    def _iterSelex(self):
        """Yield sequences from an MSA file in Stockholm/SELEX format."""

        aligned = self._aligned
        lenseq = self._lenseq
        readlines = self._readlines

        lines = readlines(NUMLINES)
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
                    raise IOError('sequence for {0} does not have '
                                  'expected length {1}'
                                  .format(label, lenseq))
                yield seq, label
            lines = readlines(NUMLINES)

    _itermap = {
        FASTA: _iterFasta,
        SELEX: _iterSelex,
        STOCKHOLM: _iterSelex,
    }

    def getTitle(self):
        """Returns title of the instance."""

        return self._title

    def setTitle(self, title):
        """Set title of the instance."""

        self._title = str(title)

    def getFilename(self):
        """Returns filename, or **None** if instance is handling a stream."""

        return self._filename

    def getFormat(self):
        """Returns file format."""

        return self._format

    def getFilter(self):
        """Returns function used for filtering sequences."""

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
        """Returns object used to slice sequences."""

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
                seq[slice]
            except Exception:
                arr = fromstring(seq, '|S1')
                try:
                    arr[slice]
                except Exception:
                    raise TypeError('invalid slice: ' + repr(slice))
                else:
                    self._slice = slice
                    self._slicer = lambda seq, slc=slice: fromstring(seq,
                                                        '|S1')[slc].tostring()
            else:
                self._slice = slice
                self._slicer = lambda seq, slc=slice: seq[slc]

    def write(self, seq):
        """Write *seq*, an :class:`.Sequence` instance, into the MSA file."""

        if self._closed:
            raise ValueError('I/O operation on closed file')
        write = self._write

        try:
            label, sequence = seq.getLabel(True), str(seq)
        except AttributeError:
            raise TypeError('seq must be a Sequence instance')

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
    """Returns an :class:`.MSA` instance that stores multiple sequence alignment
    and sequence labels parsed from Stockholm, SELEX, CLUSTAL, PIR, or FASTA format
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
    title = split(title)[1]
    aligned = kwargs.get('aligned', True)
    if (ext.lower() == '.gz' or 'filter' in kwargs or 'slice' in kwargs or
            not aligned):
        if ext.lower() == '.gz':
            title = splitext(title)[0]
        msa = MSAFile(filename, split=False, **kwargs)
        seqlist = []
        sappend = seqlist.append
        labels = []
        lappend = labels.append
        mapping = {}
        maxlen = 0
        for i, seq in enumerate(msa):
            label = seq.getLabel(True)
            lappend(label)
            if aligned:
                sappend(seq._array)
            else:
                if len(seq) > maxlen:
                    maxlen = len(seq)
                sappend(seq)
            key = splitSeqLabel(label)[0]
            if key in mapping:
                try:
                    mapping[key].append(i)
                except AttributeError:
                    mapping[key] = [mapping[key], i]
            else:
                mapping[key] = i
        if not seqlist:
            LOGGER.warn('No sequences were parsed from {0}.'.format(filename))
            return
        if aligned:
            msaarr = array(seqlist, '|S1')
        else:
            msaarr = array(seqlist, '|S' + str(maxlen))
    else:
        filesize = getsize(filename)
        format = MSAEXTMAP.get(splitext(filename)[1])

        if format == FASTA:
            from .msaio import parseFasta as parser
            msaarr = empty(filesize, '|S1')
        elif format == SELEX or format == STOCKHOLM:
            from .msaio import parseSelex as parser
            msaarr = empty(filesize, '|S1')
        elif format == CLUSTAL:
            parser = parseClustal
            msaarr = []
        elif format == PIR:
            parser = parsePIR
            msaarr = []
        else:
            raise IOError('MSA file format is not recognized from the '
                          'extension')
        msaarr, labels, mapping, lcount = parser(filename, msaarr)
        if lcount != len(msaarr):
            LOGGER.warn('Failed to parse {0} sequence labels.'
                        .format(len(msaarr) - lcount))

    msa = MSA(msa=msaarr, title=title, labels=labels, mapping=mapping,
              aligned=aligned)

    if aligned:
        LOGGER.report('{0} sequence(s) with {1} residues were parsed in '
                      '%.2fs.'.format(*msaarr.shape), '_parsemsa')
    else:
        LOGGER.report('{0} sequence(s) were parsed in %.2fs.'
                      .format(*msaarr.shape), '_parsemsa')
    return msa

def parseClustal(filename, msaarr):
    """
    Parses a CLUSTAL format (:file:`.aln`) alignment file.
    """
    msafile = open(filename,'r')
    lines = msafile.readlines()
    msafile.close()

    msa_dict = {}
    keys = []
    for line in lines:
        foundBadItem = False
        try:
            key = line.strip().split()[0]
            seq = line.strip().split()[1]
        except:
            continue
        for badItem in ['*', ' ', ':', 'CLUSTAL', '.']:
            if badItem in key:
                foundBadItem = True
                continue
        if foundBadItem:
            continue
        if not key in keys:
            keys.append(key)
            msa_dict[key] = seq
        else:
            msa_dict[key] = msa_dict[key] + seq

    for key in keys:
        msaarr.append(list(msa_dict[key]))

    return array(msaarr), keys, None, len(keys)

def parsePIR(filename, msaarr):
    msafile = open(filename,'r')
    lines = msafile.readlines()
    msafile.close()

    labels = []
    i = -1

    for line in lines:
        if line.startswith('>P1;'):
            labels.append(line.strip()[len('>P1;'):])
            i += 1
            msaarr.append([])
        elif line.startswith('s') or line.strip() == '':
            pass
        else:
            msaarr[i].append(line.strip())

    for i in range(len(msaarr)):
        msaarr[i] = list(''.join(msaarr[i]))

    return array(msaarr), labels, None, len(labels)

def writeClustal(filename, msa):
    """A simple writer for CLUSTAL format alignments.

    This lacks the characters showing degree of conservation
    but otherwise conforms to the CLUSTAL format standards."""

    msafile = open(filename, 'w')

    msafile.write('CLUSTALW file written by ProDy\n\n')

    for j in range(msa.numResidues()/60):
        for i in range(msa.numSequences()):
            sequence = str(msa[i])
            msafile.write(msa.getLabel(i) + ' '*(16-len(msa.getLabel(i))))
            msafile.write(sequence[j*60:(j+1)*60])
            msafile.write('\n')
        msafile.write('\n\n')

    for i in range(msa.numSequences()):
        sequence = str(msa[i])
        msafile.write(msa.getLabel(i) + ' '*(16-len(msa.getLabel(i))))
        msafile.write(sequence[(j+1)*60:])
        msafile.write('\n')

    msafile.close()
    return

def writePIR(filename, msa, **kwargs):
    """A function to write PIR format alignments for use with MODELLER.

    :arg filename: The name of the file to be written including .ali
    :type filename: str

    :arg msa: a multiple sequence alignment in :class:`MSA` format
    :type msa: :class:`MSA` instance

    :arg chain_sep: chain separation character or list of them
        default is '/'
    :type chain_sep: str, list

    :arg types: a list of strings for field 1, PIR types (Sequence or StructureX)
        default is all Sequence
    :type types: list

    :arg labels: a list of strings for field 2, sequence labels
        default is to take them from msa
    :type labels: list

    :arg first_resnums: contents for field 3, residue number for the first residue.
        This should be a list of strings each having length 5, 
        default is all 'FIRST'
    :type first_resnums: list

    :arg first_chains: contents for field 4, chain ID for the first residue
        This should be a list of strings each having length 1, 
        default is all '@'
    :type first_chains: list

    :arg last_resnums: contents for field 5, residue number for the last residue.
        This should be a list of strings each having length 5, 
        default is all 'LAST '
    :type last_resnums: list

    :arg last_chains: contents for field 6, chain ID for the last residue
        This should be a list of strings each having length 1, 
        default is all ' '
    :type first_chains: list

    :arg protein_names: list of strings for field 7
        default is all ''
    :type protein_names: list

    :arg protein_sources: list of strings for field 8
        default is all ''
    :type protein_sources: list

    :arg resolutions: list of strings for field 9
        default is all ''
    :type resolutions: list

    :arg r_factors: list of strings for field 10
        default is all ''
    :type r_factors: list
    """
    msafile = open(filename, 'w')

    chain_sep = kwargs.get('types', '/')
    if isinstance(chain_sep, str): 
        chain_sep = [chain_sep] * msa.numSequences()
    elif isinstance(chain_sep, list) and isinstance(chain_sep[0], str):
        if len(chain_sep) != msa.numSequences():
            raise ValueError('There should be an entry in chain_sep list for each sequence in msa')
    else:
        raise TypeError('chain_sep should be a string or list of strings')

    types = kwargs.get('types', 'Sequence')
    if types is None: 
        types = [types] * msa.numSequences()
    elif isinstance(types, list) and isinstance(types[0], str):
        if len(types) != msa.numSequences():
            raise ValueError('There should be an entry in types list for each sequence in msa')
    else:
        raise TypeError('types should be a string or list of strings')

    labels = kwargs.get('labels', None)
    if labels is None: 
        labels = []
        for sequence in msa:
            labels.append(sequence.getLabel())
    elif isinstance(labels, list) and isinstance(labels[0], str):
        if len(labels) != msa.numSequences():
            raise ValueError('There should be an entry in labels list for each sequence in msa')
    else:
        raise TypeError('labels should be a string or list of strings')

    first_resnums = kwargs.get('first_resnums', 'FIRST')
    if isinstance(first_resnums, str) and len(first_resnums) == 5: 
        first_resnums = [first_resnums] * msa.numSequences()
    elif isinstance(first_resnums, list) and isinstance(first_resnums, str):
        if len(first_resnums) != msa.numSequences():
            raise ValueError('There should be an entry in first_resnums list for each sequence in msa')
    else:
        raise TypeError('first_resnums should be a string of length 5 or list of them')

    first_chains = kwargs.get('first_chains', '@')
    if isinstance(first_chains, str) and len(first_chains) == 1: 
        first_chains = [first_chains] * msa.numSequences()
    elif isinstance(first_chains, list) and isinstance(first_chains, str):
        if len(first_chains) != msa.numSequences():
            raise ValueError('There should be an entry in first_chains list for each sequence in msa')
    else:
        raise TypeError('first_chains should be a string of length 1 or list of them')

    last_resnums = kwargs.get('last_resnums', 'LAST ')
    if isinstance(last_resnums, str) and len(last_resnums) == 5: 
        last_resnums = [last_resnums] * msa.numSequences()
    elif isinstance(last_resnums, list) and isinstance(last_resnums, str):
        if len(last_resnums) != msa.numSequences():
            raise ValueError('There should be an entry in last_resnums list for each sequence in msa')
    else:
        raise TypeError('last_resnums should be a string of length 5 or list of them')

    last_chains = kwargs.get('last_chains', ' ')
    if isinstance(last_chains, str) and len(last_chains) == 1: 
        last_chains = [last_chains] * msa.numSequences()
    elif isinstance(last_chains, list) and isinstance(last_chains, str):
        if len(last_chains) != msa.numSequences():
            raise ValueError('There should be an entry in last_chains list for each sequence in msa')
    else:
        raise TypeError('last_chains should be a string of length 1 or list of them')

    protein_names = kwargs.get('protein_names', '')
    if isinstance(protein_names, str): 
        protein_names = [protein_names] * msa.numSequences()
    elif isinstance(protein_names, list) and isinstance(protein_names, str):
        if len(protein_names) != msa.numSequences():
            raise ValueError('There should be an entry in protein_names list for each sequence in msa')
    else:
        raise TypeError('protein_names should be a string or list of strings')

    protein_sources = kwargs.get('protein_sources', '')
    if isinstance(protein_sources, str): 
        protein_sources = [protein_sources] * msa.numSequences()
    elif isinstance(protein_sources, list) and isinstance(protein_sources, str):
        if len(protein_sources) != msa.numSequences():
            raise ValueError('There should be an entry in protein_sources list for each sequence in msa')
    else:
        raise TypeError('protein_sources should be a string or list of strings')

    resolutions = kwargs.get('resolutions', '')
    if isinstance(resolutions, str): 
        resolutions = [resolutions] * msa.numSequences()
    elif isinstance(resolutions, list) and isinstance(resolutions, str):
        if len(resolutions) != msa.numSequences():
            raise ValueError('There should be an entry in resolutions list for each sequence in msa')
    else:
        raise TypeError('resolutions should be a string or list of strings')

    r_factors = kwargs.get('r_factors', '')
    if isinstance(r_factors, str): 
        r_factors = [r_factors] * msa.numSequences()
    elif isinstance(r_factors, list) and isinstance(r_factors, str):
        if len(r_factors) != msa.numSequences():
            raise ValueError('There should be an entry in r_factors list for each sequence in msa')
    else:
        raise TypeError('r_factors should be a string or list of strings')

    for i, sequence in enumerate(msa):
        sequence = str(sequence).replace(chain_sep[i],'/')
        msafile.write('>P1;' + labels[i] + '\n')
        msafile.write(types[i] + ':' + labels[i] + ':')
        msafile.write(first_resnums[i] + ':' + first_chains[i] + ':')
        msafile.write(last_resnums[i] + ':' + last_chains[i] + ':')
        msafile.write(protein_names[i] + ':' + protein_sources[i] + ':')
        msafile.write(resolutions[i] + ':' + r_factors[i])
        msafile.write('\n')

        for j in range(len(sequence)/60):
            msafile.write(sequence[j*60:(j+1)*60] + '\n')
        msafile.write(sequence[(j+1)*60:] + '*\n\n')

    msafile.close()
    return

def writeMSA(filename, msa, **kwargs):
    """Returns *filename* containing *msa*, a :class:`.MSA` or :class:`.MSAFile`
    instance, in the specified *format*, which can be *SELEX*, *Stockholm*, or
    *FASTA*.  If *compressed* is **True** or *filename* ends with :file:`.gz`,
    a compressed file will be written.  :class:`.MSA` instances will be written
    using C function into uncompressed files.
    
    Can also write *CLUSTAL* or *PIR* format files using Python functions."""

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
            raise ValueError('format {0} is not recognized'
                             .format(repr(format)))
    else:
        if ext == '.gz':
            ext = splitext(fntemp)[1].lower()

        try:
            format = MSAEXTMAP[ext]
        except KeyError:
            raise ValueError('format is not specified, and file extension '
                             '{0} is not recognized'.format(repr(ext)))

    fast = False
    try:
        seqarr, _, _ = msa._getArray(), msa.numResidues(), msa.numSequences()
    except AttributeError:
        try:
            msa.getFormat(), msa.getFilename(), msa.getFilter()
        except AttributeError:
            raise ValueError('msa must be an MSA or MSAFile instance, not '
                             .format(type(msa).__name__))
        else:
            seqiter = msa

    else:
        seqiter = msa
        fast = True

    if not fast or compressed:
        with MSAFile(filename, 'w', format=format) as out:
            write = out.write
            [write(seq) for seq in seqiter]
    else:
        from prody.utilities import backupFile
        backupFile(filename)
        if format == FASTA:
            from .msaio import writeFasta
            writeFasta(filename, msa._labels, seqarr,
                       kwargs.get('line_length', LEN_FASTA_LINE))
        elif format == CLUSTAL:
            writeClustal(filename, msa)
        elif format == PIR:
            writePIR(filename, msa, **kwargs)
        else:
            from .msaio import writeSelex
            writeSelex(filename, msa._labels, seqarr,
                       stockholm=format != SELEX,
                       label_length=kwargs.get('label_length',
                                               LEN_SELEX_LABEL))
    return filename
