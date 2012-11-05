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

"""This module defines a function for fetching a Multiple Sequence Alignment
(MSA) from Pfam"""

__author__ = 'Anindita Dutta, Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Anindita Dutta, Ahmet Bakan'

__all__ = ['fetchPfamMSA', 'MSAFile', 'MSA', 'parseMSA',
           'calcShannonEntropy']

FASTA = 'fasta'
SELEX = 'selex'
STOCKHOLM = 'stockholm'

DOWNLOAD_FORMATS = set(['seed', 'full', 'ncbi', 'metagenomics'])
FORMAT_OPTIONS = ({'format': set([FASTA, SELEX, STOCKHOLM]),
                  'order': set(['tree', 'alphabetical']),
                  'inserts': set(['lower', 'upper']),
                  'gaps': set(['mixed', 'dots', 'dashes', 'none'])})
WSJOIN = ' '.join
ESJOIN = ''.join

import re
from os.path import join, isfile, getsize, splitext, split

from numpy import zeros, dtype

from prody import LOGGER
from prody.utilities import makePath, openURL, gunzip, openFile


def fetchPfamMSA(acc, alignment='full', folder='.', compressed=False, 
                 **kwargs):
    """Return a path to the downloaded Pfam MSA file.
        
    :arg acc: Pfam ID or Accession Code
    :type acc: str
    
    :arg alignment: alignment type, one of ``'full'`` (default), ``'seed'``,
        ``'ncbi'``, or ``'metagenomics'`` 
    
    :arg folder: output folder, default is ``'.'``
    
    :arg compressed: gzip the downloaded MSA file, default is **False**
    
    *Alignment Options*
    
    :arg format: a Pfam supported MSA file format, one of ``'selex'``,  
        (default), ``'stockholm'`` or ``'fasta'``
    
    :arg order: ordering of sequences, ``'tree'`` (default) or 
        ``'alphabetical'`` 
    
    :arg inserts: letter case for inserts, ``'upper'`` (default) or ``'lower'``
    
    :arg gaps: gap character, one of ``'dashes'`` (default), ``'dots'``, 
        ``'mixed'`` or **None** for unaligned
    
    *Other Options*
    
    :arg timeout: timeout for blocking connection attempt in seconds, default 
        is 5

    :arg outname: out filename, default is input ``'acc_alignment.format'``"""   
   
    getAccUrl = 'http://pfam.sanger.ac.uk/family/acc?id=' + acc
    handle = openURL(getAccUrl)
    orig_acc = acc
    acc = handle.readline().strip()
    url_flag = False
    
    if not re.search('(?<=PF)[0-9]{5}$', acc):
        raise ValueError('No such family: check Pfam ID or Accession Code')
        
    
    if alignment not in DOWNLOAD_FORMATS:
        raise ValueError('alignment must be one of full, seed, ncbi or'
                         ' metagenomics')
    if alignment == 'ncbi' or alignment == 'metagenomics':
        url = ('http://pfam.sanger.ac.uk/family/' + acc + '/alignment/' +
               alignment + '/gzipped')
        url_flag = True
        extension = '.sth'
    else:
        if not kwargs:
            url = ('http://pfam.sanger.ac.uk/family/' + acc + '/alignment/' +
                   alignment + '/gzipped')
            url_flag = True
            extension = '.sth'
        else:
            align_format = kwargs.get('format', 'selex').lower()
            
            if align_format not in FORMAT_OPTIONS['format']:
                raise ValueError('alignment format must be of type selex'
                                 ' stockholm or fasta. MSF not supported')
            
            if align_format == SELEX:
                align_format, extension = 'pfam', '.slx'
            elif align_format == FASTA:
                extension = '.fasta'
            else:
                extension = '.sth'
            
            gaps = str(kwargs.get('gaps', 'dashes')).lower()
            if gaps not in FORMAT_OPTIONS['gaps']:
                raise ValueError('gaps must be of type mixed, dots, dashes, '
                                 'or None')
            
            inserts = kwargs.get('inserts', 'upper').lower()
            if(inserts not in FORMAT_OPTIONS['inserts']):
                raise ValueError('inserts must be of type lower or upper')
            
            order = kwargs.get('order', 'tree').lower()
            if order not in FORMAT_OPTIONS['order']:
                raise ValueError('order must be of type tree or alphabetical')
            
            url = ('http://pfam.sanger.ac.uk/family/' + acc + '/alignment/'
                   + alignment + '/format?format=' + align_format + 
                   '&alnType=' + alignment + '&order=' + order[0] +
                   '&case=' + inserts[0] + '&gaps=' + gaps + '&download=1')
        
    
    response =  openURL(url, timeout=int(kwargs.get('timeout', 5)))
    outname = kwargs.get('outname', None)
    if not outname:
        outname = orig_acc
    filepath = join(makePath(folder), outname + '_' + alignment + extension)
    if compressed:
        filepath = filepath + '.gz'
        if url_flag:
            f_out = open(filepath, 'wb')
        else:
            f_out = openFile(filepath, 'wb')
        f_out.write(response.read())
        f_out.close()
    else:        
        if url_flag:
            gunzip(response.read(), filepath)
        else:
            with open(filepath, 'wb') as f_out:
                f_out.write(response.read())
                
    LOGGER.info('Pfam MSA for {0} is written as {1}.'
                .format(orig_acc, filepath))
    
    return filepath

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
    
    >>> msafile = fetchPfamMSA('piwi', alignment='seed')
    >>> msa = MSAFile(msafile, filter=lambda tpl: 'ARATH' in tpl[0])
    >>> for label, seq, rns, rne in msa: 
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
        slice = self._slice
        numseq = 0
        for record in self._bio:
            label, seq = record.id, str(record.seq)
            if aligned and lenseq != len(seq):
                raise IOError('sequence for {0:s} does not have '
                              'expected length {1:d}'
                              .format(tpl[0], lenseq))
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
                                      .format(tpl[0], lenseq))
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
                                  .format(tpl[0], lenseq))
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
        
        if not iscallable(filter):
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
    
    >>> msa = parseMSA('piwi', alignment='seed')
    >>> msa[0]
    >>> msa[0,0]
    >>> msa[:10,]
    >>> msa[:10,20:40]
    >>> msa['GTHB2_ONCKE']"""
    
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
            row, col = index
        except (ValueError, TypeError):
            try:
                index = self._mapping.get(index, index)
            except TypeError:
                pass
        else:
            try:
                index = self._mapping.get(row, row), col
            except TypeError:
                pass
            
        result = self._msa[index]
        
        try:
            shape, ndim = result.shape, result.ndim
        except AttributeError:
            return result
        else:
            if ndim < 2:
                return result.tostring()
            else:
                msa = MSA(result)
                return result # return an MSA object here 
               
    def __iter__(self):
        
        for i, label in enumerate(self._labels):
            label, start, end = splitLabel(label)
            yield label, self._msa[i].tostring(), start, end
                
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
            return splitLabel(label)[0]
                
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
        ext = splitext(msa)[1] 
        if ext == '.gz' or 'filter' in kwargs or 'slice' in kwargs:
            return MSA(msa, **kwargs)
        else:
            filename = msa
    else:    
        try:
            filename = fetchPfamMSA(msa, **kwargs)
        except IOError:
            raise ValueError('msa must be an MSA filename or a Pfam accession')
        else:
            if compressed in kwargs:
                return MSA(filename, **kwargs)

    msafile = MSAFile(msa)
    title = splitext(split(msa)[1])[0]
    format = msafile.format
    lenseq = msafile.numResidues()
    numseq = getsize(msa) / (lenseq + 10)
    
    if format == FASTA:
        from .msatools import parseFasta
        msaarr = zeros((numseq, lenseq), '|S1')
        labels, mapping = parseFasta(msa, msaarr)
    elif format == SELEX or format == STOCKHOLM:
        from .msatools import parseSelex
        msaarr = zeros((numseq, lenseq), '|S1') 
        labels, mapping = parseSelex(msa, msaarr)

    return MSA(msa=msaarr[:len(labels)], labels=labels, title=title, 
               mapping=mapping)
    

def calcShannonEntropy(msa, dividend=False):
    """Return Shannon entropy array calculated for *msa*, which may be 
    an :class:`MSA` instance or a 2D Numpy character array.  The function
    is case insensitive and handles ambiguous amino acid characters as 
    follows:
    
      * **B** (Asx) count is shared by *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is shared by *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is shared by *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is shared by each of standard amino acids
      
    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.
    
    Gaps, which may be any non-alphabet characters, are handled in two ways:  
      
      * as a distinct character with its own probability, by default 
      * non-existent, the probability of observing amino acids in a given
        column is adjusted, when *dividend* is **True**."""
    
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
    calcShannonEntropy(msa, entropy, dividend=bool(dividend))
    return entropy
    
    
if __name__ == '__main__':

    #filepath = fetchPfamMSA('PF00497',alignment='seed',compressed=False,
    #                         format='stockholm',gaps='dashes')
    #print filepath
    #results = list(iterSequences(filepath))
    
    filepath1 = fetchPfamMSA('PF00007',alignment='seed',compressed=True, 
                             timeout=5)
    filepath2 = fetchPfamMSA('PF00007',alignment='seed',compressed=True, 
                             format='selex')
    filepath3 = fetchPfamMSA('PF00007',alignment='seed',compressed=False, 
                             format='fasta', outname='mymsa')
    results_sth = list(MSAFile(filepath1))
    results_slx = list(MSAFile(filepath2))
    results_fasta = list(MSAFile(filepath3))
    from Bio import AlignIO
    alignment = AlignIO.read(filepath3,'fasta')
    results_obj = list(MSAFile(alignment))
    import numpy
    numpy.testing.assert_equal(results_fasta,results_obj)
