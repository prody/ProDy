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

__all__ = ['fetchPfamMSA', 'MSAFile', 'MSA']

DOWNLOAD_FORMATS = set(['seed', 'full', 'ncbi', 'metagenomics'])
FORMAT_OPTIONS = ({'format': set(['selex', 'stockholm', 'fasta']),
                  'order': set(['tree', 'alphabetical']),
                  'inserts': set(['lower', 'upper']),
                  'gaps': set(['mixed', 'dots', 'dashes', 'none'])})
import re
from os.path import join, isfile

from numpy import zeros

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
            
            if align_format == 'selex':
                align_format, extension = 'pfam', '.slx'
            elif align_format == 'fasta':
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

def parseTitle(title):
    """Parses Pfam title associated with each sequence"""
    
    split_title = re.split('/*-*', title)
    if split_title:
        if (len(split_title) == 3):
            id, start, end = split_title[0], split_title[1], split_title[2]
            if start.isdigit(): 
                start = int(start)
            else:
                Logger.warn('starting residue number could not be parsed from'
                            ' {0}'.format(id))
                start = None
            if end.isdigit():
                end = int(end)
            else:
                Logger.warn('ending residue number could not be parsed from'
                            ' {0}'.format(id))
                end = None
        else:
            LOGGER.warn('missing start and/or end indicies {0}'.format(id))
            id, start, end =  split_title[0], None, None
    else:
        raise ValueError('Cannot parse MSA ID. Check MSA file!')
    return (id, start, end)
    


class MSAFile(object):
    
    """Yield tuples containing sequence id, sequence, start index, end index, 
    from an MSA file or object."""
    
    def __init__(self, msa, aligned=True):
        """*msa* may be an MSA file in fasta, Stockholm, or Selex format or a
        Biopython MSA object.  If *aligned* is **True**, length of sequences
        will be expected to be the same, and unaligned sequences will cause
        an :exc:`IOError` exception."""
        
        self._msa = None
        self._bio = None
        self._aligned = bool(aligned)
        self._lenseq = None
        self._numseq = None
        if isfile(str(msa)):    
            self._msa = msa
            with openFile(msa) as msa: 
                line = msa.readline()
                if line[0] == '>':
                    LOGGER.info('Parsing MSA in fasta format')
                    self._iter = self._iterFasta
                elif line[0] == '#' and 'STOCKHOLM' in line:
                    LOGGER.info('Parsing MSA in Stockholm format')
                    self._iter = self._iterStockholm
                else:
                    LOGGER.info('Parsing MSA in Selex format')
                    self._iter = self._iterStockholm
        else:    
            try:
                bio = iter(msa)
            except TypeError:
                raise TypeError('msa must be a filename or Bio object')
            else:
                record = bio.next()
                try:
                    title, seq = record.id, str(record.seq)
                except AttributeError:
                    raise TypeError('msa must be a filename or Bio object')
                else:
                    self._bio = iter(msa)
                    self._iter = self._iterBio

       
    def __del__(self):
        
        if self._msa:
            self._msa.close()
            
    def __iter__(self):
        
        return self._iter()
            
    def _iterBio(self):
        """Yield sequences from a Biopython MSA object."""            
        
        numseq = 0
        for record in self._bio:
            title, seq = record.id, str(record.seq)
            id, start, end = parseTitle(title)
            numseq += 1
            yield (id, seq, start, end)
        self._numseq = numseq

    def _iterFasta(self):
        """Yield sequences from an MSA file in fasta format."""

        aligned = self._aligned
        lenseq = self._lenseq
        temp = []
        numseq = 0
        with openFile(self._msa) as msa: 
            title = msa.readline()[1:]
            id, start, end = parseTitle(title.strip())
            for line in msa:
                if line[0] == '>':
                    seq = ''.join(temp)
                    if aligned:
                        if lenseq:
                            if lenseq != len(seq):
                                raise IOError('sequence for {0:s} does not '
                                              'have expected length {1:d}'
                                              .format(title, self._lenseq))
                        else: 
                            self._lenseq = lenseq = len(seq)
                    numseq += 1
                    yield (id, seq, start, end)
                    temp = []
                    title = line[1:].strip()
                    id, start, end = parseTitle(title)
                else:
                    temp.append(line.strip())
            numseq += 1
            yield (id, ''.join(temp), start, end)
        self._numseq = numseq
    
    def _iterStockholm(self):
        """Yield sequences from an MSA file in Stockholm format."""

        aligned = self._aligned
        lenseq = self._lenseq
        numseq = 0

        with openFile(self._msa) as msa:
            for line in msa:
                if line[0] == '#' or line[0] == '/':
                    continue
                items = line.split()
                title = ' '.join(items[:-1])
                seq = items[-1]
                if aligned:
                    if lenseq:
                        if lenseq != len(seq):
                            raise IOError('sequence for {0:s} does not have '
                                          'expected length {1:d}'
                                          .format(title, self._lenseq))
                    else: 
                        self._lenseq = lenseq = len(seq)

                id, start, end = parseTitle(title)
                numseq += 1
                yield (id, seq, start, end)
        self._numseq = numseq
        
    def numSequences(self):
        """Return number of sequences."""
        
        if self._numseq is None:
            for i in self:
                pass
        return self._numseq
        
    def numResidues(self):
        """Return number of residues (or columns in the MSA)."""
        
        if self._lenseq is None:
            iter(self).next()
        return self._lenseq
    
    
class MSA(object):
    
    """Store and manipulate multiple sequence alignments.
    
    >>> msa = parseMSA('piwi', alignment='seed')
    >>> msa[0]
    >>> msa[0,0]
    >>> msa[:10,]
    >>> msa[:10,20:40]
    >>> m['GTHB2_ONCKE']"""
    
    def __init__(self, msa):
        """*msa* may be an :class:`MSAFile` instance or an MSA file in a 
        supported format."""
        
        
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
        self._titles = []
        titles = self._titles.append
        self._resnums = []
        resnums = self._resnums.append
        self._tmap = tmap = {}
        
        for i, (title, seq, start, end) in enumerate(msa):
            resnums((start, end))
            titles(title)
            tmap[title] = i
            msaarray[i] = list(seq)
        
    def __getitem__(self, index):
        
        try:
            row, col = index
        except (ValueError, TypeError):
            try:
                index = self._tmap.get(index, index)
            except TypeError:
                pass
        else:
            try:
                index = self._tmap.get(row, row), col
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
                return result # return an MSA object here 
                
    def numSequences(self):
        """Return number of sequences."""
        
        return self._msa.shape[0]

    def numResidues(self):
        """Return number of residues (or columns in the MSA)."""
        
        return self._msa.shape[1]
    
    def getTitle(self, index):
        pass
        
    def getResnums(self, index):
        pass

    def getArray(self):
        """Return a copy of the MSA character array."""
        
        return self._msa.copy
    
    
def parseMSA(msa, **kwargs):
    """Should return an :class:`MSA`. A Pfam MSA sould be fetched if needed."""
    
    msa = str(msa)
    if isfile(msa):
        if 'filter' in kwargs or splitext(msa)[1] == '.gz':
            return MSA(msa)
        else:
            pass
    

    
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
