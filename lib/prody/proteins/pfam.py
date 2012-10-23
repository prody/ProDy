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

__all__ = ['fetchPfamMSA', 'iterSequences']

DOWNLOAD_FORMATS = set(['seed', 'full', 'ncbi', 'metagenomics'])
FORMAT_OPTIONS = ({'format': set(['selex', 'stockholm', 'fasta']),
                  'order': set(['tree', 'alphabetical']),
                  'inserts': set(['lower', 'upper']),
                  'gaps': set(['mixed', 'dots', 'dashes', 'none'])})
import re
from os.path import join, isfile
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
    
    :arg format: a Pfam supported MSA file format, one of ``'stockholm'`` 
        (default), ``'selex'``, or ``'fasta'``
    
    :arg order: ordering of sequences, ``'tree'`` (default) or 
        ``'alphabetical'`` 
    
    :arg inserts: letter case for inserts, ``'lower'`` (default) or ``'upper'``
    
    :arg gaps: gap character, one of ``'mixed'`` (default), ``'dots'``, 
        ``'dashes'`` or **None** for unaligned
    
    
    *Out parameter*
    
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
            align_format = kwargs.get('format', 'stockholm').lower()
            
            if align_format not in FORMAT_OPTIONS['format']:
                raise ValueError('alignment format must be of type selex'
                                 ' stockholm or fasta. MSF not supported')
            
            if align_format == 'selex':
                align_format, extension = 'pfam', '.slx'
            elif align_format == 'fasta':
                extension = '.fasta'
            else:
                extension = '.sth'
            
            gaps = str(kwargs.get('gaps', 'mixed')).lower()
            if gaps not in FORMAT_OPTIONS['gaps']:
                raise ValueError('gaps must be of type mixed, dots, dashes, '
                                 'or None')
            
            inserts = kwargs.get('inserts', 'lower').lower()
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
    


def iterSequences(msa):
    """Yield tuples containing sequence id, sequence, start index, end index, 
    from an MSA file or object.
    
    :arg msa: MSA filepath or Biopython MSA object"""

    if isfile(str(msa)):
        msa_file = openFile(msa)
        firstline = msa_file.readline()
        if firstline[0] == '>':
            LOGGER.info('Parsing fasta file')
            id, start, end = parseTitle(firstline[1:].strip())
            temp = []
            for line in msa_file:
                if line[0] == '>':
                    yield (id, ''.join(temp), start, end)
                    temp = []
                    id, start, end = parseTitle(line[1:].strip())
                else:
                    temp.append(line.strip())
            yield (id, ''.join(temp), start, end)
    
        elif firstline[0] == '#' and 'STOCKHOLM' in firstline:
            LOGGER.info('Parsing stockholm file')
            for line in msa_file:
                if line[0] != '#' and line[:2] != '//':
                    items = line.split()
                    if len(items) == 2:
                        id, start, end = parseTitle(items[0])
                        yield (id, items[1], start, end)

        else:
            items = firstline.split()
            if len(items) != 2:
                raise ValueError('file format could not be recognized')
            LOGGER.info('Parsing selex format')
            
            id, start, end = parseTitle(items[0])
            yield (id, items[1], start, end)
            for line in msa_file:
                if line[0] == '#':
                    continue
                items = line.split() 
                if len(items) == 2:
                    id, start, end = parseTitle(items[0])
                    yield (id, items[1], start, end)
        msa_file.close()

    else:
        try:
            msa = iter(msa)
        except TypeError:
            raise TypeError('msa must be a filename or Bio object')
        else:
            record = msa.next()
            try:
                title, seq = record.id, str(record.seq)
            except AttributeError:
                raise TypeError('msa must be a filename or Bio object')
            
            id, start, end = parseTitle(title)
            yield (id, seq, start, end)
            for record in msa:
                title, seq = record.id, str(record.seq)
                id, start, end = parseTitle(title)     
                yield (id, seq, start, end)

    
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
    results_sth = list(iterSequences(filepath1))
    results_slx = list(iterSequences(filepath2))
    results_fasta = list(iterSequences(filepath3))
    from Bio import AlignIO
    alignment = AlignIO.read(filepath3,'fasta')
    results_obj = list(iterSequences(alignment))
    import numpy
    numpy.testing.assert_equal(results_fasta,results_obj)
