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
    (MSA) from PFAM"""

__author__ = 'Anindita Dutta, Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Anindita Dutta, Ahmet Bakan'

__all__ = ['fetchPfamMSA', 'iterSequences']

DOWNLOAD_FORMATS = set(['seed', 'full', 'ncbi', 'metagenomics'])
FORMAT_OPTIONS = ({'format': set(['selex', 'stockholm', 'fasta']),
                  'order': set(['tree', 'alphabetical']),
                  'inserts': set(['lower', 'upper']),
                  'gaps': set(['mixed', 'dots', 'dashes'])})
import re
import os.path
from prody import LOGGER
from prody.utilities import makePath, openURL, gunzip, openFile


def fetchPfamMSA(acc, alignment='full', folder='.', compressed=False, **kwargs):
    """Returns a file path name  where the PFAM MSA has been downloaded to.
        
    :arg acc: PFAM ID or Accession Code
    :type acc: str
    
    :arg alignment: type of alignment that are available for download at PFAM, 
                    default is 'full' 
    :type alignment: str
    
    :arg folder: folder to download the MSA file to, default is '.'
    :type folder: str
    
    :arg compressed: gzip the downloaded MSA file or not, default is 'False'
    
   *Alignment Parameters*
    
    :arg format: PFAM supported MSA file formats, default is 'stockholm'
                 Other supported formats are 'selex' and 'fasta'.
    :type format: str
    
    :arg order: PFAM supported ordering of sequences, default is 'tree'
    :type order: str
    
    :arg inserts: PFAM supported case for inserts, default is 'lower'
    :type inserts: str
    
    :arg gaps: PFAM supported gap fomats, default is 'mixed'
    :type gaps: str
    
    *URL timeout*
    
    :arg timeout: timeout for url, default is 5
    :type timeout: int
    
    *Out parameter*
    
    :arg outname: out file name, default is input 'acc'
    :type outname: str"""   
   
    getAccUrl = 'http://pfam.sanger.ac.uk/family/acc?id=' + acc
    handle = openURL(getAccUrl)
    orig_acc = acc
    acc = handle.readline().strip()
    url_flag = False
    
    if not re.search('(?<=PF)[0-9]{5}$', acc):
        raise ValueError('No such family: check PFAM ID or Accession Code')
        
    
    if alignment not in DOWNLOAD_FORMATS:
        raise ValueError('alignment must be one of full, seed, ncbi or'
                         ' metagenomics')
    if(alignment == 'ncbi' or alignment == 'metagenomics'):
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
                raise ValueError('Alignment format must be of type selex'
                                 ' stockholm or fasta. MSF not supported')
            if align_format == 'selex':
                align_format, extension = 'pfam', '.slx'
            elif align_format == 'fasta':
                extension = '.fasta'
            else:
                extension = '.sth'
            gaps = kwargs.get('gaps', 'mixed').lower()
            if gaps not in FORMAT_OPTIONS['gaps']:
                raise ValueError('Gaps must be of type mixed, dots or dashes')
            inserts = kwargs.get('inserts', 'lower').lower()
            if(inserts not in FORMAT_OPTIONS['inserts']):
                raise ValueError('inserts must be of type lower or upper')
            order = kwargs.get('order', 'tree').lower()
            if order not in FORMAT_OPTIONS['order']:
                raise ValueError('Order must be of type tree or alphabetical')
            url = ('http://pfam.sanger.ac.uk/family/' + acc + '/alignment/'
                   + alignment + '/format?format=' + align_format + 
                   '&alnType=' + alignment + '&order=' + order[0] +
                   '&case=' + inserts[0] + '&gaps=' + gaps + '&download=1')
        
    
    response =  openURL(url, timeout = int(kwargs.get('timeout', 5))) 
    outname = kwargs.get('outname', None)
    if not outname:
        outname = orig_acc
    filepath = os.path.join(makePath(folder), (outname + '_' + alignment 
                            + extension))
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
                
    LOGGER.info('Pfam MSA for {0} is written as {1}'.format(orig_acc, filepath))
    
    return filepath

def parseTitle(title):
    """Parses PFAM title associated with each sequence"""
    
    split_title = re.split('/*-*', title)
    if split_title:
        if (len(split_title) == 3):
            id, start, end = split_title[0], split_title[1], split_title[2]
            if start.isdigit(): 
                start = int(start)
            else:
                Logger.warn('starting residue number could not be parsed from'
                            ' {0}, returning none'.format(id))
                start = None
            if end.isdigit():
                end = int(end)
            else:
                Logger.warn('ending residue number could not be parsed from'
                            ' {0}, returning none'.format(id))
                end = None
        else:
            LOGGER.warn('Missing start or end or both for'
                        ' {0}. Appending id.'.format(id))
            id, start, end =  split_title[0], None, None 
    else:
        raise ValueError('Cannot parse MSA ID. Check MSA file!')
    return (id, start, end)
    


def iterSequences(msa):
    """Returns a list of tuples containing sequence id, sequence,
    start index, end index, from an MSA file or object
    
    :arg msa: MSA filepath or Biopython MSA object
    :type msa: str or object"""

    if os.path.isfile(str(msa)):
        msa_file = openFile(msa)
        firstline = msa_file.readline()
        if firstline[0] == '>':
            LOGGER.info('Parsing fasta format')
            id, start, end = parseTitle(firstline[1:].strip())
            temp =[]
            for line in msa_file:
                if line[0] == '>':
                    yield (id, ''.join(temp), start, end)
                    temp = []
                    id, start, end = parseTitle(line[1:].strip())
                else:
                    temp.append(line.strip())
            yield (id, ''.join(temp), start, end)
        elif firstline[0] == '#':
            if firstline.split()[1] == 'STOCKHOLM':
                LOGGER.info('Parsing stockholm format')
            else:
                LOGGER.info('Parsing selex format')
            for line in msa_file:
                if(line[0] != '#' and line[:2] != '//'):
                    if(len(line.split()) == 2):
                        id, start, end = parseTitle(line.split()[0].strip())
                        yield (id, line.split()[1].strip(), start, end)
                    else:
                        raise ValueError('Possible incorrect file format.'
                                         ' Cannot parse file.')              
        else:
            LOGGER.info('Parsing selex format')
            if(len(firstline.split()) == 2):
                id, start, end = parseTitle(firstline.split()[0].strip())
                yield (id, firstline.split()[1].strip(), start, end)
                for line in msa_file:
                    if(line[0] != '#'):
                        if(len(line.split()) == 2):
                            id, start, end = parseTitle(line.split()[0].strip())
                            yield (id, line.split()[1].strip(), start, end)
                        else:
                            raise ValueError('Some lines in file are not the '
                                             'right format. Cannot parse file')
            else:
                raise ValueError('Possible incorrect file format.'
                                ' Cannot parse file.')
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

    #filepath = fetchPfamMSA('PF00497',alignment='seed',compressed=False,format='stockholm',gaps='dashes')
    #print filepath
    #results = list(iterSequences(filepath))
    
    filepath1 = fetchPfamMSA('PF00007',alignment='seed',compressed=True,timeout=5)
    filepath2 = fetchPfamMSA('PF00007',alignment='seed',compressed=True, format='selex')
    filepath3 = fetchPfamMSA('PF00007',alignment='seed',compressed=False, format='fasta', outname='mymsa')
    results_sth = list(iterSequences(filepath1))
    results_slx = list(iterSequences(filepath2))
    results_fasta = list(iterSequences(filepath3))
    from Bio import AlignIO
    alignment = AlignIO.read(filepath3,'fasta')
    results_obj = list(iterSequences(alignment))
    numpy.testing.assert_equal(results_fasta,results_obj)
    

    
        

