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

"""This module defines functions for interfacing Pfam database."""

__author__ = 'Anindita Dutta, Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Anindita Dutta, Ahmet Bakan'

import re
from os.path import join, isfile, getsize, splitext, split
import xml.etree.cElementTree as ET
from collections import defaultdict

from prody import LOGGER
from prody.utilities import makePath, openURL, gunzip, openFile
from prody.proteins import MSAFile


__all__ = ['searchPfam', 'fetchPfamMSA']

FASTA = 'fasta'
SELEX = 'selex'
STOCKHOLM = 'stockholm'

DOWNLOAD_FORMATS = set(['seed', 'full', 'ncbi', 'metagenomics'])
FORMAT_OPTIONS = ({'format': set([FASTA, SELEX, STOCKHOLM]),
                  'order': set(['tree', 'alphabetical']),
                  'inserts': set(['lower', 'upper']),
                  'gaps': set(['mixed', 'dots', 'dashes', 'none'])})



SPLITACCESSION = re.compile('\.*').split 

def searchPfam(seq, timeout=15, delay=5, **kwargs):
    
    """Returns a dictionary which contains Pfam accession ID as keys
    and evalue, alignment start and end properties as values.
        
    :arg seq: UniProt ID, protein sequence, or a protein sequence file
    :type seq: str
    
    :arg timeout: timeout for blocking connection attempt in seconds, default 
                  is 15
    :type timeout: int
    
    :arg delay: a delay required for sequence search on Pfam, default 
                  is 5, may use larger for longer sequences
    :type delay: int
    
    *Sequence Search Options*
    
    :arg ga: Gathering threshold value, either 1 (default) or 0
    :type ga: int
    
    :arg evalue: user specified evalue, must be smaller than < 10.0
    :type evalue: float 
    
    :arg searchBs: Search pfamB families when **True**
    :type searchBs: bool
    
    :arg skipAs: Do not search pfamA families when **True**
    :type skipAs: bool"""
    
    prefix = '{http://pfam.sanger.ac.uk/}'
    if isfile(str(seq)):
        try:
            seq_list = list(MSAFile(str(seq)))
        except:
            LOGGER.warn('MSAFile could not parse given file')
            seq_list = []        
        if not seq_list:
            try:
                with openFile(seq) as seq: 
                    seq = seq.read()
                    seq = ''.join(seq.split())
            except:
                raise ValueError('Could not open/read input file')
        else:
            if len(seq_list) > 1:
                LOGGER.warn('More than one sequence. Using first sequence')
            seq = seq_list[0][1]
    else:
        try:
            seq = ''.join(str(seq).split())   
        except:
            raise TypeError('sequence must be a string or file')
    if len(seq) > 10:
        if not seq.isalpha():
            raise ValueError(repr(seq) + ' is not a valid sequence')
        urlextension = ''
        if kwargs:
            ga = int(kwargs.get('ga', 1))
            if not (ga == 1 or ga == 0):
                raise ValueError('Must be either 0 or 1')
            evalue = kwargs.get('evalue', None)
            if evalue:
                if not float(evalue) <= 10.0:
                    raise ValueError('Must be a valid float < 10.0')
                LOGGER.info('Using evalue, ignoring ga value given, if any')
                urlextension = urlextension + '&evalue=' + str(evalue)
            else:
                urlextension = urlextension + '&ga=' + str(ga)            
            searchBs = int(bool(kwargs.get('searchBs', False)))          
            skipAs = int(bool(kwargs.get('skipAs', False)))
            if skipAs == 1:
                LOGGER.info('SkipAs is True. Setting searchBs to True.')
                searchBs = 1
            urlextension = urlextension + '&searchBs=' + str(searchBs)
            urlextension = urlextension + '&skipAs=' + str(skipAs)
        url1 = ('http://pfam.sanger.ac.uk/search/sequence?seq=' + str(seq)
               + urlextension + '&output=xml')
        url = ''
        try: 
            root = ET.XML(openURL(url1, timeout=int(timeout)).read())   
            for item in list(root.iter()):
                if item.tag[len(prefix):] == 'result_url':
                    url = item.text
            import time
            time.sleep(int(delay))    
        except:
            raise ValueError('Could not find result url for input sequence')  
        if not url:
            raise ValueError('Could not find result url for input sequence')        
    else:   
        url  = 'http://pfam.sanger.ac.uk/protein/' + seq + '?output=xml'
    
    try:
        root =  ET.XML(openURL(url).read())
    except:
        raise ValueError('Could not parse url output. Check input arguments'
                         ' or set larger delay for bigger sequence searches')
    else:
        xml_result = list(root.iter())
        matches = defaultdict(list)
        for i in range(len(xml_result)):
            if xml_result[i].tag[len(prefix):] == 'matches':
                result = xml_result[i]
                if len(result) == 0:
                    raise ValueError('No Pfam matches found')
                children = list(result.getchildren())
                if children[0].tag[len(prefix):] == 'protein':
                    LOGGER.info('Parsing sequence matches')
                    subchild = list(children[0].getchildren())
                    if (subchild[0].tag[len(prefix):] == 'database' 
                        and len(subchild) == 1):
                        result = subchild[0]
                    else:
                        raise ValueError('Found matches but cannot parse XML')
                elif children[0].tag[len(prefix):] == 'match':
                    LOGGER.info('Parsing id matches')
                else:
                    raise ValueError('Found matches but cannot parse XML')
                if len(result) > 1:
                    LOGGER.info('More than one match found for given input.')
                for child in result:
                    if child.tag[len(prefix):] == 'match':
                        acc = [item for item in child.items()
                                     if (item[0] == 'accession')]
                        if not acc:
                            raise ValueError('XML error. No accession found.')
                        accession = acc[0][1]
                        accession = SPLITACCESSION(accession)[0]
                        if not re.search('(?<=PF)[0-9]{5}$', accession):
                            raise ValueError('{0} does not match pfam accession'
                                             ' format!'.format(accession))
                        if len(child) > 1:
                            LOGGER.info('More than one location match found for'
                                        ' Pfam family {0}'.format(accession))
                        for subchild in child:
                            if(matches.has_key(accession)):
                                LOGGER.info('Appending multiple matches for same'
                                            ' pfam accession: {0:s}'
                                            .format(accession))
                            locations = sorted(subchild.items(),
                                               key=lambda x: x[0], reverse=True)    
                            matches[accession].append([item for
                                                        item in locations  if 
                                                        (item[0] == 'ali_start' 
                                                        or item[0] == 'ali_end' 
                                                        or item[0] == 'evalue')])
    if not matches:
        LOGGER.info('Found no matches for given sequence. '
                    'Returning empty dictionary')
    return matches

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


if __name__ == '__main__':

    from prody import *
    #filepath = fetchPfamMSA('PF00497',alignment='seed',compressed=False,
    #                         format='stockholm',gaps='dashes')
    #print filepath
    #results = list(iterSequences(filepath))
    
    #filepath1 = fetchPfamMSA('PF00007',alignment='seed',compressed=True, 
    #                         timeout=5)
    #filepath2 = fetchPfamMSA('PF00007',alignment='seed',compressed=True, 
    #                         format='selex')
    #filepath3 = fetchPfamMSA('PF00007',alignment='seed',compressed=False, 
    #                         format='fasta', outname='mymsa')
    #results_sth = list(MSAFile(filepath1))
    #results_slx = list(MSAFile(filepath2))
    #results_fasta = list(MSAFile(filepath3))
    #from Bio import AlignIO
    #alignment = AlignIO.read(filepath3,'fasta')
    #results_obj = list(MSAFile(alignment))
    #import numpy
    #numpy.testing.assert_equal(results_fasta,results_obj)

    matches1 = searchPfam('P12821')
    #matches2 = searchPfam('test.seq')
    matches3 = searchPfam('PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGSSEPCALCSLHSIGKIGGAQNRSYSKLLC\
GLLAERLRISPDRVYINYYDMNAANVGWNNSTFA', evalue=2, skipAs=True)
    matches3 = searchPfam('NSIQIGGLFPRGADQEYSAFRVGMVQFSTSEFRLTPHIDNLEVANSFAVTNAFCSQFSRGVYAIFGFYDKKSVNTITSFC\
GTLHVSFITPSFPTDGTHPFVIQMRPDLKGALLSLIEYYQWDKFAYLYDSDRGLSTLQAVLDSAAEKKWQVTAINVGNINNDKKDETYRSLFQDLELKKERRVILDCERDKVNDIVDQVITIGKHVKGYHYIIANLGFTDGDLLKIQFGGAEVSGFQIVD\
YDDSLVSKFIERWSTLEEKEYPGAHTATIKYTSALTYDAVQVMTEAFRNLRKQRIEISRRGNAGDCLANPAVPWGQGVEI\
ERALKQVQVEGLSGNIKFDQNGKRINYTINIMELKTNGPRKIGYWSEVDKMVLTEDDTSGLEQKTVVVTTILESPYVMMK\
ANHAALAGNERYEGYCVDLAAEIAKHCGFKYKLTIVGDGKYGARDADTKIWNGMVGELVYGKADIAIAPLTITLVREEVI\
DFSKPFMSLGISIMIKKPQKSKPGVFSFLDPLAYEIWMCIVFAYIGVSVVLFLVSRFSPYEWHTEEFEDGRETQSSESTN\
EFGIFNSLWFSLGAFMQQGADISPRSLSGRIVGGVWWFFTLIIISSYTANLAAFLTVERMVSPIESAEDLSKQTEIAYGT\
LDSGSTKEFFRRSKIAVFDKMWTYMRSAEPSVFVRTTAEGVARVRKSKGKYAYLLESTMNEYIEQRKPCDTMKVGGNLDS\
KGYGIATPKGSSLGTPVNLAVLKLSEQGLLDKLKNKWWYDKGECGAKDSGSKEKTSALSLSNVAGVFYILVGGLGLAMLV\
ALIEFCYKSRAEAKRMKGLVPRG', delay=10, evalue=2, searchBs=True)  