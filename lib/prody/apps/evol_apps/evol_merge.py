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

"""Merge multiple MSAs based on common labels."""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..apptools import DevelApp

__all__ = ['evol_merge']

APP = DevelApp('merge', 'merge multiple MSAs based on common labels')

APP.setExample(
"""This application merges a number of MSAs into one large MSA. The merging of
seqences in done based on common labels appearing across all the input MSAs The
following example show how to merge two MSAs:

    $ evol merge piwi.slx -l GTHB2_ONCKE""", [])


APP.addArgument('msas',
    nargs='+',
    help='MSA filenames to be merged')

APP.addGroup('output', 'output options')
APP.addArgument('-o', '--outname',
    dest='outname',
    help='output filename, default is first input MSA filename_merged suffix',
    type=str,
    metavar='STR',
    group='output')
    
APP.addArgument('-f', '--format', 
    dest='format', 
    type=str,
    metavar='STR', 
    help='output MSA file format, default is same as first input MSA', 
    group='output')
    
APP.addArgument('-z', '--compressed', 
    dest='compressed', 
    action='store_true', 
    help='gzip merged MSA output',
    group='output')

def evol_merge(*msas, **kwargs):
    
    import prody
    from prody import parseMSA, mergeMSA, LOGGER, writeMSA
    from os.path import splitext
    if len(msas) < 2:
        raise ValueError('There has to more than one msa for merging.')
    msaobj = []
    for i, msa in enumerate(list(msas)):
        try:
            msaobj.append(parseMSA(str(msa)))
        except:
            LOGGER.info('Could not parse msa file {0}. Ignoring that input'
                        .format(msa))
        else:
            if i == 0:
                format = kwargs.get('format', None)
                if format is None:
                    title,  ext = splitext(msa)
                    if ext.lower() == '.gz':
                        title, ext = splitext(outname)
                    if ext == 'fasta':
                        format = 'FASTA'
                    elif ext == 'sth':
                        format = 'Stockholm'
                    else:
                        format = 'SELEX'
                    kwargs['format'] = format
                    outname = kwargs.get('outname')
                    if outname is None:
                        outname = title + '_merged' + ext
                        
    writeMSA(outname, mergeMSA(*msaobj), **kwargs)    
    LOGGER.info('Wrote merged file {0}'.format(outname))
    
APP.setFunction(evol_merge)

