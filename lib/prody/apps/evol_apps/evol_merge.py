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


APP.addArgument('msa',
    nargs='+',
    help='MSA filenames to be merged')

APP.addGroup('output', 'output options')
APP.addArgument('-o', '--outname',
    dest='outname',
    help='output filename, default is first input filename '
         'with _merged suffix',
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

def evol_merge(*msa, **kwargs):
    
    import prody
    from prody import parseMSA, mergeMSA, LOGGER, writeMSA, MSAFile
    from prody.sequence.msafile import MSAEXTMAP
    from os.path import splitext
    if len(msa) < 2:
        raise ValueError('multiple msa filenames must be specified')
    msaobj = []
    try:
        msaobj = [parseMSA(fn) for fn in msa]
    except:
        raise IOError('failed to parse {0}'.format(fn))
    
    msafile = MSAFile(msa[0])

    format = kwargs.get('format') or msafile.format
    outname = kwargs.get('outname') or (msafile.getTitle() + '_merged' + 
                                        MSAEXTMAP[msafile.format])
    writeMSA(outname, mergeMSA(*msaobj), **kwargs)    
    LOGGER.info('Merged MSA is saved as: {0}'.format(outname))
    
APP.setFunction(evol_merge)

