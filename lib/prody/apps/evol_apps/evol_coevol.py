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

"""MSA residue coevolution calculation application."""

__author__ = 'Anindita Dutta, Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from ..apptools import DevelApp

__all__ = ['evol_coevol']

APP = DevelApp('coevol', 'calculate coevolution using mutual information')

APP.setExample(
"""This application calculates mutual information between MSA postions for a \
refined multiple sequence alignment.  Following example will save coevolution \
data and plot using default options:

    $ evol coevol piwi_refined.slx -S""", [])


APP.addArgument('msa', 
    help='refined MSA file')

APP.addGroup('calc', 'calculation options')
APP.addArgument('-n', '--no-ambiguity',
    dest='ambiguity',
    help='treat amino acids characters B, Z, J, and X as non-ambiguous',
    default=True,
    action='store_false',
    group='calc')

APP.addArgument('-g', '--gaps',
    dest='omitgaps',
    help='do not omit gap characters',
    default=True,
    action='store_false',
    group='calc')

APP.addArgument('-c', '--correction',
    dest='correction',
    help='correction to apply to mutual information matrix',
    choices=['apc', 'asc'],
    metavar='STR',
    type=str,
    group='calc')

APP.addGroup('output', 'output options')
APP.addArgument('-p', '--prefix',
    dest='prefix',
    help='output filename prefix, default is '
         'msa filename with _coevol suffix',
    type=str,
    metavar='STR',
    group='output')
    
APP.addArgument('-f', '--number-format', 
    dest='numformat', type=str, default='%12g', 
    metavar='STR', help='number output format', group='output')

APP.addFigarg('-L', '--cmin',
    dest='cmin',
    help='apply lower limits for figure plot',
    type=float,
    metavar='FLOAT')

APP.addFigarg('-U', '--cmax',
    dest='cmax',
    help='apply upper limits for figure plot',
    type=float,
    metavar='FLOAT')
        
APP.addFigure('-S', '--save-plot', 
    dest='figcoevol', 
    action='store_true', 
    help='save coevolution plot')


def evol_coevol(msa, **kwargs):
    
    import prody
    from prody import parseMSA, buildMutinfoMatrix, showMutualInfo
    from prody import applyMICorrection
    from prody import writeArray
    from os.path import splitext

    prefix = kwargs.get('prefix')
    if prefix is None:
        prefix, _ = splitext(msa)
        if _.lower() == '.gz':
            prefix, _ = splitext(prefix)
        prefix += '_coevol'
    msa = parseMSA(msa)
    mutinfo = buildMutinfoMatrix(msa, **kwargs)
    if kwargs.get('correction') == 'apc':
        mutinfo = applyMICorrection(mutinfo, correction='prod')
    if kwargs.get('correction') == 'asc':
        mutinfo = applyMICorrection(mutinfo, correction='sum')
    writeArray(prefix + '.txt', 
               mutinfo, format=kwargs.get('numformat', '%12g'))
    if kwargs.get('figcoevol'):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warn('Matplotlib could not be imported, '
                        'figures are not saved.')
        else:
            cmin = kwargs.get('cmin', mutinfo.min())
            cmax = kwargs.get('cmax', mutinfo.max())
            prody.SETTINGS['auto_show'] = False
            width = kwargs.get('figwidth', 8)
            height = kwargs.get('figheight', 6)
            figargs = kwargs.get('figargs', ())
            figure = plt.figure(figsize=(width, height))
            show = showMutualInfo(mutinfo, msa=msa, clim=(cmin, cmax), *figargs)
            format = kwargs.get('figformat', 'pdf')
            figure.savefig(prefix + '.' + format, format=format,
                        dpi=kwargs.get('figdpi', 300))         
    

APP.setFunction(evol_coevol)

