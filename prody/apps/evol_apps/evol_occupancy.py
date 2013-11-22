"""MSA residue coevolution calculation application."""

__author__ = 'Anindita Dutta, Ahmet Bakan'

from ..apptools import DevelApp

__all__ = ['evol_occupancy']

APP = DevelApp('occupancy',
               help='calculate occupancy of rows and columns in MSA')

APP.setExample(
"""Calculate occupancy of rows and columns of the msa. Any msa file in
supported formats 'FASTA', 'SELEX', 'Stockholm' can be given as input.
These plots can be used to identify whether MSA needs refinement to reduce
the number of gaps in order to avoid false positives in later analysis.

Following example will save row occupancy data and plot using default options:

  $ evol occupancy piwi_refined.slx

Following example will save both row and column occupancy and plot using
specified figure parameters:

  $ evol occupancy piwi_refined.slx -o both -X Residue -Y Occupancy""")


APP.addArgument('msa',
    help='MSA file')

APP.addGroup('calc', 'calculation options')
APP.addArgument('-o', '--occ-axis',
    dest='occaxis',
    help='calculate row or column occupancy or both.',
    choices=['row', 'col', 'both'],
    metavar='STR',
    type=str,
    default='row',
    group='calc')

APP.addGroup('output', 'output options')

APP.addArgument('-p', '--prefix',
    dest='prefix',
    help='output filename prefix, default is '
         'msa filename with _occupancy suffix',
    type=str,
    metavar='STR',
    group='output')

APP.addArgument('-l', '--label',
    dest='label',
    help='index for column based on msa label',
    type=str,
    metavar='STR',
    group='output')

APP.addArgument('-f', '--number-format',
    dest='numformat', type=str, default='%12g',
    metavar='STR', help='number output format', group='output')

APP.addFigarg('-X', '--xlabel',
    dest='xlabel',
    help='specify xlabel',
    type=str,
    metavar='STR',
    default=None)
APP.addFigarg('-Y', '--ylabel',
    dest='ylabel',
    help='specify ylabel',
    type=str,
    metavar='STR',
    default=None)
APP.addFigarg('-T', '--title',
    dest='title',
    help='figure title',
    type=str,
    metavar='STR',
    default=None)
APP.addFigure('-S', '--save-plot',
    dest='figocc',
    action='store_true',
    help='save occupancy plot/s')


def evol_occupancy(msa, **kwargs):

    from numpy import arange

    import prody
    from prody import parseMSA, calcMSAOccupancy, showMSAOccupancy, writeArray
    from os.path import splitext

    prefix = kwargs.get('prefix')
    if prefix is None:
        prefix, _ = splitext(msa)
        if _.lower() == '.gz':
            prefix, _ = splitext(prefix)
        prefix += '_occupancy'

    msa = parseMSA(msa)

    numformat = kwargs.get('numformat', '%12g')
    occupancy , suffix = [], []
    occaxis = kwargs.get('occaxis', 'row')
    if occaxis == 'both':
        suffix = ['_row', '_col']
        occupancy.append(calcMSAOccupancy(msa, occ='row'))
        occupancy.append(calcMSAOccupancy(msa, occ='col'))
    else:
        suffix = '_' + occaxis
        occupancy.append(calcMSAOccupancy(msa, occ=occaxis))

    for i, occ in enumerate(occupancy):
        writeArray((prefix + suffix[i] + '.txt'), occ, format=numformat)

    for i, occ in enumerate(occupancy):
        if kwargs.get('figocc'):
            try:
                import matplotlib.pyplot as plt
            except ImportError:
                LOGGER.warn('Matplotlib could not be imported, '
                            'figures are not saved.')
            else:
                prody.SETTINGS['auto_show'] = False
                width = kwargs.get('figwidth', 8)
                height = kwargs.get('figheight', 6)
                xlabel = kwargs.get('xlabel')
                title = kwargs.get('title')
                figure = plt.figure(figsize=(width, height))
                label = kwargs.get('label')
                show = showMSAOccupancy(msa=msa, occ=occ, label=label,
                                         xlabel=xlabel, title=title)
                format = kwargs.get('figformat', 'pdf')
                figure.savefig(prefix + suffix[i] + '.' + format, format=format,
                            dpi=kwargs.get('figdpi', 300))


APP.setFunction(evol_occupancy)



