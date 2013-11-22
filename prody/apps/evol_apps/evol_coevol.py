"""MSA residue coevolution calculation application."""

__author__ = 'Anindita Dutta, Ahmet Bakan'

from ..apptools import DevelApp

__all__ = ['evol_coevol']

APP = DevelApp('coevol',
               help='analyze co-evolution using mutual information')

APP.setExample(
"""Analyze coevolution by performing mutual information calculation between
MSA positions.  A refined MSA without gaps should be used.

Following example will save coevolution data and plot using default options:

  $ evol coevol piwi_refined.slx -S

Following example will save coevolution data and plot for all correction and \
normalizations:

  $ evol coevol piwi_refined.slx -S -c apc -c asc -m sument -m minent \
-m maxent -m mincon -m maxcon -m joint""", [])


APP.addArgument('msa',
    help='refined MSA file')

APP.addGroup('calc', 'calculation options')
APP.addArgument('-n', '--no-ambiguity',
    dest='ambiguity',
    help='treat amino acids characters B, Z, J, and X as non-ambiguous',
    default=True,
    action='store_false',
    group='calc')

APP.addArgument('-c', '--correction',
    dest='correction',
    help='also save corrected mutual information matrix data and plot',
    choices=['apc', 'asc'],
    metavar='STR',
    type=str,
    action='append',
    group='calc')

APP.addArgument('-m', '--normalization',
    dest='normalization',
    help='also save normalized mutual information matrix data and plot',
    choices='sument minent maxent mincon maxcon joint'.split(),
    metavar='STR',
    type=str,
    action='append',
    group='calc')

APP.addGroup('output', 'output options')

APP.addArgument('-t', '--heatmap',
    dest='heatmap',
    help='save heatmap files for all mutual information matrices',
    default=False,
    action='store_true',
    group='output')

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
APP.addFigarg('-X', '--xlabel',
    dest='xlabel',
    help='specify xlabel, by default will be applied on ylabel',
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
    dest='figcoevol',
    action='store_true',
    help='save coevolution plot')


def evol_coevol(msa, **kwargs):

    from numpy import arange

    import prody
    from prody import parseMSA, buildMutinfoMatrix, showMutinfoMatrix
    from prody import applyMutinfoCorr, calcShannonEntropy
    from prody import writeArray, LOGGER, applyMutinfoNorm, writeHeatmap
    from os.path import splitext

    prefix = kwargs.get('prefix')
    if prefix is None:
        prefix, _ = splitext(msa)
        if _.lower() == '.gz':
            prefix, _ = splitext(prefix)
        prefix += '_mutinfo'

    msa = parseMSA(msa)
    mutinfo = buildMutinfoMatrix(msa, **kwargs)
    numformat = kwargs.get('numformat', '%12g')
    heatmap = kwargs.get('heatmap', False)
    #writeArray(prefix + '.txt', mutinfo, format=numformat)
    if heatmap:
        hmargs = {
                  'xlabel': 'Residue', 'ylabel': 'Residue',
                  'xorigin': 1, 'xstep': 1,
                  'residue': arange(msa.numResidues())}

    todo = [(None, None)]
    norm = kwargs.get('normalization', [])
    corr = kwargs.get('correction', [])
    if norm is not None:
        if 'joint' in norm:
            todo.append(('norm', 'joint'))
        for which in norm:
            if which == 'join': continue
            todo.append(('norm', which))
    if corr is not None:
        for which in corr:
            todo.append(('corr', which))
    entropy = None

    for what, which in todo:
        if what is None:
            matrix = mutinfo
            suffix = ''
            tuffix = ' Mutual Information'
        elif which == 'joint':
            LOGGER.info('Applying {0} normalization.'.format(repr(which)))
            matrix = buildMutinfoMatrix(msa, norm=True, **kwargs)
            suffix = '_norm_joint'
            tuffix = ' MI - Normalization: ' + which
        elif what == 'norm':
            LOGGER.info('Applying {0} normalization.'.format(repr(which)))
            if entropy is None:
                entropy = calcShannonEntropy(msa, **kwargs)
            matrix = applyMutinfoNorm(mutinfo, entropy, norm=which)
            suffix = '_norm_' + which
            tuffix = ' MI - Normalization: ' + which
        else:
            LOGGER.info('Applying {0} correction.'.format(repr(which)))
            matrix = applyMutinfoCorr(mutinfo, which)
            suffix = '_corr_' + which
            tuffix = ' MI - Correction: ' + which

        writeArray(prefix + suffix + '.txt',
                   matrix, format=kwargs.get('numformat', '%12g'))

        if heatmap:
            writeHeatmap(prefix + suffix + '.hm', matrix,
                         title = msa.getTitle() + tuffix, **hmargs)

        if kwargs.get('figcoevol'):
            try:
                import matplotlib.pyplot as plt
            except ImportError:
                LOGGER.warn('Matplotlib could not be imported, '
                            'figures are not saved.')
            else:
                cmin = kwargs.get('cmin', matrix.min())
                cmax = kwargs.get('cmax', matrix.max())
                prody.SETTINGS['auto_show'] = False
                width = kwargs.get('figwidth', 8)
                height = kwargs.get('figheight', 6)
                xlabel = kwargs.get('xlabel')
                title = kwargs.get('title')
                figure = plt.figure(figsize=(width, height))
                show = showMutinfoMatrix(matrix, msa=msa, clim=(cmin, cmax),
                                         xlabel=xlabel, title=title)

                format = kwargs.get('figformat', 'pdf')
                figure.savefig(prefix + suffix + '.' + format, format=format,
                            dpi=kwargs.get('figdpi', 300))


APP.setFunction(evol_coevol)

