"""Calculate conservation in an MSA using Shannon entropy."""

from ..apptools import DevelApp

__all__ = ['evol_conserv']

APP = DevelApp('conserv',
               help='analyze conservation using Shannon entropy')

APP.setExample(
"""This application calculates conservation using Shannon entropy for a \
refined multiple sequence alignment.  Following example will save entropy \
data and plot using default options:

    $ evol conserv piwi_refined.slx -S""", [])


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

APP.addGroup('output', 'output options')
APP.addArgument('-p', '--prefix',
    dest='prefix',
    help='output filename prefix, default is '
         'msa filename with _conserv suffix',
    type=str,
    metavar='STR',
    group='output')

APP.addArgument('-f', '--number-format',
    dest='numformat', type=str, default='%12g',
    metavar='STR', help='number output format', group='output')

APP.addFigure('-S', '--save-plot',
    dest='figent',
    action='store_true',
    help='save conservation plot')


def evol_conserv(msa, **kwargs):

    import prody
    from prody import parseMSA, calcShannonEntropy, showShannonEntropy
    from prody import writeArray
    from os.path import splitext

    prefix = kwargs.get('prefix')
    if prefix is None:
        prefix, _ = splitext(msa)
        if _.lower() == '.gz':
            prefix, _ = splitext(prefix)
        prefix += '_conserv'
    msa = parseMSA(msa)
    entropy = calcShannonEntropy(msa, **kwargs)

    writeArray(prefix + '.txt',
               entropy, format=kwargs.get('numformat', '%12g'))

    if kwargs.get('figent'):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warn('Matplotlib could not be imported, '
                        'figures are not saved.')
        else:
            prody.SETTINGS['auto_show'] = False
            width = kwargs.get('figwidth', 8)
            height = kwargs.get('figheight', 6)
            figargs = kwargs.get('figargs', ())
            figure = plt.figure(figsize=(width, height))
            show = showShannonEntropy(entropy, msa=msa, *figargs)
            format = kwargs.get('figformat', 'pdf')
            figure.savefig(prefix + '.' + format, format=format,
                        dpi=kwargs.get('figdpi', 300))


APP.setFunction(evol_conserv)
