"""Refine MSA application."""

from ..apptools import DevelApp

__all__ = ['evol_filter']

APP = DevelApp('filter', 'filter an MSA using sequence labels')

APP.setExample(
"""Filter sequences in an MSA based on label data.

Following example will filter human sequences:

  $ evol filter piwi_seed.slx HUMAN -e""", [])


APP.addArgument('msa',
    help='MSA filename to be filtered')

APP.addArgument('word',
    help='word to be compared to sequence label', nargs='+')

APP.addGroup('filter', 'filtering method (required)', True, True)
APP.addArgument('-s', '--startswith',
    dest='startswith',
    help='sequence label starts with given words',
    action='store_true',
    group='filter')

APP.addArgument('-e', '--endswith',
    dest='endswith',
    help='sequence label ends with given words',
    action='store_true',
    group='filter')

APP.addArgument('-c', '--contains',
    dest='contains',
    help='sequence label contains with given words',
    action='store_true',
    group='filter')

APP.addGroup('filter2', 'filter option')
APP.addArgument('-F', '--full-label',
    dest='filter_full',
    help='compare full label with word(s)',
    action='store_true',
    group='filter2')


APP.addGroup('output', 'output options')
APP.addArgument('-o', '--outname',
    dest='outname',
    help='output filename, default is msa filename with _refined suffix',
    type=str,
    metavar='STR',
    group='output')

APP.addArgument('-f', '--format',
    dest='format',
    type=str,
    metavar='STR',
    help='output MSA file format, default is same as input',
    group='output')

APP.addArgument('-z', '--compressed',
    dest='compressed',
    action='store_true',
    help='gzip refined MSA output',
    group='output')


def evol_filter(msa, *word, **kwargs):

    import prody
    from prody import MSAFile, writeMSA, LOGGER
    from os.path import splitext

    outname = kwargs.get('outname')
    if outname is None:
        outname, ext = splitext(msa)
        if ext.lower() == '.gz':
            outname, _ = splitext(msa)
        outname += '_filtered' + ext

    single = len(word) == 1
    if single:
        word = word[0]

    if kwargs.get('startswith', False):
        if single:
            filter = lambda label, seq, word=word: label.startswith(word)

    elif kwargs.get('endswith', False):
        if single:
            filter = lambda label, seq, word=word: label.endswith(word)

    elif kwargs.get('contains', False):
        if single:
            filter = lambda label, seq, word=word: word in label

    elif kwargs.get('equals', False):
        if single:
            filter = lambda label, seq, word=word: word == label
        else:
            filter = lambda label, seq, word=set(word): label in word
    else:
        raise TypeError('one of startswith, endswith, contains, or equals '
                        'must be specified')

    msa = MSAFile(msa, filter=filter,
                  filter_full=kwargs.get('filter_full', False))


    LOGGER.info('Filtered MSA is written in file: ' +
                writeMSA(outname, msa, **kwargs))


APP.setFunction(evol_filter)
