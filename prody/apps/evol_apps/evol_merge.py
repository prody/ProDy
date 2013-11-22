"""Merge multiple MSAs based on common labels."""

__author__ = 'Ahmet Bakan, Anindita Dutta'

from ..apptools import DevelApp

__all__ = ['evol_merge']

APP = DevelApp('merge', 'merge multiple MSAs based on common labels')

APP.setExample(
"""Merges multiple MSAs into one MSA based on common labels
appearing across all the input MSAs. The following example shows
how to merge two MSA of a multi-domain protein.

  $ evol search 3KG2A
  $ evol fetch ANF_receptor
  $ evol fetch Lig_chan
  $ evol merge ANF_receptor_full.slx Lig_chan_full.slx
""", [])


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

