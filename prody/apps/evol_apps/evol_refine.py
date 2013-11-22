"""Refine MSA application."""

from ..apptools import DevelApp

__all__ = ['evol_refine']

APP = DevelApp('refine', 'refine an MSA by removing gapped rows/colums')

APP.setExample(
"""Refines MSA by removing gapped columns (residue positions)
and rows (sequences).

In the following, columns that are gaps in sequence with label
GTHB2_ONCKE will be removed from MSA:

  $ evol refine piwi.slx -l GTHB2_ONCKE""", [])


APP.addArgument('msa',
    help='MSA filename to be refined')

APP.addGroup('refine', 'refinement options')
APP.addArgument('-l', '--label',
    dest='label',
    help='sequence label, UniProt ID code, or PDB and chain identifier',
    type=str,
    metavar='STR',
    group='refine')

APP.addArgument('-s', '--seqid',
    dest='seqid',
    help='identity threshold for selecting unique sequences',
    default=None,
    type=float,
    metavar='FLOAT',
    group='refine')

APP.addArgument('-c', '--colocc',
    dest='colocc',
    help='column (residue position) occupancy',
    default=None,
    type=float,
    metavar='FLOAT',
    group='refine')

APP.addArgument('-r', '--rowocc',
    dest='rowocc',
    help='row (sequence) occupancy',
    default=None,
    type=float,
    metavar='FLOAT',
    group='refine')

APP.addArgument('-k', '--keep',
    dest='pdbres',
    help='keep columns corresponding to residues not resolved in PDB '
         'structure, applies label argument is a PDB identifier',
    default=False,
    action='store_true',
    group='refine')

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

def evol_refine(msa, **kwargs):

    import prody
    from prody import parseMSA, refineMSA, writeMSA, LOGGER
    from os.path import splitext

    outname = kwargs.get('outname')
    if outname is None:
        outname, ext = splitext(msa)
        if ext.lower() == '.gz':
            outname, _ = splitext(msa)
        outname += '_refined' + ext

    writeMSA(outname, refineMSA(parseMSA(msa), **kwargs), **kwargs)
    LOGGER.info('Refined MSA is written in file: ' + outname)


APP.setFunction(evol_refine)
