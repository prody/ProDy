"""Pfam MSA download application."""

__author__ = 'Ahmet Bakan, Anindita Dutta'

from ..apptools import DevelApp

__all__ = ['evol_fetch']

APP = DevelApp('fetch', help='fetch MSA files from Pfam')


APP.setExample(
"""Fetch MSA files from Pfam database for given Pfam ID code or
accession.

Fetch PFAM ID Cys_knot:

  $ evol fetch Cys_knot

Fetch PFAM accession with specific parameters:

  $ evol fetch PF00007 --compressed --format fasta --outname mymsa""", [])


APP.addArgument('acc',
    help='Pfam accession or ID',
    type=str)

APP.addGroup('download', 'download options')
APP.addArgument('-a', '--alignment',
    dest='alignment',
    type=str,
    default='full',
    metavar='STR',
    help='alignment type',
    choices='full seed ncbi metagenomics'.split(),
    group='download')
APP.addArgument('-f', '--format',
    dest='format',
    type=str,
    default='selex',
    metavar='STR',
    help='Pfam supported MSA format',
    choices='selex fasta stockholm'.split(),
    group='download')
APP.addArgument('-o', '--order',
    dest='order',
    type=str,
    default='tree',
    metavar='STR',
    help='ordering of sequences',
    choices='tree alphabetical'.split(),
    group='download')
APP.addArgument('-i', '--inserts',
    dest='inserts',
    type=str,
    default='upper',
    metavar='STR',
    help='letter case for inserts',
    choices='upper lower'.split(),
    group='download')
APP.addArgument('-g', '--gaps',
    dest='gaps',
    type=str,
    default='dashes',
    metavar='STR',
    help='gap character',
    choices='dashes dots mixed'.split(),
    group='download')
APP.addArgument('-t', '--timeout',
    dest='timeout',
    type=int,
    default=60,
    metavar='INT',
    help='timeout for blocking connection attempts',
    group='download')

APP.addGroup('output', 'output options')
APP.addArgument('-d', '--outdir',
    dest='folder',
    type=str,
    default='.',
    metavar='PATH',
    help='output directory',
    group='output')
APP.addArgument('-p', '--outname',
    dest='outname',
    type=str,
    default=None,
    metavar='STR',
    help='output filename, default is accession and alignment type',
    group='output')
APP.addArgument('-z', '--compressed',
    dest='compressed',
    action='store_true',
    help='gzip downloaded MSA file',
    group='output')

def evol_fetch(acc, **kwargs):

    import prody
    prody.fetchPfamMSA(acc, **kwargs)


APP.setFunction(evol_fetch)
