"""Download PDB files for given identifiers."""

import os.path

from ..apptools import *

__all__ = ['prody_fetch']

def prody_fetch(*pdb, **kwargs):
    """Fetch PDB files from PDB FTP server.

    :arg pdbs: PDB identifier(s) or filename(s)

    :arg dir: target directory for saving PDB file(s), default is ``'.'``

    :arg gzip: gzip fetched files or not, default is **True**"""

    import prody


    pdblist = pdb
    if len(pdblist) == 1 and os.path.isfile(pdblist[0]):
        from prody.utilities import openFile
        with openFile(pdblist[0]) as inp:
            for item in inp.read().strip().split():
                for pdb in item.split(','):
                    if len(pdb) == 4 and pdb.isalnum():
                        pdblist.append(pdb)

    prody.fetchPDB(*pdblist, folder=kwargs.get('folder', '.'),
                   compressed=kwargs.get('gzip', False),
                   copy=True)

def addCommand(commands):

    subparser = commands.add_parser('fetch',
        help='fetch a PDB file')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Download PDB file(s) by specifying identifiers:

  $ prody fetch 1mkp 1p38""",
    test_examples=[0]
    )

    subparser.add_argument('-d', '--dir', dest='folder', type=str,
        default='.', metavar='PATH',
        help=('target directory for saving PDB file(s)'))

    subparser.add_argument('-z', '--gzip', dest='gzip', action='store_true',
        default=False, help='write compressed PDB file(S)')

    subparser.add_argument('pdb', nargs='+',
        help='PDB identifier(s) or a file that contains them')

    subparser.set_defaults(func=lambda ns: prody_fetch(*ns.pdb, **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
