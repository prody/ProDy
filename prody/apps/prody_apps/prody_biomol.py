"""Generate biomolecule structure using the transformation from the header
section of the PDB file."""

from ..apptools import *

__all__ = ['prody_biomol']

def prody_biomol(pdbname,**kwargs):
    """Generate biomolecule coordinates.

    :arg pdb:  PDB identifier or filename

    :arg prefix: prefix for output files, default is :file:`_biomol`

    :arg biomol: index of the biomolecule, by default all are generated"""

    import prody
    LOGGER = prody.LOGGER
    prefix, biomol = kwargs.get('prefix',None), kwargs.get('biomol')
    pdb, header = prody.parsePDB(pdbname, header=True)
    if not prefix:
        prefix = pdb.getTitle()

    biomols = prody.buildBiomolecules(header, pdb, biomol=biomol)
    if not isinstance(biomols, list):
        biomols = [biomols]

    for i, biomol in enumerate(biomols):
        if isinstance(biomol, prody.Atomic):
            outfn = '{0}_biomol_{1}.pdb'.format(prefix, i+1)
            LOGGER.info('Writing {0}'.format(outfn))
            prody.writePDB(outfn, biomol)
        elif isinstance(biomol, tuple):
            for j, part in enumerate(biomol):
                outfn = ('{0}_biomol_{1}_part_{2}.pdb'
                         .format(prefix, i+1, j+1))
                LOGGER.info('Writing {0}'.format(outfn))
                prody.writePDB(outfn, part)

def addCommand(commands):

    subparser = commands.add_parser('biomol',
        help='build biomolecules')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Generate biomolecule coordinates:

  $ prody biomol 2bfu""",
    test_examples=[0]
    )

    subparser.add_argument('-p', '--prefix', dest='prefix', type=str,
        default=None, metavar='STR',
        help=('prefix for output files (default: pdb_biomol_)'))
    subparser.add_argument('-b', '--biomol', dest='biomol', type=int,
        default=None, metavar='INT',
        help='index of the biomolecule, by default all are generated')

    subparser.add_argument('pdb', help='PDB identifier or filename')

    subparser.set_defaults(func=lambda ns: prody_biomol(ns.pdb, **ns.__dict__))
    subparser.set_defaults(subparser=subparser)

