# -*- coding: utf-8 -*-
"""This module defines a routine for contact identification."""

from ..apptools import *

__all__ = ['prody_contacts']

def prody_contacts(**kwargs):
    """Identify contacts of a target structure with one or more ligands.
    Contacting atoms (or extended subset of atoms, such as residues) are
    outputted in PDB file format.

    :arg target: target PDB identifier or filename

    :arg ligand: ligand PDB identifier(s) or filename(s)

    :arg select: atom selection string for target structure

    :arg radius: contact radius (Å), default is ``4.0``

    :arg extend: output same ``'residue'``, ``'chain'``, or ``'segment'`` along
        with contacting atoms

    :arg prefix: prefix for output file, default is *target* filename

    :arg suffix: output filename suffix, default is *ligand* filename"""

    import prody
    LOGGER = prody.LOGGER

    target = prody.parsePDB(kwargs['target'])
    title = kwargs.get('prefix') or target.getTitle()
    selstr = kwargs.get('select')
    if selstr:
        target = target.select(selstr)
    contacts = prody.Contacts(target)
    suffix = kwargs.get('suffix', '_contacts')
    extend = kwargs.get('extend')
    radius = float(kwargs.get('radius', 4.0))
    ligands = kwargs.get('ligand')
    if len(ligands) > 1:
        outfn = lambda fn: title + suffix + '_' + fn + '.pdb'
    else:
        outfn = lambda fn: title + suffix + '.pdb'
    for pdb in ligands:
        ligand = prody.parsePDB(pdb)
        sel = contacts(radius, ligand)
        if sel:
            LOGGER.info('{0} atoms from {1} contact {2}.'
                        .format(len(sel), pdb, str(target)))
            if extend:
                sel = target.select('same ' + extend + ' as sel', sel=sel)
                LOGGER.info('Selection is extended to {0} atoms of the same '
                            '{1}(s).'.format(len(sel), extend))
            pdbfn = outfn(ligand.getTitle())
            LOGGER.info('Writing contacts into ' + pdbfn)
            prody.writePDB(pdbfn, sel)


def addCommand(commands):

    subparser = commands.add_parser('contacts',
        help='identify contacts between a target and ligand(s)')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """Identify contacts of a target structure with one or more ligands.

Fetch PDB structure 1zz2, save PDB files for individual ligands, and identify \
contacting residues of the target protein:

    $ prody select -o B11 "resname B11" 1zz2
    $ prody select -o BOG "resname BOG" 1zz2
    $ prody contacts -r 4.0 -t residue -s protein 1zz2 B11.pdb BOG.pdb
    """,
    test_examples=[(0,1,2)]
    )

    subparser.add_argument('-s', '--select', dest='select', type=str,
        metavar='SELSTR', help='selection string for target')

    subparser.add_argument('-r', '--radius', dest='radius', type=float,
        default=4.0, metavar='FLOAT',
        help='contact radius (default: %(default)s)')

    subparser.add_argument('-t', '--extend', dest='extend', type=str,
        metavar='STR', choices=set(['residue', 'chain', 'segment']),
        help=('output same residue, chain, or segment as contacting atoms'))

    subparser.add_argument('-p', '--prefix', dest='prefix', type=str,
        metavar='STR',
        help=('output filename prefix (default: target filename)'))

    subparser.add_argument('-x', '--suffix', dest='suffix', type=str,
        default='_contacts', metavar='STR',
        help=('output filename suffix (default: %(default)s)'))


    subparser.add_argument('target',
        help='target PDB identifier or filename')

    subparser.add_argument('ligand', nargs='+',
        help='ligand PDB identifier(s) or filename(s)')


    subparser.set_defaults(func=lambda opt: prody_contacts(**opt.__dict__))
    subparser.set_defaults(subparser=subparser)
