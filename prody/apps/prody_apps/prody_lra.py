"""Perform logistic regression analysis (LRA) calculations and output the results in plain text, NMD, and
graphical formats."""

from ..apptools import *
from .nmaoptions import *
from . import nmaoptions
from numbers import Integral

DEFAULTS = {}
HELPTEXT = {}
for key, txt, val in [
    ('aligned', 'trajectory is already aligned', False),
    ('altloc', 'alternative location identifiers for residues used in the calculations', "A"),
    ('outproj', 'write projections onto LR component', False),
    ('figproj', 'save projections onto specified subspaces, e.g. '
                '"1,2" for projections onto LRCs 1 and 2; '
                '"1,2 1,3" for projections onto LRCs 1,2 and 1,3; '
                '"1 1,2,3" for projections onto LRCs 1 and 1, 2, 3', ''),
    ('nShuffles', 'number of random shuffles to make', 10),
    ('labels', 'labels for classes to separate. If None, get them from k-medoids clustering using nmodes+1', 'None'),
    ('nmodes', 'number of modes of separation to calculate if labels is None', 1)]:

    DEFAULTS[key] = val
    HELPTEXT[key] = txt

DEFAULTS.update(nmaoptions.DEFAULTS)
HELPTEXT.update(nmaoptions.HELPTEXT)

# override nmodes to 1
key, txt, val = ('nmodes', 'number of modes of separation to calculate if labels is None', 1)
DEFAULTS[key] = val
HELPTEXT[key] = txt

DEFAULTS['prefix'] = '_lra'

__all__ = ['prody_lra']

def prody_lra(coords, **kwargs):
    """Perform LRA calculations for PDB or DCD format *coords* file.
    """

    for key in DEFAULTS:
        if not key in kwargs:
            kwargs[key] = DEFAULTS[key]

    from os.path import isdir, splitext, join
    outdir = kwargs.get('outdir')
    if not isdir(outdir):
        raise IOError('{0} is not a valid path'.format(repr(outdir)))

    import prody
    LOGGER = prody.LOGGER

    prefix = kwargs.get('prefix')
    nmodes = kwargs.get('nmodes')
    selstr = kwargs.get('select')
    altloc = kwargs.get('altloc')

    nproc = kwargs.get('nproc')
    nShuffles = kwargs.get('nShuffles')
    labels = eval(kwargs.get('labels'))

    ext = splitext(coords)[1].lower()
    if ext == '.gz':
        ext = splitext(coords[:-3])[1].lower()

    if ext == '.dcd':
        pdb = kwargs.get('psf') or kwargs.get('pdb')
        if pdb:
            if splitext(pdb)[1].lower() == '.psf':
                pdb = prody.parsePSF(pdb)
            else:
                pdb = prody.parsePDB(pdb, altloc=altloc)
        dcd = prody.DCDFile(coords)
        if prefix == '_lra':
            prefix = dcd.getTitle() + prefix

        if len(dcd) < 2:
            raise ValueError('DCD file must have multiple frames')
        if pdb:
            if pdb.numAtoms() == dcd.numAtoms():
                select = pdb.select(selstr)
                dcd.setAtoms(select)
                LOGGER.info('{0} atoms are selected for calculations.'
                            .format(len(select)))
            else:
                select = pdb.select(selstr)
                if select.numAtoms() != dcd.numAtoms():
                    raise ValueError('number of selected atoms ({0}) does '
                                     'not match number of atoms in the DCD '
                                     'file ({1})'.format(select.numAtoms(),
                                                           dcd.numAtoms()))
                if pdb.numCoordsets():
                    dcd.setCoords(select.getCoords())
        else:
            select = prody.AtomGroup()
            select.setCoords(dcd.getCoords())

        lra = prody.LRA(dcd.getTitle())
        ensemble = prody.Ensemble(dcd)
        ensemble.setCoords(dcd.getCoords())
        ensemble.addCoordset(dcd.getCoordsets())
    else:
        pdb = prody.parsePDB(coords)
        if pdb.numCoordsets() < 2:
            raise ValueError('PDB file must contain multiple models')

        if prefix == '_lra':
            prefix = pdb.getTitle() + prefix

        select = pdb.select(selstr)
        LOGGER.info('{0} atoms are selected for calculations.'
                    .format(len(select)))
        if select is None:
            raise ValueError('selection {0} do not match any atoms'
                                .format(repr(selstr)))
        LOGGER.info('{0} atoms will be used for LRA calculations.'
                    .format(len(select)))
        ensemble = prody.Ensemble(select)
        lra = prody.LRA(pdb.getTitle())
        if not kwargs.get('aligned'):
            ensemble.iterpose()

    if labels is None:
        _, labels, _ = prody.calcKmedoidClusters(ensemble, nmodes+1)

    nproc = kwargs.get('nproc')
    if nproc:
        try:
            from threadpoolctl import threadpool_limits
        except ImportError:
            raise ImportError('Please install threadpoolctl to control threads')

        with threadpool_limits(limits=nproc, user_api="blas"):
            lra.calcModes(ensemble, labels, n_shuffles=nShuffles)
    else:
        lra.calcModes(ensemble, labels, n_shuffles=nShuffles)

    LOGGER.info('Writing numerical output.')
    if kwargs.get('outnpz'):
        prody.saveModel(lra, join(outdir, prefix), 
                        matrices=kwargs.get('npzmatrices'))

    if kwargs.get('outscipion'):
        prody.writeScipionModes(outdir, lra)

    prody.writeNMD(join(outdir, prefix + '.nmd'), lra[:nmodes], select)

    extend = kwargs.get('extend')
    if extend:
        if pdb:
            if extend == 'all':
                extended = prody.extendModel(lra[:nmodes], select, pdb)
            else:
                extended = prody.extendModel(lra[:nmodes], select,
                                             select | pdb.bb)
            prody.writeNMD(join(outdir, prefix + '_extended_' +
                           extend + '.nmd'), *extended)
        else:
            prody.LOGGER.warn('Model could not be extended, provide a PDB or '
                              'PSF file.')
    outall = kwargs.get('outall')
    delim = kwargs.get('numdelim')
    ext = kwargs.get('numext')
    format = kwargs.get('numformat')

    if outall or kwargs.get('outeig'):
        prody.writeArray(join(outdir, prefix + '_evectors'+ext),
                         lra.getArray(), delimiter=delim, format=format)
        prody.writeArray(join(outdir, prefix + '_evalues'+ext),
                         lra.getEigvals(), delimiter=delim, format=format)
    if outall or kwargs.get('outcov'):
        prody.writeArray(join(outdir, prefix + '_covariance'+ext),
                         lra.getCovariance(), delimiter=delim, format=format)
    if outall or kwargs.get('outcc') or kwargs.get('outhm'):
        cc = prody.calcCrossCorr(lra)
        if outall or kwargs.get('outcc'):
            prody.writeArray(join(outdir, prefix + '_cross-correlations' +
                             ext), cc, delimiter=delim, format=format)
        if outall or kwargs.get('outhm'):
            resnums = select.getResnums()
            hmargs = {} if resnums is None else {'resnums': resnums}
            prody.writeHeatmap(join(outdir, prefix + '_cross-correlations.hm'),
                               cc, xlabel='Residue', ylabel='Residue',
                               title=lra.getTitle() + ' cross-correlations',
                               **hmargs)

    if outall or kwargs.get('outsf'):
        prody.writeArray(join(outdir, prefix + '_sqfluct'+ext),
                         prody.calcSqFlucts(lra), delimiter=delim,
                         format=format)
    if outall or kwargs.get('outproj'):
        prody.writeArray(join(outdir, prefix + '_proj'+ext),
                         prody.calcProjection(ensemble, lra), delimiter=delim,
                         format=format)

    figall = kwargs.get('figall')
    cc = kwargs.get('figcc')
    sf = kwargs.get('figsf')
    sp = kwargs.get('figproj')

    if figall or cc or sf or sp:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warning('Matplotlib could not be imported. '
                           'Figures are not saved.')
        else:
            prody.SETTINGS['auto_show'] = False
            LOGGER.info('Saving graphical output.')
            format = kwargs.get('figformat')
            width = kwargs.get('figwidth')
            height = kwargs.get('figheight')
            dpi = kwargs.get('figdpi')

            format = format.lower()
            if figall or cc:
                plt.figure(figsize=(width, height))
                prody.showCrossCorr(lra, interactive=False)
                plt.savefig(join(outdir, prefix + '_cc.'+format),
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or sf:
                plt.figure(figsize=(width, height))
                prody.showSqFlucts(lra)
                plt.savefig(join(outdir, prefix + '_sf.'+format),
                    dpi=dpi, format=format)
                plt.close('all')
            if figall or sp:
                indices = []
                for item in sp.split():
                    try:
                        if '-' in item:
                            item = item.split('-')
                            if len(item) == 2:
                                indices.append(list(range(int(item[0])-1,
                                                          int(item[1]))))
                        elif ',' in item:
                            indices.append([int(i)-1 for i in item.split(',')])
                        else:
                            indices.append(int(item)-1)
                    except:
                        pass
                for index in indices:
                        plt.figure(figsize=(width, height))
                        prody.showProjection(ensemble, lra[index])
                        if isinstance(index, Integral):
                            index = [index]
                        index = [str(i+1) for i in index]
                        plt.savefig(join(outdir, prefix + '_proj_' +
                            '_'.join(index) + '.' + format),
                            dpi=dpi, format=format)
                        plt.close('all')


_ = list(HELPTEXT)
_.sort()
for key in _:

    prody_lra.__doc__ += """
    :arg {0}: {1}, default is ``{2!r}``""".format(key, HELPTEXT[key],
                                                  DEFAULTS[key])


def addCommand(commands):

    subparser = commands.add_parser('lra',
        help='perform principal component analysis calculations')

    subparser.add_argument('--quiet', help='suppress info messages to stderr',
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """This command performs logistic regression analysis (LRA) calculations \
for given multi-model PDB structure or DCD format trajectory file and outputs \
results in NMD format.  If a PDB identifier is given, structure file will be \
downloaded from the PDB FTP server.  DCD files may be accompanied with \
PDB or PSF files to enable atoms selections.

Fetch pdb 2k39, perform LRA calculations, and output NMD file:

  $ prody lra 2k39

Fetch pdb 2k39 and perform calculations for backbone of residues up to 71, \
and save all output and figure files:

  $ prody lra 2k39 --select "backbone and resnum < 71" -a -A

Extra arguments are as follows:

-i or --shuffles sets the number of shuffles for identifying key residues
-l or --labels sets labels for classes to separate and overrides the number of modes

If labels are not set, then they are calculated by k-medoids clustering with k = nmodes + 1.
The default value of nmodes is 1 for separating 2 clusters.
""",
    test_examples=[0, 1])

    group = addNMAParameters(subparser, include_nmodes=False)

    group = addNMAOutput(subparser)

    group.add_argument('-i', '--shuffles', dest='nShuffles',
        type=int,
        default=DEFAULTS['nShuffles'], help=HELPTEXT['nShuffles'])
    
    group.add_argument('-l', '--labels', dest='labels',
        type=str,
        default=DEFAULTS['labels'], help=HELPTEXT['labels'])
    
    group.add_argument('-n', '--number-of-modes', dest='nmodes', type=int,
        default=DEFAULTS['nmodes'], metavar='INT',
        help=HELPTEXT['nmodes'] + ' (default: %(default)s)')

    group.add_argument('-j', '--projection', dest='outproj',
        action='store_true',
        default=DEFAULTS['outproj'], help=HELPTEXT['outproj'])

    group = addNMAOutputOptions(subparser, '_lra')

    group = addNMAFigures(subparser)

    group.add_argument('-J', '--projection-figure', dest='figproj', type=str,
        default=DEFAULTS['figproj'], metavar='STR',
        help=HELPTEXT['figproj'])

    group = addNMAFigureOptions(subparser)

    group = subparser.add_mutually_exclusive_group()
    group.add_argument('--psf', help='PSF filename')
    group.add_argument('--pdb', help='PDB filename')
    group.add_argument('-L', '--altloc', dest='altloc', type=str,
        metavar='INT', default=DEFAULTS['altloc'], help=HELPTEXT['altloc'])

    subparser.add_argument('--aligned', dest='aligned', action='store_true',
        default=DEFAULTS['aligned'], help=HELPTEXT['aligned'])

    subparser.add_argument('dcd', help='file in DCD or PDB format')

    subparser.set_defaults(func=lambda ns: prody_lra(ns.dcd, **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
