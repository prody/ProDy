# -*- coding: utf-8 -*-
"""Perform RTB calculations and output the results in plain text, NMD, and
graphical formats."""

from ..apptools import *
from .nmaoptions import *
from . import nmaoptions

__all__ = ['prody_rtb']

BLOCKS_FROM_RES = 1
BLOCKS_FROM_SECSTR = 2

DEFAULTS = {}
HELPTEXT = {}
for key, txt, val in [
    ('model', 'index of model that will be used in the calculations', 1),
    ('altloc', 'alternative location identifiers for residues used in the calculations', "A"),
    ('cutoff', 'cutoff distance (A)', '15.'),
    ('gamma', 'spring constant', '1.'),
    ('sparse', 'use sparse matrices', False),
    ('kdtree', 'use kdtree for Hessian', False),
    ('zeros', 'calculate zero modes', False),
    ('turbo', 'use memory-intensive turbo option for modes', False),

    ('outbeta', 'write beta-factors calculated from RTB modes', False),
    ('hessian', 'write Hessian matrix', False),
    ('kirchhoff', 'write Kirchhoff matrix', False),
    ('figcmap', 'save contact map (Kirchhoff matrix) figure', False),
    ('figbeta', 'save beta-factors figure', False),
    ('figmode', 'save mode shape figures for specified modes, '
                'e.g. "1-3 5" for modes 1, 2, 3 and 5', ''),
    ('blockInputType', 'type of input for blocks (1 or 2), '
                       'where 1 is number of residues and 2 is secstr elements', 2),
    ('res_per_block', 'Number of residues per block', 10),
    ('shortest_block', 'Number of residues in shortest block', 4),
    ('longest_block', 'Number of residues in longest block (-1 means whole protein)', 20),
    ('min_dist_cutoff', 'Minimum distance cutoff for including in same block', 20.0)]:

    DEFAULTS[key] = val
    HELPTEXT[key] = txt

DEFAULTS.update(nmaoptions.DEFAULTS)
HELPTEXT.update(nmaoptions.HELPTEXT)

DEFAULTS['prefix'] = '_rtb'


def prody_rtb(pdb, **kwargs):
    """Perform RTB calculations for *pdb*.

    """

    for key in DEFAULTS:
        if key not in kwargs:
            kwargs[key] = DEFAULTS[key]

    from os.path import isdir, join, exists
    outdir = kwargs.get('outdir')
    if not isdir(outdir):
        raise IOError('{0} is not a valid path'.format(repr(outdir)))

    import numpy as np
    import prody
    LOGGER = prody.LOGGER

    selstr = kwargs.get('select')
    prefix = kwargs.get('prefix')
    sparse = kwargs.get('sparse')
    kdtree = kwargs.get('kdtree')
    nmodes = kwargs.get('nmodes')
    model = kwargs.get('model')
    altloc = kwargs.get('altloc')
    zeros = kwargs.get('zeros')
    turbo = kwargs.get('turbo')
    membrane = kwargs.get('membrane')
    blockInputType = kwargs.get('blockInputType')
    res_per_block = kwargs.get('res_per_block')
    shortest_block = kwargs.get('shortest_block')
    longest_block = kwargs.get('longest_block')
    min_dist_cutoff = kwargs.get('min_dist_cutoff')

    if membrane and not exists(pdb):
        pdb = prody.fetchPDBfromOPM(pdb)

    pdb = prody.parsePDB(pdb, model=model, altloc=altloc, secondary=True)
    if prefix == '_rtb':
        prefix = pdb.getTitle() + '_rtb'

    select = pdb.select(selstr)
    if select is None:
        LOGGER.warn('Selection {0} did not match any atoms.'
                    .format(repr(selstr)))
        return
    LOGGER.info('{0} atoms will be used for RTB calculations.'
                .format(len(select)))

    rtb = prody.RTB(pdb.getTitle())
    try:
        gamma = float(kwargs.get('gamma'))
        LOGGER.info("Using gamma {0}".format(gamma))
    except ValueError:
        try:
            gamma = eval('prody.' + kwargs.get('gamma'))
            gamma = gamma(select)
            LOGGER.info("Using gamma {0}".format(gamma))
        except NameError:
            raise NameError("Please provide gamma as a float or ProDy Gamma class")
        except TypeError:
            raise TypeError("Please provide gamma as a float or ProDy Gamma class")
        
    try:
        cutoff = float(kwargs.get('cutoff'))
        LOGGER.info("Using cutoff {0}".format(cutoff))
    except ValueError:
        try:
            import math
            cutoff = eval(kwargs.get('cutoff'))
            LOGGER.info("Using cutoff {0}".format(cutoff))
        except NameError:
            raise NameError("Please provide cutoff as a float or equation using math")
        except TypeError:
            raise TypeError("Please provide cutoff as a float or equation using math")

    if blockInputType == BLOCKS_FROM_RES:
        LOGGER.info('Assigning blocks using number of residues...')
        blocks, amap = prody.assignBlocks(select, res_per_block=res_per_block,
                                          shortest_block=shortest_block,
                                          longest_block=longest_block,
                                          min_dist_cutoff=min_dist_cutoff)
    else:
        try:
            LOGGER.info('Trying assigning blocks using secondary structures...')
            blocks, amap = prody.assignBlocks(select, secstr=True,
                                              shortest_block=shortest_block,
                                              longest_block=longest_block,
                                              min_dist_cutoff=min_dist_cutoff)
        except OSError:
            LOGGER.warn('No secondary structures available.')
            LOGGER.info('Assigning blocks using number of residues...')
            blocks, amap = prody.assignBlocks(select, res_per_block=res_per_block,
                                              shortest_block=shortest_block,
                                              longest_block=longest_block,
                                              min_dist_cutoff=min_dist_cutoff)
    LOGGER.info('Assigning blocks done.')

    nproc = kwargs.get('nproc')

    LOGGER.info('Building Hessian matrix (regular then block Hessian).')
    rtb.buildHessian(amap, blocks, cutoff, gamma, sparse=sparse, kdtree=kdtree)
    LOGGER.info('Calculating modes.')
    rtb.calcModes(nmodes, zeros=zeros, turbo=turbo, nproc=nproc)

    LOGGER.info('Writing numerical output.')

    if kwargs.get('outnpz'):
        prody.saveModel(rtb, join(outdir, prefix), 
                        matrices=kwargs.get('npzmatrices'))

    if kwargs.get('outscipion'):
        prody.writeScipionModes(outdir, rtb)

    prody.writeNMD(join(outdir, prefix + '.nmd'), rtb, amap)

    extend = kwargs.get('extend')
    if extend:
        if extend == 'all':
            extended = prody.extendModel(rtb, amap, pdb)
        else:
            extended = prody.extendModel(rtb, amap, amap | pdb.bb)
        prody.writeNMD(join(outdir, prefix + '_extended_' +
                       extend + '.nmd'), *extended)

    outall = kwargs.get('outall')
    delim = kwargs.get('numdelim')
    ext = kwargs.get('numext')
    format_ = kwargs.get('numformat')


    if outall or kwargs.get('outeig'):
        prody.writeArray(join(outdir, prefix + '_evectors'+ext),
                         rtb.getArray(), delimiter=delim, format=format_)
        prody.writeArray(join(outdir, prefix + '_evalues'+ext),
                         rtb.getEigvals(), delimiter=delim, format=format_)

    if outall or kwargs.get('outbeta'):
        from prody.utilities import openFile
        fout = openFile(prefix + '_beta'+ext, 'w', folder=outdir)
        fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4s} {0[3]:5s} {0[4]:5s}\n'
                       .format(['C', 'RES', '####', 'Exp.', 'The.']))
        for data in zip(amap.getChids(), amap.getResnames(),
                        amap.getResnums(), amap.getBetas(),
                        prody.calcTempFactors(rtb, amap)):
            fout.write('{0[0]:1s} {0[1]:4s} {0[2]:4d} {0[3]:5.2f} {0[4]:5.2f}\n'
                       .format(data))
        fout.close()

    if outall or kwargs.get('outcov'):
        prody.writeArray(join(outdir, prefix + '_covariance' + ext),
                         rtb.getCovariance(), delimiter=delim, format=format_)

    if outall or kwargs.get('outcc') or kwargs.get('outhm'):
        cc = prody.calcCrossCorr(rtb)
        if outall or kwargs.get('outcc'):
            prody.writeArray(join(outdir, prefix +
                             '_cross-correlations' + ext),
                             cc, delimiter=delim,  format=format_)
        if outall or kwargs.get('outhm'):
            prody.writeHeatmap(join(outdir, prefix + '_cross-correlations.hm'),
                               cc, resnum=amap.getResnums(),
                               xlabel='Residue', ylabel='Residue',
                               title=rtb.getTitle() + ' cross-correlations')

    if outall or kwargs.get('hessian'):
        prody.writeArray(join(outdir, prefix + '_hessian'+ext),
                         rtb.getHessian(), delimiter=delim, format=format_)

    if outall or kwargs.get('kirchhoff'):
        prody.writeArray(join(outdir, prefix + '_kirchhoff'+ext),
                         rtb.getKirchhoff(), delimiter=delim, format=format_)

    if outall or kwargs.get('outsf'):
        prody.writeArray(join(outdir, prefix + '_sqflucts'+ext),
                         prody.calcSqFlucts(rtb), delimiter=delim,
                         format=format_)

    figall = kwargs.get('figall')
    cc = kwargs.get('figcc')
    sf = kwargs.get('figsf')
    bf = kwargs.get('figbeta')
    cm = kwargs.get('figcmap')


    if figall or cc or sf or bf or cm:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            LOGGER.warning('Matplotlib could not be imported. '
                           'Figures are not saved.')
        else:
            prody.SETTINGS['auto_show'] = False
            LOGGER.info('Saving graphical output.')
            format_ = kwargs.get('figformat')
            width = kwargs.get('figwidth')
            height = kwargs.get('figheight')
            dpi = kwargs.get('figdpi')
            format_ = format_.lower()

            if figall or cc:
                plt.figure(figsize=(width, height))
                prody.showCrossCorr(rtb, interactive=False)
                plt.savefig(join(outdir, prefix + '_cc.'+format_),
                    dpi=dpi, format=format_)
                plt.close('all')

            if figall or cm:
                plt.figure(figsize=(width, height))
                prody.showContactMap(rtb, interactive=False)
                plt.savefig(join(outdir, prefix + '_cm.'+format_),
                    dpi=dpi, format=format_)
                plt.close('all')

            if figall or sf:
                plt.figure(figsize=(width, height))
                prody.showSqFlucts(rtb)
                plt.savefig(join(outdir, prefix + '_sf.'+format_),
                    dpi=dpi, format=format_)
                plt.close('all')

            if figall or bf:
                plt.figure(figsize=(width, height))
                bexp = amap.getBetas()
                bcal = prody.calcTempFactors(rtb, amap)
                plt.plot(bexp, label='Experimental')
                plt.plot(bcal, label=('Theoretical (R={0:.2f})'
                                        .format(np.corrcoef(bcal, bexp)[0,1])))
                plt.legend(prop={'size': 10})
                plt.xlabel('Node index')
                plt.ylabel('Experimental B-factors')
                plt.title(pdb.getTitle() + ' B-factors')
                plt.savefig(join(outdir, prefix + '_bf.'+format_),
                    dpi=dpi, format=format_)
                plt.close('all')


_ = list(HELPTEXT)
_.sort()
for key in _:

    prody_rtb.__doc__ += """
    :arg {0}: {1}, default is ``{2!r}``""".format(key, HELPTEXT[key],
                                                  DEFAULTS[key])

def addCommand(commands):

    subparser = commands.add_parser('rtb',
        help='perform rotating and translating blocks normal mode analysis calculations')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')
    subparser.set_defaults(usage_example=
"""Perform RTB calculations for given PDB structure and output results in NMD
format.  If an identifier is passed, structure file will be downloaded from
the PDB FTP server.

Fetch PDB 1p38, run RTB calculations using default parameters, and write
NMD file:

  $ prody rtb 1p38

Fetch PDB 1aar, run RTB calculations using default parameters for chain A
carbon alpha atoms with residue numbers less than 70, and save all of the
graphical output files:

  $ prody rtb 1aar -s "calpha and chain A and resnum < 70" -A
""",
  test_examples=[0, 1])

    group = addNMAParameters(subparser)

    group.add_argument('-c', '--cutoff', dest='cutoff', type=str,
        default=DEFAULTS['cutoff'], metavar='FLOAT',
        help=HELPTEXT['cutoff'] + ' (default: %(default)s)')

    group.add_argument('-g', '--gamma', dest='gamma', type=str,
        default=DEFAULTS['gamma'], metavar='STR',
        help=HELPTEXT['gamma'] + ' (default: %(default)s)')

    group.add_argument('-C', '--sparse-hessian', dest='sparse', action='store_true',
        default=DEFAULTS['sparse'],
        help=HELPTEXT['sparse'] + ' (default: %(default)s)')

    group.add_argument('-G', '--use-kdtree', dest='kdtree', action='store_true',
        default=DEFAULTS['kdtree'],
        help=HELPTEXT['kdtree'] + ' (default: %(default)s)')

    group.add_argument('-y', '--turbo', dest='turbo', action='store_true',
        default=DEFAULTS['turbo'],
        help=HELPTEXT['turbo'] + ' (default: %(default)s)')

    group.add_argument('-m', '--model', dest='model', type=int,
        metavar='INT', default=DEFAULTS['model'],
        help=HELPTEXT['model'] + ' (default: %(default)s)')

    group.add_argument('-L', '--altloc', dest='altloc', type=str,
        metavar='STR', default=DEFAULTS['altloc'],
        help=HELPTEXT['altloc'] + ' (default: %(default)s)')

    group.add_argument('-w', '--zero-modes', dest='zeros', action='store_true',
        default=DEFAULTS['zeros'], help=HELPTEXT['zeros'] + ' (default: %(default)s)')

    group.add_argument('-i', '--block-input-type', dest='blockInputType', 
        type=int, metavar='INT', default=DEFAULTS['blockInputType'], 
        help=HELPTEXT['blockInputType'] + ' (default: %(default)s)')

    group.add_argument('-M', '--res-per-block', dest='res_per_block', 
        type=int, metavar='INT', default=DEFAULTS['res_per_block'], 
        help=HELPTEXT['res_per_block'] + ' (default: %(default)s)')

    group.add_argument('-j', '--shortest-block', dest='shortest_block', 
        type=int, metavar='INT', default=DEFAULTS['shortest_block'], 
        help=HELPTEXT['shortest_block'] + ' (default: %(default)s)')

    group.add_argument('-U', '--longest-block', dest='longest_block', 
        type=int, metavar='INT', default=DEFAULTS['longest_block'], 
        help=HELPTEXT['longest_block'] + ' (default: %(default)s)')

    group.add_argument('-V', '--min-block-dist-cutoff', dest='min_dist_cutoff',
        type=float, metavar='FLOAT', default=DEFAULTS['min_dist_cutoff'],
        help=HELPTEXT['min_dist_cutoff'] + ' (default: %(default)s)')

    group = addNMAOutput(subparser)

    group.add_argument('-b', '--beta-factors', dest='outbeta',
        action='store_true', default=DEFAULTS['outbeta'],
        help=HELPTEXT['outbeta'])

    group.add_argument('-l', '--hessian', dest='hessian',
        action='store_true',
        default=DEFAULTS['hessian'], help=HELPTEXT['hessian'])

    group.add_argument('-k', '--kirchhoff', dest='kirchhoff',
        action='store_true',
        default=DEFAULTS['kirchhoff'], help=HELPTEXT['kirchhoff'])


    group = addNMAOutputOptions(subparser, '_rtb')

    group = addNMAFigures(subparser)

    group.add_argument('-B', '--beta-factors-figure', dest='figbeta',
        action='store_true',
        default=DEFAULTS['figbeta'], help=HELPTEXT['figbeta'])

    group.add_argument('-K', '--contact-map', dest='figcmap',
        action='store_true',
        default=DEFAULTS['figcmap'], help=HELPTEXT['figcmap'])

    group = addNMAFigureOptions(subparser)

    subparser.add_argument('pdb', help='PDB identifier or filename')

    subparser.set_defaults(func=lambda ns: prody_rtb(ns.__dict__.pop('pdb'),
                                                     **ns.__dict__))

    subparser.set_defaults(subparser=subparser)
