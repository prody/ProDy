
# -*- coding: utf-8 -*-
"""Run the ClustENM(D) hybrid simulation method, 
combining clustering, ENM NMA and MD."""

from ..apptools import *
from .nmaoptions import *
from . import nmaoptions

__all__ = ['prody_clustenm']

DEFAULTS = {}
HELPTEXT = {}
for key, txt, val in [
    ('outdir', 'output directory', '.'),
    ('model', 'index of model that will be used in the simulations', 1),
    ('altloc', 'alternative location identifiers for residues used in the simulations', "A"),
    ('cutoff', 'cutoff distance (A)', '15.'),
    ('gamma', 'spring constant', '1.'),
    ('sparse', 'use sparse matrices', False),
    ('kdtree', 'use kdtree for Hessian', False),
    ('turbo', 'use memory-intensive turbo option for modes', False),

    ('ngens', 'number of generations', 5),
    ('nconfs', 'number of new conformers from each one from previous generation', 50),
    ('rmsd', 'average RMSD of new conformers from previous ones, can be tuple of floats', '1.'),
    ('maxclust', 'maximum number of clusters for each generation, can be tuple of floats', 'None'),
    ('threshold', 'RMSD threshold to apply when forming clusters, can be tuple of floats', 'None'),
    ('no_sim', 'whether a short MD simulation is not performed after energy minimization (otherwise it is)', False),
    ('parallel', 'whether conformer generation will be parallelized', False),
    ('v1', 'whether to use original sampling method with complete enumeration of ANM modes', False),
    ('no_outlier', 'whether to not exclude outliers in each generation when using implicit solvent (always False for explicit)', False),
    ('mzscore', 'modified z-score threshold to label conformers as outliers', 3.5),
    ('tolerance', 'energy tolerance to which the system should be minimized in kJ/mole', 10.),
    ('maxIterations', 'maximum number of iterations of energy minimization, 0 means until convergence', 0),
    ('temp', 'temperature at which simulations are conducted', 303.15),
    ('ionicStrength', 'total concentration (M) of ions (positive and negative), excluding those to neutralize', 0.),
    ('padding', 'padding distance to use for solvation (nm)', 1.),
    ('solvent', 'solvent model to be used, either imp for implicit or exp for explicit', 'imp'),
    ('forcefield', 'Alternative force field pair tuple for protein then solvent from openmm', 'None'),
    ('t_steps_i', 'number of 2.0 fs MD time steps for initial structure', 1000),
    ('t_steps_g', 'number of 2.0 fs MD time steps in each generation, can be tuple of floats', '7500'),
    ('multiple', 'whether each conformer will be saved as a separate PDB file', False),
    ('write_params', 'whether to write parameters', False)]:

    DEFAULTS[key] = val
    HELPTEXT[key] = txt

DEFAULTS.update(nmaoptions.DEFAULTS)
HELPTEXT.update(nmaoptions.HELPTEXT)

DEFAULTS['prefix'] = '_clustenm'


def prody_clustenm(pdb, **kwargs):
    """Run ClustENM(D) for *pdb*.
    """
    for key in DEFAULTS:
        if key not in kwargs:
            kwargs[key] = DEFAULTS[key]

    from os.path import isdir, join
    outdir = kwargs.pop('outdir')
    if not isdir(outdir):
        raise IOError('{0} is not a valid path'.format(repr(outdir)))

    import prody
    LOGGER = prody.LOGGER

    selstr = kwargs.pop('select')
    prefix = kwargs.pop('prefix')
    sparse = kwargs.pop('sparse')
    kdtree = kwargs.pop('kdtree')
    nmodes = kwargs.pop('nmodes')
    model = kwargs.pop('model')
    altloc = kwargs.pop('altloc')
    turbo = kwargs.pop('turbo')

    ngens = kwargs.pop('ngens')
    nconfs = kwargs.pop('nconfs')
    rmsd = kwargs.pop('rmsd')
    maxclust = kwargs.pop('maxclust')
    threshold = kwargs.pop('threshold')
    sim = not kwargs.pop('no_sim')
    temp = kwargs.pop('temp')
    parallel = kwargs.pop('parallel')
    write_params = kwargs.pop("write_params")

    solvent = kwargs.pop('solvent')
    forcefield = kwargs.pop('forcefield')
    t_steps_i = kwargs.pop('t_steps_i')
    t_steps_g = kwargs.pop('t_steps_g')
    outlier = not kwargs.pop('no_outlier')
    mzscore = kwargs.pop('mzscore')

    pdb = prody.parsePDB(pdb, model=model, altloc=altloc)
    if prefix == '_clustenm':
        prefix = pdb.getTitle() + '_clustenm'

    select = pdb.select(selstr)
    if select is None:
        LOGGER.warn('Selection {0} did not match any atoms.'
                    .format(repr(selstr)))
        return
    LOGGER.info('{0} atoms will be used for ANM calculations.'
                .format(len(select)))
    
    try:
        gamma = float(kwargs.pop('gamma'))
        LOGGER.info("Using gamma {0}".format(gamma))
    except ValueError:
        try:
            gamma = eval('prody.' + kwargs.pop('gamma'))
            gamma = gamma(select)
            LOGGER.info("Using gamma {0}".format(gamma))
        except NameError:
            raise NameError("Please provide gamma as a float or ProDy Gamma class")
        except TypeError:
            raise TypeError("Please provide gamma as a float or ProDy Gamma class")
        
    try:
        cutoff = float(kwargs.pop('cutoff'))
        LOGGER.info("Using cutoff {0}".format(cutoff))
    except ValueError:
        try:
            import math
            cutoff = eval(kwargs.pop('cutoff'))
            LOGGER.info("Using cutoff {0}".format(cutoff))
        except NameError:
            raise NameError("Please provide cutoff as a float or equation using math")
        except TypeError:
            raise TypeError("Please provide cutoff as a float or equation using math")

    ens = prody.ClustENM(pdb.getTitle())
    ens.setAtoms(select)
    ens.run(n_gens=ngens, n_modes=nmodes,
            n_confs=nconfs, rmsd=eval(rmsd),
            cutoff=cutoff, gamma=gamma,
            maxclust=eval(maxclust), threshold=eval(threshold),
            solvent=solvent, force_field=eval(forcefield),
            sim=sim, temp=temp, t_steps_i=t_steps_i,
            t_steps_g=eval(t_steps_g),
            outlier=outlier, mzscore=mzscore,
            sparse=sparse, kdtree=kdtree, turbo=turbo,
            parallel=parallel, **kwargs)

    single = not kwargs.pop('multiple')
    outname = join(outdir, prefix)
    ens.writePDB(outname, single=single)

    prody.saveEnsemble(ens, outname)

    if write_params:
        ens.writeParameters(outname + '.txt')

_ = list(HELPTEXT)
_.sort()
for key in _:

    prody_clustenm.__doc__ += """
    :arg {0}: {1}, default is ``{2!r}``""".format(key, HELPTEXT[key],
                                                  DEFAULTS[key])

def addCommand(commands):

    subparser = commands.add_parser('clustenm',
        help='run clustenm(d) simulations')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')
    subparser.set_defaults(usage_example=
"""run ClustENM(D) simulations for given PDB structure and output results in PDB and 
NPZ format.  If an identifier is passed, structure file will be downloaded from
the PDB FTP server.

Fetch PDB 1p38, run ClustENM(D) simulations using default parameters, and write
PDB and NPZ files:

  $ prody clustenm 1p38

Fetch PDB 1aar, run ClustENM(D) simulations using default parameters for chain A
carbon alpha atoms with residue numbers less than 70, and save all of the
graphical output files:

  $ prody clustenm 1aar -s "calpha and chain A and resnum < 70" -A""",
  test_examples=[0, 1])

    group = addNMAParameters(subparser, include_nproc=False)

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
        metavar='INT', default=DEFAULTS['model'], help=HELPTEXT['model'])

    group.add_argument('-L', '--altloc', dest='altloc', type=str,
        metavar='STR', default=DEFAULTS['altloc'], help=HELPTEXT['altloc'])

    group.add_argument('-j', '--ngens', dest='ngens', type=int,
        metavar='INT', default=DEFAULTS['ngens'], 
        help=HELPTEXT['ngens'] + ' (default: %(default)s)')
    
    group.add_argument('-k', '--nconfs', dest='nconfs', type=int,
        metavar='INT', default=DEFAULTS['nconfs'],
        help=HELPTEXT['nconfs'] + ' (default: %(default)s)')
    
    group.add_argument('-l', '--rmsd', dest='rmsd', type=str,
        default=DEFAULTS['rmsd'], metavar='STR',
        help=HELPTEXT['rmsd'] + ' (default: %(default)s)')
    
    group.add_argument('-J', '--maxclust', dest='maxclust', type=str,
        default=DEFAULTS['maxclust'], metavar='STR',
        help=HELPTEXT['maxclust'] + ' (default: %(default)s)')
    
    group.add_argument('-K', '--threshold', dest='threshold', type=str,
        default=DEFAULTS['threshold'], metavar='STR',
        help=HELPTEXT['threshold'] + ' (default: %(default)s)')
    
    group.add_argument('-b', '--no-sim', dest='no_sim',
        action='store_true',
        default=DEFAULTS['no_sim'], help=HELPTEXT['no_sim'])
    
    group.add_argument('-B', '--parallel', dest='parallel',
        action='store_true',
        default=DEFAULTS['parallel'], help=HELPTEXT['parallel'])
    
    group.add_argument('-E', '--write_params', dest='write_params',
        action='store_true',
        default=DEFAULTS['write_params'], help=HELPTEXT['write_params'])
    
    group.add_argument('-S', '--solvent', dest='solvent', type=str,
        default=DEFAULTS['solvent'], metavar='STR',
        help=HELPTEXT['solvent'] + ' (default: %(default)s)')
    
    group.add_argument('-I', '--ionicStrength', dest='ionicStrength', type=float,
        default=DEFAULTS['ionicStrength'], metavar='FLOAT',
        help=HELPTEXT['ionicStrength'] + ' (default: %(default)s)')
    
    group.add_argument('-P', '--padding', dest='padding', type=float,
        default=DEFAULTS['padding'], metavar='FLOAT',
        help=HELPTEXT['padding'] + ' (default: %(default)s)')

    group.add_argument('-f', '--force_field', dest='forcefield', type=str,
        default=DEFAULTS['forcefield'], metavar='STR',
        help=HELPTEXT['forcefield'] + ' (default: %(default)s)')

    group.add_argument('-D', '--maxIterations', dest='maxIterations', type=int,
        metavar='INT', default=DEFAULTS['maxIterations'], 
        help=HELPTEXT['maxIterations'] + ' (default: %(default)s)')
    
    subparser.add_argument('-o', '--output-dir', dest='outdir', type=str,
        default=DEFAULTS['outdir'], metavar='PATH',
        help=HELPTEXT['outdir'] + ' (default: %(default)s)')
    
    group.add_argument('-t', '--temp', dest='temp', type=float,
        default=DEFAULTS['temp'], metavar='FLOAT',
        help=HELPTEXT['temp'] + ' (default: %(default)s)')
    
    group.add_argument('-i', '--t_steps_i', dest='t_steps_i', type=int,
        default=DEFAULTS['t_steps_i'], metavar='INT',
        help=HELPTEXT['t_steps_i'] + ' (default: %(default)s)')
        
    group.add_argument('-e', '--t_steps_g', dest='t_steps_g', type=str,
        default=DEFAULTS['t_steps_g'], metavar='INT',
        help=HELPTEXT['t_steps_g'] + ' (default: %(default)s)')
    
    group.add_argument('-a', '--no-outlier', dest='no_outlier',
        action='store_true',
        default=DEFAULTS['no_outlier'], help=HELPTEXT['no_outlier'])
    
    group.add_argument('-w', '--v1', dest='v1',
        action='store_true',
        default=DEFAULTS['v1'], help=HELPTEXT['v1'])
    
    group.add_argument('-A', '--mzscore', dest='mzscore', type=float,
        default=DEFAULTS['mzscore'], metavar='FLOAT',
        help=HELPTEXT['mzscore'] + ' (default: %(default)s)')
    
    group.add_argument('-T', '--tolerance', dest='tolerance', type=float,
        default=DEFAULTS['tolerance'], metavar='FLOAT',
        help=HELPTEXT['tolerance'] + ' (default: %(default)s)')

    group.add_argument('-W', '--multiple', dest='multiple',
        action='store_true',
        default=DEFAULTS['multiple'], help=HELPTEXT['multiple'])    

    group.add_argument('-p', '--file-prefix', dest='prefix', type=str,
        default=DEFAULTS['prefix'], metavar='STR',
        help=HELPTEXT['prefix'] + ' (default: pdb%(default)s)')

    subparser.add_argument('pdb', help='PDB identifier or filename')

    subparser.set_defaults(func=lambda ns: prody_clustenm(ns.__dict__.pop('pdb'),
                                                          **ns.__dict__))

    subparser.set_defaults(subparser=subparser)