"""Extract a selection of atoms from a PDB file."""

from ..apptools import *

__all__ = ['prody_select']

def prody_energy(*pdbs, **kwargs):
    """Write out potential energy in kJ/mol from OpenMM using ClustENM wrapper.

    :arg pdbs: PDB identifier(s) or filename(s)

    :arg output: output filename, default is :file:`pdb_selected.pdb`

    :arg prefix: prefix for output file, default is PDB filename

    :arg suffix: output filename suffix, default is :file:`_selected`
    
    :arg solvent: "imp" for implicit solvent (default) or "exp" for explicit solvent
    
    :arg padding: padding in nanometers
    
    :arg force_field_protein: name of force field for protein in OpenMM
        default depends on solvent type as in ClustENM
    
    :arg force_field_sol: name of force field for solvent in OpenMM 
        default depends on solvent type as in ClustENM
        
    :arg model: model to analyse. Default is use all
    
    :arg minimise: whether to energy minimise
    
    :arg select: atom selection string, default is "all"
    """

    from os.path import isfile
    from prody import LOGGER, parsePDB, ClustENM

    if not pdbs:
        raise ValueError('pdb argument must be provided')

    prefix = kwargs.get('prefix', None)
    suffix = kwargs.get('suffix', '_energy')
    output = kwargs.get('output', None)
    altloc = kwargs.get('altloc', None)
    
    sol = kwargs.get('solvent', 'imp')
    if sol not in ["imp", "exp"]:
        LOGGER.warn('Solvent {0} did not match "imp" or "exp".'
                    .format(repr(sol)))
        return
        
    padding = kwargs.get('padding', 1.0)
    
    force_field_protein = kwargs.get('force_field_protein', None)
    force_field_solvent = kwargs.get('force_field_solvent', None)
    
    model = kwargs.get('model', 0)
    if model == 0:
        model = None
        
    minimise = kwargs.get('minimise', False)
    selstr = kwargs.get('select', 'all')
    
    force_field = (force_field_protein, force_field_solvent)
    if force_field == (None, None):
        force_field = None

    for pdb in pdbs:
        ag = parsePDB(pdb, model=model)

        outname = output or ((prefix or ag.getTitle()) + suffix)
        if outname[-4] != '.':
            outname += '.txt'        
        f = open(outname, 'w')
        
        for i in range(ag.numCoordsets()):
            
            model_num = i + 1
            
            LOGGER.info('\nAnalysing model {0} from {1}'.format(model_num, pdb))
            
            ag = parsePDB(pdb, model=model_num)

            clu = ClustENM()
            clu.setAtoms(ag)
            
            clu._sol = sol
            
            if clu._sol == 'imp':
                clu._force_field = ('amber99sbildn.xml', 'amber99_obc.xml') if force_field is None else force_field
            if clu._sol == 'exp':
                clu._force_field = ('amber14-all.xml', 'amber14/tip3pfb.xml') if force_field is None else force_field
                
            clu._padding = padding
            
            simulation = clu._prep_sim(clu._atoms.getCoords())
            
            if minimise:
                from openmm.unit import kilojoule_per_mole
                simulation.minimizeEnergy(tolerance=10.0 * kilojoule_per_mole, maxIterations=0)
                
            state = simulation.context.getState(getEnergy=True)
            energy = state.getPotentialEnergy()._value
            
            f.write(str(energy) + " kJ/mol\n")
        
        f.close()
        LOGGER.info('\nEnergy is written into: ' + outname + '\n')



def addCommand(commands):

    subparser = commands.add_parser('energy',
    help='fix missing atoms, solvate and write a file containing energy')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)

    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')
    
    subparser.add_argument('-m', '--model', dest='model', type=int,
        default=0, metavar='INT',
        help=('index of model that will be used in the calculations (default: all of them)'))

    subparser.set_defaults(usage_example=
    """This command selects specified atoms and writes them in a PDB file.

Fetch PDB files 1p38 and 1r39 and write backbone atoms in a file:

  $ prody select backbone 1p38 1r39""",
    test_examples=[0])


    group = subparser.add_argument_group('output options')

    group.add_argument('-o', '--output', dest='output', metavar='STR',
        type=str, help='output PDB filename (default: pdb_energy.txt)')

    group.add_argument('-p', '--prefix', dest='prefix', metavar='STR',
        type=str, help=('output filename prefix (default: PDB filename)'))

    group.add_argument('-x', '--suffix', dest='suffix', metavar='STR',
        type=str, default='_energy',
        help=('output filename suffix (default: %(default)s)'))

    subparser.add_argument('pdb', nargs='+',
        help='PDB identifier(s) or filename(s)')

    subparser.set_defaults(func=lambda ns: prody_energy(*ns.pdb,
                                                        **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
    
    group_energy = subparser.add_argument_group('energy options')

    group_energy.add_argument('-s', '--solvent', 
                              dest='solvent', metavar='STR',
                              type=str, default="imp",
            help=('name of force field for protein in OpenMM (default: %(default)s)'))
    
    group_energy.add_argument('-P', '--force_field_protein', 
                              dest='force_field_protein', metavar='STR',
                              type=str, default=None,
            help=('name of force field for protein in OpenMM (default: ClustENM default)'))

    group_energy.add_argument('-S', '--force_field_sol', 
                              dest='force_field_sol', metavar='STR',
                              type=str, default=None,
            help=('name of force field for solvent in OpenMM (default: ClustENM default)'))

    group_energy.add_argument('-M', '--minimise', dest='minimise', action='store_true',
        default=False, help=('whether to energy minimise (default: %(default)s)'))
    