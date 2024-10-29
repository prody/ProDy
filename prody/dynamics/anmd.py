# -*- coding: utf-8 -*-
'''
Copyright (c) 2020-2022 Burak Kaynak, Pemra Doruker.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

__author__ = 'Anupam Banerjee'
__credits__ = ['James Krieger']
__email__ = ['anupam.banerjee@stonybrook.edu', 'jamesmkrieger@gmail.com']

from numbers import Number
import os

from prody import LOGGER
from prody.atomic.atomic import Atomic
from prody.ensemble.ensemble import Ensemble
from prody.proteins.pdbfile import parsePDB, writePDB

from prody.dynamics.anm import ANM
from prody.dynamics.clustenm import ClustENM
from prody.dynamics.editing import extendModel
from prody.dynamics.sampling import traverseMode


__all__ = ['runANMD']

def runANMD(atoms, num_modes=2, max_rmsd=2., num_steps=2,
            tolerance=10.0, solvent='imp', force_field=None):
    """Runs the ANMD hybrid simulation method ([CM22]_), which generates conformations
    along single modes and minimises them. The first non-zero mode is scaled to *max_rmsd*
    and the remaining modes are scaled accordingly using their eigenvalues.

    :arg atoms: an object with atom and coordinate data
    :type atoms: :class:`.Atomic`

    :arg num_modes: number of modes to analyse
    :type num_modes: int

    :arg max_rmsd: maximum rmsd for non-zero mode 1
    :type max_rmsd: float

    :arg num_steps: number of conformers in each direction for each mode
    :type num_steps: int

    :arg tolerance: tolerance for energy minimisation in OpenMM 
        in kilojoule_per_mole. Default is 10
    :type tolerance: float

    :arg solvent: Solvent model to be used. If it is set to 'imp' (default),
        implicit solvent model will be used, whereas 'exp' stands for explicit solvent model.
        Warning: In the case of nucleotide chains, explicit solvent model is automatically set.
    :type solvent: str

    :arg force_field: Implicit solvent force field is ('amber99sbildn.xml', 'amber99_obc.xml'). 
        Explicit solvent force field is ('amber14-all.xml', 'amber14/tip3pfb.xml').
        Experimental feature: Forcefields already implemented in OpenMM can be used. 
    :type force_field: tuple of strings

    .. [CM22] Mary Hongying Cheng, James M Krieger, Anupam Banerjee, Yufei Xiang, 
            Burak Kaynak, Yi Shi, Moshe Arditi, Ivet Bahar. 
            Impact of new variants on SARS-CoV-2 infectivity and neutralization: 
            A molecular assessment of the alterations in the spike-host protein 
            interactions. *iScience* **2022** 25(3):103939.
    """
    try:
        from simtk.openmm.app import PDBFile, ForceField, \
            Simulation, HBonds, NoCutoff
        from simtk.openmm import LangevinIntegrator
        from simtk.unit import nanometer, kelvin, picosecond, picoseconds, \
            angstrom
    except ImportError:
        raise ImportError('Please install PDBFixer and OpenMM to use ANMD')

    if not isinstance(atoms, Atomic):
        raise TypeError('atoms should be an Atomic object')

    if not isinstance(num_modes, int):
        raise TypeError('num_modes should be an integer')
    
    if not isinstance(max_rmsd, Number):
        raise TypeError('max_rmsd should be a float')

    if not isinstance(tolerance, Number):
        raise TypeError('tolerance should be a float')
    
    if not isinstance(solvent, str):
        raise TypeError('solvent should be a string')
    if solvent not in ['imp', 'exp']:
        raise ValueError('solvent should be "imp" or "exp"')
    
    if not isinstance(force_field, (type(None), tuple)):
        raise TypeError('force_field should be None or a tuple')
    if force_field is not None:
        if len(force_field) != 2:
            raise ValueError('force_field tuple should have 2 elements')
        if not isinstance(force_field[0], str) or not isinstance(force_field[1], str):
            raise TypeError('force_field tuple entries should be strings')

    pdb_name=atoms.getTitle().replace(' ', '_')

    fix_name = pdb_name + '_fixed.pdb'
    if os.path.exists(fix_name):
        LOGGER.info('\nFixed structure found')
    else:
        clustenm=ClustENM()
        clustenm.setAtoms(atoms)
        clustenm.writePDBFixed()
    pdb_fix = PDBFile(fix_name)

    fixmin_name=pdb_name + '_fixedmin.pdb'
    if os.path.exists(fixmin_name):
        LOGGER.info('\nMinimised fixed structure found')
    else:
        LOGGER.info('\nMinimising fixed structure ...')
        LOGGER.timeit('_anmd_min')

        forcefield = ForceField("amber99sbildn.xml", "amber99_obc.xml")
        system = forcefield.createSystem(pdb_fix.topology, 
                                         nonbondedMethod=NoCutoff,
                                         nonbondedCutoff=1*nanometer, 
                                         constraints=HBonds)
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 
                                        0.003*picoseconds)
        simulation = Simulation(pdb_fix.topology, system, integrator)
        simulation.context.setPositions(pdb_fix.positions)
        simulation.minimizeEnergy(tolerance=tolerance)
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open(fixmin_name, 'w'))
        LOGGER.report('The fixed structure was minimised in %.2fs.\n',
                        label='_anmd_min')

    pdb_fixed=parsePDB(fixmin_name, compressed=False)
    calphas=pdb_fixed.select('calpha')
    anm=ANM()
    anm.buildHessian(calphas)
    anm.calcModes(n_modes=num_modes)
    anm_ex, atoms_all = extendModel(anm, calphas, pdb_fixed)
    eval_0=anm[0].getEigval()

    ensembles = []
    for i in range(num_modes):
        ip1 = i+1
        num_confs = 2*num_steps+1

        LOGGER.info('\nGenerating {0} conformers for mode {1} ...'.format(num_confs, ip1))
        eval_i=anm[i].getEigval()
        sc_rmsd=((1/eval_i)**0.5/(1/eval_0)**0.5)*max_rmsd
        traj_aa=traverseMode(anm_ex[i], atoms_all, n_steps=num_steps, rmsd=sc_rmsd)
        traj_aa.setAtoms(atoms_all)
       
        target_ensemble = Ensemble('mode {0} ensemble'.format(ip1))
        target_ensemble.setAtoms(atoms_all)
        target_ensemble.setCoords(atoms_all)
        
        for j in range(num_confs):
            jp1 = j+1
            writePDB('temp1.pdb', traj_aa[j])
            pdb = PDBFile('temp1.pdb')
            os.remove("temp1.pdb")

            LOGGER.info('\nMinimising structure {0} along mode {1} ...'.format(jp1, ip1))
            LOGGER.timeit('_anmd_min')
            forcefield = ForceField("amber99sbildn.xml", "amber99_obc.xml")
            system = forcefield.createSystem(pdb.topology, 
                                             nonbondedMethod=NoCutoff,
                                             nonbondedCutoff=1*nanometer, 
                                             constraints=HBonds)
            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, \
                                            0.003*picoseconds)
            simulation = Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy(tolerance=tolerance)
            positions = simulation.context.getState(getPositions=True).getPositions(
                asNumpy=True).value_in_unit(angstrom)[:pdb.topology.getNumAtoms()]
            
            target_ensemble.addCoordset(positions)

            LOGGER.report('The structure was minimised in %.2fs.',
                            label='_anmd_min')
            
        ensembles.append(target_ensemble)

    os.remove(fix_name)
    os.remove(fixmin_name)

    return ensembles


if __name__=='__main__':
    import sys
    from prody.tests.datafiles import pathDatafile
    pdb_name=sys.argv[1] if len(sys.argv) > 1 else pathDatafile('1ubi')
    num_modes=int(sys.argv[2]) if len(sys.argv) > 2 else 2
    max_rmsd=float(sys.argv[3]) if len(sys.argv) > 3 else 2.
    num_steps = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    tol=float(sys.argv[4]) if len(sys.argv) > 4 else 10.

    pdb_name_ext = pdb_name
    if pdb_name_ext.endswith('.pdb'):
        pdb_name = pdb_name_ext[:-4]
    else:
        pdb_name_ext += '.pdb'

    pdb = parsePDB(pdb_name_ext, compressed=False)

    x = runANMD(pdb, num_modes, max_rmsd, num_steps, tol)

    pdb_basename = os.path.basename(pdb_name)
    for ens in x:
        filename = pdb_basename + '_' + ens.getTitle().replace(' ', '_')
        LOGGER.info('writing PDB file {0}'.format(filename))
        writePDB(filename, ens)
