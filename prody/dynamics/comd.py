from prody.proteins.pdbfile import parsePDB, writePDBStream, parsePDBStream
from prody.measure.transform import calcRMSD, calcTransformation, getRMSD, applyTransformation
from prody.measure.measure import buildDistMatrix, calcDeformVector
from prody.ensemble.ensemble import Ensemble
from prody.utilities import createStringIO, importLA, checkCoords, copy

norm = importLA().norm

from prody import LOGGER

from numpy import *
import numpy as np
from multiprocessing import cpu_count, Pool
from random import random
import os.path
import sys
from decimal import Decimal, ROUND_HALF_UP

from .adaptive import ONEWAY, ALTERNATING, SERIAL, DEFAULT
from .hybrid import Hybrid
from .nma import NMA
from .gnm import GNM, ZERO
from .anm import ANM
from .sampling import traverseMode

__all__ = ['calcANMMC', 'CoMD']

def calcANMMC(initial, final, **kwargs):
    """Perform ANM-MC calculations from CoMD.
    """
    devi = kwargs.pop('devi', 0.5)
    stepcutoff = kwargs.pop('stepcutoff', 2.)
    acceptance_ratio = kwargs.pop('acceptance_ratio', 0.9) # f in the paper

    cutoff = kwargs.pop('cutoff', 15.)
    anm_cut = kwargs.pop('anm_cut', cutoff)
    log = kwargs.pop('log', True)

    N = kwargs.pop('N', 10000)
    usePseudoatoms = kwargs.pop('usePseudoatoms', False)

    initial_pdb_id = kwargs.pop('initial_pdb_id', initial.getTitle())
    original_initial_pdb = kwargs.pop('original_initial_pdb', initial_pdb_id)

    original_final_pdb = kwargs.pop('original_final_pdb', final.getTitle())

    if usePseudoatoms:
        initial_ca = initial
        final_ca = final
    else:
        initial_ca = initial.select('name CA or name BB')
        final_ca = final.select('name CA or name BB')

    n_modes = kwargs.pop('n_modes', 20)
    # ANM calculation based on current
    pdb_anm = ANM('pdb ca')
    pdb_anm.buildHessian(initial_ca, cutoff=anm_cut, **kwargs)
    pdb_anm.calcModes(n_modes=n_modes, **kwargs)

    # Cumulative sum vector preparation for metropolis sampling
    eigs = 1/sqrt(pdb_anm.getEigvals())
    eigs_n = zeros(eigs.shape)
    eigs_n = eigs / sum(eigs)
    eigscumsum = eigs_n.cumsum()
    U = pdb_anm.getEigvecs()

    # Take a step along mode 1 (ID 0) to calculate the scale factor
    pdb_ca = initial_ca
    pdb_ca_temp = pdb_ca.copy()
    ID = 0
    direction = 1.
    coords_temp = pdb_ca_temp.getCoords()
    coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID] * eigs[ID]
    coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] * eigs[ID]
    coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] * eigs[ID]
    pdb_ca_temp.setCoords(coords_temp)
    pdb_ca = pdb_ca_temp.copy()
    biggest_rmsd = calcRMSD(pdb_ca.getCoords(), initial_ca.getCoords())
    scale_factor = devi/biggest_rmsd # This means that devi is the maximum deviation in RMSD for any step

    # counts for metropolis sampling
    count1 = 0 # Up-hill moves
    count2 = 0 # Accepted up-hill moves
    count3 = 0 # Down-hill moves

    # read MC parameter from file
    if (os.path.isfile(initial_pdb_id + '_ratio.dat') and
            os.stat(initial_pdb_id + '_ratio.dat').st_size != 0):
        MCpara = loadtxt(initial_pdb_id + '_ratio.dat')
        accept_para = MCpara[4]
        if MCpara[1] > acceptance_ratio + 0.05:
            accept_para *= 1.5
        elif MCpara[1] < acceptance_ratio - 0.05:
            accept_para /= 1.5
        else:
            savetxt(initial_pdb_id + '_status.dat',[1])
    else:
        accept_para = 0.1

    # MC parameter 1 is the acceptance ratio, f, which should converge on
    # the selected value with a tolerance of 0.05 either side
    # and accept_para, gamma, is adjusted to help bring it within these limits.
    # This also happens every 5 steps during the run.

    if original_initial_pdb != original_final_pdb:
        # difference from the target structure is defined as the energy and the minimum is zero. 
        
        native_dist = buildDistMatrix(final_ca)
        dist = buildDistMatrix(initial_ca)
        Ep = sum((native_dist - dist)**2)

    # Reset pdb_ca (the current structure whole the steps back to the original)
    pdb_ca = initial_ca

    step_count = 0
    check_step_counts = [0]

    if log:
        sys.stdout.write(' '*2 + 'rmsd' + ' '*2 + 'rand' + ' '*2 + 'ID' + ' '*3 + 'step'
                         + ' '*2 + 'accept_para' + ' '*5 + 'f' + '\n')

    # MC Loop 
    for k in range(N):
        pdb_ca_temp = pdb_ca.copy()
        rand = random()
        ID = argmax(rand<eigscumsum)
        direction = 2*(random()>0.5)-1

        coords_temp = pdb_ca_temp.getCoords()
        coords_temp[0:,0] = coords_temp[0:,0] + direction * U[range(0,len(U),3),ID] * eigs[ID] * scale_factor
        coords_temp[0:,1] = coords_temp[0:,1] + direction * U[range(1,len(U),3),ID] * eigs[ID] * scale_factor
        coords_temp[0:,2] = coords_temp[0:,2] + direction * U[range(2,len(U),3),ID] * eigs[ID] * scale_factor
        pdb_ca_temp.setCoords(coords_temp)

        if original_initial_pdb != original_final_pdb:   
            dist = buildDistMatrix(pdb_ca_temp)
            En = sum((native_dist - dist)**2)

            # Check whether you are heading the right way and accept uphill moves 
            # depending on the Metropolis criterion. Classically this depends on RT 
            # but this is subsumed by the unknown units from having a uniform 
            # spring constant that is set to 1.
            if Ep > En:
                count3 += 1
                pdb_ca = pdb_ca_temp.copy()
                Ep = En
                accepted = 1

            elif exp(-(En-Ep) * accept_para) > random():
                pdb_ca = pdb_ca_temp.copy()
                count1 += 1
                count2 += 1
                Ep = En
                accepted = 1

            else:
                count1 += 1
                accepted = 0

            if count1 == 0:
                f = 1.
            else:
                f = float(count2)/float(count1)

            if (mod(k,5)==0 and not(k==0)):
                # Update of the accept_para to keep the MC para reasonable
                # See comment lines above. 
                if f > acceptance_ratio + 0.05:
                    accept_para /= 1.5
                elif f < acceptance_ratio - 0.05:
                    accept_para *= 1.5

            if accept_para < 0.001: accept_para = 0.001

        else:
            # for exploration based on one structure
            # all moves are uphill but will be accepted anyway
            pdb_ca = pdb_ca_temp.copy()
            count3 += 1
            accepted = 1
            f = 1.

        rmsd = calcRMSD(pdb_ca.getCoords(), initial_ca.getCoords())

        if log:
            sys.stdout.write('{:6.2f}'.format(rmsd) + ' ' + '{:5.2f}'.format(rand) +
                             '{:4d}'.format(ID) + '{:7d}'.format(k) + ' '*2 + str(accepted) + ' '*2 +
                             '{:5.4e}'.format(accept_para) + ' '*2 + '{:5.4f}'.format(f) + '\n')

        if rmsd > stepcutoff:
            break
        
    # Build an ensemble for writing the final structure to a dcd file
    ensemble_final = Ensemble()
    ensemble_final.setAtoms(initial_ca)
    ensemble_final.setCoords(initial_ca)
    ensemble_final.addCoordset(pdb_ca.getCoords())

    return ensemble_final, count1, count2, count3, k, accept_para, rmsd


class CoMD(Hybrid):
    '''
    This is a new version of the Collective Molecular Dynamics (CoMD) sampling algorithm [MG13]_. 
    
    It requires PDBFixer and OpenMM for performing energy minimization and MD simulations 
    in implicit/explicit solvent.

    Instantiate a ClustENM-like hybrid method object for CoMD.

    .. [MG13] Mert Gur, Jeffry D Madura, Ivet Bahar. Global transitions of proteins 
            explored by a multiscale hybrid methodology: application to adenylate kinase. 
            *Biophys J* **2013** 105:1643-1652.
    '''
    def __init__(self, title, **kwargs):
        super().__init__(title=title)
        self._atomsB = None
        self._coordsB = None
        self._indicesB = None
        
        self._defvecs = []
        self._rmsds = []
        self._traj_rmsds = []
        self._cg_ensA = Ensemble(title=title)
        self._cg_ensB = Ensemble(title=title)
        self._targeted = True

        self._tmdk = kwargs.get('tmdk', 2000) 
        # 20,000 (like paper) works for forces on CA but leads to unrealistic deformations
        # 200 (like NAMD website) is too weak for such small conformational changes

    def _sample(self, **kwargs):
        log = kwargs.pop('log', False)

        conf, conf2 = self._conformers[-2], self._conformers[-1]
        
        tmp = self._atoms.copy()
        tmpB = self._atomsB.copy()

        if self._direction == 1:
            tmp.setCoords(conf)
            tmpB.setCoords(conf2)
        else:
            tmpB.setCoords(conf)
            tmp.setCoords(conf2)

        cg = tmp[self._idx_cg]
        cgB = tmpB[self._idx_cg]

        anm_cg = self._buildANM(cg)
        if not self._checkANM(anm_cg):
            return None

        anm_cgB = self._buildANM(cgB)
        if not self._checkANM(anm_cgB):
            return None

        self._direction_mode = kwargs.pop('mode', DEFAULT)
        rmsd = self._rmsd[self._cycle]

        if self._direction_mode == ONEWAY:
            LOGGER.info('\nStarting cycle with structure A')
            self._cg_ensA, _, _, _, _, _, rmsd = calcANMMC(cg, cgB, log=log,
                                                           stepcutoff=rmsd,
                                                           n_modes=self._n_modes,
                                                           **kwargs)
            cg_ens = self._cg_ensA

        elif self._direction_mode == ALTERNATING:
            if self._direction == 1:
                LOGGER.info('\nStarting cycle with structure A')
                self._cg_ensA, _, _, _, _, _, rmsd = calcANMMC(cg, cgB, log=log,
                                                               stepcutoff=rmsd,
                                                               n_modes=self._n_modes,
                                                               **kwargs)
                cg_ens = self._cg_ensA

            else:
                LOGGER.info('\nStarting cycle with structure B')
                self._cg_ensB, _, _, _, _, _, rmsd = calcANMMC(cgB, cg, log=log,
                                                               stepcutoff=rmsd,
                                                               n_modes=self._n_modes,
                                                               **kwargs)
                cg_ens = self._cg_ensB

        elif self._direction_mode == SERIAL:
            if self._direction == 1:
                LOGGER.info('\nStarting cycle with structure A')
                self._cg_ensA, _, _, _, _, _, rmsd = calcANMMC(cg, cgB, log=log,
                                                               stepcutoff=rmsd,
                                                               n_modes=self._n_modes,
                                                               **kwargs)
                cg_ens = self._cg_ensA

            else:
                LOGGER.info('\nStarting cycle with structure B')
                self._cg_ensB, _, _, _, _, _, rmsd = calcANMMC(cgB, cg, log=log,
                                                               stepcutoff=rmsd,
                                                               n_modes=self._n_modes,
                                                               **kwargs)
                cg_ens = self._cg_ensB

        else:
            raise ValueError('unknown aANM mode: %d' % self._direction_mode)
        
        if self._direction == 1:
            defvec = calcDeformVector(cg, cg_ens.getCoordsets()[-1])
            model = NMA()
            model.setEigens(defvec.getArray().reshape((defvec.getArray().shape[0], 1)))
            model_ex = self._extendModel(model, cg, tmp)
            def_ens = traverseMode(model_ex[0], tmp, 1, rmsd)
            coordsets = [def_ens.getCoordsets()[-1]]
            
            if self._direction_mode == ALTERNATING:
                self._direction = 2
        else:
            defvec = calcDeformVector(cgB, cg_ens.getCoordsets()[-1])
            model = NMA()
            model.setEigens(defvec.getArray().reshape((defvec.getArray().shape[0], 1)))
            model_ex = self._extendModel(model, cgB, tmpB)
            def_ens = traverseMode(model_ex[0], tmpB, 1, rmsd)
            coordsets = [def_ens.getCoordsets()[-1]]
            
            if self._direction_mode == ALTERNATING:
                self._direction = 1

        if self._targeted:
            if self._parallel:
                with Pool(cpu_count()) as p:
                    pot_conf = p.map(self._multi_targeted_sim,
                                     [(conf, coords) for coords in coordsets])
            else:
                pot_conf = [self._multi_targeted_sim((conf, coords)) for coords in coordsets]

            pots, poses = list(zip(*pot_conf))

            idx = np.logical_not(np.isnan(pots))
            coordsets = np.array(poses)[idx]

            LOGGER.debug('%d/%d sets of coordinates were moved to the target' % (len(coordsets), len(poses)))

        return coordsets


    def _multi_targeted_sim(self, args):

        conf = args[0]
        coords = args[1]

        return self._targeted_sim(conf, coords, tmdk=self._tmdk)


    def _targeted_sim(self, coords0, coords1, tmdk=15.,
                      d_steps=100, n_max_steps=10000, ddtol=1e-3, n_conv=5):

        try:
            from simtk.openmm import CustomExternalForce
            from simtk.openmm.app import StateDataReporter
            from simtk.unit import nanometer, angstrom, kilocalorie_per_mole, kilojoule_per_mole
        except ImportError:
            raise ImportError('Please install PDBFixer and OpenMM in order to use Hybrid.')

        tmdk *= kilocalorie_per_mole/angstrom**2
        n_atoms = coords0.shape[0]
        tmdk /= n_atoms

        pos1 = coords1 * angstrom

        force = CustomExternalForce("tmdk*periodicdistance(x, y, z, x0, y0, z0)^2")
        force.addGlobalParameter('tmdk', 0) 
        force.addPerParticleParameter('x0')
        force.addPerParticleParameter('y0')
        force.addPerParticleParameter('z0')
        force.setForceGroup(1)

        for i, atm_idx in enumerate(np.arange(n_atoms)):
            pars = pos1[i, :].value_in_unit(nanometer)
            force.addParticle(int(atm_idx), pars)

        simulation = self._prep_sim(coords0, external_forces=[force])

        # automatic conversion into nanometer will be carried out.
        simulation.context.setPositions(coords0 * angstrom)

        dist = dist0 = calcRMSD(coords0, coords1)
        m_conv = 0
        n_steps = 0
        try:
            simulation.minimizeEnergy(tolerance=self._tolerance*kilojoule_per_mole,
                                      maxIterations=self._maxIterations)

            # update parameters
            while n_steps < n_max_steps:
                simulation.context.setParameter('tmdk', tmdk)
                force.updateParametersInContext(simulation.context)

                simulation.step(d_steps)
                n_steps += d_steps

                # evaluate distance to destination
                pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(angstrom)
                d = calcRMSD(pos, coords1)
                dd = np.abs(dist - d)

                dist = d

                if dd < ddtol:
                    m_conv += 1

                if m_conv >= n_conv:
                    break

            LOGGER.debug('RMSD: %4.2f -> %4.2f' % (dist0, dist))

            simulation.context.setParameter('tmdk', 0.0)
            simulation.minimizeEnergy(tolerance=self._tolerance*kilojoule_per_mole,
                                      maxIterations=self._maxIterations)

            pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(angstrom)
            pot = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)

            return pot, pos

        except BaseException as be:
            LOGGER.warning('OpenMM exception: ' + be.__str__() + ' so the corresponding conformer will be discarded!')

            return np.nan, np.full_like(coords0, np.nan)

    def setAtoms(self, atomsA, atomsB=None, pH=7.0, **kwargs):
        aligned = kwargs.get('aligned', False)
        if not aligned and atomsB is not None:
            T = calcTransformation(atomsA.ca, atomsB.ca, weights=atomsA.ca.getFlags("mapped"))
            _ = applyTransformation(T, atomsA)

        if self._isBuilt():
            super(Hybrid, self).setAtoms(atomsA)
            self._atomsB = atomsB
        else:
            A_B_dict = {0: 'A', 1: 'B'}
            for i, atoms in enumerate([atomsA, atomsB]):
                if i == 1 and atomsB is None:
                    break

                atoms = atoms.select('not hetatm')

                self._nuc = atoms.select('nucleotide')

                if self._nuc is not None:

                    idx_p = []
                    for c in self._nuc.getChids():
                        tmp = self._nuc[c].iterAtoms()
                        for a in tmp:
                            if a.getName() in ['P', 'OP1', 'OP2', 'OP3']:
                                idx_p.append(a.getIndex())

                    if idx_p:
                        nsel = 'not index ' + ' '.join([str(i) for i in idx_p])
                        atoms = atoms.select(nsel)

                LOGGER.info('Fixing structure {0}...'.format(A_B_dict[i]))
                LOGGER.timeit('_clustenm_fix')
                self._ph = pH
                self._fix(atoms, i)
                LOGGER.report('The structure was fixed in %.2fs.',
                            label='_clustenm_fix')

                if self._nuc is None:
                    self._idx_cg = self._atoms.ca.getIndices()
                    self._n_cg = self._atoms.ca.numAtoms()
                else:
                    self._idx_cg = self._atoms.select("name CA C2 C4' P").getIndices()
                    self._n_cg = self._atoms.select("name CA C2 C4' P").numAtoms()

                self._n_atoms = self._atoms.numAtoms()
                self._indices = None

            if i == 1:
                self._cg_ensA.setAtoms(self._atoms[self._idx_cg])
            else:
                self._cg_ensB.setAtoms(self._atomsB[self._idx_cg])

    def _fix(self, atoms, i):
        try:
            from pdbfixer import PDBFixer
            from simtk.openmm.app import PDBFile
        except ImportError:
            raise ImportError('Please install PDBFixer and OpenMM in order to use ClustENM and related hybrid methods.')

        stream = createStringIO()
        title = atoms.getTitle()
        writePDBStream(stream, atoms)
        stream.seek(0)
        fixed = PDBFixer(pdbfile=stream)
        stream.close()

        fixed.missingResidues = {}
        fixed.findNonstandardResidues()
        fixed.replaceNonstandardResidues()
        fixed.removeHeterogens(False)
        fixed.findMissingAtoms()
        fixed.addMissingAtoms()
        fixed.addMissingHydrogens(self._ph)

        stream = createStringIO()
        PDBFile.writeFile(fixed.topology, fixed.positions,
                          stream, keepIds=True)
        stream.seek(0)
        if i == 0:
            self._atoms = parsePDBStream(stream)
            self._atoms.setTitle(title)
        else:
            self._atomsB = parsePDBStream(stream)
            self._atomsB.setTitle(title)            
        stream.close()

        self._topology = fixed.topology
        self._positions = fixed.positions

    def getAtomsA(self, selected=True):
        'Returns atoms for structure A (main atoms).'
        return super(CoMD, self).getAtoms(selected)

    def getAtomsB(self, selected=True):
        'Returns atoms for structure B.'
        if self._atomsB is None:
            return None
        if self._indicesB is None or not selected:
            return self._atomsB
        return self._atomsB[self._indicesB]

    def getRMSDsB(self):
        if self._confs is None or self._coords is None:
            return None

        indices = self._indicesB
        if indices is None:
            indices = np.arange(self._confs.shape[1])
        
        weights = self._weights[indices] if self._weights is not None else None

        return calcRMSD(self._atomsB, self._confs[:, indices], weights)

    def getConvergenceRMSDs(self):
        if self._confs is None or self._coords is None:
            return None

        indices = self._indices
        if indices is None:
            indices = np.arange(self._confs.shape[1])
        
        weights = self._weights[indices] if self._weights is not None else None

        n_confs = self.numConfs()
        n_confsA = int(Decimal(n_confs/2).to_integral(rounding=ROUND_HALF_UP))

        confsA = self._confs[:n_confsA]
        if n_confs % 2:
            confsB = self._confs[n_confsA:]
        else:
            confsB = self._confs[n_confsA:]

        RMSDs = []
        for i in range(n_confsA):
            for j in range(2):
                if i + j > n_confsA - 1:
                    break
                RMSDs.append(getRMSD(confsA[i+j], confsB[n_confsA-(i+1)], weights=weights))

        return np.array(RMSDs)

    def run(self, cutoff=15., n_modes=20, gamma=1., n_confs=50, rmsd=1.0,
            n_gens=5, solvent='imp', sim=False, force_field=None, temp=303.15,
            t_steps_i=1000, t_steps_g=7500,
            outlier=True, mzscore=3.5, **kwargs):

        '''
        Performs a ClustENM-like run.

        :arg cutoff: Cutoff distance (A) for pairwise interactions used in ANM
            computations, default is 15.0 A.
        :type cutoff: float

        :arg gamma: Spring constant of ANM, default is 1.0.
        :type gamma: float

        :arg n_modes: Number of non-zero eigenvalues/vectors to calculate.
        :type n_modes: int

        :arg n_confs: Number of new conformers to be generated based on any conformer
            from the previous generation, default is 50.
        :type n_confs: int
            
        :arg rmsd: Average RMSD of the new conformers with respect to the conformer
            from which they are generated, default is 1.0 A.
            A tuple of floats can be given, e.g. (1.0, 1.5, 1.5) for subsequent generations.
            Note: In the case of ClustENMv1, this value is the maximum rmsd, not the average.
        :type rmsd: float, tuple of floats

        :arg n_gens: Number of generations.
        :type n_gens: int

        :arg solvent: Solvent model to be used. If it is set to 'imp' (default),
            implicit solvent model will be used, whereas 'exp' stands for explicit solvent model.
            Warning: In the case of nucleotide chains, explicit solvent model is automatically set.
        :type solvent: str

        :arg padding: Padding distance to use for solvation. Default is 1.0 nm.
        :type padding: float

        :arg ionicStrength: Total concentration of ions (both positive and negative) to add.
            This does not include ions that are added to neutralize the system.
            Default concentration is 0.0 molar.
        :type ionicStrength: float

        :arg force_field: Implicit solvent force field is ('amber99sbildn.xml', 'amber99_obc.xml'). 
            Explicit solvent force field is ('amber14-all.xml', 'amber14/tip3pfb.xml').
            Experimental feature: Forcefields already implemented in OpenMM can be used. 
        :type force_field: a tuple of str
        
        :arg tolerance: Energy tolerance to which the system should be minimized, default is 10.0 kJ/mole.
        :type tolerance: float
        
        :arg maxIterations: Maximum number of iterations to perform during energy minimization.
            If this is 0 (default), minimization is continued until the results converge without
            regard to how many iterations it takes.
        :type maxIterations: int

        :arg sim: If it is True (default), a short MD simulation is performed after energy minimization.
            Note: There is also a heating-up phase until the desired temperature is reached.
        :type sim: bool

        :arg temp: Temperature at which the simulations are conducted, default is 303.15 K.
        :type temp: float

        :arg t_steps_i: Duration of MD simulation (number of time steps) for the starting structure
            following the heating-up phase, default is 1000. Each time step is 2.0 fs.
            Note: Default value reduces possible drift from the starting structure. 
        :type t_steps_i : int

        :arg t_steps_g: Duration of MD simulations (number of time steps) to run for each conformer
            following the heating-up phase, default is 7500. Each time step is 2.0 fs.
            A tuple of int's can be given, e.g. (3000, 5000, 7000) for subsequent generations.
        :type t_steps_g: int or tuple of int's

        :arg outlier: Exclusion of conformers detected as outliers in each generation.
            Default is True for implicit solvent. Outliers, if any, are detected by
            the modified z-scores of the conformers' potential energies over a generation.
            Note: It is automatically set to False when explicit solvent model is being used.
        :type outlier: bool

        :arg mzscore: Modified z-score threshold to label conformers as outliers. Default is 3.5.
        :type mzscore: float

        :arg v1: Original sampling method with complete enumeration of desired ANM modes is used.
            Default is False. Maximum number of modes should not exceed 5 for efficiency.
        :type v1: bool

        :arg platform: Architecture on which the OpenMM part runs, default is None.
            It can be chosen as 'CUDA', 'OpenCL' or 'CPU'.
            For efficiency, 'CUDA' or 'OpenCL' is recommended.
        :type platform: str

        :arg parallel: If it is True (default is False), conformer generation will be parallelized.
        :type parallel: bool

        :arg min_rmsd_diff: Difference between *rmsds* from previous step to current for checking 
            convergence. Default is 0.6 A
        :type min_rmsd_diff: float

        :arg target_rmsd: Target RMSD for stopping.
            Default is 0.125 A
        :type target_rmsd: float
        '''

        if self._isBuilt():
            raise ValueError('CoMD ensemble has been built; please start a new instance')

        # set up parameters
        self._cutoff = cutoff
        self._n_modes0 = self._n_modes = n_modes
        self._gamma = gamma
        self._n_confs = n_confs
        self._rmsd = (0.,) + rmsd if isinstance(rmsd, tuple) else (0.,) + (rmsd,) * n_gens
        self._n_gens = n_gens
        self._platform = kwargs.pop('platform', None)
        self._parallel = kwargs.pop('parallel', False)
        self._targeted = kwargs.pop('targeted', self._targeted)
        self._tmdk = kwargs.pop('tmdk', self._tmdk)

        self._target_rmsd = kwargs.pop('target_rmsd', 0.6)
        self._min_rmsd_diff = kwargs.pop('min_rmsd_diff', 0.125)

        self._direction_mode = kwargs.get('mode', DEFAULT)
        if self._direction_mode == ALTERNATING and self._n_gens % 2:
            raise ValueError('ALTERNATING modes needs even n_gens')

        self._direction = 1

        self._sol = solvent if self._nuc is None else 'exp'
        self._padding = kwargs.pop('padding', 1.0)
        self._ionicStrength = kwargs.pop('ionicStrength', 0.0)
        if self._sol == 'imp':
            self._force_field = ('amber99sbildn.xml', 'amber99_obc.xml') if force_field is None else force_field
        if self._sol == 'exp':
            self._force_field = ('amber14-all.xml', 'amber14/tip3pfb.xml') if force_field is None else force_field
        self._tolerance = kwargs.pop('tolerance', 10.0)
        self._maxIterations = kwargs.pop('maxIterations', 0)
        self._sim = sim
        self._temp = temp

        if self._sim:
            if isinstance(t_steps_g, tuple):
                self._t_steps = (t_steps_i,) + t_steps_g
            else:
                self._t_steps = (t_steps_i,) + (t_steps_g,) * n_gens

        self._outlier = False if self._sol == 'exp' else outlier
        self._mzscore = mzscore
        self._v1 = kwargs.pop('v1', False)

        self._cycle = 0

        # check for discontinuity in the structure A
        gnm = GNM()
        gnm.buildKirchhoff(self._atoms[self._idx_cg], cutoff=self._cutoff)
        K = gnm.getKirchhoff()
        rank_diff = (len(K) - 1
                     - np.linalg.matrix_rank(K, tol=ZERO, hermitian=True))
        if rank_diff != 0:
            raise ValueError('atoms has disconnected parts; please check the structure')

        LOGGER.timeit('_hybrid_overall')

        LOGGER.info('Generation 0 for structure A ...')

        if self._sim:
            if self._t_steps[0] != 0:
                LOGGER.info('Minimization, heating-up & simulation in generation 0 for structure A ...')
            else:
                LOGGER.info('Minimization & heating-up in generation 0 for structure A ...')
        else:
            LOGGER.info('Minimization in generation 0 for structure A ...')
        LOGGER.timeit('_clustenm_min')
        potential, conformer = self._min_sim(self._atoms.getCoords())
        if np.isnan(potential):
            raise ValueError('Initial structure A could not be minimized. Try again and/or check your structure.')

        LOGGER.report(label='_clustenm_min')

        LOGGER.info('#' + '-' * 19 + '/*\\' + '-' * 19 + '#')

        self.setCoords(conformer)

        potentials = [potential]
        sizes = [1]
        new_shape = [1]
        for s in conformer.shape:
            new_shape.append(s)
        conf = conformer.reshape(new_shape)
        self._conformers = start_confs = conf
        keys = [(0, 0)]

        # check for discontinuity in the structure B
        gnmB = GNM()
        gnmB.buildKirchhoff(self._atomsB[self._idx_cg], cutoff=self._cutoff)
        KB = gnmB.getKirchhoff()
        rank_diffB = (len(KB) - 1
                     - np.linalg.matrix_rank(KB, tol=ZERO, hermitian=True))
        if rank_diffB != 0:
            raise ValueError('atoms B has disconnected parts; please check the structure')

        LOGGER.info('Generation 0 for structure B ...')

        if self._sim:
            if self._t_steps[0] != 0:
                LOGGER.info('Minimization, heating-up & simulation in generation 0 for structure B ...')
            else:
                LOGGER.info('Minimization & heating-up in generation 0 for structure B ...')
        else:
            LOGGER.info('Minimization in generation 0 for structure B ...')
        LOGGER.timeit('_clustenm_min')
        potentialB, conformerB = self._min_sim(self._atomsB.getCoords())
        if np.isnan(potentialB):
            raise ValueError('Initial structure B could not be minimized. Try again and/or check your structure.')

        LOGGER.report(label='_clustenm_min')

        LOGGER.info('#' + '-' * 19 + '/*\\' + '-' * 19 + '#')

        self.setCoordsB(conformerB)

        potentials.extend([potentialB])
        sizes.extend([1])
        new_shape = [1]
        for s in conformerB.shape:
            new_shape.append(s)
        confB = conformerB.reshape(new_shape)
        start_confsB = confB
        self._conformers = np.vstack((self._conformers, start_confsB))
        keys.extend([(0, 0)])
        
        curr_rmsd = calcRMSD(self._conformers[-1], self._conformers[-2])
        self._traj_rmsds.append(curr_rmsd)

        for i in range(1, self._n_gens+1):
            self._cycle += 1
            LOGGER.info('Generation %d ...' % i)

            confs, weights = self._generate(**kwargs)

            LOGGER.timeit('_clustenm_min_sim')

            pot_conf = [self._min_sim(conf) for conf in confs]

            LOGGER.report('Structures were sampled in %.2fs.',
                        label='_clustenm_min_sim')
            LOGGER.info('#' + '-' * 19 + '/*\\' + '-' * 19 + '#')

            pots, confs = list(zip(*pot_conf))
            idx = np.logical_not(np.isnan(pots))
            weights = np.array(weights)[idx]
            pots = np.array(pots)[idx]
            confs = np.array(confs)[idx]

            if self._outlier:
                idx = np.logical_not(self._outliers(pots))
            else:
                idx = np.full(pots.size, True, dtype=bool)

            sizes.extend(weights[idx])
            potentials.extend(pots[idx])

            start_confs = self._superpose_cg(confs[idx])

            for j in range(start_confs.shape[0]):
                keys.append((i, j))
            self._conformers = np.vstack((self._conformers, start_confs))

            curr_rmsd = calcRMSD(self._conformers[-1], self._conformers[-2])
            self._traj_rmsds.append(curr_rmsd)

            diff_rmsds = self._traj_rmsds[-2] - self._traj_rmsds[-1]

            if self._cycle == 1:
                g = "generation"
            else:
                g = "generations"

            LOGGER.info('\nRan {:d} {:s}, RMSD {:5.2f}, change in RMSD {:5.2f}\n'.format(self._cycle, g,
                                                                                         curr_rmsd,
                                                                                         diff_rmsds))

            if curr_rmsd < self._target_rmsd or diff_rmsds < self._min_rmsd_diff:
                if self._direction_mode == 2 and self._direction == 1:
                    self._direction = 2
                else:
                    LOGGER.report('Transition converged in %.2fs.', '_hybrid_overall')
                    break

        LOGGER.timeit('_hybrid_ens')
        LOGGER.info('Creating an ensemble of conformers ...')

        self._build(self._conformers, keys, potentials, sizes)
        LOGGER.report('Ensemble was created in %.2fs.', label='_hybrid_ens')

        self._time = LOGGER.timing(label='_hybrid_overall')
        LOGGER.report('All completed in %.2fs.', label='_hybrid_overall')

    def _generate(self, **kwargs):
        LOGGER.info('Sampling conformers in generation %d ...' % self._cycle)
        LOGGER.timeit('_clustenm_gen')

        tmp = [self._sample(**kwargs)]

        tmp = [r for r in tmp if r is not None]

        confs_ex = np.concatenate(tmp)

        return confs_ex, [1]

    def setCoordsB(self, coords):
        """Set *coords* as the ensemble reference coordinate set.  *coords*
        may be an array with suitable data type, shape, and dimensionality, or
        an object with :meth:`getCoords` method."""

        atoms = coords
        try:
            if isinstance(coords, Ensemble):
                coords = copy(coords._coords)
            else:
                coords = coords.getCoords()
        except AttributeError:
            pass
        finally:
            if coords is None:
                raise ValueError('coordinates of {0} are not set'
                                 .format(str(atoms)))

        try:
            checkCoords(coords, natoms=self._n_atoms)
        except TypeError:
            raise TypeError('coords must be a numpy array or an object '
                            'with `getCoords` method')

        if coords.shape == self._coords.shape:
            self._coordsB = coords
            self._n_atomsB = coords.shape[0]

            if isinstance(atoms, Ensemble):
                self._indicesB = atoms._indices
                self._atomsB = atoms._atoms
        else:
            raise ValueError('coordsB must have the same shape as main coords')

    def getCoordsB(self, selected=True):
        """Returns a copy of reference coordinates for selected atoms."""

        if self._coordsB is None:
            return None
        if self._indices is None or not selected:
            return self._coordsB.copy()
        return self._coordsB[self._indices].copy()


if __name__ == '__main__':

    from prody import *
    from numpy import *
    import time

    time.sleep(10)
    ar = []
    for arg in sys.argv:
        ar.append(arg)

    if len(ar) > 1:
        initial_pdbn=ar[1]
    else:
        raise ValueError('Please provide at least 1 argument (a PDB filename)')

    if len(ar) > 2:
        final_pdbn=ar[2]
    else:
        final_pdbn = initial_pdbn
        
    initial_pdb_id = initial_pdbn[:initial_pdbn.rfind('.')]
    final_pdb_id = final_pdbn[:final_pdbn.rfind('.')]

    if len(ar) > 3 and ar[3].strip() != '0':
        original_initial_pdb = ar[3]
    else:
        original_initial_pdb = initial_pdb_id

    if len(ar) > 4 and ar[4].strip() != '0':
        original_final_pdb = ar[4]
    else:
        original_final_pdb = final_pdb_id

    if len(ar) > 5 and ar[5].strip() != '0':
        comd_cycle_number = ar[5]
    else:
        comd_cycle_number = 1

    if len(ar) > 6 and ar[6].strip() != '0':
        devi = float(ar[6])
    else:
        devi = 0.5

    if len(ar) > 7 and ar[7].strip() != '0':
        stepcutoff=float(ar[7])
    else:
        stepcutoff=2.

    if len(ar) > 8 and ar[8].strip() != '0':
        acceptance_ratio = float(ar[8])
    else:
        acceptance_ratio = 0.9

    if len(ar) > 9 and ar[9].strip() != '0':
        anm_cut=float(ar[9])
    else:
        anm_cut=15

    if len(ar) > 10 and ar[10].strip() != '0':
        N=int(ar[10])
    else:
        N=10000

    if len(ar) > 11 and ar[11].strip() != '0':
        final_structure_dcd_name = ar[11]
    else:
        final_structure_dcd_name = 'cycle_{0}_'.format(int(comd_cycle_number)) + \
                                    initial_pdb_id + '_' + final_pdb_id + '_final_structure.dcd'

    if len(ar) > 12 and ar[12].strip() != '0':
        usePseudoatoms = int(ar[12])
    else:
        usePseudoatoms = 0

    initial_pdb = parsePDB(initial_pdbn)
    final_pdb = parsePDB(final_pdbn)

    ensemble_final, count1, count2, count3, k, accept_para, rmsd = calcANMMC(initial_pdb, final_pdb,
                                                                             initial_pdb_id=initial_pdb_id,
                                                                             original_initial_pdb=original_initial_pdb,
                                                                             original_final_pdb=original_final_pdb,
                                                                             comd_cycle_number=comd_cycle_number,
                                                                             devi=devi, stepcutoff=stepcutoff,
                                                                             acceptance_ratio=acceptance_ratio,
                                                                             anm_cut=anm_cut, N=N,
                                                                             usePseudoatoms=usePseudoatoms)
    writeDCD(final_structure_dcd_name, ensemble_final)

    ratios = [count2/N, count2/count1 if count1 != 0 else 0, count2, k, accept_para ]
    savetxt(initial_pdb_id + '_ratio.dat', ratios, fmt='%.2e')

