"""
Copyright (c) 2020 Burak Kaynak, Pemra Doruker.

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
"""

__author__ = 'Burak Kaynak'
__credits__ = ['Pemra Doruker', 'She Zhang']
__email__ = ['burak.kaynak@pitt.edu', 'doruker@pitt.edu', 'shz66@pitt.edu']

from itertools import product
from multiprocessing import cpu_count, Pool
from collections import OrderedDict
from os import chdir, mkdir
from os.path import isdir
import pickle
from shutil import rmtree
from sys import stdout
#from time import perf_counter

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.stats import median_absolute_deviation

from prody import LOGGER
from .anm import ANM
from .gnm import GNM, ZERO
from .editing import extendModel
from .sampling import sampleModes
from prody.atomic import AtomGroup
from prody.measure import calcTransformation, applyTransformation
from prody.ensemble import Ensemble
from prody.proteins import writePDB, parsePDB, writePDBStream, parsePDBStream
from prody.utilities import createStringIO

__all__ = ['ClustENM']

class ClustENM(Ensemble):
    '''
    ClustENMv2 is the new version of ClustENM(v1) conformation sampling algorithm [KZ16]_.
    This ANM-based hybrid algorithm requires PDBFixer and OpenMM for performing energy minimization and MD simulations in implicit solvent.
    It is Python 3.7 compatible and has been only tested on Linux machines.

    .. [KZ16] Kurkcuoglu Z, Bahar I, Doruker P. ClustENM: ENM-based sampling of essential conformational space at full atomic resolution. 
       *J Chem* **2016** 12(9):4549-4562.
    
    .. [PE17] Eastman P, Swails J, Chodera JD, McGibbon RT, Zhao Y, Beauchamp KA, Wang LP, Simmonett AC, Harrigan MP, Stern CD, Wiewiora RP, 
       Brooks BR, Pande VS. OpenMM 7: Rapid Development of High Performance Algorithms for Molecular Dynamics. *PLoS Comput Biol* **2017** 
       13:e1005659.

    Instantiate a ClustENM object.

    Parameters
    ----------
    pdb : str
        pdb name without .pdb to download the structure or with .pdb to read it from disk.
        This will be used for the initial conformer.
    chain : str (optional)
            Chain Ids. If None, all chains in the PDB file are parsed.
            Otherwise, a set of chains is parsed, e. g. 'AC'.
    pH : float
        the pH based on which to select protonation states, default is 7.0.
    cutoff : float
            Cutoff distance (A) for pairwise interactions used in ANM computations, default is 15.0 A.
    n_modes : int
            Number of non-zero eigenvalues/vectors to calculate.
    n_confs : int
            Number of new conformations to be generated based on each conformer coming from the previous generation, default is 50.
    rmsd : float, tuple of floats
        Average RMSD of the new conformers with respect to the original conformation, default is 1.0 A.
        A tuple of floats can be given, e.g. (1.0, 1.5, 1.5) for subsequent generations.
        Note: In the case of ClustENMv1, this value is the maximum rmsd, not the average.
    n_gens : int
            Number of generations.
    maxclust : int or a tuple of int's.
            The maximum number of clusters for each generation. Default in None.
            A tuple of int's can be given, e.g. (10, 30 ,50) for subsequent generations.
            Warning: Either threshold or maxclust should be given! For large number of generations and/or structures,
            specifying maxclust is more efficient.
    threshold : float or tuple of floats.
                If it is true, a short MD simulation will be performed after energy minimization. Default is True.
                Note: There is a heating-up phase until the desired temperature is reached before MD simulation starts.
    sim : bool
        If it is true, a short MD simulation will be performed after energy minimization. Default is True.
        Note: There is a heating-up phase until the desired temperature is reached before MD simulation starts.
    temp : float.
        Temperature at which the simulations are conducted. Default is 300.0 K.
    t_steps_i : int
                Duration of MD simulation (number of time steps) for the initial conformer, default is 1000.
                Note: Each time step is 2.0 fs.
    t_steps_g : int or tuple of int's
                Duration of MD simulations (number of time steps) to run for each conformer, default is 7500.
                A tuple of int's for subsequent generations, e.g. (3000, 5000, 7000).
                Note: Each time step is 2.0 fs.
    outlier : bool
            Exclusion of conformers detected as outliers in each generation,
            based on their potential energy using their modified z-scores over the potentials in that generation,         
            default is True.
    mzscore : float
            Modified z-score threshold to label conformers as outliers. Default is 3.5.
    v1 : bool
        The sampling method used in the original article is utilized, a complete enumeration of desired ANM Modes.
        Note: Maximum number of modes should not exceed 5 for efficiency.
    platform: str
            The architecture on which the OpenMM part runs, default is None. It can be chosen as 'OpenCL' or 'CPU'.
            For efficiency, 'CUDA' or 'OpenCL' is recommended.
    missing_residues : bool
        whether fixes missing residues. Default is **True**.
    '''

    def __init__(self, title=None):

        self._atoms = None
        
        self._ph = 7.0

        self._cutoff = 15.
        self._n_modes = 3
        self._n_confs = 50
        self._rmsd = (0.,) 
        self._n_gens = 5

        self._maxclust = None
        self._threshold = None

        self._sim = True
        self._temp = 300
        self._t_steps = None

        self._outlier = True
        self._mzscore = 3.5
        self._v1 = False
        self._platform = None 

        self._topology = None
        self._positions = None
        self._idx_ca = None
        self._n_ca = None
        self._cycle = 0
        # self._counts = OrderedDict()   # possibly deprecated
        # self._potentials = OrderedDict()
        # self._conformers = OrderedDict()
        self._indexer = None

        super(ClustENM, self).__init__('Unknown') # dummy title; will be replaced in next line
        self._title = title

    def __getitem__(self, index):
        if isinstance(index, tuple):
            I = self._slice(index)
            if len(I) == 0:
                raise IndexError('index out of range (%d, %d)'%index)
            index = I
        
        return super(ClustENM, self).__getitem__(index)

    def getAtoms(self):
        return self._atoms

    def isBuilt(self):
        return self._confs is not None

    def setAtoms(self, atoms, pH=7.0, fix_missing_residues=True):
        if self.isBuilt():
            super(ClustENM, self).setAtoms(atoms)
        else:
            LOGGER.info('Fixing the structure ...')
            LOGGER.timeit('_clustenm_fix')
            self._ph = pH
            self._fix(atoms, fix_missing_residues)
            LOGGER.report('The structure was fixed in %.2fs.', label='_clustenm_fix')

            # self._atoms must be an AtomGroup instance because of self._fix().
            self._idx_ca = self._atoms.ca.getIndices()
            self._n_ca = self._atoms.ca.numAtoms()
            self._n_atoms = self._atoms.numAtoms()
            self._indices = None

            # check for discontinuity in the structure
            gnm = GNM()
            gnm.buildKirchhoff(self._atoms.ca)
            K = gnm.getKirchhoff()
            rank_diff = (len(K) - 1
                        - np.linalg.matrix_rank(K, tol=ZERO, hermitian=True))
            if rank_diff != 0:
                raise ValueError('atoms has disconnected parts; please check the structure')

    def getTitle(self):
        title = 'Unknown'
        if self._title is None:
            atoms = self.getAtoms()
            if atoms is not None:
                title = atoms.getTitle() + '_clustenm'
        else:
            title = self._title

        return title

    def setTitle(self, value):
        if not isinstance(value, str) and value is not None:
            raise TypeError('title must be either str or None')
        self._title = value

    def _fix(self, atoms, fix_missing_residues=True):

        try:
            from pdbfixer import PDBFixer
            from simtk.openmm.app import PDBFile
        except ImportError:
            raise ImportError('Please install PDBFixer and OpenMM in order to use ClustENM.')

        stream = createStringIO()
        title = atoms.getTitle()
        writePDBStream(stream, atoms)
        stream.seek(0)
        fixed = PDBFixer(pdbfile=stream)
        stream.close()

        if fix_missing_residues:
            fixed.findMissingResidues()
        else:
            # skipping modeling of missing residues neither at the chain ends
            # nor in the middle, since sometimes it models unnaturally
            fixed.missingResidues = {}

        fixed.findNonstandardResidues()
        fixed.replaceNonstandardResidues()
        fixed.removeHeterogens(False)
        fixed.findMissingAtoms()
        fixed.addMissingAtoms()
        fixed.addMissingHydrogens(self._ph)

        stream = createStringIO()
        PDBFile.writeFile(fixed.topology, fixed.positions, stream, keepIds=True)
        stream.seek(0)
        self._atoms = parsePDBStream(stream)
        self._atoms.setTitle(title)
        stream.close()

        self._topology = fixed.topology
        self._positions = fixed.positions

    def _min_sim(self, arg):

        # arg: coordset   (numAtoms, 3) in Angstrom, which should be converted into nanometer

        # we are not using self._positions!
        # arg will be set as positions

        try:
            from simtk.openmm import Platform, LangevinIntegrator
            from simtk.openmm.app import Modeller, ForceField, CutoffNonPeriodic, Simulation, \
                                         HBonds, StateDataReporter
            from simtk.unit import nanometers, picosecond, kelvin, angstrom, kilojoule_per_mole, \
                                   MOLAR_GAS_CONSTANT_R
        except ImportError:
            raise ImportError('Please install PDBFixer and OpenMM in order to use ClustENM.')

        modeller = Modeller(self._topology,
                            self._positions)
        forcefield = ForceField('amber99sbildn.xml',
                                'amber99_obc.xml')

        system = forcefield.createSystem(modeller.topology,
                                         nonbondedMethod=CutoffNonPeriodic,
                                         nonbondedCutoff=1.0*nanometers,
                                         constraints=HBonds)

        integrator = LangevinIntegrator(self._temp*kelvin,
                                        1/picosecond,
                                        0.002*picosecond)

        # precision could be mixed, but single is okay.
        platform = self._platform if self._platform is None else Platform.getPlatformByName(self._platform)
        properties = None

        if self._platform is None:
            properties = {'Precision': 'single'}
        elif self._platform in ['CUDA', 'OpenCL']:
            properties = {'Precision': 'single'}
            
        simulation = Simulation(modeller.topology, system, integrator,
                                platform, properties)

        # automatic conversion into nanometer will be carried out.

        simulation.context.setPositions(arg * angstrom)

        try:
            if self._sim:
                simulation.minimizeEnergy()
                # heating-up the system incrementally
                sdr = StateDataReporter(stdout, 1, step=True, temperature=True)
                sdr._initializeConstants(simulation)
                temp = 0.0

                # instantaneous temperature could be obtained by openmmtools module
                # but its installation using conda may have problem due to repository freezing,
                # therefore, we are evaluating by hand.

                while temp < self._temp:
                    simulation.step(1)
                    ke = simulation.context.getState(getEnergy=True).getKineticEnergy()
                    temp = (2 * ke / (sdr._dof * MOLAR_GAS_CONSTANT_R)).value_in_unit(kelvin)

                simulation.step(self._t_steps[self._cycle])
            else:
                simulation.minimizeEnergy()
            pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(angstrom)
            pot = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilojoule_per_mole)

            return pot, pos

        except BaseException as be:
            LOGGER.warning('OpenMM exception: ' + be.__str__() + ' so the corresponding conformer will be discarded!')

            return np.nan, np.full_like(arg, np.nan)

    def _sample_v1(self, conf):

        # arg: conf idx

        tmp = self._atoms.copy()
        tmp.setCoords(conf)
        ca = tmp.ca

        anm_ca = ANM()
        anm_ca.buildHessian(ca, cutoff=self._cutoff)

        # 1e-6 is the same value of prody's ZERO parameter
        rank_diff = (3 * self._n_ca - 6
                     - np.linalg.matrix_rank(anm_ca.getHessian(),
                                             tol=ZERO, hermitian=True))
        if rank_diff != 0:
            # taking care cases with more than 6 zeros
            # maybe an exception can be raised in debug mode
            LOGGER.warn('Abnormal number of zeros modes detected (%d detected, 6 expected).'%(6 + rank_diff))
            return None

        anm_ca.calcModes(self._n_modes)

        anm_ex, _ = extendModel(anm_ca, ca, tmp, norm=True)
        a = np.array(list(product([-1, 0, 1], repeat=self._n_modes)))

        nv = (anm_ex.getEigvecs() / np.sqrt(anm_ex.getEigvals())) @ a.T

        nvn = nv / np.linalg.norm(nv, axis=0).max()

        d = (self._rmsd[self._cycle] * np.sqrt(tmp.numAtoms()) * nvn).T
        d = d.reshape(d.shape[0], -1, 3)

        r0 = tmp.getCoords()
        r = r0 + d

        return r

    def _sample(self, conf):
        tmp = self._atoms.copy()
        tmp.setCoords(conf)
        ca = tmp.ca

        anm_ca = ANM()
        anm_ca.buildHessian(ca, cutoff=self._cutoff)

        # use prody's ZERO parameter
        rank_diff = (3 * self._n_ca - 6
                     - np.linalg.matrix_rank(anm_ca.getHessian(),
                                             tol=ZERO, hermitian=True))
        if rank_diff != 0:
            # taking care cases with more than 6 zeros
            # maybe an exception can be raised in debug mode
            LOGGER.warn('Abnormal number of zeros modes detected (%d detected, 6 expected).'%(6 + rank_diff))
            return None

        anm_ca.calcModes(self._n_modes)

        anm_ex, _ = extendModel(anm_ca, ca, tmp, norm=True)
        ens_ex = sampleModes(anm_ex, atoms=tmp,
                                n_confs=self._n_confs,
                                rmsd=self._rmsd[self._cycle])

        return ens_ex.getCoordsets()

    def _rmsds(self, arg):

        # as long as there is no need for superposing conformations
        # only anm modes are used for perturbation
        # so no translation or rotation would involve

        # arg: coords (n_conf, n_ca, 3)

        tmp = arg.reshape(-1, 3 * self._n_ca)

        return pdist(tmp) / np.sqrt(self._n_ca)

    def _hc(self, arg):

        # arg: coords   (n_conf, n_ca, 3)

        rmsds = self._rmsds(arg)
        # optimal_ordering=True can be slow, particularly on large datasets.
        link = linkage(rmsds, method='average')

        # fcluster gives cluster labels starting from 1

        if self._threshold is not None:
            hcl = fcluster(link, t=self._threshold[self._cycle],
                           criterion='distance') - 1

        if self._maxclust is not None:
            hcl = fcluster(link, t=self._maxclust[self._cycle],
                           criterion='maxclust') - 1

        return hcl

    def _centroid(self, arg):

        # arg: coords   (n_conf_clust, n_ca, 3)

        if arg.shape[0] > 2:
            rmsds = self._rmsds(arg)
            sim = np.exp(- squareform(rmsds) / rmsds.std())
            idx = sim.sum(1).argmax()
            return idx
        else:
            return 0   # or np.random.randint(low=0, high=arg.shape[0])

    def _centers(self, *args):

        # args[0]: coords   (n_conf, n_ca, 3)
        # args[1]: labels

        nl = np.unique(args[1])
        idx = OrderedDict()
        for i in nl:
            idx[i] = np.where(args[1] == i)[0]

        # Dictionary order is guaranteed to be insertion order by Python 3.7!
        wei = [idx[k].size for k in idx.keys()]
        centers = np.empty(nl.size, dtype=int)
        for i in nl:
            tmp = self._centroid(args[0][idx[i]])
            centers[i] = idx[i][tmp]

        return centers, wei

    def _generate(self, confs):

        # arg: previous generation no

        LOGGER.info('Sampling conformers in generation %d ...'%self._cycle)
        LOGGER.timeit('_clustenm_gen')
        tmp = []
        for conf in confs:
            if not self._v1:
                ret = self._sample(conf)
                if ret is not None:
                    tmp.append(ret)
                else:
                    # we may raise an exception in debug mode
                    LOGGER.info('more than 6 zero eigenvalues!')
            else:
                ret = self._sample_v1(conf)
                if ret is not None:
                    tmp.append(ret)
                else:
                    LOGGER.info('more than 6 zero eigenvalues!')

        confs_ex = np.concatenate(tmp)

        confs_ca = confs_ex[:, self._idx_ca]

        LOGGER.info('Clustering in generation %d ...'%self._cycle)
        label_ca = self._hc(confs_ca)
        centers, wei = self._centers(confs_ca, label_ca)
        LOGGER.report('Centroids were generated in %.2fs.', label='_clustenm_gen')

        return confs_ex[centers], wei

    def _outliers(self, arg):

        # arg : potential energies
        # outliers are detected by modified z_score.

        tmp = 0.6745 * (arg - np.median(arg)) / median_absolute_deviation(arg)
        # here the assumption is that there is not just one conformer

        return tmp > 3.5

    def _superpose_ca(self, confs):
        tmp0 = self._getCoords()
        n = confs.shape[0]
        tmp1 = []
        for i in range(n):
            tmp2 = calcTransformation(confs[i, self._idx_ca],
                                         tmp0[self._idx_ca])
            tmp1.append(applyTransformation(tmp2, confs[i]))

        return np.array(tmp1)

    def _build(self, conformers, keys, potentials, counts):

        #coords = conformers[0][0]
        #self.setCoords(coords)
        self.addCoordset(conformers)
        self.setData('count', counts)
        self.setData('key', keys)
        self.setData('potential', potentials)

    def addCoordset(self, coords):
        self._indexer = None
        super(ClustENM, self).addCoordset(coords)

    def getData(self, key, gen=None):
        keys = super(ClustENM, self)._getData('key')
        data = super(ClustENM, self).getData(key)

        if gen is not None:
            data_ = []
            for k, d in zip(keys, data):
                g, _ = k
                if g == gen:
                    data_.append(d)
            data = np.array(data_)
        return data

    def getKeys(self, gen=None):
        return self.getData('key', gen)

    def getLabels(self, gen=None):
        keys = self.getKeys(gen)
        labels = ['%d_%d'%tuple(k) for k in keys]
        return labels

    def getPotentials(self, gen=None):
        return self.getData('potential', gen)

    def getSizes(self, gen=None):
        return self.getData('size', gen)

    def numGenerations(self):
        return self._n_gens

    def numConfs(self, gen=None):
        if gen is None:
            return super(ClustENM, self).numConfs()
        
        keys = self._getData('key')
        n_confs = 0
        for g, _ in keys:
            if g == gen:
                n_confs += 1
        return n_confs

    def _slice(self, indices):
        if len(indices) == 0:
            raise ValueError('indices (tuple) cannot be empty')
        
        if self._indexer is None:
            keys = self._getData('key')
            entries = [[] for _ in range(self.numGenerations() + 1)]
            for i, (gen, _) in enumerate(keys):
                entries[gen].append(i)
            
            n_conf_per_gen = np.max([len(entry) for entry in entries])
            for entry in entries:
                for i in range(len(entry), n_conf_per_gen):
                    entry.append(-1)

            indexer = self._indexer = np.array(entries)
        else:
            indexer = self._indexer
        full_serials = indexer[indices].flatten()
        
        serials = []
        for s in full_serials:
            if s != -1:
                serials.append(s)

        return np.array(serials)

    def _getCoordsets(self, indices=None, selected=True):
        """Returns the coordinate set(s) at given *indices*, which may be
        an integer, a list of integers, a tuple of (generation, index), or **None**. 
        **None** returns all coordinate sets. For reference coordinates, use :meth:`getCoords`
        method."""

        if isinstance(indices, tuple):
            I = self._slice(indices)
            if len(I) == 0:
                raise IndexError('index out of range (%d, %d)'%indices)
            return super(ClustENM, self)._getCoordsets(I, selected)
        else:
            return super(ClustENM, self)._getCoordsets(indices, selected)

    def writePDB(self, folder='.', single=True, **kwargs):

        # single -> True, save as a single pdb file with each conformer as a model
        # otherwise, each conformer is saved as a separate pdb file
        # in the directory pdbs_pdbname

        title = self.getTitle()

        LOGGER.timeit('t0')
        if single:
            LOGGER.info('Saving %s.pdb ...'%title)
            writePDB(folder + '/' + title, self)
        else:
            direc = folder + '/' + title
            if isdir(direc):
                LOGGER.warn('%s is not empty; will be flooded'%direc)
            else:
                mkdir(direc)

            LOGGER.info('Saving in %s ...'%direc)
            for i, lab in enumerate(self.getLabels()):
                filename = '%s/%s'%(direc, lab)
                writePDB(filename, self, csets=i)
        LOGGER.report(label='t0')

    def run(self, cutoff=15., n_modes=3, n_confs=50, rmsd=1.0,
            n_gens=5, maxclust=None, threshold=None,
            sim=True, temp=300, t_steps_i=1000, t_steps_g=7500,
            outlier=True, mzscore=3.5, **kwargs):

        if self.isBuilt():
            raise ValueError('ClustENM ensemble has been built; please start a new instance')

        # set up parameters
        self._cutoff = cutoff
        self._n_modes = n_modes
        self._n_confs = n_confs
        self._rmsd = (0.,) + rmsd if isinstance(rmsd, tuple) else (0.,) + (rmsd,) * n_gens
        self._n_gens = n_gens
        self._v1 = kwargs.pop('v1', False) 
        self._platform = kwargs.pop('platform', None) 

        if maxclust is None:
            self._maxclust = None
        else:
            if isinstance(maxclust, tuple):
                self._maxclust = (0,) + maxclust
            else:
                self._maxclust = (0,) + (maxclust,) * n_gens

        if threshold is None:
            self._threshold = None
        else:
            if isinstance(threshold, tuple):
                self._threshold = (0,) + threshold
            else:
                self._threshold = (0,) + (threshold,) * n_gens

        self._sim = sim
        self._temp = temp

        if self._sim:
            if isinstance(t_steps_g, tuple):
                self._t_steps = (t_steps_i,) + t_steps_g
            else:
                self._t_steps = (t_steps_i,) + (t_steps_g,) * n_gens

        self._outlier = outlier
        self._mzscore = mzscore

        self._cycle = 0

        # t0 = perf_counter()
        # prody.logger.timeit doesn't give the correct overal time,
        # that's why, perf_counter is being used!
        # but prody.logger.timeit is still here to see
        # how much it differs
        LOGGER.timeit('_clustenm_overall')

        LOGGER.info('Generation 0 ...')

        if self._sim:
            LOGGER.info('Minimization, heating-up & simulation in generation 0 ...')
        else:
            LOGGER.info('Minimization in generation 0 ...')
        LOGGER.timeit('_clustenm_min')
        potential, conformer = self._min_sim(self._atoms.getCoords())
        if np.isnan(potential):
            raise ValueError('Initial structure could not be minimized. Try again and/or check your structure.')

        LOGGER.report('Structure was energetically minizied in %.2fs.', label='_clustenm_min')
        LOGGER.info('#' + '-' * 19 + '/*\\' + '-' * 19 + '#')
        # self._potentials[0] = np.array([potential])
        # self._conformers[0] = conformer.reshape(1, *conformer.shape)
        self.setCoords(conformer)

        potentials = [potential]
        sizes = [1]
        conf = conformer.reshape((1, *conformer.shape))
        conformers = start_confs = conf
        keys = [(0, 0)]

        for i in range(1, self._n_gens+1):
            self._cycle += 1
            LOGGER.info('Generation %d ...'%i)
            confs, weights = self._generate(start_confs)
            if self._sim:
                LOGGER.info('Minimization, heating-up & simulation in generation %d ...'%i)
            else:
                LOGGER.info('Minimization in generation %d ...'%i)
            LOGGER.timeit('_clustenm_min_sim')

            # <simtk.openmm.openmm.Platform;
            # proxy of a <Swig Object of type 'OpenMM::Platform'>>
            # which isn't picklable because it is a SwigPyObject
            # that's why, the loop is serial when a platform is specified.
            # ThreadPool can overcome this issue, however,
            # threading can slow calculations down
            # since calculations here are CPU-bound.
            # we may totally discard Pool and just use serial version!

            if self._platform is None:
                with Pool(cpu_count()) as p:
                    pot_conf = p.map(self._min_sim, confs)
            else:
                pot_conf = [self._min_sim(conf) for conf in confs]

            LOGGER.report('Structures were sampled in %.2fs.', label='_clustenm_min_sim')
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
            start_confs = self._superpose_ca(confs[idx])
            #conformers[i] = self._superpose_ca(conformers[idx])
            for j in range(start_confs.shape[0]):
                keys.append((i, j))
            conformers = np.vstack((conformers, start_confs))

        LOGGER.timeit('_clustenm_ens')
        LOGGER.info('Creating an ensemble of conformers ...')

        self._build(conformers, keys, potentials, sizes)
        LOGGER.report('Ensemble was created in %.2fs.', label='_clustenm_ens')

        LOGGER.report('All completed in %.2fs.', label='_clustenm_overall')
        # t1 = perf_counter()
        # t10 = round(t1 - t0, 2)
        # LOGGER.info('All completed in %.2fs.'%t10)

    def writeParameters(self, filename=None):

        title = self.getTitle()
        if filename is None:
            filename = '%s_parameters.txt'%title

        with open(filename, 'w') as f:
            f.write(f'title = {title}\n')
            f.write(f'pH = {self._ph}\n')
            f.write(f'cutoff = {self._cutoff}\n')
            f.write(f'n_modes = {self._n_modes}\n')
            if not self._v1:
                f.write(f'n_confs = {self._n_confs}\n')
            f.write(f'rmsd = {self._rmsd[1:]}\n')
            f.write(f'n_gens = {self._n_gens}\n')
            if self._threshold is not None:
                f.write(f'threshold = {self._threshold[1:]}\n')
            if self._maxclust is not None:
                f.write(f'maxclust = {self._maxclust[1:]}\n')
            if self._sim:
                f.write(f'temp = {self._temp}\n')
                f.write(f't_steps = {self._t_steps}\n')
            if self._outlier:
                f.write(f'outlier = {self._outlier}\n')
                f.write(f'mzscore = {self._mzscore}\n')
            if self._v1:
                f.write(f'v1 = {self._v1}\n')
            if self._platform is not None:
                f.write(f'platform = %s\n'%self._platform)
            else:
                f.write(f'platform = Default\n')
