from os import chdir, listdir, mkdir, system
from os.path import isdir
from pickle import dump
from re import findall
from numpy import argsort, array, c_, count_nonzero, hstack, mean, median, quantile, save
from scipy.stats import zscore, median_absolute_deviation
import matplotlib.pyplot as plt
from prody import LOGGER
from prody.atomic.functions import extendAtomicData
from .anm import ANM
from .gnm import GNM
from prody.proteins import parsePDB, writePDB
from .editing import reduceModel
from .plotting import showAtomicLines
from .signature import ModeEnsemble, saveModeEnsemble
from prody.utilities import which


__all__ = ['ESSA']


class ESSA:
    '''docstring'''

    def __init__(self, pdb, lig=None):

        self._atoms = pdb
        self._title = pdb.getTitle()
        self._lig = lig
        self._heavy = pdb.select('protein and heavy and not hetatm')
        self._ca = self._heavy.ca
        self._n_modes = None
        self._enm = None
        self._cutoff = None
        self._lig = lig
        self._dist = None
        self._ensemble = None
        self._labels = None
        self._zscore = None
        if lig:
            self._lig_idx = {}

    def scanResidues(self, n_modes=10, enm='gnm', cutoff=None, dist=4.5):

        # pdb : a pdb parsed by ProDy
        # n_modes : n_modes of GNM object
        # lig : string of ligands' chainIDs and resSeqs (resnum) separated by a whitespace, e.g., 'A 300 B 301'
        # dist : the protein residues within a distance of ligands
        # enm: the type of enm
        # cutoff: enm cutoff

        self._n_modes = n_modes
        self._enm = enm
        if self._lig is not None:
            self._dist = dist

        self._ensemble = ModeEnsemble(f'{self._title}')
        self._ensemble.setAtoms(self._ca)
        self._labels = ['ref']

        # --- reference model --- #

        if self._enm == 'gnm':
            ca_enm = GNM('ca')
            if cutoff is not None:
                self._cutoff = cutoff
                ca_enm.buildKirchhoff(self._ca, cutoff=self._cutoff)
            else:
                ca_enm.buildKirchhoff(self._ca)
                self._cutoff = ca_enm.getCutoff()

        if self._enm == 'anm':
            ca_enm = ANM('ca')
            if cutoff is not None:
                self._cutoff = cutoff
                ca_enm.buildHessian(self._ca, cutoff=self._cutoff)
            else:
                ca_enm.buildHessian(self._ca)
                self._cutoff = ca_enm.getCutoff()

        ca_enm.calcModes(n_modes=n_modes)
        self._ensemble.addModeSet(ca_enm[:])

        # --- perturbed models --- #

        LOGGER.progress(msg='', steps=(self._ca.numAtoms()))
        for i in self._ca.getResindices():
            LOGGER.update(step=i+1, msg=f'scanning residue {i+1}')
            sel = f'calpha or resindex {i}'
            tmp = self._heavy.select(sel)

            if self._enm == 'gnm':
                tmp_enm = GNM(f'res_{i}')
                tmp_enm.buildKirchhoff(tmp, cutoff=self._cutoff)

            if self._enm == 'anm':
                tmp_enm = ANM(f'res_{i}')
                tmp_enm.buildHessian(tmp, cutoff=self._cutoff)

            tmp_enm_red, _ = reduceModel(tmp_enm, tmp, self._ca)
            tmp_enm_red.calcModes(n_modes=self._n_modes)

            self._ensemble.addModeSet(tmp_enm_red[:])
            self._labels.append(tmp_enm.getTitle())

        self._ensemble.setLabels(self._labels)
        self._ensemble.match()

        # --- ESSA computation part --- #

        denom = self._ensemble[0].getEigvals()
        num = self._ensemble[1:].getEigvals() - denom

        eig_diff = num / denom * 100
        eig_diff_mean = mean(eig_diff, axis=1)

        self._zscore = zscore(eig_diff_mean)

        # --- residue indices of protein residues that are within dist (4.5 A) of ligands --- #

        if self._lig:
            ligs = self._lig.split()
            ligs = list(zip(ligs[::2], ligs[1::2]))
            for chid, resnum in ligs:
                key = ''.join(chid + str(resnum))
                sel_lig = 'calpha and not hetatm and (same residue as ' \
                          f'exwithin {self._dist} of (chain {chid} and resnum {resnum}))'
                self._lig_idx[key] = self._atoms.select(sel_lig).getResindices()

    def getESSAZscores(self):

        return self._zscore

    def getESSAEnsemble(self):

        return self._ensemble[:]
    
    def saveESSAEnsemble(self):

        saveModeEnsemble(self._ensemble, filename=f'{self._title}_{self._enm}')

    def saveESSAZscores(self):

        save(f'{self._title}_{self._enm}_zs', self._zscore)

    def writeESSAZscoresToPDB(self):

        writePDB(f'{self._title}_{self._enm}_zs', self._heavy,
                 beta=extendAtomicData(self._zscore, self._ca, self._heavy)[0])

    def saveLigandIndices(self):

        if self._lig:
            dump(self._lig_idx, open(f'{self._title}_ligand_resindices.pkl', 'wb'))
        else:
            LOGGER.warning('No ligand provided.')

    def showESSAProfile(self, quant=.75):

        showAtomicLines(self._zscore, atoms=self._ca, c='k', linewidth=1.)

        if self._lig:
            zs_lig = {k: self._zscore[v] for k, v in self._lig_idx.items()}
            for k in self._lig_idx.keys():
                plt.scatter(self._lig_idx[k], zs_lig[k], label=k)
            plt.legend()
        plt.hlines(quantile(self._zscore, q=quant),
                   xmin=0., xmax=self._ca.numAtoms(),
                   linestyle='--', color='c')

        plt.xlabel('Residue')
        plt.ylabel('Z-Score')

        plt.tight_layout()

    def scanPockets(self):

        fpocket = which('fpocket')

        if fpocket is None:
            LOGGER.info('Fpocket 3.0 was not found, please install it.')
            return None

        try:
            from pandas import Index, DataFrame
        except ImportError as ie:
            LOGGER.info(ie.__str__() + ' was found, please install it.')
            return None

        rcr = {(i, j): k for i, j, k in zip(self._ca.getChids(),
                                            self._ca.getResnums(),
                                            self._ca.getResindices())}

        writePDB(f'{self._title}_pro', self._heavy)

        direc = f'{self._title}_pro_out'
        if not isdir(direc):
            system(f'fpocket -f {self._title}_pro.pdb')

        chdir(direc + '/pockets')
        l = [x for x in listdir('.') if x.endswith('.pdb')]
        l.sort(key=lambda x:int(x.partition('_')[0][6:]))

        ps = []
        for x in l:
            with open(x, 'r') as f:
                tmp0 = f.read()
                tmp1 = [float(x[1]) for x in findall(r'(.+:\s+)(-*[\d.]+)(\n)', tmp0)]
            ps.append(tmp1)
        pdbs = parsePDB(l)
        chdir('../..')

        # ----- # ----- #

        ps = array(ps)

        pcn = {int(pdb.getTitle().partition('_')[0][6:]):
               set(zip(pdb.getChids().tolist(),
                       pdb.getResnums().tolist())) for pdb in pdbs}
        pi = {p: [rcr[x] for x in crn] for p, crn in pcn.items()}

        pzs_max = {k: max(self._zscore[v]) for k, v in pi.items()}
        pzs_med = {k: median(self._zscore[v]) for k, v in pi.items()}

        # ----- # ----- #

        indices = Index(range(1, ps.shape[0] + 1), name='Pocket')

        columns = Index(['Pocket Score',
                         'Drug Score',
                         'Number of alpha spheres',
                         'Mean alpha-sphere radius',
                         'Mean alpha-sphere Solvent Acc.',
                         'Mean B-factor of pocket residues',
                         'Hydrophobicity Score',
                         'Polarity Score',
                         'Amino Acid based volume Score',
                         'Pocket volume (Monte Carlo)',
                         'Pocket volume (convex hull)',
                         'Charge Score',
                         'Local hydrophobic density Score',
                         'Number of apolar alpha sphere',
                         'Proportion of apolar alpha sphere'],
                        name='Feature')

        self._df = DataFrame(index=indices, columns=columns, data=ps)

        # ----- # ----- #

        columns_zs = Index(['Pocket Z-score',
                            'Drug Z-score',
                            'Number of alpha spheres',
                            'Mean alpha-sphere radius',
                            'Mean alpha-sphere Solvent Acc.',
                            'Mean B-factor of pocket residues',
                            'Hydrophobicity Z-score',
                            'Polarity Z-score',
                            'Amino Acid based volume Z-score',
                            'Pocket volume (Monte Carlo)',
                            'Pocket volume (convex hull)',
                            'Charge Z-score',
                            'Local hydrophobic density Z-score',
                            'Number of apolar alpha sphere',
                            'Proportion of apolar alpha sphere',
                            'Maximum ESSA Z-score of pocket residues',
                            'Median ESSA Z-score of pocket residues'],
                           name='Feature')

        zps = zscore(ps, axis=0)
        zps = hstack((zps, c_[list(pzs_max.values())]))
        zps = hstack((zps, c_[list(pzs_med.values())]))

        self._df_zs = DataFrame(index=indices, columns=columns_zs, data=zps)

    def rankPockets(self):

        from pandas import DataFrame

        pzs_max = self._df_zs['Maximum ESSA Z-score of pocket residues']
        pzs_med = self._df_zs['Median ESSA Z-score of pocket residues']
        self._idx_pzs_max = argsort(list(pzs_max))[::-1] + 1
        self._idx_pzs_med = argsort(list(pzs_med))[::-1] + 1

        lhd = self._df_zs.loc[:, 'Local hydrophobic density Z-score']
        n = count_nonzero(lhd >= 0.)
        q = quantile(lhd, q=.85)

        # ----- # ------ #

        s_max = ['Maximum ESSA Z-score of pocket residues',
                 'Local hydrophobic density Z-score']

        zf_max = self._df_zs[s_max].copy()

        if n >= lhd.size // 4:
            f_max = zf_max.iloc[:, 1] >= 0.
        else:
            f_max = zf_max.iloc[:, 1] >= q

        zf_max = zf_max[f_max]

        zf_max.iloc[:, 0] = zf_max.iloc[:, 0].round(1)
        zf_max.iloc[:, 1] = zf_max.iloc[:, 1].round(2)

        self._idx_max = zf_max.sort_values(s_max, ascending=False).index

        # ----- # ----- #

        s_med = ['Median ESSA Z-score of pocket residues',
                 'Local hydrophobic density Z-score']

        zf_med = self._df_zs[s_med].copy()

        if n >= lhd.size // 4:
            f_med = zf_med.iloc[:, 1] >= 0.
        else:
            f_med = zf_med.iloc[:, 1] >= q

        zf_med = zf_med[f_med]

        zf_med.iloc[:, 0] = zf_med.iloc[:, 0].round(1)
        zf_med.iloc[:, 1] = zf_med.iloc[:, 1].round(2)

        self._idx_med = zf_med.sort_values(s_med, ascending=False).index

        self._pocket_ranks = DataFrame(columns=['ESSA_max', 'ESSA_med',
                                                'ESSA_max_loc_hydro',
                                                'ESSA_med_loc_hydro'])

        self._pocket_ranks.iloc[:, 0] = self._idx_pzs_max
        self._pocket_ranks.iloc[:, 1] = self._idx_pzs_med
        self._pocket_ranks.iloc[:self._idx_max.size, 2] = self._idx_max
        self._pocket_ranks.iloc[:self._idx_med.size, 3] = self._idx_med

    def getPocketRanks(self):

        return self._pocket_ranks

    def getPocketFeatures(self):

        return self._df

    def getPocketZscores(self):

        return self._df_zs

    def showPocketZscores(self):

        self._df_zs[['Maximum ESSA Z-score of pocket residues',
                     'Median ESSA Z-score of pocket residues',
                     'Local hydrophobic density Z-score']].plot.bar(figsize=(25, 10))
        plt.xticks(rotation=0)
        plt.ylabel('Z-score')
        plt.tight_layout()

    def savePocketFeatures(self):

        self._df.to_pickle(f'{self._title}_pocket_features.pkl')

    def savePocketZscores(self):

        self._df_zs.to_pickle(f'{self._title}_pocket_zscores.pkl')

    def savePocketRanks(self):

        save(f'{self._title}_{self._enm}_pocket_ranks_wrt_ESSA_max',
             self._idx_pzs_max)
        save(f'{self._title}_{self._enm}_pocket_ranks_wrt_ESSA_med',
             self._idx_pzs_med)
        save(f'{self._title}_{self._enm}_pocket_ranks_wrt_ESSA_max_loc_hydro',
             self._idx_max)
        save(f'{self._title}_{self._enm}_pocket_ranks_wrt_ESSA_med_loc_hydro',
             self._idx_med)
        
    def writePocketRanksToCSV(self):

        self._pocket_ranks.to_csv(f'{self._title}_{self._enm}_pocket_ranks.csv', index=False)
