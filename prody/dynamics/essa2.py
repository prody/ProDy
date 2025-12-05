"""
ESSA-2 implementation for ProDy integration

"""
import sys
import subprocess

from prody import *
from pylab import *
from typing import *
from collections import defaultdict
from Bio.PDB import PDBIO, PDBParser, Atom
import numpy as np
from numpy import array, median, c_, hstack, savez, save
import matplotlib.pyplot as plt       
from tqdm import tqdm
import os
import contextlib
import sys
import tempfile
from scipy.ndimage import gaussian_filter
from scipy.stats import rankdata
import logging
import shutil
import subprocess
from pathlib import Path
import re

try:
    from pandas import DataFrame, Index
except Exception:
    DataFrame = Index = None

class ESSA2:

    def __init__(self):
        self._system_ready = False
        confProDy(verbosity='none')
        self._zscores = None
        self._msf_matrix = None
        self._ref_msfs = None
        self._ligres_idx = {}
        self._ligres_code = {}
        self.TRP_ROTAMERS = [                      
                {"chi1": -60, "chi2":  90},
                #{"chi1": -60, "chi2": -90},
                #{"chi1":  60, "chi2":  90},
                #{"chi1":  60, "chi2": -90},
                #{"chi1": 180, "chi2":  90},
                #{"chi1": 180, "chi2": -90},
            ]

        self.TRP_TEMPLATE = {
                    "CG":  [1.50,   0.00,   0.00],
                    "CD1": [2.10,  -1.20,   0.00],
                    "NE1": [3.40,  -1.20,   0.00],
                    "CE2": [4.00,   0.00,   0.00],
                    "CD2": [2.80,   1.20,   0.00],
                    "CE3": [2.80,   2.50,   0.00],
                    "CZ3": [4.00,   2.50,   0.00],
                    "CH2": [4.60,   1.20,   0.00],
                    "CZ2": [5.20,   0.00,   0.00],
                }


    def setSystem(self, pdb_input, lig: Optional[str] = None, cutoff: float = 10.0,
                target_variance: Optional[float] = None, num_modes: Optional[int] = None,
                lig_dist: float = 4.5):

        """
        Set protein system, compute reference GNM, and identify ligand-binding residues.

        Parameters
        ----------
        pdb_input : str or AtomGroup
            PDB file path or ProDy AtomGroup.
        lig : str, optional
            Ligand specification, e.g., 'A 300 A 301'.
        cutoff : float
            Cutoff for GNM Kirchhoff matrix (Å).
        target_variance : float, optional
            Fraction of variance (%) to determine number of modes.
        num_modes : int, optional
            Number of modes to compute; overrides default 10 if target_variance is None.
        lig_dist : float, default=4.5
            Distance (Å) to define residues near the ligand.
        """

        if isinstance(pdb_input, str):
            self._pdb_path = pdb_input
            self._atoms = parsePDB(pdb_input)

        elif hasattr(pdb_input, "getTitle"):  # AtomGroup
            self._atoms = pdb_input
            self._pdb_path = os.path.join(tempfile.gettempdir(), "essa_input.pdb")
            writePDB(self._pdb_path, self._atoms)

        else:
            raise TypeError("pdb_input must be a file path or ProDy AtomGroup.")

        self._title = self._atoms.getTitle()
        self._cutoff = cutoff
        self._lig = lig
        self._lig_dist = lig_dist

        self._heavy = self._atoms.select('protein and heavy')
        self._ca = self._heavy.ca

        gnm = GNM()
        gnm.buildKirchhoff(self._ca, cutoff=self._cutoff)

        # Determine number of modes
        if target_variance is not None:
            gnm.calcModes('all')
            cumulative = 0.0
            count = 0
            for mode in gnm:
                cumulative += calcFractVariance(mode) * 100
                count += 1
                if cumulative >= target_variance:
                    break
            self.num_modes = max(1, count)
            print(f"Calculating {self.num_modes} modes for {target_variance}% target variance.")
        else:
            self.num_modes = 10 if num_modes is None else num_modes
            print(f"Calculating {self.num_modes} modes.")

        gnm.calcModes(self.num_modes)
        self._gnm = gnm
        self._ref_msfs = calcSqFlucts(gnm)

        # Protein residues that are within dist Å of ligands
        self._lig = lig
        self._lig_dist = lig_dist
        self._ligres_idx = {}
        self._ligres_code = {}

        if self._lig:
            ligs = self._lig.split()
            ligs = list(zip(ligs[::2], ligs[1::2]))  # (chain, resnum) pairs

            for chid, resnum in ligs:
                key = ''.join(chid + str(resnum))

                sel_lig = (
                    'calpha and not hetatm and (same residue as '
                    'exwithin {} of (chain {} and resnum {}))'
                ).format(self._lig_dist, chid, resnum)

                
                self._ligres_idx[key] = self._atoms.select(sel_lig).getResindices()

            for k, v in self._ligres_idx.items():
                atoms = self._ca.select('resindex ' + ' '.join([str(i) for i in v]))
                tmp0 = defaultdict(list)

                for ch, rn in zip(atoms.getChids(), atoms.getResnums()):
                    tmp0[ch].append(str(rn))

                tmp1 = {ch: ' '.join(rn) for ch, rn in tmp0.items()}

                self._ligres_code[k] = [
                    'chain {} and resnum {}'.format(ch, rn)
                    for ch, rn in tmp1.items()
                ]
        self._system_ready = True


    def refMSFs(self) -> np.ndarray:
          
        if not self._system_ready:
            raise RuntimeError("System not set. Call setSystem() first.")        
        #print(f"Computing reference MSFs using {self.num_modes} modes.")

        return self._ref_msfs

    def _write_trp_mutation_pdb(self, res_idx, chi, tmp_pdb_path):
        """
        Write a PDB file with a single residue mutated to TRP and 
        rotated according to chi angles.

        Handles missing CB atoms (e.g., Glycine) by estimating CB position.
        """

        if self._pdb_path is None:
            raise ValueError("PDB path is not set. Call setSystem first.")

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("WT", self._pdb_path)
        model = next(structure.get_models())  # assume first model

        # Identify residue by index (assuming 0-based indexing)
        res = next((r for r in model.get_residues() if r.get_id()[1] == res_idx), None)
        if res is None:
            raise ValueError(f"Residue index {res_idx} not found in PDB.")

        CA = res['CA'].get_coord()

        # Estimate CB if missing (for glycine or incomplete PDB)
        try:
            CB = res['CB'].get_coord()
        except KeyError:
            N = res['N'].get_coord()
            C = res['C'].get_coord()
            b = CA - N
            c = CA - C
            n = np.cross(b, c)
            n = n / np.linalg.norm(n)
            CB = CA + 1.53 * n  # typical CA-CB bond length

        # Prepare new TRP atoms
        trp_atoms = []
        for name, xyz_local in self.TRP_TEMPLATE.items():
            local = np.array(xyz_local)
            global_xyz = CA + local 
            # Create a new Bio.PDB Atom
            atom = Atom.Atom(
                name=name,
                coord=global_xyz,
                bfactor=0.0,
                occupancy=1.0,
                altloc=' ',
                fullname=name.ljust(4),
                serial_number=0,
                element=name[0] 
            )
            trp_atoms.append(atom)

        for atom in list(res):
            if atom.get_name() != 'CA':
                res.detach_child(atom.get_id())

        for atom in trp_atoms:
            res.add(atom)

        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_pdb_path)

    def scanResidues(self, residue_indices=None, tmp_dir=None):
        if not self._system_ready:
            raise RuntimeError("Call setSystem() first.")

        print(f"Running mutation scan with {self.num_modes} GNM modes.")

        if residue_indices is None:
            residue_indices = self._ca.getResnums()

        if tmp_dir is None:
            tmp_dir = tempfile.gettempdir()

        msf_matrix = np.zeros((len(residue_indices), len(self._ca)))

        print(f"Scanning {len(residue_indices)} residues.")

        for i, res_idx in enumerate(tqdm(residue_indices, desc="Residues", unit="res")):
            rotamer_flucts = []

            for rot in self.TRP_ROTAMERS:
                tmp_pdb = os.path.join(tmp_dir, f"tmp_res{res_idx}.pdb")
                self._write_trp_mutation_pdb(res_idx, rot, tmp_pdb)

                ag = parsePDB(tmp_pdb)
                ca = ag.select('calpha')
                heavy = ag.select('protein and heavy')
                cwd = heavy.select(f'calpha or resnum {res_idx}')

                gnm = GNM()
                gnm.buildKirchhoff(cwd, cutoff=self._cutoff)
                gnm.calcModes(self.num_modes)

                gnm_red, _ = reduceModel(gnm, cwd, ca)
                gnm_red.calcModes(self.num_modes)

                _, matched = matchModes(self._gnm, gnm_red)
                msf = calcSqFlucts(matched[:self.num_modes])

                rotamer_flucts.append(msf)

                try: os.remove(tmp_pdb)
                except: pass

            msf_matrix[i,:] = np.mean(rotamer_flucts, axis=0)
           
        self._msf_matrix = msf_matrix
        print(f"Stored MSF matrix of shape: {msf_matrix.shape}")
        return msf_matrix

    def zscores(self, ref_msf: np.ndarray = None, msf_perturbed: np.ndarray = None, rank_norm: bool = False) -> np.ndarray:
        """
        Compute Z-scores from MSF differences.

        Parameters
        ----------
        ref_msf : np.ndarray, optional
            Reference MSF vector. If None, uses self._ref_msfs.
        msf_perturbed : np.ndarray, optional
            Perturbed MSF matrix. If None, uses self._msf_matrix.
        rank_norm : bool
            If True, compute Z-scores based on rank normalization.

        Returns
        -------
        np.ndarray
            Computed Z-scores.
        """
        
        if ref_msf is None:
            if hasattr(self, "_ref_msfs"):
                ref_msf = self._ref_msfs
            else:
                raise ValueError("Reference MSF must be computed first.")

        if msf_perturbed is None:
            if hasattr(self, "_msf_matrix"):
                msf_perturbed = self._msf_matrix
            else:
                raise ValueError("Perturbed MSF matrix must be computed first.")

        ref_msf = ref_msf[:msf_perturbed.shape[1]]
        msf_diff = ref_msf - msf_perturbed
        avg_diff = np.mean(msf_diff, axis=1)

        if rank_norm:
            # Rank-based normalization
            ranks = rankdata(avg_diff, method="average")
            mean = ranks.mean()
            std = ranks.std()
            zscores = (ranks - mean) / std if std > 0 else ranks
            self._zscores = zscores
        else:
            mean = avg_diff.mean()
            std = avg_diff.std()
            zscores = (avg_diff - mean) / std if std > 0 else avg_diff
            self._zscores = zscores

        return zscores

    def saveMSFs(self, fn: Optional[str] = None):
        if self._ref_msfs is None:
            raise RuntimeError("Reference MSFs not computed. Run setSystem() first.")
        if fn is None:
            fn = f"{self._title}_msfs.npy"
        np.save(fn, self._ref_msfs)
        print(f"Saved reference MSFs to {fn}")

    def saveMSFmatrix(self, fn: Optional[str] = None):
        if self._msf_matrix is None:
            raise RuntimeError("MSF matrix not available. Run msfscan() first.")
        if fn is None:
            fn = f"{self._title}_msf_matrix.npy"
        np.save(fn, self._msf_matrix)
        print(f"Saved MSF matrix to {fn}")

    def saveESSAZscores(self, fn: Optional[str] = None):
        if self._zscores is None:
            raise RuntimeError("ESSA Z-scores not computed. Run zscores() first.")
        if fn is None:
            fn = f"{self._title}_zs.npy"
        np.save(fn, self._zscores)
        print(f"Saved ESSA Z-scores to {fn}")

    def loadESSAZscores(self, fn: Optional[str] = None) -> np.ndarray:
        if fn is None:
            fn = f"{self._title}_zs.npy"
        if not os.path.exists(fn):
            raise FileNotFoundError(f"Z-scores file not found: {fn}")
        self._zscores = np.load(fn)
        return self._zscores

    def saveLigResidues(self, fn="lig_resindices.npy"):
        all_res = sorted({idx for group in self._ligres_idx.values() for idx in group})
        np.save(fn, np.array(all_res, dtype=int))
        print(f"Saved ligand-binding residue indices to {fn}")

    def saveLigResidueZscores(self, fn="lig_zscores.npy"):
        if self._zscores is None:
            raise RuntimeError("ESSA Z-scores not computed. Run zscores() first.")
        all_res = sorted({idx for group in self._ligres_idx.values() for idx in group})
        np.save(fn, self._zscores[all_res])
        print(f"Saved ligand-binding Z-scores to {fn}")

    def plotESSAZscores(self, zscores: Optional[np.ndarray] = None,
                        title: str = "ESSA Profile",
                        figsize: Tuple[int,int] = (6, 5)):
        if zscores is None:
            if self._zscores is None:
                raise RuntimeError("No ESSA Z-scores available. Compute or load them first.")
            zscores = self._zscores

        p75 = np.percentile(zscores, 75)
        plt.figure(figsize=figsize)
        plt.plot(zscores, lw=1.8, color='gray', label="ESSA Z-scores")
        plt.axhline(p75, linestyle='--', color='black', alpha=0.6,
                    label=f"75th percentile = {p75:.2f}")

        if self._ligres_idx:
            lig_indices = sorted({idx for group in self._ligres_idx.values() for idx in group})
            plt.scatter(lig_indices, zscores[lig_indices], color='red', s=45,
                        label="Ligand-binding residues")

        plt.xlabel("Residue index")
        plt.ylabel("ESSA Z-score")
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        plt.show()

    def writeESSAZscoresToPDB(self, zscores: Optional[np.ndarray] = None, fn: Optional[str] = None):

        if zscores is None:
            if self._zscores is None:
                raise RuntimeError("No ESSA Z-scores available. Compute or load them first.")
            zscores = self._zscores

        if not hasattr(self, "_heavy") or not hasattr(self, "_ca"):
            raise RuntimeError("System not set. Run setSystem() first.")
      
        if fn is None:
            fn = f"{self._title}_zs.pdb"

        from prody.atomic import extendAtomicData
        beta_column, _ = extendAtomicData(zscores, self._ca, self._heavy)       
        writePDB(fn, self._heavy, beta=beta_column)
        print(f"Saved ESSA Z-scores to PDB file: {fn}")

    def plotMSFchangeMap(self, msf_vals: Optional[np.ndarray] = None, ref_msf: Optional[np.ndarray] = None,
                       highlight_indices: Optional[Sequence[int]] = None,
                       sigma: float = 0.8, figsize: tuple = (5, 4), cmap: str = 'bwr_r',
                        percentile_clip: float = 98):

        if msf_vals is None:
            if not hasattr(self, "_msf_matrix"):
                raise RuntimeError("MSF matrix not computed. Run msfscan() first.")
            msf_vals = self._msf_matrix

        if ref_msf is None:
            if not hasattr(self, "_ref_msfs"):
                raise RuntimeError("Reference MSF not computed. Run setSystem() first.")
            ref_msf = self._ref_msfs

        ref_msf = ref_msf[:msf_vals.shape[1]]

        msf_diff = ref_msf - msf_vals
        smoothed = gaussian_filter(msf_diff, sigma=sigma)
        avg_diff = np.mean(msf_diff, axis=1)
        residue_indices = np.arange(len(avg_diff))

        if highlight_indices is None:
            if hasattr(self, "_ligres_idx") and self._ligres_idx:
                highlight_indices = sorted({idx for group in self._ligres_idx.values() for idx in group})
            else:
                highlight_indices = []

        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_axes([0.06, 0.1, 0.64, 0.8])
        v = np.percentile(np.abs(smoothed), percentile_clip)
        im = ax1.imshow(smoothed, aspect='auto', origin='lower', cmap=cmap, vmin=-v, vmax=v)
        ax1.set_xlabel("Residues (j)", fontsize=16)
        ax1.set_ylabel("Residues (i)", fontsize=16)
        ax1.tick_params(labelsize=12)

        ax2 = fig.add_axes([0.71, 0.1, 0.10, 0.8], sharey=ax1)
        ax2.plot(avg_diff, residue_indices, color='gray', lw=1.5, zorder=1)
        if len(highlight_indices) > 0:
            ax2.scatter(avg_diff[highlight_indices], highlight_indices,
                        color='red', edgecolors='black', s=20)
        ax2.set_xticks([])
        ax2.tick_params(left=False, labelleft=False)
        for spine in ax2.spines.values():
            spine.set_visible(False)

        cax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label(r"$\Delta \mathrm{MSF}_{j|i}$ ", fontsize=16)
        cbar.ax.tick_params(labelsize=12)

        plt.show()
        plt.close(fig)

    def scanPockets(self):
        """
        Generates ESSA z-scores for pockets and parses pocket features.
        Requires Fpocket 3.0 and Pandas.
        """
        import os
        import numpy as np
        from re import findall
        from pandas import Index, DataFrame
        from shutil import which
        from scipy.stats import zscore

        fpocket = which('fpocket')
        if fpocket is None:
            print("Fpocket (version >= 3.0) not found. Please install it.")
            return None

        rcr = {(ch, rn): idx for ch, rn, idx in zip(
            self._ca.getChids(),
            self._ca.getResnums(),
            self._ca.getResindices()
        )}

        pdb_name = f"{self._title}_pro.pdb"
        writePDB(pdb_name, self._heavy)

        out_dir = f"{self._title}_pro_out"
        if not os.path.isdir(out_dir):
            os.system(f"fpocket -f {pdb_name}")

        pockets_dir = os.path.join(out_dir, "pockets")
        if not os.path.isdir(pockets_dir):
            print("Fpocket output directory not found.")
            return None

        cwd0 = os.getcwd()
        os.chdir(pockets_dir)

        pocket_files = [x for x in os.listdir('.') if x.endswith('.pdb')]
        pocket_files.sort(key=lambda x: int(x.partition('_')[0][6:]))

        ps = []
        fea = None
        pdbs = []

        for pdb_file in pocket_files:
            with open(pdb_file, 'r') as f:
                content = f.read()
                matches = findall(r'(.+?):\s*([\d.-]+)', content)
                if matches:
                    tmp_fea, tmp_score = zip(*matches)
                    tmp_fea = [f.strip() for f in tmp_fea]
                    tmp_score = [float(s) for s in tmp_score]
                    ps.append(tmp_score)
                    if fea is None:
                        fea = tmp_fea  
            pdbs.append(parsePDB(pdb_file))

        os.chdir(cwd0)

        max_len = max(len(row) for row in ps)
        ps_fixed = [row + [0.0] * (max_len - len(row)) for row in ps]
        ps_arr = np.array(ps_fixed)

        pcn = {int(pdb.getTitle().partition('_')[0][6:]):
            set(zip(pdb.getChids().tolist(), pdb.getResnums().tolist()))
            for pdb in pdbs}
        pi = {p: [rcr[x] for x in crn] for p, crn in pcn.items()}

        pzs_max = {k: max(self._zscores[v]) for k, v in pi.items()}
        pzs_med = {k: np.median(self._zscores[v]) for k, v in pi.items()}

        indices = Index(range(1, ps_arr.shape[0] + 1), name='Pocket #')
        columns = Index(fea, name='Feature')
        self._df = DataFrame(index=indices, columns=columns, data=ps_arr)

        lhd_col = [c for c in self._df.columns if 'Local hydrophobic density Score' in c]
        if lhd_col:
            lhd_vals = self._df[[lhd_col[0]]].values
            lhd_z = zscore(lhd_vals, axis=0, ddof=0)
        else:
            print("Warning: Local hydrophobic density feature not found. LHD set to 0.")
            lhd_z = np.zeros((ps_arr.shape[0], 1))

        essa_max_arr = np.array(list(pzs_max.values())).reshape(-1, 1)
        essa_med_arr = np.array(list(pzs_med.values())).reshape(-1, 1)
        zps = np.hstack((essa_max_arr, essa_med_arr, lhd_z))

        columns_zs = Index(['ESSA_max', 'ESSA_med', 'LHD'], name='Z-score')
        self._df_zs = DataFrame(index=indices, columns=columns_zs, data=zps)

        print("Pocket features and ESSA z-scores successfully computed.")

    def rankPockets(self):
        """
        Ranks pockets in terms of their allosteric potential based on ESSA z-scores
        (max/median) and local hydrophobic density (LHD) screening.
        """
        import numpy as np
        import pandas as pd

        if not hasattr(self, "_df_zs") or self._df_zs is None:
            raise RuntimeError("Pocket Z-scores not available. Run scanPockets() first.")

        lhd = self._df_zs.loc[:, "LHD"]
        n = np.count_nonzero(lhd >= 0.0)
        q = np.quantile(lhd, 0.85)

        s_max = ["ESSA_max", "LHD"]
        zf_max = self._df_zs[s_max].copy()
        if n >= lhd.size // 4:
            f_max = zf_max.iloc[:, 1] >= 0.0
        else:
            f_max = zf_max.iloc[:, 1] >= q
        zf_max = zf_max[f_max]
        zf_max.iloc[:, 0] = zf_max.iloc[:, 0].round(1)
        zf_max.iloc[:, 1] = zf_max.iloc[:, 1].round(2)
        self._idx_max = zf_max.sort_values(s_max, ascending=False).index

        s_med = ["ESSA_med", "LHD"]
        zf_med = self._df_zs[s_med].copy()
        if n >= lhd.size // 4:
            f_med = zf_med.iloc[:, 1] >= 0.0
        else:
            f_med = zf_med.iloc[:, 1] >= q
        zf_med = zf_med[f_med]
        zf_med.iloc[:, 0] = zf_med.iloc[:, 0].round(1)
        zf_med.iloc[:, 1] = zf_med.iloc[:, 1].round(2)
        self._idx_med = zf_med.sort_values(s_med, ascending=False).index

        ranks = pd.Index(range(1, len(self._idx_max) + 1), name="Rank")
        columns_ranks = pd.Index(
            ["Pocket # (ESSA_max & LHD)", "Pocket # (ESSA_med & LHD)"]
        )
        self._pocket_ranks = pd.DataFrame(index=ranks, columns=columns_ranks)
        self._pocket_ranks.iloc[:, 0] = self._idx_max
        self._pocket_ranks.iloc[:, 1] = self._idx_med

    def getPocketFeatures(self):
        """Return pocket features as a Pandas DataFrame."""
        return getattr(self, "_df", None)

    def getPocketZscores(self):
        """Return ESSA + LHD z-scores for pockets as a Pandas DataFrame."""
        return getattr(self, "_df_zs", None)

    def getPocketRanks(self):
        """Return pocket ranks (allosteric potential) as a Pandas DataFrame."""
        return getattr(self, "_pocket_ranks", None)

    def showPocketZscores(self):
        """Plot max/median ESSA and LHD z-scores for pockets."""
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import style, xticks, xlabel, ylabel, tight_layout

        if not hasattr(self, "_df_zs") or self._df_zs is None:
            raise RuntimeError("Pocket Z-scores not available. Run scanPockets() first.")

        with style.context(
            {
                "xtick.major.size": 10,
                "xtick.labelsize": 30,
                "ytick.major.size": 10,
                "ytick.labelsize": 30,
                "axes.labelsize": 35,
                "legend.fontsize": 25,
                "legend.title_fontsize": 0,
            }
        ):
            self._df_zs[["ESSA_max", "ESSA_med", "LHD"]].plot.bar(figsize=(25, 10))
            xticks(rotation=0)
            xlabel("Pocket")
            ylabel("Z-score")
            tight_layout()
            plt.show()

    def savePocketFeatures(self):
        """Save pocket features to a pickle `.pkl` file."""
        import pandas as pd

        if not hasattr(self, "_df") or self._df is None:
            raise RuntimeError("Pocket features not available. Run scanPockets() first.")
        self._df.to_pickle(f"{self._title}_pocket_features.pkl")

    def savePocketZscores(self):
        """Save ESSA + LHD z-scores for pockets to a pickle `.pkl` file."""
        if not hasattr(self, "_df_zs") or self._df_zs is None:
            raise RuntimeError("Pocket z-scores not available. Run scanPockets() first.")
        self._df_zs.to_pickle(f"{self._title}_pocket_zscores.pkl")

    def savePocketRanks(self):
        """Save pocket ranks to `.npy` files."""
        import numpy as np

        if not hasattr(self, "_idx_max") or not hasattr(self, "_idx_med"):
            raise RuntimeError("Pocket ranks not available. Run rankPockets() first.")

        np.save(f"{self._title}_pocket_ranks_wrt_ESSAmax_LHD.npy", self._idx_max)
        np.save(f"{self._title}_pocket_ranks_wrt_ESSAmed_LHD.npy", self._idx_med)

    def writePocketRanksToCSV(self):
        """Save pocket ranks to a CSV file."""
        if not hasattr(self, "_pocket_ranks") or self._pocket_ranks is None:
            raise RuntimeError("Pocket ranks not available. Run rankPockets() first.")

        self._pocket_ranks.to_csv(f"{self._title}_pocket_ranks.csv", index=False)

    