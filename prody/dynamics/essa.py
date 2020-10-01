import pickle
from os import chdir, listdir, mkdir
from os.path import isdir
from shutil import rmtree
from numpy import array, mean, median, quantile, save
from scipy.stats import zscore, median_absolute_deviation
from matplotlib import rc_context, pyplot

from prody import LOGGER
from .anm import ANM
from .gnm import GNM
# from prody.proteins import parsePDB
from prody.proteins import writePDB
# from .compare import matchModes
from .editing import reduceModel, extendVector
# from .functions import loadModel, saveModel
from .mode import Vector
from .plotting import showDomainBar
from .signature import ModeEnsemble, loadModeEnsemble, saveModeEnsemble 

__all__ = ['ESSA', 'showESSAprofile']


def ESSA(pdb, n_modes=20, s_modes=10, lig=None,
         dist=4.5, enm='gnm', cutoff=None):

    # pdb : a pdb parsed by ProDy
    # n_modes : n_modes of GNM object
    # s_modes : the number of modes to calculate the eigenvalue difference, s_modes <= n_modes
    # lig : string of ligands' chainIDs and resSeqs (resnum) separated by a whitespace, e.g., 'A 300 B 301'
    # dist : the protein residues within a distance of ligands
    # enm: the type of enm
    # cutoff: enm cutoff

    # --- ESSA scanning part --- #

    heavy = pdb.select('protein and heavy and not hetatm')
    ca = heavy.ca

    # saving pdb with a new name for Fpocket

    writePDB(f'{pdb.getTitle()}_pro', heavy)

    # --- reference model --- #

    if enm == 'gnm':
        ca_enm = GNM('ca')
        if cutoff is not None:
            ca_enm.buildKirchhoff(ca, cutoff=cutoff)
        else:
            ca_enm.buildKirchhoff(ca)

    if enm == 'anm':
        ca_enm = ANM('ca')
        if cutoff is not None:
            ca_enm.buildHessian(ca, cutoff=cutoff)
        else:
            ca_enm.buildHessian(ca)

    ca_enm.calcModes(n_modes=n_modes)

    ens = ModeEnsemble(f'{pdb.getTitle()}')
    ens.setAtoms(ca)
    ens.addModeSet(ca_enm[:])
    labels = ['ref']

    # --- perturbed models --- #

    LOGGER.progress(msg='', steps=(ca.numAtoms()))
    for i in ca.getResindices():
        LOGGER.update(step=i+1, msg=f'scanning residue {i+1}')
        sel = f'calpha or resindex {i}'
        tmp = heavy.select(sel)

        if enm == 'gnm':
            tmp_enm = GNM(f'res_{i}')
            if cutoff is not None:
                tmp_enm.buildKirchhoff(tmp, cutoff=cutoff)
            else:
                tmp_enm.buildKirchhoff(tmp)

        if enm == 'anm':
            tmp_enm = ANM(f'res_{i}')
            if cutoff is not None:
                tmp_enm.buildHessian(tmp, cutoff=cutoff)
            else:
                tmp_enm.buildHessian(tmp)

        tmp_enm_red, _ = reduceModel(tmp_enm, tmp, ca)
        tmp_enm_red.calcModes(n_modes=n_modes)

        ens.addModeSet(tmp_enm_red[:])
        labels.append(tmp_enm.getTitle())

    ens.setLabels(labels)
    ens.match()

    # --- ESSA computation part --- #

    def beta_heavy(*args):

        # args[0] : vector
        # args[1] : nodes
        # args[2] : atoms

        tmp0 = Vector(args[0], is3d=False)
        tmp1, _ = extendVector(tmp0, args[1], args[2])

        return tmp1.getArray()

    denom = ens[0].getEigvals()
    num = ens[1:].getEigvals() - denom

    eig_diff = num / denom * 100
    eig_diff_mean = mean(eig_diff, axis=1)

    #-----------#

    zs = zscore(eig_diff_mean)

    # automatically save zscores both as a numpy array and as a pdb
    # the numpy array of zscores will be used by pocket function later

    save(f'{pdb.getTitle()}_{enm}_zs', zs)
    writePDB(f'{pdb.getTitle()}_{enm}_zs', heavy,
             beta=beta_heavy(zs, ca, heavy))

    # --- residue indices of protein residues that are within dist (4.5 A) of ligands --- #

    if lig is not None:
        ligs = lig.split()
        ligs = list(zip(ligs[::2], ligs[1::2]))
        idx_lig = {}
        for chid, resnum in ligs:
            key = ''.join(chid + str(resnum))
            sel_lig = 'calpha and not hetatm and (same residue as ' \
                      f'exwithin {dist} of (chain {chid} and resnum {resnum}))'
            idx_lig[key] = pdb.select(sel_lig).getResindices()

        pickle.dump(idx_lig,
                    open(f'{pdb.getTitle()}_resindices_ligs.pkl', 'wb'))

    #-----------#

    if lig:
        return zs, idx_lig
    else:
        return zs


def showESSAprofile(pdb, zscore, lig=None, quant=.75, enm='gnm', save=False):

    # pdb : pdb parsed by ProDy, it should be the same as the one used in ESSA
    # zscore: zscore as a numpy array
    # lig: dictionary containing resindices of protein residues within 4.5 A of ligands, e.g.
    #      {('A', 300): np.array([0, 15, 78, 79, 83]), ('B', 201): np.array([135, 137, 139])}
    # quantile: drawing quantule hlines at that specified value

    ca = pdb.select('calpha and not hetatm')
    can = ca.numAtoms()

    # plotting confs can also be parameters
    # let's not change rcParams globaly
    with rc_context({'axes.labelsize': 'xx-large',
                     'xtick.labelsize': 'x-large',
                     'ytick.labelsize': 'x-large',
                     'legend.fontsize': 'large',
                     'figure.figsize': (15, 7),
                     'figure.dpi': 600}):

        pyplot.figure()
        pyplot.plot(zscore, 'k', linewidth=1.)

        if lig:
            zs_lig = {k: zscore[v] for k, v in lig.items()}
            for k in lig.keys():
                pyplot.scatter(lig[k], zs_lig[k], s=100, label=k)
            pyplot.legend(fontsize='large')

        pyplot.hlines(quantile(zscore, q=quant),
                      xmin=0., xmax=can,
                      linewidth=3., linestyle='--', color='c')

        # if there is more than one chain, then let's put a domainbar
        if len(set(ca.getChids())) > 1:
            showDomainBar(ca.getChids(), fontdict={'size': 'x-large'})

        # resindices as labels, not (chanID, resSeq),
        # but we have a domainbar if the number of chains > 1
        pyplot.xlabel('Residue')
        pyplot.ylabel('Z-Score')

        pyplot.tight_layout()

        if save:
            pyplot.savefig(f'{pdb.getTitle()}_{enm}_zs')
            pyplot.close()
        else:
            pyplot.show()
