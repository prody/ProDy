import pickle
from os import chdir, listdir, mkdir
from os.path import isdir
from shutil import rmtree
import numpy as np
from scipy.stats import zscore, median_absolute_deviation
import matplotlib as mpl
import prody as pr


def ESSA(pdb, n_modes=20, s_modes=10, lig=None, cutoff_l=4.5, enm='gnm', cutoff=None):

    # pdb : a pdb parsed by ProDy
    # n_modes : n_modes of GNM object
    # s_modes : the number of modes to calculate the eigenvalue difference, s_modes <= n_modes
    # lig : list of ligands as a tuple of chainID and resSeq (resnum), e.g., [('A', 300), ('B', 301)]
    # cutoff_l : the protein residues within cutoff_l of ligands
    # enm: the type of enm
    # cutoff: enm cutoff


    # --- ESSA scanning part --- #

    heavy = pdb.select('protein and heavy and not hetatm')
    ca = heavy.ca

    # saveing pdb with a new name for Fpocket
    pr.writePDB(f'{pdb.getTitle()}_pro', heavy)

    # --- creating a directory to store models --- #

    direc = f'{pdb.getTitle()}_{enm}_data'
    if isdir(direc):
        rmtree(direc)
        mkdir(direc)
        chdir(direc)
    else:
        mkdir(direc)
        chdir(direc)

    # --- reference model --- #

    if enm == 'gnm':
        ca_enm = pr.GNM('ca')
        if cutoff is not None:
            ca_enm.buildKirchhoff(ca, cutoff=cutoff)
        else:
            ca_enm.buildKirchhoff(ca)

    if enm == 'anm':
        ca_enm = pr.ANM('ca')
        if cutoff is not None:
            ca_enm.buildHessian(ca, cutoff=cutoff)
        else:
            ca_enm.buildHessian(ca)

    ca_enm.calcModes(n_modes=n_modes)
    pr.saveModel(ca_enm, matrices=False, filename=pdb.getTitle() + '_ca_20')

    # --- perturbed models --- #

    pr.LOGGER.progress(msg='', steps=(ca.numAtoms()))
    for i in ca.getResindices():
        pr.LOGGER.update(step=i+1, msg=f'scanning residue {i+1}')
        sel = f'resindex {i} or (calpha and resindex != {i})'
        tmp = heavy.select(sel)

        if enm == 'gnm':
            tmp_enm = pr.GNM(f'res_{i}')
            if cutoff is not None:
                tmp_enm.buildKirchhoff(tmp, cutoff=cutoff)
            else:
                tmp_enm.buildKirchhoff(tmp)

        if enm == 'anm':
            tmp_enm = pr.ANM(f'res_{i}')
            if cutoff is not None:
                tmp_enm.buildHessian(tmp, cutoff=cutoff)
            else:
                tmp_enm.buildHessian(tmp)

        tmp_enm_red, _ = pr.reduceModel(tmp_enm, tmp, ca)
        tmp_enm_red.calcModes(n_modes=n_modes)
        pr.saveModel(tmp_enm_red, matrices=False, filename=pdb.getTitle() + f'_red_20_{i}')

    chdir('..')
    print('\n')

    # --- ESSA computation part --- #

    def beta_heavy(*args):

        # args[0] : vector
        # args[1] : nodes
        # args[2] : atoms

        tmp0 = pr.Vector(args[0], is3d=False)
        tmp1, _ = pr.extendVector(tmp0, args[1], args[2])

        return tmp1.getArray()

    def modified_zscore(arg):

        # possbily be deprecated

        return 0.6745 * (arg - np.median(arg)) / median_absolute_deviation(arg)


    #-----------#

    chdir(f'{pdb.getTitle()}_{enm}_data')

    ml = [x for x in listdir('.') if x.endswith(f'.{enm}.npz')]
    ml.sort()

    ref = pr.loadModel(ml[0])

    redm = ml[1:]
    redm.sort(key=lambda x: int(x.partition('.')[0].rpartition('_')[-1]))

    red = [pr.loadModel(x) for x in redm]

    chdir('..')

    #-----------#

    mn = len(red)
    red_mm = []
    for i in range(mn):
        tmp = pr.matchModes(ref[:s_modes], red[i][:s_modes])
        red_mm.append(tmp)

    #-----------#

    diff_red_ca_mm = np.array([(red_mm[i][1].getEigvals() - red_mm[0][0].getEigvals()) / red_mm[0][0].getEigvals() * 100 for i in range(mn)])

    diff_red_ca_mm_m = diff_red_ca_mm.mean(1)

    #-----------#

    zs = zscore(diff_red_ca_mm_m)

    # automatically save zscores both as a numpy array and as a pdb
    # the numpy array of zscores will be used by pocket function later

    np.save(f'{pdb.getTitle()}_{enm}_zs', zs)
    pr.writePDB(f'{pdb.getTitle()}_{enm}_zs', heavy, beta=beta_heavy(zs, ca, heavy))

    # possibly be deprecated
    zsm = modified_zscore(diff_red_ca_mm_m)
    np.save(f'{pdb.getTitle()}_{enm}_mzs', zsm)

    # --- residue indices of protein residues that are within 4.5 A of ligands --- #

    if lig is not None:
        idx_lig = {''.join(chid + str(resnum)): pdb.select(f'calpha and not hetatm and (same residue as exwithin 4.5 of (chain {chid} and resnum {resnum}))').getResindices()
                for chid, resnum in lig}

        pickle.dump(idx_lig, open(f'{pdb.getTitle()}_resindices_ligs.pkl', 'wb'))

    #-----------#

    if lig:
        return zs, idx_lig
    else:
        return zs


def ESSA_profile(pdb, zscore, lig=None, quantile=.75, enm='gnm', save=False):

    # pdb : pdb parsed by ProDy, it should be the same as the one used in ESSA
    # zscore: zscore as a numpy array
    # lig: dictionary containing resindices of protein residues within 4.5 A of ligands, e.g.
    #      {('A', 300): np.array([0, 15, 78, 79, 83]), ('B', 201): np.array([135, 137, 139])}
    # quantile: drawing quantule hlines at that specified value

    ca = pdb.select('calpha and not hetatm')
    can = ca.numAtoms()

    # plotting confs can also be parameters
    # let's not change rcParams globaly
    with mpl.rc_context({'axes.labelsize': 'xx-large', 'xtick.labelsize': 'x-large',
                         'ytick.labelsize': 'x-large', 'legend.fontsize': 'large',
                         'figure.figsize': (15, 7), 'figure.dpi': 600}):

        mpl.pyplot.figure()
        mpl.pyplot.plot(zscore, 'k', linewidth=1.)

        if lig:
            zs_lig = {k: zscore[v] for k, v in lig.items()}
            for k in lig.keys():
                mpl.pyplot.scatter(lig[k], zs_lig[k], s=100, label=k)
            # legends : (chainID & resSeq)   resSeq (resnum), not resindex
            mpl.pyplot.legend(fontsize='large')

        mpl.pyplot.hlines(np.quantile(zscore, q=quantile),
                         xmin=0., xmax=can, linewidth=3., linestyle='--', color='c')

        # if there is more than one chain, then let's put a domainbar
        if len(set(ca.getChids())) > 1:
            pr.showDomainBar(ca.getChids(), fontdict={'size': 'x-large'})

        # resindices as labels, not (chanID, resSeq),
        # but we have a domainbar if the number of chains > 1
        mpl.pyplot.xlabel('Residue')
        mpl.pyplot.ylabel('Z-Score')

        mpl.pyplot.tight_layout()

        if save:
            mpl.pyplot.savefig(f'{pdb.getTitle()}_{enm}_zs')
            mpl.pyplot.close()
        else:
            mpl.pyplot.show()