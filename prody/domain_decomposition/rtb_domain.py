import numpy as np
from scipy import sparse
from sklearn.neighbors import radius_neighbors_graph
from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
from prody import ANM, RTB, calcSqFlucts

__all__ = ['rtb_domain', 'plot_cc', 'plot_msf']


def rtb_domain(pdb,
               ndomains_l,
               ndomains_u,
               msf_other=None,
               n_modes=10,
               radius=10.,
               affinity=None,
               n_init=10,
               n_jobs=-1,
               method='discretize'):

    assert n_modes is None or isinstance(n_modes, int) and n_modes > 0, \
            'n_modes must be a positive integer or None'

    coo = pdb.getCoords()
    bfact = pdb.getBetas()
    labels = {}

    if msf_other is None:
        anm = ANM()
        anm.buildHessian(pdb)
        anm.calcModes(n_modes=n_modes)
        msf_other = calcSqFlucts(anm)

    if affinity is None:
        affinity = radius_neighbors_graph(coo, radius)
    else:
        affinity = affinity

    for n in range(ndomains_l, ndomains_u + 1):
        sc_pre = SpectralClustering(
            n_clusters=n, affinity='precomputed', n_init=n_init, n_jobs=n_jobs, assign_labels=method)
        sc_pre_labels = sc_pre.fit(affinity).labels_
        labels[n] = sc_pre_labels

    msf_rtb = {}
    cc_other_rtb = []
    cc_bfact_rtb = []

    rtb = RTB()

    for n in range(ndomains_l, ndomains_u + 1):
        rtb.buildHessian(coords=coo, blocks=labels[n])
        rtb.calcModes(n_modes=None)
        _n_modes = rtb._dof if n_modes is None else min(rtb._dof, n_modes) 
        msf_rtb[n] = calcSqFlucts(rtb[:_n_modes])
        cc_other_rtb.append(np.corrcoef(msf_other, msf_rtb[n])[0, 1])
        cc_bfact_rtb.append(np.corrcoef(bfact, msf_rtb[n])[0, 1])

    return np.array(cc_other_rtb), np.array(cc_bfact_rtb), msf_rtb, labels


def plot_cc(rtb, n_domains_l, n_domains_u):

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), dpi=100)
    ax[0].grid()
    ax[0].plot(np.arange(n_domains_l, n_domains_u + 1), rtb[0], '.-')
    ax[1].grid()
    ax[1].plot(np.arange(n_domains_l, n_domains_u + 1), rtb[1], '.-')
    ax[0].set_xlabel('# of domains')
    ax[0].set_ylabel("CC between Other's and RTBs' MSFs")

    ax[1].set_xlabel('# of domains')
    ax[1].set_ylabel('CC between Bfactors and MSFs of RTBs')

    plt.tight_layout()
    plt.show()


def plot_msf(other, rtb, n_domains_l, n_domains_u):

    plt.figure(figsize=(10, 5), dpi=100)
    plt.plot(other, label='Other MSFs')
    for n in range(n_domains_l, n_domains_u + 1):
        plt.plot(
            rtb[n] * np.mean(other) / np.mean(rtb[n]), label='{}'.format(n))
    plt.legend()
    plt.xlabel('# of res')
    plt.ylabel('MSFs of both other and RTBs')
    plt.tight_layout()
    plt.show()
