# -*- coding: utf-8 -*-
"""This module defines :func:`.calcNeoAdaptiveANM`, which generates a
conformational transition path between two structures in which every accepted
step is physical.

Like adaptive ANM (Yang, Majek and Bahar, *PLoS Comput Biol* 2009;
:func:`.calcAdaptiveANM`), it deforms one endpoint
along the soft modes of an elastic network rebuilt at every step. NeoAdaptiveANM
adds, all on by default and each switchable with one keyword:

* a physical network: the generalized ANM (:class:`.genANM`, ``model='genanm'``)
  and stiff contact springs between approaching non-bonded residues;
* soft modes first: a mode pool that opens with progress, modes taken softest
  first, and a stiffness filter that damps stiffer modes;
* physical steps: a bond guard, a steric test, lower limits on virtual bonds and on
  residues two apart, and geometric correction after every accepted step;
* recovery: a backtracking line search and trajectory rewind.

:func:`.calcAdaptiveANM` is unchanged; this module adds a separate entry point.
"""

__author__ = 'Anupam Banerjee'
__credits__ = ['Anupam Banerjee']
__email__ = ['anupam.banerjee@stonybrook.edu']

from numbers import Integral
import numpy as np
from datetime import datetime

from prody import LOGGER
import logging

from prody.atomic import Atomic, AtomMap
from prody.ensemble import Ensemble
from prody.measure import calcRMSD, superpose
from prody.utilities import getCoords, importLA

from .functions import calcENM
from .generalized_anm import GENANM_MODEL_NAMES
from .adaptive import AANM_ONEWAY, AANM_ALTERNATING, AANM_BOTHWAYS


def _rmsd_to(coords, target, weights=None, superposed=False):
    """RMSD of *coords* to *target*. With *superposed* the RMSD is taken after the optimal superposition of a copy of
    *coords* onto *target*, so it measures shape only; the coordinates themselves are never moved. ANM modes carry no
    rigid-body motion, so a rotation that builds up between the moving structure and the target cannot be removed by
    later steps; measured without superposition it would keep the RMSD above the target and stall the run."""
    if superposed:
        coords = superpose(np.array(coords, dtype=float), np.asarray(target, dtype=float), weights)[0]
    return calcRMSD(coords, target, weights)

# Defaults: these reproduce the published NeoAdaptiveANM configuration.
DEFAULT_ADJ_MAX = 6.0      # bond guard: longest allowed sequential Ca-Ca distance (A); native breaks exempt
DEFAULT_F = 0.2            # step factor
DEFAULT_MIN_RMSD_DIFF = 0.0       # stall test on the last RMSD changes (0 = off)
DEFAULT_MIN_RMSD_STEP = 0.0       # smallest accepted RMSD step, A
DEFAULT_MIN_RMSD_IMPROVE = 5e-3   # stop when a step improves the RMSD by less than this, A
DEFAULT_GEOM_BOND_LO = 3.4   # geometric correction and step test: shortest allowed virtual bond (A), native-aware
DEFAULT_GEOM_D13 = 4.4       # shortest allowed distance between residues i and i+2 (A), native-aware
DEFAULT_CLASH_MIN = 3.0     # steric test on by default; pass clash_min=0 to turn it off
DEFAULT_GEOM_BOND_HI = 4.2  # geometric correction on by default; pass geom_bond_hi=0 and geom_clash=0 to turn it off
DEFAULT_GEOM_CLASH = 3.2
DEFAULT_CONTACT_ON = 4.5    # contact-aware network on by default; pass contact_on=0 to turn it off
DEFAULT_CONTACT_K = 10.0
DEFAULT_SOFT_FILTER = 3.0         # stiffness filter on by default (soft modes first); soft_filter=0 gives the standard driver
DEFAULT_SOFT_FILTER_ADAPT = True  # relax the filter on steps that cannot progress without stiffer contributions

__all__ = ['calcNeoAdaptiveANM']


norm = importLA().norm

# -------------------------------------------------------------------
# helpers: geometry / adjacency
# -------------------------------------------------------------------

def sequential_edge_lengths(coords):
    r = np.asarray(coords, float)
    if r.shape[0] < 2:
        return np.zeros(0, dtype=float)
    diffs = r[1:] - r[:-1]
    return np.sqrt(np.sum(diffs * diffs, axis=1))


def build_adj_exempt(native_coords, adj_max=DEFAULT_ADJ_MAX, tol=0.0):
    """
    Return boolean mask for sequential edges that are already longer than
    adj_max (+tol) in the native structure. These edges are exempt from
    adjacency checks.
    """
    elen_native = sequential_edge_lengths(native_coords)
    if elen_native.size == 0:
        return None
    threshold = float(adj_max) + float(tol)
    return (elen_native > threshold)


def compute_slack_info(coords, adj_max=DEFAULT_ADJ_MAX, risk_topk=8, adj_exempt=None):
    """Return (slack array, min_slack, worst_edge_indices[:topk])."""
    elen  = sequential_edge_lengths(coords)
    slack = float(adj_max) - elen
    if slack.size == 0:
        return slack, np.inf, np.array([], dtype=int)
    if adj_exempt is not None and adj_exempt.shape == slack.shape:
        slack = slack.copy()
        slack[adj_exempt] = np.inf
    order = np.argsort(slack)
    worst = order[:min(risk_topk, order.size)]
    return slack, float(slack.min()), worst


def checkDisconnection(coords, cutoff, adj_max=DEFAULT_ADJ_MAX, adj_exempt=None):
    """
    True if there is global isolation (>cutoff) or any non-exempt
    sequential edge > adj_max.
    """
    r = np.asarray(coords, float)
    N = r.shape[0]
    if N < 2:
        return False

    all_dists = np.sqrt(((r[:, None, :] - r[None, :, :]) ** 2).sum(axis=2))
    nn = np.empty(N, dtype=float)
    for i in range(N):
        row  = all_dists[i]
        nn[i] = min(
            np.min(row[:i])     if i > 0     else np.inf,
            np.min(row[i + 1:]) if i + 1 < N else np.inf
        )
    isolated = nn > float(cutoff)

    broken_adj = False
    if N >= 2:
        edge_d = sequential_edge_lengths(coords)
        if adj_exempt is not None and adj_exempt.shape == edge_d.shape:
            broken_adj = np.any((edge_d > float(adj_max)) & (~adj_exempt))
        else:
            broken_adj = np.any(edge_d > float(adj_max))

    return bool(np.any(isolated) or broken_adj)


# -------------------------------------------------------------------
# steric (excluded-volume) test on non-bonded Calpha pairs
# -------------------------------------------------------------------
def _close_pairs(coords, r, sep):
    """Index pairs (i < j, j - i >= sep) with Calpha-Calpha distance < r."""
    X = np.asarray(coords, float)
    try:
        from scipy.spatial import cKDTree
        P = np.asarray(sorted(cKDTree(X).query_pairs(float(r))), dtype=int).reshape(-1, 2)
    except ImportError:
        i, j = np.triu_indices(len(X), k=1)
        m = np.linalg.norm(X[i] - X[j], axis=1) < float(r)
        P = np.stack([i[m], j[m]], axis=1)
    if P.size == 0:
        return P.reshape(0, 2)
    P = np.sort(P, axis=1)
    return P[(P[:, 1] - P[:, 0]) >= int(sep)]





def _opt(kwargs, name, default):
    """Option value from kwargs; a missing key or None means the module default."""
    v = kwargs.get(name, None)
    return default if v is None else v


def _resolve_defaults(kwargs):
    """Write the effective steric, geometric-correction and contact settings into kwargs."""
    kwargs['clash_min'] = float(_opt(kwargs, 'clash_min', DEFAULT_CLASH_MIN))
    kwargs['geom_bond_hi'] = float(_opt(kwargs, 'geom_bond_hi', DEFAULT_GEOM_BOND_HI))
    kwargs['geom_clash'] = float(_opt(kwargs, 'geom_clash', DEFAULT_GEOM_CLASH))
    kwargs['geom_bond_lo'] = float(_opt(kwargs, 'geom_bond_lo', DEFAULT_GEOM_BOND_LO))
    kwargs['geom_d13'] = float(_opt(kwargs, 'geom_d13', DEFAULT_GEOM_D13))
    kwargs['contact_on'] = float(_opt(kwargs, 'contact_on', DEFAULT_CONTACT_ON))
    kwargs['contact_k'] = float(_opt(kwargs, 'contact_k', DEFAULT_CONTACT_K))
    kwargs['soft_filter'] = float(_opt(kwargs, 'soft_filter', DEFAULT_SOFT_FILTER))
    kwargs['soft_filter_adapt'] = bool(_opt(kwargs, 'soft_filter_adapt', DEFAULT_SOFT_FILTER_ADAPT))
    return kwargs


def _clash_min_of(kwargs):
    """Steric threshold from kwargs, defaulting to DEFAULT_CLASH_MIN (None also means default)."""
    v = kwargs.get('clash_min', None)
    return DEFAULT_CLASH_MIN if v is None else float(v)


def build_chain_floors(native_coords_list, bond_lo=DEFAULT_GEOM_BOND_LO, d13=DEFAULT_GEOM_D13):
    """
    Native-aware lower limits for the chain: for every virtual bond (i, i+1) the floor is the smaller of bond_lo and
    its length in any native structure, and likewise for every pair (i, i+2) with d13. Native geometry (for example a
    cis peptide, about 2.9 A) is therefore never treated as a violation. Returns (bond_floor, d13_floor); either is
    None when its limit is switched off (<= 0).
    """
    X = [np.asarray(c, float) for c in native_coords_list]
    bf = df = None
    if bond_lo and bond_lo > 0:
        bf = np.minimum(bond_lo, np.min([np.linalg.norm(x[1:] - x[:-1], axis=1) for x in X], axis=0))
    if d13 and d13 > 0 and len(X[0]) > 2:
        df = np.minimum(d13, np.min([np.linalg.norm(x[2:] - x[:-2], axis=1) for x in X], axis=0))
    return bf, df


def build_bond_ceiling(native_coords_list, bond_hi=4.2):
    """
    Native-aware upper limit for every virtual bond: the larger of bond_hi and its length in any native structure,
    so native chain breaks and long native bonds are never treated as violations. None when bond_hi <= 0.
    """
    if not bond_hi or bond_hi <= 0: return None
    X = [np.asarray(c, float) for c in native_coords_list]
    return np.maximum(bond_hi, np.max([np.linalg.norm(x[1:] - x[:-1], axis=1) for x in X], axis=0))


def checkChainGeometry(coords, bond_floor=None, d13_floor=None, tol=0.05, bond_ceiling=None, exempt_floor=None):
    """
    Chain test. Returns (violation, kind, amount): violation is True if any virtual bond is shorter than its floor or
    longer than its ceiling, any pair (i, i+2) closer than its floor, or any steric-exempt pair closer than its
    native distance (exempt_floor, from build_exempt_floors), by more than tol (A).
    """
    X = np.asarray(coords, float)
    for kind, floor, k in (("bond", bond_floor, 1), ("i,i+2", d13_floor, 2)):
        if floor is None: continue
        short = floor - np.linalg.norm(X[k:] - X[:-k], axis=1)
        if short.max() > tol:
            return True, kind, float(short.max())
    if bond_ceiling is not None:
        over = np.linalg.norm(X[1:] - X[:-1], axis=1) - bond_ceiling
        if over.max() > tol:
            return True, "bond above its ceiling", float(over.max())
    if exempt_floor is not None:
        P, fl = exempt_floor
        short = fl - np.linalg.norm(X[P[:, 0]] - X[P[:, 1]], axis=1)
        if short.max() > tol:
            return True, "native close pair below its native distance", float(short.max())
    return False, None, 0.0


def build_clash_exempt(native_coords_list, clash_min=3.0, clash_sep=3):
    """
    Non-bonded pairs already closer than clash_min in any native structure.
    These are exempt from the steric test, so genuine tight contacts present
    in the endpoints are not treated as clashes (mirrors build_adj_exempt).
    """
    ex = set()
    for X in native_coords_list:
        for i, j in _close_pairs(X, clash_min, clash_sep):
            ex.add((int(i), int(j)))
    return ex


def build_exempt_floors(native_coords_list, exempt):
    """
    Native floor for each steric-exempt pair: its shortest distance in any native structure. An exempt pair is a
    genuine tight contact, so it is allowed to stay that close but not to come closer. Returns (pairs, floors) as
    arrays, or None when there are no exempt pairs.
    """
    if not exempt: return None
    P = np.array(sorted(exempt), dtype=int).reshape(-1, 2)
    X = [np.asarray(c, float) for c in native_coords_list]
    return P, np.min([np.linalg.norm(x[P[:, 0]] - x[P[:, 1]], axis=1) for x in X], axis=0)


def contact_pairs(coords, contact_on, contact_sep=3, contact_exempt=None):
    """Non-bonded pairs (sequence separation >= contact_sep) currently closer than
    contact_on, excluding pairs already that close in a native structure."""
    P = _close_pairs(coords, contact_on, contact_sep)
    if contact_exempt:
        P = np.array([p for p in P if (int(p[0]), int(p[1])) not in contact_exempt], dtype=int).reshape(-1, 2)
    return P


def _calcENM_contact(coords, contact_on=0.0, contact_k=10.0, contact_sep=3,
                     contact_exempt=None, **kw):
    """
    Contact-aware network. Builds the model with calcENM, then, if contact_on > 0,
    adds a stiff spring (stiffness contact_k, rest length = current distance) between
    every non-bonded pair that has come within contact_on, and recomputes the modes.
    Motions that would push these residues together become stiff and leave the soft
    subspace, while sliding past each other stays soft, so steps built from the soft
    modes route around the obstacle. Each spring has the ANM form k u u^T, so the six
    rigid-body zero modes are preserved.
    """
    enm, atoms = calcENM(coords, **kw)
    if not contact_on or float(contact_on) <= 0.0 or kw.get('select') is not None:
        return enm, atoms
    X = np.asarray(coords, float)
    P = contact_pairs(X, float(contact_on), int(contact_sep), contact_exempt)
    if P.size == 0:
        return enm, atoms
    H = enm.getHessian()
    H = np.array(H.toarray() if hasattr(H, 'toarray') else H, dtype=float, copy=True)
    k = float(contact_k)
    for i, j in P:
        u = X[j] - X[i]; u = u / (np.linalg.norm(u) + 1e-12)
        B = k * np.outer(u, u); a, b = 3 * int(i), 3 * int(j)
        H[a:a+3, a:a+3] += B; H[b:b+3, b:b+3] += B
        H[a:a+3, b:b+3] -= B; H[b:b+3, a:a+3] -= B
    enm.setHessian(H)
    enm.calcModes(n_modes=kw.get('n_modes'), zeros=False)
    LOGGER.info("contact: %d pair springs (k=%.1f, within %.2f A) added before diagonalization"
                % (len(P), k, float(contact_on)))
    return enm, atoms


def geometric_correct(X, bond_hi=4.2, clash_d=3.2, clash_sep=3, adj_exempt=None,
                      clash_exempt=None, iters=200, tol=1e-3, relax=0.5, bond_floor=None, d13_floor=None,
                      exempt_floor=None):
    """
    Position-based constraint satisfaction applied to a candidate step. Iteratively:
      - virtual bonds (i, i+1) longer than bond_hi are pulled back to bond_hi
        (native chain breaks, flagged in adj_exempt, are left alone);
      - non-bonded pairs (sequence separation >= clash_sep) closer than clash_d are
        pushed apart to clash_d (pairs close in a native structure are exempt);
      - virtual bonds shorter than bond_floor, and pairs (i, i+2) closer than d13_floor, are pushed
        out to their floors (native-aware, see build_chain_floors).
    Each correction moves both atoms along their connecting vector; the result is recentred so
    the centre of mass is exactly unchanged. Returns
    (corrected coords, rms correction, remaining worst violation).
    Setting bond_hi <= 0 or clash_d <= 0 disables that part.
    """
    Y = np.array(X, dtype=float, copy=True); N = len(Y)
    worst = 0.0
    for _ in range(int(iters)):
        dY = np.zeros_like(Y); cnt = np.zeros(N); worst = 0.0
        if bond_hi and bond_hi > 0:
            v = Y[1:] - Y[:-1]; L = np.linalg.norm(v, axis=1)
            bad = L > bond_hi
            if adj_exempt is not None and np.shape(adj_exempt) == L.shape:
                bad &= ~np.asarray(adj_exempt, bool)
            if bad.any():
                k = np.where(bad)[0]; ex = (L[k] - bond_hi)
                worst = max(worst, float(ex.max()))
                u = v[k] / L[k][:, None]; step = 0.5 * ex[:, None] * u
                np.add.at(dY, k, step); np.add.at(dY, k + 1, -step)
                np.add.at(cnt, k, 1); np.add.at(cnt, k + 1, 1)
        if exempt_floor is not None:                              # native close pairs: no closer than native
            Pe, fe = exempt_floor
            v = Y[Pe[:, 1]] - Y[Pe[:, 0]]; L = np.linalg.norm(v, axis=1) + 1e-9
            sh = fe - L; k = np.where(sh > 0)[0]
            if k.size:
                worst = max(worst, float(sh[k].max()))
                u = v[k] / L[k][:, None]; step = 0.5 * sh[k][:, None] * u
                np.add.at(dY, Pe[k, 0], -step); np.add.at(dY, Pe[k, 1], step)
                np.add.at(cnt, Pe[k, 0], 1); np.add.at(cnt, Pe[k, 1], 1)
        for floor, kk in ((bond_floor, 1), (d13_floor, 2)):       # lower limits along the chain
            if floor is None: continue
            v = Y[kk:] - Y[:-kk]; L = np.linalg.norm(v, axis=1) + 1e-9
            sh = floor - L; k = np.where(sh > 0)[0]
            if k.size:
                worst = max(worst, float(sh[k].max()))
                u = v[k] / L[k][:, None]; step = 0.5 * sh[k][:, None] * u
                np.add.at(dY, k, -step); np.add.at(dY, k + kk, step)
                np.add.at(cnt, k, 1); np.add.at(cnt, k + kk, 1)
        if clash_d and clash_d > 0:
            P = _close_pairs(Y, clash_d, clash_sep)
            if clash_exempt:
                P = np.array([p for p in P if (int(p[0]), int(p[1])) not in clash_exempt], dtype=int).reshape(-1, 2)
            if P.size:
                v = Y[P[:, 1]] - Y[P[:, 0]]; L = np.linalg.norm(v, axis=1) + 1e-9
                sh = clash_d - L; worst = max(worst, float(sh.max()))
                u = v / L[:, None]; step = 0.5 * sh[:, None] * u
                np.add.at(dY, P[:, 0], -step); np.add.at(dY, P[:, 1], step)
                np.add.at(cnt, P[:, 0], 1); np.add.at(cnt, P[:, 1], 1)
        if worst < tol:
            break
        m = cnt > 0
        Y[m] += relax * dY[m] / cnt[m][:, None]
    Y -= Y.mean(axis=0) - np.asarray(X, float).mean(axis=0)     # keep the centre of mass fixed
    rms = float(np.sqrt(((Y - np.asarray(X, float)) ** 2).sum(axis=1).mean()))
    return Y, rms, worst


def checkClash(coords, clash_min=3.0, clash_sep=3, clash_exempt=None):
    """
    Steric test. Returns (clash, d_min): clash is True if any non-exempt
    non-bonded Calpha pair (sequence separation >= clash_sep) is closer than
    clash_min; d_min is the shortest such offending distance (or None).
    """
    P = _close_pairs(coords, clash_min, clash_sep)
    if clash_exempt:
        P = np.array([p for p in P if (int(p[0]), int(p[1])) not in clash_exempt], dtype=int).reshape(-1, 2)
    if P.size == 0:
        return False, None
    X = np.asarray(coords, float)
    d = np.linalg.norm(X[P[:, 0]] - X[P[:, 1]], axis=1)
    return True, float(d.min())

# -------------------------------------------------------------------
# helpers: input / scheduling / convergence
# -------------------------------------------------------------------

def checkInput(a, b, **kwargs):
    coordsA = getCoords(a)
    if isinstance(a, Atomic):
        title = a.getTitle()
        atoms = a
    else:
        title = None
        atoms = None

    coordsB = getCoords(b)
    if title is None:
        if isinstance(b, Atomic):
            title = b.getTitle()
            atoms = b
        else:
            title = 'Unknown'
            atoms = None

    maskA = a.getFlags("mapped") if isinstance(a, AtomMap) else 1.0
    maskB = b.getFlags("mapped") if isinstance(b, AtomMap) else 1.0
    weights = maskA * maskB

    if np.isscalar(weights):
        weights = None
    if np.isscalar(maskA):
        maskA = None
    if np.isscalar(maskB):
        maskB = None

    if not kwargs.get('aligned', False):
        coordsA, _ = superpose(coordsA, coordsB, weights)

    rmsd = _rmsd_to(coordsA, coordsB, weights, bool(kwargs.get('superposed_rmsd', False)))
    LOGGER.info('Initialized Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd))

    return coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd


def getTitle(structure, def_title='structure'):
    return structure.getTitle() if isinstance(structure, Atomic) else def_title


def computeAdaptiveF(defvecs, f_max=1.0, f_min=0.2, gamma=1.0):
    """Optional global f scheduler (kept for compatibility, f_auto=True only)."""
    if not defvecs:
        return float(f_max)
    d0, dk   = defvecs[0], defvecs[-1]
    n0, nk   = norm(d0), norm(dk)
    if n0 <= 1e-12:
        return float(f_min)
    progress = 1.0 - np.sqrt(max(nk, 0.0) / max(n0, 1e-12))
    progress = float(np.clip(progress, 0.0, 1.0))
    return float(f_min + (f_max - f_min) * (1.0 - (progress ** float(gamma))))


def computeProgressiveNmax(defvecs, n_max_modes, frac0=0.2, eta=1.5):
    """
    Progressive cap for high-frequency modes (legacy, soft_pool=False only).
    Cap starts at n_max_modes * frac0 and grows toward n_max_modes as the
    defvec shrinks.
    """
    n_max_modes = int(max(1, n_max_modes))
    if n_max_modes <= 1:
        return n_max_modes
    if not defvecs:
        return int(max(1, int(np.floor(n_max_modes * float(frac0)))))
    d0, dk = defvecs[0], defvecs[-1]
    n0, nk = norm(d0), norm(dk)
    if n0 <= 1e-12:
        progress = 1.0
    else:
        progress = 1.0 - np.sqrt(max(nk, 0.0) / max(n0, 1e-12))
    progress = float(np.clip(progress, 0.0, 1.0))
    frac     = float(frac0) + (1.0 - float(frac0)) * (progress ** float(eta))
    eff_n    = int(np.floor(n_max_modes * frac))
    return int(max(1, min(n_max_modes, eff_n)))


def checkConvergence(rmsds, coords, **kwargs):
    min_rmsd_diff = kwargs.get('min_rmsd_diff', DEFAULT_MIN_RMSD_DIFF)
    target_rmsd   = kwargs.get('target_rmsd', 1.0)
    cutoff        = kwargs.get('cutoff', 15.0)
    adj_max       = kwargs.get('adj_max', DEFAULT_ADJ_MAX)
    adj_exempt    = kwargs.get('adj_exempt', None)

    if len(rmsds) > 4 and np.all(np.abs(np.diff(rmsds[-4:])) < min_rmsd_diff):
        LOGGER.warning(f'RMSD change fell below {min_rmsd_diff}')
        return True

    if rmsds[-1] < target_rmsd:
        LOGGER.info(
            "Convergence reached: RMSD %.3f is below target %.3f." %
            (rmsds[-1], target_rmsd)
        )
        return True

    if checkDisconnection(coords, cutoff, adj_max=adj_max, adj_exempt=adj_exempt):
        LOGGER.warning('Disconnection found in the structure.')
        return True

    return False

# -------------------------------------------------------------------
# RMSD-based Fmin schedule
# -------------------------------------------------------------------

def _compute_Fmin_from_RMSD(rmsd_current, rmsd0, Fmin_min=0.4, Fmin_max=0.9):
    """
    - At RMSD = rmsd0: Fmin ~ Fmin_min
    - At RMSD -> 0:    Fmin -> Fmin_max
    """
    if rmsd0 <= 1e-8:
        progress = 1.0
    else:
        progress = 1.0 - float(rmsd_current / rmsd0)
    progress = float(np.clip(progress, 0.0, 1.0))
    Fmin     = float(Fmin_min + (Fmin_max - Fmin_min) * progress)
    return float(min(Fmin, Fmin_max))

# -------------------------------------------------------------------
# soft_pool: spectral coverage diagnostic
# -------------------------------------------------------------------

def _soft_pool_select(eigvecs, d_loc, n_pool, coverage_floor):
    """
    Select softest modes sufficient to cover coverage_floor of ||defvec||^2.

    Modes are taken in index order (soft first). Returns the smallest prefix
    of modes whose cumulative normalised overlap^2 >= coverage_floor.

    Returns
    -------
    sel_idx    : int array of selected column indices into eigvecs
    coverage   : float, actual coverage achieved
    sufficient : bool, True if coverage_floor met within n_pool
    """
    norm_d = norm(d_loc)
    if norm_d < 1e-12:
        return np.array([], dtype=int), 0.0, False

    ov         = np.dot(d_loc, eigvecs[:, :n_pool])     # (n_pool,)
    cumcov     = np.cumsum((ov / norm_d) ** 2)

    meets = np.where(cumcov >= coverage_floor)[0]
    if meets.size > 0:
        k          = int(meets[0]) + 1
        sel_idx    = np.arange(k, dtype=int)
        coverage   = float(cumcov[k - 1])
        sufficient = True
    else:
        sel_idx    = np.arange(n_pool, dtype=int)
        coverage   = float(cumcov[-1]) if cumcov.size > 0 else 0.0
        sufficient = False

    return sel_idx, coverage, sufficient

# -------------------------------------------------------------------
# first-step diagnostic capture + auto-calibration
# -------------------------------------------------------------------

def _capture_diagnostic(eigvecs, d_loc, n_modes_avail, kwargs):
    """
    Called from calcStep on the first step using already-computed eigvecs.
    Computes soft-first cumulative overlap spectrum, logs diagnostics,
    and auto-calibrates parameters not already set by the user.
    Sets kwargs['_diag_done'] = True so it only runs once.
    """
    norm_d = norm(d_loc)
    if norm_d < 1e-8:
        kwargs['_diag_done'] = True
        return

    # Cumulative overlap spectrum in soft-first (index) order
    ov     = np.dot(d_loc, eigvecs[:, :n_modes_avail])
    ov_sq  = (ov / norm_d) ** 2
    cumcov = np.cumsum(ov_sq)

    # Coverage thresholds
    thresholds   = [0.25, 0.50, 0.75, 0.90]
    modes_needed = {}
    for thresh in thresholds:
        hits = np.where(cumcov >= thresh)[0]
        modes_needed[thresh] = int(hits[0]) + 1 if hits.size > 0 else n_modes_avail

    n_soft5        = max(1, int(np.floor(0.05 * n_modes_avail)))
    n_soft20       = max(1, int(np.floor(0.20 * n_modes_avail)))
    soft5_coverage = float(cumcov[n_soft5 - 1])
    soft20_coverage= float(cumcov[n_soft20 - 1])
    peak_mode      = int(np.argmax(ov_sq)) + 1
    hard_dominated = soft5_coverage < 0.10

    LOGGER.info("=" * 60)
    LOGGER.info("Transition diagnostic (from first-step ENM, %d modes)" % n_modes_avail)
    LOGGER.info("  Peak mode         = %d (coverage=%.3f)" % (peak_mode, float(ov_sq[peak_mode-1])))
    LOGGER.info("  Soft 5%%  coverage = %.3f (%d modes)" % (soft5_coverage, n_soft5))
    LOGGER.info("  Soft 20%% coverage = %.3f (%d modes)" % (soft20_coverage, n_soft20))
    for thresh in thresholds:
        LOGGER.info("  Modes for %.0f%% coverage = %d" % (thresh * 100, modes_needed[thresh]))
    if hard_dominated:
        LOGGER.warning(
            "Transition appears hard-mode dominated (soft 5%% coverage=%.3f < 0.10). "
            "ANM may not describe this transition well." % soft5_coverage
        )
    LOGGER.info("=" * 60)

    auto_calibrate = bool(kwargs.get('auto_calibrate', True))
    if auto_calibrate:
        rmsd0 = float(kwargs.get('_rmsd0', 1.0))

        # coverage_floor: based on overlap spectrum character
        if 'coverage_floor' not in kwargs:
            if modes_needed[0.50] <= n_soft5:
                cf = 0.70   # very soft-mode dominated
            elif modes_needed[0.50] <= n_soft20:
                cf = 0.50   # moderately soft
            else:
                cf = 0.35   # hard-mode dominated — be conservative
            kwargs['coverage_floor'] = cf
            LOGGER.info("auto_calibrate: coverage_floor = %.2f" % cf)

        # convergence_zone: 25% of rmsd0
        if 'convergence_zone' not in kwargs:
            cz = max(0.75, 0.25 * rmsd0)
            kwargs['convergence_zone'] = cz
            LOGGER.info("auto_calibrate: convergence_zone = %.3f A" % cz)

        # Fmin_min: based on soft5_coverage
        # clip(soft5*2, 0.15, 0.40): low start gives schedule room to ramp
        if 'Fmin_min' not in kwargs:
            fmin_min = float(np.clip(soft5_coverage * 2.0, 0.15, 0.40))
            kwargs['Fmin_min'] = fmin_min
            LOGGER.info("auto_calibrate: Fmin_min = %.3f" % fmin_min)

        # Fmin_max: cap at actual achievable coverage from available modes
        if 'Fmin_max' not in kwargs:
            total_cov = float(cumcov[-1])
            fmin_max = float(np.clip(total_cov * 0.95, 0.5, 0.9))
            kwargs['Fmin_max'] = fmin_max
            LOGGER.info("auto_calibrate: Fmin_max = %.3f (total_coverage=%.3f)" %
                        (fmin_max, total_cov))

        # pool_exp: smaller for large transitions so pool opens faster
        if 'pool_exp' not in kwargs:
            pe = float(np.clip(10.0 / rmsd0, 0.5, 2.0))
            kwargs['pool_exp'] = pe
            LOGGER.info("auto_calibrate: pool_exp = %.2f" % pe)

        # f: conservative for large transitions
        if 'f' not in kwargs:
            f_auto = float(np.clip(1.0 / rmsd0, 0.02, 0.5))
            kwargs['f'] = f_auto
            LOGGER.info("auto_calibrate: f = %.3f" % f_auto)

        # max_rewinds: more rewinds for larger transitions
        if 'max_rewinds' not in kwargs:
            mr = max(3, int(rmsd0 / 3.0))
            kwargs['max_rewinds'] = mr
            LOGGER.info("auto_calibrate: max_rewinds = %d" % mr)

        # target_rmsd: fixed 1.0 A — structurally meaningful regardless of transition size
        if 'target_rmsd' not in kwargs:
            kwargs['target_rmsd'] = 1.0
            LOGGER.info("auto_calibrate: target_rmsd = 1.000 A")

    kwargs['_diag_done'] = True
    kwargs['_diag'] = dict(
        soft5_coverage   = soft5_coverage,
        soft20_coverage  = soft20_coverage,
        modes_needed     = modes_needed,
        peak_mode        = peak_mode,
        hard_dominated   = hard_dominated,
        n_modes_avail    = n_modes_avail,
    )


# -------------------------------------------------------------------
# core single step
# -------------------------------------------------------------------

def calcStep(initial, target, n_modes, ensemble, defvecs, rmsds, Fmin,
             mask=None, callback_func=None, runtime_state=None, **kwargs):
    """
    Single adaptive step (v14).

    Revised controller support:
    - L1 uses a conservative soft-step scale (f_soft).
    - If L1 fails or only yields non-meaningful progress, L2 rescue is tried
      immediately at the same geometry using a separate rescue scale
      (f_rescue * rescue_boost).
    - Only meaningful steps are accepted; tiny positive steps do not advance the
      outer cycle.
    - Near the convergence zone, soft-mode insufficiency no longer short-circuits
      before L2 gets one last chance at the same geometry.
    """

    if runtime_state is None:
        runtime_state = {}

    # ---- core parameters ----
    f_soft         = float(kwargs.get('f_soft', kwargs.get('f', DEFAULT_F)))
    f_rescue_factor = float(kwargs.get('f_rescue_factor', 1.0))
    f_rescue_base  = float(kwargs.get('f_rescue', max(f_soft, f_soft * f_rescue_factor)))
    rescue_boost   = float(kwargs.get('rescue_boost', 1.0))
    f_rescue       = float(max(f_rescue_base * rescue_boost, f_soft))
    Fmin_max       = float(kwargs.get('Fmin_max', 0.9))
    cutoff         = float(kwargs.get('cutoff', 15.0))
    adj_max        = float(kwargs.get('adj_max', DEFAULT_ADJ_MAX))
    max_backtracks = int(kwargs.get('max_backtracks', 5))
    min_f_eff      = float(kwargs.get('min_f_eff', 1e-3))
    min_rmsd_step  = float(kwargs.get('min_rmsd_step', DEFAULT_MIN_RMSD_STEP))
    adj_exempt     = kwargs.get('adj_exempt', None)
    clash_min      = _clash_min_of(kwargs)    # default ON (DEFAULT_CLASH_MIN); 0 = off
    clash_sep      = int(kwargs.get('clash_sep', 3))
    clash_exempt   = kwargs.get('clash_exempt', None)
    _resolve_defaults(kwargs)                                     # effective defaults (see DEFAULT_* above)
    geom_bond_hi   = kwargs['geom_bond_hi']                       # 0 = geometric correction off
    geom_clash     = kwargs['geom_clash']
    geom_iters     = int(kwargs.get('geom_iters', 200))
    bond_floor     = kwargs.get('bond_floor', None)              # native-aware chain floors (build_chain_floors)
    d13_floor      = kwargs.get('d13_floor', None)
    bond_ceiling   = kwargs.get('bond_ceiling', None)            # native-aware upper bond limit (build_bond_ceiling)
    exempt_floor   = kwargs.get('exempt_floor', None)            # native distances of steric-exempt pairs
    risk_topk      = int(kwargs.get('risk_topk', 8))
    model          = kwargs.get('model', 'anm')
    target_rmsd    = float(kwargs.get('target_rmsd', 1.0))
    superposed_rmsd = bool(kwargs.get('superposed_rmsd', False))   # measure RMSD to the target after superposition
    meaningful_rmsd_improve = float(
        kwargs.get('meaningful_rmsd_improve', max(0.05, 0.05 * target_rmsd))
    )

    # ---- soft_pool parameters ----
    use_soft_pool      = bool(kwargs.get('soft_pool', True))
    coverage_floor     = float(kwargs.get('coverage_floor', 0.5))
    soft_filter        = float(_opt(kwargs, 'soft_filter', DEFAULT_SOFT_FILTER))   # 0 = off: standard driver
    soft_ladder        = bool(kwargs.get('soft_ladder', False))          # widen the soft prefix before any L2 rescue
    soft_filter_adapt  = bool(_opt(kwargs, 'soft_filter_adapt', DEFAULT_SOFT_FILTER_ADAPT))  # relax on stalled steps
    mu_step = soft_filter                                                # filter strength used for this step
    convergence_zone   = float(kwargs.get('convergence_zone', 1.5 * target_rmsd))
    defvec_noise_floor = float(kwargs.get('defvec_noise_floor', 0.05))

    # ---- legacy parameters (soft_pool=False only) ----
    progressive_nmax_flag = (
        bool(kwargs.get('progressive_nmax', False)) and not use_soft_pool
    )
    nmax_frac0 = float(kwargs.get('nmax_frac0', 0.2))
    nmax_eta   = float(kwargs.get('nmax_eta', 1.5))

    f_hist           = kwargs.setdefault('f_hist', [])
    min_slack_hist   = kwargs.setdefault('min_slack_hist', [])
    worst_edges_hist = kwargs.setdefault('worst_edges_hist', [])

    weights = ensemble.getWeights()
    if weights is not None:
        weights = weights.flatten()

    coords_init = initial
    coords_tar  = target

    # Number of trivial (rigid-body) modes: six for the classical ANM and for
    # the constrained ANM at kappa=0, three when a curvature weight is applied
    # (the curvature term preserves translations but lifts the three rotations).
    _kappa = float(kwargs.get('kappa', 0.0))
    _model_name = str(model).lower() if isinstance(model, str) else ''
    _is_genanm = _model_name in GENANM_MODEL_NAMES
    n_trivial = 3 if (_is_genanm and _kappa != 0.0) else 6

    dof = coords_init.size - n_trivial
    if dof <= 0:
        LOGGER.warning("Not enough DOF; returning without update.")
        return n_modes, 0.0, [], None

    raw_n_max_modes = kwargs.get('n_max_modes', dof)
    if isinstance(raw_n_max_modes, float) and 0 < raw_n_max_modes <= 1:
        n_max_modes = max(1, int(raw_n_max_modes * dof))
    else:
        n_max_modes = min(dof, int(raw_n_max_modes))

    defvec = coords_tar - coords_init
    d      = defvec.flatten()
    if weights is not None:
        d *= weights.repeat(3)

    if kwargs.get('f_auto', False):
        defvecs.append(d)
        f_soft = computeAdaptiveF(
            defvecs,
            f_max=kwargs.get('f_max', 1.0),
            f_min=kwargs.get('f_min_sched', 0.2),
            gamma=kwargs.get('f_gamma', 1.0)
        )
        defvecs.pop()
        f_rescue_base = float(max(kwargs.get('f_rescue', f_soft), f_soft * f_rescue_factor))
        f_rescue = float(max(f_rescue_base * rescue_boost, f_soft))
        kwargs['f'] = f_soft
        kwargs['f_soft'] = f_soft
        kwargs['f_rescue'] = f_rescue_base

    d_loc  = defvec.reshape(-1).copy()
    if weights is not None:
        d_loc *= weights.repeat(3)
    norm_d = norm(d_loc)

    Fmin = float(min(Fmin, Fmin_max))

    enm_kwargs = dict(kwargs)
    for key in ('model', 'n_modes', 'trim', 'superposed_rmsd'):
        enm_kwargs.pop(key, None)

    def _try_step(sel_idx_, sel_vecs_, ov_, f_seed_, rel_stiff_=None, mu_=None):
        n_sel = int(sel_idx_.size)
        if n_sel == 0:
            return None, None, None, 0.0
        mu_eff_ = mu_step if mu_ is None else float(mu_)
        if mu_eff_ > 0.0 and rel_stiff_ is not None:
            # stiffness filter: mode j contributes ov_j / (1 + mu (lambda_j / lambda_1 - 1)), so the softest mode
            # passes in full and stiffer modes are damped smoothly (mu = soft_filter)
            w_ = 1.0 / (1.0 + mu_eff_ * np.maximum(np.asarray(rel_stiff_, float) - 1.0, 0.0))
            v_ = sel_vecs_.dot(w_ * ov_)
        else:
            v_ = sel_vecs_.dot(ov_)
        denom_ = float(np.dot(v_, v_))
        if denom_ == 0.0:
            LOGGER.warning("Degenerate step direction (v·v == 0); returning without update.")
            return None, None, None, 0.0
        s_base_ = float(np.dot(v_, d_loc) / denom_)
        v3_     = v_.reshape(-1, 3)
        best_imp, best_cand_, best_f_, best_k_ = 0.0, None, None, None
        for k_ in range(int(max_backtracks) + 1):
            f_eff_ = float(f_seed_) / (2 ** k_)
            if f_eff_ < min_f_eff:
                break
            cand_ = coords_init + f_eff_ * s_base_ * v3_
            if geom_bond_hi > 0.0 or geom_clash > 0.0 or bond_floor is not None or d13_floor is not None or exempt_floor is not None:
                cand_, corr_rms_, corr_left_ = geometric_correct(
                    cand_, bond_hi=geom_bond_hi, clash_d=geom_clash, clash_sep=clash_sep,
                    adj_exempt=adj_exempt, clash_exempt=clash_exempt, iters=geom_iters,
                    bond_floor=bond_floor, d13_floor=d13_floor, exempt_floor=exempt_floor)
                step_rms_ = float(np.sqrt(((f_eff_ * s_base_ * v3_) ** 2).sum(axis=1).mean()))
                LOGGER.info("geom: f_eff=%.4f mode step %.3f A rms, correction %.3f A rms, residual violation %.3f A"
                            % (f_eff_, step_rms_, corr_rms_, corr_left_))
            if checkDisconnection(cand_, cutoff, adj_max=adj_max,
                                  adj_exempt=adj_exempt):
                continue
            if bond_floor is not None or d13_floor is not None or bond_ceiling is not None or exempt_floor is not None:
                # test the endpoint and the midpoint of the step, as for the steric test
                bad_, kind_, amt_ = checkChainGeometry(cand_, bond_floor, d13_floor, bond_ceiling=bond_ceiling,
                                                       exempt_floor=exempt_floor)
                if not bad_:
                    bad_, kind_, amt_ = checkChainGeometry(0.5 * (coords_init + cand_), bond_floor, d13_floor,
                                                           bond_ceiling=bond_ceiling, exempt_floor=exempt_floor)
                if bad_:
                    LOGGER.info("chain: rejected step at f_eff=%.4f (%s, by %.2f A)" % (f_eff_, kind_, amt_))
                    continue
            if clash_min > 0.0:
                # test the endpoint and the midpoint of the step, so that two
                # residues passing through each other within one step are caught
                hit_, dmin_ = checkClash(cand_, clash_min, clash_sep, clash_exempt)
                if not hit_:
                    mid_ = 0.5 * (coords_init + cand_)
                    hit_, dmin_ = checkClash(mid_, clash_min, clash_sep, clash_exempt)
                if hit_:
                    LOGGER.info("steric: rejected step at f_eff=%.4f (non-bonded Ca-Ca %.2f A < %.2f A)"
                                % (f_eff_, dmin_, clash_min))
                    continue
            imp_ = rmsd_before_step - _rmsd_to(cand_, coords_tar, weights, superposed_rmsd)
            if imp_ <= 0.0:
                continue
            if imp_ > best_imp:
                best_imp, best_cand_, best_f_, best_k_ = imp_, cand_, f_eff_, k_
        return best_cand_, best_f_, best_k_, best_imp

    _soft_pool_log    = None
    _first_step_diag  = None
    rmsd_before_step  = rmsds[-1]

    best_cand = None
    best_f_eff = None
    best_k = None
    best_improvement = 0.0
    sel_idx = np.array([], dtype=int)
    sel_vecs = None
    ov = None
    n_used = 0

    if use_soft_pool:
        if norm_d < defvec_noise_floor:
            LOGGER.info(
                "soft_pool: ||defvec||=%.4g < noise_floor=%.4g; signalling convergence." %
                (norm_d, defvec_noise_floor)
            )
            return 0, 0.0, [], None

        rmsd0_ref     = float(rmsds[0]) if rmsds[0] > 1e-8 else 1.0
        rmsd_progress = float(np.clip(1.0 - rmsds[-1] / rmsd0_ref, 0.0, 1.0))
        pool_exp      = float(kwargs.get('pool_exp', 2.0))
        if len(rmsds) <= 1:
            pool_floor = max(1, int(np.floor(0.20 * n_max_modes)))
        else:
            pool_floor = max(1, int(np.floor(0.05 * n_max_modes)))
        n_pool = int(max(pool_floor, int(np.floor(n_max_modes * (rmsd_progress ** pool_exp)))))
        n_pool = min(n_pool, n_max_modes)

        _cache       = runtime_state.get('_enm_cache', None)
        _cache_valid = (
            _cache is not None and
            _cache['n_pool'] == n_pool and
            np.array_equal(_cache['coords'], coords_init)
        )
        if _cache_valid:
            eigvecs = _cache['eigvecs']; eigvals = _cache['eigvals']
        else:
            anm_h, _ = _calcENM_contact(
                coords_init, select=mask, mask=mask,
                model=model, trim='trim', n_modes=n_pool, **enm_kwargs
            )
            if mask is not None:
                anm_h.masked = False
            eigvecs = anm_h.getEigvecs()[:, :n_pool]; eigvals = anm_h.getEigvals()[:n_pool]
            runtime_state['_enm_cache'] = dict(
                coords=coords_init.copy(), n_pool=n_pool, eigvecs=eigvecs, eigvals=eigvals
            )

        if len(defvecs) == 0 and '_diag_done' not in kwargs:
            _first_step_diag = (eigvecs.copy(), d_loc.copy(), n_pool)

        sel_idx_l1, coverage, sufficient = _soft_pool_select(
            eigvecs, d_loc, n_pool, coverage_floor
        )
        soft_exhausted_near_target = (not sufficient and rmsds[-1] <= convergence_zone)

        if not sufficient:
            if soft_exhausted_near_target:
                LOGGER.info(
                    "soft_pool: coverage=%.3f < floor=%.3f within pool=%d and RMSD=%.3f <= convergence_zone=%.3f; "
                    "giving L2 one final rescue attempt before signalling exhaustion." %
                    (coverage, coverage_floor, n_pool, rmsds[-1], convergence_zone)
                )
            else:
                LOGGER.info(
                    "soft_pool: coverage=%.3f < floor=%.3f within pool=%d; "
                    "trying full soft pool, then L2 rescue if needed." %
                    (coverage, coverage_floor, n_pool)
                )

        ov_soft = np.dot(d_loc, eigvecs[:, sel_idx_l1])
        l1_cand, l1_f_eff, l1_k, l1_improvement = _try_step(
            sel_idx_l1, eigvecs[:, sel_idx_l1], ov_soft, f_soft, eigvals[sel_idx_l1] / eigvals[0]
        )

        if mu_step > 0.0 and soft_filter_adapt and not (l1_cand is not None and l1_improvement >= meaningful_rmsd_improve):
            # adaptive filter: on a stalled filtered step, retry the same soft modes with the filter halved
            # (mu, mu/2, ... down to mu/8, then off); stiff contributions enter only when the softer step cannot progress
            mu_try = mu_step
            while mu_try > 0.0:
                mu_try = mu_try / 2.0 if mu_try / 2.0 >= soft_filter / 8.0 else 0.0
                c2 = _try_step(sel_idx_l1, eigvecs[:, sel_idx_l1], ov_soft, f_soft,
                               eigvals[sel_idx_l1] / eigvals[0], mu_=mu_try)
                if c2[0] is not None and c2[3] >= meaningful_rmsd_improve:
                    LOGGER.info("soft_filter: relaxed to mu = %.3g for this step" % mu_try)
                    l1_cand, l1_f_eff, l1_k, l1_improvement = c2
                    break
        if soft_ladder and not (l1_cand is not None and l1_improvement >= meaningful_rmsd_improve):
            # soft ladder: before any overlap-ranked (stiffness-blind) rescue, widen the soft prefix in order of
            # softness, doubling the coverage target each time, and keep the first widening that makes real progress
            cf_ = coverage_floor
            while cf_ < 0.95 and sel_idx_l1.size < n_pool:
                cf_ = min(0.95, 2.0 * cf_)
                sel2, cov2, _ = _soft_pool_select(eigvecs, d_loc, n_pool, cf_)
                if sel2.size <= sel_idx_l1.size:
                    continue
                ov2 = np.dot(d_loc, eigvecs[:, sel2])
                c2 = _try_step(sel2, eigvecs[:, sel2], ov2, f_soft, eigvals[sel2] / eigvals[0])
                if c2[0] is not None and c2[3] >= meaningful_rmsd_improve:
                    LOGGER.info("soft_ladder: widened the soft prefix from %d to %d modes (coverage %.3f)"
                                % (sel_idx_l1.size, sel2.size, cov2))
                    sel_idx_l1, coverage, ov_soft = sel2, cov2, ov2
                    l1_cand, l1_f_eff, l1_k, l1_improvement = c2
                    break
                sel_idx_l1, coverage, ov_soft = sel2, cov2, ov2            # keep widening from here
        if l1_cand is not None and l1_improvement >= meaningful_rmsd_improve:
            best_cand, best_f_eff, best_k, best_improvement = (
                l1_cand, l1_f_eff, l1_k, l1_improvement
            )
            sel_idx  = sel_idx_l1
            sel_vecs = eigvecs[:, sel_idx_l1]
            ov       = ov_soft
            n_used   = int(sel_idx.size)
            _soft_pool_log = (n_pool, n_used, coverage, coverage_floor, 'soft-prefix')
            runtime_state['_last_strategy'] = 'soft-prefix'
            runtime_state['_last_channel'] = 'soft'
        else:
            if l1_cand is not None and l1_improvement > 0.0:
                LOGGER.info(
                    "soft_pool: L1 produced only marginal progress (ΔRMSD=%.4g < meaningful=%.4g); trying L2 rescue at same geometry." %
                    (l1_improvement, meaningful_rmsd_improve)
                )

            n_l2_modes   = dof
            _cache_l2    = runtime_state.get('_enm_cache_l2', None)
            _l2_valid    = (
                _cache_l2 is not None and
                _cache_l2['n_pool'] == n_l2_modes and
                np.array_equal(_cache_l2['coords'], coords_init)
            )
            if _l2_valid:
                eigvecs_full = _cache_l2['eigvecs']; eigvals_full = _cache_l2['eigvals']
            else:
                anm_l2, _ = _calcENM_contact(
                    coords_init, select=mask, mask=mask,
                    model=model, trim='trim', n_modes=n_l2_modes, **enm_kwargs
                )
                if mask is not None:
                    anm_l2.masked = False
                eigvecs_full = anm_l2.getEigvecs()[:, :n_l2_modes]; eigvals_full = anm_l2.getEigvals()[:n_l2_modes]
                runtime_state['_enm_cache_l2'] = dict(
                    coords=coords_init.copy(), n_pool=n_l2_modes,
                    eigvecs=eigvecs_full, eigvals=eigvals_full
                )

            norm_d_safe = norm_d or 1.0
            ov_all      = np.dot(d_loc, eigvecs_full)
            norm_ov     = ov_all / norm_d_safe
            sort_idx    = np.argsort(-np.abs(norm_ov))
            c_sq        = np.cumsum(norm_ov[sort_idx] ** 2)
            pick_mask   = c_sq <= Fmin
            if not np.any(pick_mask):
                pick_mask[0] = True
            sel_idx_ov = sort_idx[pick_mask]
            ov_ov      = ov_all[sel_idx_ov]

            l2_cand, l2_f_eff, l2_k, l2_improvement = _try_step(
                sel_idx_ov, eigvecs_full[:, sel_idx_ov], ov_ov, f_rescue, eigvals_full[sel_idx_ov] / eigvals_full[0]
            )

            if l2_cand is not None and l2_improvement >= meaningful_rmsd_improve:
                best_cand, best_f_eff, best_k, best_improvement = (
                    l2_cand, l2_f_eff, l2_k, l2_improvement
                )
                sel_idx  = sel_idx_ov
                sel_vecs = eigvecs_full[:, sel_idx_ov]
                ov       = ov_ov
                n_used   = int(sel_idx.size)
                _soft_pool_log = (n_l2_modes, n_used, coverage, coverage_floor,
                                  'overlap-ranked-full')
                runtime_state['_last_strategy'] = 'overlap-ranked-full'
                runtime_state['_last_channel'] = 'rescue'
            else:
                if l2_cand is not None and l2_improvement > 0.0:
                    LOGGER.info(
                        "soft_pool: L2 produced only marginal progress (ΔRMSD=%.4g < meaningful=%.4g); treating this geometry as stalled." %
                        (l2_improvement, meaningful_rmsd_improve)
                    )
                if soft_exhausted_near_target:
                    LOGGER.info(
                        "soft_pool: both L1 and L2 are exhausted near target (RMSD=%.3f, target=%.3f); signalling completion/exhaustion." %
                        (rmsds[-1], target_rmsd)
                    )
                    return 0, 0.0, [], None
                LOGGER.warning('All step attempts failed or were non-meaningful; exiting step with no update.')
                return n_modes, 0.0, [], None

    else:
        if progressive_nmax_flag:
            eff_n_max_modes = computeProgressiveNmax(
                defvecs, n_max_modes, frac0=nmax_frac0, eta=nmax_eta
            )
        else:
            eff_n_max_modes = n_max_modes
        eff_n_max_modes = int(max(1, min(n_max_modes, eff_n_max_modes)))
        n_modes = min(max(1, int(n_modes)), eff_n_max_modes)

        _cache       = runtime_state.get('_enm_cache', None)
        _cache_valid = (
            _cache is not None and
            _cache['n_pool'] == eff_n_max_modes and
            np.array_equal(_cache['coords'], coords_init)
        )
        if _cache_valid:
            eigvecs = _cache['eigvecs']; eigvals = _cache['eigvals']
        else:
            anm_h, _ = _calcENM_contact(
                coords_init, select=mask, mask=mask,
                model=model, trim='trim', n_modes=eff_n_max_modes, **enm_kwargs
            )
            if mask is not None:
                anm_h.masked = False
            eigvecs = anm_h.getEigvecs()[:, :eff_n_max_modes]; eigvals = anm_h.getEigvals()[:eff_n_max_modes]
            runtime_state['_enm_cache'] = dict(
                coords=coords_init.copy(), n_pool=eff_n_max_modes, eigvecs=eigvecs, eigvals=eigvals
            )

        if len(defvecs) == 0 and '_diag_done' not in kwargs:
            _first_step_diag = (eigvecs.copy(), d_loc.copy(), eff_n_max_modes)

        norm_d_safe = norm_d or 1.0
        ov_all      = np.dot(d_loc, eigvecs)
        norm_ov     = ov_all / norm_d_safe
        sort_idx    = np.argsort(-np.abs(norm_ov))
        c_sq        = np.cumsum(norm_ov[sort_idx] ** 2)
        pick_mask   = c_sq <= Fmin
        if not np.any(pick_mask):
            pick_mask[0] = True

        sel_idx  = sort_idx[pick_mask]
        ov       = ov_all[sel_idx]
        sel_vecs = eigvecs[:, sel_idx]
        n_used   = int(sel_idx.size)

        LOGGER.info(
            "original: eff_n_max=%d, n_used=%d, nmax_frac=%.3f" %
            (eff_n_max_modes, n_used,
             float(eff_n_max_modes) / float(n_max_modes) if n_max_modes > 0 else 1.0)
        )

        runtime_state['_last_strategy'] = 'original'
        runtime_state['_last_channel'] = 'original'

        best_cand, best_f_eff, best_k, best_improvement = _try_step(
            sel_idx, sel_vecs, ov, f_soft, (eigvals[sel_idx] / eigvals[0]) if sel_idx.size else None
        )

        if best_cand is None or best_improvement < meaningful_rmsd_improve:
            if best_cand is not None and best_improvement > 0.0:
                LOGGER.info(
                    "original: candidate ΔRMSD=%.4g < meaningful=%.4g; treating geometry as stalled." %
                    (best_improvement, meaningful_rmsd_improve)
                )
            LOGGER.warning('All step attempts failed or were non-meaningful; exiting step with no update.')
            return n_modes, 0.0, [], None

    if n_used == 0:
        LOGGER.warning("No modes selected (n_used == 0); returning without update.")
        return n_modes, 0.0, [], None

    coords_updated = best_cand
    accepted_f_eff = best_f_eff
    used_h         = best_k
    runtime_state['_last_h'] = int(used_h)

    mid = 0.5 * (coords_init + coords_updated)
    ensemble.addCoordset(mid.copy())
    ensemble.addCoordset(coords_updated.copy())

    f_hist.append(float(accepted_f_eff))
    f_hist.append(float(accepted_f_eff))

    s_mid, ms_mid, worst_mid = compute_slack_info(
        mid, adj_max=adj_max, risk_topk=risk_topk, adj_exempt=adj_exempt
    )
    min_slack_hist.append(ms_mid)
    worst_edges_hist.append(worst_mid.tolist())

    s_end, ms_end, worst_end = compute_slack_info(
        coords_updated, adj_max=adj_max, risk_topk=risk_topk, adj_exempt=adj_exempt
    )
    min_slack_hist.append(ms_end)
    worst_edges_hist.append(worst_end.tolist())

    initial[:]       = coords_updated
    rmsd_after_step  = _rmsd_to(initial, coords_tar, weights, superposed_rmsd)
    rmsds.append(rmsd_after_step)
    rmsd_improvement = rmsd_before_step - rmsd_after_step

    defvecs.append(d)

    h         = int(used_h)
    scale_str = f"1/2^{h}"

    if use_soft_pool:
        eff_n_used  = int(sel_idx.size)
        nmax_frac   = float(eff_n_used) / float(n_max_modes) if n_max_modes > 0 else 1.0
    else:
        nmax_frac   = float(eff_n_max_modes) / float(n_max_modes) if n_max_modes > 0 else 1.0

    highest_mode = int(sel_idx.max()) + 1 if sel_idx.size > 0 else 0

    if use_soft_pool and _soft_pool_log is not None:
        pool_, k_, cov_, floor_, strategy_ = _soft_pool_log
        if strategy_ == 'overlap-ranked-full':
            LOGGER.info(
                "soft_pool: L1 (soft-prefix) was insufficient; L2 (overlap-ranked full-basis) provided the accepted move."
            )
        LOGGER.info(
            "soft_pool: pool=%d, k_selected=%d, coverage=%.3f "
            "(floor=%.3f, strategy=%s)" %
            (pool_, k_, cov_, floor_, strategy_)
        )

    step_line = (
        "Step successful [h=%d, scale=%s] "
        "(Fmin=%.3f, f_eff=%.3f, highest_mode=%d, ΔRMSD=%.4g)" %
        (used_h, scale_str, Fmin, accepted_f_eff, highest_mode, rmsd_improvement)
    )
    summary_line = (
        "Step summary (Fmin=%.3f, n_used=%d, nmax_frac=%.3f)" %
        (Fmin, n_used, nmax_frac)
    )
    rmsd_line = 'Current RMSD is %.6f\n' % rmsd_after_step

    LOGGER.info(step_line)
    LOGGER.info(summary_line)
    LOGGER.info(rmsd_line)

    step_log_lines = [step_line, summary_line, rmsd_line]

    if checkConvergence(rmsds, initial, **kwargs):
        n_modes_out = 0
    else:
        n_modes_out = max(1, min(int(n_modes), n_used))

    return n_modes_out, rmsd_improvement, step_log_lines, _first_step_diag

# -------------------------------------------------------------------
# trajectory history: snapshot management for rewind
# -------------------------------------------------------------------

class TrajectoryHistory:
    """
    Maintains a rolling window of accepted trajectory snapshots
    to support trajectory-level backtracking.

    Each snapshot stores everything needed to fully restore driver
    state: coords, rmsds, defvecs, ensemble coordset count, hists,
    rmsd_improvement, cycle_idx, f_global, and step strategy
    ('soft-prefix' or 'overlap-ranked').

    Rewind target selection:
    - Uses a single mixed candidate pool. soft-prefix snapshots are kept
      slightly preferred as a tie-breaker, but there is no hard first-pass
      filter anymore.
    - L2 / overlap-ranked snapshots remain eligible for rewind selection,
      but are not explicitly promoted just because they were large rescue
      steps.
    - Tracks already-tried targets so successive rewinds go progressively
      further back in history.
    """

    def __init__(self, max_depth=5):
        self.max_depth          = int(max_depth)
        self._history           = []
        self._tried_ids         = set()
        self._uid_counter       = 0
        self._last_snapshot_added = False  # flag for post-step update

    def snapshot(self, coords, rmsds, defvecs, n_coordsets,
                 f_hist, min_slack_hist, worst_edges_hist,
                 rmsd_improvement, cycle_idx, f_global, strategy='soft-prefix'):
        rmsd_val = float(rmsds[-1])

        # Deduplicate: same RMSD within 0.01A = same geometry, don't add again
        for existing in self._history:
            if abs(existing['rmsd_val'] - rmsd_val) < 0.01:
                if (rmsd_improvement > existing['rmsd_improvement'] and
                        existing['uid'] not in self._tried_ids):
                    existing['rmsd_improvement'] = float(rmsd_improvement)
                self._last_snapshot_added = False
                return

        uid = self._uid_counter
        self._uid_counter += 1
        self._last_snapshot_added = True
        entry = dict(
            coords           = coords.copy(),
            rmsds            = list(rmsds),
            defvecs          = list(defvecs),
            n_coordsets      = int(n_coordsets),
            f_hist           = list(f_hist),
            min_slack_hist   = list(min_slack_hist),
            worst_edges_hist = list(worst_edges_hist),
            rmsd_improvement = float(rmsd_improvement),
            cycle_idx        = int(cycle_idx),
            f_global         = float(f_global),
            strategy         = str(strategy),
            uid              = uid,
            rmsd_val         = rmsd_val,
            h_used           = 0,   # updated after step
            steps_after      = 0,   # how many successful steps followed this snapshot
        )
        # Increment steps_after for all existing snapshots — they all have
        # one more followable step than before
        for existing in self._history:
            existing['steps_after'] += 1
        self._history.append(entry)
        # Soft cap: drop highest-RMSD untried entry when over limit
        if len(self._history) > self.max_depth:
            worst = max(
                (e for e in self._history if e['uid'] not in self._tried_ids),
                key=lambda e: e['rmsd_val'], default=None
            )
            if worst is not None and worst['rmsd_val'] > rmsd_val:
                # remove by identity: entries hold NumPy arrays, so equality-based
                # list.remove() is ambiguous and raises
                self._history = [e for e in self._history if e is not worst]
            else:
                self._history.pop(0)

    def find_rewind_target(self, healthy_threshold, rewind_f_factor=0.5):
        """
        Pick the best untried snapshot from a single mixed pool.

        soft-prefix snapshots are *not* given a hard first-pass priority anymore.
        overlap-ranked / L2 snapshots are also eligible, but they are not
        explicitly promoted either. A very small strategy bias keeps soft-prefix
        snapshots slightly preferred when all else is nearly equal, while still
        allowing strong L2 snapshots to win on merit.
        """
        LOGGER.info("rewind_search: history=%d, tried=%d, threshold=%.4f" % (
            len(self._history), len(self._tried_ids), healthy_threshold))
        for e in sorted(self._history, key=lambda x: x['rmsd_val']):
            tried = e['uid'] in self._tried_ids
            LOGGER.info("  uid=%d rmsd=%.3f impr=%.4f h=%d strat=%s tried=%s" % (
                e['uid'], e['rmsd_val'], e['rmsd_improvement'],
                e.get('h_used', 0), e['strategy'], tried))

        def _score(e, rewind_f_factor):
            """
            Score a snapshot by expected post-rewind progress.

            expected_progress = improvement * rewind_factor / 2^h_used

            Divided by rmsd_val to prefer closer geometries when yield is similar.
            Multiplied by log(1 + steps_after) as a continuability bonus —
            a geometry followed by many successful steps is more likely to be
            followable again after a rewind with halved f.

            Strategy handling:
            - no hard soft-prefix pass
            - no explicit promotion of large L2 rescue snapshots
            - only a mild soft-prefix preference as a tie-breaker
            """
            h = e.get('h_used', 0)
            expected = e['rmsd_improvement'] * rewind_f_factor / max(2 ** h, 1)
            continuability = np.log1p(e.get('steps_after', 0) + 1)
            base = expected * continuability / max(e['rmsd_val'], 1.0)

            strategy = e.get('strategy', 'unknown')
            if strategy == 'soft-prefix':
                strategy_weight = 1.00
            elif strategy in ('overlap-ranked-full', 'overlap-ranked'):
                strategy_weight = 0.97
            else:
                strategy_weight = 0.99

            return -(base * strategy_weight)

        candidates = [
            e for e in self._history
            if (e['rmsd_improvement'] >= healthy_threshold and
                e['uid'] not in self._tried_ids and
                e['strategy'] in ('soft-prefix', 'overlap-ranked-full',
                                  'overlap-ranked', 'original', 'unknown'))
        ]
        if candidates:
            return min(candidates, key=lambda e: _score(e, rewind_f_factor))

        return None

    def mark_tried(self, cycle_idx_or_uid, rmsd_tolerance=0.05):
        """Mark uid as tried, plus any other snapshots within rmsd_tolerance."""
        uid = int(cycle_idx_or_uid)
        self._tried_ids.add(uid)
        # Find the RMSD of this target
        target_rmsd = None
        for e in self._history:
            if e['uid'] == uid:
                target_rmsd = e['rmsd_val']
                break
        # Mark nearby snapshots too — prevents revisiting same geometry region
        if target_rmsd is not None:
            for e in self._history:
                if abs(e['rmsd_val'] - target_rmsd) <= rmsd_tolerance:
                    self._tried_ids.add(e['uid'])

    def clear(self):
        self._history.clear()
        self._tried_ids.clear()
        self._uid_counter = 0
        self._last_snapshot_added = False

    def keep_before(self, target_cycle_idx):
        """No-op — all snapshots are preserved across rewinds."""
        pass

    def seed_with(self, entry):
        """
        Legacy: clear history and seed with a single entry.
        Prefer keep_before for progressive rewind.
        """
        self._history.clear()
        self._history.append(entry)

    def __len__(self):
        return len(self._history)


# -------------------------------------------------------------------
# temporary verbose logger for trajectory backtracking events
# -------------------------------------------------------------------

class _TempLog:
    """
    Writes verbose backtracking debug lines to a separate .rewind.log
    file alongside the main log if log_path is provided, otherwise
    falls back to LOGGER.info so rewind events are always visible.
    """
    def __init__(self, base_log_path):
        if base_log_path:
            stem = base_log_path.rsplit('.log', 1)[0]
            self.path = stem + '.rewind.log'
            self._fh  = open(self.path, 'w')
            self._fh.write(
                "# Trajectory rewind log — %s\n" %
                datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            )
            self._fh.flush()
        else:
            self.path = None
            self._fh  = None

    def write(self, msg):
        if self._fh is not None:
            self._fh.write(msg + '\n')
            self._fh.flush()
        else:
            # Fallback: emit via LOGGER so rewind events are never silent
            LOGGER.info('[REWIND] ' + msg)

    def close(self):
        if self._fh is not None:
            self._fh.close()


# -------------------------------------------------------------------
# final clean log writer
# -------------------------------------------------------------------

def _write_final_log(log_path, title, rmsd0, step_log_buffer):
    """
    Write a clean sequential final log directly to log_path.
    The raw verbose output (with rewinds) goes to basename.temp.log,
    arranged by the calling script. This writes only surviving steps
    with sequential cycle numbers — same format as original log so
    parse_logs.py works unchanged on the final .log file.
    """
    if not log_path:
        return

    with open(log_path, 'w') as fh:
        fh.write('@> Initialized Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd0))
        fh.write('\n')
        for seq_cycle, (cs, orig_cycle, lines) in enumerate(step_log_buffer):
            if seq_cycle > 0:
                fh.write('@> \n')
                fh.write('@> Starting cycle %d with initial structure %s\n' %
                         (seq_cycle, title))
            for line in lines:
                fh.write('@> %s\n' % line)


# -------------------------------------------------------------------
# dispatcher
# -------------------------------------------------------------------

def calcNeoAdaptiveANM(a, b, n_steps=100, mode=AANM_ONEWAY, **kwargs):
    """Generate a transition path from *a* to *b* with NeoAdaptiveANM and return
    it as an :class:`.Ensemble`, one coordinate set per accepted step.

    With no keyword arguments the published configuration is used, except the
    network: pass ``model='genanm'`` for the generalized ANM, which is
    recommended for pathways.

    :arg a: the structure to start from
    :type a: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg b: the structure to target, with atoms matching *a*
    :type b: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg n_steps: maximum number of steps (per direction for *AANM_BOTHWAYS*),
        default 100
    :type n_steps: int

    :arg mode: *AANM_ONEWAY* (default, from *a* to *b*), *AANM_ALTERNATING*, or
        *AANM_BOTHWAYS*
    :type mode: int

    **Network**

    :kwarg model: ``'anm'`` (default) or ``'genanm'``
    :kwarg k_bond: genANM bond stiffness, default 20.0 (also *k_theta* 10.0,
        *k_phi* 1.0, *kappa* 0.0; keep *kappa* at 0 for pathways)
    :kwarg cutoff: contact cutoff in A, default 15.0
    :kwarg contact_on: add stiff springs between non-bonded residues closer
        than this, in A, default 4.5 (0 turns them off)
    :kwarg contact_k: stiffness of those springs relative to the ANM spring
        constant, default 10.0
    :kwarg contact_sep: minimum sequence separation for a contact spring,
        default 3

    **Soft modes first**

    :kwarg soft_filter: stiffness filter strength mu, default 3.0. Each mode is
        scaled by 1 / (1 + mu (s - 1)), where s is its stiffness relative to the
        softest mode. 0 gives the standard (unfiltered) driver.
    :kwarg soft_filter_adapt: on a step that cannot progress, halve mu for that
        step (down to mu/8, then off), default **True**
    :kwarg soft_filter_fallback: if the run would stop short of the target,
        weaken the filter instead of stopping, default **False** (one-way only)
    :kwarg soft_ladder: on a step that cannot progress, widen the set of softest
        modes (coverage doubled each time) before the rescue that ranks modes by
        overlap alone, default **False**
    :kwarg coverage_floor: share of the remaining deformation the selected soft
        modes must cover, set automatically from the overlap spectrum by default

    **Physical steps**

    :kwarg adj_max: bond guard, the longest allowed sequential Calpha-Calpha
        distance in A, default 6.0; bonds already longer in the native structure
        (chain breaks) are exempt
    :kwarg clash_min: steric test, reject steps bringing non-bonded residues
        closer than this, in A, at the step end or midpoint, default 3.0
        (0 turns it off); pairs already closer in a native structure are exempt,
        but may not come closer than their native distance (by more than 0.05 A)
    :kwarg clash_sep: minimum sequence separation for the steric test, default 3
    :kwarg geom_bond_hi: geometric correction pulls bonds longer than this back,
        in A, default 4.2, and steps that still leave a bond above it (or above its native length, where longer)
        by more than 0.05 A are rejected (0 turns both off)
    :kwarg geom_clash: geometric correction pushes non-bonded pairs closer than
        this apart, in A, default 3.2 (0 turns it off)
    :kwarg geom_iters: maximum correction iterations per step, default 200
    :kwarg geom_bond_lo: shortest allowed virtual bond, in A, default 3.4; bonds already shorter in a native
        structure (for example cis peptides) keep their native length as the floor (0 turns it off)
    :kwarg geom_d13: shortest allowed distance between residues i and i+2, in A, default 4.4, native-aware in the
        same way (0 turns it off). Both limits are applied by the geometric correction and tested at the end and
        the midpoint of every step.

    **Step size, convergence and recovery**

    :kwarg f: step factor, default 0.2
    :kwarg target_rmsd: RMSD at which the path counts as converged, default 1.0
    :kwarg superposed_rmsd: measure every RMSD to the target after superposition
        (shape only), default **False**. ANM modes carry no rigid-body motion,
        so a rotation that builds up along the path is never removed; without
        superposition it keeps the RMSD above *target_rmsd* although the shapes
        have met. The RMSD drives the stop tests, the mode budget (Fmin), the
        acceptance of each step and the rewinds, so the option changes the
        path; the coordinates themselves are never superposed.
    :kwarg min_rmsd_diff: stop when the last RMSD changes are all smaller than
        this, default 0.0 (off)
    :kwarg min_rmsd_improve: stop when a step improves the RMSD by less than
        this, in A, default 0.005
    :kwarg min_rmsd_step: smallest accepted RMSD step, in A, default 0.0
    :kwarg max_backtracks: how many times a failing step is shortened, default 5
    :kwarg max_rewinds: how many times the path may return to an earlier
        accepted state, default 3

    **Output**

    :kwarg log_path: write a per-step log (modes used, step size, corrections,
        rejections) to this file, default **None**

    Every accepted step keeps virtual bonds within *adj_max* and non-bonded
    residues apart by at least *clash_min*. For large rearrangements that the
    soft-first driver cannot complete, run again with ``soft_filter=0`` and keep
    the path that ends closer to the target.

    Each geometric correction is logged as a share of the mode step. Report this
    share with the path; if it exceeds about 30%, the constraints rather than the
    modes drive the motion.
    """

    if not isinstance(n_steps, Integral):
        raise TypeError('n_steps must be an integer')

    if mode == AANM_ONEWAY:
        return calcOneWayAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == AANM_ALTERNATING:
        return calcAlternatingAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == AANM_BOTHWAYS:
        return calcBothWaysAdaptiveANM(a, b, n_steps, **kwargs)
    else:
        raise ValueError('unknown mode: %d' % mode)

# -------------------------------------------------------------------
# one-way driver
# -------------------------------------------------------------------

def calcOneWayAdaptiveANM(a, b, n_steps, **kwargs):
    """
    One-way adaptive ANM (v14, revised generic controller).

    Controller design:
    - L1 soft-prefix first with conservative f_soft.
    - If L1 is non-meaningful, try L2 rescue immediately at the same geometry,
      using a separate rescue scale (f_rescue_factor * rescue_boost).
    - If neither yields meaningful progress, perform a *short* same-geometry
      escalation: relaxed rescue, then max rescue.
    - Only after the short escalation is exhausted do we rewind.
    - After a strong L2 rescue, keep the trajectory in that basin briefly before
      allowing rewind again.
    """
    n_modes = kwargs.pop('n_modes', 20)
    n_modes_init = int(n_modes)

    coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd0 = checkInput(
        a, b, **kwargs
    )
    coordsA = coordsA.copy()

    LOGGER.timeit('_prody_calcAdaptiveANM')
    kwargs['_rmsd0'] = float(rmsd0)

    defvecs = []
    rmsds   = [rmsd0]

    adj_max = kwargs.get('adj_max', DEFAULT_ADJ_MAX)
    adj_tol = kwargs.get('adj_tol', 0.0)
    if 'adj_exempt' not in kwargs or kwargs['adj_exempt'] is None:
        kwargs['adj_exempt'] = build_adj_exempt(coordsA, adj_max=adj_max, tol=adj_tol)
    _resolve_defaults(kwargs)
    LOGGER.info("settings: steric test %s, geometric correction %s, contact-aware network %s" % (
        "on (%.1f A)" % kwargs['clash_min'] if kwargs['clash_min'] > 0 else "off",
        "on (bonds <= %.1f A, pairs >= %.1f A)" % (kwargs['geom_bond_hi'], kwargs['geom_clash'])
            if (kwargs['geom_bond_hi'] > 0 or kwargs['geom_clash'] > 0) else "off",
        "on (within %.1f A, k = %.0f)" % (kwargs['contact_on'], kwargs['contact_k']) if kwargs['contact_on'] > 0 else "off"))
    LOGGER.info("settings: stiffness filter %s" % (
        ("mu = %.3g, %s" % (kwargs['soft_filter'], "adaptive" if kwargs['soft_filter_adapt'] else "fixed"))
        if kwargs['soft_filter'] > 0 else "off (standard driver)"))
    if kwargs.get('bond_floor') is None and kwargs.get('d13_floor') is None:
        kwargs['bond_floor'], kwargs['d13_floor'] = build_chain_floors(
            [coordsA, coordsB], kwargs['geom_bond_lo'], kwargs['geom_d13'])
        LOGGER.info("chain limits: virtual bonds >= %s, residues i,i+2 >= %s (native-aware)" % (
            ("%.2f A" % kwargs['geom_bond_lo']) if kwargs['geom_bond_lo'] > 0 else "off",
            ("%.2f A" % kwargs['geom_d13']) if kwargs['geom_d13'] > 0 else "off"))
    if kwargs.get('bond_ceiling') is None:
        kwargs['bond_ceiling'] = build_bond_ceiling([coordsA, coordsB], kwargs['geom_bond_hi'])
    if kwargs['clash_min'] > 0.0 and kwargs.get('clash_exempt') is None:
        kwargs['clash_exempt'] = build_clash_exempt(
            [coordsA, coordsB], float(kwargs['clash_min']), int(kwargs.get('clash_sep', 3)))
        LOGGER.info("steric test on: clash_min=%.2f A, sequence separation >= %d, "
                    "%d native close pairs exempt" % (float(kwargs['clash_min']),
                    int(kwargs.get('clash_sep', 3)), len(kwargs['clash_exempt'])))
    if kwargs.get('clash_exempt') and kwargs.get('exempt_floor') is None:
        kwargs['exempt_floor'] = build_exempt_floors([coordsA, coordsB], kwargs['clash_exempt'])
    if float(kwargs.get('contact_on', 0.0) or 0.0) > 0.0 and kwargs.get('contact_exempt') is None:
        kwargs['contact_exempt'] = build_clash_exempt(
            [coordsA, coordsB], float(kwargs['contact_on']), int(kwargs.get('contact_sep', 3)))
        LOGGER.info("contact-aware network on: springs k=%.1f for non-bonded pairs within %.2f A "
                    "(sequence separation >= %d); %d native close pairs exempt" % (
                    float(kwargs.get('contact_k', 10.0)), float(kwargs['contact_on']),
                    int(kwargs.get('contact_sep', 3)), len(kwargs['contact_exempt'])))

    f0 = float(kwargs.get('f_soft', kwargs.get('f', DEFAULT_F)))
    kwargs.setdefault('f_hist', [float(f0)])
    kwargs.setdefault('min_slack_hist', [])
    kwargs.setdefault('worst_edges_hist', [])

    s0, ms0, worst0 = compute_slack_info(
        coordsA, adj_max=adj_max,
        risk_topk=kwargs.get('risk_topk', 8),
        adj_exempt=kwargs.get('adj_exempt', None)
    )
    kwargs['min_slack_hist'].append(ms0)
    kwargs['worst_edges_hist'].append(worst0.tolist())

    ensemble = Ensemble(title + '_aANM')
    ensemble.setAtoms(atoms)
    ensemble.setCoords(coordsB)
    ensemble.setWeights(weights)
    ensemble.addCoordset(coordsA.copy())

    Fmin_min            = kwargs.get('Fmin_min', 0.4)
    Fmin_max            = kwargs.get('Fmin_max', 0.9)
    min_rmsd_improve    = float(kwargs.get('min_rmsd_improve', DEFAULT_MIN_RMSD_IMPROVE))
    use_dynamic_fmin    = bool(kwargs.get('dynamic_fmin', True))
    stall_steps         = int(kwargs.get('stall_steps', 3))
    target_rmsd         = float(kwargs.get('target_rmsd', 1.0))
    base_meaningful     = float(kwargs.get('meaningful_rmsd_improve', max(0.05, 0.05 * target_rmsd)))
    meaningful_fraction = float(kwargs.get('meaningful_rmsd_fraction', 0.01))
    relaxed_fmin_gap    = float(kwargs.get('relaxed_fmin_gap', 0.05))

    f_rescue_factor     = float(kwargs.get('f_rescue_factor', 1.0))
    rescue_boost_relaxed= float(kwargs.get('rescue_boost_relaxed', 1.25))
    rescue_boost_max    = float(kwargs.get('rescue_boost_max', 1.50))
    strong_l2_abs       = float(kwargs.get('strong_l2_abs', 0.5))
    strong_l2_frac      = float(kwargs.get('strong_l2_frac', 0.05))
    rescue_commit_window = int(kwargs.get('rescue_commit_window', 2))
    rescue_commit_grace_rounds = int(kwargs.get('rescue_commit_grace_rounds', 1))

    max_rewinds          = int(kwargs.get('max_rewinds', 3))
    rewind_history_depth = int(kwargs.get('rewind_history_depth',
                               max(50, int(kwargs.get('max_rewinds', 3)) * 10)))
    rewind_f_factor      = float(kwargs.get('rewind_f_factor', 0.75))
    f_recovery_rate      = float(kwargs.get('f_recovery_rate', 2.0))
    base_log_path        = kwargs.get('log_path', None)

    traj_history = TrajectoryHistory(max_depth=rewind_history_depth)
    tmp_log      = _TempLog(base_log_path)
    rewind_count = 0
    f_current    = float(f0)
    rescue_boost = 1.0
    runtime_state = {}

    step_log_buffer = []
    stall_counter   = 0
    fmin_offset     = 0.0
    terminated_on_fmin_cap = False
    commit_cycles_remaining = 0
    commit_grace_remaining  = 0

    n = 0
    attempt_idx = 0
    max_attempts_per_phase = int(kwargs.get('max_attempts_per_phase', n_steps * 10))

    def _current_convergence_zone():
        return float(kwargs.get('convergence_zone', max(0.75, 0.25 * rmsd0)))

    def _meaningful_threshold(current_rmsd):
        zone = _current_convergence_zone()
        if current_rmsd > zone:
            return max(base_meaningful, meaningful_fraction * float(current_rmsd))
        return base_meaningful

    def _healthy_threshold(current_rmsd):
        current_meaningful = _meaningful_threshold(current_rmsd)
        return max(min_rmsd_improve, 0.5 * current_meaningful, float(current_rmsd) * 0.002)

    def _strong_l2_threshold(current_rmsd):
        return max(strong_l2_abs, strong_l2_frac * float(current_rmsd))

    def _trim_to_coordsets(n_keep):
        while ensemble.numCoordsets() > n_keep:
            ensemble.delCoordset(-1)
        step_log_buffer[:] = [
            (cs, ci, ll) for cs, ci, ll in step_log_buffer if cs <= n_keep
        ]

    def _restore_snapshot(target, reason_label):
        nonlocal rewind_count, f_current, fmin_offset, stall_counter
        nonlocal n, n_modes, attempt_idx, runtime_state
        nonlocal rescue_boost, commit_cycles_remaining, commit_grace_remaining

        rewind_count += 1
        traj_history.mark_tried(target['uid'])
        old_f = f_current
        target_f = float(target.get('f_global', f_current))
        new_f = max(1e-6, target_f * rewind_f_factor)

        tmp_log.write(
            "[REWIND %d/%d] %s Rewinding to cycle %d (RMSD=%.4f, improvement was %.4f). "
            "Reducing f_soft: %.4f -> %.4f. Resetting rescue state and Fmin escalation." % (
                rewind_count, max_rewinds, reason_label,
                target['cycle_idx'], target['rmsds'][-1], target['rmsd_improvement'],
                old_f, new_f,
            )
        )

        coordsA[:]                     = target['coords']
        rmsds[:]                       = target['rmsds']
        defvecs[:]                     = target['defvecs']
        kwargs['f_hist'][:]            = target['f_hist']
        kwargs['min_slack_hist'][:]    = target['min_slack_hist']
        kwargs['worst_edges_hist'][:]  = target['worst_edges_hist']
        _trim_to_coordsets(target['n_coordsets'])

        runtime_state.clear()
        f_current     = new_f
        kwargs['f']   = f_current
        kwargs['f_soft'] = f_current
        kwargs.pop('f_rescue', None)
        fmin_offset   = 0.0
        stall_counter = 0
        rescue_boost  = 1.0
        commit_cycles_remaining = 0
        commit_grace_remaining  = 0
        n            = int(target['cycle_idx'])
        n_modes      = n_modes_init
        attempt_idx  = 0
        traj_history.keep_before(target['cycle_idx'])

        tmp_log.write(
            "[REWIND %d/%d] State restored. Resuming from RMSD=%.4f with f_soft=%.4f." % (
                rewind_count, max_rewinds, rmsds[-1], f_current
            )
        )

    def _filter_fallback(reason):
        """Filter fallback: instead of ending a filtered run short of the target, weaken the filter (mu/2 ... mu/8,
        then off) and continue. mu lives in kwargs, so rewinds (which clear runtime_state) cannot reset it."""
        mu = float(kwargs.get('soft_filter', 0.0) or 0.0)
        if not kwargs.get('soft_filter_fallback') or mu <= 0.0 or rmsds[-1] <= target_rmsd:
            return False
        mu0 = float(kwargs.setdefault('_soft_filter0', mu))
        new = mu / 2.0 if mu / 2.0 >= mu0 / 8.0 else 0.0
        kwargs['soft_filter'] = new
        LOGGER.info("soft_filter: fallback instead of stopping (%s) at RMSD %.3f; mu %.3g -> %.3g"
                    % (reason, rmsds[-1], mu, new))
        return True

    while n < n_steps:
        if attempt_idx >= max_attempts_per_phase:
            LOGGER.warning(
                "Max attempts per phase (%d) reached without convergence." %
                max_attempts_per_phase
            )
            break

        rmsd_current = rmsds[-1]
        current_meaningful = _meaningful_threshold(rmsd_current)
        kwargs['meaningful_rmsd_improve'] = current_meaningful
        kwargs['f'] = f_current
        kwargs['f_soft'] = f_current
        kwargs.pop('f_rescue', None)
        kwargs['f_rescue_factor'] = f_rescue_factor
        kwargs['rescue_boost'] = rescue_boost

        _Fmin_min = float(kwargs.get('Fmin_min', Fmin_min))
        _Fmin_max = float(kwargs.get('Fmin_max', Fmin_max))
        base_Fmin = _compute_Fmin_from_RMSD(
            rmsd_current, rmsd0, Fmin_min=_Fmin_min, Fmin_max=_Fmin_max
        )
        Fmin_step = min(base_Fmin + fmin_offset, _Fmin_max)

        max_offset = max(0.0, _Fmin_max - _Fmin_min)
        relaxed_abs = max(_Fmin_min, min(_Fmin_max, _Fmin_max - relaxed_fmin_gap))
        relaxed_offset = min(max_offset, max(0.0, relaxed_abs - base_Fmin))
        max_rescue_offset = min(max_offset, max(0.0, _Fmin_max - base_Fmin))

        prestep_idx = None
        if stall_counter == 0 and fmin_offset == 0.0:
            traj_history.snapshot(
                coords           = coordsA.copy(),
                rmsds            = list(rmsds),
                defvecs          = list(defvecs),
                n_coordsets      = ensemble.numCoordsets(),
                f_hist           = list(kwargs['f_hist']),
                min_slack_hist   = list(kwargs['min_slack_hist']),
                worst_edges_hist = list(kwargs['worst_edges_hist']),
                rmsd_improvement = 0.0,
                cycle_idx        = n,
                f_global         = f_current,
                strategy         = 'unknown',
            )
            if traj_history._last_snapshot_added and traj_history._history:
                prestep_idx = len(traj_history._history) - 1

        n_modes_out, rmsd_improvement, _step_log_lines, _diag_info = calcStep(
            coordsA, coordsB, n_modes, ensemble, defvecs, rmsds,
            Fmin=Fmin_step, mask=maskA, runtime_state=runtime_state, **kwargs
        )
        attempt_idx += 1

        if _diag_info is not None and '_diag_done' not in kwargs:
            _capture_diagnostic(_diag_info[0], _diag_info[1], _diag_info[2], kwargs)
            Fmin_min    = float(kwargs.get('Fmin_min', Fmin_min))
            Fmin_max    = float(kwargs.get('Fmin_max', Fmin_max))
            max_rewinds = int(kwargs.get('max_rewinds', max_rewinds))
            target_rmsd = float(kwargs.get('target_rmsd', target_rmsd))
            base_meaningful = float(kwargs.get('meaningful_rmsd_improve', max(0.05, 0.05 * target_rmsd)))
            f0 = float(kwargs.get('f_soft', kwargs.get('f', f0)))
            if n == 0:
                f_current = f0
            new_depth = max(50, max_rewinds * 10)
            traj_history.max_depth = new_depth
            LOGGER.info("auto_calibrate: rewind_history_depth = %d" % new_depth)

        if rmsd_improvement > 0.0 and _step_log_lines:
            step_log_buffer.append((ensemble.numCoordsets(), n, _step_log_lines))

        if n_modes_out == 0:
            if rmsd_improvement > 0.0:
                if prestep_idx is not None and prestep_idx < len(traj_history._history):
                    traj_history._history[prestep_idx]['rmsd_improvement'] = float(rmsd_improvement)
                    traj_history._history[prestep_idx]['strategy'] = runtime_state.get('_last_strategy', 'unknown')
                    traj_history._history[prestep_idx]['h_used'] = int(runtime_state.get('_last_h', 0))
                LOGGER.report('One-way Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
                break

            if rmsds[-1] <= target_rmsd:
                LOGGER.report('One-way Adaptive ANM converged in %.2fs.', '_prody_calcAdaptiveANM')
                break

            healthy_thresh = _healthy_threshold(rmsds[-1])
            target = None
            if rewind_count < int(kwargs.get('max_rewinds', max_rewinds)) and commit_cycles_remaining <= 0:
                target = traj_history.find_rewind_target(
                    healthy_thresh,
                    rewind_f_factor=rewind_f_factor,
                )
            if target is not None:
                _restore_snapshot(
                    target,
                    "Soft modes and rescue exhausted above target."
                )
                continue

            if _filter_fallback("soft modes and rescue exhausted"):
                continue
            tmp_log.write(
                "[EXHAUSTED] L1/L2 exhausted near RMSD=%.4f but no healthy rewind target remained. Accepting termination." %
                rmsds[-1]
            )
            break

        if rmsd_improvement > 0.0:
            if prestep_idx is not None and prestep_idx < len(traj_history._history):
                traj_history._history[prestep_idx]['rmsd_improvement'] = float(rmsd_improvement)
                traj_history._history[prestep_idx]['strategy'] = runtime_state.get('_last_strategy', 'unknown')
                traj_history._history[prestep_idx]['h_used'] = int(runtime_state.get('_last_h', 0))

            last_strategy = runtime_state.get('_last_strategy', 'unknown')
            strong_l2 = (
                last_strategy == 'overlap-ranked-full' and
                rmsd_improvement >= _strong_l2_threshold(rmsd_current)
            )
            if strong_l2:
                commit_cycles_remaining = max(commit_cycles_remaining, rescue_commit_window)
                commit_grace_remaining = max(commit_grace_remaining, rescue_commit_grace_rounds)
                tmp_log.write(
                    "[COMMIT] Strong L2 rescue detected at RMSD=%.4f (ΔRMSD=%.4f). Protecting current basin for %d succeeding cycles with %d grace round(s)." %
                    (rmsd_current, rmsd_improvement, commit_cycles_remaining, commit_grace_remaining)
                )
            elif commit_cycles_remaining > 0:
                commit_cycles_remaining -= 1
                if commit_cycles_remaining == 0:
                    commit_grace_remaining = 0

            should_capture_postjump = (
                rmsd_improvement > current_meaningful * 5.0 and
                last_strategy != 'overlap-ranked-full'
            )
            if should_capture_postjump:
                traj_history.snapshot(
                    coords           = coordsA.copy(),
                    rmsds            = list(rmsds),
                    defvecs          = list(defvecs),
                    n_coordsets      = ensemble.numCoordsets(),
                    f_hist           = list(kwargs['f_hist']),
                    min_slack_hist   = list(kwargs['min_slack_hist']),
                    worst_edges_hist = list(kwargs['worst_edges_hist']),
                    rmsd_improvement = float(rmsd_improvement),
                    cycle_idx        = n,
                    f_global         = f_current,
                    strategy         = last_strategy,
                )
                if traj_history._last_snapshot_added and traj_history._history:
                    traj_history._history[-1]['h_used'] = int(runtime_state.get('_last_h', 0))

            stall_counter = 0
            fmin_offset   = 0.0
            rescue_boost  = 1.0
            n += 1
            n_modes = max(1, int(n_modes_out))

            if rewind_count > 0 and f_current < f0:
                old_f = f_current
                f_current = min(f_current * f_recovery_rate, f0)
                kwargs['f'] = f_current
                kwargs['f_soft'] = f_current
                if f_current > old_f:
                    tmp_log.write(
                        "[F-RECOVERY] f_soft recovered %.4f -> %.4f (target=%.4f)" %
                        (old_f, f_current, f0)
                    )

            if n < n_steps:
                LOGGER.info('')
                LOGGER.info(f"Starting cycle {n + 1} with initial structure {title}")
            continue

        if not use_dynamic_fmin:
            LOGGER.warning(
                "No meaningful progress (threshold=%.4g) with dynamic_fmin disabled; stopping." %
                current_meaningful
            )
            break

        stall_counter += 1
        LOGGER.info(
            "No meaningful progress at current geometry (threshold=%.4g). Stall counter: %d/%d" %
            (current_meaningful, stall_counter, stall_steps)
        )

        if stall_counter < stall_steps:
            continue

        stall_counter = 0
        old_offset = fmin_offset
        stage_label = None
        if old_offset + 1e-12 < relaxed_offset:
            fmin_offset = relaxed_offset
            rescue_boost = max(rescue_boost, rescue_boost_relaxed)
            stage_label = 'relaxed-rescue'
        elif old_offset + 1e-12 < max_rescue_offset:
            fmin_offset = max_rescue_offset
            rescue_boost = max(rescue_boost, rescue_boost_max)
            stage_label = 'max-rescue'

        if stage_label is not None and fmin_offset > old_offset + 1e-12:
            LOGGER.warning(
                "Same-geometry stall persisted; escalating to %s (Fmin offset=%.3f, rescue_boost=%.2f)." %
                (stage_label, fmin_offset, rescue_boost)
            )
            continue

        if commit_cycles_remaining > 0 and commit_grace_remaining > 0:
            commit_grace_remaining -= 1
            fmin_offset = max_rescue_offset
            rescue_boost = max(rescue_boost, rescue_boost_max)
            LOGGER.warning(
                "Post-rescue commitment active; granting one extra max-rescue round before rewind (remaining grace=%d)." %
                commit_grace_remaining
            )
            continue

        healthy_thresh = _healthy_threshold(rmsds[-1])
        target = None
        if rewind_count < int(kwargs.get('max_rewinds', max_rewinds)):
            target = traj_history.find_rewind_target(
                healthy_thresh,
                rewind_f_factor=rewind_f_factor,
            )
        if target is not None:
            _restore_snapshot(
                target,
                "Short same-geometry rescue ladder exhausted."
            )
            continue

        if _filter_fallback("rescue ladder exhausted"):
            stall_counter = 0
            continue
        LOGGER.info(
            "Short same-geometry rescue ladder exhausted at maximum Fmin offset (%.3f); ending simulation gracefully." %
            fmin_offset
        )
        terminated_on_fmin_cap = True
        break

    if terminated_on_fmin_cap:
        LOGGER.info(
            "One-way Adaptive ANM terminated after exhausting the short rescue ladder at Fmin offset %.3f." %
            fmin_offset
        )

    LOGGER.report('One-way Adaptive ANM finished in %.2fs.', '_prody_calcAdaptiveANM')
    tmp_log.write("[DONE] Simulation finished. Total rewinds: %d/%d." % (
        rewind_count, max_rewinds
    ))
    tmp_log.close()

    if base_log_path:
        _write_final_log(base_log_path, title, rmsd0, step_log_buffer)
    else:
        LOGGER.warning(
            "No log_path provided; final clean log not written. Pass log_path=<path> to generate a rewind-free sequential log."
        )

    return ensemble

# -------------------------------------------------------------------
# alternating driver
# -------------------------------------------------------------------

def calcAlternatingAdaptiveANM(a, b, n_steps, **kwargs):
    """
    Alternating A->B and B->A driver (v14).
    n_pool is shared across directions (same protein, same mode space).
    Fmin-bump stall counters are independent per direction.
    """
    n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd0 = checkInput(
        a, b, **kwargs
    )
    coordsA, coordsB = coordsA.copy(), coordsB.copy()

    LOGGER.timeit('_prody_calcAdaptiveANM')

    # Pre-simulation diagnostic and auto-calibration
    kwargs['_rmsd0'] = float(rmsd0)  # for auto-calibration in _capture_diagnostic

    n                  = 0
    defvecsA, defvecsB = [], []
    rmsds              = [rmsd0]

    adj_max = kwargs.get('adj_max', DEFAULT_ADJ_MAX)
    adj_tol = kwargs.get('adj_tol', 0.0)

    adj_exempt_common = kwargs.get('adj_exempt', None)
    adj_exempt_A = kwargs.get('adj_exempt_A') or adj_exempt_common or \
                   build_adj_exempt(coordsA, adj_max=adj_max, tol=adj_tol)
    adj_exempt_B = kwargs.get('adj_exempt_B') or adj_exempt_common or \
                   build_adj_exempt(coordsB, adj_max=adj_max, tol=adj_tol)
    _resolve_defaults(kwargs)
    LOGGER.info("settings: steric test %s, geometric correction %s, contact-aware network %s" % (
        "on (%.1f A)" % kwargs['clash_min'] if kwargs['clash_min'] > 0 else "off",
        "on (bonds <= %.1f A, pairs >= %.1f A)" % (kwargs['geom_bond_hi'], kwargs['geom_clash'])
            if (kwargs['geom_bond_hi'] > 0 or kwargs['geom_clash'] > 0) else "off",
        "on (within %.1f A, k = %.0f)" % (kwargs['contact_on'], kwargs['contact_k']) if kwargs['contact_on'] > 0 else "off"))
    LOGGER.info("settings: stiffness filter %s" % (
        ("mu = %.3g, %s" % (kwargs['soft_filter'], "adaptive" if kwargs['soft_filter_adapt'] else "fixed"))
        if kwargs['soft_filter'] > 0 else "off (standard driver)"))
    if kwargs.get('bond_floor') is None and kwargs.get('d13_floor') is None:
        kwargs['bond_floor'], kwargs['d13_floor'] = build_chain_floors(
            [coordsA, coordsB], kwargs['geom_bond_lo'], kwargs['geom_d13'])
        LOGGER.info("chain limits: virtual bonds >= %s, residues i,i+2 >= %s (native-aware)" % (
            ("%.2f A" % kwargs['geom_bond_lo']) if kwargs['geom_bond_lo'] > 0 else "off",
            ("%.2f A" % kwargs['geom_d13']) if kwargs['geom_d13'] > 0 else "off"))
    if kwargs.get('bond_ceiling') is None:
        kwargs['bond_ceiling'] = build_bond_ceiling([coordsA, coordsB], kwargs['geom_bond_hi'])
    if kwargs['clash_min'] > 0.0 and kwargs.get('clash_exempt') is None:
        kwargs['clash_exempt'] = build_clash_exempt(
            [coordsA, coordsB], float(kwargs['clash_min']), int(kwargs.get('clash_sep', 3)))
        LOGGER.info("steric test on: clash_min=%.2f A, sequence separation >= %d, "
                    "%d native close pairs exempt" % (float(kwargs['clash_min']),
                    int(kwargs.get('clash_sep', 3)), len(kwargs['clash_exempt'])))
    if kwargs.get('clash_exempt') and kwargs.get('exempt_floor') is None:
        kwargs['exempt_floor'] = build_exempt_floors([coordsA, coordsB], kwargs['clash_exempt'])
    if float(kwargs.get('contact_on', 0.0) or 0.0) > 0.0 and kwargs.get('contact_exempt') is None:
        kwargs['contact_exempt'] = build_clash_exempt(
            [coordsA, coordsB], float(kwargs['contact_on']), int(kwargs.get('contact_sep', 3)))
        LOGGER.info("contact-aware network on: springs k=%.1f for non-bonded pairs within %.2f A "
                    "(sequence separation >= %d); %d native close pairs exempt" % (
                    float(kwargs.get('contact_k', 10.0)), float(kwargs['contact_on']),
                    int(kwargs.get('contact_sep', 3)), len(kwargs['contact_exempt'])))

    f0 = kwargs.get('f', 1.0)
    kwargs.setdefault('f_histA', [float(f0)])
    kwargs.setdefault('f_histB', [float(f0)])
    kwargs.setdefault('min_slack_histA', [])
    kwargs.setdefault('min_slack_histB', [])
    kwargs.setdefault('worst_edges_histA', [])
    kwargs.setdefault('worst_edges_histB', [])

    for coords, hist_ms, hist_we, adj_ex in [
        (coordsA, 'min_slack_histA', 'worst_edges_histA', adj_exempt_A),
        (coordsB, 'min_slack_histB', 'worst_edges_histB', adj_exempt_B),
    ]:
        s, ms, worst = compute_slack_info(
            coords, adj_max=adj_max,
            risk_topk=kwargs.get('risk_topk', 8), adj_exempt=adj_ex
        )
        kwargs[hist_ms].append(ms)
        kwargs[hist_we].append(worst.tolist())

    ensA = Ensemble(getTitle(a, 'structureA') + '_aANM')
    ensA.setAtoms(atoms); ensA.setCoords(coordsB)
    ensA.setWeights(weights); ensA.addCoordset(coordsA.copy())

    ensB = Ensemble(getTitle(b, 'structureB') + '_aANM')
    ensB.setAtoms(atoms); ensB.setCoords(coordsA)
    ensB.setWeights(weights); ensB.addCoordset(coordsB.copy())

    Fmin_min            = kwargs.get('Fmin_min', 0.4)
    Fmin_max            = kwargs.get('Fmin_max', 0.9)
    min_rmsd_improve    = float(kwargs.get('min_rmsd_improve', DEFAULT_MIN_RMSD_IMPROVE))
    use_dynamic_fmin    = bool(kwargs.get('dynamic_fmin', True))
    stall_steps         = int(kwargs.get('stall_steps', 3))
    Fmin_bump_increment = float(kwargs.get('Fmin_bump_increment', 0.05))

    # n_pool is recomputed fresh each step from RMSD progress — no initialisation needed.

    stateA = {'stall_counter': 0, 'fmin_offset': 0.0}
    stateB = {'stall_counter': 0, 'fmin_offset': 0.0}
    runtime_stateA = {}
    runtime_stateB = {}
    terminated_on_fmin_cap = False

    base_log_path = kwargs.get('log_path', None)
    step_log_buffer = []  # (n_coordsets_after, cycle_idx, log_lines)

    while n < n_steps:

        # ---------- A -> B ----------
        rmsd_current = rmsds[-1]
        Fmin_stepA   = min(
            _compute_Fmin_from_RMSD(rmsd_current, rmsd0, Fmin_min, Fmin_max)
            + stateA['fmin_offset'], Fmin_max
        )

        LOGGER.info('')
        LOGGER.info(
            f"Starting cycle {n + 1} with initial structure "
            f"{getTitle(a, 'structure A')}"
        )

        n_modes_A, rmsd_improvement_A, _step_log_A, _diag_info_A = calcStep(
            coordsA, coordsB, n_modes, ensA, defvecsA, rmsds,
            Fmin=Fmin_stepA, mask=maskA,
            f_hist=kwargs['f_histA'],
            min_slack_hist=kwargs['min_slack_histA'],
            worst_edges_hist=kwargs['worst_edges_histA'],
            adj_exempt=adj_exempt_A,
            runtime_state=runtime_stateA,
            **kwargs
        )
        if _diag_info_A is not None and '_diag_done' not in kwargs:
            _capture_diagnostic(_diag_info_A[0], _diag_info_A[1], _diag_info_A[2], kwargs)
            Fmin_min = float(kwargs.get('Fmin_min', Fmin_min))
            Fmin_max = float(kwargs.get('Fmin_max', Fmin_max))
        if rmsd_improvement_A > 0.0 and _step_log_A:
            step_log_buffer.append((ensA.numCoordsets(), n, _step_log_A))

        if use_dynamic_fmin:
            if rmsd_improvement_A < min_rmsd_improve and n_modes_A > 0:
                stateA['stall_counter'] += 1
            else:
                stateA['stall_counter'] = 0
            if stateA['stall_counter'] >= stall_steps:
                max_offset = max(0.0, Fmin_max - Fmin_min)
                old        = stateA['fmin_offset']
                stateA['fmin_offset'] = min(old + Fmin_bump_increment, max_offset)
                if stateA['fmin_offset'] <= old:
                    LOGGER.info(
                        "A->B Fmin offset has reached its maximum (%.3f) after "
                        "%d small-improvement cycles; ending simulation "
                        "gracefully." % (stateA['fmin_offset'], stall_steps)
                    )
                    terminated_on_fmin_cap = True
                    break
                LOGGER.warning(
                    "A->B stalled. Increasing Fmin offset additively to "
                    f"{stateA['fmin_offset']:.3f}."
                )
                stateA['stall_counter'] = 0
        else:
            if rmsd_improvement_A < min_rmsd_improve and n_modes_A > 0:
                LOGGER.warning(
                    "A->B ΔRMSD=%.2e < min_rmsd_improve=%.2e. Stopping." %
                    (rmsd_improvement_A, min_rmsd_improve)
                )
                break

        if n_modes_A == 0:
            LOGGER.report('Alternating Adaptive ANM converged in %.2fs.',
                          '_prody_calcAdaptiveANM')
            break

        # ---------- B -> A ----------
        rmsd_current = rmsds[-1]
        Fmin_stepB   = min(
            _compute_Fmin_from_RMSD(rmsd_current, rmsd0, Fmin_min, Fmin_max)
            + stateB['fmin_offset'], Fmin_max
        )

        LOGGER.info('')
        LOGGER.info(
            f"Starting cycle {n + 1} with initial structure "
            f"{getTitle(b, 'structure B')}"
        )

        n_modes_B, rmsd_improvement_B, _step_log_B, _diag_info_B = calcStep(
            coordsB, coordsA, n_modes, ensB, defvecsB, rmsds,
            Fmin=Fmin_stepB, mask=maskB,
            f_hist=kwargs['f_histB'],
            min_slack_hist=kwargs['min_slack_histB'],
            worst_edges_hist=kwargs['worst_edges_histB'],
            adj_exempt=adj_exempt_B,
            runtime_state=runtime_stateB,
            **kwargs
        )
        if rmsd_improvement_B > 0.0 and _step_log_B:
            step_log_buffer.append((ensB.numCoordsets(), n, _step_log_B))

        if use_dynamic_fmin:
            if rmsd_improvement_B < min_rmsd_improve and n_modes_B > 0:
                stateB['stall_counter'] += 1
            else:
                stateB['stall_counter'] = 0
            if stateB['stall_counter'] >= stall_steps:
                max_offset = max(0.0, Fmin_max - Fmin_min)
                old        = stateB['fmin_offset']
                stateB['fmin_offset'] = min(old + Fmin_bump_increment, max_offset)
                if stateB['fmin_offset'] <= old:
                    LOGGER.info(
                        "B->A Fmin offset has reached its maximum (%.3f) after "
                        "%d small-improvement cycles; ending simulation "
                        "gracefully." % (stateB['fmin_offset'], stall_steps)
                    )
                    terminated_on_fmin_cap = True
                    break
                LOGGER.warning(
                    "B->A stalled. Increasing Fmin offset additively to "
                    f"{stateB['fmin_offset']:.3f}."
                )
                stateB['stall_counter'] = 0
        else:
            if rmsd_improvement_B < min_rmsd_improve and n_modes_B > 0:
                LOGGER.warning(
                    "B->A ΔRMSD=%.2e < min_rmsd_improve=%.2e. Stopping." %
                    (rmsd_improvement_B, min_rmsd_improve)
                )
                break

        n += 1
        n_modes = n_modes_B
        if n_modes == 0:
            LOGGER.report('Alternating Adaptive ANM converged in %.2fs.',
                          '_prody_calcAdaptiveANM')
            break

    if terminated_on_fmin_cap:
        LOGGER.info(
            "Alternating Adaptive ANM terminated because Fmin offset reached "
            "its maximum value (A: %.3f, B: %.3f)." %
            (stateA['fmin_offset'], stateB['fmin_offset'])
        )

    LOGGER.report('Alternating Adaptive ANM finished in %.2fs.',
                  '_prody_calcAdaptiveANM')

    _write_final_log(base_log_path, title, rmsd0, step_log_buffer)

    ensemble = ensA + ensB[::-1]
    return ensemble

# -------------------------------------------------------------------
# both-ways driver
# -------------------------------------------------------------------

def calcBothWaysAdaptiveANM(a, b, n_steps, **kwargs):
    """
    Run one-way A->B then one-way B->A and stitch into a single Ensemble.
    Each one-way pass writes its own final log if log_path is provided,
    with '_fwd' and '_rev' suffixes before '.log'.
    """
    base_log_path = kwargs.get('log_path', None)

    LOGGER.info('Running forward pass: A -> B')
    fwd_kwargs = dict(kwargs)
    if base_log_path:
        stem = base_log_path.rsplit('.log', 1)[0]
        fwd_kwargs['log_path'] = stem + '_fwd.log'
    forward_ensemble = calcOneWayAdaptiveANM(a, b, n_steps, **fwd_kwargs)

    LOGGER.info('Running reverse pass: B -> A')
    rev_kwargs = dict(kwargs)
    if base_log_path:
        rev_kwargs['log_path'] = stem + '_rev.log'
    reverse_ensemble = calcOneWayAdaptiveANM(b, a, n_steps, **rev_kwargs)

    full_ensemble = Ensemble(forward_ensemble.getTitle() + '_bothways')
    full_ensemble.setCoords(forward_ensemble.getCoords())
    full_ensemble.setAtoms(forward_ensemble.getAtoms())
    full_ensemble.setWeights(forward_ensemble.getWeights())

    for i in range(forward_ensemble.numCoordsets()):
        full_ensemble.addCoordset(forward_ensemble.getCoordsets()[i])
    for i in range(1, reverse_ensemble.numCoordsets()):
        full_ensemble.addCoordset(reverse_ensemble.getCoordsets()[-i])

    LOGGER.info('Both-way Adaptive ANM completed successfully.')
    return full_ensemble
