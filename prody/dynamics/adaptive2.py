# -*- coding: utf-8 -*-

"""
Adaptive ANM (adaptive2)
========================

This module implements an advanced, physically guided variant of Adaptive
Anisotropic Network Modeling (Adaptive ANM). It generates a realistic
transition pathway between two conformations by iteratively deforming one
structure toward the other using selected ENM/ANM modes and progressively
refining the deformation as RMSD decreases.

Compared to the legacy adaptive implementation, *adaptive2* introduces:

• **Sorted mode selection (overlap–ranked)**
  Modes are first ranked by absolute overlap with the deformation vector
  before applying the cumulative-overlap threshold (Fmin). This ensures that
  the most relevant modes are always considered first, greatly improving
  stability and physical relevance.

• **Adaptive Fmin scheduling (automatic mode-selection threshold)**
  Fmin is dynamically updated according to deformation progress.
  Large deformations → low Fmin → broad exploratory mode usage.
  Small deformations → higher Fmin → fine-grained, localized refinement.
  This eliminates the need to guess a good Fmin and avoids runaway mode
  accumulation.

• **Adaptive step-size scheduling ("adaptive f") [optional]**
  When `f_auto=True`, the global step size *f* is automatically adjusted
  according to deformation magnitude:
      – Far from target → larger steps
      – Near target → smaller, safer steps
  This complements backtracking and improves convergence smoothness.

• **Stringent, physically meaningful adjacency checks**
  Sequential distances are continuously monitored to prevent unphysical chain
  stretching or breakage. Naturally long native edges can be exempted, but all
  other adjacency violations are treated strictly. Any step that causes
  disconnection or unrealistic stretching is automatically rejected.

• **Backtracking-based step acceptance**
  If a trial deformation does not improve RMSD or violates adjacency
  constraints, the step size is halved repeatedly until a valid improvement is
  found—ensuring that every accepted update is physically reasonable.

• **Optional progressive mode cap (progressive_nmax)**
  The maximum number of allowable modes can increase as the deformation
  shrinks, stabilizing large-scale early steps while enabling high-resolution
  adjustments later in the trajectory.

Driver Modes
------------
The calculation can be performed in three modes:

- ``AANM_ONEWAY``  
    Deform structure A → B in a single direction.

- ``AANM_ALTERNATING``  
    Alternate between A → B and B → A, maintaining independent mode selection
    and adjacency checks for each direction.

- ``AANM_BOTHWAYS``  
    Perform a full one-way run A → B followed by B → A, stitching both
    trajectories into a single ensemble.

Core Functions
--------------
- ``calcAdaptiveANM``  
    Main dispatcher selecting one-way, alternating, or both-way execution.

- ``calcOneWayAdaptiveANM``  
    Forward adaptive deformation using adaptive Fmin, optional adaptive f,
    backtracking, and adjacency filtering.

- ``calcAlternatingAdaptiveANM``  
    Alternating-direction adaptive ANM with separate Fmin/f histories.

- ``calcBothWaysAdaptiveANM``  
    Performs two one-way transitions and merges the path.

- ``calcStep``  
    Core update routine performing ENM construction, overlap sorting,
    adaptive Fmin selection, adaptive/initial step-size handling,
    candidate-step evaluation, adjacency checks, and RMSD-based acceptance.

Parameter Summary (with Defaults)
---------------------------------
Below are the most commonly used parameters and their defaults
(as accepted by `calcStep`, drivers, or both):

General deformation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``n_modes`` (int, default: 20)  
    Starting number of modes to compute.

- ``n_max_modes`` (int or float, default: DOF)  
    Maximum number of allowable modes.  
    If ``0 < n_max_modes <= 1``, interpreted as a fraction of DOF.

- ``model`` (str, default: "anm")  
    ENM model used in mode computation.

Step-size parameters
~~~~~~~~~~~~~~~~~~~~
- ``f`` (float, default: 1.0 or 0.2 depending on driver)  
    Global step-size multiplier.

- ``f_auto`` (bool, default: False)  
    Enable adaptive step-size scheduling.

- ``f_max`` (float, default: 1.0)  
    Maximum adaptive f when ``f_auto=True``.

- ``f_min_sched`` (float, default: 0.2)  
    Minimum adaptive f when ``f_auto=True``.

- ``f_gamma`` (float, default: 1.0)  
    Controls curvature of adaptive f schedule.

Fmin (mode-selection) parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``Fmin`` (float or None, default: None)  
    If None, Fmin is computed adaptively.

- ``Fmin_min`` (float, default: 0.4)  
    Minimum allowable Fmin.

- ``Fmin_max`` (float, default: 0.9)  
    Maximum allowable Fmin.

- ``dynamic_fmin`` (bool, default: True)  
    Enable dynamic adjustment of Fmin based on stall detection.

- ``Fmin_bump_increment`` (float, default: 0.05)  
    Amount to raise Fmin when stall occurs.

- ``stall_steps`` (int, default: 3)  
    Number of zero/small-improvement steps required to trigger an Fmin bump.

Adjacency and physical constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``adj_max`` (float, default: 4.5)  
    Maximum allowed sequential distance.

- ``adj_tol`` (float, default: 0.0)  
    Additional tolerance for identifying exempt native edges.

- ``adj_exempt`` (array or None, default: None)  
    Boolean mask specifying which native edges are exempt from adjacency checks.

- ``cutoff`` (float, default: 15.0)  
    Distance threshold for non-bonded disconnection checks.

Backtracking and convergence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- ``max_backtracks`` (int, default: 5)  
    Maximum number of backtracking levels tried per step.

- ``min_f_eff`` (float, default: 1e-3)  
    Smallest acceptable effective step size.

- ``min_rmsd_step`` (float, default: 1e-3)  
    Minimum RMSD improvement required to accept a step unless stalling.

- ``min_rmsd_diff`` (float, default: 0.05)  
    Convergence criterion for RMSD change over multiple steps.

- ``target_rmsd`` (float, default: 1.0)  
    Stop when RMSD falls below this threshold.

Miscellaneous
~~~~~~~~~~~~~
- ``risk_topk`` (int, default: 8)  
    Number of highest-risk adjacency edges to log.

- ``mask`` (array or None)  
    Atom mask for mode computation.

Notes
-----
This module is fully compatible with ProDy’s ENM/ANM framework and supports
all ENM-related keyword arguments (e.g. cutoff, gamma, mass weighting).
Logging reports details of mode usage, Fmin evolution, adaptive f values,
adjacency slack, RMSD progression, and step acceptance or rejection.

"""

__author__ = 'Anupam Banerjee'
__credits__ = ['Anthony Bogetti']
__email__ = ['anupam.banerjee@stonybrook.edu', 'anthony.bogetti@stonybrook.edu']

from numbers import Integral
import numpy as np

from prody import LOGGER
import logging

class DropWarnings(logging.Filter):
    def filter(self, record):
        return record.levelno != logging.WARNING

LOGGER._logger.addFilter(DropWarnings())

from prody.atomic import Atomic, AtomMap
from prody.ensemble import Ensemble
from prody.measure import calcRMSD, superpose
from prody.utilities import getCoords, importLA

# Adjust import according to your project layout
try:
    from .functions import calcENM
except Exception:  # pragma: no cover
    from functions import calcENM

__all__ = [
    'calcAdaptiveANM', 'AANM_ONEWAY', 'AANM_ALTERNATING',
    'AANM_BOTHWAYS', 'AANM_DEFAULT'
]

AANM_ALTERNATING = 0
AANM_ONEWAY = 1
AANM_BOTHWAYS = 2
AANM_DEFAULT = AANM_ALTERNATING

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


def build_adj_exempt(native_coords, adj_max=4.5, tol=0.0):
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


def compute_slack_info(coords, adj_max=4.5, risk_topk=8, adj_exempt=None):
    """Return (slack array, min_slack, worst_edge_indices[:topk])."""
    elen = sequential_edge_lengths(coords)
    slack = float(adj_max) - elen
    if slack.size == 0:
        return slack, np.inf, np.array([], dtype=int)
    if adj_exempt is not None and adj_exempt.shape == slack.shape:
        slack = slack.copy()
        slack[adj_exempt] = np.inf
    order = np.argsort(slack)  # ascending: most risky first
    worst = order[:min(risk_topk, order.size)]
    return slack, float(slack.min()), worst


def checkDisconnection(coords, cutoff, adj_max=4.5, adj_exempt=None):
    """
    True if there is global isolation (>cutoff) or any non-exempt
    sequential edge > adj_max.
    """
    r = np.asarray(coords, float)
    N = r.shape[0]
    if N < 2:
        return False

    # Global isolation
    all_dists = np.sqrt(((r[:, None, :] - r[None, :, :]) ** 2).sum(axis=2))
    nn = np.empty(N, dtype=float)
    for i in range(N):
        row = all_dists[i]
        nn[i] = min(
            np.min(row[:i]) if i > 0 else np.inf,
            np.min(row[i + 1:]) if i + 1 < N else np.inf
        )
    isolated = nn > float(cutoff)

    # Sequential adjacency
    broken_adj = False
    if N >= 2:
        edge_d = sequential_edge_lengths(coords)
        if adj_exempt is not None and adj_exempt.shape == edge_d.shape:
            broken_adj = np.any((edge_d > float(adj_max)) & (~adj_exempt))
        else:
            broken_adj = np.any(edge_d > float(adj_max))

    return bool(np.any(isolated) or broken_adj)

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

    rmsd = calcRMSD(coordsA, coordsB, weights)
    LOGGER.info('Initialized Adaptive ANM with RMSD {:4.3f}\n'.format(rmsd))

    return coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd


def getTitle(structure, def_title='structure'):
    return structure.getTitle() if isinstance(structure, Atomic) else def_title


def computeAdaptiveF(defvecs, f_max=1.0, f_min=0.2, gamma=1.0):
    """
    Optional global f scheduler (kept for compatibility).
    Used only if f_auto=True.
    """
    if not defvecs:
        return float(f_max)

    d0, dk = defvecs[0], defvecs[-1]
    n0, nk = norm(d0), norm(dk)
    if n0 <= 1e-12:
        return float(f_min)

    progress = 1.0 - np.sqrt(max(nk, 0.0) / max(n0, 1e-12))
    progress = float(np.clip(progress, 0.0, 1.0))
    return float(f_min + (f_max - f_min) * (1.0 - (progress ** float(gamma))))


def computeProgressiveNmax(defvecs, n_max_modes, frac0=0.2, eta=1.5):
    """
    Progressive cap for high-frequency modes.

    If progressive_nmax=True, the effective cap on the number of modes
    starts as n_max_modes * frac0 and grows smoothly toward n_max_modes
    as the defvec shrinks.
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

    frac = float(frac0) + (1.0 - float(frac0)) * (progress ** float(eta))
    eff_n = int(np.floor(n_max_modes * frac))
    return int(max(1, min(n_max_modes, eff_n)))


def checkConvergence(rmsds, coords, **kwargs):
    min_rmsd_diff = kwargs.get('min_rmsd_diff', 0.05)
    target_rmsd = kwargs.get('target_rmsd', 1.0)
    cutoff = kwargs.get('cutoff', 15.0)
    adj_max = kwargs.get('adj_max', 4.5)
    adj_exempt = kwargs.get('adj_exempt', None)

    if len(rmsds) > 4 and np.all(np.abs(np.diff(rmsds[-4:])) < min_rmsd_diff):
        LOGGER.warn(f'RMSD change fell below {min_rmsd_diff}')
        return True

    if rmsds[-1] < target_rmsd:
        LOGGER.info(
            "Convergence reached: RMSD %.3f is below target %.3f." %
            (rmsds[-1], target_rmsd)
        )
        return True

    if checkDisconnection(coords, cutoff, adj_max=adj_max, adj_exempt=adj_exempt):
        LOGGER.warn('Disconnection found in the structure.')
        return True

    return False

# -------------------------------------------------------------------
# RMSD-based Fmin schedule
# -------------------------------------------------------------------

def _compute_Fmin_from_RMSD(rmsd_current, rmsd0, Fmin_min=0.4, Fmin_max=0.9):
    """
    Very lightweight Fmin schedule (v13):

    - At RMSD = rmsd0: Fmin ~ Fmin_min
    - At RMSD → 0:     Fmin → Fmin_max
    """
    if rmsd0 <= 1e-8:
        progress = 1.0
    else:
        progress = 1.0 - float(rmsd_current / rmsd0)
    progress = float(np.clip(progress, 0.0, 1.0))

    Fmin = float(Fmin_min + (Fmin_max - Fmin_min) * progress)
    return float(min(Fmin, Fmin_max))

# -------------------------------------------------------------------
# core single step: best-improvement backtracking, progressive_nmax
# -------------------------------------------------------------------

def calcStep(initial, target, n_modes, ensemble, defvecs, rmsds, Fmin,
             mask=None, callback_func=None, **kwargs):
    """
    Single adaptive step (v13 core, v12b-style logging):

    - Uses Fmin passed from the driver (RMSD-based + dynamic bump).
    - Optional progressive_nmax.
    - Simple best-improvement backtracking (no Fmin tweaks inside).
    - No culprit-aware rewind or deep history scanning.
    """

    # Basic parameters / defaults
    f_global = float(kwargs.get('f', 0.2))
    Fmin_max = float(kwargs.get('Fmin_max', 0.9))
    cutoff = float(kwargs.get('cutoff', 15.0))
    adj_max = float(kwargs.get('adj_max', 4.5))
    max_backtracks = int(kwargs.get('max_backtracks', 5))
    min_f_eff = float(kwargs.get('min_f_eff', 1e-3))
    min_rmsd_step = float(kwargs.get('min_rmsd_step', 1e-3))
    adj_exempt = kwargs.get('adj_exempt', None)
    progressive_nmax_flag = bool(kwargs.get('progressive_nmax', False))
    nmax_frac0 = float(kwargs.get('nmax_frac0', 0.2))
    nmax_eta = float(kwargs.get('nmax_eta', 1.5))

    # kept for compatibility / logging
    risk_topk = int(kwargs.get('risk_topk', 8))

    f_hist = kwargs.setdefault('f_hist', [])
    min_slack_hist = kwargs.setdefault('min_slack_hist', [])
    worst_edges_hist = kwargs.setdefault('worst_edges_hist', [])

    weights = ensemble.getWeights()
    if weights is not None:
        weights = weights.flatten()

    coords_init = initial
    coords_tar = target

    dof = coords_init.size - 6
    if dof <= 0:
        LOGGER.warn("Not enough DOF; returning without update.")
        return n_modes, 0.0

    # base cap for modes
    raw_n_max_modes = kwargs.get('n_max_modes', dof)
    if isinstance(raw_n_max_modes, float) and 0 < raw_n_max_modes <= 1:
        n_max_modes = max(1, int(raw_n_max_modes * dof))
    else:
        n_max_modes = min(dof, int(raw_n_max_modes))

    # defvec (for progressive_nmax and f_auto)
    defvec = coords_tar - coords_init
    d = defvec.flatten()
    if weights is not None:
        d *= weights.repeat(3)
    defvecs.append(d)

    if kwargs.get('f_auto', False):
        f_global = computeAdaptiveF(
            defvecs,
            f_max=kwargs.get('f_max', 1.0),
            f_min=kwargs.get('f_min_sched', 0.2),
            gamma=kwargs.get('f_gamma', 1.0)
        )
        kwargs['f'] = f_global

    if progressive_nmax_flag:
        eff_n_max_modes = computeProgressiveNmax(
            defvecs, n_max_modes,
            frac0=nmax_frac0,
            eta=nmax_eta
        )
    else:
        eff_n_max_modes = n_max_modes
    eff_n_max_modes = int(max(1, min(n_max_modes, eff_n_max_modes)))

    n_modes = min(max(1, int(n_modes)), eff_n_max_modes)
    model = kwargs.get('model', 'anm')

    # cap Fmin
    Fmin = float(min(Fmin, Fmin_max))

    # ENM calculation
    enm_kwargs = dict(kwargs)
    enm_kwargs.pop('model', None)
    enm_kwargs.pop('n_modes', None)
    enm_kwargs.pop('trim', None)

    anm_h, _ = calcENM(
        coords_init,
        select=mask,
        mask=mask,
        model=model,
        trim='trim',
        n_modes=eff_n_max_modes,
        **enm_kwargs
    )
    if mask is not None:
        anm_h.masked = False

    all_eigvecs = anm_h.getEigvecs()[:, :eff_n_max_modes]

    # project defvec onto modes
    d_loc = (coords_tar - coords_init).reshape(-1)
    if weights is not None:
        d_loc *= weights.repeat(3)

    norm_d = norm(d_loc) or 1.0
    ov = np.dot(d_loc, all_eigvecs)
    norm_ov = ov / norm_d

    sort_idx = np.argsort(-np.abs(norm_ov))
    c_sq_sorted = np.cumsum(norm_ov[sort_idx] ** 2)
    pick_mask = c_sq_sorted <= Fmin
    if not np.any(pick_mask):
        pick_mask[0] = True

    sel_idx = sort_idx[pick_mask]
    sel_ov = ov[sel_idx]
    sel_vecs = all_eigvecs[:, sel_idx]

    n_used = int(sel_idx.size)
    if n_used == 0:
        LOGGER.warn("No modes selected (n_used == 0); returning without update.")
        return n_modes, 0.0

    v = sel_vecs.dot(sel_ov)
    denom = float(np.dot(v, v))
    if denom == 0.0:
        LOGGER.warn("Degenerate step direction (v·v == 0); returning without update.")
        return n_modes, 0.0

    s_base = float(np.dot(v, d_loc) / denom)
    v3 = v.reshape(-1, 3)

    rmsd_before_step = rmsds[-1]

    # ---- best-improvement backtracking on f ----
    accepted = False
    coords_updated = coords_init.copy()
    accepted_f_eff = 0.0
    used_h = 0

    best_improvement = 0.0
    best_cand = None
    best_f_eff = None
    best_k = None

    for k in range(int(max_backtracks) + 1):
        f_eff = float(f_global) / (2 ** k)
        if f_eff < min_f_eff:
            break

        s_step = f_eff * s_base
        cand = coords_init + s_step * v3

        if checkDisconnection(cand, cutoff, adj_max=adj_max, adj_exempt=adj_exempt):
            continue

        trial_rmsd = calcRMSD(cand, coords_tar, weights)
        improvement = rmsd_before_step - trial_rmsd
        if improvement <= 0.0:
            continue

        if improvement > best_improvement:
            best_improvement = improvement
            best_cand = cand
            best_f_eff = f_eff
            best_k = k

    if best_cand is not None:
        if best_improvement < min_rmsd_step:
            LOGGER.info(
                "Accepting sub-threshold improvement ΔRMSD=%.4g "
                "(min_rmsd_step=%.2e) to avoid stall." %
                (best_improvement, min_rmsd_step)
            )
        accepted = True
        coords_updated = best_cand
        accepted_f_eff = best_f_eff
        used_h = best_k

    if not accepted:
        LOGGER.warn('All step attempts failed; exiting step with no update.')
        return n_modes, 0.0

    # ---- update coords, ensemble, histories, logging ----
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

    initial[:] = coords_updated
    rmsd_after_step = calcRMSD(initial, coords_tar, weights)
    rmsds.append(rmsd_after_step)
    rmsd_improvement = rmsd_before_step - rmsd_after_step

    h = int(used_h)
    scale_str = f"1/2^{h}"

    if sel_idx.size > 0:
        highest_mode = int(sel_idx.max()) + 1
    else:
        highest_mode = 0

    LOGGER.info(
        "Step successful [h=%d, scale=%s] "
        "(Fmin=%.3f, f_eff=%.3f, highest_mode=%d, ΔRMSD=%.4g)" %
        (h, scale_str, Fmin, accepted_f_eff, highest_mode, rmsd_improvement)
    )

    if n_max_modes and n_max_modes > 0:
        nmax_frac = float(eff_n_max_modes) / float(n_max_modes)
    else:
        nmax_frac = 1.0

    LOGGER.info(
        "Step summary (Fmin=%.3f, n_used=%d, nmax_frac=%.3f)" %
        (Fmin, n_used, nmax_frac)
    )
    LOGGER.info('Current RMSD is %.6f\n' % rmsd_after_step)

    if checkConvergence(rmsds, initial, **kwargs):
        n_modes_out = 0
    else:
        n_modes_out = max(1, min(int(n_modes), n_used))

    return n_modes_out, rmsd_improvement

# -------------------------------------------------------------------
# dispatcher
# -------------------------------------------------------------------

def calcAdaptiveANM(a, b, n_steps, mode=AANM_DEFAULT, **kwargs):
    if not isinstance(n_steps, Integral):
        raise TypeError('n_steps must be an integer')

    if mode == AANM_ONEWAY:
        return calcOneWayAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == AANM_ALTERNATING:
        return calcAlternatingAdaptiveANM(a, b, n_steps, **kwargs)
    elif mode == AANM_BOTHWAYS:
        return calcBothWaysAdaptiveANM(a, b, n_steps, **kwargs)
    else:
        raise ValueError('unknown aANM mode: %d' % mode)

# -------------------------------------------------------------------
# one-way driver with dynamic Fmin bumping
# -------------------------------------------------------------------

def calcOneWayAdaptiveANM(a, b, n_steps, **kwargs):
    """
    One-way adaptive ANM with:
    - v13 RMSD-based Fmin ramp.
    - optional progressive_nmax.
    - v12b-style log lines.
    - additive dynamic Fmin bumping (+Fmin_bump_increment).
    """
    n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd0 = checkInput(
        a, b, **kwargs
    )
    coordsA = coordsA.copy()

    LOGGER.timeit('_prody_calcAdaptiveANM')

    n = 0
    defvecs = []
    rmsds = [rmsd0]

    adj_max = kwargs.get('adj_max', 4.5)
    adj_tol = kwargs.get('adj_tol', 0.0)
    if 'adj_exempt' not in kwargs or kwargs['adj_exempt'] is None:
        kwargs['adj_exempt'] = build_adj_exempt(coordsA, adj_max=adj_max, tol=adj_tol)

    f0 = kwargs.get('f', 1.0)
    kwargs.setdefault('f_hist', [float(f0)])
    kwargs.setdefault('min_slack_hist', [])
    kwargs.setdefault('worst_edges_hist', [])

    s0, ms0, worst0 = compute_slack_info(
        coordsA,
        adj_max=adj_max,
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

    Fmin_min = kwargs.get('Fmin_min', 0.4)
    Fmin_max = kwargs.get('Fmin_max', 0.9)
    min_rmsd_improve = float(kwargs.get('min_rmsd_improve', 1e-4))

    use_dynamic_fmin = bool(kwargs.get('dynamic_fmin', True))
    stall_steps = int(kwargs.get('stall_steps', 3))
    Fmin_bump_increment = float(kwargs.get('Fmin_bump_increment', 0.05))

    stall_counter = 0
    fmin_offset = 0.0
    
    # Track whether we stopped because Fmin offset saturated
    terminated_on_fmin_cap = False

    # n = user-visible cycle index (only successful/terminal cycles)
    n = 0
    # attempt_idx = raw attempts (always increments, limits total tries)
    attempt_idx = 0

    while attempt_idx < n_steps:
        rmsd_current = rmsds[-1]
        base_Fmin = _compute_Fmin_from_RMSD(
            rmsd_current, rmsd0,
            Fmin_min=Fmin_min, Fmin_max=Fmin_max
        )
        Fmin_step = min(base_Fmin + fmin_offset, Fmin_max)

        # internal attempt (no "Starting cycle" log yet)
        n_modes, rmsd_improvement = calcStep(
            coordsA, coordsB, n_modes, ensemble, defvecs, rmsds,
            Fmin=Fmin_step, mask=maskA, **kwargs
        )

        attempt_idx += 1  # always advance attempt count

        # If this attempt did nothing and there are still modes, treat as failed.
        if rmsd_improvement == 0.0 and n_modes > 0:
            if use_dynamic_fmin:
                stall_counter += 1
                if stall_counter >= stall_steps:
                    max_offset = max(0.0, Fmin_max - Fmin_min)
                    old_offset = fmin_offset
                    fmin_offset = min(fmin_offset + Fmin_bump_increment,
                                      max_offset)
                    if fmin_offset <= old_offset:
                        LOGGER.info(
                            "Fmin offset has reached its maximum (%.3f) after "
                            "%d zero-improvement attempts; ending simulation "
                            "gracefully." % (fmin_offset, stall_steps)
                        )
                        terminated_on_fmin_cap = True
                        break

                    LOGGER.warn(
                        f"Stalled for {stall_steps} zero-improvement attempts. "
                        f"Increasing Fmin offset additively to "
                        f"{fmin_offset:.3f}."
                    )
                    stall_counter = 0
            # No logging of cycle index; just try another attempt
            continue

        # From here: either made progress (rmsd_improvement != 0) or n_modes == 0

        # If the step signaled convergence (e.g. RMSD < target, or disconnection),
        # stop *before* starting a new cycle.
        if n_modes == 0:
            LOGGER.report(
                'One-way Adaptive ANM converged in %.2fs.',
                '_prody_calcAdaptiveANM'
            )
            break

        # Now we know this is a "real" non-terminal cycle → log it with cycle index n+1
        LOGGER.info('')
        LOGGER.info(f"Starting cycle {n + 1} with initial structure {title}")

        if use_dynamic_fmin:
            if rmsd_improvement < min_rmsd_improve and n_modes > 0:
                stall_counter += 1
                LOGGER.info(
                    f"Stall detected (improvement={rmsd_improvement:.2e}). "
                    f"Counter: {stall_counter}/{stall_steps}"
                )
            else:
                stall_counter = 0

            if stall_counter >= stall_steps:
                max_offset = max(0.0, Fmin_max - Fmin_min)
                old_offset = fmin_offset
                fmin_offset = min(fmin_offset + Fmin_bump_increment,
                                  max_offset)
                if fmin_offset <= old_offset:
                    LOGGER.info(
                        "Fmin offset has reached its maximum (%.3f) after "
                        "%d small-improvement cycles; ending simulation "
                        "gracefully." % (fmin_offset, stall_steps)
                    )
                    terminated_on_fmin_cap = True
                    break

                LOGGER.warn(
                    f"Stalled for {stall_steps} small-improvement cycles. "
                    f"Increasing Fmin offset additively to {fmin_offset:.3f}."
                )
                stall_counter = 0
        else:
            if rmsd_improvement < min_rmsd_improve and n_modes > 0:
                LOGGER.warn(
                    "RMSD improvement %.2e < min_rmsd_improve=%.2e. Stopping." %
                    (rmsd_improvement, min_rmsd_improve)
                )
                # Count this as the last cycle and exit
                n += 1
                break

        # Successful / terminal cycle → advance visible cycle index
        n += 1
        
    if terminated_on_fmin_cap:
        LOGGER.info(
            "One-way Adaptive ANM terminated because Fmin offset reached its "
            "maximum value (%.3f)." % fmin_offset
        )
        
    LOGGER.report(
        'One-way Adaptive ANM finished in %.2fs.',
        '_prody_calcAdaptiveANM'
    )
    return ensemble

# -------------------------------------------------------------------
# alternating driver with same schedule + dynamic bump
# -------------------------------------------------------------------

def calcAlternatingAdaptiveANM(a, b, n_steps, **kwargs):
    """
    Alternating A->B and B->A driver:

    - Uses same RMSD-based Fmin schedule for both directions.
    - Optional progressive_nmax.
    - Optional dynamic Fmin bumping (+Fmin_bump_increment) per direction.
    - Uses native adjacency exemptions for both A and B separately.
    """
    n_modes = kwargs.pop('n_modes', 20)

    coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd0 = checkInput(
        a, b, **kwargs
    )
    coordsA, coordsB = coordsA.copy(), coordsB.copy()

    LOGGER.timeit('_prody_calcAdaptiveANM')

    n = 0
    defvecsA, defvecsB = [], []
    rmsds = [rmsd0]

    adj_max = kwargs.get('adj_max', 4.5)
    adj_tol = kwargs.get('adj_tol', 0.0)

    # Build adjacency exemptions separately for A and B, while
    # remaining backward-compatible with a single adj_exempt mask.
    adj_exempt_common = kwargs.get('adj_exempt', None)

    if 'adj_exempt_A' in kwargs:
        adj_exempt_A = kwargs['adj_exempt_A']
    elif adj_exempt_common is not None:
        adj_exempt_A = adj_exempt_common
    else:
        adj_exempt_A = build_adj_exempt(coordsA, adj_max=adj_max, tol=adj_tol)

    if 'adj_exempt_B' in kwargs:
        adj_exempt_B = kwargs['adj_exempt_B']
    elif adj_exempt_common is not None:
        adj_exempt_B = adj_exempt_common
    else:
        adj_exempt_B = build_adj_exempt(coordsB, adj_max=adj_max, tol=adj_tol)

    f0 = kwargs.get('f', 1.0)
    kwargs.setdefault('f_histA', [float(f0)])
    kwargs.setdefault('f_histB', [float(f0)])
    kwargs.setdefault('min_slack_histA', [])
    kwargs.setdefault('min_slack_histB', [])
    kwargs.setdefault('worst_edges_histA', [])
    kwargs.setdefault('worst_edges_histB', [])

    sA0, msA0, worstA0 = compute_slack_info(
        coordsA, adj_max=adj_max,
        risk_topk=kwargs.get('risk_topk', 8),
        adj_exempt=adj_exempt_A
    )
    kwargs['min_slack_histA'].append(msA0)
    kwargs['worst_edges_histA'].append(worstA0.tolist())

    sB0, msB0, worstB0 = compute_slack_info(
        coordsB, adj_max=adj_max,
        risk_topk=kwargs.get('risk_topk', 8),
        adj_exempt=adj_exempt_B
    )
    kwargs['min_slack_histB'].append(msB0)
    kwargs['worst_edges_histB'].append(worstB0.tolist())

    ensA = Ensemble(getTitle(a, 'structureA') + '_aANM')
    ensA.setAtoms(atoms)
    ensA.setCoords(coordsB)
    ensA.setWeights(weights)
    ensA.addCoordset(coordsA.copy())

    ensB = Ensemble(getTitle(b, 'structureB') + '_aANM')
    ensB.setAtoms(atoms)
    ensB.setCoords(coordsA)
    ensB.setWeights(weights)
    ensB.addCoordset(coordsB.copy())

    Fmin_min = kwargs.get('Fmin_min', 0.4)
    Fmin_max = kwargs.get('Fmin_max', 0.9)
    min_rmsd_improve = float(kwargs.get('min_rmsd_improve', 1e-4))

    use_dynamic_fmin = bool(kwargs.get('dynamic_fmin', True))
    stall_steps = int(kwargs.get('stall_steps', 3))
    Fmin_bump_increment = float(kwargs.get('Fmin_bump_increment', 0.05))

    stateA = {'stall_counter': 0, 'fmin_offset': 0.0}
    stateB = {'stall_counter': 0, 'fmin_offset': 0.0}

    terminated_on_fmin_cap = False

    while n < n_steps:
        # ---------- A -> B ----------
        rmsd_current = rmsds[-1]
        base_FminA = _compute_Fmin_from_RMSD(
            rmsd_current, rmsd0,
            Fmin_min=Fmin_min, Fmin_max=Fmin_max
        )
        Fmin_stepA = min(base_FminA + stateA['fmin_offset'], Fmin_max)

        LOGGER.info('')
        LOGGER.info(
            f"Starting cycle {n + 1} with initial structure "
            f"{getTitle(a, 'structure A')}"
        )

        n_modes_A, rmsd_improvement_A = calcStep(
            coordsA, coordsB, n_modes, ensA, defvecsA, rmsds,
            Fmin=Fmin_stepA, mask=maskA,
            f_hist=kwargs['f_histA'],
            min_slack_hist=kwargs['min_slack_histA'],
            worst_edges_hist=kwargs['worst_edges_histA'],
            adj_exempt=adj_exempt_A,
            **kwargs
        )

        if use_dynamic_fmin:
            if rmsd_improvement_A < min_rmsd_improve and n_modes_A > 0:
                stateA['stall_counter'] += 1
            else:
                stateA['stall_counter'] = 0

            if stateA['stall_counter'] >= stall_steps:
                max_offset = max(0.0, Fmin_max - Fmin_min)
                old_offset = stateA['fmin_offset']
                stateA['fmin_offset'] = min(
                    stateA['fmin_offset'] + Fmin_bump_increment,
                    max_offset
                )
                if stateA['fmin_offset'] <= old_offset:
                    LOGGER.info(
                        "A->B Fmin offset has reached its maximum (%.3f) after "
                        "%d small-improvement cycles; ending simulation "
                        "gracefully." % (stateA['fmin_offset'], stall_steps)
                    )
                    terminated_on_fmin_cap = True
                    break

                LOGGER.warn(
                    "A->B stalled. Increasing Fmin offset additively to "
                    f"{stateA['fmin_offset']:.3f}."
                )
                stateA['stall_counter'] = 0
        else:
            if rmsd_improvement_A < min_rmsd_improve and n_modes_A > 0:
                LOGGER.warn(
                    "A->B ΔRMSD=%.2e < min_rmsd_improve=%.2e. Stopping." %
                    (rmsd_improvement_A, min_rmsd_improve)
                )
                break

        if n_modes_A == 0:
            LOGGER.report(
                'Alternating Adaptive ANM converged in %.2fs.',
                '_prody_calcAdaptiveANM'
            )
            break

        # ---------- B -> A ----------
        rmsd_current = rmsds[-1]
        base_FminB = _compute_Fmin_from_RMSD(
            rmsd_current, rmsd0,
            Fmin_min=Fmin_min, Fmin_max=Fmin_max
        )
        Fmin_stepB = min(base_FminB + stateB['fmin_offset'], Fmin_max)

        LOGGER.info('')
        LOGGER.info(
            f"Starting cycle {n + 1} with initial structure "
            f"{getTitle(b, 'structure B')}"
        )

        n_modes_B, rmsd_improvement_B = calcStep(
            coordsB, coordsA, n_modes, ensB, defvecsB, rmsds,
            Fmin=Fmin_stepB, mask=maskB,
            f_hist=kwargs['f_histB'],
            min_slack_hist=kwargs['min_slack_histB'],
            worst_edges_hist=kwargs['worst_edges_histB'],
            adj_exempt=adj_exempt_B,
            **kwargs
        )

        if use_dynamic_fmin:
            if rmsd_improvement_B < min_rmsd_improve and n_modes_B > 0:
                stateB['stall_counter'] += 1
            else:
                stateB['stall_counter'] = 0

            if stateB['stall_counter'] >= stall_steps:
                max_offset = max(0.0, Fmin_max - Fmin_min)
                old_offset = stateB['fmin_offset']
                stateB['fmin_offset'] = min(
                    stateB['fmin_offset'] + Fmin_bump_increment,
                    max_offset
                )
                if stateB['fmin_offset'] <= old_offset:
                    LOGGER.info(
                        "B->A Fmin offset has reached its maximum (%.3f) after "
                        "%d small-improvement cycles; ending simulation "
                        "gracefully." % (stateB['fmin_offset'], stall_steps)
                    )
                    terminated_on_fmin_cap = True
                    break

                LOGGER.warn(
                    "B->A stalled. Increasing Fmin offset additively to "
                    f"{stateB['fmin_offset']:.3f}."
                )
                stateB['stall_counter'] = 0
        else:
            if rmsd_improvement_B < min_rmsd_improve and n_modes_B > 0:
                LOGGER.warn(
                    "B->A ΔRMSD=%.2e < min_rmsd_improve=%.2e. Stopping." %
                    (rmsd_improvement_B, min_rmsd_improve)
                )
                break

        n += 1
        n_modes = n_modes_B
        if n_modes == 0:
            LOGGER.report(
                'Alternating Adaptive ANM converged in %.2fs.',
                '_prody_calcAdaptiveANM'
            )
            break

    if terminated_on_fmin_cap:
        LOGGER.info(
            "Alternating Adaptive ANM terminated because Fmin offset reached "
            "its maximum value (A: %.3f, B: %.3f)." %
            (stateA['fmin_offset'], stateB['fmin_offset'])
        )

    LOGGER.report(
        'Alternating Adaptive ANM finished in %.2fs.',
        '_prody_calcAdaptiveANM'
    )
    ensemble = ensA + ensB[::-1]
    return ensemble

# -------------------------------------------------------------------
# both-ways driver (forward + reverse one-way, stitched)
# -------------------------------------------------------------------

def calcBothWaysAdaptiveANM(a, b, n_steps, **kwargs):
    """
    Convenience driver:
    - run one-way v12e A->B,
    - then one-way v12e B->A,
    - stitch trajectories into a single Ensemble.
    """
    LOGGER.info('Running forward pass: A -> B')
    forward_ensemble = calcOneWayAdaptiveANM(a, b, n_steps, **kwargs)

    LOGGER.info('Running reverse pass: B -> A')
    reverse_ensemble = calcOneWayAdaptiveANM(b, a, n_steps, **kwargs)

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
