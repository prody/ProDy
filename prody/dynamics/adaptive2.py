# -*- coding: utf-8 -*-

import logging
from numbers import Integral

import numpy as np
from prody import LOGGER
from prody.atomic import Atomic, AtomMap
from prody.ensemble import Ensemble
from prody.measure import calcRMSD, superpose
from prody.utilities import getCoords, importLA

try:
    # future: when this file lives inside prody/dynamics
    from .functions import calcENM
except Exception:  # pragma: no cover
    # current: use installed ProDy
    from prody.dynamics.functions import calcENM

AANM_ALTERNATING = 0
AANM_ONEWAY = 1
AANM_BOTHWAYS = 2
AANM_DEFAULT = AANM_ALTERNATING

__all__ = [
    "AdaptiveANM",
    "AANM_ONEWAY",
    "AANM_ALTERNATING",
    "AANM_BOTHWAYS",
    "AANM_DEFAULT",
]


norm = importLA().norm


class DropWarnings(logging.Filter):
    def filter(self, record):
        return record.levelno != logging.WARNING


LOGGER._logger.addFilter(DropWarnings())


class AdaptiveANM:
    """Class-based implementation of the adaptive anisotropic network model (adaptive2).

    This class generates a transition pathway between two conformations (A and B)
    by iteratively deforming one structure toward the other using ENM/ANM modes.
    It encapsulates all state (coordinates, masks, weights, RMSD history, Fmin and
    step-size schedules, adjacency slack, and logging) that was previously spread
    across module-level functions.

    Parameters
    ----------
    a, b : Atomic or array-like
        Input structures or coordinate arrays representing the start (A) and
        target (B) conformations. If `aligned` is False, A is first superposed
        onto B and the initial RMSD is logged.
    aligned : bool, optional
        If True, structures are assumed to be pre-aligned and no superposition
        is performed.
    n_modes : int, optional
        Initial number of ANM modes to consider in the deformation.
    model : str, optional
        ENM model type passed through to the normal-mode calculator (e.g. "anm").
    n_max_modes : int or float, optional
        Maximum number of modes allowed. If a float with 0 < value <= 1, it is
        interpreted as a fraction of the total degrees of freedom.
    f : float, optional
        Global step-size multiplier used in the backtracking line search.
    f_auto : bool, optional
        If True, enables an adaptive schedule for `f` based on deformation
        progress.
    f_max, f_min_sched, f_gamma : float, optional
        Parameters controlling the adaptive step-size schedule when `f_auto`
        is enabled.
    Fmin : float or None, optional
        Mode-selection cumulative-overlap threshold. If None, a lightweight
        RMSD-based schedule between `Fmin_min` and `Fmin_max` is used.
    Fmin_min, Fmin_max : float, optional
        Lower and upper bounds for the adaptive Fmin schedule.
    dynamic_fmin : bool, optional
        If True, enables additive Fmin bumping in response to stalls in RMSD
        improvement.
    Fmin_bump_increment : float, optional
        Amount added to the effective Fmin offset each time a stall condition
        is triggered.
    stall_steps : int, optional
        Number of consecutive zero- or small-improvement attempts required to
        trigger a dynamic Fmin bump.
    adj_max : float, optional
        Maximum allowed distance for sequential (i,i+1) edges; longer edges
        are treated as adjacency violations unless exempted.
    adj_tol : float, optional
        Extra tolerance for defining native sequential edges that are exempt
        from later adjacency checks.
    adj_exempt : array-like or None, optional
        Common boolean mask for exempting specific native edges from adjacency
        checks. If None, a mask is built from the native coordinates.
    cutoff : float, optional
        Maximum nearest-neighbor distance allowed before an atom is considered
        isolated, used in global disconnection checks.
    max_backtracks : int, optional
        Maximum number of halving steps in the backtracking line search per
        update cycle.
    min_f_eff : float, optional
        Minimum effective step size `f / 2^k` allowed before abandoning a
        trial update.
    min_rmsd_step : float, optional
        Minimum acceptable RMSD improvement in a single step. Smaller but
        positive improvements may still be accepted to avoid stalling.
    min_rmsd_diff : float, optional
        Convergence threshold on the change in RMSD over several recent steps.
    target_rmsd : float, optional
        RMSD target; if the running RMSD drops below this value, the run is
        considered converged.
    min_rmsd_improve : float, optional
        Minimum RMSD improvement per cycle used by the drivers to detect
        stalling when `dynamic_fmin` is disabled.
    progressive_nmax : bool, optional
        If True, progressively increases the effective cap on the number of
        modes used as the deformation shrinks.
    nmax_frac0 : float, optional
        Initial fraction of `n_max_modes` used when the deformation is large
        in progressive-nmax scheduling.
    nmax_eta : float, optional
        Exponent controlling how quickly the progressive-nmax cap approaches
        `n_max_modes` as the deformation shrinks.
    risk_topk : int, optional
        Number of highest-risk adjacency edges to track in the slack history.
    mask, maskA, maskB : array-like or None, optional
        Optional atom masks for ENM construction and projection for the
        combined system or separately for A and B.
    adj_exempt_A, adj_exempt_B : array-like or None, optional
        Separate adjacency exemption masks for the native A and B structures.
    **enm_kwargs :
        Additional keyword arguments forwarded to the ENM/ANM construction
        routine.

    Attributes
    ----------
    coordsA, coordsB : ndarray
        Working coordinate arrays for the current A and B structures.
    atoms : Atomic or None
        Underlying atomic object, if available, for building ensembles.
    weights : ndarray or None
        Per-atom weights derived from mapped flags on A and B, if present.
    rmsds : list of float
        History of RMSD values between current coordinates and the target.
    defvecsA, defvecsB : list of ndarray
        Deformation vectors accumulated for A->B and B->A runs.
    f_hist, f_histA, f_histB : list of float
        Accepted effective step sizes recorded globally and per direction.
    min_slack_hist* : list of float
        Minimum adjacency slack values per accepted update.
    worst_edges_hist* : list of list[int]
        Indices of the highest-risk sequential edges per accepted update.

    Methods
    -------
    run(n_steps, mode=AANM_DEFAULT)
        Execute the adaptive ANM deformation for a given number of steps
        in one-way, alternating, or both-ways mode, returning a ProDy
        Ensemble containing the generated trajectory.
    """
    def __init__(
        self,
        a,
        b,
        *,
        aligned=False,
        n_modes=20,
        model="anm",
        n_max_modes=None,
        # step-size control
        f=0.2,
        f_auto=False,
        f_max=1.0,
        f_min_sched=0.2,
        f_gamma=1.0,
        # Fmin control
        Fmin=None,
        Fmin_min=0.4,
        Fmin_max=0.9,
        dynamic_fmin=True,
        Fmin_bump_increment=0.05,
        stall_steps=3,
        # adjacency / physical constraints
        adj_max=6.0,
        adj_tol=0.0,
        adj_exempt=None,
        cutoff=15.0,
        # backtracking / convergence
        max_backtracks=5,
        min_f_eff=1e-3,
        min_rmsd_step=1e-3,
        min_rmsd_diff=0.05,
        target_rmsd=1.0,
        min_rmsd_improve=1e-4,
        # progressive_nmax and logging
        progressive_nmax=False,
        nmax_frac0=0.2,
        nmax_eta=1.5,
        risk_topk=8,
        # masks
        mask=None,
        maskA=None,
        maskB=None,
        adj_exempt_A=None,
        adj_exempt_B=None,
        **enm_kwargs,
    ):
        # basic config
        self.a = a
        self.b = b
        self.aligned = aligned
        self.n_modes = int(n_modes)
        self.model = model
        self.n_max_modes = n_max_modes
        self.f = float(f)
        self.f_auto = bool(f_auto)
        self.f_max = float(f_max)
        self.f_min_sched = float(f_min_sched)
        self.f_gamma = float(f_gamma)

        self.Fmin = Fmin
        self.Fmin_min = float(Fmin_min)
        self.Fmin_max = float(Fmin_max)
        self.dynamic_fmin = bool(dynamic_fmin)
        self.Fmin_bump_increment = float(Fmin_bump_increment)
        self.stall_steps = int(stall_steps)

        self.adj_max = float(adj_max)
        self.adj_tol = float(adj_tol)
        self.adj_exempt = adj_exempt
        self.cutoff = float(cutoff)

        self.max_backtracks = int(max_backtracks)
        self.min_f_eff = float(min_f_eff)
        self.min_rmsd_step = float(min_rmsd_step)
        self.min_rmsd_diff = float(min_rmsd_diff)
        self.target_rmsd = float(target_rmsd)
        self.min_rmsd_improve = float(min_rmsd_improve)

        self.progressive_nmax = bool(progressive_nmax)
        self.nmax_frac0 = float(nmax_frac0)
        self.nmax_eta = float(nmax_eta)
        self.risk_topk = int(risk_topk)

        self.mask = mask
        self.maskA = maskA
        self.maskB = maskB
        self.adj_exempt_A = adj_exempt_A
        self.adj_exempt_B = adj_exempt_B

        self.enm_kwargs = dict(enm_kwargs)

        (
            self.coordsA,
            self.coordsB,
            self.title,
            self.atoms,
            self.weights,
            self.maskA,
            self.maskB,
            rmsd0,
        ) = self._check_input(a, b, aligned=self.aligned)

        self.rmsd0 = rmsd0
        self.rmsds = [rmsd0]
        self.defvecsA = []
        self.defvecsB = []

        # adjacency exemption default building if not provided
        if self.adj_exempt is None:
            self.adj_exempt = None
        if self.adj_exempt_A is None:
            self.adj_exempt_A = self._build_adj_exempt(
                self.coordsA, adj_max=self.adj_max, tol=self.adj_tol
            )
        if self.adj_exempt_B is None:
            self.adj_exempt_B = self._build_adj_exempt(
                self.coordsB, adj_max=self.adj_max, tol=self.adj_tol
            )

        # histories
        self.f_hist = [float(self.f)]
        self.f_histA = [float(self.f)]
        self.f_histB = [float(self.f)]
        self.min_slack_hist = []
        self.min_slack_histA = []
        self.min_slack_histB = []
        self.worst_edges_hist = []
        self.worst_edges_histA = []
        self.worst_edges_histB = []

    # ---------- public API ----------

    def run(self, n_steps, mode=AANM_DEFAULT):
        if not isinstance(n_steps, Integral):
            raise TypeError("n_steps must be an integer")

        if mode == AANM_ONEWAY:
            return self._run_oneway(n_steps)
        elif mode == AANM_ALTERNATING:
            return self._run_alternating(n_steps)
        elif mode == AANM_BOTHWAYS:
            return self._run_bothways(n_steps)
        else:
            raise ValueError(f"unknown aANM mode: {mode}")

    # ---------- core drivers ----------

    def _run_oneway(self, n_steps):
        n_modes = self.n_modes
        coordsA = self.coordsA.copy()
        coordsB = self.coordsB.copy()

        LOGGER.timeit("_prody_calcAdaptiveANM")

        _, ms0, worst0 = self._compute_slack_info(
            coordsA,
            adj_max=self.adj_max,
            risk_topk=self.risk_topk,
            adj_exempt=self.adj_exempt or self.adj_exempt_A,
        )
        self.min_slack_hist.append(ms0)
        self.worst_edges_hist.append(worst0.tolist())

        ensemble = Ensemble(self.title + "_aANM")
        ensemble.setAtoms(self.atoms)
        ensemble.setCoords(coordsB)
        ensemble.setWeights(self.weights)
        ensemble.addCoordset(coordsA.copy())

        stall_counter = 0
        fmin_offset = 0.0
        terminated_on_fmin_cap = False
        n = 0
        attempt_idx = 0

        while attempt_idx < n_steps:
            rmsd_current = self.rmsds[-1]
            base_Fmin = self._compute_Fmin_from_RMSD(
                rmsd_current, self.rmsd0, self.Fmin_min, self.Fmin_max
            )
            Fmin_step = min(base_Fmin + fmin_offset, self.Fmin_max)

            n_modes, rmsd_improvement = self._step(
                coordsA,
                coordsB,
                n_modes,
                ensemble,
                self.defvecsA,
                self.rmsds,
                Fmin=Fmin_step,
                mask=self.maskA or self.mask,
                adj_exempt=self.adj_exempt or self.adj_exempt_A,
                f_hist=self.f_hist,
                min_slack_hist=self.min_slack_hist,
                worst_edges_hist=self.worst_edges_hist,
            )

            attempt_idx += 1

            if rmsd_improvement == 0.0 and n_modes > 0:
                if self.dynamic_fmin:
                    stall_counter += 1
                    if stall_counter >= self.stall_steps:
                        max_offset = max(0.0, self.Fmin_max - self.Fmin_min)
                        old_offset = fmin_offset
                        fmin_offset = min(
                            fmin_offset + self.Fmin_bump_increment, max_offset
                        )
                        if fmin_offset <= old_offset:
                            LOGGER.info(
                                "Fmin offset has reached its maximum "
                                f"({fmin_offset:.3f}); ending simulation gracefully."
                            )
                            terminated_on_fmin_cap = True
                            break
                        LOGGER.warn(
                            f"Stalled for {self.stall_steps} zero-improvement attempts. "
                            f"Increasing Fmin offset to {fmin_offset:.3f}."
                        )
                        stall_counter = 0
                continue

            if n_modes == 0:
                LOGGER.report(
                    "One-way Adaptive ANM converged in %.2fs.",
                    "_prody_calcAdaptiveANM",
                )
                break

            LOGGER.info("")
            LOGGER.info(f"Starting cycle {n + 1} with initial structure {self.title}")

            if self.dynamic_fmin:
                if rmsd_improvement < self.min_rmsd_improve and n_modes > 0:
                    stall_counter += 1
                    LOGGER.info(
                        "Stall detected (improvement=%.2e). Counter: %d/%d"
                        % (rmsd_improvement, stall_counter, self.stall_steps)
                    )
                else:
                    stall_counter = 0

                if stall_counter >= self.stall_steps:
                    max_offset = max(0.0, self.Fmin_max - self.Fmin_min)
                    old_offset = fmin_offset
                    fmin_offset = min(
                        fmin_offset + self.Fmin_bump_increment, max_offset
                    )
                    if fmin_offset <= old_offset:
                        LOGGER.info(
                            "Fmin offset has reached its maximum "
                            f"({fmin_offset:.3f}); ending simulation gracefully."
                        )
                        terminated_on_fmin_cap = True
                        break
                    LOGGER.warn(
                        f"Stalled for {self.stall_steps} small-improvement cycles. "
                        f"Increasing Fmin offset additively to {fmin_offset:.3f}."
                    )
                    stall_counter = 0
            else:
                if rmsd_improvement < self.min_rmsd_improve and n_modes > 0:
                    LOGGER.warn(
                        "RMSD improvement %.2e < min_rmsd_improve=%.2e. Stopping."
                        % (rmsd_improvement, self.min_rmsd_improve)
                    )
                    n += 1
                    break

            n += 1

        if terminated_on_fmin_cap:
            LOGGER.info(
                "One-way Adaptive ANM terminated because Fmin offset reached its "
                f"maximum value ({fmin_offset:.3f})."
            )

        LOGGER.report(
            "One-way Adaptive ANM finished in %.2fs.", "_prody_calcAdaptiveANM"
        )
        return ensemble

    def _run_alternating(self, n_steps):
        n_modes = self.n_modes
        coordsA = self.coordsA.copy()
        coordsB = self.coordsB.copy()

        LOGGER.timeit("_prody_calcAdaptiveANM")

        _, msA0, worstA0 = self._compute_slack_info(
            coordsA,
            adj_max=self.adj_max,
            risk_topk=self.risk_topk,
            adj_exempt=self.adj_exempt_A,
        )
        self.min_slack_histA.append(msA0)
        self.worst_edges_histA.append(worstA0.tolist())

        _, msB0, worstB0 = self._compute_slack_info(
            coordsB,
            adj_max=self.adj_max,
            risk_topk=self.risk_topk,
            adj_exempt=self.adj_exempt_B,
        )
        self.min_slack_histB.append(msB0)
        self.worst_edges_histB.append(worstB0.tolist())

        ensA = Ensemble(self._get_title(self.a, "structureA") + "_aANM")
        ensA.setAtoms(self.atoms)
        ensA.setCoords(coordsB)
        ensA.setWeights(self.weights)
        ensA.addCoordset(coordsA.copy())

        ensB = Ensemble(self._get_title(self.b, "structureB") + "_aANM")
        ensB.setAtoms(self.atoms)
        ensB.setCoords(coordsA)
        ensB.setWeights(self.weights)
        ensB.addCoordset(coordsB.copy())

        stateA = {"stall_counter": 0, "fmin_offset": 0.0}
        stateB = {"stall_counter": 0, "fmin_offset": 0.0}
        terminated_on_fmin_cap = False

        n = 0
        while n < n_steps:
            # A -> B
            rmsd_current = self.rmsds[-1]
            base_FminA = self._compute_Fmin_from_RMSD(
                rmsd_current, self.rmsd0, self.Fmin_min, self.Fmin_max
            )
            Fmin_stepA = min(base_FminA + stateA["fmin_offset"], self.Fmin_max)

            LOGGER.info("")
            LOGGER.info(
                "Starting cycle %d with initial structure %s"
                % (n + 1, self._get_title(self.a, "structure A"))
            )

            n_modes_A, rmsd_improvement_A = self._step(
                coordsA,
                coordsB,
                n_modes,
                ensA,
                self.defvecsA,
                self.rmsds,
                Fmin=Fmin_stepA,
                mask=self.maskA or self.mask,
                adj_exempt=self.adj_exempt_A,
                f_hist=self.f_histA,
                min_slack_hist=self.min_slack_histA,
                worst_edges_hist=self.worst_edges_histA,
            )

            if self.dynamic_fmin:
                if rmsd_improvement_A < self.min_rmsd_improve and n_modes_A > 0:
                    stateA["stall_counter"] += 1
                else:
                    stateA["stall_counter"] = 0

                if stateA["stall_counter"] >= self.stall_steps:
                    max_offset = max(0.0, self.Fmin_max - self.Fmin_min)
                    old_offset = stateA["fmin_offset"]
                    stateA["fmin_offset"] = min(
                        stateA["fmin_offset"] + self.Fmin_bump_increment, max_offset
                    )
                    if stateA["fmin_offset"] <= old_offset:
                        LOGGER.info(
                            "A->B Fmin offset has reached its maximum (%.3f); "
                            "ending simulation gracefully."
                            % stateA["fmin_offset"]
                        )
                        terminated_on_fmin_cap = True
                        break
                    LOGGER.warn(
                        "A->B stalled. Increasing Fmin offset additively to %.3f."
                        % stateA["fmin_offset"]
                    )
                    stateA["stall_counter"] = 0
            else:
                if rmsd_improvement_A < self.min_rmsd_improve and n_modes_A > 0:
                    LOGGER.warn(
                        "A->B ΔRMSD=%.2e < min_rmsd_improve=%.2e. Stopping."
                        % (rmsd_improvement_A, self.min_rmsd_improve)
                    )
                    break

            if n_modes_A == 0:
                LOGGER.report(
                    "Alternating Adaptive ANM converged in %.2fs.",
                    "_prody_calcAdaptiveANM",
                )
                break

            # B -> A
            rmsd_current = self.rmsds[-1]
            base_FminB = self._compute_Fmin_from_RMSD(
                rmsd_current, self.rmsd0, self.Fmin_min, self.Fmin_max
            )
            Fmin_stepB = min(base_FminB + stateB["fmin_offset"], self.Fmin_max)

            LOGGER.info("")
            LOGGER.info(
                "Starting cycle %d with initial structure %s"
                % (n + 1, self._get_title(self.b, "structure B"))
            )

            n_modes_B, rmsd_improvement_B = self._step(
                coordsB,
                coordsA,
                n_modes,
                ensB,
                self.defvecsB,
                self.rmsds,
                Fmin=Fmin_stepB,
                mask=self.maskB or self.mask,
                adj_exempt=self.adj_exempt_B,
                f_hist=self.f_histB,
                min_slack_hist=self.min_slack_histB,
                worst_edges_hist=self.worst_edges_histB,
            )

            if self.dynamic_fmin:
                if rmsd_improvement_B < self.min_rmsd_improve and n_modes_B > 0:
                    stateB["stall_counter"] += 1
                else:
                    stateB["stall_counter"] = 0

                if stateB["stall_counter"] >= self.stall_steps:
                    max_offset = max(0.0, self.Fmin_max - self.Fmin_min)
                    old_offset = stateB["fmin_offset"]
                    stateB["fmin_offset"] = min(
                        stateB["fmin_offset"] + self.Fmin_bump_increment, max_offset
                    )
                    if stateB["fmin_offset"] <= old_offset:
                        LOGGER.info(
                            "B->A Fmin offset has reached its maximum (%.3f); "
                            "ending simulation gracefully."
                            % stateB["fmin_offset"]
                        )
                        terminated_on_fmin_cap = True
                        break
                    LOGGER.warn(
                        "B->A stalled. Increasing Fmin offset additively to %.3f."
                        % stateB["fmin_offset"]
                    )
                    stateB["stall_counter"] = 0
            else:
                if rmsd_improvement_B < self.min_rmsd_improve and n_modes_B > 0:
                    LOGGER.warn(
                        "B->A ΔRMSD=%.2e < min_rmsd_improve=%.2e. Stopping."
                        % (rmsd_improvement_B, self.min_rmsd_improve)
                    )
                    break

            n += 1
            n_modes = n_modes_B
            if n_modes == 0:
                LOGGER.report(
                    "Alternating Adaptive ANM converged in %.2fs.",
                    "_prody_calcAdaptiveANM",
                )
                break

        if terminated_on_fmin_cap:
            LOGGER.info(
                "Alternating Adaptive ANM terminated because Fmin offset reached "
                "its maximum value (A: %.3f, B: %.3f)."
                % (stateA["fmin_offset"], stateB["fmin_offset"])
            )

        LOGGER.report(
            "Alternating Adaptive ANM finished in %.2fs.", "_prody_calcAdaptiveANM"
        )
        ensemble = ensA + ensB[::-1]
        return ensemble

    def _run_bothways(self, n_steps):
        LOGGER.info("Running forward pass: A -> B")
        forward = self._run_oneway(n_steps)

        LOGGER.info("Running reverse pass: B -> A")
        reverse_runner = AdaptiveANM(
            self.b,
            self.a,
            aligned=True,
            n_modes=self.n_modes,
            model=self.model,
            n_max_modes=self.n_max_modes,
            f=self.f,
            f_auto=self.f_auto,
            f_max=self.f_max,
            f_min_sched=self.f_min_sched,
            f_gamma=self.f_gamma,
            Fmin=self.Fmin,
            Fmin_min=self.Fmin_min,
            Fmin_max=self.Fmin_max,
            dynamic_fmin=self.dynamic_fmin,
            Fmin_bump_increment=self.Fmin_bump_increment,
            stall_steps=self.stall_steps,
            adj_max=self.adj_max,
            adj_tol=self.adj_tol,
            cutoff=self.cutoff,
            max_backtracks=self.max_backtracks,
            min_f_eff=self.min_f_eff,
            min_rmsd_step=self.min_rmsd_step,
            min_rmsd_diff=self.min_rmsd_diff,
            target_rmsd=self.target_rmsd,
            min_rmsd_improve=self.min_rmsd_improve,
            progressive_nmax=self.progressive_nmax,
            nmax_frac0=self.nmax_frac0,
            nmax_eta=self.nmax_eta,
            risk_topk=self.risk_topk,
            mask=self.mask,
            **self.enm_kwargs,
        )
        reverse = reverse_runner._run_oneway(n_steps)

        full = Ensemble(forward.getTitle() + "_bothways")
        full.setCoords(forward.getCoords())
        full.setAtoms(forward.getAtoms())
        full.setWeights(forward.getWeights())
        for i in range(forward.numCoordsets()):
            full.addCoordset(forward.getCoordsets()[i])
        for i in range(1, reverse.numCoordsets()):
            full.addCoordset(reverse.getCoordsets()[-i])
        LOGGER.info("Both-way Adaptive ANM completed successfully.")
        return full

    # ---------- single step (was calcStep) ----------

    def _step(
        self,
        initial,
        target,
        n_modes,
        ensemble,
        defvecs,
        rmsds,
        *,
        Fmin,
        mask=None,
        adj_exempt=None,
        f_hist,
        min_slack_hist,
        worst_edges_hist,
    ):
        f_global = float(self.f)
        Fmin_max = self.Fmin_max
        cutoff = self.cutoff
        adj_max = self.adj_max
        max_backtracks = self.max_backtracks
        min_f_eff = self.min_f_eff
        min_rmsd_step = self.min_rmsd_step

        weights = ensemble.getWeights()
        if weights is not None:
            weights = weights.flatten()

        coords_init = initial
        coords_tar = target
        dof = coords_init.size - 6
        if dof <= 0:
            LOGGER.warn("Not enough DOF; returning without update.")
            return n_modes, 0.0

        raw_n_max_modes = self.n_max_modes if self.n_max_modes is not None else dof
        if isinstance(raw_n_max_modes, float) and 0 < raw_n_max_modes <= 1:
            n_max_modes = max(1, int(raw_n_max_modes * dof))
        else:
            n_max_modes = min(dof, int(raw_n_max_modes))

        # defvec
        defvec = coords_tar - coords_init
        d = defvec.flatten()
        if weights is not None:
            d *= weights.repeat(3)
        defvecs.append(d)

        if self.f_auto:
            f_global = self._compute_adaptive_f(
                defvecs, f_max=self.f_max, f_min=self.f_min_sched, gamma=self.f_gamma
            )
            self.f = f_global

        if self.progressive_nmax:
            eff_n_max_modes = self._compute_progressive_nmax(
                defvecs,
                n_max_modes,
                frac0=self.nmax_frac0,
                eta=self.nmax_eta,
            )
        else:
            eff_n_max_modes = n_max_modes
        eff_n_max_modes = int(max(1, min(n_max_modes, eff_n_max_modes)))

        n_modes = min(max(1, int(n_modes)), eff_n_max_modes)
        model = self.model

        Fmin = float(min(Fmin, Fmin_max))

        enm_kwargs = dict(self.enm_kwargs)
        anm_h, _ = calcENM(
            coords_init,
            select=mask,
            mask=mask,
            model=model,
            trim="trim",
            n_modes=eff_n_max_modes,
            **enm_kwargs,
        )
        if mask is not None:
            anm_h.masked = False
        all_eigvecs = anm_h.getEigvecs()[:, :eff_n_max_modes]

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

        rmsd_before = rmsds[-1]

        accepted = False
        coords_updated = coords_init.copy()
        accepted_f_eff = 0.0
        used_h = 0
        best_improvement = 0.0
        best_cand = None
        best_f_eff = None
        best_k = None

        for k in range(int(max_backtracks) + 1):
            f_eff = float(f_global) / (2**k)
            if f_eff < min_f_eff:
                break
            s_step = f_eff * s_base
            cand = coords_init + s_step * v3
            if self._check_disconnection(
                cand, cutoff, adj_max=adj_max, adj_exempt=adj_exempt
            ):
                continue
            trial_rmsd = calcRMSD(cand, coords_tar, weights)
            improvement = rmsd_before - trial_rmsd
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
                    "(min_rmsd_step=%.2e) to avoid stall."
                    % (best_improvement, min_rmsd_step)
                )
            accepted = True
            coords_updated = best_cand
            accepted_f_eff = best_f_eff
            used_h = best_k

        if not accepted:
            LOGGER.warn("All step attempts failed; exiting step with no update.")
            return n_modes, 0.0

        mid = 0.5 * (coords_init + coords_updated)
        ensemble.addCoordset(mid.copy())
        ensemble.addCoordset(coords_updated.copy())

        f_hist.append(float(accepted_f_eff))
        f_hist.append(float(accepted_f_eff))

        s_mid, ms_mid, worst_mid = self._compute_slack_info(
            mid,
            adj_max=self.adj_max,
            risk_topk=self.risk_topk,
            adj_exempt=adj_exempt,
        )
        min_slack_hist.append(ms_mid)
        worst_edges_hist.append(worst_mid.tolist())

        s_end, ms_end, worst_end = self._compute_slack_info(
            coords_updated,
            adj_max=self.adj_max,
            risk_topk=self.risk_topk,
            adj_exempt=adj_exempt,
        )
        min_slack_hist.append(ms_end)
        worst_edges_hist.append(worst_end.tolist())

        initial[:] = coords_updated
        rmsd_after = calcRMSD(initial, coords_tar, weights)
        rmsds.append(rmsd_after)
        rmsd_improvement = rmsd_before - rmsd_after

        h = int(used_h)
        scale_str = f"1/2^{h}"
        highest_mode = int(sel_idx.max()) + 1 if sel_idx.size > 0 else 0

        LOGGER.info(
            "Step successful [h=%d, scale=%s] (Fmin=%.3f, f_eff=%.3f, "
            "highest_mode=%d, ΔRMSD=%.4g)"
            % (h, scale_str, Fmin, accepted_f_eff, highest_mode, rmsd_improvement)
        )

        if n_max_modes and n_max_modes > 0:
            nmax_frac = float(eff_n_max_modes) / float(n_max_modes)
        else:
            nmax_frac = 1.0
        LOGGER.info(
            "Step summary (Fmin=%.3f, n_used=%d, nmax_frac=%.3f)"
            % (Fmin, n_used, nmax_frac)
        )
        LOGGER.info("Current RMSD is %.6f\n" % rmsd_after)

        if self._check_convergence(rmsds, initial, adj_exempt=adj_exempt):
            n_modes_out = 0
        else:
            n_modes_out = max(1, min(int(n_modes), n_used))
        return n_modes_out, rmsd_improvement

    # ---------- helper methods ----------

    @staticmethod
    def _sequential_edge_lengths(coords):
        r = np.asarray(coords, float)
        if r.shape[0] < 2:
            return np.zeros(0, dtype=float)
        diffs = r[1:] - r[:-1]
        return np.sqrt(np.sum(diffs * diffs, axis=1))

    @classmethod
    def _build_adj_exempt(cls, native_coords, adj_max=4.5, tol=0.0):
        elen_native = cls._sequential_edge_lengths(native_coords)
        if elen_native.size == 0:
            return None
        threshold = float(adj_max) + float(tol)
        return elen_native > threshold

    @classmethod
    def _compute_slack_info(cls, coords, adj_max=4.5, risk_topk=8, adj_exempt=None):
        elen = cls._sequential_edge_lengths(coords)
        slack = float(adj_max) - elen
        if slack.size == 0:
            return slack, np.inf, np.array([], dtype=int)
        if adj_exempt is not None and adj_exempt.shape == slack.shape:
            slack = slack.copy()
            slack[adj_exempt] = np.inf
        order = np.argsort(slack)
        worst = order[: min(risk_topk, order.size)]
        return slack, float(slack.min()), worst

    @classmethod
    def _check_disconnection(cls, coords, cutoff, adj_max=4.5, adj_exempt=None):
        r = np.asarray(coords, float)
        N = r.shape[0]
        if N < 2:
            return False
        all_dists = np.sqrt(((r[:, None, :] - r[None, :, :]) ** 2).sum(axis=2))
        nn = np.empty(N, dtype=float)
        for i in range(N):
            row = all_dists[i]
            nn[i] = min(
                np.min(row[:i]) if i > 0 else np.inf,
                np.min(row[i + 1 :]) if i + 1 < N else np.inf,
            )
        isolated = nn > float(cutoff)

        broken_adj = False
        if N >= 2:
            edge_d = cls._sequential_edge_lengths(coords)
            if adj_exempt is not None and adj_exempt.shape == edge_d.shape:
                broken_adj = np.any((edge_d > float(adj_max)) & (~adj_exempt))
            else:
                broken_adj = np.any(edge_d > float(adj_max))
        return bool(np.any(isolated) or broken_adj)

    @staticmethod
    def _compute_Fmin_from_RMSD(rmsd_current, rmsd0, Fmin_min=0.4, Fmin_max=0.9):
        if rmsd0 <= 1e-8:
            progress = 1.0
        else:
            progress = 1.0 - float(rmsd_current / rmsd0)
        progress = float(np.clip(progress, 0.0, 1.0))
        Fmin = float(Fmin_min + (Fmin_max - Fmin_min) * progress)
        return float(min(Fmin, Fmin_max))

    @staticmethod
    def _compute_adaptive_f(defvecs, f_max=1.0, f_min=0.2, gamma=1.0):
        if not defvecs:
            return float(f_max)
        d0, dk = defvecs[0], defvecs[-1]
        n0, nk = norm(d0), norm(dk)
        if n0 <= 1e-12:
            return float(f_min)
        progress = 1.0 - np.sqrt(max(nk, 0.0) / max(n0, 1e-12))
        progress = float(np.clip(progress, 0.0, 1.0))
        return float(f_min + (f_max - f_min) * (1.0 - (progress ** float(gamma))))

    @staticmethod
    def _compute_progressive_nmax(defvecs, n_max_modes, frac0=0.2, eta=1.5):
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
            progress = 1.0 - np.sqrt(max(nk, 0.0) / max(n0, 1.0e-12))
        progress = float(np.clip(progress, 0.0, 1.0))
        frac = float(frac0) + (1.0 - float(frac0)) * (progress ** float(eta))
        eff_n = int(np.floor(n_max_modes * frac))
        return int(max(1, min(n_max_modes, eff_n)))

    def _check_convergence(self, rmsds, coords, adj_exempt=None):
        if len(rmsds) > 4 and np.all(
            np.abs(np.diff(rmsds[-4:])) < self.min_rmsd_diff
        ):
            LOGGER.warn(f"RMSD change fell below {self.min_rmsd_diff}")
            return True
        if rmsds[-1] < self.target_rmsd:
            LOGGER.info(
                "Convergence reached: RMSD %.3f is below target %.3f."
                % (rmsds[-1], self.target_rmsd)
            )
            return True
        if self._check_disconnection(
            coords,
            self.cutoff,
            adj_max=self.adj_max,
            adj_exempt=adj_exempt if adj_exempt is not None else self.adj_exempt,
        ):
            LOGGER.warn("Disconnection found in the structure.")
            return True
        return False

    @staticmethod
    def _get_title(structure, def_title="structure"):
        return structure.getTitle() if isinstance(structure, Atomic) else def_title

    @staticmethod
    def _check_input(a, b, **kwargs):
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
                title = "Unknown"
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

        if not kwargs.get("aligned", False):
            coordsA, _ = superpose(coordsA, coordsB, weights)
        rmsd = calcRMSD(coordsA, coordsB, weights)
        LOGGER.info("Initialized Adaptive ANM with RMSD {:4.3f}\n".format(rmsd))
        return coordsA, coordsB, title, atoms, weights, maskA, maskB, rmsd

