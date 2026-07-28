# -*- coding: utf-8 -*-
"""Additively-weighted (Apollonius) Voronoi tessellation for :func:`.calcChannels`
``diagram="weighted"``, built on the third-party ``vorpy`` package.

vorpy computes the *true* additively-weighted (a.k.a. Apollonius) Voronoi network,
which bakes van der Waals radii into the diagram exactly instead of approximating it
by homogenising atoms into uniform balls. This module:

1. Ships a compiled (numba) drop-in for vorpy's ``calc_vert`` (the 4-sphere Apollonius
   vertex solve). vorpy's own ``calc_vert`` solves a quadratic with ``numpy.roots``
   (a companion-matrix eigenvalue decomposition) and wraps every call in numpy
   ``seterr`` context managers, costing ~116 us/call and ~80% of the total run time.
   The compiled version uses the closed-form quadratic root and is bit-equivalent on
   the common (non-degenerate) path, ~37x faster, deferring to the original only on
   the rare singular (``F ~= 0``) branch.
2. Drives ``vorpy.Network`` in additively-weighted mode and adapts its vertex network
   into the ``(simplices, neighbors, vertices)`` arrays the channel finder consumes,
   plus an exact per-vertex clearance.

The heavy imports (``vorpy``, ``numba``) happen at import time of *this* module, which
is imported lazily from the ``diagram="weighted"`` branch of :func:`.calcChannels`, so
they are not a hard dependency of ProDy or of the default channel path.
"""

import os
import io
import re
import sys
import time
import hashlib
import contextlib
import tempfile
from math import sqrt

import numpy as np

__all__ = ['buildAwTessellation', 'accelerateVorpy', 'resolveCachePath']

# On-disk cache of the raw additively-weighted network. Building it with vorpy is the
# expensive step (minutes on a full structure); the downstream channel search is cheap
# and often re-run with different r2/bottleneck/sparsity/start_point, none of which
# change the diagram. Persisting the (simplices, neighbors, vertices, clearances)
# arrays lets those re-runs -- and repeated debugging passes -- skip vorpy entirely.
_CACHE_VERSION = 1   # bump whenever the cached array layout/semantics change

_EPS_F = 1e-12       # |F| below this => degenerate spatial matrix, defer to vorpy
_TOL_IMAG = 1e-12    # matches vorpy _real_roots_quadratic imaginary tolerance


def _importNumbaKernel():
    """Build and return the numba-compiled AW vertex core. Kept in a function so the
    numba import/compile only happens when the weighted path is actually used."""
    from numba import njit

    @njit(cache=True, fastmath=False)
    def _awCore(b0x, b0y, b0z, b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z,
                r0, r1, r2, r3):
        # The four ball centres as 12 scalars followed by the four radii. Taking
        # scalars (not (4,3)/(4,) arrays) lets calcVertFast skip building a numpy
        # array per call -- np.ascontiguousarray was ~3 s of the accelerated run.
        # Returns (degenerate, code, ax,ay,az,ar, bx,by,bz,br).
        #   degenerate=1 -> caller falls back to vorpy's original calc_vert (F ~= 0).
        #   code: 0 no vertex, 1 loc only, 2 loc+loc2 (ordered by |R| ascending).
        l0x, l0y, l0z = b0x, b0y, b0z
        l1x, l1y, l1z = b1x - b0x, b1y - b0y, b1z - b0z
        l2x, l2y, l2z = b2x - b0x, b2y - b0y, b2z - b0z
        l3x, l3y, l3z = b3x - b0x, b3y - b0y, b3z - b0z
        r0_2 = r0 * r0

        a1, b1, c1, d1 = 2*l1x, 2*l1y, 2*l1z, 2*(r1 - r0)
        f1 = r0_2 - r1*r1 + l1x*l1x + l1y*l1y + l1z*l1z
        a2, b2, c2, d2 = 2*l2x, 2*l2y, 2*l2z, 2*(r2 - r0)
        f2 = r0_2 - r2*r2 + l2x*l2x + l2y*l2y + l2z*l2z
        a3, b3, c3, d3 = 2*l3x, 2*l3y, 2*l3z, 2*(r3 - r0)
        f3 = r0_2 - r3*r3 + l3x*l3x + l3y*l3y + l3z*l3z

        F = a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1
        if F > -_EPS_F and F < _EPS_F:
            return 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        F_2 = F * F
        F10 =  b1*c2*f3 - b1*c3*f2 - b2*c1*f3 + b2*c3*f1 + b3*c1*f2 - b3*c2*f1
        F11 = -b1*c2*d3 + b1*c3*d2 + b2*c1*d3 - b2*c3*d1 - b3*c1*d2 + b3*c2*d1
        F20 = -a1*c2*f3 + a1*c3*f2 + a2*c1*f3 - a2*c3*f1 - a3*c1*f2 + a3*c2*f1
        F21 =  a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1
        F30 =  a1*b2*f3 - a1*b3*f2 - a2*b1*f3 + a2*b3*f1 + a3*b1*f2 - a3*b2*f1
        F31 = -a1*b2*d3 + a1*b3*d2 + a2*b1*d3 - a2*b3*d1 - a3*b1*d2 + a3*b2*d1

        qa = (F11*F11 + F21*F21 + F31*F31) - F_2
        qb = 2.0 * ((F10*F11 + F20*F21 + F30*F31) - r0*F_2)
        qc = (F10*F10 + F20*F20 + F30*F30) - r0_2*F_2

        R_a = 0.0; R_b = 0.0; nR = 0
        if qa > -1e-300 and qa < 1e-300:
            if not (qb > -1e-300 and qb < 1e-300):
                R_a = -qc / qb; nR = 1
        else:
            disc = qb*qb - 4.0*qa*qc
            if disc >= 0.0:
                sq = np.sqrt(disc)
                R_a = (-qb + sq) / (2.0*qa)
                R_b = (-qb - sq) / (2.0*qa)
                nR = 2
            elif disc > -_TOL_IMAG:
                R_a = -qb / (2.0*qa); nR = 1

        if nR == 0:
            return 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        if nR == 1:
            x = (F10 + R_a*F11)/F + l0x
            y = (F20 + R_a*F21)/F + l0y
            z = (F30 + R_a*F31)/F + l0z
            return 0, 1, x, y, z, R_a, 0.0, 0.0, 0.0, 0.0

        if abs(R_a) > abs(R_b):
            R_a, R_b = R_b, R_a
        x0 = (F10 + R_a*F11)/F + l0x
        y0 = (F20 + R_a*F21)/F + l0y
        z0 = (F30 + R_a*F31)/F + l0z
        x1 = (F10 + R_b*F11)/F + l0x
        y1 = (F20 + R_b*F21)/F + l0y
        z1 = (F30 + R_b*F31)/F + l0z
        max_ball_rad = max(max(r0, r1), max(r2, r3))

        if R_a < 0.0 or R_b < 0.0:
            if R_a > 0.0 or abs(R_a) < max_ball_rad:
                if R_b > 0.0 or abs(R_b) < max_ball_rad:
                    return 0, 2, x0, y0, z0, R_a, x1, y1, z1, R_b
                return 0, 1, x0, y0, z0, R_a, 0.0, 0.0, 0.0, 0.0
            elif R_b > 0.0 or abs(R_b) < max_ball_rad:
                return 0, 1, x1, y1, z1, R_b, 0.0, 0.0, 0.0, 0.0
            return 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        return 0, 2, x0, y0, z0, R_a, x1, y1, z1, R_b

    return _awCore


_ACCEL_APPLIED = False


def _patchEverywhere(name, orig, replacement):
    """Rebind ``name`` to ``replacement`` in every already-imported vorpy module that
    currently holds ``orig`` (vorpy does ``from ... import name`` widely, so each such
    binding must be replaced). Returns the number of bindings swapped."""
    n = 0
    for modname, mod in list(sys.modules.items()):
        if mod is None or not modname.startswith('vorpy'):
            continue
        if getattr(mod, name, None) is orig:
            setattr(mod, name, replacement)
            n += 1
    return n


def accelerateVorpy():
    """Monkeypatch vorpy's two pure-Python hot functions with faster drop-ins,
    everywhere they are bound: ``calc_vert`` (the 4-sphere Apollonius solve, replaced
    by the compiled closed-form kernel) and ``calc_dist`` (a 3-D Euclidean distance
    that vorpy implements with per-call numpy array allocation + a Python ``sum``,
    called millions of times). Idempotent. Returns the number of bindings replaced."""
    global _ACCEL_APPLIED
    if _ACCEL_APPLIED:
        return 0

    # calc_vert -> compiled closed-form kernel
    import vorpy.src.calculations.vert as _vmod
    vert_orig = _vmod.calc_vert
    core = _importNumbaKernel()

    def calcVertFast(locs, rads):
        # locs is a 4-list of ball centres, rads a 4-list of radii (vorpy always
        # calls with [locs[i] for i in balls]). Hand the kernel their 16 scalars
        # directly instead of materialising a numpy array per call.
        l0, l1, l2, l3 = locs
        deg, code, ax, ay, az, ar, bx, by, bz, br = core(
            l0[0], l0[1], l0[2], l1[0], l1[1], l1[2],
            l2[0], l2[1], l2[2], l3[0], l3[1], l3[2],
            rads[0], rads[1], rads[2], rads[3])
        if deg:
            return vert_orig(locs, rads)
        if code == 0:
            return None, None, None, None
        if code == 1:
            return [ax, ay, az], ar, None, None
        return [ax, ay, az], ar, [bx, by, bz], br

    # calc_dist -> allocation-free 3-D distance (falls back for non-3-D inputs)
    import vorpy.src.calculations.calcs as _cmod
    dist_orig = _cmod.calc_dist

    def calcDistFast(l0, l1):
        if len(l0) == 3:
            d0 = l0[0] - l1[0]
            d1 = l0[1] - l1[1]
            d2 = l0[2] - l1[2]
            return sqrt(d0 * d0 + d1 * d1 + d2 * d2)
        return dist_orig(l0, l1)

    # write_verts -> no-op: vorpy dumps every vertex to a .txt at the end of
    # find_net_verts via pandas iterrows (slow over tens of thousands of rows); we
    # build our own arrays from net.verts and never read that file.
    import vorpy.src.output as _omod
    wv_orig = _omod.write_verts

    def writeVertsNoop(net):
        return None

    n = _patchEverywhere('calc_vert', vert_orig, calcVertFast)
    n += _patchEverywhere('calc_dist', dist_orig, calcDistFast)
    n += _patchEverywhere('write_verts', wv_orig, writeVertsNoop)
    _ACCEL_APPLIED = True
    return n


# vorpy prints a "\r...finding vertices: <n> verts - <pct> %" line to stdout as it
# traces the network; we parse it to drive a ProDy progress bar and otherwise swallow
# vorpy's chatter.
_PROG_RE = re.compile(r'finding vertices:\s*(\d+)\s*verts\s*-\s*([\d.]+)\s*%')


class _ProgressSink(io.StringIO):
    """stdout replacement that discards vorpy's output but extracts its vertex-finding
    percentage and forwards it to ``on_progress(n_verts, percent)``."""

    def __init__(self, on_progress):
        super().__init__()
        self._cb = on_progress

    def write(self, s):
        if '%' not in s:      # skip the regex on vorpy's non-progress chatter
            return len(s)
        m = _PROG_RE.search(s)
        if m is not None:
            try:
                self._cb(int(m.group(1)), float(m.group(2)))
            except Exception:
                pass
        return len(s)


@contextlib.contextmanager
def _captureVorpy(on_progress=None):
    """Contain vorpy's on-disk vertex writes in a throwaway temp directory (cwd
    restored afterwards) and capture its stdout: parsed for progress when
    ``on_progress`` is given, silently discarded otherwise."""
    cwd = os.getcwd()
    tmp = tempfile.mkdtemp(prefix='prody_awvor_')
    sink = _ProgressSink(on_progress) if on_progress is not None else io.StringIO()
    try:
        os.chdir(tmp)
        with contextlib.redirect_stdout(sink):
            yield
    finally:
        os.chdir(cwd)


def _cacheKey(coords, vdw_radii, max_vert):
    """Content hash identifying a tessellation input: the ball centres, their radii
    and the ``max_vert`` cutoff (the only inputs that change the diagram). A cache is
    reused only when all three match, so editing the structure -- or, via ``r1``,
    ``max_vert`` -- transparently forces a recompute."""
    h = hashlib.sha1()
    h.update(np.ascontiguousarray(coords, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(vdw_radii, dtype=np.float64).tobytes())
    h.update(np.float64(max_vert).tobytes())
    return 'v{0}:{1}'.format(_CACHE_VERSION, h.hexdigest())


def _loadCache(path, key):
    """Return the cached ``(simplices, neighbors, vertices, clearances)`` at ``path``
    if it exists and was built for ``key``, else ``None`` (caller recomputes). Any
    read/format error is swallowed and treated as a miss."""
    if not path or not os.path.isfile(path):
        return None
    try:
        with np.load(path, allow_pickle=False) as d:
            if str(d['key']) != key:
                return None
            return (d['simplices'], d['neighbors'], d['vertices'], d['clearances'])
    except Exception:
        return None


def _saveCache(path, key, result):
    """Persist a tessellation ``result`` under ``key`` to ``path`` (compressed npz),
    writing to a temporary file and renaming so a crashed write can't leave a
    half-written cache. Failures warn but do not abort the calculation."""
    simplices, neighbors, vertices, clearances = result
    tmp = path + '.tmp.npz'
    try:
        d = os.path.dirname(path)
        if d and not os.path.isdir(d):
            os.makedirs(d)
        np.savez_compressed(tmp, key=np.array(key), simplices=simplices,
                            neighbors=neighbors, vertices=vertices,
                            clearances=clearances)
        os.replace(tmp, path)
    except Exception as exc:
        from prody import LOGGER
        LOGGER.warn('Could not write additively-weighted Voronoi cache to {0}: {1}'
                    .format(path, exc))
        try:
            if os.path.isfile(tmp):
                os.remove(tmp)
        except OSError:
            pass


def resolveCachePath(cache, output_path=None, title=None):
    """Resolve the ``weighted_cache`` argument of :func:`.calcChannels` to a concrete
    ``.npz`` path or ``None``. ``False``/``None`` (and ``""``) disable caching; a
    string is used verbatim; ``True`` auto-derives a path next to ``output_path``
    (e.g. ``aw.pqr`` -> ``aw.awvoronoi.npz``) or, lacking an ``output_path``, from the
    structure ``title`` in the current working directory."""
    if not cache:
        return None
    if isinstance(cache, str):
        return cache
    if output_path:
        return os.path.splitext(output_path)[0] + '.awvoronoi.npz'
    title = (title or 'structure').replace(' ', '_')
    return os.path.join(os.getcwd(), title + '.awvoronoi.npz')


def buildAwTessellation(coords, vdw_radii, max_vert=8.0, accelerate=True,
                        target_atoms=500, cache=None):
    """Build the additively-weighted (Apollonius) Voronoi network of the balls
    ``(coords, vdw_radii)`` with vorpy and adapt it to the channel finder's arrays.

    Large *and spatially extended* structures are split into local boxes (each padded
    with a correctness halo), tessellated one after another and merged by global 4-ball
    signature. Splitting only engages where a box core would be comfortably larger than
    its halo (see ``_planBoxes``); a compact globular structure -- even a few thousand
    atoms -- stays on the single-pass path, since there the halo would re-cover the
    whole structure and decomposition would be a net loss.

    :arg coords: (N,3) ball centres (atom coordinates), in Angstrom.
    :arg vdw_radii: (N,) van der Waals radii aligned with ``coords``.
    :arg max_vert: probe distance / maximum empty-sphere radius. Vertices whose
        clearance exceeds this are dropped. ~8 A keeps every channel-relevant vertex
        while pruning irrelevant large voids. Default 8.0.
    :arg accelerate: apply the compiled ``calc_vert`` kernel first. Default True.
    :arg target_atoms: approximate number of balls per box core; also the threshold
        below which the structure is tessellated in a single pass. Default 500.
    :arg cache: optional path to an ``.npz`` cache file. When given and it already
        holds a diagram built for this exact ``(coords, vdw_radii, max_vert)``, it is
        loaded and vorpy is skipped; otherwise the freshly built diagram is written
        there for next time. A mismatching or unreadable file is ignored and
        overwritten. ``None`` (default) disables caching.

    :returns: ``(simplices, neighbors, vertices, clearances)`` where
        ``simplices`` is (M,4) int (the 4 tangent atoms per vertex),
        ``neighbors`` is (M,4) int (face-adjacent vertices sharing 3 atoms, -1 padded),
        ``vertices`` is (M,3) float (Voronoi-vertex positions), and
        ``clearances`` is (M,) float (exact additively-weighted clearance per vertex).
    """
    coords = np.asarray(coords, dtype=float)
    vdw_radii = np.asarray(vdw_radii, dtype=float)
    N = len(coords)
    from prody import LOGGER

    # Reuse a previously saved diagram for this exact input if one is cached.
    key = _cacheKey(coords, vdw_radii, max_vert)
    if cache:
        cached = _loadCache(cache, key)
        if cached is not None:
            LOGGER.info('Reusing additively-weighted Voronoi diagram cached at {0} '
                        '({1} vertices); skipping vorpy.'.format(cache, len(cached[2])))
            return cached

    # Each box must be padded with a correctness halo (see _buildAwDecomposed). A box
    # is only worth splitting out if its core is meaningfully larger than that halo,
    # otherwise core+halo re-tessellates almost the whole structure at a net loss. So
    # decomposition only engages for large/elongated systems; a compact globular
    # protein (even a few thousand atoms) stays on the single-pass path.
    halo = _haloWidth(max_vert, vdw_radii)
    boxes = _planBoxes(coords, target_atoms, halo) if N > target_atoms else None
    if boxes and len(boxes) > 1:
        # Only decompose if it meaningfully shrinks the per-box work. Every box is
        # padded with the halo, so for a compact (globular) structure each box's
        # core+halo re-covers most of the atoms and splitting is a net loss (extra
        # overlap work, no gain). In that case fall back to the single pass, which is
        # both faster here and shows a fine-grained progress bar.
        loads = _boxLoads(coords, boxes, halo)
        worst = max(ch for _, ch in loads)
        if worst > 0.6 * N:
            LOGGER.info('Weighted tessellation: skipping decomposition (largest of {0} '
                        'boxes would still span {1}/{2} atoms once haloed); using a '
                        'single pass.'.format(len(boxes), worst, N))
            boxes = None
    if not boxes or len(boxes) == 1:
        if accelerate:
            accelerateVorpy()
        prog = _MonoProgress()
        verts = _runVorpy(coords, vdw_radii, max_vert, accelerate, on_progress=prog)
        prog.finish()
        balls = list(verts['balls'])
        result = _assemble(balls, list(verts['loc']), list(verts['rad']),
                           list(verts['dub']) if 'dub' in verts else None)
    else:
        result = _buildAwDecomposed(coords, vdw_radii, max_vert, accelerate,
                                    boxes, halo)

    if cache:
        _saveCache(cache, key, result)
        LOGGER.info('Cached additively-weighted Voronoi diagram ({0} vertices) to {1}.'
                    .format(len(result[2]), cache))
    return result


def _haloWidth(max_vert, vdw_radii):
    """Halo padding a box needs so that every vertex it owns (clearance <= max_vert)
    has all its defining/blocking balls present: clearance + ball radius + margin."""
    return float(max_vert) + float(np.max(vdw_radii)) + 1.0


def _boxLoads(coords, boxes, halo):
    """For each box return ``(core_count, core_plus_halo_count)``: how many atoms it
    owns and how many it must actually tessellate (core plus halo)."""
    loads = []
    for lo, hi in boxes:
        core = int(np.all((coords >= lo) & (coords < hi), axis=1).sum())
        near = int(np.all((coords >= lo - halo) & (coords < hi + halo), axis=1).sum())
        loads.append((core, near))
    return loads


class _MonoProgress(object):
    """Throttled adapter turning vorpy's vertex-finding percentage into periodic ProDy
    log lines for the single-pass (mono) path. Uses plain ``LOGGER.info`` (not the
    ``progress``/``update`` bar) so it neither assumes a TTY nor mutates the logger's
    verbosity state, and so it can never suppress later messages."""

    _msg = 'Additively-weighted Voronoi diagram'
    _interval = 3.0     # seconds between log lines

    def __init__(self):
        from prody import LOGGER
        self._logger = LOGGER
        self._last_t = 0.0

    def __call__(self, n_verts, percent):
        # Throttle by wall time (vorpy calls this per edge). Report the vertex COUNT,
        # which is exact and monotonic; vorpy's percentage is an estimate that saturates
        # at 100% while work continues, so it is only a secondary hint. Time-based
        # throttling keeps the line moving even past that bogus 100%.
        now = time.perf_counter()
        if now - self._last_t >= self._interval:
            self._last_t = now
            self._logger.info('{0}: {1} vertices ({2}%)'.format(
                self._msg, n_verts, min(percent, 100)))

    def finish(self):
        # Nothing to restore: LOGGER.info leaves no state. The caller's LOGGER.report
        # ('... tessellation constructed in %.2fs') closes out the stage.
        pass


def _runVorpy(coords, vdw_radii, max_vert, accelerate, on_progress=None):
    """Run vorpy's additively-weighted vertex finder on the given balls and return
    its ``verts`` DataFrame (columns: balls, loc, rad, dub). When ``on_progress`` is
    given it is called with ``(n_verts, percent)`` as vorpy traces the network."""
    if accelerate:
        accelerateVorpy()

    from vorpy.src.network import Network

    coords = np.asarray(coords, dtype=float)
    vdw_radii = np.asarray(vdw_radii, dtype=float)
    settings = {'surf_res': 0.2, 'surf_col': 'plasma', 'surf_scheme': 'mean',
                'max_vert': float(max_vert), 'box_size': 1.5, 'net_type': 'aw',
                'build_type': 'all', 'num_splits': None, 'print_metrics': False,
                'ball_type': 'mol', 'sys_dir': os.getcwd(), 'foam_box': None,
                'atom_rad': None, 'scheme_factor': 'log'}
    locs = [row.copy() for row in coords]
    rads = [float(x) for x in vdw_radii]

    with _captureVorpy(on_progress):
        net = Network(locs=locs, rads=rads, settings=settings, sort_balls=True)
        net.find_verts()

    return net.verts


_BIG = 1.0e6    # sentinel bound: outer leaf faces extend to +/- this so that every
                # Voronoi vertex (including those outside the atom bounding box) is
                # owned by exactly one box core.
_SPLIT_HALO_FACTOR = 5.0   # only split an axis whose extent exceeds this * halo, so
                           # each child core is comfortably wider than its halo shell.


def _planBoxes(coords, target, halo):
    """Partition space into axis-aligned leaf boxes each owning <= ``target`` centres,
    by recursive median bisection along the widest data axis. An axis is only split
    when its extent exceeds ``_SPLIT_HALO_FACTOR * halo``; recursion stops otherwise,
    so compact structures collapse to a single box (and take the single-pass path).

    Returns a list of ``(lo, hi)`` half-open bound pairs (float (3,) arrays) that tile
    *all* of space: internal cut planes are finite, while any leaf face on the global
    boundary is pushed out to +/-``_BIG`` so exterior vertices are still owned. Every
    centre falls in exactly one box."""
    lo0 = coords.min(axis=0)
    hi0 = coords.max(axis=0)
    min_extent = _SPLIT_HALO_FACTOR * halo
    boxes = []

    def rec(idx, lo, hi, lo_out, hi_out):
        # Effective bounds: outer faces open out to the +/-_BIG sentinel.
        eff_lo = np.where(lo_out, -_BIG, lo)
        eff_hi = np.where(hi_out, _BIG, hi)
        extent = hi - lo
        # Widest axis that is both long enough to be worth splitting and actually
        # spans the points in this cell.
        order = np.argsort(extent)[::-1]
        axis = None
        for a in order:
            if extent[a] > min_extent:
                axis = int(a)
                break
        if len(idx) <= target or axis is None:
            boxes.append((eff_lo, eff_hi))
            return
        split = float(np.median(coords[idx, axis]))
        left = idx[coords[idx, axis] < split]
        right = idx[coords[idx, axis] >= split]
        # guard against a degenerate median that fails to divide the points
        if split <= lo[axis] or split >= hi[axis] or len(left) == 0 or len(right) == 0:
            boxes.append((eff_lo, eff_hi))
            return
        hl = hi.copy(); hl[axis] = split
        lr = lo.copy(); lr[axis] = split
        lo_out_r = lo_out.copy(); lo_out_r[axis] = False
        hi_out_l = hi_out.copy(); hi_out_l[axis] = False
        rec(left, lo, hl, lo_out, hi_out_l)
        rec(right, lr, hi, lo_out_r, hi_out)

    ones = np.ones(3, dtype=bool)
    rec(np.arange(len(coords)), lo0, hi0, ones.copy(), ones.copy())
    return boxes


def _buildAwDecomposed(coords, vdw_radii, max_vert, accelerate, boxes, halo):
    """Tessellate each box (core + halo) independently, keep only the vertices whose
    location falls in the box core, and merge them by global 4-ball signature.

    The halo width is chosen so a ball can invalidate a kept vertex (clearance <=
    max_vert) only if its centre lies within clearance + its radius <= max_vert +
    max_radius of the vertex; since the vertex sits in the core, that halo around the
    core contains every ball a core vertex depends on (+1 A margin)."""
    from prody import LOGGER

    tasks, cores, haloed = [], [], []
    for lo, hi in boxes:
        core = np.all((coords >= lo) & (coords < hi), axis=1)
        near = np.all((coords >= lo - halo) & (coords < hi + halo), axis=1)
        if not core.any():
            continue
        gidx = np.nonzero(near)[0]
        tasks.append((coords[gidx].copy(), vdw_radii[gidx].copy(), gidx,
                      lo.copy(), hi.copy(), float(max_vert), bool(accelerate)))
        cores.append(int(core.sum()))
        haloed.append(len(gidx))

    LOGGER.info('Weighted tessellation: {0} boxes, {1}-{2} atoms each ({3}-{4} with '
                'halo).'.format(len(tasks), min(cores), max(cores),
                min(haloed), max(haloed)))

    results = []
    for t in tasks:
        results.append(_awBoxTask(t))
        LOGGER.info('Additively-weighted Voronoi diagram: box {0}/{1} done'
                    .format(len(results), len(tasks)))

    # Merge box outputs, deduplicating by global 4-ball signature. Ownership-by-core-
    # location already assigns each vertex to a single box, so collisions are rare and
    # identical; keep the first seen.
    seen = {}
    balls, locs, rads = [], [], []
    for out in results:
        for sig, x, y, z, rad in out:
            if sig in seen:
                continue
            seen[sig] = True
            balls.append(sig)
            locs.append((x, y, z))
            rads.append(rad)

    return _assemble(balls, locs, rads, None)


def _awBoxTask(task):
    """Tessellate one box's balls and return the core-owned vertices as
    ``(global_signature, x, y, z, clearance)`` records."""
    coords, vdw_radii, gidx, lo, hi, max_vert, accelerate = task
    verts = _runVorpy(coords, vdw_radii, max_vert, accelerate)
    out = []
    for bl, loc, rad in zip(verts['balls'], verts['loc'], verts['rad']):
        x, y, z = float(loc[0]), float(loc[1]), float(loc[2])
        if not (lo[0] <= x < hi[0] and lo[1] <= y < hi[1] and lo[2] <= z < hi[2]):
            continue
        sig = tuple(sorted(int(gidx[int(b)]) for b in bl))
        out.append((sig, x, y, z, float(rad)))
    return out


def _assemble(balls, locs, rads, dubs=None):
    """Build ``(simplices, neighbors, vertices, clearances)`` from per-vertex lists of
    4-ball tuples, locations and additively-weighted clearances. ``balls`` indices must
    already be in the final (global) numbering."""
    from prody import LOGGER

    if dubs is not None:
        ndub = int(sum(1 for d in dubs if d != 0))
        if ndub:
            LOGGER.warn('vorpy AW network has {0} doublet vertices; their face '
                        'adjacency may be incomplete.'.format(ndub))

    M = len(balls)
    simplices = np.array([sorted(int(b) for b in bl) for bl in balls],
                         dtype=np.int64).reshape(M, 4) if M else \
        np.zeros((0, 4), dtype=np.int64)
    vertices = np.array([np.asarray(l, dtype=float) for l in locs],
                        dtype=float).reshape(M, 3) if M else np.zeros((0, 3), float)
    clearances = np.array([float(r) for r in rads], dtype=float)

    # Face adjacency: two vertices are neighbours iff they share 3 of their 4 atoms.
    # Map each triangular face (sorted 3-atom tuple) -> the vertices carrying it;
    # a face shared by exactly two vertices links them.
    face_to_verts = {}
    for i in range(M):
        b0, b1, b2, b3 = simplices[i]
        for face in ((b0, b1, b2), (b0, b1, b3), (b0, b2, b3), (b1, b2, b3)):
            face_to_verts.setdefault(face, []).append(i)

    neighbors = np.full((M, 4), -1, dtype=np.int64)
    for i in range(M):
        b0, b1, b2, b3 = simplices[i]
        for k, face in enumerate(((b0, b1, b2), (b0, b1, b3),
                                  (b0, b2, b3), (b1, b2, b3))):
            owners = face_to_verts[face]
            if len(owners) == 2:
                neighbors[i, k] = owners[0] if owners[1] == i else owners[1]
            # len 1 -> boundary face (-1); len >2 -> degenerate, left as -1

    return simplices, neighbors, vertices, clearances
