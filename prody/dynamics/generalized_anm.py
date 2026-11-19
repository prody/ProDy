# -*- coding: utf-8 -*-
"""This module defines the generalized anisotropic network model,
:class:`.genANM`, that augments the classical ANM pairwise Hessian with
virtual-bond angle and signed virtual-bond dihedral stiffness along the Calpha
trace, plus an optional Cartesian curvature regularizer.

The angle and torsion terms are assembled as positive-semidefinite rank-one
gradient outer products built from rigid-body-invariant internal coordinates, so
they do not lift the six rigid-body zero modes.  The optional curvature term is a
Cartesian second-difference penalty; it preserves translations but not rotations,
so it reduces the exact zero-mode count from six to three.

Consequences for the zero-mode count of a connected, non-degenerate structure:

  * ``kappa == 0``  ->  6 exact zero modes (3 translations + 3 rotations)
  * ``kappa != 0``  ->  3 exact zero modes (3 translations)

:class:`.genANM` records the ``kappa`` used when the Hessian is built and reports
the corresponding number of trivial modes through :meth:`.genANM.numZeros`, and its
:meth:`.genANM.calcModes` passes that count to the eigensolver so that the correct
number of trivial modes is skipped regardless of ``kappa``.

When ``k_theta == k_phi == kappa == 0`` the model reduces exactly to the classical
ANM built with the same cutoff, spring constant, and sequential-pair convention.
"""

import math

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.utilities import checkCoords, solveEig
from .anm import ANM
from .gnm import checkENMParameters

__author__ = 'Anupam Banerjee'
__credits__ = ['Anupam Banerjee', 'Anthony Bogetti', 'Ivet Bahar']
__email__ = ['anupam.banerjee@stonybrook.edu']

__all__ = ['genANM', 'calcGenANM', 'buildGeneralizedHessian']

#: Model names that select the generalized ANM, shared by :func:`.calcENM`
#: and :func:`.calcNeoAdaptiveANM` so that both always agree.
GENANM_MODEL_NAMES = ('genanm', 'generalized', 'generalized_anm')

#: Keyword arguments understood by :meth:`.genANM.buildHessian` (besides
#: *coords* and *gamma*). :func:`.calcENM` forwards only these, so it accepts
#: and ignores unrelated keywords exactly as it does for the classical ANM.
GENANM_KEYWORDS = ('cutoff', 'k_bond', 'k_theta', 'k_phi', 'kappa', 'include_sequential',
                   'symmetrize', 'torsion_fd_eps', 'enforce_rigid', 'tol')


# ---------------------------------------------------------------------------
# Geometry utilities
# ---------------------------------------------------------------------------

def _wrap_to_pi(angle):
    """Wrap *angle* to ``(-pi, pi]``."""

    return math.atan2(math.sin(angle), math.cos(angle))


def _rigid_body_basis(coords_subset, tol=1e-12):
    """Return an orthonormal basis for the rigid-body motions of a local subset
    of atoms with coordinates *coords_subset* (shape ``(m, 3)``).  The returned
    array has shape ``(3m, r)`` with ``0 <= r <= 6`` depending on the local
    geometric rank (collinear or coincident subsets give ``r < 6``)."""

    coords_subset = np.asarray(coords_subset, dtype=float)
    m = coords_subset.shape[0]
    xc = coords_subset - coords_subset.mean(axis=0, keepdims=True)

    basis = []
    for e in np.eye(3):
        basis.append(np.tile(e, m))
    for e in np.eye(3):
        basis.append(np.cross(np.tile(e, (m, 1)), xc).reshape(-1))

    B = np.column_stack(basis)
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    if s.size == 0 or s[0] <= tol:
        return np.zeros((3 * m, 0), dtype=float)
    keep = s > (tol * s[0])
    return U[:, keep]


def _project_out_rigid_components(grad_flat, coords_subset, tol=1e-12):
    """Project the local gradient *grad_flat* onto the complement of the local
    rigid-body subspace.  Exact internal-coordinate gradients already satisfy
    this; the projection removes finite-difference and roundoff contamination."""

    grad_flat = np.asarray(grad_flat, dtype=float).reshape(-1)
    Q = _rigid_body_basis(coords_subset, tol=tol)
    if Q.shape[1] == 0:
        return grad_flat.copy()
    return grad_flat - Q @ (Q.T @ grad_flat)


def _bonded_mask(atoms, n_atoms, coords=None, break_dist=4.6):
    """Return a boolean array *bonded* of length ``n_atoms - 1`` where
    ``bonded[i]`` is **True** if residues *i* and *i+1* are sequentially bonded,
    i.e. they belong to the same chain and segment and have consecutive residue
    numbers.  Triples and quadruples that span a **False** entry are excluded
    from the angle, torsion, and forced-sequential-pair terms so that no spurious
    bonded stiffness is introduced across chain breaks, missing residues, or
    insertion discontinuities.

    If *atoms* carries no residue/chain information (for example a bare
    :class:`~numpy.ndarray`), a single contiguous chain is assumed and every
    entry is **True**.  In that case the caller is responsible for supplying a
    gap-free selection.

    As a distance backstop, if *coords* is supplied, any sequential pair whose
    Ca-Ca separation exceeds *break_dist* (default 4.6 A) in the reference
    structure is treated as a break regardless of annotation. A covalent virtual
    bond is ~3.8 A; a pair sitting tens of angstroms apart is a chain break even
    when the PDB labels it as one continuous chain with consecutive residue
    numbers (blank chain IDs, renumbered multi-chain complexes). This makes
    genANM robust to un-annotated gaps, matching the ANM's intrinsic tolerance
    of them, so no stiff covalent term is ever placed across a physical gap."""

    if n_atoms < 2:
        return np.zeros(max(n_atoms - 1, 0), dtype=bool)

    bonded = np.ones(n_atoms - 1, dtype=bool)

    # distance backstop, applied regardless of whether annotation is available:
    # a sequential pair physically far apart in the reference is a chain break.
    ref = coords
    if ref is None and not isinstance(atoms, Atomic):
        ref = atoms                      # a bare ndarray was passed as `atoms`
    if ref is not None:
        C = np.asarray(ref, dtype=float)
        if C.ndim == 2 and C.shape[0] == n_atoms:
            dseq = np.linalg.norm(C[1:] - C[:-1], axis=1)
            bonded = bonded & (dseq <= float(break_dist))

    if not isinstance(atoms, Atomic):
        return bonded

    try:
        resnums = atoms.getResnums()
    except Exception:
        resnums = None
    try:
        chids = atoms.getChids()
    except Exception:
        chids = None
    try:
        segs = atoms.getSegnames()
    except Exception:
        segs = None

    for i in range(n_atoms - 1):
        ok = True
        if chids is not None and chids[i] != chids[i + 1]:
            ok = False
        if segs is not None and segs[i] != segs[i + 1]:
            ok = False
        if resnums is not None and int(resnums[i + 1]) - int(resnums[i]) != 1:
            ok = False
        # combine with the distance backstop above; assigning ``ok`` alone would
        # discard it whenever annotation is present
        bonded[i] = bonded[i] and ok

    n_breaks = int(np.count_nonzero(~bonded))
    if n_breaks:
        LOGGER.info('genANM: detected {0} sequential chain break(s) (by '
                    'annotation and/or Ca-Ca distance > {1:.1f} A); angle, '
                    'torsion, and forced-sequential-pair terms will not cross '
                    'them.'.format(n_breaks, break_dist))
    return bonded


# ---------------------------------------------------------------------------
# Hessian contributions
# ---------------------------------------------------------------------------

def _pairwise_blocks(coords, cutoff, k_pair, include_sequential, bonded):
    """Classical ANM pairwise Hessian blocks.  Sequential neighbors are forced
    connected only across genuine bonds when *include_sequential* is **True**."""

    N = coords.shape[0]
    H = np.zeros((3 * N, 3 * N), dtype=float)

    diff = coords[:, None, :] - coords[None, :, :]
    D = np.linalg.norm(diff, axis=-1)

    for i in range(N):
        for j in range(i + 1, N):
            connect = (D[i, j] <= cutoff)
            if include_sequential and (j - i == 1) and bonded[i]:
                connect = True
            if not connect:
                continue

            rij = coords[j] - coords[i]
            r = D[i, j]
            if r < 1e-12:
                continue

            u = rij / r
            block = -k_pair * np.outer(u, u)

            ii = slice(3 * i, 3 * i + 3)
            jj = slice(3 * j, 3 * j + 3)
            H[ii, jj] += block
            H[jj, ii] += block
            H[ii, ii] -= block
            H[jj, jj] -= block

    return H


def _bond_blocks(coords, k_bond, bonded):
    """Virtual Ca-Ca bond-stretching Hessian blocks.

    This is formally the i, i+1 special case of the pairwise term, but carries
    its own force constant so that the covalent bond can be made orders of
    magnitude stiffer than a transient contact. It is the term that prevents the
    backbone from stretching when the structure is DEFORMED along the modes, as
    in adaptive ANM; without it a virtual bond is as soft as any other spring.

    Depends only on the internal coordinate d_i = |r_{i+1} - r_i|, so it
    preserves all six rigid-body zero modes exactly.
    """

    N = coords.shape[0]
    H = np.zeros((3 * N, 3 * N), dtype=float)
    if k_bond == 0.0:
        return H

    for i in range(N - 1):
        if not bonded[i]:
            continue
        rij = coords[i + 1] - coords[i]
        r = np.linalg.norm(rij)
        if r < 1e-12:
            continue

        u = rij / r
        block = -k_bond * np.outer(u, u)

        ii = slice(3 * i, 3 * i + 3)
        jj = slice(3 * (i + 1), 3 * (i + 1) + 3)
        H[ii, jj] += block
        H[jj, ii] += block
        H[ii, ii] -= block
        H[jj, jj] -= block

    return H


def _angle_blocks(coords, k_theta, bonded, enforce_rigid=True, tol=1e-12):
    """Virtual-bond-angle contribution, ``sum_j k_theta * J_j J_j^T``.  A triple
    ``(j-1, j, j+1)`` is included only if both its sequential edges are bonded."""

    N = coords.shape[0]
    H = np.zeros((3 * N, 3 * N), dtype=float)
    if k_theta == 0.0:
        return H

    for j in range(1, N - 1):
        if not (bonded[j - 1] and bonded[j]):
            continue
        i, k = j - 1, j + 1
        b = coords[i] - coords[j]
        c = coords[k] - coords[j]
        nb = np.linalg.norm(b)
        nc = np.linalg.norm(c)
        if nb < tol or nc < tol:
            continue

        bh = b / nb
        ch = c / nc
        ct = float(np.dot(bh, ch))
        st = math.sqrt(max(1.0 - ct * ct, tol))

        dti = (bh * ct - ch) / (nb * st)
        dtk = (ch * ct - bh) / (nc * st)
        dtj = -(dti + dtk)

        g = np.concatenate((dti, dtj, dtk))
        if enforce_rigid:
            g = _project_out_rigid_components(g, coords[[i, j, k]], tol=tol)

        # rank-one update k * g g^T on the 9 x 9 block of the three atoms only (identical values to adding the outer
        # product of the full-length gradient, which is zero elsewhere, at O(1) instead of O(N^2) cost per angle)
        idx = np.r_[3 * i:3 * i + 3, 3 * j:3 * j + 3, 3 * k:3 * k + 3]
        H[np.ix_(idx, idx)] += k_theta * np.outer(g, g)

    return H


def _signed_dihedral(p0, p1, p2, p3, tol=1e-12):
    """Signed dihedral (radians) from a stable atan2 formula, or **None** for a
    degenerate four-atom geometry."""

    b1 = p2 - p1
    n1 = np.linalg.norm(b1)
    if n1 < tol:
        return None
    b1h = b1 / n1

    v = (p0 - p1) - np.dot(p0 - p1, b1h) * b1h
    w = (p3 - p2) - np.dot(p3 - p2, b1h) * b1h
    if np.linalg.norm(v) < tol or np.linalg.norm(w) < tol:
        return None

    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1h, v), w))
    return math.atan2(y, x)


def _torsion_gradient_fd(coords4, fd_eps=1e-6, enforce_rigid=True, tol=1e-12):
    """Central-difference gradient of one signed dihedral with respect to the 12
    Cartesian coordinates of the four atoms, projected onto the complement of the
    local rigid-body subspace.  Returns **None** for degenerate geometries."""

    coords4 = np.asarray(coords4, dtype=float)
    phi0 = _signed_dihedral(coords4[0], coords4[1], coords4[2], coords4[3], tol=tol)
    if phi0 is None:
        return None

    local_scale = max(1.0,
                      float(np.linalg.norm(coords4[1] - coords4[0])),
                      float(np.linalg.norm(coords4[2] - coords4[1])),
                      float(np.linalg.norm(coords4[3] - coords4[2])))
    eps = float(fd_eps) * local_scale

    work = coords4.copy()
    g = np.zeros((4, 3), dtype=float)
    for a in range(4):
        for d in range(3):
            work[a, d] += eps
            phi_p = _signed_dihedral(work[0], work[1], work[2], work[3], tol=tol)
            work[a, d] -= 2.0 * eps
            phi_m = _signed_dihedral(work[0], work[1], work[2], work[3], tol=tol)
            work[a, d] += eps
            if phi_p is None or phi_m is None:
                return None
            g[a, d] = (_wrap_to_pi(phi_p - phi0) - _wrap_to_pi(phi_m - phi0)) / (2.0 * eps)

    g_flat = g.reshape(-1)
    if enforce_rigid:
        g_flat = _project_out_rigid_components(g_flat, coords4, tol=tol)
    return g_flat.reshape(4, 3)


def _torsion_blocks(coords, k_phi, bonded, fd_eps=1e-6, enforce_rigid=True, tol=1e-12):
    """Signed-dihedral contribution, ``sum_q k_phi * J_q J_q^T``.  A quadruple
    ``(q, q+1, q+2, q+3)`` is included only if all three of its sequential edges
    are bonded."""

    N = coords.shape[0]
    H = np.zeros((3 * N, 3 * N), dtype=float)
    if k_phi == 0.0:
        return H

    for i in range(N - 3):
        if not (bonded[i] and bonded[i + 1] and bonded[i + 2]):
            continue
        g4 = _torsion_gradient_fd(coords[i:i + 4], fd_eps=fd_eps,
                                  enforce_rigid=enforce_rigid, tol=tol)
        if g4 is None:
            continue
        # rank-one update on the contiguous 12 x 12 block of the four atoms only (see _angle_blocks)
        g = g4.reshape(-1)
        H[3 * i:3 * (i + 4), 3 * i:3 * (i + 4)] += k_phi * np.outer(g, g)

    return H


def _curvature_blocks(coords, kappa, bonded):
    """Cartesian second-difference penalty with stencil ``(1, -2, 1)``.  A stencil
    centered at *i* is applied only if both its edges are bonded.

    Warning
    -------
    This term preserves translations but not rotations, so a nonzero *kappa*
    reduces the exact zero-mode count from six to three."""

    N = coords.shape[0]
    H = np.zeros((3 * N, 3 * N), dtype=float)
    if kappa == 0.0:
        return H
    I3 = np.eye(3)

    for i in range(1, N - 1):
        if not (bonded[i - 1] and bonded[i]):
            continue
        a, b, c = i - 1, i, i + 1
        aa = slice(3 * a, 3 * a + 3)
        bb = slice(3 * b, 3 * b + 3)
        cc = slice(3 * c, 3 * c + 3)
        H[aa, aa] += kappa * I3
        H[bb, bb] += 4.0 * kappa * I3
        H[cc, cc] += kappa * I3
        H[aa, bb] += -2.0 * kappa * I3
        H[bb, aa] += -2.0 * kappa * I3
        H[bb, cc] += -2.0 * kappa * I3
        H[cc, bb] += -2.0 * kappa * I3
        H[aa, cc] += kappa * I3
        H[cc, aa] += kappa * I3

    return H


def buildGeneralizedHessian(coords, cutoff=15.0, gamma=1.0, k_bond=20.0,
                            k_theta=10.0, k_phi=1.0, kappa=0.0,
                            include_sequential=True,
                            symmetrize=True, bonded=None, torsion_fd_eps=1e-6,
                            enforce_rigid=True, tol=1e-12):
    """Assemble and return the constrained (generalized) ANM Hessian for an
    ``(N, 3)`` coordinate array.  *bonded* is a boolean array of length ``N-1``
    marking sequential bonds; if **None**, a single contiguous chain is assumed.

    With ``kappa == 0`` the six rigid-body zero modes are preserved; with
    ``kappa != 0`` only the three translational zero modes are exact.

    Default constants (gamma = 1, cutoff = 15 A):
    ``k_bond = 20, k_theta = 10, k_phi = 1, kappa = 0``."""

    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError('coords must have shape (N, 3)')
    N = coords.shape[0]

    if bonded is None:
        # no explicit mask: still guard against physical chain breaks by distance,
        # so a bare-array call never bonds across a gap.
        bonded = np.ones(max(N - 1, 0), dtype=bool)
        if N >= 2:
            dseq = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
            bonded = dseq <= 4.6
    else:
        bonded = np.asarray(bonded, dtype=bool)
        if bonded.shape[0] != N - 1:
            raise ValueError('bonded must have length N - 1')

    H = _pairwise_blocks(coords, cutoff, gamma, include_sequential, bonded)
    H += _bond_blocks(coords, k_bond, bonded)
    H += _angle_blocks(coords, k_theta, bonded, enforce_rigid=enforce_rigid, tol=tol)
    H += _torsion_blocks(coords, k_phi, bonded, fd_eps=torsion_fd_eps,
                         enforce_rigid=enforce_rigid, tol=tol)
    if kappa != 0.0:
        H += _curvature_blocks(coords, kappa, bonded)

    if symmetrize:
        H = 0.5 * (H + H.T)
    return H


# ---------------------------------------------------------------------------
# genANM model class
# ---------------------------------------------------------------------------

class genANM(ANM):
    """A constrained (generalized) ANM model.  Subclass of :class:`.ANM`, so it
    is compatible with all downstream ProDy analysis, masking, slicing, and
    trimming.  Build the Hessian with :meth:`.buildHessian` and diagonalize with
    :meth:`.calcModes`."""

    def __init__(self, name='genANM'):
        super(genANM, self).__init__(name)
        self._k_bond = 20.0
        self._k_theta = 10.0
        self._k_phi = 1.0
        self._kappa = 0.0
        self._n_zeros = 6

    def numZeros(self):
        """Return the number of trivial (rigid-body) zero modes expected for the
        Hessian that was built: 6 when ``kappa == 0`` and 3 otherwise."""

        return self._n_zeros

    def getKappa(self):
        """Return the curvature weight ``kappa`` used to build the Hessian."""

        return self._kappa

    def getForceConstants(self):
        """Return ``(k_bond, k_theta, k_phi, kappa)`` used to build the Hessian.

        The default set is ``(20.0, 10.0, 1.0, 0.0)`` at ``gamma = 1`` and
        ``cutoff = 15 A``."""

        return (self._k_bond, self._k_theta, self._k_phi, self._kappa)

    def buildHessian(self, coords, cutoff=15.0, gamma=1.0, **kwargs):
        """Build the constrained (generalized) ANM Hessian.

        :arg coords: a coordinate set or an object with ``getCoords``.  When an
            :class:`.Atomic` instance is passed, its residue numbers, chain ids,
            and segment names are used to detect chain breaks so that angle,
            torsion, and forced-sequential-pair terms never cross them.
        :arg cutoff: pairwise interaction cutoff in Angstroms, default 15.0
        :arg gamma: pairwise spring constant, default 1.0

        Keyword arguments (defaults shown):

        ============= ======= ======================================================
        Keyword       Default Meaning
        ============= ======= ======================================================
        k_bond        20.0    virtual Ca-Ca bond-stretching stiffness
        k_theta       10.0    virtual-bond-angle stiffness
        k_phi         1.0     signed virtual-bond-dihedral stiffness
        kappa         0.0     Cartesian curvature weight; nonzero lifts rotations
        include_seq.  True    force sequential neighbors connected across bonds
        symmetrize    True    symmetrize the assembled Hessian
        ============= ======= ======================================================

        With the default ``kappa=0.0`` the six rigid-body zero modes are
        preserved.  Set ``kappa`` (for example 20.0) to add curvature damping for
        static mobility analyses, in which case only three zero modes remain and
        :meth:`.calcModes` skips three trivial modes accordingly."""

        try:
            coords_arr = (coords._getCoords() if hasattr(coords, '_getCoords')
                          else coords.getCoords())
        except AttributeError:
            coords_arr = coords
            try:
                checkCoords(coords_arr)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object '
                                'with a getCoords method')

        cutoff, g, _ = checkENMParameters(cutoff, gamma)

        k_bond = float(kwargs.pop('k_bond', 20.0))
        k_theta = float(kwargs.pop('k_theta', 10.0))
        k_phi = float(kwargs.pop('k_phi', 1.0))
        kappa = float(kwargs.pop('kappa', 0.0))
        include_sequential = bool(kwargs.pop('include_sequential', True))
        symmetrize = bool(kwargs.pop('symmetrize', True))
        torsion_fd_eps = float(kwargs.pop('torsion_fd_eps', 1e-6))
        enforce_rigid = bool(kwargs.pop('enforce_rigid', True))
        tol = float(kwargs.pop('tol', 1e-12))

        if kwargs.pop('sparse', False):
            raise NotImplementedError('sparse=True is not supported by genANM; '
                                      'the reference implementation is dense.')
        if kwargs:
            raise TypeError('unrecognized keyword arguments: {0}'
                            .format(', '.join(sorted(kwargs))))

        coords_arr = np.asarray(coords_arr, dtype=float)
        if coords_arr.ndim != 2 or coords_arr.shape[1] != 3:
            raise ValueError('coords must have shape (N, 3)')
        n_atoms = coords_arr.shape[0]

        bonded = _bonded_mask(coords, n_atoms, coords=coords_arr)

        LOGGER.timeit('_genanm_build')
        H = buildGeneralizedHessian(coords_arr, cutoff=cutoff, gamma=g,
                                    k_bond=k_bond, k_theta=k_theta,
                                    k_phi=k_phi, kappa=kappa,
                                    include_sequential=include_sequential,
                                    symmetrize=symmetrize, bonded=bonded,
                                    torsion_fd_eps=torsion_fd_eps,
                                    enforce_rigid=enforce_rigid, tol=tol)

        self._reset()
        self._cutoff = cutoff
        self._gamma = g
        self._k_bond = k_bond
        self._k_theta = k_theta
        self._k_phi = k_phi
        self._kappa = kappa
        self._n_zeros = 6 if kappa == 0.0 else 3
        self.setHessian(H)
        self._n_atoms = n_atoms
        self._dof = 3 * n_atoms
        LOGGER.report('genANM Hessian was built in %.2fs.', label='_genanm_build')
        return H

    def calcModes(self, n_modes=20, zeros=False, turbo=True, **kwargs):
        """Calculate normal modes, skipping the correct number of trivial modes
        for the ``kappa`` that was used (six for ``kappa == 0``, three otherwise).
        See :meth:`.ANM.calcModes` for the meaning of the arguments."""

        if self._hessian is None:
            raise ValueError('Hessian matrix is not built or set')
        if str(n_modes).lower() == 'all':
            n_modes = None
        assert n_modes is None or isinstance(n_modes, int) and n_modes > 0, \
            'n_modes must be a positive integer'
        assert isinstance(zeros, bool), 'zeros must be a boolean'
        assert isinstance(turbo, bool), 'turbo must be a boolean'

        self._clear()
        LOGGER.timeit('_genanm_calc_modes')
        values, vectors, vars = solveEig(self._hessian, n_modes=n_modes,
                                         zeros=zeros, turbo=turbo,
                                         expct_n_zeros=self._n_zeros, **kwargs)
        self._eigvals = values
        self._array = vectors
        self._vars = vars
        self._trace = self._vars.sum()
        self._n_modes = len(self._eigvals)
        LOGGER.report('{0} modes were calculated in %.2fs.'
                      .format(self._n_modes), label='_genanm_calc_modes')




# ---------------------------------------------------------------------------
# Convenience driver
# ---------------------------------------------------------------------------

def calcGenANM(pdb, selstr='calpha', cutoff=15.0, gamma=1.0, n_modes=20,
             zeros=False, title=None, **kwargs):
    """Parse *pdb*, build a :class:`.genANM` on the *selstr* selection, calculate
    *n_modes* modes, and return ``(genANM, selection)``.  Extra keyword arguments
    (``k_bond``, ``k_theta``, ``k_phi``, ``kappa``, ``include_sequential``,
    ``symmetrize``) are passed to :meth:`.genANM.buildHessian`."""

    from prody.proteins import parsePDB

    if isinstance(pdb, np.ndarray):
        model = genANM(title or 'Unknown')
        model.setHessian(pdb)
        model.calcModes(n_modes, zeros)
        return model

    if isinstance(pdb, str):
        ag = parsePDB(pdb)
        title = title or ag.getTitle()
    elif isinstance(pdb, Atomic):
        ag = pdb
        if title is None:
            title = (ag.getTitle() if isinstance(pdb, AtomGroup)
                     else ag.getAtomGroup().getTitle())
    else:
        raise TypeError('pdb must be an atomic class, not {0}'.format(type(pdb)))

    model = genANM(title)
    sel = ag.select(selstr)
    model.buildHessian(sel, cutoff, gamma, **kwargs)
    model.calcModes(n_modes, zeros)
    return model, sel


# Earlier name of buildGeneralizedHessian, kept so that existing imports keep working.
# Note that the defaults are now k_bond=20, k_theta=10, k_phi=1 and kappa=0.
build_generalized_hessian = buildGeneralizedHessian
