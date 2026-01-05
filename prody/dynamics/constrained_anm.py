# -*- coding: utf-8 -*-
"""Constrained ANM (cANM)
A ProDy-compatible ANM subclass that builds a constrained Hessian with:
 - pairwise distance springs (ANM-style)
 - 3-body bond-angle terms
 - 4-body torsional/dihedral terms
 - backbone curvature (discrete second-difference) term

This module provides:
 - class cANM(ANM)
 - calcCANM(...) convenience function, analogous to calcANM in anm.py

Notes
-----
This implementation uses dense NumPy arrays (no sparse support). For very large
systems, consider implementing a sparse version.
"""

from __future__ import annotations

__author__ = 'Anupam Banerjee'
__credits__ = ['Anthony Bogetti']
__email__ = ['anupam.banerjee@stonybrook.edu', 'anthony.bogetti@stonybrook.edu']

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.utilities import checkCoords
from prody.dynamics.anm import ANM
from prody.dynamics.gnm import checkENMParameters

__all__ = ['cANM', 'calcCANM']

# ---------------------------------------------------------------------------
# Internal helpers to assemble constrained Hessian
# ---------------------------------------------------------------------------

def _pairwise_blocks(coords: np.ndarray, cutoff: float, k_pair: float, include_sequential: bool) -> np.ndarray:
    N = coords.shape[0]
    H = np.zeros((3*N, 3*N), dtype=float)
    diff = coords[:, None, :] - coords[None, :, :]
    D = np.linalg.norm(diff, axis=-1)
    for i in range(N):
        for j in range(i+1, N):
            connect = (D[i, j] <= cutoff)
            if include_sequential and abs(i - j) == 1:
                connect = True
            if not connect:
                continue
            rij = coords[j] - coords[i]
            r = np.linalg.norm(rij)
            if r < 1e-12:
                continue
            u = rij / r
            block = -k_pair * np.outer(u, u)
            ii = slice(3*i, 3*i+3)
            jj = slice(3*j, 3*j+3)
            H[ii, jj] += block
            H[jj, ii] += block
            H[ii, ii] -= block
            H[jj, jj] -= block
    return H


def _angle_blocks(coords: np.ndarray, k_theta: float) -> np.ndarray:
    N = coords.shape[0]
    H = np.zeros((3*N, 3*N), dtype=float)
    for j in range(1, N-1):
        i, k = j - 1, j + 1
        b = coords[i] - coords[j]
        c = coords[k] - coords[j]
        nb = np.linalg.norm(b)
        nc = np.linalg.norm(c)
        if nb < 1e-12 or nc < 1e-12:
            continue
        bh = b / nb
        ch = c / nc
        ct = float(np.dot(bh, ch))
        st = np.sqrt(max(1.0 - ct*ct, 1e-12))
        dti = (bh * ct - ch) / (nb * st)
        dtk = (ch * ct - bh) / (nc * st)
        dtj = -(dti + dtk)
        J = np.zeros((1, 3*N), dtype=float)
        J[0, 3*i:3*i+3] = dti
        J[0, 3*j:3*j+3] = dtj
        J[0, 3*k:3*k+3] = dtk
        H += k_theta * (J.T @ J)
    return H


def _torsion_blocks(coords: np.ndarray, k_phi: float) -> np.ndarray:
    N = coords.shape[0]
    H = np.zeros((3*N, 3*N), dtype=float)
    for i in range(N-3):
        j, k, l = i+1, i+2, i+3
        b1 = coords[j] - coords[i]
        b2 = coords[k] - coords[j]
        b3 = coords[l] - coords[k]
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        A1 = float(np.dot(n1, n1))
        A2 = float(np.dot(n2, n2))
        B = float(np.linalg.norm(b2))
        if A1 < 1e-10 or A2 < 1e-10 or B < 1e-10:
            continue
        dpi = -(B / A1) * n1
        dpl =  (B / A2) * n2
        dpj = (np.dot(b1, b2) / (B*B)) * dpi - (np.dot(b3, b2) / (B*B)) * dpl
        dpk = -(dpi + dpj + dpl)
        J = np.zeros((1, 3*N), dtype=float)
        J[0, 3*i:3*i+3] = dpi
        J[0, 3*j:3*j+3] = dpj
        J[0, 3*k:3*k+3] = dpk
        J[0, 3*l:3*l+3] = dpl
        H += k_phi * (J.T @ J)
    return H


def _backbone_curvature_blocks(coords: np.ndarray, kappa: float) -> np.ndarray:
    N = coords.shape[0]
    H = np.zeros((3*N, 3*N), dtype=float)
    I3 = np.eye(3)
    for i in range(1, N-1):
        a, b, c = i - 1, i, i + 1
        aa = slice(3*a, 3*a+3)
        bb = slice(3*b, 3*b+3)
        cc = slice(3*c, 3*c+3)
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


def build_constrained_hessian(
    coords: np.ndarray,
    cutoff: float = 10.0,
    gamma: float = 1.0,
    k_theta: float = 0.5,
    k_phi: float = 0.2,
    kappa: float = 0.8,
    include_sequential: bool = True,
    symmetrize: bool = True,
) -> np.ndarray:
    """Assemble the constrained Hessian from an (N,3) coordinate array."""
    H = (
        _pairwise_blocks(coords, cutoff=cutoff, k_pair=gamma, include_sequential=include_sequential)
        + _angle_blocks(coords, k_theta=k_theta)
        + _torsion_blocks(coords, k_phi=k_phi)
        + _backbone_curvature_blocks(coords, kappa=kappa)
    )
    if symmetrize:
        H = 0.5 * (H + H.T)
    return H


# ---------------------------------------------------------------------------
# ProDy-compatible subclass
# ---------------------------------------------------------------------------

class cANM(ANM):
    """Constrained ANM that behaves like a ProDy ANM/NMA object."""

    def __init__(self, name: str = 'cANM') -> None:
        super().__init__(name)

    def buildHessian(self, coords, cutoff=15., gamma=1., **kwargs):
        """Build constrained Hessian from coordinates or an object exposing getCoords().

        Additional keyword arguments (defaults shown):
          k_theta=0.5
          k_phi=0.2
          kappa=0.8
          include_sequential=True
          symmetrize=True

        Note: this implementation uses dense arrays (no 'sparse' support).
        """
        try:
            coords = (coords._getCoords() if hasattr(coords, '_getCoords') else coords.getCoords())
        except AttributeError:
            try:
                checkCoords(coords)
            except TypeError:
                raise TypeError('coords must be a Numpy array or an object with `getCoords` method')

        #cutoff, g, gamma_checked = checkENMParameters(cutoff, gamma)
        cutoff, g, gamma_func = checkENMParameters(cutoff, gamma)

        # constrained-specific parameters
        k_theta = float(kwargs.pop('k_theta', 0.5))
        k_phi = float(kwargs.pop('k_phi', 0.2))
        kappa = float(kwargs.pop('kappa', 0.8))
        include_sequential = bool(kwargs.pop('include_sequential', True))
        symmetrize = bool(kwargs.pop('symmetrize', True))

        # Do not support sparse construction in this reference implementation
        if kwargs.pop('sparse', False):
            raise NotImplementedError('sparse=True is not supported by constrained ANM (dense NumPy used).')

        coords = np.asarray(coords, dtype=float)
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError('coords must have shape (N, 3)')

        self._reset()
        self._cutoff = cutoff
        self._gamma = g
        n_atoms = coords.shape[0]
        #H = build_constrained_hessian(
        #    coords,
        #    cutoff=cutoff,
        #    gamma=gamma_checked,
        #    k_theta=k_theta,
        #    k_phi=k_phi,
        #    kappa=kappa,
        #    include_sequential=include_sequential,
        #    symmetrize=symmetrize,
        #)

        H = build_constrained_hessian(
            coords,
            cutoff=cutoff,
            gamma=g,  # <--- Use 'g' (the value) instead of the function
            k_theta=k_theta,
            k_phi=k_phi,
            kappa=kappa,
            include_sequential=include_sequential,
            symmetrize=symmetrize,
        )

        # setHessian performs validation and sets internal sizes
        self.setHessian(H)
        self._n_atoms = n_atoms
        self._dof = n_atoms * 3
        return H


# ---------------------------------------------------------------------------
# Convenience functions
# ---------------------------------------------------------------------------

def calcCANM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20,
             zeros=False, title=None, **kwargs):
    """Perform cANM calculation and return (cANM, selection) similar to calcANM.

    Accepts the same kinds of inputs as calcANM:
      - numpy.ndarray (Hessian or coords) => if array is Hessian, accepted directly
      - string PDB code => parsed with parsePDB
      - Atomic instance => selection on selstr is taken
    """
    import numpy as _np
    from prody.proteins import parsePDB

    if isinstance(pdb, _np.ndarray):
        H = pdb
        if title is None:
            title = 'Unknown'
        model = cANM(title)
        model.setHessian(H)
        model.calcModes(n_modes, zeros)
        return model

    else:
        if isinstance(pdb, str):
            ag = parsePDB(pdb)
            if title is None:
                title = ag.getTitle()
        elif isinstance(pdb, Atomic):
            ag = pdb
            if title is None:
                if isinstance(pdb, AtomGroup):
                    title = ag.getTitle()
                else:
                    title = ag.getAtomGroup().getTitle()
        else:
            raise TypeError('pdb must be an atomic class, not {0}'.format(type(pdb)))

        model = cANM(title)
        sel = ag.select(selstr)
        model.buildHessian(sel, cutoff, gamma, **kwargs)
        model.calcModes(n_modes, zeros)
        return model, sel
