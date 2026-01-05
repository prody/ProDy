import numpy as np
import pytest

from prody.dynamics.constrained_anm import cANM, calcCANM

def _sample_coords():
    return np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        [0.0, 0.0, 1.5],
    ], dtype=float)

def test_build_hessian_symmetry_and_shape():
    coords = _sample_coords()
    model = cANM('test')
    H = model.buildHessian(coords, cutoff=5.0, gamma=1.0)
    assert H.shape == (12, 12)
    assert np.allclose(H, H.T, atol=1e-8)

def test_calc_modes_and_calcCANM_with_array():
    coords = _sample_coords()
    model = cANM('test2')
    H = model.buildHessian(coords, cutoff=5.0)
    model.calcModes(n_modes=6, zeros=False)
    eigvals = model.getEigvals()
    assert eigvals is not None and eigvals.size > 0

    # calcCANM accepting Hessian array should return an ANM-like model
    m2 = calcCANM(H, n_modes=6, zeros=False)
    assert hasattr(m2, 'getHessian')
    # calcCANM sets a symmetric copy via setHessian in our implementation
    assert np.allclose(m2.getHessian(), 0.5*(H + H.T), atol=1e-8)
