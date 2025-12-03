import numpy as np
import pytest

from prody.dynamics import adaptive2
from prody.ensemble import Ensemble

def test_sequential_edge_lengths():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 1.0, 0.0],
    ], dtype=float)
    elen = adaptive2.sequential_edge_lengths(coords)
    assert elen.shape == (3,)
    assert np.allclose(elen, np.array([1.0, 1.0, 1.0]))

def test_compute_progressive_nmax_no_defvecs():
    # When no defvecs are provided, progressive cap should be at least 1
    out = adaptive2.computeProgressiveNmax([], n_max_modes=10, frac0=0.2, eta=1.5)
    assert isinstance(out, int)
    assert 1 <= out <= 10

def test_calcAdaptiveANM_oneway_zero_steps_smoke():
    # Smoke test: request zero steps -> should return an Ensemble quickly
    a = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 1.5, 0.0],
    ], dtype=float)
    # small perturbation target
    b = a + np.array([0.1, 0.0, 0.0])
    ensemble = adaptive2.calcAdaptiveANM(a, b, n_steps=0, mode=adaptive2.AANM_ONEWAY)
    assert isinstance(ensemble, Ensemble)
    # initial coordset should have been added
    assert ensemble.numCoordsets() >= 1