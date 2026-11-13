"""Prédiction entière exacte."""
import numpy as np

from predict import predict, reconstruct, residual


def test_prediction_is_rounded_barycentre():
    pts = np.array([[0, 0, 0], [3, 0, 1], [0, 3, 2]])
    assert predict(pts).tolist() == [1, 1, 1]


def test_reconstruction_is_exact():
    rng = np.random.default_rng(0)
    for d in range(3, 13):
        pts = rng.integers(0, 4096, size=(d, 3))
        v = rng.integers(0, 4096, size=3)
        assert np.array_equal(reconstruct(pts, residual(v, pts)), v)
