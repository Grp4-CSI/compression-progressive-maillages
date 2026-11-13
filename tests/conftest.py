"""Données de test partagées : petites formes générées et dossier des modèles."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, "tools")]

import make_test_meshes as shapes  # noqa: E402


def as_arrays(shape):
    v, f = shape
    return np.asarray(v, dtype=np.float64), np.asarray(f, dtype=np.int64)


@pytest.fixture
def tetra():
    return as_arrays(shapes.tetrahedron())


@pytest.fixture
def icosa():
    return as_arrays(shapes.icosahedron())


@pytest.fixture
def sphere():
    return as_arrays(shapes.icosphere(3))


@pytest.fixture
def torus():
    return as_arrays(shapes.torus())


@pytest.fixture
def grid():
    return as_arrays(shapes.grid())


@pytest.fixture
def models_dir():
    return os.path.join(ROOT, "models")
