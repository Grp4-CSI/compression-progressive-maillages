"""Prédiction géométrique d'un sommet retiré à partir du contour de son patch.

Tous les calculs sont faits sur des entiers : codeur et décodeur obtiennent
exactement la même prédiction, sans dérive d'un niveau à l'autre.
"""
from __future__ import annotations

import numpy as np


def predict(boundary_positions: np.ndarray) -> np.ndarray:
    """Barycentre du contour, arrondi à l'entier le plus proche : floor((2s + d) / 2d)."""
    pts = np.asarray(boundary_positions, dtype=np.int64)
    d = len(pts)
    return np.floor_divide(2 * pts.sum(axis=0) + d, 2 * d)


def residual(position: np.ndarray, boundary_positions: np.ndarray) -> np.ndarray:
    """Côté codeur : position réelle moins prédiction."""
    return np.asarray(position, dtype=np.int64) - predict(boundary_positions)


def reconstruct(boundary_positions: np.ndarray, res: np.ndarray) -> np.ndarray:
    """Côté décodeur : prédiction plus résidu."""
    return predict(boundary_positions) + np.asarray(res, dtype=np.int64)
