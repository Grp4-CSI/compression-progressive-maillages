"""Préparation des données du visualiseur .pmc (sans fenêtre OpenGL)."""
import os
import sys

import bitstream
from levels import encode

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))
from view_pmc import build_snapshots, triangle_arrays  # noqa: E402


def test_snapshots_follow_the_levels(sphere):
    pm = encode(*sphere)
    data, stats = bitstream.encode(pm)
    snaps = build_snapshots(data)
    assert len(snaps) == len(pm.levels) + 1
    assert [s.bits for s in snaps] == stats.cumulative_bits()
    assert len(snaps[0].vertices) == len(pm.base_positions) and len(snaps[-1].vertices) == len(sphere[0])
    assert not snaps[0].new.any() and all(s.new.any() for s in snaps[1:])
    pts, normals, colors = triangle_arrays(snaps[-1], highlight=True)
    assert pts.shape == normals.shape == colors.shape == (3 * len(snaps[-1].faces), 3)
