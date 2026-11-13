"""2-coloriage : zigzags, règles de lecture de la bande, décodage sans perte."""
import numpy as np
import pytest

from coloring import strip_end_edges
from decimation import _layouts, decimate_level
from io_utils import Quantizer
from mesh import Mesh


def regular_polygon(d, radius=1000):
    angles = 2 * np.pi * np.arange(d) / d
    return np.stack([2048 + radius * np.cos(angles), 2048 + radius * np.sin(angles), np.zeros(d)], axis=1).astype(np.int64)


@pytest.mark.parametrize("d", range(5, 11))
def test_strip_rules_recover_the_end_triangles(d):
    pos = regular_polygon(d)
    for lay in _layouts(d, 2):
        middles = [lay.tris[t] for t in lay.middles]
        expected = {frozenset(lay.chords[0]), frozenset(lay.chords[-1])}
        assert set(strip_end_edges(middles, pos)) == expected


@pytest.mark.parametrize("d", range(5, 11))
def test_zigzag_counts_and_orientation(d):
    for lay in _layouts(d, 2):
        assert len(lay.tris) == d - 2 and len(lay.middles) == d - 4 and len(lay.chords) == d - 3
        # les triangles d'un trou triangulé forment un disque orienté valide
        Mesh.from_faces(lay.tris).boundary_cycle(range(d - 2))


def test_ones_of_different_patches_never_touch(sphere):
    v, f = sphere
    q = Quantizer.fit(v).quantize(v)
    removals = decimate_level(Mesh.from_faces(f), q, scheme=2)
    owner = {}
    for i, r in enumerate(removals):
        for face, lab in zip(r.patch, r.labels):
            if lab == 1:
                for e in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
                    owner.setdefault(frozenset(e), set()).add(i)
    assert all(len(o) == 1 for o in owner.values())
    assert all(sum(r.labels) == len(r.ring) - 4 for r in removals)
    assert min(len(r.ring) for r in removals) >= 5
