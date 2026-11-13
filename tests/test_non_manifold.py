"""Maillages non-variétés : faces mises de côté, sommets verrouillés, sans perte."""
import numpy as np
import pytest

import bitstream
from decoder import Decoder
from levels import HierarchyParams, encode
from mesh import Mesh, canonical_face


def damaged_sphere(sphere):
    """Sphère + une face « aileron » sur une arête déjà pleine (non-variété)
    + une face retournée qui réutilise une demi-arête existante."""
    v, f = sphere
    a, b, _ = f[0]
    fin = np.array([[a, b, len(v)]])                  # 3e face sur l'arête a-b
    flipped = np.array([f[5]])                        # même face, même orientation : doublon
    extra_vertex = (v[a] + v[b]) / 2 + 0.3
    return np.vstack([v, extra_vertex]), np.vstack([f, fin, flipped])


def test_tolerant_loading_sets_bad_faces_aside(sphere):
    v, f = damaged_sphere(sphere)
    with pytest.raises(ValueError):
        Mesh.from_faces(f)
    mesh, fixed = Mesh.from_faces_tolerant(f)
    assert len(fixed) == 2 and mesh.num_faces == len(f) - 2


@pytest.mark.parametrize("coloring", ["4", "2", "auto"])
def test_lossless_with_fixed_faces(sphere, coloring):
    v, f = damaged_sphere(sphere)
    pm = encode(v, f, params=HierarchyParams(coloring=coloring))
    assert len(pm.fixed_faces) == 2 and pm.levels
    decoded = bitstream.decode(bitstream.encode(pm)[0])
    decoded.vertex_order = pm.vertex_order
    pos, faces = Decoder(decoded).run().current()
    order = pm.vertex_order
    assert np.array_equal(pos, pm.quantizer.quantize(v)[order])
    got = sorted(canonical_face(tuple(int(order[x]) for x in t)) for t in faces)
    assert got == sorted(canonical_face(tuple(int(x) for x in t)) for t in f)
