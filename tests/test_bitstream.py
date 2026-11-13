"""Flux binaire : Huffman, rangement des étiquettes, aller-retour complet."""
import numpy as np
import pytest

import bitstream
from coloring import pack_labels, unpack_labels
from levels import HierarchyParams, encode


def test_huffman_roundtrip():
    rng = np.random.default_rng(0)
    for values in (rng.integers(-40, 40, 3000), np.zeros(7, dtype=int), np.array([5]), np.array([], dtype=int)):
        assert np.array_equal(bitstream.huffman_decode(bitstream.huffman_encode(values)), values)


@pytest.mark.parametrize("scheme, top", [(4, 4), (2, 2)])
def test_label_packing_roundtrip(scheme, top):
    c = np.random.default_rng(1).integers(0, top, 1001).astype(np.uint8)
    assert np.array_equal(unpack_labels(pack_labels(c, scheme), len(c), scheme), c)


@pytest.mark.parametrize("coloring", ["4", "2", "auto"])
def test_stream_roundtrip(torus, coloring):
    pm = encode(*torus, params=HierarchyParams(coloring=coloring))
    data, stats = bitstream.encode(pm)
    back = bitstream.decode(data)
    assert stats.total_bytes == len(data)
    assert back.quantizer == pm.quantizer and back.num_vertices == pm.num_vertices
    assert np.array_equal(back.base_positions, pm.base_positions)
    assert np.array_equal(back.base_faces, pm.base_faces)
    for a, b in zip(back.levels, pm.levels):
        assert a.scheme == b.scheme
        assert np.array_equal(a.labels, b.labels) and np.array_equal(a.residuals, b.residuals)
