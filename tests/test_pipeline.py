"""Bout en bout : sans perte pour chaque schéma, progressivité, rejeu de l'OBJA."""
import os

import numpy as np
import pytest

import bitstream
from decoder import Decoder
from io_utils import read_obj
from levels import HierarchyParams, encode
from mesh import canonical_face
from obja_writer import replay_obja, write_obja


def assert_lossless(pm, v, f):
    pos, faces = Decoder(pm).run().current()
    order = pm.vertex_order
    assert np.array_equal(pos, pm.quantizer.quantize(v)[order])
    got = {canonical_face(tuple(int(order[x]) for x in t)) for t in faces}
    assert got == {canonical_face(tuple(int(x) for x in t)) for t in f}


@pytest.mark.parametrize("coloring", ["4", "2", "auto"])
@pytest.mark.parametrize("name", ["sphere", "torus", "grid", "icosa", "tetra"])
def test_lossless_on_test_shapes(name, coloring, request):
    v, f = request.getfixturevalue(name)
    pm = encode(v, f, params=HierarchyParams(coloring=coloring))
    decoded = bitstream.decode(bitstream.encode(pm)[0])
    decoded.vertex_order = pm.vertex_order
    assert_lossless(decoded, v, f)


@pytest.mark.parametrize("coloring", ["4", "auto"])
@pytest.mark.parametrize("path", ["fandisk.obj", "edge_cases/cow.obj"])
def test_lossless_on_real_models(path, coloring, models_dir):
    full = os.path.join(models_dir, path)
    if not os.path.exists(full):
        pytest.skip("modèle absent")
    v, f = read_obj(full)
    assert_lossless(encode(v, f, params=HierarchyParams(coloring=coloring)), v, f)


def test_auto_mixes_schemes(sphere):
    pm = encode(*sphere, params=HierarchyParams(coloring="auto"))
    assert {level.scheme for level in pm.levels} == {2, 4}


def test_levels_are_progressive(sphere):
    pm = encode(*sphere)
    decoder = Decoder(pm)
    counts = [decoder.count]
    while not decoder.done:
        decoder.step()
        counts.append(decoder.count)
        assert len(decoder.mesh.faces) == 2 * decoder.count - 4      # sphère fermée à chaque niveau
    assert counts == sorted(counts) and counts[-1] == len(sphere[0])


def test_obja_replay_matches_decoder(torus, tmp_path):
    pm = encode(*torus)
    path = str(tmp_path / "torus.obja")
    write_obja(path, pm)
    *_, (size, vertices, faces) = replay_obja(path)
    pos, ref_faces = Decoder(pm).run().current()
    assert size == os.path.getsize(path)
    assert np.allclose(vertices, pm.quantizer.dequantize(pos), atol=pm.quantizer.step)
    assert {canonical_face(tuple(t)) for t in faces.tolist()} == {canonical_face(tuple(t)) for t in ref_faces.tolist()}


def test_obja_declares_sizes_and_colors(torus, tmp_path):
    pm = encode(*torus)
    _, stats = bitstream.encode(pm)
    sizes = [b // 8 for b in stats.cumulative_bits()]
    path = str(tmp_path / "torus.obja")
    write_obja(path, pm, sizes, color_levels=True)
    snaps = list(replay_obja(path, at_declared_sizes=True))
    assert [s for s, _, _ in snaps] == sizes                     # une ligne s après M0 et après chaque niveau
    assert len(snaps[-1][1]) == len(torus[0])                    # la dernière correspond au maillage complet
    assert sum(1 for line in open(path) if line.startswith("fc ")) > 0
