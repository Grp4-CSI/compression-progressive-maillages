"""Un niveau de décimation : indépendance, bilan des faces, validité, couleurs."""
import pytest

from decimation import decimate_level
from io_utils import Quantizer
from mesh import Mesh


def decimate(shape, scheme=4):
    v, f = shape
    q = Quantizer.fit(v).quantize(v)
    mesh = Mesh.from_faces(f)
    before = mesh.copy()
    return before, mesh, decimate_level(mesh, q, scheme=scheme)


@pytest.mark.parametrize("scheme", [4, 2])
def test_removed_vertices_form_an_independent_set(sphere, scheme):
    before, _, removals = decimate(sphere, scheme)
    removed = {r.vertex for r in removals}
    assert removals
    for v in removed:
        assert not (before.neighbors(v) & removed)


@pytest.mark.parametrize("scheme", [4, 2])
def test_each_removal_removes_one_vertex_and_two_faces(sphere, torus, scheme):
    for shape in (sphere, torus):
        before, after, removals = decimate(shape, scheme)
        assert after.num_vertices == before.num_vertices - len(removals)
        assert after.num_faces == before.num_faces - 2 * len(removals)
        assert all(len(r.patch) == len(r.ring) - 2 == len(r.labels) for r in removals)


@pytest.mark.parametrize("scheme", [4, 2])
def test_result_is_a_valid_manifold(sphere, grid, scheme):
    for shape in (sphere, grid):
        _, after, _ = decimate(shape, scheme)
        Mesh.from_faces(after.face_array())        # lève une erreur si non-variété
        assert all(after.valence(v) >= 3 for v in after.vertices() if not after.is_boundary_vertex(v))


def test_adjacent_patches_have_different_colors(sphere):
    _, _, removals = decimate(sphere)
    owners = {}
    for i, r in enumerate(removals):
        for f in r.patch:
            for e in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
                owners.setdefault(frozenset(e), set()).add(i)
    shared = [o for o in owners.values() if len(o) == 2]
    assert shared, "le maillage de test doit contenir des patchs adjacents"
    for i, j in shared:
        assert removals[i].labels[0] != removals[j].labels[0]
    assert {r.labels[0] for r in removals} <= {1, 2, 3}


def test_boundary_vertices_are_kept(grid):
    before, _, removals = decimate(grid)
    assert all(not before.is_boundary_vertex(r.vertex) for r in removals)
