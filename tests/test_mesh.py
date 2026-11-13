"""Structure de demi-arêtes : Euler, anneau, bord, garde-fous."""
import pytest

from mesh import Mesh, TopologyError, canonical_face


def euler(mesh):
    edges = {frozenset(e) for e in mesh.halfedges}
    return mesh.num_vertices - len(edges) + mesh.num_faces


def test_closed_meshes_have_expected_euler_characteristic(tetra, icosa, torus):
    assert euler(Mesh.from_faces(tetra[1])) == 2
    assert euler(Mesh.from_faces(icosa[1])) == 2
    assert euler(Mesh.from_faces(torus[1])) == 0


def test_ring_is_ordered_and_complete(icosa):
    mesh = Mesh.from_faces(icosa[1])
    for v in mesh.vertices():
        ring = mesh.ring(v)
        assert len(ring) == 5 == mesh.valence(v)
        for j in range(5):
            # deux voisins consécutifs forment une face (v, w_j, w_{j+1})
            assert mesh.halfedges[(ring[j], ring[(j + 1) % 5])] == mesh.halfedges[(v, ring[j])]


def test_boundary_vertices_have_no_ring(grid):
    mesh = Mesh.from_faces(grid[1])
    assert mesh.ring(0) is None and mesh.is_boundary_vertex(0)
    assert mesh.ring(9) is not None


def test_non_manifold_insertion_is_rejected(tetra):
    mesh = Mesh.from_faces(tetra[1])
    with pytest.raises(TopologyError):
        mesh.add_face(*tetra[1][0])
    with pytest.raises(TopologyError):
        mesh.add_face(0, 0, 1)


def test_remove_then_add_restores_state(icosa):
    mesh = Mesh.from_faces(icosa[1])
    before = (dict(mesh.halfedges), mesh.num_vertices)
    face = mesh.remove_face(3)
    mesh.add_face(*face, fid=3)
    assert (mesh.halfedges, mesh.num_vertices) == before


def test_boundary_cycle_of_vertex_fan(icosa):
    mesh = Mesh.from_faces(icosa[1])
    ring = mesh.ring(0)
    cycle = mesh.boundary_cycle(mesh.faces_around(0))
    assert sorted(cycle) == sorted(ring)
    # le contour suit l'orientation des faces : c'est une rotation de l'anneau
    k = cycle.index(ring[0])
    assert cycle[k:] + cycle[:k] == ring


def test_canonical_face_keeps_orientation():
    assert canonical_face((5, 2, 9)) == (2, 9, 5)
    assert canonical_face((9, 5, 2)) == (2, 9, 5)
