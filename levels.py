"""Hiérarchie multirésolution : décimation (du fin au grossier), puis
sérialisation des niveaux de raffinement dans l'ordre de transmission
(du grossier au fin).

Les sommets sont renumérotés dans l'ordre de transmission : d'abord ceux du
maillage de base, puis les sommets réinsérés dans l'ordre où le décodeur les
reconstruit. Codeur et décodeur partagent ainsi la même numérotation, qui
définit l'ordre canonique de lecture des étiquettes.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from coloring import canonical_order, find_level_patches
from decimation import DecimationParams, Removal, decimate_level
from io_utils import Quantizer
from mesh import Mesh, canonical_face
from predict import residual

COLORINGS = ("4", "2", "auto")


@dataclass
class Level:
    scheme: int               # 4 ou 2 couleurs
    labels: np.ndarray        # une étiquette par face du maillage plus grossier, ordre canonique
    residuals: np.ndarray     # (k, 3) entiers, un par patch, ordre de découverte


@dataclass
class ProgressiveMesh:
    quantizer: Quantizer
    base_positions: np.ndarray            # (n0, 3) positions quantifiées de M0
    base_faces: np.ndarray                # (m0, 3) faces de M0, ordre canonique
    levels: list[Level]                   # du plus grossier au plus fin
    num_vertices: int                     # nombre total de sommets à la fin
    fixed_faces: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=np.int64))
    # faces non insérables dans la structure (non-variétés, mal orientées) : transmises
    # telles quelles avec M0, leurs sommets ne sont jamais retirés
    vertex_order: np.ndarray | None = field(default=None, compare=False)  # nouveau -> origine (codeur)


@dataclass(frozen=True)
class HierarchyParams:
    decimation: DecimationParams = DecimationParams()
    coloring: str = "4"                   # "4", "2" ou "auto" (choix par niveau)
    min_vertices: int = 16
    max_levels: int = 64
    min_removed_fraction: float = 0.01


def insert_vertex(mesh: Mesh, patch: list[int], cycle: list[int], v: int) -> list[tuple[int, tuple[int, int, int], bool]]:
    """Remplace un patch par l'éventail de triangles autour du nouveau sommet v.

    Les d-2 numéros de faces du patch sont réutilisés ; deux numéros neufs sont créés.
    Renvoie (numéro de face, face, est_nouvelle) pour chaque face de l'éventail.
    """
    for fid in patch:
        mesh.remove_face(fid)
    d = len(cycle)
    events = []
    for j in range(d):
        face = (v, cycle[j], cycle[(j + 1) % d])
        if j < len(patch):
            events.append((mesh.add_face(*face, fid=patch[j]), face, False))
        else:
            events.append((mesh.add_face(*face), face, True))
    return events


def connectivity_cost(scheme: int, faces_after: int, removed: int) -> float:
    """Bits bruts d'étiquettes par sommet retiré pour un niveau."""
    return (2 if scheme == 4 else 1) * faces_after / removed if removed else float("inf")


def _decimate(mesh: Mesh, pos: np.ndarray, params: HierarchyParams,
              locked: frozenset[int]) -> tuple[Mesh, int, list[Removal]]:
    """Un niveau, avec le schéma imposé ou, en mode auto, le moins coûteux des deux."""
    if params.coloring != "auto":
        scheme = int(params.coloring)
        return mesh, scheme, decimate_level(mesh, pos, params.decimation, scheme, locked=locked)
    best = None
    for scheme in (2, 4):
        trial = mesh.copy()
        removals = decimate_level(trial, pos, params.decimation, scheme, locked=locked)
        cost = connectivity_cost(scheme, trial.num_faces, len(removals))
        if best is None or cost < best[0]:
            best = (cost, trial, scheme, removals)
    return best[1], best[2], best[3]


def build_hierarchy(mesh: Mesh, pos: np.ndarray, params: HierarchyParams = HierarchyParams(),
                    locked: frozenset[int] = frozenset()) -> tuple[Mesh, list[tuple[int, list[Removal]]]]:
    """Décime jusqu'au maillage de base. Renvoie (M0, [(schéma, retraits)] du plus fin au plus grossier)."""
    if params.coloring not in COLORINGS:
        raise ValueError(f"coloriage inconnu : {params.coloring}")
    records: list[tuple[int, list[Removal]]] = []
    while len(records) < params.max_levels and mesh.num_vertices > params.min_vertices:
        n = mesh.num_vertices
        mesh, scheme, removals = _decimate(mesh, pos, params, locked)
        if not removals:
            break
        records.append((scheme, removals))
        if len(removals) < params.min_removed_fraction * n:
            break
    return mesh, records


def serialize(base: Mesh, records: list[tuple[int, list[Removal]]], q: np.ndarray, quantizer: Quantizer,
              fixed: list[tuple[int, int, int]] = ()) -> ProgressiveMesh:
    """Rejoue les niveaux du plus grossier au plus fin, exactement comme le décodeur,
    pour produire les étiquettes et les résidus dans l'ordre où il les lira."""
    # les sommets des faces mises de côté font partie du maillage de base
    base_order = sorted(set(base.vertices()) | {v for f in fixed for v in f})
    order = list(base_order)
    new_id = {v: i for i, v in enumerate(order)}
    base_faces = sorted(canonical_face(tuple(new_id[x] for x in f)) for f in base.faces.values())
    current = Mesh.from_faces(base_faces)
    positions = [q[v] for v in order]

    levels = []
    for scheme, removals in reversed(records):
        owner, label_of = {}, {}
        for idx, r in enumerate(removals):
            for f, lab in zip(r.patch, r.labels):
                key = canonical_face(tuple(new_id[x] for x in f))
                owner[key], label_of[key] = idx, lab

        face_order = canonical_order(current)
        labels = np.array([label_of.get(canonical_face(current.faces[fid]), 0) for fid in face_order], dtype=np.uint8)
        patches = find_level_patches(current, face_order, labels, scheme, positions)
        if len(patches) != len(removals):
            raise RuntimeError("les patchs retrouvés ne correspondent pas aux retraits")

        res = np.zeros((len(patches), 3), dtype=np.int64)
        for i, patch in enumerate(patches):
            keys = {canonical_face(current.faces[fid]) for fid in patch}
            r = removals[owner[next(iter(keys))]]
            if keys != {canonical_face(tuple(new_id[x] for x in f)) for f in r.patch}:
                raise RuntimeError("patch retrouvé incomplet")
            cycle = current.boundary_cycle(patch)
            v = len(order)
            new_id[r.vertex] = v
            order.append(r.vertex)
            res[i] = residual(q[r.vertex], np.array([positions[w] for w in cycle]))
            positions.append(q[r.vertex])
            insert_vertex(current, patch, cycle, v)
        levels.append(Level(scheme, labels, res))

    return ProgressiveMesh(
        quantizer=quantizer,
        base_positions=q[base_order].astype(np.int64).reshape(-1, 3),
        base_faces=np.array(base_faces, dtype=np.int64).reshape(-1, 3),
        levels=levels,
        num_vertices=len(order),
        fixed_faces=np.array([[new_id[v] for v in f] for f in fixed], dtype=np.int64).reshape(-1, 3),
        vertex_order=np.array(order, dtype=np.int64),
    )


def encode(vertices: np.ndarray, faces: np.ndarray, bits: int = 12,
           params: HierarchyParams = HierarchyParams()) -> ProgressiveMesh:
    """Codeur complet : quantification, hiérarchie de décimation, sérialisation."""
    quantizer = Quantizer.fit(vertices, bits)
    q = quantizer.quantize(vertices)
    mesh, fixed = Mesh.from_faces_tolerant(faces)
    locked = frozenset(v for f in fixed for v in f)
    base, records = build_hierarchy(mesh, q, params, locked)
    return serialize(base, records, q, quantizer, fixed)
