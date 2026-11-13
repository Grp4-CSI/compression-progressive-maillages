"""Un niveau de décimation : ensemble indépendant, retrait des sommets,
retriangulation des trous et étiquetage des faces (4 ou 2 couleurs)."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from coloring import choose_color, strip_end_edges
from mesh import Mesh
from predict import residual


@dataclass(frozen=True)
class DecimationParams:
    max_valence: int = 12
    min_cos: float = 0.2      # cosinus minimal entre la normale d'un triangle du patch et celle du sommet


@dataclass
class Removal:
    vertex: int
    ring: list[int]                       # voisins, ordre trigonométrique
    patch: list[tuple[int, int, int]]     # les d-2 faces qui rebouchent le trou
    labels: list[int]                     # étiquette de chaque face du patch (couleur ou bit)


@dataclass(frozen=True)
class Layout:
    """Une triangulation possible du trou, exprimée en positions dans l'anneau."""
    tris: tuple[tuple[int, int, int], ...]
    chords: tuple[tuple[int, int], ...]   # arêtes intérieures créées
    middles: tuple[int, ...]              # 2 couleurs : indices des triangles à 1
    middle_edges: tuple[tuple[int, int], ...]  # arête de contour de chaque triangle à 1


def _edge(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


@lru_cache(maxsize=None)
def _layouts(d: int, scheme: int) -> tuple[Layout, ...]:
    """Les d triangulations candidates d'un trou à d côtés.

    4 couleurs : éventails depuis chaque sommet k.
    2 couleurs : zigzags k, k+1, k-1, k+2, k-2, ... ; les triangles sont orientés
    dans l'ordre cyclique de l'anneau, donc comme le reste du maillage.
    """
    layouts = []
    for k in range(d):
        if scheme == 4:
            tris = tuple((k, (k + j) % d, (k + j + 1) % d) for j in range(1, d - 1))
            chords = tuple((k, (k + j) % d) for j in range(2, d - 1))
            layouts.append(Layout(tris, chords, (), ()))
            continue
        seq = [k] + [(k + (t + 1) // 2) % d if t % 2 else (k - t // 2) % d for t in range(1, d)]
        tris = tuple(tuple(sorted(seq[i:i + 3], key=lambda x: (x - k) % d)) for i in range(d - 2))
        chords = tuple((seq[j], seq[j + 1]) for j in range(1, d - 2))
        middles = tuple(range(1, d - 3))
        middle_edges = tuple((seq[i], seq[i + 2]) for i in middles)
        layouts.append(Layout(tris, chords, middles, middle_edges))
    return tuple(layouts)


def _topology_ok(mesh: Mesh, ring: list[int], valence: list[int], layout: Layout) -> bool:
    """Aucune corde ne doit déjà exister, et chaque sommet du contour garde une valence >= 3."""
    if any(mesh.has_edge(ring[a], ring[b]) for a, b in layout.chords):
        return False
    gain = Counter()
    for a, b in layout.chords:
        gain[a] += 1
        gain[b] += 1
    return all(valence[i] - 1 + gain[i] >= 3 for i in range(len(ring)))


def _decodable(ring: list[int], pos: np.ndarray, layout: Layout) -> bool:
    """2 couleurs : le décodeur doit retrouver exactement les deux triangles du bout."""
    middles = [tuple(ring[i] for i in layout.tris[t]) for t in layout.middles]
    ends = strip_end_edges(middles, pos)
    first, last = layout.chords[0], layout.chords[-1]
    expected = {frozenset((ring[first[0]], ring[first[1]])), frozenset((ring[last[0]], ring[last[1]]))}
    return ends is not None and set(ends) == expected


def _best_layout(mesh: Mesh, pos: np.ndarray, v: int, ring: list[int], scheme: int,
                 params: DecimationParams, ones_edges: set) -> Layout | None:
    """Meilleure triangulation valide du trou laissé par v, ou None.

    Un candidat est classé par le pire cosinus entre la normale de ses triangles et
    la normale du sommet (un triangle retourné ou plat a un cosinus faible).
    """
    d = len(ring)
    p = pos[ring].astype(np.float64)
    rel = p - pos[v]
    normal = np.cross(rel, np.roll(rel, -1, axis=0)).sum(axis=0)
    norm = float(np.linalg.norm(normal))
    if norm < 1e-12:
        return None
    normal /= norm

    layouts = _layouts(d, scheme)
    idx = np.array([lay.tris for lay in layouts])            # (d, d-2, 3)
    a, b, c = p[idx[..., 0]], p[idx[..., 1]], p[idx[..., 2]]
    n = np.cross(b - a, c - a)
    length = np.linalg.norm(n, axis=-1)
    cos = np.where(length > 1e-12, (n @ normal) / np.maximum(length, 1e-12), -1.0)
    score = cos.min(axis=1)

    valence = [mesh.valence(w) for w in ring]
    for k in sorted(range(d), key=lambda k: (-score[k], k)):
        if score[k] <= params.min_cos:
            break
        lay = layouts[k]
        if not _topology_ok(mesh, ring, valence, lay):
            continue
        if scheme == 2:
            # deux triangles à 1 de patchs différents ne doivent jamais se toucher
            if any(_edge(ring[i], ring[j]) in ones_edges for i, j in lay.middle_edges):
                continue
            if not _decodable(ring, pos, lay):
                continue
        return lay
    return None


def decimate_level(mesh: Mesh, pos: np.ndarray, params: DecimationParams = DecimationParams(),
                   scheme: int = 4, stats: Counter | None = None,
                   locked: frozenset[int] = frozenset()) -> list[Removal]:
    """Retire un ensemble indépendant maximal de sommets de `mesh` (sur place).

    Les candidats sont traités par erreur de prédiction croissante : les sommets
    bien prédits (zones plates) partent d'abord, ceux qui portent la forme restent.
    En 2 couleurs, seules les valences >= 5 sont codables. Les sommets de `locked`
    (ceux des faces mises de côté au chargement) ne sont jamais retirés.
    """
    stats = Counter() if stats is None else stats
    min_valence = 5 if scheme == 2 else 3
    candidates = []
    for v in mesh.vertices():
        if v in locked:
            stats["sommet verrouillé"] += 1
            continue
        ring = mesh.ring(v)
        if ring is None:
            stats["bord ou non-variété"] += 1
            continue
        if not min_valence <= len(ring) <= params.max_valence:
            stats["valence hors limites"] += 1
            continue
        err = residual(pos[v], pos[ring])
        candidates.append((int(err @ err), len(ring), v, ring))
    candidates.sort()

    blocked: set[int] = set()
    edge_color: dict[tuple[int, int], int] = {}    # 4 couleurs : couleur du patch par arête de contour
    ones_edges: set[tuple[int, int]] = set()       # 2 couleurs : arêtes de contour bordées par un triangle à 1
    removals: list[Removal] = []
    for _, d, v, ring in candidates:
        if v in blocked:
            stats["voisin déjà retiré"] += 1
            continue
        ring_edges = [_edge(ring[j], ring[(j + 1) % d]) for j in range(d)]
        color = 1
        if scheme == 4:
            color = choose_color(edge_color[e] for e in ring_edges if e in edge_color)
            if color is None:
                stats["pas de couleur libre"] += 1
                continue
        lay = _best_layout(mesh, pos, v, ring, scheme, params, ones_edges)
        if lay is None:
            stats["pas de retriangulation valide"] += 1
            continue

        for fid in mesh.faces_around(v):
            mesh.remove_face(fid)
        patch = [tuple(ring[i] for i in t) for t in lay.tris]
        for face in patch:
            mesh.add_face(*face)
        if scheme == 4:
            labels = [color] * len(patch)
            for e in ring_edges:
                edge_color[e] = color
        else:
            labels = [1 if t in lay.middles else 0 for t in range(len(patch))]
            ones_edges.update(_edge(ring[i], ring[j]) for i, j in lay.middle_edges)
        blocked.add(v)
        blocked.update(ring)
        removals.append(Removal(v, ring, patch, labels))
    return removals
