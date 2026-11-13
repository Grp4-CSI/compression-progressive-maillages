"""Codage de la connectivité par coloriage des faces.

Deux schémas du papier sont disponibles pour chaque niveau :

* 4 couleurs (2 bits par face) : 0 hors patch, 1..3 pour les faces d'un patch,
  deux patchs adjacents ayant toujours des couleurs différentes.
* 2 couleurs (1 bit par face) : chaque patch est triangulé en zigzag ; ses
  triangles du milieu valent 1, ses deux triangles du bout valent 0. Le
  décodeur retrouve les bouts à partir de la forme de la bande de 1.
"""
from __future__ import annotations

from collections import defaultdict, deque
from typing import Iterable, Sequence

import numpy as np

from mesh import Mesh, canonical_face

NUM_COLORS = 3


# ------------------------------------------------------------------ communs
def canonical_order(mesh: Mesh) -> list[int]:
    """Numéros de faces triés par triplet canonique : identique chez le codeur et
    le décodeur, puisqu'ils partagent la même numérotation des sommets."""
    return sorted(mesh.faces, key=lambda fid: canonical_face(mesh.faces[fid]))


def _components(mesh: Mesh, order: list[int], label_of: dict[int, int], keep) -> list[list[int]]:
    """Composantes connexes (par arêtes) des faces de même étiquette retenue,
    dans l'ordre de découverte (parcours en largeur depuis l'ordre canonique)."""
    visited: set[int] = set()
    comps = []
    for seed in order:
        c = label_of[seed]
        if not keep(c) or seed in visited:
            continue
        visited.add(seed)
        comp, queue = [], deque([seed])
        while queue:
            fid = queue.popleft()
            comp.append(fid)
            for u, w in Mesh.face_halfedges(mesh.faces[fid]):
                other = mesh.halfedges.get((w, u))
                if other is not None and other not in visited and label_of[other] == c:
                    visited.add(other)
                    queue.append(other)
        comps.append(comp)
    return comps


# ------------------------------------------------------------- 4 couleurs
def choose_color(neighbor_colors: Iterable[int]) -> int | None:
    """Plus petite couleur de 1..3 non utilisée par un patch voisin, ou None."""
    used = set(neighbor_colors)
    for c in range(1, NUM_COLORS + 1):
        if c not in used:
            return c
    return None


def find_patches(mesh: Mesh, order: list[int], colors: np.ndarray) -> list[list[int]]:
    """4 couleurs : un patch = faces de même couleur non nulle reliées par une arête."""
    if len(order) != len(colors):
        raise ValueError(f"{len(colors)} couleurs pour {len(order)} faces")
    label_of = dict(zip(order, (int(c) for c in colors)))
    return _components(mesh, order, label_of, lambda c: c != 0)


# ------------------------------------------------------------- 2 couleurs
def _sq_length(pos, a: int, b: int) -> int:
    d = np.asarray(pos[a], dtype=np.int64) - np.asarray(pos[b], dtype=np.int64)
    return int(d @ d)


def strip_end_edges(middles: Sequence[tuple[int, int, int]], pos) -> list[frozenset] | None:
    """Arêtes (non orientées) par lesquelles les deux triangles du bout s'attachent
    à la bande des triangles du milieu. None si la forme est ambiguë.

    Règles (identiques au codeur et au décodeur) :
    * 1 triangle du milieu (pentagone) : les bouts touchent ses deux arêtes autres
      que sa plus courte, qui doit être strictement la plus courte ;
    * 2 triangles (hexagone) : le Z n'est jamais en miroir ; le bout est attaché par
      la demi-arête qui entre dans le sommet opposé à l'arête commune ;
    * 3 triangles ou plus : les diagonales d'un zigzag alternent, ce qui fixe le bout.
    """
    m = len(middles)
    if m == 1:
        a, b, c = middles[0]
        edges = [(a, b), (b, c), (c, a)]
        lengths = [_sq_length(pos, u, w) for u, w in edges]
        shortest = sorted(lengths)
        if shortest[0] == shortest[1]:
            return None
        i = lengths.index(shortest[0])
        return [frozenset(e) for j, e in enumerate(edges) if j != i]

    # graphe d'adjacence des triangles du milieu
    owners = defaultdict(list)
    for t, f in enumerate(middles):
        for e in Mesh.face_halfedges(f):
            owners[frozenset(e)].append(t)
    links = defaultdict(list)                       # triangle -> [(voisin, arête commune)]
    for e, ts in owners.items():
        if len(ts) == 2:
            links[ts[0]].append((ts[1], e))
            links[ts[1]].append((ts[0], e))
    ends = [t for t in range(m) if len(links[t]) == 1]
    if len(ends) != 2 or any(len(links[t]) > 2 for t in range(m)):
        return None

    result = []
    for t in ends:
        nxt, shared = links[t][0]
        apex = (set(middles[t]) - shared).pop()    # sommet opposé à l'arête commune
        if m == 2:
            entering = [(u, w) for u, w in Mesh.face_halfedges(middles[t]) if w == apex][0]
            result.append(frozenset(entering))
        else:
            following = [e for k, e in links[nxt] if k != t][0]
            common = shared & following
            if len(common) != 1:
                return None
            other = next(iter(shared - common))
            result.append(frozenset((apex, other)))
    return result


def find_patches_2(mesh: Mesh, order: list[int], bits: np.ndarray, pos) -> list[list[int]]:
    """2 couleurs : chaque bande de faces à 1 plus ses deux faces du bout forme un patch."""
    if len(order) != len(bits):
        raise ValueError(f"{len(bits)} bits pour {len(order)} faces")
    label_of = dict(zip(order, (int(b) for b in bits)))
    patches, used = [], set()
    for comp in _components(mesh, order, label_of, lambda b: b == 1):
        inside = set(comp)
        ends = strip_end_edges([mesh.faces[f] for f in comp], pos)
        if ends is None:
            raise ValueError("bande de triangles ambiguë")
        end_faces = []
        for edge in ends:
            u, w = tuple(edge)
            # la demi-arête appartient à une face de la bande ; le bout est de l'autre côté
            inner = (u, w) if mesh.halfedges.get((u, w)) in inside else (w, u)
            end = mesh.halfedges.get((inner[1], inner[0]))
            if end is None or end in inside or end in used or label_of[end] != 0:
                raise ValueError("triangle du bout introuvable")
            end_faces.append(end)
        used.update(end_faces)
        patches.append([end_faces[0], *comp, end_faces[1]])
    return patches


def find_level_patches(mesh: Mesh, order: list[int], labels: np.ndarray, scheme: int, pos) -> list[list[int]]:
    """Recherche des patchs d'un niveau selon son schéma (4 ou 2 couleurs)."""
    if scheme == 4:
        return find_patches(mesh, order, labels)
    return find_patches_2(mesh, order, labels, pos)


# ------------------------------------------------------- rangement en octets
def pack_labels(labels: np.ndarray, scheme: int) -> bytes:
    """4 couleurs : 2 bits par face (4 par octet) ; 2 couleurs : 1 bit par face."""
    c = np.asarray(labels, dtype=np.uint8)
    if scheme == 2:
        return np.packbits(c).tobytes()
    c = np.concatenate([c, np.zeros((-len(c)) % 4, dtype=np.uint8)]).reshape(-1, 4)
    return (c[:, 0] | (c[:, 1] << 2) | (c[:, 2] << 4) | (c[:, 3] << 6)).astype(np.uint8).tobytes()


def unpack_labels(data: bytes, count: int, scheme: int) -> np.ndarray:
    b = np.frombuffer(data, dtype=np.uint8)
    if scheme == 2:
        return np.unpackbits(b, count=count).astype(np.uint8)
    c = np.stack([b & 3, (b >> 2) & 3, (b >> 4) & 3, (b >> 6) & 3], axis=1).reshape(-1)
    return c[:count].astype(np.uint8)
