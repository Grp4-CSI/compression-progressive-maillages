"""Maillage triangulaire en demi-arêtes.

Les demi-arêtes sont stockées implicitement comme arêtes orientées dans un
dictionnaire : (a, b) -> numéro de la face qui contient la demi-arête a->b.
La jumelle de (a, b) est (b, a) ; la suivante se lit dans l'ordre cyclique de
la face. Une modification topologique revient donc à insérer ou effacer trois
clés par face, sans pointeurs à recâbler.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Iterator


class TopologyError(ValueError):
    """Levée quand une opération casserait la structure de variété orientée."""


class Mesh:
    def __init__(self) -> None:
        self.faces: dict[int, tuple[int, int, int]] = {}          # numéro -> (a, b, c)
        self.halfedges: dict[tuple[int, int], int] = {}           # (a, b) -> numéro de face
        self._out: defaultdict[int, set[int]] = defaultdict(set)  # v -> {w : v->w existe}
        self._in: defaultdict[int, set[int]] = defaultdict(set)   # v -> {u : u->v existe}
        self._next_fid = 0                                         # prochain numéro de face libre

    @classmethod
    def from_faces(cls, faces: Iterable[Iterable[int]]) -> "Mesh":
        mesh = cls()
        for a, b, c in faces:
            mesh.add_face(int(a), int(b), int(c))
        return mesh

    @classmethod
    def from_faces_tolerant(cls, faces: Iterable[Iterable[int]]) -> tuple["Mesh", list[tuple[int, int, int]]]:
        """Comme from_faces, mais les faces qui casseraient la variété (arête partagée
        par plus de deux faces, face mal orientée) sont mises de côté au lieu de lever
        une erreur. Renvoie (maillage, faces mises de côté)."""
        mesh, rejected = cls(), []
        for a, b, c in faces:
            try:
                mesh.add_face(int(a), int(b), int(c))
            except TopologyError:
                rejected.append((int(a), int(b), int(c)))
        return mesh, rejected

    # ------------------------------------------------------------------ faces
    def add_face(self, a: int, b: int, c: int, fid: int | None = None) -> int:
        """Ajoute la face orientée (a, b, c), éventuellement sous un numéro imposé."""
        if a == b or b == c or c == a:
            raise TopologyError(f"face dégénérée {(a, b, c)}")
        hes = ((a, b), (b, c), (c, a))
        for he in hes:
            # une demi-arête déjà utilisée signale une arête non-variété ou une face retournée
            if he in self.halfedges:
                raise TopologyError(f"demi-arête {he} déjà utilisée")
        if fid is None:
            fid = self._next_fid
        elif fid in self.faces:
            raise TopologyError(f"numéro de face {fid} déjà utilisé")
        self._next_fid = max(self._next_fid, fid + 1)
        self.faces[fid] = (a, b, c)
        for u, v in hes:
            self.halfedges[(u, v)] = fid
            self._out[u].add(v)
            self._in[v].add(u)
        return fid

    def remove_face(self, fid: int) -> tuple[int, int, int]:
        face = self.faces.pop(fid)
        a, b, c = face
        for u, v in ((a, b), (b, c), (c, a)):
            del self.halfedges[(u, v)]
            self._out[u].discard(v)
            self._in[v].discard(u)
        for x in face:
            if not self._out[x] and not self._in[x]:
                del self._out[x], self._in[x]
        return face

    @staticmethod
    def face_halfedges(face: tuple[int, int, int]) -> tuple[tuple[int, int], ...]:
        a, b, c = face
        return (a, b), (b, c), (c, a)

    # -------------------------------------------------------------- requêtes
    def vertices(self) -> Iterator[int]:
        return iter(self._out.keys() | self._in.keys())

    @property
    def num_vertices(self) -> int:
        return len(self._out.keys() | self._in.keys())

    @property
    def num_faces(self) -> int:
        return len(self.faces)

    def neighbors(self, v: int) -> set[int]:
        return self._out.get(v, set()) | self._in.get(v, set())

    def valence(self, v: int) -> int:
        return len(self.neighbors(v))

    def has_edge(self, a: int, b: int) -> bool:
        return (a, b) in self.halfedges or (b, a) in self.halfedges

    def faces_around(self, v: int) -> list[int]:
        """Éventail de v : une face par demi-arête sortante."""
        return [self.halfedges[(v, w)] for w in self._out.get(v, ())]

    def is_boundary_vertex(self, v: int) -> bool:
        # sommet intérieur : chaque voisin est à la fois sortant et entrant
        return self._out.get(v, set()) != self._in.get(v, set())

    def ring(self, v: int) -> list[int] | None:
        """Voisins de v dans l'ordre trigonométrique, ou None si v est sur le bord
        ou non-variété (plusieurs éventails qui se touchent en v)."""
        out = self._out.get(v)
        if not out or self.is_boundary_vertex(v):
            return None
        start = min(out)
        ring = [start]
        w = start
        for _ in range(len(out)):
            a, b, c = self.faces[self.halfedges[(v, w)]]
            w = a + b + c - v - w          # troisième sommet de la face (v, w, x)
            if w == start:
                break
            ring.append(w)
        else:
            return None
        return ring if len(ring) == len(out) else None

    def boundary_cycle(self, fids: Iterable[int]) -> list[int]:
        """Contour orienté d'un ensemble de faces formant un disque."""
        fids = list(fids)
        inside = set(fids)
        succ: dict[int, int] = {}
        for fid in fids:
            for u, w in self.face_halfedges(self.faces[fid]):
                # demi-arête de contour : sa jumelle est hors du patch (ou absente)
                if self.halfedges.get((w, u)) not in inside:
                    if u in succ:
                        raise TopologyError("le contour du patch n'est pas une boucle simple")
                    succ[u] = w
        start = self.faces[fids[0]][0]
        if start not in succ:
            start = next(iter(succ))
        cycle = [start]
        while (w := succ[cycle[-1]]) != start:
            cycle.append(w)
            if len(cycle) > len(succ):
                raise TopologyError("le contour du patch n'est pas une boucle simple")
        if len(cycle) != len(succ):
            raise TopologyError("le contour du patch comporte plusieurs boucles")
        return cycle

    def face_array(self) -> list[tuple[int, int, int]]:
        return list(self.faces.values())

    def copy(self) -> "Mesh":
        other = Mesh()
        other.faces = dict(self.faces)
        other.halfedges = dict(self.halfedges)
        other._out = defaultdict(set, {k: set(s) for k, s in self._out.items()})
        other._in = defaultdict(set, {k: set(s) for k, s in self._in.items()})
        other._next_fid = self._next_fid
        return other


def canonical_face(face: tuple[int, int, int]) -> tuple[int, int, int]:
    """Rotation de la face qui place le plus petit sommet en tête (orientation conservée)."""
    a, b, c = face
    if a < b and a < c:
        return a, b, c
    if b < c:
        return b, c, a
    return c, a, b
