"""Sortie OBJA : maillage de base, puis une réinsertion de sommet par patch,
niveau par niveau.

La classe `Output` du fichier obja.py du cours est utilisée si elle est présente ;
`ObjaOutput` offre la même interface et la même numérotation (identifiants
quelconques traduits en indices à partir de 1, dans l'ordre de création).
"""
from __future__ import annotations

import math
from collections import namedtuple
from typing import IO, Iterator

import numpy as np

from decoder import Decoder
from levels import ProgressiveMesh

try:
    import obja as _obja
except ImportError:
    _obja = None

Face = namedtuple("Face", "a b c")


class ObjaOutput:
    def __init__(self, stream: IO[str]) -> None:
        self.stream = stream
        self.vertex_mapping: dict[int, int] = {}   # notre numéro -> rang de création
        self.face_mapping: dict[int, int] = {}

    def add_vertex(self, index: int, vertex) -> None:
        self.vertex_mapping[index] = len(self.vertex_mapping)
        print("v {} {} {}".format(*vertex), file=self.stream)

    def add_face(self, index: int, face: Face) -> None:
        self.face_mapping[index] = len(self.face_mapping)
        print("f {} {} {}".format(*self._vertices(face)), file=self.stream)

    def edit_face(self, index: int, face: Face) -> None:
        print("ef {} {} {} {}".format(self.face_mapping[index] + 1, *self._vertices(face)), file=self.stream)

    def delete_face(self, index: int) -> None:
        print("df {}".format(self.face_mapping[index] + 1), file=self.stream)

    def _vertices(self, face: Face) -> tuple[int, int, int]:
        m = self.vertex_mapping
        return m[face.a] + 1, m[face.b] + 1, m[face.c] + 1


def make_output(stream: IO[str]):
    """obja.Output du cours si disponible, sinon notre implémentation équivalente."""
    if _obja is not None and all(hasattr(getattr(_obja, "Output", None), m) for m in ("add_vertex", "add_face", "edit_face")):
        return _obja.Output(stream)
    return ObjaOutput(stream)


LEVEL_COLORS = [(0.85, 0.85, 0.85), (0.95, 0.75, 0.1), (0.9, 0.45, 0.65), (0.15, 0.65, 0.7),
                (0.45, 0.7, 0.25), (0.55, 0.45, 0.85), (0.95, 0.55, 0.3)]


def write_obja(path: str, pm: ProgressiveMesh, declared_sizes: list[int] | None = None,
               color_levels: bool = False) -> None:
    """Pilote un décodeur et traduit chaque événement en ligne OBJA :
    une ligne v, d-2 lignes ef (faces du patch réutilisées) et 2 lignes f par sommet.

    declared_sizes : tailles cumulées (en octets) à déclarer par une ligne « s » après
    M0 puis après chaque niveau (instruction du format OBJA servant à l'évaluation du
    débit) ; None pour ne rien déclarer.
    color_levels : colore les faces créées ou modifiées à chaque niveau (lignes « fc »).
    """
    q = pm.quantizer
    digits = max(0, math.ceil(-math.log10(q.step)) + 1)   # juste assez de décimales pour le pas

    def coords(p: np.ndarray) -> tuple[float, float, float]:
        return tuple(round(float(x), digits) for x in q.dequantize(p))

    decoder = Decoder(pm)
    with open(path, "w") as stream:
        out = make_output(stream)
        mapping = getattr(out, "face_mapping", None)        # numéro OBJA de chaque face

        def color(fid, level: int) -> None:
            if color_levels and mapping is not None:
                r, g, b = LEVEL_COLORS[level % len(LEVEL_COLORS)] if level else LEVEL_COLORS[0]
                print(f"fc {mapping[fid] + 1} {r} {g} {b}", file=stream)

        def declare(level: int) -> None:
            if declared_sizes is not None:
                print(f"s {declared_sizes[level]}", file=stream)

        for v, p in enumerate(pm.base_positions):
            out.add_vertex(v, coords(p))
        for fid, face in decoder.mesh.faces.items():
            out.add_face(fid, Face(*face))
            color(fid, 0)
        for i, face in enumerate(pm.fixed_faces):          # faces mises de côté, jamais modifiées
            out.add_face(("fixe", i), Face(*(int(x) for x in face)))
            color(("fixe", i), 0)
        declare(0)
        while not decoder.done:
            level = decoder.level + 1
            for ref in decoder.step():
                out.add_vertex(ref.vertex, coords(ref.position))
                for fid, face, is_new in ref.faces:
                    (out.add_face if is_new else out.edit_face)(fid, Face(*face))
                    color(fid, level)
            declare(level)


def replay_obja(path: str, snapshots: int = 0, at_declared_sizes: bool = False) -> Iterator[tuple[int, np.ndarray, np.ndarray]]:
    """Exécute un fichier OBJA. Fournit (octets lus, sommets, faces visibles) à
    `snapshots` instants régulièrement espacés, puis une fois à la fin.
    Avec at_declared_sizes, fournit plutôt (taille déclarée, sommets, faces) à chaque
    ligne « s » du fichier."""
    with open(path, "rb") as f:
        lines = f.read().splitlines(keepends=True)
    marks = {round(i * len(lines) / snapshots) for i in range(1, snapshots)} if snapshots > 1 else set()
    vertices: list[list[float]] = []
    faces: list[list[int] | None] = []
    read = 0

    def state() -> tuple[int, np.ndarray, np.ndarray]:
        n = len(vertices)
        # une face qui cite un sommet pas encore créé est ignorée dans l'instantané
        f = np.array([x for x in faces if x is not None and max(x) < n], dtype=np.int64).reshape(-1, 3)
        return read, np.array(vertices, dtype=np.float64).reshape(-1, 3), f

    for i, raw in enumerate(lines, 1):
        read += len(raw)
        parts = raw.split()
        if parts:
            op = parts[0]
            if op == b"s" and at_declared_sizes:
                n = len(vertices)
                f = np.array([x for x in faces if x is not None and max(x) < n], dtype=np.int64).reshape(-1, 3)
                yield int(parts[1]), np.array(vertices, dtype=np.float64).reshape(-1, 3), f
                continue
            if op == b"v":
                vertices.append([float(x) for x in parts[1:4]])
            elif op == b"f":
                faces.append([int(x.split(b"/")[0]) - 1 for x in parts[1:4]])
            elif op == b"ef":
                faces[int(parts[1]) - 1] = [int(x) - 1 for x in parts[2:5]]
            elif op == b"df":
                faces[int(parts[1]) - 1] = None
            elif op == b"ev":
                vertices[int(parts[1]) - 1] = [float(x) for x in parts[2:5]]
            elif op == b"tv":
                idx = int(parts[1]) - 1
                vertices[idx] = [a + float(b) for a, b in zip(vertices[idx], parts[2:5])]
        if i in marks and not at_declared_sizes:
            yield state()
    if not at_declared_sizes:
        yield state()
