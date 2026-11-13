"""Lecture des fichiers OBJ et quantification des coordonnées."""
from __future__ import annotations

import struct
from dataclasses import dataclass

import numpy as np


def read_obj(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Lit un fichier OBJ et renvoie (sommets, faces).

    Les polygones sont découpés en éventail, les indices de texture et de
    normale sont ignorés, et les sommets non utilisés sont supprimés (les faces
    sont renumérotées en conséquence).
    """
    vertices, faces = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "v":
                vertices.append([float(x) for x in parts[1:4]])
            elif parts[0] == "f":
                idx = [int(p.split("/")[0]) for p in parts[1:]]
                # indices OBJ à partir de 1 ; un indice négatif compte depuis la fin
                idx = [i - 1 if i > 0 else len(vertices) + i for i in idx]
                for k in range(1, len(idx) - 1):
                    faces.append((idx[0], idx[k], idx[k + 1]))

    vertices = np.asarray(vertices, dtype=np.float64).reshape(-1, 3)
    faces = np.asarray(faces, dtype=np.int64).reshape(-1, 3)
    # suppression des faces dégénérées (deux fois le même sommet)
    faces = faces[(faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 2] != faces[:, 0])]

    # suppression des sommets inutilisés et renumérotation
    used = np.unique(faces)
    remap = np.full(len(vertices), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    return vertices[used], remap[faces]


def write_obj(path: str, vertices: np.ndarray, faces: np.ndarray) -> None:
    with open(path, "w") as f:
        for x, y, z in vertices:
            f.write(f"v {x:.7g} {y:.7g} {z:.7g}\n")
        for a, b, c in faces:
            f.write(f"f {a + 1} {b + 1} {c + 1}\n")


@dataclass(frozen=True)
class Quantizer:
    """Quantification uniforme des coordonnées sur `bits` bits, sur la boîte englobante.

    Le même pas est utilisé pour les trois axes (calculé sur la plus grande
    dimension) afin de ne pas déformer l'objet.
    """

    origin: tuple[float, float, float]
    step: float
    bits: int

    # en-tête binaire : 1 octet (bits) + 4 flottants 64 bits (origine, pas) = 33 octets
    HEADER = struct.Struct("<B4d")

    @classmethod
    def fit(cls, vertices: np.ndarray, bits: int = 12) -> "Quantizer":
        if not 1 <= bits <= 16:
            raise ValueError("bits doit être compris entre 1 et 16")
        lo = vertices.min(axis=0)
        extent = float((vertices.max(axis=0) - lo).max())
        step = extent / ((1 << bits) - 1) if extent > 0 else 1.0
        return cls(tuple(float(x) for x in lo), step, bits)

    def quantize(self, vertices: np.ndarray) -> np.ndarray:
        q = np.rint((vertices - np.asarray(self.origin)) / self.step).astype(np.int64)
        return np.clip(q, 0, (1 << self.bits) - 1)

    def dequantize(self, q: np.ndarray) -> np.ndarray:
        return np.asarray(self.origin) + np.asarray(q, dtype=np.float64) * self.step

    def to_bytes(self) -> bytes:
        return self.HEADER.pack(self.bits, *self.origin, self.step)

    @classmethod
    def from_bytes(cls, data: bytes) -> "Quantizer":
        bits, ox, oy, oz, step = cls.HEADER.unpack(data)
        return cls((ox, oy, oz), step, bits)
