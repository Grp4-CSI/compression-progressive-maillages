"""Décodeur progressif : reconstruit le maillage niveau par niveau."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from coloring import canonical_order, find_level_patches
from levels import ProgressiveMesh, insert_vertex
from mesh import Mesh
from predict import reconstruct


@dataclass
class Refinement:
    vertex: int
    position: np.ndarray                                   # quantifiée
    faces: list[tuple[int, tuple[int, int, int], bool]]     # (numéro de face, face, est_nouvelle)


class Decoder:
    def __init__(self, pm: ProgressiveMesh) -> None:
        self.pm = pm
        self.mesh = Mesh.from_faces(pm.base_faces)
        self.positions = np.zeros((pm.num_vertices, 3), dtype=np.int64)
        self.count = len(pm.base_positions)
        self.positions[: self.count] = pm.base_positions
        self.level = 0

    @property
    def done(self) -> bool:
        return self.level >= len(self.pm.levels)

    def step(self) -> list[Refinement]:
        """Applique le niveau suivant : retrouve les patchs, reconstruit et réinsère les sommets."""
        level = self.pm.levels[self.level]
        order = canonical_order(self.mesh)
        patches = find_level_patches(self.mesh, order, level.labels, level.scheme, self.positions)
        if len(patches) != len(level.residuals):
            raise ValueError(f"niveau {self.level} : {len(patches)} patchs pour {len(level.residuals)} résidus")
        refinements = []
        for patch, res in zip(patches, level.residuals):
            cycle = self.mesh.boundary_cycle(patch)
            v = self.count
            self.count += 1
            self.positions[v] = reconstruct(self.positions[cycle], res)
            events = insert_vertex(self.mesh, patch, cycle, v)
            refinements.append(Refinement(v, self.positions[v].copy(), events))
        self.level += 1
        return refinements

    def run(self) -> "Decoder":
        while not self.done:
            self.step()
        return self

    def current(self) -> tuple[np.ndarray, np.ndarray]:
        """Maillage courant : (positions quantifiées, faces), faces mises de côté comprises."""
        faces = np.array(self.mesh.face_array(), dtype=np.int64).reshape(-1, 3)
        return self.positions[: self.count], np.vstack([faces, self.pm.fixed_faces.reshape(-1, 3)])
