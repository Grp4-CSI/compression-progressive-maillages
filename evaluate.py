"""Mesures de distorsion et courbes débit/distorsion."""
from __future__ import annotations

import csv

import numpy as np
from scipy.spatial import cKDTree

from bitstream import StreamStats
from decoder import Decoder
from levels import ProgressiveMesh
from obja_writer import replay_obja


def sample_surface(vertices: np.ndarray, faces: np.ndarray, n: int, seed: int = 0) -> np.ndarray:
    """Échantillonnage uniforme d'un maillage (tirage des triangles proportionnel à leur aire)."""
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    if area.sum() <= 0:
        return vertices.copy()
    rng = np.random.default_rng(seed)
    tri = rng.choice(len(faces), size=n, p=area / area.sum())
    u, v = rng.random(n), rng.random(n)
    flip = u + v > 1                     # repli des points tombés hors du triangle
    u[flip], v[flip] = 1 - u[flip], 1 - v[flip]
    return a[tri] + u[:, None] * (b[tri] - a[tri]) + v[:, None] * (c[tri] - a[tri])


def closest_point_on_triangles(p: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Point du triangle (a, b, c) le plus proche de p, vectorisé (Ericson, RTCD §5.1.5) :
    on teste si la projection tombe près d'un sommet, d'une arête ou à l'intérieur."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = (ab * ap).sum(-1), (ac * ap).sum(-1)
    bp = p - b
    d3, d4 = (ab * bp).sum(-1), (ac * bp).sum(-1)
    cp = p - c
    d5, d6 = (ab * cp).sum(-1), (ac * cp).sum(-1)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2

    denom = va + vb + vc
    denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    v, w = vb / denom, vc / denom
    result = a + ab * v[..., None] + ac * w[..., None]

    def blend(mask, value):             # remplace le résultat là où la condition est vraie
        nonlocal result
        result = np.where(mask[..., None], value, result)

    safe = lambda x: np.where(np.abs(x) < 1e-30, 1e-30, x)
    t_bc = (d4 - d3) / safe((d4 - d3) + (d5 - d6))
    blend((va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0), b + (c - b) * t_bc[..., None])
    t_ac = d2 / safe(d2 - d6)
    blend((vb <= 0) & (d2 >= 0) & (d6 <= 0), a + ac * t_ac[..., None])
    t_ab = d1 / safe(d1 - d3)
    blend((vc <= 0) & (d1 >= 0) & (d3 <= 0), a + ab * t_ab[..., None])
    blend((d6 >= 0) & (d5 <= d6), c)
    blend((d3 >= 0) & (d4 <= d3), b)
    blend((d1 <= 0) & (d2 <= 0), a)
    return result


class SurfaceDistance:
    """Distance point-surface : distance exacte aux k triangles dont les centres sont les plus proches."""

    def __init__(self, vertices: np.ndarray, faces: np.ndarray, k: int = 16) -> None:
        self.tri = vertices[faces]
        self.tree = cKDTree(self.tri.mean(axis=1))
        self.k = min(k, len(faces))

    def __call__(self, points: np.ndarray, chunk: int = 4096) -> np.ndarray:
        out = np.empty(len(points))
        for s in range(0, len(points), chunk):
            p = points[s:s + chunk]
            _, idx = self.tree.query(p, k=self.k)
            idx = idx.reshape(len(p), -1)
            t = self.tri[idx]
            q = closest_point_on_triangles(p[:, None, :], t[..., 0, :], t[..., 1, :], t[..., 2, :])
            out[s:s + chunk] = np.linalg.norm(q - p[:, None, :], axis=-1).min(axis=1)
        return out


class Reference:
    """Surface originale, échantillonnée une fois, à laquelle on compare les approximations."""

    def __init__(self, vertices: np.ndarray, faces: np.ndarray, samples: int = 30000, seed: int = 0) -> None:
        self.points = sample_surface(vertices, faces, samples, seed)
        self.surface = SurfaceDistance(vertices, faces)
        self.diagonal = float(np.linalg.norm(vertices.max(axis=0) - vertices.min(axis=0)))
        self.samples, self.seed = samples, seed

    def distortion(self, vertices: np.ndarray, faces: np.ndarray) -> tuple[float, float]:
        """Distances RMS et de Hausdorff symétriques, divisées par la diagonale de la boîte englobante."""
        if len(faces) == 0:
            return float("nan"), float("nan")
        pts = sample_surface(vertices, faces, self.samples, self.seed + 1)
        d_ab = SurfaceDistance(vertices, faces)(self.points)   # original -> approximation
        d_ba = self.surface(pts)                               # approximation -> original
        rms = np.sqrt((np.mean(d_ab ** 2) + np.mean(d_ba ** 2)) / 2)
        haus = max(d_ab.max(), d_ba.max())
        return rms / self.diagonal, haus / self.diagonal


def stream_curve(pm: ProgressiveMesh, stats: StreamStats, ref: Reference) -> list[dict]:
    """Un point par préfixe du flux compressé (maillage de base, puis chaque niveau)."""
    decoder = Decoder(pm)
    q = pm.quantizer
    rows = []
    for level, bits in enumerate(stats.cumulative_bits()):
        if level > 0:
            decoder.step()
        pos, faces = decoder.current()
        rms, haus = ref.distortion(q.dequantize(pos), faces)
        rows.append({"level": level, "vertices": decoder.count, "bits": bits,
                     "bpv": bits / pm.num_vertices, "rms": rms, "hausdorff": haus})
    return rows


def obja_curve(path: str, ref: Reference, num_vertices: int, points: int = 20) -> list[dict]:
    """Distorsion en fonction des octets lus d'un fichier OBJA, ou des tailles déclarées
    par ses lignes « s » s'il en contient (comme le site d'évaluation)."""
    with open(path, "rb") as f:
        declared = any(line.startswith(b"s ") for line in f)
    rows = []
    for read, vertices, faces in replay_obja(path, snapshots=points, at_declared_sizes=declared):
        rms, haus = ref.distortion(vertices, faces)
        rows.append({"bytes": read, "vertices": len(vertices), "bits": 8 * read,
                     "bpv": 8 * read / num_vertices, "rms": rms, "hausdorff": haus})
    return rows


def write_csv(path: str, rows: list[dict]) -> None:
    """Enregistre les points d'une courbe."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_curves(path: str, curves: dict[str, list[dict]], title: str = "") -> None:
    """Trace RMS et Hausdorff en fonction des bits par sommet (axes logarithmiques)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for name, rows in curves.items():
        x = [r["bpv"] for r in rows]
        axes[0].plot(x, [r["rms"] for r in rows], marker="o", ms=3, label=name)
        axes[1].plot(x, [r["hausdorff"] for r in rows], marker="o", ms=3, label=name)
    for ax, label in zip(axes, ("distance RMS / diagonale", "distance de Hausdorff / diagonale")):
        ax.set_xlabel("bits par sommet (maillage original), échelle log")
        ax.set_xscale("log")
        ax.set_ylabel(label)
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
