"""Évalue le programme sur un ensemble de modèles : tailles, temps, et distorsion à
plusieurs débits, pour les schémas 4 couleurs et auto. Produit un tableau Markdown
et une figure par modèle.

    python tools/benchmark.py modeles/*.obj -o out/benchmark
"""
from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import bitstream  # noqa: E402
from evaluate import Reference, plot_curves, stream_curve  # noqa: E402
from io_utils import read_obj  # noqa: E402
from levels import HierarchyParams, encode  # noqa: E402


def at_fraction(rows: list[dict], fraction: float, key: str) -> float:
    """Distorsion au premier point dont le débit atteint `fraction` du débit total."""
    total = rows[-1]["bits"]
    return next(r[key] for r in rows if r["bits"] >= fraction * total)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("models", nargs="+")
    parser.add_argument("-o", "--outdir", default="out/benchmark")
    parser.add_argument("--samples", type=int, default=15000)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lines = ["| modèle | sommets | schéma | niveaux | M0 | connect. | géom. | total (b/s) | Hausdorff à 10 % | à 25 % | à 50 % | codage |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for path in args.models:
        name = os.path.splitext(os.path.basename(path))[0]
        vertices, faces = read_obj(path)
        n = len(vertices)
        ref = Reference(vertices, faces, samples=args.samples)
        curves = {}
        for coloring in ("4", "auto"):
            t0 = time.perf_counter()
            pm = encode(vertices, faces, params=HierarchyParams(coloring=coloring))
            dt = time.perf_counter() - t0
            _, s = bitstream.encode(pm)
            rows = stream_curve(pm, s, ref)
            curves[f"coloriage {coloring}"] = rows
            h = [at_fraction(rows, f, "hausdorff") for f in (0.10, 0.25, 0.50)]
            lines.append(f"| {name} | {n} | {coloring} | {len(pm.levels)} | {len(pm.base_positions)} | "
                         f"{8 * sum(s.connectivity_bytes) / n:.2f} | {8 * sum(s.geometry_bytes) / n:.2f} | "
                         f"{8 * s.total_bytes / n:.2f} | {h[0]:.4f} | {h[1]:.4f} | {h[2]:.4f} | {dt:.1f} s |")
            print(lines[-1], flush=True)
        plot_curves(os.path.join(args.outdir, f"{name}.png"), curves, title=name)
    table = "\n".join(lines)
    with open(os.path.join(args.outdir, "resultats.md"), "w") as f:
        f.write("Distorsions relatives à la diagonale de la boîte englobante ; "
                "débits en bits par sommet du modèle original (flux .pmc).\n\n" + table + "\n")
    print(f"\ntableau : {os.path.join(args.outdir, 'resultats.md')}")


if __name__ == "__main__":
    main()
