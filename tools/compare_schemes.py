"""Compare les schémas de coloriage (4, 2, auto) sur plusieurs modèles.

    python tools/compare_schemes.py models/fandisk.obj models/edge_cases/cow.obj
    python tools/compare_schemes.py models/fandisk.obj --plot out/schemes.png
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
from levels import COLORINGS, HierarchyParams, encode  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("models", nargs="+")
    parser.add_argument("--plot", help="trace les courbes du flux pour chaque schéma (premier modèle)")
    args = parser.parse_args()

    print(f"{'modèle':16s} {'schéma':6s} {'niveaux':>7s} {'M0':>5s} {'connect.':>9s} {'géom.':>7s} {'total':>7s} {'codage':>7s}")
    curves = {}
    for i, path in enumerate(args.models):
        vertices, faces = read_obj(path)
        n = len(vertices)
        ref = Reference(vertices, faces) if args.plot and i == 0 else None
        for coloring in COLORINGS:
            t0 = time.perf_counter()
            pm = encode(vertices, faces, params=HierarchyParams(coloring=coloring))
            dt = time.perf_counter() - t0
            _, s = bitstream.encode(pm)
            conn, geom = 8 * sum(s.connectivity_bytes) / n, 8 * sum(s.geometry_bytes) / n
            name = os.path.splitext(os.path.basename(path))[0]
            print(f"{name:16s} {coloring:6s} {len(pm.levels):7d} {len(pm.base_positions):5d} "
                  f"{conn:9.2f} {geom:7.2f} {8 * s.total_bytes / n:7.2f} {dt:6.1f}s")
            if ref is not None:
                curves[f"coloriage {coloring}"] = stream_curve(pm, s, ref)
    if curves:
        plot_curves(args.plot, curves, title="flux compressé selon le coloriage")
        print(f"courbes : {args.plot}")


if __name__ == "__main__":
    main()
