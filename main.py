"""Compression progressive d'un maillage triangulaire (Cohen-Or, Levin & Remez, 1999).

    python main.py models/fandisk.obj                      # 4 couleurs (par défaut)
    python main.py models/fandisk.obj --coloring auto      # 2 ou 4 couleurs selon le niveau
    python main.py models/fandisk.obj --compare naive=out/fandisk_naive.obja
    python main.py models/fandisk.obj --color-levels --export-levels out/niveaux

L'OBJA contient une ligne « s » (taille cumulée en octets) après M0 puis après chaque
niveau : par défaut, la taille du flux compressé .pmc, c'est-à-dire le vrai débit.
"""
from __future__ import annotations

import argparse
import os
import time
from collections import Counter

import numpy as np

import bitstream
from decimation import DecimationParams
from decoder import Decoder
from evaluate import Reference, obja_curve, plot_curves, stream_curve, write_csv
from io_utils import read_obj
from levels import COLORINGS, HierarchyParams, encode
from mesh import canonical_face
from io_utils import write_obj
from obja_writer import replay_obja, write_obja


def check_lossless(pm, vertices: np.ndarray, faces: np.ndarray) -> None:
    """Décode entièrement et compare à l'entrée quantifiée (positions et connectivité)."""
    pos, out_faces = Decoder(pm).run().current()
    order = pm.vertex_order
    if not np.array_equal(pos, pm.quantizer.quantize(vertices)[order]):
        raise AssertionError("positions décodées différentes de l'entrée quantifiée")
    got = {canonical_face(tuple(int(order[x]) for x in f)) for f in out_faces}
    expected = {canonical_face(tuple(int(x) for x in f)) for f in faces}
    if got != expected:
        raise AssertionError("connectivité décodée différente de l'entrée")


def check_obja(path: str, pm) -> None:
    """Rejoue le fichier OBJA écrit et vérifie qu'il redonne le maillage décodé."""
    *_, (_, vertices, faces) = replay_obja(path)
    pos, ref_faces = Decoder(pm).run().current()
    ok = (len(vertices) == len(pos) and np.allclose(vertices, pm.quantizer.dequantize(pos), atol=pm.quantizer.step)
          and {canonical_face(tuple(t)) for t in faces.tolist()} == {canonical_face(tuple(t)) for t in ref_faces.tolist()})
    if not ok:
        raise AssertionError("le rejeu du fichier OBJA ne redonne pas le maillage décodé")


def export_levels(directory: str, name: str, pm) -> None:
    """Écrit un OBJ par niveau (M0, M1, ...) pour les visualiser avec n'importe quel logiciel."""
    os.makedirs(directory, exist_ok=True)
    decoder = Decoder(pm)
    while True:
        pos, faces = decoder.current()
        write_obj(os.path.join(directory, f"{name}_niveau{decoder.level:02d}.obj"), pm.quantizer.dequantize(pos), faces)
        if decoder.done:
            break
        decoder.step()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="fichier .obj d'entrée")
    parser.add_argument("-o", "--outdir", default="out", help="dossier de sortie")
    parser.add_argument("--coloring", choices=COLORINGS, default="4",
                        help="schéma de coloriage : 4 couleurs, 2 couleurs, ou choix par niveau (auto)")
    parser.add_argument("--bits", type=int, default=12, help="bits de quantification par coordonnée (<= 16)")
    parser.add_argument("--max-valence", type=int, default=12)
    parser.add_argument("--min-cos", type=float, default=0.2)
    parser.add_argument("--min-vertices", type=int, default=16, help="taille cible du maillage de base")
    parser.add_argument("--declared-size", choices=("pmc", "none"), default="pmc",
                        help="taille déclarée par les lignes « s » de l'OBJA : celle du flux compressé (.pmc) ou aucune")
    parser.add_argument("--color-levels", action="store_true", help="colore les faces de chaque niveau dans l'OBJA (lignes fc)")
    parser.add_argument("--export-levels", metavar="DOSSIER", help="écrit un OBJ par niveau dans ce dossier")
    parser.add_argument("--no-eval", action="store_true", help="sans courbes débit/distorsion")
    parser.add_argument("--compare", action="append", default=[], metavar="NOM=FICHIER.obja",
                        help="ajoute un autre fichier OBJA aux courbes")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    name = os.path.splitext(os.path.basename(args.input))[0]
    out = lambda ext: os.path.join(args.outdir, name + ext)

    vertices, faces = read_obj(args.input)
    print(f"{name} : {len(vertices)} sommets, {len(faces)} faces, coloriage {args.coloring}")

    t0 = time.perf_counter()
    params = HierarchyParams(DecimationParams(args.max_valence, args.min_cos), args.coloring,
                             min_vertices=args.min_vertices)
    pm = encode(vertices, faces, args.bits, params)
    t_enc = time.perf_counter() - t0

    data, stats = bitstream.encode(pm)
    with open(out(".pmc"), "wb") as f:
        f.write(data)

    # preuve du sans-perte : on relit le flux écrit, on décode et on compare
    t0 = time.perf_counter()
    decoded = bitstream.decode(data)
    Decoder(decoded).run()
    t_dec = time.perf_counter() - t0
    decoded.vertex_order = pm.vertex_order
    check_lossless(decoded, vertices, faces)
    sizes = [bits // 8 for bits in stats.cumulative_bits()] if args.declared_size == "pmc" else None
    write_obja(out(".obja"), pm, sizes, args.color_levels)
    check_obja(out(".obja"), pm)
    if args.export_levels:
        export_levels(args.export_levels, name, pm)

    n = len(vertices)
    conn, geom = sum(stats.connectivity_bytes), sum(stats.geometry_bytes)
    schemes = Counter(level.scheme for level in pm.levels)
    print(f"niveaux : {len(pm.levels)} ({schemes[4]} en 4 couleurs, {schemes[2]} en 2 couleurs)   "
          f"maillage de base : {len(pm.base_positions)} sommets, {len(pm.base_faces)} faces")
    print(f"codage {t_enc:.2f} s   décodage {t_dec:.2f} s   sans perte : OK   rejeu OBJA : OK")
    print(f"flux {stats.total_bytes} o = {8 * stats.total_bytes / n:.2f} bits/sommet "
          f"(connectivité {8 * conn / n:.2f}, géométrie {8 * geom / n:.2f}, "
          f"en-tête + base {8 * (stats.header_bytes + stats.base_bytes) / n:.2f})")
    print(f"obja   {os.path.getsize(out('.obja'))} o   obj {os.path.getsize(args.input)} o")

    if args.no_eval:
        return
    ref = Reference(vertices, faces)
    curves = {"flux compressé (.pmc)": stream_curve(pm, stats, ref)}
    if args.declared_size == "none":
        curves["fichier .obja"] = obja_curve(out(".obja"), ref, n)
    for spec in args.compare:
        label, path = spec.split("=", 1)
        curves[label] = obja_curve(path, ref, n)
    for i, rows in enumerate(curves.values()):
        write_csv(out(f"_rd_{i}.csv"), rows)
    plot_curves(out("_rd.png"), curves, title=f"{name} ({args.coloring} couleurs)")
    print(f"débit/distorsion : {out('_rd.png')}")


if __name__ == "__main__":
    main()
