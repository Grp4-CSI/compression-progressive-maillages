"""Génère de petits maillages de test (OBJ triangulaires, orientés de façon cohérente).

Usage : python make_test_meshes.py   -> écrit dans ../models/test/
"""
import math
import os

OUT = os.path.join(os.path.dirname(__file__), "..", "models", "test")


def write_obj(name, verts, faces):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w") as f:
        f.write(f"# {name} : {len(verts)} sommets, {len(faces)} faces\n")
        for x, y, z in verts:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        for a, b, c in faces:
            f.write(f"f {a + 1} {b + 1} {c + 1}\n")   # OBJ : indices à partir de 1


def tetrahedron():
    v = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    f = [(0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2)]
    return v, f


def octahedron():
    v = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    f = [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
         (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)]
    return v, f


def icosahedron():
    t = (1 + 5 ** 0.5) / 2
    v = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
         (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
         (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    f = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
         (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
         (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
         (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    return v, f


def subdivide_sphere(verts, faces):
    """Coupe chaque triangle en 4 et projette les nouveaux sommets sur la sphère."""
    verts = list(verts)
    cache = {}

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in cache:
            p = [(verts[a][i] + verts[b][i]) / 2 for i in range(3)]
            n = math.sqrt(sum(c * c for c in p))
            verts.append(tuple(c / n for c in p))
            cache[key] = len(verts) - 1
        return cache[key]

    new = []
    for a, b, c in faces:
        ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
        new += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
    return verts, new


def icosphere(level):
    v, f = icosahedron()
    n = [math.sqrt(sum(c * c for c in p)) for p in v]
    v = [tuple(c / n[i] for c in p) for i, p in enumerate(v)]
    for _ in range(level):
        v, f = subdivide_sphere(v, f)
    return v, f


def grid(n=8, m=8):
    """Grille plate n x m avec bord (légère bosse en z pour que la prédiction ait du sens)."""
    v = [(i, j, 0.3 * math.sin(i / 2) * math.cos(j / 2)) for j in range(m) for i in range(n)]
    f = []
    for j in range(m - 1):
        for i in range(n - 1):
            a, b, c, d = j * n + i, j * n + i + 1, (j + 1) * n + i, (j + 1) * n + i + 1
            f += [(a, b, d), (a, d, c)]
    return v, f


def torus(n=24, m=12, R=2.0, r=0.7):
    v = []
    for i in range(n):
        for j in range(m):
            u, w = 2 * math.pi * i / n, 2 * math.pi * j / m
            v.append(((R + r * math.cos(w)) * math.cos(u),
                      (R + r * math.cos(w)) * math.sin(u),
                      r * math.sin(w)))
    f = []
    for i in range(n):
        for j in range(m):
            a = i * m + j
            b = ((i + 1) % n) * m + j
            c = ((i + 1) % n) * m + (j + 1) % m
            d = i * m + (j + 1) % m
            f += [(a, b, c), (a, c, d)]
    return v, f


if __name__ == "__main__":
    write_obj("tetrahedron.obj", *tetrahedron())
    write_obj("octahedron.obj", *octahedron())
    write_obj("icosahedron.obj", *icosahedron())
    for lvl in (1, 2, 3):
        write_obj(f"icosphere{lvl}.obj", *icosphere(lvl))
    write_obj("grid8.obj", *grid())
    write_obj("torus.obj", *torus())
    print("Maillages écrits dans", os.path.abspath(OUT))
