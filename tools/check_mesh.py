"""Diagnostic d'un maillage OBJ : est-il utilisable pour notre méthode ?

Usage : python check_mesh.py models/cow.obj [autres.obj ...]
"""
import sys
from collections import Counter, defaultdict


def read_obj(path):
    """Lit un OBJ : renvoie (sommets, faces). Ignore normales/textures, garde les polygones tels quels."""
    verts, faces = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "v":
                verts.append(tuple(float(x) for x in parts[1:4]))
            elif parts[0] == "f":
                # "12/5/3" -> 12 ; indices OBJ à partir de 1 -> on passe à 0
                faces.append(tuple(int(p.split("/")[0]) - 1 for p in parts[1:]))
    return verts, faces


def diagnose(path):
    verts, faces = read_obj(path)
    sizes = Counter(len(f) for f in faces)
    tris = [f for f in faces if len(f) == 3]

    # arêtes orientées et non orientées
    edge_faces = defaultdict(int)    # {a,b} -> nb de faces
    directed = Counter()             # (a,b) -> nb d'occurrences
    for a, b, c in tris:
        for u, v in ((a, b), (b, c), (c, a)):
            edge_faces[frozenset((u, v))] += 1
            directed[(u, v)] += 1

    used = {v for f in tris for v in f}
    boundary = sum(1 for n in edge_faces.values() if n == 1)
    nonmanifold = sum(1 for n in edge_faces.values() if n > 2)
    bad_orient = sum(1 for n in directed.values() if n > 1)
    degenerate = sum(1 for f in tris if len(set(f)) < 3)

    # sommets non-variétés : l'anneau de triangles autour du sommet doit être d'un seul tenant
    around = defaultdict(list)
    for t in tris:
        for v in t:
            around[v].append(t)
    pinched = 0
    for v, ts in around.items():
        # composantes connexes des triangles autour de v (liés s'ils partagent une arête contenant v)
        parent = list(range(len(ts)))
        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        by_edge = defaultdict(list)
        for i, t in enumerate(ts):
            for w in t:
                if w != v:
                    by_edge[w].append(i)
        for idx in by_edge.values():
            for j in idx[1:]:
                parent[find(j)] = find(idx[0])
        if len({find(i) for i in range(len(ts))}) > 1:
            pinched += 1

    valence = Counter()
    for e in edge_faces:
        for v in e:
            valence[v] += 1
    V, E, F = len(used), len(edge_faces), len(tris)

    print(f"== {path}")
    print(f"  sommets : {len(verts)} (utilisés : {V})   faces : {len(faces)}  tailles : {dict(sizes)}")
    print(f"  arêtes : {E}   Euler V-E+F = {V - E + F}")
    if valence:
        vals = list(valence.values())
        print(f"  valence : moyenne {sum(vals)/len(vals):.2f}, min {min(vals)}, max {max(vals)}")
    print(f"  arêtes de bord : {boundary}   non-variétés (>2 faces) : {nonmanifold}")
    print(f"  orientation incohérente : {bad_orient}   triangles dégénérés : {degenerate}"
          f"   sommets non-variétés : {pinched}")
    clean = sizes.keys() <= {3} and nonmanifold == bad_orient == degenerate == pinched == 0
    if clean and boundary == 0:
        print("  => fermé et propre : idéal pour commencer")
    elif clean:
        print("  => propre mais avec des bords : bon test des cas limites")
    else:
        print("  => à nettoyer avant usage")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        diagnose(p)
