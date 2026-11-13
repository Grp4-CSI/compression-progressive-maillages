# Modèles de test

## models/test/ — petits maillages générés (tools/make_test_meshes.py)
| Fichier | Sommets | Usage |
|---|---|---|
| tetrahedron.obj | 4 | tests unitaires minimaux (valence 3) |
| octahedron.obj | 6 | valence 4 partout (cas où le 2-coloriage échoue) |
| icosahedron.obj | 12 | valence 5 partout |
| icosphere1/2/3.obj | 42 / 162 / 642 | sphères régulières : premiers essais multi-niveaux |
| torus.obj | 288 | fermé, genre 1, valence 6 partout |
| grid8.obj | 64 | surface avec bord (test des sommets de bord) |

## models/ — modèles réels fermés et propres
| Fichier | Sommets | Remarque |
|---|---|---|
| fandisk.obj | 6 475 | pièce mécanique, arêtes vives |
| armadillo.obj | 49 990 | gros modèle pour les courbes finales |
| stanford-bunny.obj | 35 947 | 223 arêtes de bord (trous), 1 113 sommets inutilisés |

## models/edge_cases/ — à ne PAS utiliser au début
| Fichier | Problème |
|---|---|
| cow.obj | 1 sommet non-variété (deux éventails qui se touchent en un point) |
| teapot.obj | 1 036 arêtes de bord, 38 sommets non-variétés |

Source des modèles réels : github.com/alecjacobson/common-3d-test-models.
Diagnostic : `python tools/check_mesh.py models/xxx.obj`
