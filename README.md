# Compression progressive de maillages triangulaires

Implémentation de Cohen-Or, Levin & Remez, *Progressive Compression of Arbitrary
Triangular Meshes* (IEEE Visualization 1999), dans le cadre du projet OBJA
(fork de [tforgione/obja](https://git.tforgione.fr/tforgione/obja)).

Le programme transforme un `.obj` en `.obja` progressif et produit le flux binaire
compressé du papier (`.pmc`). Le décodage est vérifié sans perte à chaque exécution,
et le fichier OBJA écrit est rejoué pour vérifier qu'il redonne le même maillage.

## Installation
```
git clone <url de ce dépôt>
cd <dossier>
pip install -r requirements.txt
```

## Utilisation
```
python main.py example/suzanne.obj                    # -> out/suzanne.obja, out/suzanne.pmc, courbes
python main.py example/suzanne.obj --no-eval          # sans courbes (plus rapide)
python main.py example/suzanne.obj --coloring auto    # 2 ou 4 couleurs selon le niveau
python main.py example/suzanne.obj --color-levels     # faces colorées par niveau (démo)
python main.py example/suzanne.obj --export-levels out/niveaux   # un OBJ par niveau
```
Visualisation du streaming (serveur fourni par le projet OBJA) :
```
./server.py            # puis ouvrir http://localhost:8000/?out/suzanne.obja
```

| Option | Défaut | Rôle |
|---|---|---|
| `--coloring` | `4` | schéma de coloriage : `4`, `2` ou `auto` |
| `--declared-size` | `pmc` | taille déclarée par les lignes `s` de l'OBJA : celle du flux compressé, ou `none` |
| `--color-levels` | — | lignes `fc` : une couleur par niveau |
| `--export-levels` | — | écrit M0, M1, … en OBJ dans un dossier |
| `--bits` | `12` | bits de quantification par coordonnée |
| `--min-vertices` | `16` | taille cible du maillage de base |
| `--no-eval` | — | sans courbes débit/distorsion |
| `--compare NOM=F.obja` | — | ajoute un autre OBJA aux courbes |

Le fichier OBJA contient une ligne `s` après le maillage de base puis après chaque
niveau : elle déclare la taille cumulée du flux compressé, c'est-à-dire le vrai débit.

## Outils
```
python -m pytest tests                              # 66 tests
python tools/check_mesh.py example/suzanne.obj        # diagnostic d'un maillage
python tools/compare_schemes.py models/fandisk.obj  # 4 / 2 / auto
python tools/benchmark.py example/*.obj             # tableau et courbes par modèle
```
Résultats sur les modèles du site : `docs/benchmark/resultats.md`.

## Organisation
| Fichier | Rôle |
|---|---|
| `obja.py`, `server.py`, `decimate.py`, `index.html`, `js/`, `src/` | fournis par le projet OBJA |
| `io_utils.py` | lecture OBJ, quantification |
| `mesh.py` | maillage en demi-arêtes |
| `predict.py` | prédiction barycentrique entière |
| `coloring.py` | 4 et 2 couleurs : étiquettes, ordre canonique, recherche des patchs |
| `decimation.py` | un niveau : ensemble indépendant, retrait, retriangulation |
| `levels.py` | hiérarchie, choix du schéma, sérialisation, `encode()` |
| `decoder.py` | reconstruction niveau par niveau |
| `bitstream.py` | flux binaire : zlib (étiquettes) + Huffman (résidus) |
| `obja_writer.py` | écriture (lignes `v`, `f`, `ef`, `s`, `fc`) et rejeu d'OBJA |
| `evaluate.py` | distances point-surface, courbes débit/distorsion |
| `main.py` | pipeline en ligne de commande |
| `tools/`, `tests/`, `models/`, `docs/` | outils, tests, modèles de test, documentation |

## Usage de l'IA

On va le compléter 