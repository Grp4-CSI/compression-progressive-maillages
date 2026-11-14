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

### Visualisation niveau par niveau (sans serveur)

Le visualiseur OpenGL lit directement le flux compressé `.pmc` et affiche le maillage
après M0 puis après chaque niveau.

```
pip install glfw PyOpenGL
python main.py example/suzanne.obj --no-eval       # produit out/suzanne.pmc
python tools/view_pmc.py out/suzanne.pmc
```

| Touche | Effet |
|---|---|
| → / ← | niveau suivant / précédent |
| Début / Fin | maillage de base / maillage complet |
| C | surbrillance des faces ajoutées au niveau |
| W | afficher / masquer les arêtes |
| R | vue initiale |
| souris | clic gauche : tourner, clic droit : déplacer, molette : zoomer |
| Échap | quitter |

Le titre de la fenêtre indique le niveau, le nombre de sommets et de faces,
et la quantité du flux reçue (octets et bits par sommet).

## Outils
```
python -m pytest tests                              # 67 tests
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

## Utilisation de l’intelligence artificielle

L’IA a été utilisée comme outil d’aide à la compréhension, à la conception, au développement et à la présentation du projet. Ses propositions ont été examinées, testées et adaptées par le groupe.

### 1. Compréhension de l’article

L’IA nous a aidés à comprendre l’article de Cohen-Or, Levin et Remez (1999), notamment les étapes de compression progressive, les méthodes de coloration et la reconstruction du maillage. Elle nous a également permis d’identifier les choix laissés à notre implémentation, ainsi que les limites et hypothèses de la méthode.

**Exemple de prompt utilisé :**

> Nous travaillons sur une implémentation de l’article « Progressive Compression of Arbitrary Triangular Meshes » de Cohen-Or, Levin et Remez (1999). Aide-nous à comprendre précisément la méthode, étape par étape.
>
> Pour chaque mécanisme, distingue ce qui est explicitement décrit dans l’article de ce qui relève d’une hypothèse ou d’un choix d’implémentation. Explique également les contraintes topologiques, les limites de la méthode et les difficultés possibles lors de son implémentation.
>
> Ne présente pas les informations manquantes comme des faits : signale les ambiguïtés et indique les points qui doivent être vérifiés directement dans l’article.

### 2. Réflexion sur les choix techniques

L’IA nous a aidés à comparer et à questionner plusieurs choix d’implémentation :
- **Sélection des sommets indépendants :** réflexion sur les stratégies de parcours et la prise en compte de la valence des sommets.
- **Ordre de reconstruction :** comparaison entre un parcours en largeur et un tri canonique, notamment pour garantir un ordre déterministe.
- **Coloration des patches :** comparaison des méthodes à 2 couleurs, 4 couleurs et du choix automatique. Les résultats obtenus nous ont conduits à privilégier la méthode à 4 couleurs.
- **Robustesse :** proposition de vérifications supplémentaires. Certaines étant trop restrictives sur les modèles du site, nous avons préféré signaler les défauts et mettre à l’écart les faces non-manifold sans bloquer systématiquement le traitement.

**Exemple de prompt utilisé :**

> Nous implémentons une méthode de compression progressive de maillages triangulaires. Nous devons faire plusieurs choix qui ne sont pas nécessairement imposés par l’article : stratégie de parcours des sommets indépendants, ordre déterministe de reconstruction et méthode de coloration des patches (2 couleurs, 4 couleurs ou choix automatique).
>
> Compare les différentes possibilités en expliquant leurs avantages, leurs inconvénients et leurs conséquences sur la compression, la reconstruction et la robustesse.
>
> Ne choisis pas automatiquement une solution à notre place. Présente les critères permettant de trancher et les expériences à réaliser pour vérifier les hypothèses. Distingue les conclusions théoriques des points qui nécessitent des mesures sur nos propres modèles.


### 3. Tests

Une partie des tests unitaires et d’intégration a été générée avec l’aide de l’IA. Au début du projet, nous avons utilisé des modèles OBJ de référence trouvés en ligne, les modèles officiels du site n’étant pas encore disponibles.

Les tests proposés ont été relus, vérifiés, exécutés et adaptés par le groupe. La validation finale a ensuite été réalisée sur les huit modèles fournis par le site.

**Exemple de prompt utilisé :**

> Nous développons un encodeur et un décodeur pour une méthode de compression progressive de maillages triangulaires. Propose une série de tests unitaires et d’intégration couvrant la sélection de sommets indépendants, l’identification des patches, la coloration, la suppression et la réinsertion des sommets, la conservation de la connectivité et la lecture du flux compressé.
>
> Pour chaque test, précise le comportement attendu, les propriétés à vérifier et les cas limites à envisager.
>
> N’invente pas de résultats d’exécution. Signale les hypothèses nécessaires et privilégie des tests reproductibles sur de petits modèles OBJ dont les résultats peuvent être vérifiés. Le code proposé doit être adapté à notre architecture existante.

### 4. Aide au développement

L’IA nous a aidés à améliorer la structure du code, notamment les signatures de fonctions, les annotations de types, les commentaires et la répartition des responsabilités entre modules.

Elle nous a également orientés vers certaines bibliothèques existantes. Nous avons notamment étudié OpenMesh, sans le retenir, car notre structure half-edge interne répondait déjà à nos besoins. Des bibliothèques spécialisées en compression ont aussi été envisagées.

### 5. Visualisation OpenGL

L’IA nous a accompagnés dans la création du visualiseur OpenGL, notamment pour la prise en main des buffers, des shaders, de la caméra et du rendu des triangles.

Le visualiseur exploite les états réellement produits par notre décodeur afin de montrer la reconstruction progressive du maillage, les faces ajoutées et l’évolution du nombre d’octets reçus.

**Exemple de prompt utilisé :**

> Nous souhaitons développer un visualiseur OpenGL pour notre projet de compression progressive de maillages. Il doit afficher les états successifs réellement produits par notre décodeur, afin de visualiser l’ajout progressif des faces et l’évolution du nombre d’octets reçus.
>
> Aide-nous à comprendre les éléments nécessaires en OpenGL (buffers, shaders, caméra, transformations et rendu des triangles) et propose une architecture adaptée.
>
> Le visualiseur ne doit pas implémenter sa propre méthode de simplification ou de reconstruction : il doit exploiter les données produites par notre décodeur. Tiens compte de notre architecture existante et explique les modifications nécessaires avant de proposer du code.

### 6. Illustrations de la présentation

Plusieurs schémas de la présentation ont été générés avec l’aide de l’IA pour illustrer les étapes de décimation, les patches, les méthodes de coloration, la prédiction géométrique et la réinsertion des sommets.

### 7. Vérification des propositions

Les réponses de l’IA ont servi de pistes de réflexion et non de références systématiques. Les choix finaux ont été confrontés à l’article, à notre implémentation, aux tests et aux résultats expérimentaux.
