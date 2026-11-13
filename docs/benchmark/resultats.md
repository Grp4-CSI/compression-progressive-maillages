# Résultats sur les 8 modèles du site d'évaluation

Flux compressé `.pmc`, quantification 12 bits. Débits en bits par sommet du modèle
original ; distorsions (Hausdorff) relatives à la diagonale de la boîte englobante,
mesurées après avoir reçu 10 %, 25 % et 50 % du flux. Tous les modèles sont décodés
sans perte (à la quantification près). Produit par `tools/benchmark.py`.

| modèle | sommets | schéma | niveaux | M0 | connect. | géom. | total | Hausdorff 10 % | 25 % | 50 % | codage |
|---|---|---|---|---|---|---|---|---|---|---|---|
| icosphere | 642 | 4 | 12 | 13 | 12.52 | 49.07 | 64.65 | 0.0490 | 0.0177 | 0.0126 | 0.1 s |
| icosphere | 642 | auto | 15 | 12 | 9.61 | 54.77 | 67.30 | 0.0913 | 0.0335 | 0.0178 | 0.3 s |
| sphere | 2342 | 4 | 17 | 13 | 12.75 | 32.53 | 46.12 | 0.0216 | 0.0118 | 0.0049 | 0.5 s |
| sphere | 2342 | auto | 19 | 15 | 8.41 | 38.36 | 47.69 | 0.0407 | 0.0250 | 0.0092 | 1.2 s |
| bunny | 2503 | 4 | 24 | 43 | 12.67 | 34.14 | 48.44 | 0.0775 | 0.0396 | 0.0118 | 0.5 s |
| bunny | 2503 | auto | 31 | 44 | 9.98 | 40.46 | 52.12 | 0.0897 | 0.0303 | 0.0249 | 1.4 s |
| cow | 2903 | 4 | 22 | 15 | 12.66 | 29.14 | 42.54 | 0.0847 | 0.0313 | 0.0181 | 0.6 s |
| cow | 2903 | auto | 24 | 13 | 9.57 | 32.93 | 43.18 | 0.0605 | 0.0238 | 0.0237 | 1.7 s |
| suzanne | 7830 | 4 | 22 | 182 | 11.84 | 21.72 | 35.18 | 0.0465 | 0.0288 | 0.0072 | 1.7 s |
| suzanne | 7830 | auto | 30 | 181 | 8.02 | 24.57 | 34.18 | 0.0939 | 0.0380 | 0.0149 | 4.0 s |
| fandisk | 9864 | 4 | 22 | 13 | 11.31 | 20.34 | 31.85 | 0.0370 | 0.0139 | 0.0092 | 2.3 s |
| fandisk | 9864 | auto | 27 | 16 | 7.68 | 23.09 | 30.99 | 0.0414 | 0.0208 | 0.0119 | 5.3 s |
| pokemon | 25926 | 4 | 22 | 794 | 10.52 | 20.97 | 34.22 | 0.0731 | 0.0264 | 0.0070 | 8.1 s |
| pokemon | 25926 | auto | 24 | 767 | 8.36 | 22.66 | 33.64 | 0.0775 | 0.0443 | 0.0135 | 17.7 s |
| hippo | 32144 | 4 | 26 | 252 | 8.54 | 11.09 | 20.26 | 0.0240 | 0.0107 | 0.0054 | 9.2 s |
| hippo | 32144 | auto | 31 | 244 | 6.91 | 12.13 | 19.66 | 0.0394 | 0.0229 | 0.0096 | 18.6 s |

Lecture rapide :
- le mode auto baisse la connectivité de 19 à 34 % sur tous les modèles, mais augmente
  la géométrie de 8 à 19 % ; le total baisse sur les 4 plus gros modèles (−2 à −3 %)
  et augmente sur les 4 petits (+2 à +8 %) ;
- à fraction de flux égale, le 4-coloriage a la plus faible distorsion dans 21 mesures
  sur 24 (exceptions : la vache à 10 et 25 %, le lapin à 25 %) : c'est le meilleur choix
  pour la progressivité, d'où le choix par défaut ;
- plus le modèle est gros et régulier, moins il coûte par sommet (hippo : 20 b/s ;
  icosphere : 65 b/s, où l'en-tête et les grands écarts des niveaux grossiers pèsent lourd).
