#!/bin/bash

# Fichiers utilisés
SHAPES_LIST="shapes_list.txt"
OUTPUT_FILE="results.txt"

# Vérification que les fichiers requis existent
if [[ ! -f "$SHAPES_LIST" ]]; then
  echo "The file $SHAPES_LIST is not found."
  exit 1
fi

# Initialisation du fichier de résultats
echo "SHAPE NAME, RESULT, TIME(s)" > "$OUTPUT_FILE"

# Boucle sur chaque nom de figure dans le fichier shapes_list.txt
while IFS= read -r FIGURE; do
  # Ignorer les lignes vides
  [[ -z "$FIGURE" ]] && continue

  # Chemins vers les fichiers JSON et la base de données des figures (modifier si nécessaire)
  JSON_FILE="../../mctsblock/tst/data/params.json"
  DATABASE_FILE="/Users/paulbourmaud/Projects/blocking-learning-data/input/"+FIGURE+".vtk"

  # Mesure du temps d'exécutiong
  START_TIME=$(date +%s.%N)

  # Exécution du programme C++
  ./main_blocker.cpp "$JSON_FILE" "$DATABASE_FILE" "$FIGURE"
  RETURN_CODE=$?

  # Calcul du temps écoulé
  END_TIME=$(date +%s.%N)
  ELAPSED_TIME=$(echo "$END_TIME - $START_TIME" | bc)

  # Traduction du code de retour en texte
  case $RETURN_CODE in
    0) RESULT="Defeat" ;;
    1) RESULT="Draw" ;;
    2) RESULT="Win" ;;
    *) RESULT="Inconnu" ;;
  esac

  # Écriture des résultats dans le fichier
  echo "$FIGURE, $RESULT, $ELAPSED_TIME" >> "$OUTPUT_FILE"

done < "$SHAPES_LIST"

echo "Traitement terminé. Résultats sauvegardés dans $OUTPUT_FILE."