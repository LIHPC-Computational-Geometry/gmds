#!/bin/bash

# Fichiers utilisés
SHAPES_LIST="shapes_list.txt"
OUTPUT_FILE="results.txt"

# Vérification que les fichiers requis existent
if [[ ! -f "$SHAPES_LIST" ]]; then
  echo "Le fichier $SHAPES_LIST est introuvable."
  exit 1
fi

# Initialisation du fichier de résultats
echo "Figure, Résultat, Temps (s)" > "$OUTPUT_FILE"

# Boucle sur chaque nom de figure dans le fichier shapes_list.txt
while IFS= read -r FIGURE; do
  # Ignorer les lignes vides
  [[ -z "$FIGURE" ]] && continue

  # Chemins vers les fichiers JSON et la base de données des figures (modifier si nécessaire)
  JSON_FILE="path/to/json_file.json"
  DATABASE_FILE="path/to/database.db"

  # Mesure du temps d'exécution
  START_TIME=$(date +%s.%N)

  # Exécution du programme C++
  ./main "$JSON_FILE" "$DATABASE_FILE" "$FIGURE"
  RETURN_CODE=$?

  # Calcul du temps écoulé
  END_TIME=$(date +%s.%N)
  ELAPSED_TIME=$(echo "$END_TIME - $START_TIME" | bc)

  # Traduction du code de retour en texte
  case $RETURN_CODE in
    0) RESULT="Défaite" ;;
    1) RESULT="Égalité" ;;
    2) RESULT="Victoire" ;;
    *) RESULT="Inconnu" ;;
  esac

  # Écriture des résultats dans le fichier
  echo "$FIGURE, $RESULT, $ELAPSED_TIME" >> "$OUTPUT_FILE"

done < "$SHAPES_LIST"

echo "Traitement terminé. Résultats sauvegardés dans $OUTPUT_FILE."