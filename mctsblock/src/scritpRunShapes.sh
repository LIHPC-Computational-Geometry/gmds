#!/bin/bash

# Files used
SHAPES_LIST="shapes_list.txt"
OUTPUT_FILE="results.txt"
CPP_FILE="main_blocker.cpp"
EXECUTABLE="main_blocker"

# Check if the C++ source file exists
if [[ ! -f "$CPP_FILE" ]]; then
  echo "The source file $CPP_FILE was not found."
  exit 1
fi

# Compile the C++ program
echo "Compiling $CPP_FILE..."
g++ -o "$EXECUTABLE" "$CPP_FILE" -std=c++17
if [[ $? -ne 0 ]]; then
  echo "Error while compiling $CPP_FILE."
  exit 1
fi
echo "Compilation successful."

# Check if the shapes list file exists
if [[ ! -f "$SHAPES_LIST" ]]; then
  echo "The file $SHAPES_LIST was not found."
  exit 1
fi

# Initialize the results file
echo "SHAPE NAME, RESULT, TIME(s)" > "$OUTPUT_FILE"

# Loop through each shape name in the shapes_list.txt file
while IFS= read -r SHAPE; do
  # Skip empty lines
  [[ -z "$SHAPE" ]] && continue

  # Paths to the JSON file and the shape database
  JSON_FILE="../../mctsblock/tst/data/params.json"
  DATABASE_FILE="/Users/paulbourmaud/Projects/blocking-learning-data/input/${SHAPE}.vtk"

  # Measure execution time
  START_TIME=$(date +%s.%N)

  # Execute the compiled program
  ./"$EXECUTABLE" "$JSON_FILE" "$DATABASE_FILE" "$SHAPE"
  RETURN_CODE=$?

  # Calculate elapsed time
  END_TIME=$(date +%s.%N)
  ELAPSED_TIME=$(echo "$END_TIME - $START_TIME" | bc)

  # Translate the return code into a result
  case $RETURN_CODE in
    0) RESULT="Defeat" ;;
    1) RESULT="Draw" ;;
    2) RESULT="Win" ;;
    *) RESULT="Unknown" ;;
  esac

  # Write the results to the file
  echo "$SHAPE, $RESULT, $ELAPSED_TIME" >> "$OUTPUT_FILE"

done < "$SHAPES_LIST"

echo "Processing completed. Results saved in $OUTPUT_FILE."