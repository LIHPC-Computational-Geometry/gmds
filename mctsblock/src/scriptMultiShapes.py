import os
import subprocess
import time
import csv

# Files used
SHAPES_LIST = "shapes_list.txt"
OUTPUT_FILE = "results.txt"
CPP_FILE = "main_blocker.cpp"
EXECUTABLE = "main_blocker"

# Function to compile the C++ program
def compile_cpp_program():
    if not os.path.exists(CPP_FILE):
        print(f"The source file {CPP_FILE} was not found.")
        return False

    print(f"Compiling {CPP_FILE}...")
    compile_command = ["g++", "-o", EXECUTABLE, CPP_FILE, "-std=c++17"]
    result = subprocess.run(compile_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        print(f"Error while compiling {CPP_FILE}:\n{result.stderr.decode()}")
        return False

    print("Compilation successful.")
    return True

# Function to process shapes
def process_shapes():
    if not os.path.exists(SHAPES_LIST):
        print(f"The file {SHAPES_LIST} was not found.")
        return

    # Initialize the results file
    with open(OUTPUT_FILE, mode="w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(["SHAPE NAME", "RESULT", "TIME(s)"])  # CSV header

        # Read shapes from the file
        with open(SHAPES_LIST, mode="r") as shapes_file:
            for line in shapes_file:
                shape = line.strip()
                if not shape:  # Skip empty lines
                    continue

                # Paths to required files
                json_file = "../../mctsblock/tst/data/params.json"
                database_file = f"/Users/paulbourmaud/Projects/blocking-learning-data/input/{shape}.vtk"

                # Measure execution time
                start_time = time.time()

                # Execute the compiled program
                command = [f"./{EXECUTABLE}", json_file, database_file, shape]
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                elapsed_time = time.time() - start_time

                # Map return codes to results
                if result.returncode == 0:
                    result_text = "Defeat"
                elif result.returncode == 1:
                    result_text = "Draw"
                elif result.returncode == 2:
                    result_text = "Win"
                else:
                    result_text = "Unknown"

                # Write result to the CSV file
                writer.writerow([shape, result_text, round(elapsed_time, 3)])

    print(f"Processing completed. Results saved in {OUTPUT_FILE}.")

# Main script logic
if __name__ == "__main__":
    if compile_cpp_program():
        process_shapes()