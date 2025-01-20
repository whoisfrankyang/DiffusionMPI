import os
import re
import numpy as np
import argparse

def integrate_states(input_folder, output_file):
    """
    Integrates all state files (including initial_state.txt and final_state.txt)
    into a single file.

    Parameters:
        input_folder (str): Path to the folder containing state files.
        output_file (str): Path to the output integrated file.
    """
    files = sorted(
        [f for f in os.listdir(input_folder) if re.match(r"(initial_state|state_\d{3}|final_state)\.txt", f)]
    )
    if not files:
        raise FileNotFoundError(f"No state files found in {input_folder}")

    with open(output_file, "w") as out_f:
        # Write header
        first_file_path = os.path.join(input_folder, files[0])
        first_state = np.loadtxt(first_file_path)
        nx, ny = first_state.shape
        num_states = len(files)
        out_f.write(f"# Grid dimensions: {nx} x {ny}\n")
        out_f.write(f"# Number of timesteps: {num_states}\n\n")

        # Process each file
        for idx, file in enumerate(files):
            file_path = os.path.join(input_folder, file)
            state_data = np.loadtxt(file_path)

            # Write timestep and data to the output file
            out_f.write(f"# Timestep: {idx}\n")  # Use index for timestep
            np.savetxt(out_f, state_data, fmt="%.6f")
            out_f.write("\n")  # Separate timesteps with a blank line

    print(f"All states (including initial and final) combined into {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Integrate state files into a single file.")
    parser.add_argument("input_folder", type=str, help="Path to the folder containing state files.")
    parser.add_argument("output_file", type=str, help="Path to the output file.")
    args = parser.parse_args()

    integrate_states(args.input_folder, args.output_file)
