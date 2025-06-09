import subprocess
import re

# Number of spatial cells
num_cells = 233

# Path to WSL PumpKin directory and input.txt file
pumpkin_dir = '/mnt/d/OneDrive/Water Worlds/PumpKin/src/Examples/AIOLOS_New'
input_file_path = f'{pumpkin_dir}/input.txt'

# Path to executable
#command_base = f'"cd \\"/mnt/d/OneDrive/Water Worlds/PumpKin/src\\" && ./PumpKin Examples/AIOLOS_New/"'
command = 'cd "/mnt/d/OneDrive/Water Worlds/PumpKin/src" && ./PumpKin Examples/AIOLOS_New/'

# Input species list (1 to 26) + end flag -1
input_data = "\n".join(str(i) for i in range(1, 27)) + "\n-1\n"

def modify_input_file(cell_index):
    """Updates t_init and t_end in input.txt for a given cell."""
    with open(input_file_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if line.strip().startswith("t_init"):
            new_lines.append(f"t_init = {cell_index}\n")
        elif line.strip().startswith("t_end"):
            new_lines.append(f"t_end = {cell_index + 1}\n")
        else:
            new_lines.append(line)

    with open(input_file_path, "w") as f:
        f.writelines(new_lines)

# Loop through all cells
for cell_index in range(num_cells):
    print(f"Running PumpKin for cell {cell_index}...")

    # Step 1: Update input.txt with current cell
    modify_input_file(cell_index)

    # Step 2: Run PumpKin with current input
    result = subprocess.run(command, input=input_data, capture_output=True, text=True, shell=True)
    if result.returncode != 0:
        print(f"PumpKin exited with error code {result.returncode} at cell {cell_index}")

    # Combine stdout and stderr for debugging
    output_combined = result.stdout + "\n--- STDERR ---\n" + result.stderr

    # Step 3: Save unique output
    output_filename = f"pumpkin_output_cell_{cell_index:03}.txt"
    with open(output_filename, "w") as f:
        f.write(output_combined)

print("All cells processed and outputs saved.")
