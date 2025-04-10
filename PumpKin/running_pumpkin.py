"""
This script automates the execution of the PumpKin C program via WSL (Windows Subsystem for Linux)
and passes a predefined list of interested species as input.

Steps performed:
1. Executes the PumpKin executable located in a specified WSL directory using a shell command.
2. Provides input to the program representing interested species (numbers 1 to 26, followed by -1 to signal the end).
3. Captures the standard output of the PumpKin run.
4. Saves the output to a text file named 'pumpkin_output.txt' for further use or inspection.

Note:
- Ensure WSL is installed and properly configured on your system.
- The path to the PumpKin executable and example directory must be correctly mapped to the WSL file system (/mnt/...).
- Adjust the `command` and input range if your use case differs.
"""
import subprocess

# Full shell command to run in WSL
command = 'wsl bash -c "cd \\"/mnt/d/OneDrive/Water Worlds/PumpKin/src\\" && ./PumpKin Examples/AIOLOS_New/"'

# The input we want to pass into PumpKin (1 to 26, then -1), i.e. all of the species as interested species and -1 to finish.
input_data = "\n".join(str(i) for i in range(1, 27)) + "\n-1\n"

# Run the command and pass input
result = subprocess.run(command, input=input_data, capture_output=True, text=True, shell=True)

# Save output to a .txt file
with open("pumpkin_output.txt", "w") as f:
    f.write(result.stdout)

# Print message to show code has run to completion
print("Output saved to pumpkin_output.txt")
