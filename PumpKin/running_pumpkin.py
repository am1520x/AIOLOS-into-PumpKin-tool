import subprocess

# Full shell command to run in WSL
command = 'wsl bash -c "cd \\"/mnt/d/OneDrive/Water Worlds/PumpKin/src\\" && ./PumpKin Examples/AIOLOS/"'

# Input you want to pass into the C program (1 to 22, then -1)
input_data = "\n".join(str(i) for i in range(1, 23)) + "\n-1\n"

# Run the command and pass input
result = subprocess.run(command, input=input_data, capture_output=True, text=True, shell=True)

# Save output to a .txt file
with open("pumpkin_output.txt", "w") as f:
    f.write(result.stdout)

print("Output saved to pumpkin_output.txt")
