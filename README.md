# AIOLOS-into-PumpKin-tool
This tool takes the simulation outputs from the AIOLOS code and converts it into inputs for use on the chemical reaction network analysis code PumpKin.

To set up the tool do the following:
1) Fork and Clone from this github page.
2) Move to the aiolos_pumpkin_tool dir.
3) setup virtual environment, e.g.
python -m venv .venv
source .venv/bin/activate # Windows: .venv\Scripts\activate or .venv\Scripts\activate.bat
4) pip install -e .
Now it should be ready to go, test with;
python -m aiolos_pumpkin_tool.main --help


Full Instructions to run:
Only seems to work on linux or wsl because of the syntax.
Go to the root directory of the code base. Then setup virtual environment. 
/mnt/d/OneDrive/Water Worlds/AIOLOS-into-PumpKin-tool$ 
python -m venv .venv

.venv/Scripts/activate

pip install -e .

python -m aiolos_pumpkin_tool.main --process-data --logfile log_HD40307_10Fearth_transcells2_long.log --chemfile chemistry_HD40307_10Fearth_transcells2_long_t19.dat --simulation-dir "../ozone/" --simulation-name HD40307_10Fearth_transcells2_long --timestep 19 

python -m aiolos_pumpkin_tool.main --run-pumpkin --pumpkin-dir "/mnt/d/OneDrive/Water Worlds/PumpKin/src" --pumpkin_in-dir "/mnt/d/OneDrive/Water Worlds/PumpKin/src/Examples/AIOLOS_New"

 python -m aiolos_pumpkin_tool.main --generate-plots --pumpkin-dir "/mnt/d/OneDrive/Water Worlds/PumpKin/src" --pumpkin_in-dir "/mnt/d/OneDrive/Water Worlds/PumpKin/src/Examples/AIOLOS_New"

 If doing running all will need all of the arguments, so becomes a long prompt. 

