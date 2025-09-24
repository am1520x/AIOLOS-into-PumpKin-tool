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