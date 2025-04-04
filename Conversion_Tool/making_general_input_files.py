

import re
import numpy as np
from processing_aiolos_reac_file import parse_reaction_data, process_reaction_line

try:
  with open('steamfull_step3.reac', 'r') as f:  # Replace 'your_file_path.py' with the actual path
    reac_text = f.read()
except FileNotFoundError:
  print("File not found.")
  reac_text = None
  
if reac_text:
    reaction_list, stoich_list, species_list = parse_reaction_data(reac_text)

    # --- File Output ---
    # Write species list
    with open("qt_species_list.txt", "w") as f:
        for i, sp in enumerate(species_list, start=1):
            f.write(f"{i} {sp}\n")

    # Write reactions list (modified format)
    with open("qt_reactions_list.txt", "w") as f:
        for idx, reaction in enumerate(reaction_list, start=1):
            modified = process_reaction_line(reaction)
            f.write(f"{idx} {modified}\n")

    # Write stoichiometry matrix
    with open("qt_matrix.txt", "w") as f:
        for sp in species_list:
            row = [str(int(coeff)) if isinstance(coeff, (int, float)) and coeff.is_integer() else str(coeff) for coeff in (r.get(sp, 0) for r in stoich_list)]
            f.write("\t".join(row) + "\n")

    print("Files generated: qt_species_list.txt, qt_reactions_list.txt, qt_matrix.txt")