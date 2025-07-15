from .processing_aiolos_reac_file import (
    parse_reaction_data,
    process_reaction_line,
    transform_species_set,
)
from .making_densities_file import process_timesteps
from .making_rates_file import make_rates

try:
    with open(
        "thermo.reac", "r"
    ) as f:  # Replace 'your_file_path.py' with the actual path
        reac_text = f.read()
except FileNotFoundError:
    print("File not found.")
    reac_text = None

if reac_text:
    reaction_list, stoich_list, species_list = parse_reaction_data(reac_text)
    species_list = transform_species_set(species_list)

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
            row = [
                (
                    str(int(coeff))
                    if isinstance(coeff, (int, float)) and coeff.is_integer()
                    else str(coeff)
                )
                for coeff in (r.get(sp, 0) for r in stoich_list)
            ]
            f.write("\t".join(row) + "\n")

    print("Files generated: qt_species_list.txt, qt_reactions_list.txt, qt_matrix.txt")

# Define your species (example format)
species = [
    ["H", 1, 1, r"$\rm H$"],
    ["H-", 1, -1, r"$\rm H^-$"],
    ["H2", 2, 0, r"$\rm H_2$"],
    ["H2O", 18, 0, r"$\rm H_2O$"],
    ["H2O2", 34, 0, r"$\rm H_2O_2$"],
    ["H2Op", 18, +1, r"$\rm H_2O^+$"],
    ["H2p", 2, +1, r"$\rm H_2^+$"],
    ["H3Op", 19, +1, r"$\rm H_3O^+$"],
    ["H3p", 3, +1, r"$\rm H_3^+$"],
    ["O", 16, 0, r"$\rm O$"],
    ["O-", 16, -1, r"$\rm O^-$"],
    ["O2", 32, 0, r"$\rm O_2$"],
    ["O2-", 32, -1, r"$\rm O_2^-$"],
    ["O2H", 33, 0, r"$\rm O_2H$"],
    ["O2Hp", 33, +1, r"$\rm O_2H^+$"],
    ["O2p", 32, +1, r"$\rm O_2^+$"],
    ["OH", 17, 0, r"$\rm OH$"],
    ["OH-", 17, -1, r"$\rm OH^-$"],
    ["OHp", 17, +1, r"$\rm OH^+$"],
    ["Op", 16, +1, r"$\rm O^+$"],
    ["S1", 1, +1, r"$\rm H$"],
    ["e-", 5e-4, -1, r"$\rm e^-$"],
]

# Define timesteps to process
timesteps = range(13, 17, 1)  # From 0 to 100 in steps of 5

# Makes the densities file and also stores the avgerage temperature and \n# species number densities at a given radial index to calculate rates
avg_T_at_index, number_densities = process_timesteps(
    directory="../dynamic_cond0_data/",
    sim="dynamic_ignoreleectrontest-cell200-newsol2-h2e+2-long-cond0",
    timesteps=timesteps,
    species=species,
    rplanet=1.37e8,  # Planet radius in cm
    index=10,  # Extract data at radial index 10 can be changed
    output_file="species_densities_at_r10.txt",
)

# Make the rates file
make_rates(
    output_file="rates.txt",
    avg_T_at_index=avg_T_at_index,
    number_densities=number_densities,
    timesteps=timesteps,
)
