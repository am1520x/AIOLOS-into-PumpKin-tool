import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import traceback
from glob import glob
from collections import defaultdict

DATA_DIR = ""  # path to directory containing all pumpkin_output_cell_*.txt files
N_CELLS = 233

deleted_contributions = defaultdict(lambda: {"production": [np.nan] * N_CELLS, "consumption": [np.nan] * N_CELLS})
dominant_pathways = defaultdict(list)
flagged_entries = []

species_list = []

def parse_deleted_percent_block(lines, cell_number):
    """
    Parse a block of table lines reporting deleted pathway contributions.

    Parameters
    ----------
    lines : list of str
        Table lines, each containing species name and percent deleted for production and consumption.
    cell_number : int
        Index of the current cell (0-based).

    Side effects
    ------------
    - Updates global `deleted_contributions`, `species_list`, and `flagged_entries`.

    Example
    -------
    >>> lines = [
    ...     "|   H2O   |   15.0%   |   3.2%   |",
    ...     "|   CO2   |   7.1%    |   12.0%  |"
    ... ]
    >>> parse_deleted_percent_block(lines, 0)
    >>> 'H2O' in species_list
    True
    >>> deleted_contributions['CO2']['production'][0] == 7.1
    True
    >>> flagged_entries[-1][1]  # last flagged species
    'CO2'
    """
    pattern = r"\|\s+(.+?)\s+\|\s+([-\deE.+nan% ]+)\s+\|\s+([-\deE.+nan% ]+)\s+\|"
    for line in lines:
        match = re.match(pattern, line)
        if match:
            species = match.group(1).strip()
            prod_str = match.group(2).strip().replace('%', '')
            cons_str = match.group(3).strip().replace('%', '')

            try:
                prod = float(prod_str)
            except:
                prod = np.nan
            try:
                cons = float(cons_str)
            except:
                cons = np.nan

            if species not in species_list:
                species_list.append(species)
                if len(deleted_contributions[species]["production"]) != N_CELLS:
                     deleted_contributions[species]["production"] = [np.nan] * N_CELLS
                if len(deleted_contributions[species]["consumption"]) != N_CELLS:
                     deleted_contributions[species]["consumption"] = [np.nan] * N_CELLS

            if 0 <= cell_number < N_CELLS:
                deleted_contributions[species]["production"][cell_number] = prod
                deleted_contributions[species]["consumption"][cell_number] = cons
            else:
                print(f"Warning: Cell number {cell_number} is outside the expected range [0, {N_CELLS-1}]. Data not stored.")

            if (prod is not np.nan and prod > 10) or (cons is not np.nan and cons > 10):
                flagged_entries.append((cell_number, species, prod, cons))


def parse_file(filepath, cell_number):
    """
    Parse a single cell's output file for deleted pathway contributions.

    Parameters
    ----------
    filepath : str
        Path to the pumpkin_output_cell_XXX.txt file.
    cell_number : int
        The index for this cell (0-based).

    Side effects
    ------------
    - Calls `parse_deleted_percent_block` to fill global state.

    Example
    -------
    >>> import tempfile
    >>> with tempfile.NamedTemporaryFile("w+", delete=False) as tf:
    ...     _ = tf.write('''
    ... some header
    ... Procent of species from deleted pathways
    ... +----------+----------+----------+
    ... |   H2O    |   15.0%  |   2.1%   |
    ... |   CO2    |   7.1%   |   12.0%  |
    ... +----------+----------+----------+
    ... ''')
    ...     tf.flush()
    ...     parse_file(tf.name, 1)
    >>> deleted_contributions['H2O']['production'][1] == 15.0
    True
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        header_idx = None
        for idx, line in enumerate(lines):
            if "Procent of species from deleted pathways" in line:
                header_idx = idx
                break

        if header_idx is None:
            print(f"Section 'Procent of species from deleted pathways' not found in {filepath}. Skipping block parsing for this file.")
            return

        table_lines = []
        i = header_idx + 1
        sep_pattern = re.compile(r"^\+[-+]+\+\s*$")
        data_pattern = re.compile(r"^\|\s*.+\s*\|\s*[-\deE.+nan% ]+\s*\|\s*[-\deE.+nan% ]+\s*\|")

        while i < len(lines):
            line = lines[i].rstrip('\n')
            if sep_pattern.match(line) or data_pattern.match(line):
                table_lines.append(line)
                i += 1
            else:
                break

        data_rows = [line for line in table_lines if data_pattern.match(line)]
        parse_deleted_percent_block(data_rows, cell_number)

    except FileNotFoundError:
        print(f"File not found: {filepath}. Skipping cell {cell_number}.")
    except Exception as e:
        print(f"Error processing file {filepath}: {e}")
        traceback.print_exc()


def plot_deleted_contributions(save_path=None, show=True):
    """
    Plot heatmaps of production/consumption % from deleted pathways for each species vs cell.

    Side effects
    ------------
    - Plots two seaborn heatmaps (production and consumption) using global state.
    - Prints shapes of matrices for debugging.

    Returns
    -------
    None

    Example
    -------
    >>> # After running main() or populating deleted_contributions/species_list
    >>> plot_deleted_contributions()  # doctest: +SKIP
    """
    species_order = sorted(species_list)

    try:
        prod_list = [deleted_contributions[sp]["production"] for sp in species_order]
        cons_list = [deleted_contributions[sp]["consumption"] for sp in species_order]
        if not all(len(lst) == N_CELLS for lst in prod_list):
            print("Error: Inconsistent list lengths for production data.")
            return
        if not all(len(lst) == N_CELLS for lst in cons_list):
            print("Error: Inconsistent list lengths for consumption data.")
            return
        prod_matrix = np.array(prod_list)
        cons_matrix = np.array(cons_list)
        print(f"Shape of prod_matrix: {prod_matrix.shape}")
        print(f"Shape of cons_matrix: {cons_matrix.shape}")
    except Exception as e:
        print(f"Error creating numpy arrays: {e}")
        traceback.print_exc()
        return

    fig, axs = plt.subplots(2, 1, figsize=(16, 12), sharex=True)
    cell_labels = np.arange(N_CELLS)
    sns.heatmap(prod_matrix, xticklabels=10, yticklabels=species_order, cmap="Reds", ax=axs[0], cbar_kws={'label': 'Production %'})
    axs[0].set_title("Production % from Deleted Pathways")
    axs[0].set_ylabel("Species")
    axs[0].set_xlim(0, N_CELLS)
    sns.heatmap(cons_matrix, xticklabels=10, yticklabels=species_order, cmap="Blues", ax=axs[1], cbar_kws={'label': 'Consumption %'})
    axs[1].set_title("Consumption % from Deleted Pathways")
    axs[1].set_xlabel("Radial Cell Index")
    axs[1].set_ylabel("Species")
    axs[1].set_xlim(0, N_CELLS)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to {save_path}")
    if show:
        plt.show()
    else:
        plt.close()

def main():
    """
    Main driver function. Parses all cell files and plots results.
    Make sure you are in the directory folder of the data not this code, and run 'python Outputs/deleted_pathways_plots.py'.

    Returns
    -------
    None

    Example
    -------
    >>> main()  # doctest: +SKIP
    """
    print("Parsing files...")
    global species_list
    species_list = []

    for i in range(N_CELLS):
        filename = os.path.join(DATA_DIR, f"pumpkin_output_cell_{i:03d}.txt")
        parse_file(filename, i)

    print(f"\nFlagged Entries (Deleted Contribution > 10%):")
    for cell, species, prod, cons in flagged_entries:
        print(f"Cell {cell:03d} | Species: {species} | Production: {prod:.2f}% | Consumption: {cons:.2f}%")

    print("\nGenerating plots...")
    if species_list:
        plot_deleted_contributions(save_path="deleted_contributions.png", show=False)
        print("Plots generated and saved as 'deleted_contributions.png'.")
    else:
        print("No species data found to plot deleted contributions.")

if __name__ == "__main__":
    main()
