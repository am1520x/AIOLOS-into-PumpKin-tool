import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import traceback
import os
from plot_styles import (
    apply_plot_style, get_species_style, get_species_color, 
    get_species_linestyle, format_species_name
)

def plot_densities_and_mixing_ratios(
    species_file,
    densities_file,
    save_path="species_density_and_mixing_ratios.png"
):
    """
    Reads species names and densities, plots densities and mixing ratios vs. cell number on a log-log scale,
    and saves the plot.

    Args:
        species_file (str): Path to file with species names (one per line).
        densities_file (str): Path to file with densities. First col: cell, rest: densities.
        save_path (str): Output path to save the PNG plot.
    """
    try:
        # Read species names
        with open(species_file, 'r') as f:
            species = [line.strip() for line in f if line.strip()]
        expected_cols = len(species) + 1

        # Read density data (header row, then data)
        column_names = ['Cell'] + [f'Density_Sp{i+1}' for i in range(len(species))]
        dtype_dict = {'Cell': float}
        dtype_dict.update({f'Density_Sp{i+1}': float for i in range(len(species))})
        try:
            density_df = pd.read_csv(
                densities_file,
                sep='\s+',
                header=0,
                names=column_names,
                dtype=dtype_dict,
                engine='python'
            )
        except Exception as e:
            print(f"Error reading densities file: {e}")
            traceback.print_exc()
            return

        if 'Cell' not in density_df.columns:
            print("Error: 'Cell' column not found after reading.")
            return
        cells = density_df['Cell'].values

        density_cols = [col for col in density_df.columns if col.startswith('Density_Sp')]
        if len(density_cols) != len(species):
            print(f"Error: Number of density columns ({len(density_cols)}) does not match species list ({len(species)}).")
            return
        densities = density_df[density_cols].values

        # Calculate total density and mixing ratios safely
        total_density = np.sum(densities, axis=1)
        mixing_ratios = np.divide(
            densities,
            total_density[:, None],
            out=np.zeros_like(densities),
            where=total_density[:, None] > 0
        )

        # Apply consistent styling
        apply_plot_style()
        
        # Create plots
        fig, axs = plt.subplots(1, 2, figsize=(16, 6))

        # Plot densities
        ax1 = axs[0]
        for i, sp in enumerate(species):
            if np.any(densities[:, i] > 0):
                style = get_species_style(sp)
                formatted_sp = format_species_name(sp)
                ax1.plot(cells, densities[:, i], label=formatted_sp,
                        color=style['color'], linestyle=style['linestyle'], linewidth=2)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel('Cell Number')
        ax1.set_ylabel('Number Density')
        ax1.set_title('Species Number Density vs. Cell Number', fontweight='bold')
        ax1.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)

        # Plot mixing ratios
        ax2 = axs[1]
        for i, sp in enumerate(species):
            if np.any(mixing_ratios[:, i] > 0):
                style = get_species_style(sp)
                formatted_sp = format_species_name(sp)
                ax2.plot(cells, mixing_ratios[:, i], label=formatted_sp,
                        color=style['color'], linestyle=style['linestyle'], linewidth=2)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel('Cell Number')
        ax2.set_ylabel('Mixing Ratio')
        ax2.set_title('Species Mixing Ratio vs. Cell Number', fontweight='bold')
        ax2.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.3)

        plt.tight_layout()
        plt.subplots_adjust(right=0.8)
        print(f"Saving plot to {save_path}")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("Plot saved successfully.")

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        traceback.print_exc()

def main(
    species_file="qt_species_list.txt",
    densities_file="qt_densities.txt",
    output_plot="species_density_and_mixing_ratios.png"
):
    """
    Run the densities/mixing ratio plot and save to a file.

    Args:
        species_file (str): Path to the species list file.
        densities_file (str): Path to the densities file.
        output_plot (str): Output image file name.
    """
    plot_densities_and_mixing_ratios(species_file, densities_file, save_path=output_plot)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot densities and mixing ratios, save to PNG.")
    parser.add_argument("--species", type=str, default="../qt_species_list.txt", help="Species list file")
    parser.add_argument("--densities", type=str, default="../qt_densities.txt", help="Densities file")
    parser.add_argument("--out", type=str, default="./species_density_and_mixing_ratios.png", help="Output plot PNG file")
    args = parser.parse_args()
    main(
        species_file=args.species,
        densities_file=args.densities,
        output_plot=args.out
    )
