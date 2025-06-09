import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import traceback
import os

def plot_reaction_rates(
    reactions_file,
    rates_file,
    save_path="reaction_rates.png"
):
    """
    Reads reaction names and rates, plots reaction rates vs cell number, and saves the plot.

    Args:
        reactions_file (str): Path to file with reaction names (one per line).
        rates_file (str): Path to file with reaction rates. First col is cell number.
        save_path (str): Where to save the output PNG plot.
    """
    try:
        print(f"Reading reactions from: {reactions_file}")
        with open(reactions_file, 'r') as f:
            reactions = [line.strip() for line in f if line.strip()]
        print(f"Found {len(reactions)} reactions.")

        print(f"Reading rates from: {rates_file}")
        try:
            rate_df = pd.read_csv(
                rates_file,
                sep='\s+',
                header=0,
                dtype=float,
                engine='python'
            )
            print(f"Successfully read rates file using read_csv with header=0. Shape: {rate_df.shape}")

        except Exception as e:
            print(f"Error reading rates file with pandas: {e}")
            print("Please check the file format, especially the header and the data columns.")
            traceback.print_exc()
            return

        if rate_df.shape[1] == len(reactions) + 1:
            cells = rate_df.iloc[:, 0].values
            rates_data = rate_df.iloc[:, 1:].values
        elif rate_df.shape[1] == len(reactions):
            try:
                cells = rate_df.index.astype(float).values
                rates_data = rate_df.values
            except ValueError:
                print("Error: Cannot determine cell numbers. Rates file should have cell numbers in header or first column.")
                return
        else:
            print(f"Error: Unexpected number of columns in rates file ({rate_df.shape[1]}). Expected {len(reactions)} or {len(reactions)+1}.")
            return

        rates = rates_data.T
        if rates.shape[0] != len(reactions):
            print(f"Error: Number of transposed rates rows ({rates.shape[0]}) does not match the number of reactions ({len(reactions)}).")
            return

        # Create plot
        fig, axs = plt.subplots(2, 1, figsize=(14, 16), sharex=True)

        # Main (all rates) plot
        reactions_to_plot_indices = [
            i for i in range(rates.shape[0])
            if np.any(rates[i, :] > 0)
        ]

        ax = axs[0]
        for i in reactions_to_plot_indices:
            reaction_name = reactions[i]
            reaction_rates_across_cells = rates[i, :]
            ax.plot(cells, reaction_rates_across_cells, label=reaction_name, linewidth=0.8)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Cell Number')
        ax.set_ylabel('Reaction Rate')
        ax.set_title('Reaction Rate vs. Cell Number')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend(title='Reaction', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small', ncol=1)

        # Zoomed-in subplot (thresholded)
        ax_zoomed = axs[1]
        zoom_threshold = 1e-12
        for i in reactions_to_plot_indices:
            reaction_name = reactions[i]
            reaction_rates_across_cells = rates[i, :]
            rates_filtered = np.where(reaction_rates_across_cells > zoom_threshold, reaction_rates_across_cells, np.nan)
            ax_zoomed.plot(cells, rates_filtered, label=reaction_name, linewidth=0.8)
        ax_zoomed.set_xscale('log')
        ax_zoomed.set_yscale('log')
        ax_zoomed.set_xlabel('Cell Number')
        ax_zoomed.set_ylabel('Reaction Rate')
        ax_zoomed.set_title(f'Reaction Rate vs. Cell Number (Rates > {zoom_threshold:.0e})')
        ax_zoomed.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax_zoomed.set_ylim(bottom=zoom_threshold)

        plt.tight_layout()
        plt.subplots_adjust(right=0.8)  # Make space for legend

        # Save the plot
        print(f"Saving plot to {save_path}")
        plt.savefig(save_path, bbox_inches='tight')
        plt.close(fig)
        print("Plot saved successfully.")

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        traceback.print_exc()

def main(
    reactions_path="qt_reactions_list.txt",
    rates_path="qt_rates.txt",
    output_plot="reaction_rates.png"
):
    """
    Run the reaction rate plotting and save the output.

    Args:
        reactions_path (str): Path to the reactions list file.
        rates_path (str): Path to the rates file.
        output_plot (str): Output image filename.
    """
    plot_reaction_rates(reactions_path, rates_path, output_plot)

# --- If run as a script, use command-line args or defaults ---
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot reaction rates and save plot.")
    parser.add_argument("--reactions", type=str, default="../qt_reactions_list.txt", help="Path to reactions list")
    parser.add_argument("--rates", type=str, default="../qt_rates.txt", help="Path to rates file")
    parser.add_argument("--out", type=str, default="reaction_rates.png", help="Output plot image file")
    args = parser.parse_args()

    main(
        reactions_path=args.reactions,
        rates_path=args.rates,
        output_plot=args.out
    )
