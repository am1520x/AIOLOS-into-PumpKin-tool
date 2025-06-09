import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def parse_all_species_pathway_tables_multiline(text):
    """
    Parses a multi-species pathway block into a DataFrame with columns:
    species, pathway, rate, production_pct, consumption_pct, net_pct.

    Parameters
    ----------
    text : str
        File content as a single string.

    Returns
    -------
    pd.DataFrame

    Example
    -------
    >>> test_block = '''
    ... | Pathway | Production of H | Consumption of H |
    ... +---------------------------------------------------+---------------------------+---------------------------+
    ... | 1 * (H+H2O=>OH+H2)             5.2465e+08 |                       0 % |                      91 % |
    ... +---------------------------------------------------+---------------------------+---------------------------+
    ... # End block
    ... '''
    >>> df = parse_all_species_pathway_tables_multiline(test_block)
    >>> df.iloc[0]['species']
    'H'
    >>> bool(abs(df.iloc[0]['rate']) == 5.2465e+08)  # No rate in this test
    True
    """
    lines = text.splitlines()
    results = []
    i = 0

    while i < len(lines) - 2:
        if "Pathway" in lines[i] and "Production of" in lines[i] and "Consumption of" in lines[i]:
            header_line = lines[i]
            species_match = re.search(r'Production of\s+([^\|\s]+)', header_line)
            if not species_match:
                i += 1
                continue

            species = species_match.group(1).strip()
            i += 1

            pathway_buffer = []
            while i < len(lines):
                line = lines[i].strip()

                if line.startswith('#'):
                    break

                if line.startswith('+'):
                    i += 1
                    continue

                if line.startswith('|'):
                    parts = [p.strip() for p in line.strip('|').split('|')]
                    if len(parts) != 3:
                        i += 1
                        continue

                    raw_pathway, prod_col, cons_col = parts

                    # Check if this is the final line of a pathway
                    rate_match = re.search(r"([\d\.eE\+\-]+)\s*$", raw_pathway)
                    if rate_match:
                        rate = float(rate_match.group(1))
                        cleaned_pathway = re.sub(r"([\d\.eE\+\-]+)\s*$", "", raw_pathway).strip()
                        if cleaned_pathway:
                            pathway_buffer.append(cleaned_pathway)

                        def clean_step(step):
                            step = re.sub(r"^\s*\d+\s*\*\s*", "", step)  # Remove '1 *' or similar
                            step = step.strip()
                            # Clean inside parentheses if present
                            step = re.sub(r'\((.*?)\)', lambda m: f"({m.group(1).strip()})", step)
                            return step

                        full_pathway = " -> ".join(clean_step(p) for p in pathway_buffer)

                        pathway_buffer = []

                        def parse_pct(p):
                            try:
                                return float(p.replace('%', '').replace('e', 'E'))
                            except:
                                return np.nan

                        prod_pct = parse_pct(prod_col)
                        cons_pct = parse_pct(cons_col)

                        results.append({
                            "species": species,
                            "pathway": full_pathway,
                            "rate": rate,
                            "production_pct": prod_pct,
                            "consumption_pct": cons_pct,
                            "net_pct": prod_pct - cons_pct if not np.isnan(prod_pct) and not np.isnan(cons_pct) else np.nan
                        })
                    else:
                        if raw_pathway:
                            pathway_buffer.append(raw_pathway.strip())
                i += 1

            if not any(r['species'] == species for r in results):
                results.append({
                    "species": species,
                    "pathway": np.nan,
                    "rate": np.nan,
                    "production_pct": np.nan,
                    "consumption_pct": np.nan,
                    "net_pct": np.nan
                })
        else:
            i += 1

    return pd.DataFrame(results)


def load_all_cells_as_dict(data_dir, file_pattern="pumpkin_output_cell_*.txt"):
    """
    Loads all cell text files in a directory into a dict of DataFrames.

    Parameters
    ----------
    data_dir : str
        Directory containing files.
    file_pattern : str
        Pattern for files. Default is 'pumpkin_output_cell_*.txt'.

    Returns
    -------
    dict
        Maps cell number (int) to pd.DataFrame for that cell.

    Example
    -------
    >>> import tempfile, os
    >>> tmpdir = tempfile.mkdtemp()
    >>> fname = os.path.join(tmpdir, "pumpkin_output_cell_004.txt")
    >>> with open(fname, "w") as f:
    ...     _ = f.write('| Pathway | Production of H | Consumption of H |\\n+----+----+----+\\n| 1 * (A=>B) (0) | 10% | 3% |\\n#')
    >>> d = load_all_cells_as_dict(tmpdir)
    >>> 4 in d
    True
    """
    cell_data = {}
    regex_pattern = file_pattern.replace("*", r"\d+")
    for filename in sorted(os.listdir(data_dir)):
        if re.match(regex_pattern, filename):
            cell_match = re.search(r'cell_(\d+)', filename)
            if not cell_match:
                continue
            cell_number = int(cell_match.group(1))
            filepath = os.path.join(data_dir, filename)
            with open(filepath, 'r') as f:
                text = f.read()
            df = parse_all_species_pathway_tables_multiline(text)
            cell_data[cell_number] = df
    return cell_data


def plot_top_species_pathways_by_rate(cell_data_dict, species_name, top_n=3, save_path=None):
    """
    Plots and optionally saves top N production/consumption pathways for a species.

    Parameters
    ----------
    cell_data_dict : dict
        Maps cell numbers to DataFrames.
    species_name : str
    top_n : int
    save_path : str or None
        If provided, saves figure here.

    Returns
    -------
    None
    """
    prod_pathway_data = defaultdict(lambda: [0.0] * len(cell_data_dict))
    cons_pathway_data = defaultdict(lambda: [0.0] * len(cell_data_dict))
    cell_numbers = sorted(cell_data_dict.keys())

    for idx, cell in enumerate(cell_numbers):
        df = cell_data_dict[cell]
        df_species = df[df['species'] == species_name]

        for _, row in df_species.iterrows():
            path = row['pathway']
            rate = row['rate']
            if pd.isna(path) or pd.isna(rate):
                continue
            if row['production_pct'] > row['consumption_pct']:
                prod_pathway_data[path][idx] = rate
            elif row['consumption_pct'] > row['production_pct']:
                cons_pathway_data[path][idx] = abs(rate)

    # Rank pathways by total contribution
    prod_ranked = sorted(prod_pathway_data.items(), key=lambda kv: sum(kv[1]), reverse=True)[:top_n]
    cons_ranked = sorted(cons_pathway_data.items(), key=lambda kv: sum(kv[1]), reverse=True)[:top_n]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharex=True, sharey=True)
    fig.suptitle(f"Top {top_n} Pathways Producing and Consuming {species_name}", fontsize=14)

    ax = axes[0]
    for path, rates in prod_ranked:
        x_vals = cell_numbers
        y_vals = rates
        if any(r > 0 for r in y_vals):
            ax.plot(x_vals, y_vals, marker='.', label=path)
    ax.set_title("Production Pathways")
    ax.set_xlabel("Cell Number")
    ax.set_ylabel("Rate")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8)

    ax = axes[1]
    for path, rates in cons_ranked:
        x_vals = cell_numbers
        y_vals = rates
        if any(r > 0 for r in y_vals):
            ax.plot(x_vals, y_vals, marker='.', label=path)
    ax.set_title("Consumption Pathways")
    ax.set_xlabel("Cell Number")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved: {save_path}")
    plt.close(fig)

def plot_species_net_rate(cell_data_dict, species_name, save_path=None):
    """
    Plots and optionally saves the net rate (prod-cons) per cell for a species.

    Parameters
    ----------
    cell_data_dict : dict
    species_name : str
    save_path : str or None

    Returns
    -------
    None
    """
    cell_numbers = sorted(cell_data_dict.keys())
    net_rates = []

    for cell in cell_numbers:
        df = cell_data_dict[cell]
        df_species = df[df['species'] == species_name]
        prod_total = df_species.loc[df_species['production_pct'] > df_species['consumption_pct'], 'rate'].sum()
        cons_total = df_species.loc[df_species['consumption_pct'] > df_species['production_pct'], 'rate'].sum()
        net_rate = prod_total - cons_total
        net_rates.append(net_rate)

    plt.figure(figsize=(10, 6))
    plt.plot(cell_numbers, net_rates, marker='o')
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    plt.xscale('log')
    plt.yscale('symlog', linthresh=1e-2)
    plt.xlabel("Cell Number (log)")
    plt.ylabel("Net Rate (Production − Consumption)")
    plt.title(f"Net Rate of {species_name} per Cell (Log-SymLog)")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved: {save_path}")
    plt.close()

def plot_pathway_heatmap(cell_data_dict, species_name, save_path=None):
    """
    Plots and optionally saves a pathway-cell heatmap for a species.

    Parameters
    ----------
    cell_data_dict : dict
    species_name : str
    save_path : str or None

    Returns
    -------
    None
    """
    all_data = []
    for cell, df in cell_data_dict.items():
        df_species = df[df['species'] == species_name]
        for _, row in df_species.iterrows():
            all_data.append({
                "cell": cell,
                "pathway": row["pathway"],
                "rate": row["rate"]
            })

    df_all = pd.DataFrame(all_data)
    if df_all.empty:
        print(f"No data found for species {species_name}")
        return

    pivot = df_all.pivot_table(index="pathway", columns="cell", values="rate", aggfunc="sum", fill_value=0)
    plt.figure(figsize=(12, max(4, 0.5 * len(pivot))))
    sns.heatmap(np.log10(pivot + 1e-5), annot=False, cmap="viridis", cbar_kws={'label': 'log10(rate)'})
    plt.title(f"Heatmap of Pathway Rates for {species_name}")
    plt.xlabel("Cell")
    plt.ylabel("Pathway")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved: {save_path}")
    plt.close()

def plot_stacked_production_fraction(cell_data_dict, species_name, top_n=3, save_path=None):
    """
    Plots and optionally saves a stacked area chart of top N production pathways' fraction for a species.

    Parameters
    ----------
    cell_data_dict : dict
    species_name : str
    top_n : int
    save_path : str or None

    Returns
    -------
    None
    """
    cell_numbers = sorted(cell_data_dict.keys())
    contrib_data = defaultdict(lambda: [0.0] * len(cell_numbers))
    for i, cell in enumerate(cell_numbers):
        df = cell_data_dict[cell]
        df_species = df[df['species'] == species_name]
        for _, row in df_species.iterrows():
            if row["production_pct"] > row["consumption_pct"]:
                contrib_data[row["pathway"]][i] += row["rate"]

    if not contrib_data:
        print(f"No production data found for species {species_name}")
        return

    top_pathways = sorted(contrib_data.items(), key=lambda x: sum(x[1]), reverse=True)[:top_n]
    df_plot = pd.DataFrame({p: r for p, r in top_pathways}, index=cell_numbers)
    df_plot_frac = df_plot.div(df_plot.sum(axis=1), axis=0).fillna(0)

    df_plot_frac.plot(kind='area', stacked=True, figsize=(12, 6))
    plt.title(f"Top {top_n} Production Pathway Fractions for {species_name}")
    plt.xlabel("Cell Number")
    plt.ylabel("Fraction of Production Rate")
    plt.ylim(0, 1)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved: {save_path}")
    plt.close()

def compute_pathway_species_matrix(cell_data_dict):
    """
    Computes a matrix (DataFrame) with pathways as rows and species as columns.
    1 if pathway involves species, else 0.

    Parameters
    ----------
    cell_data_dict : dict

    Returns
    -------
    pd.DataFrame

    Example
    -------
    >>> # See test file for examples.
    """
    species_pathways = defaultdict(set)
    for df in cell_data_dict.values():
        for _, row in df.iterrows():
            if not pd.isna(row["pathway"]):
                species_pathways[row["pathway"]].add(row["species"])
    all_pathways = sorted(species_pathways.keys())
    all_species = sorted(set(s for species_set in species_pathways.values() for s in species_set))
    matrix = pd.DataFrame(0, index=all_pathways, columns=all_species)
    for p, s_set in species_pathways.items():
        for s in s_set:
            matrix.at[p, s] = 1
    return matrix

def plot_pathway_species_heatmap(matrix, save_path=None):
    """
    Plots and optionally saves a pathway-species matrix heatmap.

    Parameters
    ----------
    matrix : pd.DataFrame
    save_path : str or None

    Returns
    -------
    None
    """
    if matrix.empty:
        print("Pathway-species matrix is empty.")
        return
    plt.figure(figsize=(12, min(12, 0.5 * len(matrix))))
    sns.heatmap(matrix, cmap="Greys", cbar=False)
    plt.title("Pathway–Species Occurrence Matrix")
    plt.xlabel("Species")
    plt.ylabel("Pathway")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved: {save_path}")
    plt.close()

def plot_species_net_percentage_vs_temperature(cell_data_dict, temp_file_path, species_name, log_x=True, save_path=None):
    """
    Plots and optionally saves net production-consumption % vs. temperature.

    Parameters
    ----------
    cell_data_dict : dict
    temp_file_path : str
    species_name : str
    log_x : bool
    save_path : str or None

    Returns
    -------
    None
    """
    temp_df = pd.read_csv(temp_file_path, delim_whitespace=True)
    temp_map = dict(zip(temp_df['Radius_index'], temp_df['Avg_Temperature']))
    net_pct_data = []
    for cell, df in cell_data_dict.items():
        if cell not in temp_map:
            continue
        T = temp_map[cell]
        df_species = df[df['species'] == species_name]
        if not df_species.empty:
            net_pct_total = df_species['net_pct'].sum()
            net_pct_data.append((T, net_pct_total))
    df_plot = pd.DataFrame(net_pct_data, columns=['Temperature', 'NetPct'])
    df_plot = df_plot.sort_values('Temperature')
    plt.figure(figsize=(10, 6))
    plt.plot(df_plot['Temperature'], df_plot['NetPct'], marker='o')
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    if log_x:
        plt.xscale('log')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Net Production − Consumption (%)")
    plt.title(f"Net Production/Consumption % of {species_name} vs Temperature")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved: {save_path}")
    plt.close()

def main(
    data_dir,
    temp_file_path,
    species_list,
    outdir="./",
    top_n=5
):
    """
    Generate and save all pathway plots for a list of species.

    Parameters
    ----------
    data_dir : str
    temp_file_path : str
    species_list : list of str
    outdir : str
    top_n : int

    Returns
    -------
    None

    Example
    -------
    >>> # main("data/", "qt_conditions.txt", ["H", "H2"])
    """
    os.makedirs(outdir, exist_ok=True)
    cell_dict = load_all_cells_as_dict(data_dir)
    for species in species_list:
        print(f"Plotting for species: {species}")
        plot_top_species_pathways_by_rate(
            cell_dict, species, top_n=top_n,
            save_path=os.path.join(outdir, f"{species}_top_pathways.png"))
        plot_species_net_rate(
            cell_dict, species,
            save_path=os.path.join(outdir, f"{species}_net_rate.png"))
        plot_pathway_heatmap(
            cell_dict, species,
            save_path=os.path.join(outdir, f"{species}_heatmap.png"))
        plot_stacked_production_fraction(
            cell_dict, species, top_n=top_n,
            save_path=os.path.join(outdir, f"{species}_stacked_fraction.png"))
        plot_species_net_percentage_vs_temperature(
            cell_dict, temp_file_path, species,
            save_path=os.path.join(outdir, f"{species}_net_pct_vs_temp.png")
        )
    matrix = compute_pathway_species_matrix(cell_dict)
    plot_pathway_species_heatmap(
        matrix, save_path=os.path.join(outdir, "pathway_species_matrix.png")
    )
    print("All plots generated.")

# If running directly:
if __name__ == "__main__":
    # Example usage:
    main(
        data_dir="./Outputs/",
        temp_file_path="../qt_conditions.txt",
        species_list=["H", "H2", "O", "OH"],
        outdir="./Plots/",
        top_n=5
    )