import re
import os
import traceback
import matplotlib.pyplot as plt

def parse_multi_line_pathway_table(file_content):
    """
    Parses both single-line and multi-line pathway table formats.

    Cleans whitespace, removes leading coefficients like '1 * ', and reconstructs
    multi-line pathways with proper formatting.

    Parameters
    ----------
    file_content : str
        Content of the file as a string.

    Returns
    -------
    dict
        Dictionary mapping pathway strings to their rates.

    Example
    -------
    >>> test_table = '''
    ... | Some header | Rate of pathway |
    ... +------------+-----------------+
    ... | H2 + O2 (1) | 3.5e-05        |
    ... +------------+-----------------+
    ... '''
    >>> d = parse_multi_line_pathway_table(test_table)
    Parsed 1 pathways.
    >>> list(d.items())[0][0].startswith('H2 + O2')
    True
    >>> abs(list(d.items())[0][1] - 3.5e-5) < 1e-8
    True
    """
    lines = file_content.splitlines()
    pathway_dict = {}
    in_table_section = False
    past_header_separator = False

    header_content_pattern = re.compile(r"^\|.*Rate of pathway.*\|")
    separator_pattern = re.compile(r"^\+[-+]+\+$")
    end_marker_pattern = re.compile(r"^#+\s*$")
    rate_line_pattern = re.compile(r"\|\s*(.*?)\(\s*\d+\)\s*\|\s*([-\deE+\.]+)\s*\|")

    buffer = []

    def clean_step(step):
        step = re.sub(r"^\s*\d+\s*\*\s*", "", step)  # Remove '1 *' or similar
        step = step.strip()
        match = re.search(r"\((.*)\)", step)
        if match:
            return f"({match.group(1).strip()})"
        return step

    for i, line in enumerate(lines):
        if header_content_pattern.search(line):
            in_table_section = True
            continue

        if in_table_section and not past_header_separator:
            if separator_pattern.match(line):
                past_header_separator = True
            continue

        if past_header_separator:
            if end_marker_pattern.match(line):
                break

            if separator_pattern.match(line):
                continue

            rate_match = rate_line_pattern.match(line)
            if rate_match:
                buffer.append(rate_match.group(1).strip())
                steps = [clean_step(p) for p in buffer if p]
                pathway = " -> ".join(steps)
                try:
                    rate = float(rate_match.group(2).strip())
                    pathway_dict[pathway] = rate
                except ValueError:
                    print(f"Warning: Failed to parse rate: {rate_match.group(2)} on line {i}")
                buffer = []
            else:
                if "|" in line:
                    part = line.strip().strip("|").split("|")[0].strip()
                    if part:
                        buffer.append(part)

    print(f"Parsed {len(pathway_dict)} pathways.")
    return pathway_dict

def parse_single_line_pathway_table(file_content):
    """
    Extracts pathway strings and their rates from a single-line table format.

    Parameters
    ----------
    file_content : str
        Content of the file as a string.

    Returns
    -------
    dict
        Dictionary mapping pathway strings to their rates.

    Example
    -------
    >>> test_table = '''
    ... | Pathway (Index)      | Rate of pathway |
    ... +---------------------+-----------------+
    ... | H2 + O2 (1)         | 3.5e-05         |
    ... | CO + O2 (2)         | 4.1e-07         |
    ... ##############################
    ... '''
    >>> d = parse_single_line_pathway_table(test_table)
    --- Starting parse_single_line_pathway_table ---
    --- Found header content line at line 1 ---
    --- Passed header separator at line 2 ---
    --- Found data line at line 3 ---
    --- Found data line at line 4 ---
    --- Found end marker at line 5 ---
    --- Finished parse_single_line_pathway_table. Found 2 pathways. ---
    >>> abs(d["H2 + O2"] - 3.5e-5) < 1e-8
    True
    >>> abs(d["CO + O2"] - 4.1e-7) < 1e-8
    True
    """
    lines = file_content.splitlines()
    pathway_dict = {}
    in_table_section = False
    past_header_separator = False

    header_content_pattern = re.compile(r"^\|.*Rate of pathway.*\|")
    data_line_pattern = re.compile(r"\|\s*(.+?\s*\(\s*\d+\)\s*)\s*\|\s*([-\deE\+\-\.]+)\s*\|")
    separator_pattern = re.compile(r"^\+[-+]+\+\s*$")
    end_marker_pattern = re.compile(r"^#+\s*$")

    print("--- Starting parse_single_line_pathway_table ---")

    for i, line in enumerate(lines):
        if header_content_pattern.search(line):
            print(f"--- Found header content line at line {i} ---")
            in_table_section = True
            continue

        if in_table_section and not past_header_separator:
            if separator_pattern.match(line):
                print(f"--- Passed header separator at line {i} ---")
                past_header_separator = True
                continue

        if past_header_separator:
            if end_marker_pattern.match(line):
                print(f"--- Found end marker at line {i} ---")
                break
            if separator_pattern.match(line):
                continue

            data_match = data_line_pattern.match(line)
            if data_match:
                print(f"--- Found data line at line {i} ---")
                raw_pathway_part = data_match.group(1).strip()
                rate_str = data_match.group(2).strip()
                pathway_str = re.sub(r'\s*\(\s*\d+\)\s*$', '', raw_pathway_part).strip()
                try:
                    rate = float(rate_str)
                    pathway_dict[pathway_str] = rate
                except ValueError:
                    print(f"Warning: Could not parse rate '{rate_str}' in line: {line.strip()}")
                    continue

    print(f"--- Finished parse_single_line_pathway_table. Found {len(pathway_dict)} pathways. ---")
    return pathway_dict

def main(
    data_dir="./",
    num_cells=233,
    use_multi_line_parser=True,
    top_n_pathways=25,
    save_plot_path=None
):
    """
    Main driver: parses all cell files, aggregates and plots pathway data.

    Parameters
    ----------
    data_dir : str
        Directory containing input files.
    num_cells : int
        Number of cell files to process.
    use_multi_line_parser : bool
        If True, uses the multi-line parser, else uses single-line.
    top_n_pathways : int
        Number of top pathways to plot.
    save_plot_path : str or None
        If not None, saves the plot to this path.

    Returns
    -------
    None

    Example
    -------
    >>> import tempfile, os
    >>> tmpdir = tempfile.mkdtemp()
    >>> fname = os.path.join(tmpdir, "pumpkin_output_cell_000.txt")
    >>> with open(fname, "w") as f:
    ...     _ = f.write('''
    ... | Pathway (Index)      | Rate of pathway |
    ... +---------------------+-----------------+
    ... | H2 + O2 (1)         | 3.5e-05         |
    ... ##############################
    ... ''')
    >>> main(tmpdir, num_cells=1, use_multi_line_parser=False, top_n_pathways=1)
    Parsing files...
    --- Starting parse_single_line_pathway_table ---
    --- Found header content line at line 1 ---
    --- Passed header separator at line 2 ---
    --- Found data line at line 3 ---
    --- Found end marker at line 4 ---
    --- Finished parse_single_line_pathway_table. Found 1 pathways. ---
    Finished parsing 1 files.
    Found a total of 1 unique pathways.
    Pathways found (sample): 
    - H2 + O2
    Plotting top pathways (log-log)...
    <BLANKLINE>
    Top 1 Pathways ranked by total rate:
    Pathway: H2 + O2, Total Rate: 3.50e-05, Nonzero cells: 1
    <BLANKLINE>
    Summary for first 10 pathways:
    Pathway: H2 + O2, Max rate: 3.50e-05, Nonzero points: 1
    dict_keys(['H2 + O2'])
    """
    print("Parsing files...")
    all_pathways = set()
    parsed_rates_per_cell = []

    for i in range(num_cells):
        filename = os.path.join(data_dir, f"pumpkin_output_cell_{i:03d}.txt")
        if not os.path.exists(filename):
            print(f"Warning: {filename} not found, skipping.")
            parsed_rates_per_cell.append({})
            continue

        try:
            with open(filename, "r") as f:
                content = f.read()

            if use_multi_line_parser:
                rates = parse_multi_line_pathway_table(content)
            else:
                rates = parse_single_line_pathway_table(content)
            parsed_rates_per_cell.append(rates)
            all_pathways.update(rates.keys())
        except Exception as e:
            print(f"Error processing file {filename}: {e}")
            traceback.print_exc()
            parsed_rates_per_cell.append({})

    num_cells = len(parsed_rates_per_cell)
    print(f"Finished parsing {num_cells} files.")
    print(f"Found a total of {len(all_pathways)} unique pathways.")

    pathway_data = {p: [0.0] * num_cells for p in all_pathways}

    for cell_i, rates_dict in enumerate(parsed_rates_per_cell):
        for pathway_str, rate in rates_dict.items():
            if pathway_str in pathway_data:
                pathway_data[pathway_str][cell_i] = rate

    any_none = any(
        any(val is None for val in rates_list)
        for rates_list in pathway_data.values()
    )
    # print("Any None left after filling?", any_none)

    print("Pathways found (sample): ")
    for i, p in enumerate(list(all_pathways)[:10]):
        print(f"- {p}")
    if len(all_pathways) > 10:
        print("...")

    print("Plotting top pathways (log-log)...")

    ranked_pathways = sorted(
        [(p, sum(rates)) for p, rates in pathway_data.items()],
        key=lambda x: x[1],
        reverse=True
    )
    top_pathways_to_plot = ranked_pathways[:top_n_pathways]

    print(f"\nTop {len(top_pathways_to_plot)} Pathways ranked by total rate:")
    for p_str, total_rate in top_pathways_to_plot:
        nonzero_count = sum(1 for r in pathway_data.get(p_str, []) if r > 0)
        print(f"Pathway: {p_str[:80]}{'...' if len(p_str) > 80 else ''}, Total Rate: {total_rate:.2e}, Nonzero cells: {nonzero_count}")

    plt.figure(figsize=(16, 10))
    for p_str, _ in top_pathways_to_plot:
        rates = pathway_data.get(p_str, [])
        x_vals, y_vals = zip(*[(i, rate) for i, rate in enumerate(rates) if isinstance(rate, (int, float)) and rate > 0]) if rates else ([], [])
        if x_vals:
            plt.plot(x_vals, y_vals, marker='.', linestyle='-', linewidth=0.8, markersize=4, label=p_str)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Cell Number (log)")
    plt.ylabel("Pathway Rate (log)")
    plt.title(f"Top {len(top_pathways_to_plot)} Pathway Rates vs Cell Number (Log-Log)")
    plt.legend(title='Pathway', loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=7)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    if save_plot_path:
        plt.savefig(save_plot_path, dpi=300)
        print(f"Plot saved to {save_plot_path}")
    plt.show()

    # print("\nPlotting all pathways (non-ranked, log-log)...")
    # plt.figure(figsize=(16, 10))
    # for p_str, rates in pathway_data.items():
    #     x_vals, y_vals = zip(*[(i, rate) for i, rate in enumerate(rates) if isinstance(rate, (int, float)) and rate > 0]) if rates else ([], [])
    #     if x_vals:
    #         plt.plot(x_vals, y_vals, 'x', markersize=4, label=p_str)
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel("Cell Number (log)")
    # plt.ylabel("Pathway Rate (log)")
    # plt.title("All Pathway Rates vs Cell Number (Log-Log)")
    # plt.show()

    print("\nSummary for first 10 pathways:")
    for p, rates in list(pathway_data.items())[:10]:
        nonzero_count = sum(r > 0 for r in rates)
        valid_rates = [r for r in rates if isinstance(r, (int, float))]
        max_rate = max(valid_rates) if valid_rates else 0.0
        print(f"Pathway: {p[:80]}{'...' if len(p) > 80 else ''}, Max rate: {max_rate:.2e}, Nonzero points: {nonzero_count}")
    print(pathway_data.keys())

if __name__ == "__main__":
    main()
