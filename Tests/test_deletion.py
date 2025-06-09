import tempfile
import numpy as np
import os
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from Outputs.deleted_pathways_plots import parse_deleted_percent_block, parse_file, deleted_contributions, species_list, flagged_entries, N_CELLS

def test_parse_deleted_percent_block_basic():
    # Setup
    species_list.clear()
    flagged_entries.clear()
    lines = [
        "|   H2O   |   15.0%   |   3.2%   |",
        "|   CO2   |   7.1%    |   12.0%  |"
    ]
    parse_deleted_percent_block(lines, 0)

    assert 'H2O' in species_list
    assert deleted_contributions['CO2']['production'][0] == 7.1
    assert flagged_entries[-1][1] == 'H2O'

def test_parse_file_integration(tmp_path):
    # Prepare a fake file
    content = """
Procent of species from deleted pathways
+----------+----------+----------+
|   H2O    |   11.0%  |   2.1%   |
|   CO2    |   8.1%   |   18.0%  |
+----------+----------+----------+
"""
    fpath = tmp_path / "pumpkin_output_cell_000.txt"
    fpath.write_text(content)
    cell_idx = 5
    parse_file(str(fpath), cell_idx)
    assert deleted_contributions['H2O']['production'][cell_idx] == 11.0
    assert deleted_contributions['CO2']['consumption'][cell_idx] == 18.0

def test_flagged_entries_threshold():
    species_list.clear()
    flagged_entries.clear()
    lines = [
        "|   ABC   |   2.0%    |   20.1%  |",  # should be flagged due to cons > 10
        "|   XYZ   |   15.2%   |   0.2%   |",  # should be flagged due to prod > 10
        "|   LMN   |   9.0%    |   9.0%   |",  # should NOT be flagged
    ]
    parse_deleted_percent_block(lines, 1)
    flagged_species = {entry[1] for entry in flagged_entries}
    assert "ABC" in flagged_species
    assert "XYZ" in flagged_species
    assert "LMN" not in flagged_species

def test_index_range_protection(capsys):
    species_list.clear()
    flagged_entries.clear()
    lines = ["|   OOB   |   5.0%    |   9.0%   |"]
    parse_deleted_percent_block(lines, N_CELLS + 10)
    captured = capsys.readouterr()
    assert "outside the expected range" in captured.out

def test_nan_on_invalid_data():
    species_list.clear()
    flagged_entries.clear()
    lines = ["|  BAD  |  nan%  |  nan%  |"]
    parse_deleted_percent_block(lines, 2)
    val = deleted_contributions["BAD"]["production"][2]
    assert np.isnan(val)
