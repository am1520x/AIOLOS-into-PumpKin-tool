import tempfile, os, sys
import numpy as np
import pandas as pd

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from Outputs.species_specific_plots import parse_all_species_pathway_tables_multiline, load_all_cells_as_dict, compute_pathway_species_matrix 

def test_parse_all_species_pathway_tables_multiline_basic():
    text = '''
| Pathway | Production of H | Consumption of H |
+---+---+---+
| 1 * (A=>B) 1202 | 22% | 5% |
# End block
'''
    df = parse_all_species_pathway_tables_multiline(text)
    assert "species" in df.columns
    assert "H" in df["species"].values
    assert abs(df.iloc[0]["production_pct"] - 22.0) < 1e-8

def test_load_all_cells_as_dict(tmp_path):
    fpath = tmp_path / "pumpkin_output_cell_001.txt"
    text = (
        "| Pathway | Production of X | Consumption of X |\n"
        "+-+-+-+\n"
        "| 1 * (A=>B) (0) | 90% | 70% |\n#"
    )
    fpath.write_text(text)
    d = load_all_cells_as_dict(tmp_path)
    assert 1 in d
    assert d[1]["species"].iloc[0] == "X"

def test_basic_pathway_table_single_row():
    block = '''
|       Pathway                                     | Production of H           | Consumption of H          |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H+H2O=>OH+H2)             5.2465e+08 |                       0 % |                      91 % |
+---------------------------------------------------+---------------------------+---------------------------+
# End block
'''
    df = parse_all_species_pathway_tables_multiline(block)
    assert len(df) == 1
    assert df.iloc[0]['species'] == "H"
    assert abs(df.iloc[0]['rate'] - 5.2465e8) < 1e-3
    assert df.iloc[0]['consumption_pct'] == 91.0
    assert df.iloc[0]['production_pct'] == 0.0

def test_empty_pathway_table():
    block = '''
|       Pathway                                     | Production of E           | Consumption of E          |
+---------------------------------------------------+---------------------------+---------------------------+
#############################################################################################################
'''
    df = parse_all_species_pathway_tables_multiline(block)
    # Should be empty DataFrame or contain only NaNs
    assert df.empty or all(np.isnan(df['rate']))

def test_complex_multistep_pathway_table():
    block = '''
|       Pathway                                     | Production of H           | Consumption of H          |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H2+O=>OH+H)            6.52246e+08 |                      62 % |                       0 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H+H2=>H)             1.5601e+08 |                      30 % |                       0 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (ANY_NEUTRAL+H2=>ANY_NEUTRAL+H)             2.3101e+07 |                     4.4 % |                       0 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H2+O=>OH+H                         )         |                           |                           |
| 1 * (H+OH=>O+H)            1.59806e+07 |                     3.1 % |                       0 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H+O2=>OH+O                         )         |                           |                           |
| 1 * (H+OH=>O+H2)                 499051 |                       0 % |                      38 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H+H2O=>OH+H2                       )         |                           |                           |
| 1 * (H+OH=>O+H2)                 331737 |                       0 % |                      25 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H+OH=>O+H2                         )         |                           |                           |
| 1 * (H+O=>OH)                 246812 |                       0 % |                      19 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H+O2=>OH+O                         )         |                           |                           |
| 1 * (H+OH=>O+H2)                 165588 |                       0 % |                      13 % |
| 1 * (O=>O2                              )         |                           |                           |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (ANY_NEUTRAL+H=>ANY_NEUTRAL+H2)                  30365 |                       0 % |                     2.3 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (O+H2O=>OH                          )         |                           |                           |
| 2 * (H+OH=>O+H2)                  21531 |                       0 % |                     1.6 % |
+---------------------------------------------------+---------------------------+---------------------------+
#############################################################################################################
'''
    df = parse_all_species_pathway_tables_multiline(block)
    # It should parse all lines with a rate (multi-line support)
    assert (df['rate'] > 0).sum() >= 6  # Should parse several rates > 0
    assert "H" in df['species'].unique()

def test_pathway_table_missing_columns():
    block = '''
| Pathway | Production of H | Consumption of H |
+----+----+----+
| nonsense line |
+----+----+----+
# End block
'''
    df = parse_all_species_pathway_tables_multiline(block)
    # Should skip lines with missing/invalid columns
    assert df.empty or all(df['pathway'].isna())

def test_pathway_with_percent_symbols_and_weird_spacing():
    block = '''
| Pathway | Production of O | Consumption of O |
+----+----+----+
| 1 * (A=>B)     1234 |    33%    |   12 % |
+----+----+----+
# End block
'''
    df = parse_all_species_pathway_tables_multiline(block)
    assert df.iloc[0]['species'] == "O"
    assert df.iloc[0]['production_pct'] == 33.0
    assert df.iloc[0]['consumption_pct'] == 12.0
    assert df.iloc[0]['rate'] == 1234

# This one uses the full realistic block you posted:
def test_full_realistic_block():
    block = '''
|       Pathway                                     | Production of H           | Consumption of H          |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H2+OH=>H2O+H)              5.248e+08 |                      87 % |                       0 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H+H2O=>OH+H2)             5.2465e+08 |                       0 % |                      91 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (S1+H2O=>H2O^++H)            5.13142e+07 |                     8.5 % |                       0 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (ANY_NEUTRAL+H=>ANY_NEUTRAL+H2)               2.71e+07 |                       0 % |                     9.4 % |
+---------------------------------------------------+---------------------------+---------------------------+
| 1 * (H2^++H2O=>H3O^++H)            2.52915e+07 |                     4.2 % |                       0 % |
+---------------------------------------------------+---------------------------+---------------------------+
#############################################################################################################
'''
    df = parse_all_species_pathway_tables_multiline(block)
    assert set(df['species']) == {"H"}
    assert "(H2+OH=>H2O+H)" in set(df['pathway'])
    assert np.isclose(df[df['pathway'] == '(H2+OH=>H2O+H)']['rate'].iloc[0], 5.248e8)
    assert df[df['pathway'] == '(H+H2O=>OH+H2)']['consumption_pct'].iloc[0] == 91.0

def test_load_all_cells_as_dict(tmp_path):
    # Write two cell files
    cell_content = (
        "| Pathway | Production of H | Consumption of H |\n"
        "+----+----+----+\n"
        "| 1 * (A=>B) 1234 | 10% | 3% |\n"
        "#"
    )
    for i in [2, 4]:
        fname = tmp_path / f"pumpkin_output_cell_{i:03d}.txt"
        fname.write_text(cell_content)
    d = load_all_cells_as_dict(str(tmp_path))
    assert set(d.keys()) == {2, 4}
    assert d[2].iloc[0]['species'] == "H"
    assert d[2].iloc[0]['rate'] == 1234

def test_load_all_cells_empty_dir(tmp_path):
    d = load_all_cells_as_dict(str(tmp_path))
    assert d == {}

def test_load_all_cells_missing_proper_filenames(tmp_path):
    # Put in a non-matching file
    (tmp_path / "some_other_file.txt").write_text("irrelevant")
    d = load_all_cells_as_dict(str(tmp_path))
    assert d == {}

def test_compute_pathway_species_matrix_basic():
    df = pd.DataFrame([
        {"species": "H", "pathway": "(A=>B)", "rate": 10},
        {"species": "O", "pathway": "(A=>B)", "rate": 10},
        {"species": "O", "pathway": "(B=>C)", "rate": 5},
        {"species": "H", "pathway": "(C=>D)", "rate": 3},
    ])
    data = {0: df}
    mtx = compute_pathway_species_matrix(data)
    assert set(mtx.index) == {"(A=>B)", "(B=>C)", "(C=>D)"}
    assert set(mtx.columns) == {"H", "O"}
    assert mtx.loc["(A=>B)", "H"] == 1
    assert mtx.loc["(A=>B)", "O"] == 1
    assert mtx.loc["(B=>C)", "O"] == 1
    assert mtx.loc["(C=>D)", "H"] == 1
    assert mtx.loc["(B=>C)", "H"] == 0

def test_compute_pathway_species_matrix_empty():
    mtx = compute_pathway_species_matrix({})
    assert isinstance(mtx, pd.DataFrame)
    assert mtx.empty