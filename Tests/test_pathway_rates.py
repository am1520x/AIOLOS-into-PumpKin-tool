import pytest, sys, os

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from Outputs.pathway_rates_plots import (
    parse_multi_line_pathway_table,
    parse_single_line_pathway_table,
)


def test_parse_multi_line_pathway_table_basic():
    content = """
|   Pathway | Rate of pathway |
+----------------+-----------------+
| H2 + O2 (1)    | 1.2e-04         |
+----------------+-----------------+
"""
    d = parse_multi_line_pathway_table(content)
    assert any("H2 + O2" in k for k in d)
    assert abs(list(d.values())[0] - 1.2e-4) < 1e-8


def test_parse_single_line_pathway_table_basic():
    content = """
| Pathway (Index) | Rate of pathway |
+-----------------+-----------------+
| H2 + O2 (1)     | 1.2e-04         |
| CO + O2 (2)     | 3.1e-05         |
#########################
"""
    d = parse_single_line_pathway_table(content)
    assert "H2 + O2" in d
    assert abs(d["H2 + O2"] - 1.2e-4) < 1e-8
    assert "CO + O2" in d
    assert abs(d["CO + O2"] - 3.1e-5) < 1e-8


def test_parse_multi_line_pathway_table_malformed_line():
    content = """
|   Reaction Step | Rate of pathway |
+----------------+-----------------+
| Bad line       | not_a_number    |
+----------------+-----------------+
"""
    d = parse_multi_line_pathway_table(content)
    assert len(d) == 0  # Should skip line


def test_parse_single_line_pathway_table_with_index_spaces():
    content = """
| Pathway (Index) | Rate of pathway |
+-----------------+-----------------+
| H2 + O2 ( 1 )   | 9.9e-01         |
############################
"""
    d = parse_single_line_pathway_table(content)
    assert "H2 + O2" in d
    assert abs(d["H2 + O2"] - 0.99) < 1e-8


def test_parse_multi_line_pathway_table_multistep():
    content = """
| Step | Rate of pathway |
+------+-----------------+
| 1 * H2 + O2 (1) | 1.0e-02 |
+------+-----------------+
| CO2 + H2O (2)   | 5.0e-03 |
+------+-----------------+
"""
    d = parse_multi_line_pathway_table(content)
    assert any("H2 + O2" in k for k in d)
    assert abs(list(d.values())[0] - 1.0e-2) < 1e-8


def test_parse_multi_line_pathway_table_complex_block():
    content = """
+---------------------------------------------------+---------------------------+
|       Pathway                                     |       Rate of pathway     |
+---------------------------------------------------+---------------------------+
| 1 * (H2+OH=>H2O+H                        ) (     0) | 5.248e+08                 |
+---------------------------------------------------+---------------------------+
| 1 * (H+H2O=>OH+H2                        ) (     1) | 5.2465e+08                |
+---------------------------------------------------+---------------------------+
| 1 * (H+O2H=>O2+H2                        ) (    2) | 3027                      |
+---------------------------------------------------+---------------------------+
| 1 * (S1+H2O=>H2O^++H                            ) |                           |
| 1 * (H+O^+=>O+S1                         ) (    3) | 3025.98                   |
+---------------------------------------------------+---------------------------+
| 1 * (OH=>O+H                             ) (    4) | 2945.5                    |
+---------------------------------------------------+---------------------------+
| 1 * (S1+H2O=>H2O^++H                            ) |                           |
| 1 * (H+H2^+=>H2+S1                       ) (    5) | 2797.48                   |
+---------------------------------------------------+---------------------------+
#################################################################################
"""
    # Call your parser
    result = parse_multi_line_pathway_table(content)
    # Check that only the pathways with rates are present
    assert len(result) == 6  # Only 6 have rates, 2 are blank and should be skipped

    # The actual keys will look like e.g. "(H2+OH=>H2O+H)", "(H+H2O=>OH+H2)", etc.
    # The clean_step logic wraps them in parentheses as in your code.
    print(result)
    expected = {
        "(H2+OH=>H2O+H)": 5.248e8,
        "(H+H2O=>OH+H2)": 5.2465e8,
        "(H+O2H=>O2+H2)": 3027.0,
        "(S1+H2O=>H2O^++H) -> (H+O^+=>O+S1)": 3025.98,
        "(OH=>O+H)": 2945.5,
        "(S1+H2O=>H2O^++H) -> (H+H2^+=>H2+S1)": 2797.48,
    }

    for k, v in expected.items():
        assert k in result, f"Missing pathway: {k}"
        assert abs(result[k] - v) < 1e-3, f"Value mismatch for {k}: {result[k]} vs {v}"

    # Check that no blank/empty rates have been added
    for k in result:
        assert result[k] != "" and result[k] is not None
