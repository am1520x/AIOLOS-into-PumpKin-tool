"""
This script extracts chemical reaction information from the AIOLOS log file.

It identifies and parses both thermochemical and photochemical reactions, returning them
as a structured pandas DataFrame for further analysis or processing.

Functions:
- extract_reaction_blocks(content): Separates the thermochemical and photochemical blocks.
- extract_reactions(block): Extracts individual reaction strings from a block.
- extract_all_reactions(file_path): Processes a log file and returns a DataFrame of reactions.

Example usage:
    python getting_reactions_from_log.py  # Will print the first and last few parsed reactions from the given file.
"""
import re
import pandas as pd

def extract_reaction_blocks(content):
    """
    Extract thermochemical and photochemical reaction sections from the log file content.

    Args:
        content (str): The full content of a chemistry log file.

    Returns:
        tuple: (thermo_block (str), photo_block (str))

    Example:
        >>> thermo, photo = extract_reaction_blocks(\"""
        In init chemistry, reporting thermochemical reactions:
        Reaction #1: A + B -> C ...
        Reporting photoreactions:
        Reaction #2: X + photon -> Y ...
        \""")
        >>> "Reaction #1" in thermo
        True
        >>> "Reaction #2" in photo
        True
    """
    thermo_block = re.search(
        r"In init chemistry, reporting thermochemical reactions:(.*?)(?:Reporting photoreactions:|$)",
        content,
        re.DOTALL
    )
    photo_block = re.search(
        r"Reporting photoreactions:(.*?)(?:Photoreaction opacities table|$)",
        content,
        re.DOTALL
    )
    return (thermo_block.group(1) if thermo_block else ""), (photo_block.group(1) if photo_block else "")

def extract_reactions(block):
    """
    Extract individual reactions from a given text block.

    Args:
        block (str): A string block containing reaction definitions.

    Returns:
        list: A list of reaction strings.

    Example:
        >>> extract_reactions("Reaction #1: A + B -> C ...\\nAvailable flux\\nReaction #2: X -> Y ...")
        ['Reaction #1: A + B -> C', 'Reaction #2: X -> Y']
    """
    return re.findall(
        r"(Reaction\s+#\d+:\s+.+?->.+?)(?=\s+\.{2,}|Available flux|\ndG|\n|$)", 
        block,
        re.DOTALL
    )

def extract_all_reactions(file_path):
    """
    Extract all thermochemical and photochemical reactions from the log file.

    Args:
        file_path (str): Path to the log file.

    Returns:
        pd.DataFrame: A DataFrame with columns ['number', 'reaction'].

    Example:
        >>> import tempfile
        >>> content = '''In init chemistry, reporting thermochemical reactions:
        ... Reaction #1: H2 + O -> OH + H ...
        ... Reporting photoreactions:
        ... Reaction #2: CO2 + photon -> CO + O + (hv 10.5 eV ) ...'''
        >>> with tempfile.NamedTemporaryFile(mode="w+", delete=False) as f:
        ...     _ = f.write(content)
        ...     f.flush()
        ...     df = extract_all_reactions(f.name)
        >>> df.shape[0]
        2
        >>> df.iloc[0]['reaction']
        'H2 + O -> OH + H'
    """
    with open(file_path, 'r') as f:
        content = f.read()

    thermo_block, photo_block = extract_reaction_blocks(content)

    thermo_reactions = extract_reactions(thermo_block)
    photo_reactions = extract_reactions(photo_block)

    reactions = {
        "thermo": thermo_reactions,
        "photo": photo_reactions
    }

    data = []
    
    for reac in reactions["thermo"]:
        match = re.match(r'Reaction\s+#(\d+):\s+(.*)', reac)
        if match:
            number = int(match.group(1))
            reaction = match.group(2)
            data.append({'number': number, 'reaction': reaction})
    for reac in reactions["photo"]:
        match = re.match(r'Reaction\s+#(\d+):\s+(.*)', reac)
        if match:
            number = int(match.group(1))
            reaction = match.group(2).split(' + (')[0] + match.group(2).split(' + (')[1].split('eV ')[1].strip()
            data.append({'number': number, 'reaction': reaction})
        
    table = pd.DataFrame(data)
    return table


# Example usage
if __name__ == "__main__":
    file_path = r"2025_waterworlds_cheminfo_h2-2\log_dynamic_h2-2-h2odom_frac-long153cut-step8-sol2eps-4-dist003-fullspectrum2-res35nocut-cheminfo.log"  # Replace with your file path
    reactions = extract_all_reactions(file_path)
    
    print(reactions.head())
    print(reactions.tail())
