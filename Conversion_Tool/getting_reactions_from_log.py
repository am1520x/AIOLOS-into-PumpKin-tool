import re
import pandas as pd

def extract_reaction_blocks(content):
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
    # Match Reaction #XXX: <reaction> [stop at ...... or Available flux or newlines]
    return re.findall(
        r"(Reaction\s+#\d+:\s+.+?->.+?)(?=\s+\.{2,}|Available flux|\ndG|\n|$)", 
        block,
        re.DOTALL
    )

def extract_all_reactions(file_path):
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
