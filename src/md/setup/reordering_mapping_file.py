order_list = [
    'P',
    'O',
    'O1',
    'O2',
    'P1',
    'O3',
    'O4',
    'O5',
    'P2',
    'O6',
    'O7',
    'O8',
    'O9',
    'C',
    'C1',
    'O10',
    'C2',
    'O11',
    'C3',
    'O12',
    'C4',
    'N',
    'C5',
    'N1',
    'C6',
    'C7',
    'N2',
    'N3',
    'C8',
    'N4',
    'C9',
    'H',
    'H1',
    'H2',
    'H3',
    'H4',
    'H5',
    'H6',
    'H7',
    'H8',
    'H9',
    'H10',
    'H11'
]

def reorder_mapping_file(input_file, output_file=None):
    """
    Reorder lines in a mapping file based on order_list.
    Assumes the atom name to match is in a specific column.
    
    Parameters:
    -----------
    input_file : str
        Path to the input mapping file
    output_file : str, optional
        Path to save the reordered file (if None, overwrites input)
    """
    if output_file is None:
        output_file = input_file
    
    # Read all lines from the file
    lines_dict = {}
    header_lines = []
    unmatched_lines = []
    
    try:
        with open(input_file, 'r') as f:
            for line in f:
                # Keep header/comment lines
                if line.startswith('#'):
                    header_lines.append(line)
                    continue
                
                stripped = line.strip()
                if not stripped:
                    continue
                
                # Split line and get atom name (adjust column index as needed)
                parts = stripped.split()
                if len(parts) >= 2:
                    # Assuming format: "RESNAME ATOM1 ATOM2" or similar
                    # Adjust the index based on your file format
                    atom_name = parts[2].strip()  # Change index if needed
                    
                    if atom_name in order_list:
                        # Store line with its order index
                        if atom_name not in lines_dict:
                            lines_dict[atom_name] = []
                        lines_dict[atom_name].append(line)
                    else:
                        unmatched_lines.append(line)
        
        # Write reordered lines
        with open(output_file, 'w') as f:
            # Write headers first
            for header in header_lines:
                f.write(header)
            
            # Write lines in order_list order
            for atom in order_list:
                if atom in lines_dict:
                    for line in lines_dict[atom]:
                        f.write(line)
            
            # Write unmatched lines at the end
            if unmatched_lines:
                f.write("\n# Unmatched atoms:\n")
                for line in unmatched_lines:
                    f.write(line)
        
        print(f"✓ Reordered {len(lines_dict)} atom types")
        print(f"✓ Output saved to: {output_file}")
        if unmatched_lines:
            print(f"⚠ {len(unmatched_lines)} unmatched lines appended at end")
        
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found!")
    except Exception as e:
        print(f"Error: {e}")


# Example usage
if __name__ == "__main__":
    # Reorder the mapping file
    from pathlib import Path
    _here = Path(__file__).resolve().parent
    input_mapping_file = str(_here / "ATP_AF3_to_amber_mapping.txt")
    output_mapping_file = str(_here / "atom_mapping_reordered.txt")

    reorder_mapping_file(input_mapping_file, output_mapping_file)
    
    # Or overwrite the original file:
    # reorder_mapping_file("atom_mapping.txt")



