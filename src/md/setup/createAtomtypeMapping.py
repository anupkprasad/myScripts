#!/usr/bin/env python3
"""
Create a mapping between atom names in two PDB files based on coordinates (having same coordinates).
"""

import sys
import math


def read_pdb(filename):
    """Read PDB and return list of (atom_name, resname, x, y, z)."""
    atoms = []
    with open(filename) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_name, resname, x, y, z))
    return atoms

def distance(a, b):
    """Calculate distance between two atoms based on coordinates."""
    return math.sqrt((a[2]-b[2])**2 + (a[3]-b[3])**2 + (a[4]-b[4])**2)

def build_mapping(pdb1, pdb2, tol=0.5):
    """Return mapping between atom names of pdb1 and pdb2 based on coordinates."""
    mapping = []
    unmatched = []
    for a1 in pdb1:
        found = None
        for a2 in pdb2:
            if distance(a1, a2) < tol:
                found = a2
                break
        if found:
            mapping.append((a1[1], a1[0], found[0]))  # (resname, name1, name2)
        else:
            unmatched.append((a1[1], a1[0]))  # (resname, name1)
    return mapping, unmatched

def main(pdb1_file, pdb2_file, outfile="atom_mapping.txt"):
    pdb1 = read_pdb(pdb1_file)
    pdb2 = read_pdb(pdb2_file)

    mapping, unmatched = build_mapping(pdb1, pdb2)

    with open(outfile, "w") as f:
        f.write("# Residue  AtomName_file1  AtomName_file2\n")
        for resname, name1, name2 in mapping:
            f.write(f"{resname:>4}  {name1:>6}  {name2:>6}\n")
        if unmatched:
            f.write("\n# Unmatched atoms from file1:\n")
            for resname, name1 in unmatched:
                f.write(f"{resname:>4}  {name1:>6}\n")

    print(f"Mapping written to {outfile}")
    print(f"Matched atoms: {len(mapping)}, Unmatched: {len(unmatched)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python create_atomtype_mapping.py file1.pdb file2.pdb [output.txt]")
        sys.exit(1)
    pdb1_file = sys.argv[1]
    pdb2_file = sys.argv[2]
    outfile = sys.argv[3] if len(sys.argv) > 3 else "atom_mapping.txt"
    main(pdb1_file, pdb2_file, outfile)
