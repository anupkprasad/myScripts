#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-08-22 (Y/M/D) at 18:08
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

import MDAnalysis as mda
from parmed import gromacs
import numpy as np
import argparse

def createTopology(top_file, gro_file, output_file="topol.prmtop"):
    """
    Create an Amber topology file from a GROMACS topology and GRO file.
    Parameters:
    top_file (str): Path to the GROMACS topology file (.top).
    gro_file (str): Path to the GRO file (.gro).
    output_file (str): Path to the output Amber topology file (.prmtop).
    """
    g = gromacs.GromacsTopologyFile(top_file, xyz=gro_file)
    g.save(output_file, format="amber")
    print(f"Converted {top_file} and {gro_file} to {output_file}.")


# path = "/mnt/mydrive/GT_project/ruili/Q8N5D6_34_UDPGlcA/equilibration/"

# u = mda.Universe(path + "topol.prmtop", path + "pdb_ion.gro")
# protein_charge = np.sum(u.select_atoms("protein").charges)
# ligand_charge = np.sum(u.select_atoms("resname UGA").charges)
# ion_charge = np.sum(u.select_atoms("resname MN").charges)

# print("Protein charge:", protein_charge)
# print("Ligand charge:", ligand_charge)
# print("Ions charge:", ion_charge)
# print("Net system charge:", np.sum(u.atoms.charges))
# jf;lg
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create Amber topology from GROMACS files.")
    parser.add_argument("top_file", type=str, help="Path to the GROMACS topology file (.top).")
    parser.add_argument("gro_file", type=str, help="Path to the GRO file (.gro).")
    parser.add_argument("--output_file", type=str, default="topol.prmtop", help="Path to the output Amber topology file (.prmtop).")
    args = parser.parse_args()
    createTopology(args.top_file, args.gro_file, args.output_file)

    # u = mda.Universe(args.output_file, args.gro_file)
    # # Example: total charge of protein, ligand, ions
    # protein_charge = np.sum(u.select_atoms("protein").charges)
    # ligand_charge = np.sum(u.select_atoms("resname UGA").charges)
    # ion_charge = np.sum(u.select_atoms("resname MN").charges)

    # print("Protein charge:", protein_charge)
    # print("Ligand charge:", ligand_charge)
    # print("Ions charge:", ion_charge)
    # print("Net system charge:", np.sum(u.atoms.charges))


