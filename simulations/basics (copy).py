#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-08-22 (Y/M/D) at 18:08
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

import MDAnalysis as mda
from parmed import gromacs
from MDAnalysis.coordinates.GRO import GROWriter
from MDAnalysis.core.universe import Merge
import numpy as np

path = "/mnt/mydrive/GT_project/ruili/Q8N5D6_13_UDPGalNAc_GrTr2/equilibration/"

def mergeGroFiles(gro_files, output_file = "complex.gro"):
    """
    Merge multiple GRO files into a single GRO file.
    Parameters:
    gro_files (list of str): List of paths to GRO files to be merged.
    output_file (str): Path to the output merged GRO file.
    """
    universes = [mda.Universe(gro) for gro in gro_files]
    merged = Merge(*[u.atoms for u in universes])

    with GROWriter(output_file, n_atoms=merged.atoms.n_atoms) as w:
        w.write(merged.atoms)

    print(f"Merged {len(gro_files)} GRO files into {output_file} with {merged.atoms.n_atoms} atoms.")



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





# Input files (adjust filenames)
protein_gro = "protein_processed.gro"
ligand_gro  = "ligand_h_GMX.gro"
ion_gro     = "mn.gro"

# Load as separate Universes
u_protein = mda.Universe(path + protein_gro)
u_ligand  = mda.Universe(path + ligand_gro)
u_ion     = mda.Universe(path + ion_gro)


# Merge Universes (preserve original coordinates)
merged = Merge(u_protein.atoms, u_ligand.atoms, u_ion.atoms)



# Write merged GRO
with GROWriter(path + "prot_lig_mn.gro", n_atoms=merged.atoms.n_atoms) as w:
    w.write(merged.atoms)

print(f"✅ Merged system written to {path}prot_lig_mn.gro with {merged.atoms.n_atoms} atoms")


dfdf



# Load GROMACS topology + coordinates
g = gromacs.GromacsTopologyFile(path + "topol.top",  xyz= path + "pdb_ion.gro")
g.save(path + "structure.prmtop", format="amber")

u = mda.Universe(path + "structure.prmtop", path + "pdb_ion.gro")

# Example: total charge of protein, ligand, ions
protein_charge = np.sum(u.select_atoms("protein").charges)
ligand_charge = np.sum(u.select_atoms("resname LIG").charges)
ion_charge = np.sum(u.select_atoms("resname MN").charges)

print("Protein charge:", protein_charge)
print("Ligand charge:", ligand_charge)
print("Ions charge:", ion_charge)
print("Net system charge:", np.sum(u.atoms.charges))


