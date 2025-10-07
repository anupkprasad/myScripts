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




if __name__ == "__main__":
    # Example usage
    gro_files = ["protein_processed.gro", "ligand_h_GMX.gro", "ion.gro"]
    mergeGroFiles(gro_files, output_file="prot_lig_ion.gro")

