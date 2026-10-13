#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-06-24 (Y/M/D) at 17:24
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

from pymol import cmd

def mutateProtein(selection='all', mutation='ALA'):
    """
    Mutate all residues in the selection to the specified amino acid.
    Args:
        selection (str): PyMOL selection string. Default is 'all'.
        mutation (str): Three-letter code of the amino acid to mutate to. Default is 'ALA'.
    """
    model = cmd.get_model(selection)
    for atom in model.atom:
        # Keep backbone atoms
        if atom.name not in ['N', 'CA', 'C', 'O']:
            resi = atom.resi
            chain = atom.chain
            # Delete sidechain atoms except backbone
            cmd.remove(f"chain {chain} and resi {resi} and not name N+CA+C+O")
            # Rename residue to ALA
            cmd.alter(f"chain {chain} and resi {resi}", "resn= {mutation}")
    cmd.sort()

cmd.extend('mutateProtein', mutate_all_to_alanine)

