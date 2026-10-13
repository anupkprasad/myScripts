#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-08-01 (Y/M/D) at 12:33
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""

from pyrosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.pose import Pose, append_pose_to_pose
from pyrosetta import init, pose_from_pdb

def CheckComponent(pdbfile):
    """
    Load a PDB into a Pose and print residue classification:
    Protein / Ligand / Metal ion / Other
    """
    init(extra_options="-ignore_unrecognized_res")  # safe init
    pose = pose_from_pdb(pdbfile)

    print(f"=== Classification for {pdbfile} ===")
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i)
        chain = pose.pdb_info().chain(i)
        resnum = pose.pdb_info().number(i)

        # classify residue type
        if res.is_protein():
            moltype = "Protein"
        elif res.is_ligand():
            moltype = "Ligand"
        elif res.is_metal():
            moltype = "Metal ion"
        else:
            moltype = "Other"

        print(f"Residue {i:3d} {res.name():6s}  chain={chain} resnum={resnum}  → {moltype}")

def hydrogenatePose(pdb_file):
    # Load the PDB file
    pose = pose_from_pdb(pdb_file)
    # Save the fixed PDB with hydrogens
    fixed_pdb_file = pdb_file.replace(".pdb", "_H.pdb")
    pose.dump_pdb(fixed_pdb_file)
    print(f"Fixed PDB saved to: {fixed_pdb_file}")
    return


def minimizePose(pdb_file):
    # Load complex (protein + ATP ligand)
    pose = pose_from_pdb(pdb_file)

    # Use ref2015 full-atom score function
    scorefxn = get_fa_scorefxn()

    # Create a MoveMap to allow backbone and sidechains to move
    movemap = MoveMap()
    movemap.set_bb(True)
    movemap.set_chi(True)

    # Setup and run MinMover (gradient-based energy minimization)
    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)
    min_mover.min_type('dfpmin_armijo_nonmonotone')  # robust choice
    min_mover.max_iter(200)

    print("Initial Score:", scorefxn(pose))
    min_mover.apply(pose)
    print("Minimized Score:", scorefxn(pose))

    # Save minimized structure
    minimized_pdb_file = pdb_file.replace(".pdb", "_min.pdb")
    dssp = Dssp(pose)
    dssp.insert_ss_into_pose(pose)
    pose.dump_pdb(minimized_pdb_file)
    
    return pose

def computeBindingEnergy(pdb, ion=None):
    # Load the PDB file
    pose = pose_from_pdb(pdb)

    # Score full complex
    scorefxn = get_fa_scorefxn()
    bound_score = scorefxn(pose)

    # # Print chain-residue info for debugging
    # for i in range(1, pose.total_residue() + 1):
    #     res = pose.residue(i)
    #     print(f"{i}: {res.name()} chain={pose.pdb_info().chain(i)}")

    # Split: assume protein is chain A, ligand is chain L
    if ion is not None:
        if pose.num_chains() < 3:
            raise ValueError("Pose must have at least 3 chains (protein, ligand, ion) for binding energy calculation.")
        protein_pose = pose.split_by_chain(1)  # protein chain A
        ligand = pose.split_by_chain(2)   # ligand chain L
        ligand_pose = Pose()
        ligand_pose.assign(ligand)
        ion_pose = pose.split_by_chain(3)       # ion chain I      
        append_pose_to_pose(ligand_pose, ion_pose)
        if ion == 2:
            ion2_pose = pose.split_by_chain(4)  # ion2 chain M
            append_pose_to_pose(ligand_pose, ion2_pose)
            
    else:
        protein_pose = pose.split_by_chain(1)  # protein chain A
        ligand_pose = pose.split_by_chain(2)   # ligand chain L

    if protein_pose.total_atoms() == 0 or ligand_pose.total_atoms() == 0:
        raise ValueError("Failed to split pose into protein and ligand chains.")
    # Print atom numbers info for debugging
    print(f"Protein atoms: {protein_pose.total_atoms()}")
    print(f"Ligand atoms: {ligand_pose.total_atoms()}")

    # Score unbound parts
    score_protein = scorefxn(protein_pose)
    score_ligand = scorefxn(ligand_pose)

    # Compute binding energy
    delta_G_bind = bound_score - (score_protein + score_ligand)
    print(f"Binding free energy (REU): {delta_G_bind:.3f}")
    
    return delta_G_bind


if __name__ == "__main__":
    # Example usage
    # path = "/path/to/structures/"
    # delta_G_bind = computeBindingEnergy(path + "0_min.pdb")
    pass
