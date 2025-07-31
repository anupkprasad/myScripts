from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap

init()


path = "/mnt/mydrive/pseudokinase/ephb6_human_amsa_atemp/pdbs/"
# Load complex (protein + ATP ligand)
pose = pose_from_pdb(path + "8_H.pdb")

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
pose.dump_pdb(path + "minimized_proteinATP.pdb")

# Split complex into protein and ligand
# You may need to adapt chain IDs based on your PDB
protein = pose.split_by_chain(1)
ligand = pose.split_by_chain(2)

# Score individual components
score_protein = scorefxn(protein)
score_ligand = scorefxn(ligand)
score_complex = scorefxn(pose)

# Compute binding energy
delta_G_bind = score_complex - (score_protein + score_ligand)

print(f"\nBinding free energy (REU): {delta_G_bind:.3f}")
