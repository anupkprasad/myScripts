from pyrosetta import *
from pyrosetta.teaching import *

# # Initialize PyRosetta with params file
# init(extra_options='-extra_res_fa /home/anup/myScripts/machineLearning/pyrosetta/ATP.params')
init()
safsda
path = "/mnt/mydrive/pseudokinase/ephb6_human_amsa_atemp/pdbs/"

delta_G_bind_list = []
for i in range(0, 15):
    # Load the PDB file
    pdb_file = f"{path}{i}.pdb"
    pose = pose_from_pdb(pdb_file)

    # Save the fixed PDB with hydrogens
    fixed_pdb_file = f"{path}{i}_H.pdb"
    pose.dump_pdb(fixed_pdb_file)
    print(f"Fixed PDB saved to: {fixed_pdb_file}")
    # Load the fixed PDB with hydrogens
    pose = pose_from_pdb(f"{path}{i}_H.pdb")
    scorefxn = get_fa_scorefxn()
    bound_score = scorefxn(pose)
    # Print chain-residue info for debugging
    for i in range(1, pose.total_residue() + 1):
        res = pose.residue(i)
        print(f"{i}: {res.name()} chain={pose.pdb_info().chain(i)}")

    # Split: assume protein is chain A, ligand is chain X
    protein_pose = pose.split_by_chain(1)  # protein chain A
    ligand_pose = pose.split_by_chain(2)   # ligand chain X

    # Optional: relax or repack before scoring
    # relax.apply(protein_pose)
    # relax.apply(ligand_pose)

    # Score unbound parts
    score_protein = scorefxn(protein_pose)
    score_ligand = scorefxn(ligand_pose)

    # Compute binding energy
    delta_G_bind = bound_score - (score_protein + score_ligand)
    print(f"Binding free energy (REU): {delta_G_bind:.3f}")
    delta_G_bind_list.append(delta_G_bind)




hfds
# Load pose from the fixed PDB
pose = pose_from_pdb(fixed_pdb_file)
#pose.dump_pdb("with_H_pyrosetta.pdb")

# Score full complex
scorefxn = get_fa_scorefxn()
bound_score = scorefxn(pose)

# Print chain-residue info for debugging
for i in range(1, pose.total_residue() + 1):
    res = pose.residue(i)
    print(f"{i}: {res.name()} chain={pose.pdb_info().chain(i)}")

# Split: assume protein is chain A, ligand is chain X
protein_pose = pose.split_by_chain(1)  # protein chain A
ligand_pose = pose.split_by_chain(2)   # ligand chain X

# Optional: relax or repack before scoring
# relax.apply(protein_pose)
# relax.apply(ligand_pose)

# Score unbound parts
score_protein = scorefxn(protein_pose)
score_ligand = scorefxn(ligand_pose)

# Compute binding energy
delta_G_bind = bound_score - (score_protein + score_ligand)
print(f"Binding free energy (REU): {delta_G_bind:.3f}")
