from pyrosetta import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap

def clean_vmd_pdb(input_pdb, output_pdb):
    with open(input_pdb) as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 54:
                atom = line[12:16]    # Atom name (columns 13-16)
                altloc = line[16]     # AltLoc (column 17)
                resname = line[17:20] # Residue name (18-20)
                chain = line[21]      # Chain ID (column 22)
                rest = line[22:]      # Remaining line (starting from residue number)

                # Fix altloc if it's a number (like '1') → treat it as empty
                altloc = ' ' if altloc not in ['A', 'B', ' '] else altloc

                # Format properly
                new_line = (
                    line[:12] +
                    f"{atom:<4}" +      # Atom name (left-aligned)
                    f"{altloc}" +
                    f"{resname:>3}" +   # Residue name (right-aligned)
                    f"{chain}" +
                    rest
                )
                f_out.write(new_line)
            else:
                f_out.write(line)







# Initialization and example usage should be in main()
if __name__ == "__main__":
    # init(extra_options='-extra_res_fa /path/to/params.params')
    # path = "/path/to/pdbs/"
    # pose = pose_from_pdb(path + "0.pdb")
    # protein = pose.split_by_chain(1)
    pass
ligand = pose.split_by_chain(2)

# Print chain-residue info for debugging
residue_count = pose.total_residue()
chain_count = pose.conformation().num_chains()
for i in range(1, pose.total_residue() + 1):
    res = pose.residue(i)
    print(f"{i}: {res.name()} chain={pose.pdb_info().chain(i)}")






djflkdj
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
