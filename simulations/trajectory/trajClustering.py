# Python script to cluster a molecular dynamics trajectory using RMSD and hierarchical clustering
import os
import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
from scipy.cluster.hierarchy import linkage, fcluster
import numpy as np

def CalRmsd (u, sel="protein and backbone", pos_ref=None):
    """
    Calculate RMSD for a given selection in a trajectory.
    
    Parameters:
    u (MDAnalysis Universe): The universe containing the trajectory.
    sel (str): Selection string to choose atoms for RMSD calculation.
    
    Returns:
    np.ndarray: RMSD values for each frame.
    """
    # Select protein backbone for alignment and RMSD calculation
    protein = u.select_atoms(sel)

    # Compute RMSD for each frame
    n_frames = len(u.trajectory)
    rmsd_values = np.zeros(n_frames)

    for i in range(n_frames):
        u.trajectory[i]  # Move to frame i
        if pos_ref is None:
            # If a reference structure is not provided
            pos_ref = protein.positions.copy()

        pos_i = protein.positions.copy()
        rmsd_values[i] = rmsd(pos_i, pos_ref, superposition=True)
    return rmsd_values


def CalRMSDMatrix(u, sel="protein and backbone"):
    """
    Calculate RMSD matrix for a given selection in a trajectory.
    
    Parameters:
    u (MDAnalysis Universe): The universe containing the trajectory.
    sel (str): Selection string to choose atoms for RMSD calculation.
    
    Returns:
    np.ndarray: RMSD matrix.
    """
    # Select protein backbone for alignment and RMSD calculation
    protein = u.select_atoms(sel)

    # Compute RMSD matrix
    n_frames = len(u.trajectory)
    rmsd_matrix = np.zeros((n_frames, n_frames))

    for i in range(n_frames):
        u.trajectory[i]  # Move to frame i
        pos_i = protein.positions.copy()
        
        for j in range(i, n_frames):  # Only compute upper triangle
            u.trajectory[j]  # Move to frame j
            pos_j = protein.positions.copy()
            
            rmsd_matrix[i, j] = rmsd(pos_i, pos_j, superposition=True)
            rmsd_matrix[j, i] = rmsd_matrix[i, j]  # Symmetric matrix

    return rmsd_matrix


def cluster_trajectory(u, sel="protein and backbone", threshold= None, clstr_num = None):
    """
    Cluster the trajectory based on RMSD.
    
    Parameters:
    u (MDAnalysis Universe): The universe containing the trajectory.
    sel (str): Selection string to choose atoms for RMSD calculation.
    threshold (float): Distance threshold for clustering.
    
    Returns:
    np.ndarray: Cluster labels for each frame.
    """
    rmsd_matrix = CalRMSDMatrix(u, sel)
    
    # Perform hierarchical clustering
    Z = linkage(rmsd_matrix, method='ward')  # Ward’s method
    if threshold is None and clstr_num is not None:
        # If only number of clusters is specified, use maxclust criterion
        clusters = fcluster(Z, clstr_num, criterion='maxclust')
    elif threshold is not None and clstr_num is None:
        # If only distance threshold is specified, use distance criterion
        clusters = fcluster(Z, threshold, criterion='distance')
    else:
        raise ValueError("Specify either 'threshold' or 'clstr_num', not both.")
    
    print(f"Identified {len(set(clusters))} clusters.")


    # Save each frame in its corresponding cluster directory
    output_dir = os.path.join(path, "cluster_results")
    os.makedirs(output_dir, exist_ok=True)
    for i, cluster_id in enumerate(clusters):
        cluster_folder = os.path.join(output_dir, f"cluster_{cluster_id}")
        os.makedirs(cluster_folder, exist_ok=True)  # Create folder if not exists
        u.trajectory[i]  # Move to the frame
        u.atoms.write(os.path.join(cluster_folder, f"frame_{i}.pdb"))  # Save PDB
    print(f"Saved trajectory frames into cluster folders in '{output_dir}/'.")
    
    return clusters, rmsd_matrix









selection = "(backbone and around 15 resname ADE) or (resname ADE)"  # Selection string for protein backbone
path = "/home/anup/anup_data/temp/ryan/test"
# Load trajectory and topology
u = mda.Universe( path + '/em_hg.tpr', path + '/md_50k.xtc', in_memory=True, in_memory_step=100)
ref_u = mda.Universe(path + '/em_hg.tpr')  # Reference structure for alignment
ref = u.select_atoms(selection).positions.copy()  # Reference selection for RMSD calculation

ref_u1 = mda.Universe(path + '/em_wc.tpr')  # Reference structure for alignment
ref1 = u.select_atoms(selection).positions.copy() 


rmsrefs = rmsd(ref1, ref, superposition=True)





# Calculate RMSD
rmsd_values = CalRmsd(u, sel=selection, pos_ref=ref)

# Plot RMSD
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))
plt.plot(rmsd_values, label="RMSD")
plt.xlabel("Frame")
plt.ylabel("RMSD (Å)")
plt.title("RMSD over Time")
plt.legend()
plt.show()



dfdfsd
# Cluster trajectory
clusters, rmsd_matrix = cluster_trajectory(u, sel="protein and backbone", threshold=None, clstr_num=3) ## and resid 1:258