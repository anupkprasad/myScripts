#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2025-07-08 (Y/M/D) at 09:56
@author: Anup K. Prasad
email: anupkprasad121@gmail.com
"""
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np
from MDAnalysis.transformations.positionaveraging import PositionAverager
from MDAnalysis.coordinates.PDB import PDBWriter


def CalRmsf(u, outpath = None, selRMSF = "protein and name CA", selAlign = "protein and backbone and resid 1:258",
             ref: str = None, rmsfTimepoint: list = None):    
    if ref==0:  # 0 means first frame of trajectory is used as reference for alignment
        align.AlignTraj(u, reference=u, ref_frame=0, select=selAlign, in_memory=True).run()
    elif ref=="average":
        average = align.AverageStructure(u, u, ref_frame=0).run()
        alignuni = average.universe
        align.AlignTraj(u, reference=alignuni, select=selAlign, in_memory=True).run()
    elif ref:
        uref = ref
        align.AlignTraj(u, reference=uref, select=selAlign, in_memory=True).run()
    else:
        print("No alignment of trajectory")    


    rmsf_atoms = u.select_atoms(selRMSF)
    R = rms.RMSF(rmsf_atoms).run() ##start=0, stop=10, step=1
    rmsf_results = np.vstack((rmsf_atoms.resids, R.results.rmsf)).T
    if rmsfTimepoint:
        rmsfs_t = []
        avg = 10
        for t in rmsfTimepoint:
            if t == rmsfTimepoint[0]:
                Rt = rms.RMSF(rmsf_atoms).run(start=t, stop=t+avg, step=1)
            else:
                Rt = rms.RMSF(rmsf_atoms).run(start=t-avg, stop=t, step=1)
            rmsfs_t.append(Rt.results.rmsf)
        rmsfs_t = np.array(rmsfs_t)
        filename = "rmsfTimepoints" + str(rmsfTimepoint).replace(" ", "").replace("[", "").replace("]", "").replace(",", "_") + ".csv"
        np.savetxt(outpath + '/' + filename, rmsfs_t, delimiter=',')

    print (f"CA: {len(rmsf_atoms), len(rmsf_results)}")
    if outpath:
        np.savetxt(outpath + '/rmsf.csv', rmsf_results, delimiter=',')
        if ref == 0 or ref is None:
            u.trajectory[-1]
            u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
            protein = u.select_atoms(selRMSF)
            for residue, r_value in zip(protein.residues, rmsf_results[:, 1]):
                residue.atoms.tempfactors = r_value
            u.select_atoms("protein").write(outpath + '/proteinWithBetaLastFrame.pdb') #write the bfactor to the last frame of trajectory for all protein atoms
        elif ref == "average":
            alignuni.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
            protein = alignuni.select_atoms(selRMSF)
            for residue, r_value in zip(protein.residues, rmsf_results[:, 1]):
                residue.atoms.tempfactors = r_value
            AvgCA = alignuni.select_atoms("protein and backbone")
            AvgCA.write(outpath + '/proteinWithBetaAvgStr.pdb') #write the bfactor to the avg structure for CA atoms
        else:
            uref.add_TopologyAttr('tempfactors')
            protein = uref.select_atoms(selRMSF)
            for residue, r_value in zip(protein.residues, rmsf_results[:, 1]):
                residue.atoms.tempfactors = r_value
            uref.select_atoms("protein").write(outpath + '/proteinWithBetaRefStr.pdb') #write the bfactor to reference structure all protein atoms
    return rmsf_results


def avarageAtomPositions (u, sel="protein and backbone"):
    """
    Calculate the average position of atoms in a selection over the trajectory.
    
    Parameters:
    u (MDAnalysis Universe): The universe containing the trajectory.
    sel (str): Selection string to choose atoms for averaging.
    
    Returns:
    np.ndarray: Average positions of selected atoms.
    """
    # Select atoms based on the provided selection string
    atoms = u.select_atoms(sel)
    
    # Initialize an array to accumulate positions
    sum_positions = np.zeros((len(atoms), 3))
    n_frames = 0
    
    # Loop over the trajectory frames
    for ts in u.trajectory:
        sum_positions += atoms.positions
        n_frames += 1
    
    # Compute the average positionref
    mean_positions = sum_positions / n_frames

    avg_universe = mda.Merge(atoms)
    avg_universe.atoms.positions = mean_positions
    avg_universe.atoms.write(path + "/avarageAtomPositions.pdb")
    

############################## Code initialization #####################################

path = "/mnt/mydrive/alpfaFold3/dynamics/DCLK/DCLK_KD_tail/MD/dclk3_amsa_atemp"
selRMSF = "protein and name CA"
MD = True
rep = "/rep1"
rmsfTimepoint=[0, 100, 200, 300]
refuni = mda.Universe('/home/anup/myScripts/referencefiles/KinaseStructures/DCLK3_KD.pdb')



if MD:
    u = mda.Universe(path + '/equilibration/renamed_frm0.pdb', path + rep + '/mdWrap1ns.xtc', in_memory=True, in_memory_step=1)
    path = path + rep
else:
    u = mda.Universe(path + '/frm0.pdb', path + '/traj.trr', in_memory=True, in_memory_step=1)
    
res_i = np.min(u.select_atoms("protein").resids)
selAlign = f"protein and name CA and resid {res_i}:{res_i + 257}"
print(f"Str-Alignament using {selAlign}")
rmsf = CalRmsf(u, outpath = path, selRMSF=selRMSF, selAlign=selAlign, ref= refuni, rmsfTimepoint=rmsfTimepoint)  # ref = 0 means first frame of trajectory is used as reference for alignment

c_alphas = u.select_atoms(selRMSF).resids  # Residue indices for the plot 




##########plotting RMSF results##########
plotYlim = 5
fig, ax = plt.subplots(figsize=(3.5, 2.2), dpi = 600, constrained_layout=True)   
ax.plot(rmsf[:,0], rmsf[:,1], color = "black", label =path.split("/")[-2])
ax.legend(fontsize = 6, ncol=2)
from matplotlib.ticker import (MultipleLocator)
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(MultipleLocator(25))
#ax.yaxis.set_major_locator(MultipleLocator(5))
#ax.yaxis.set_minor_locator(MultipleLocator(2.5))
ax.tick_params(axis='both', which='major', labelsize=8)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.set_xlabel("Residue index", fontsize=10)
ax.set_ylabel("RMSF (\u00c5)", fontsize=10)
fig.savefig(path +  "/RMSF.png", dpi = 600)
ax.set_xlim(min(c_alphas), min(c_alphas)+258)  ## save KD only
ax.set_ylim(0, plotYlim)
fig.savefig(path +  "/RMSF_KD.png", dpi = 600)


############ multiline plot of RMSF results ##########
fig, ax = plt.subplots(figsize=(3.5, 2.2), dpi = 600, constrained_layout=True)
filename = "rmsfTimepoints" + str(rmsfTimepoint).replace(" ", "").replace("[", "").replace("]", "").replace(",", "_") + ".csv"
rmsfs_t = np.loadtxt(path + '/' + filename, delimiter=',')

# Set up color map (e.g., sequential colormap)
import matplotlib.cm as cm
cmap = cm.get_cmap('viridis', len(rmsfs_t))  # or 'plasma', 'cividis', 'inferno', etc.

for i, t in enumerate(rmsfTimepoint):
    ax.plot(c_alphas, rmsfs_t[i], label=f"t={t} ns", linewidth=0.5, color=cmap(i))
ax.legend(fontsize = 6, ncol=2)
from matplotlib.ticker import (MultipleLocator)
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(MultipleLocator(25))
#ax.yaxis.set_major_locator(MultipleLocator(5))
#ax.yaxis.set_minor_locator(MultipleLocator(2.5))
ax.tick_params(axis='both', which='major', labelsize=8)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.set_xlabel("Residue index", fontsize=10)
ax.set_ylabel("RMSF (\u00c5)", fontsize=10)
fig.savefig(path +  "/RMSF_timepoints.png", dpi = 600)
ax.set_xlim(min(c_alphas), min(c_alphas)+258)  ## save KD only
ax.set_ylim(0, plotYlim)
fig.savefig(path +  "/RMSF_KD_timepoints.png", dpi = 600)